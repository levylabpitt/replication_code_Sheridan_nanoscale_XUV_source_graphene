!! Copyright (C) 2015 N. Helbig and M. Verstraete
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

!> project out states with proper symmetry for cases which are of symmetry = unknown
subroutine X(modelmb_sym_state)(space, mesh, modelmbparticles, ncombo, young_used, &
  wf, symmetries_satisfied, tproj_1yd, nspindown_out, iyoung_out, norm)
  type(space_t),            intent(in)    :: space
  class(mesh_t),            intent(in)    :: mesh
  type(modelmb_particle_t), intent(in)    :: modelmbparticles
  integer,                  intent(in)    :: ncombo
  integer,                  intent(inout) :: young_used(:) !< (1:ncombo)
  R_TYPE,                   intent(inout) :: wf(:) !< will be antisymmetrized on output
  logical,                  intent(out)   :: symmetries_satisfied
  logical,                  intent(in)    :: tproj_1yd
  integer,                  intent(out)   :: nspindown_out(:) !< (1:modelmbparticles%ntype_of_particle)
  integer,                  intent(out)   :: iyoung_out(:) !< (1:modelmbparticles%ntype_of_particle)
  FLOAT,                    intent(out)   :: norm

  integer :: npptype
  integer :: iyoung
  integer :: itype
  integer :: ikeeppart
  integer :: nspindown, nspinup
  integer :: idiagram_combo

  type(young_t) :: young

  integer, allocatable :: dg_combo_iy(:,:)
  integer, allocatable :: dg_combo_ndown(:,:)
  integer, allocatable :: sym_ok_alltypes(:)

  R_TYPE, allocatable  :: antisymwf(:,:,:)              ! stores progressive steps of antisymmetrized wf
  R_TYPE, allocatable  :: fermicompwf(:,:,:)              ! stores progressive steps of antisymmetrized wf

  type(batch_t) :: wfbatch
  R_TYPE :: wfdotp(1,1)

  PUSH_SUB(X(modelmb_sym_state))

  symmetries_satisfied = .false.

  SAFE_ALLOCATE(sym_ok_alltypes(1:modelmbparticles%ntype_of_particle))

  !set up combinations of young diagrams (1 for each type)

  SAFE_ALLOCATE(dg_combo_ndown(1:modelmbparticles%ntype_of_particle, 1:ncombo))
  dg_combo_ndown = 0
  SAFE_ALLOCATE(dg_combo_iy(1:modelmbparticles%ntype_of_particle, 1:ncombo))
  dg_combo_iy = 0

  ! find list of all Young Diagram combinations for all particle types
  idiagram_combo = 1
  do itype = 1, modelmbparticles%ntype_of_particle
    ikeeppart = modelmbparticles%particles_of_type(1,itype)
    if (modelmbparticles%bosonfermion(ikeeppart) /= 1) then ! 1 is for fermion - might introduce a parameter in modelmb_particles
      dg_combo_ndown(itype, idiagram_combo) = 0
      dg_combo_iy(itype, idiagram_combo) = 1
      idiagram_combo = idiagram_combo + 1
      cycle
    end if
    npptype = modelmbparticles%nparticles_per_type(itype)
    do nspindown = 0, floor(npptype/M_TWO)
      nspinup = npptype - nspindown
      call young_init(young, nspinup, nspindown)
      do iyoung = 1, young%nyoung
        dg_combo_ndown(itype, idiagram_combo) = nspindown
        dg_combo_iy(itype, idiagram_combo) = iyoung
        idiagram_combo = idiagram_combo + 1
      end do
      call young_end (young)
    end do
  end do

  ! NB: this is getting ridiculous - too many copies of the wf. Do not have a suggestion, though.
  !   see other subroutines below for even more copies!!!
  SAFE_ALLOCATE(antisymwf(1:mesh%np,1,1))
  SAFE_ALLOCATE(fermicompwf(1:mesh%np,1,1))
  fermicompwf = M_ZERO

  ! index for combination of Young diagrams for all particle types
  do idiagram_combo = 1, ncombo
    antisymwf(:,1,1) = wf(1:mesh%np)
    ! skip diagram combinations already used in present degenerate subspace
    if (young_used (idiagram_combo) > 0) cycle

    call X(modelmb_sym_state_1diag)(space, mesh, modelmbparticles, dg_combo_ndown(:, idiagram_combo), &
      dg_combo_iy(:, idiagram_combo), antisymwf, sym_ok_alltypes, norm)

    ! test the overall symmetrization (no 0.0 norms for present combination of Young diagrams)
    ! check if all types of particles have been properly symmetrized
    if (sum(sym_ok_alltypes) == modelmbparticles%ntype_of_particle .and. abs(norm) > CNST(1.e-6)) then
      fermicompwf = fermicompwf + antisymwf
      symmetries_satisfied = .true.
      young_used (idiagram_combo) = 1
      ! eventually exit the combo loop
      if (tproj_1yd) then
        nspindown_out = dg_combo_ndown(:, idiagram_combo)
        iyoung_out = dg_combo_iy(:, idiagram_combo)
        exit
      end if
    end if

  end do ! idiagram_combo

  ! if we are projecting on all Fermionic YD, need to renormalize here
  if (symmetries_satisfied) then
    ! Only normalize if symmetries are satisfied, other wise the norm might be zero.
    if (mesh%parallel_in_domains) then
      call batch_init(wfbatch, 1, 1, 1, fermicompwf)
      call X(mesh_batch_dotp_self)(mesh, wfbatch, wfdotp, reduce=.true.)
      norm = TOFLOAT(wfdotp(1,1))
      call wfbatch%end()
    else
      norm = TOFLOAT(sum(R_CONJ(fermicompwf(:,1,1))*fermicompwf(:,1,1)))
      norm = norm * product(mesh%spacing(1:space%dim)) !1/units_out%length**space%dim
    end if

    wf(:) = fermicompwf(:,1,1) / sqrt(norm)
  end if

  SAFE_DEALLOCATE_A(antisymwf)
  SAFE_DEALLOCATE_A(fermicompwf)

  SAFE_DEALLOCATE_A(sym_ok_alltypes)
  SAFE_DEALLOCATE_A(dg_combo_ndown)
  SAFE_DEALLOCATE_A(dg_combo_iy)

  POP_SUB(X(modelmb_sym_state))
end subroutine X(modelmb_sym_state)


! ---------------------------------------------------------
!> project out states for a single combination of Young diagrams (1 diagram for each particle type)
subroutine X(modelmb_sym_state_1diag)(space, mesh, modelmbparticles, nspindown_in, iyoung_in, antisymwf, sym_ok_alltypes, norm)
  type(space_t),            intent(in)    :: space
  class(mesh_t),            intent(in)    :: mesh
  type(modelmb_particle_t), intent(in)    :: modelmbparticles
  integer,                  intent(in)    :: nspindown_in(:) !< (1:modelmbparticles%ntype_of_particle)
  integer,                  intent(in)    :: iyoung_in(:) !< (1:modelmbparticles%ntype_of_particle)
  R_TYPE,                   intent(inout) :: antisymwf(:,:,:) !< will be antisymmetrized on output
  integer,                  intent(out)   :: sym_ok_alltypes(:) !< (1:modelmbparticles%ntype_of_particle)
  FLOAT,                    intent(out)   :: norm

  !local vars
  integer :: ipart1, npptype
  integer :: ikeeppart, itype
  integer :: ndimmb
  integer :: ofst_so_far
  integer :: nspinup

  type(permutations_t) :: perms_up, perms_down
  type(young_t) :: young

  type(batch_t) :: antisymwfbatch

  integer, allocatable :: ofst(:)
  integer, allocatable :: p_of_type_up(:), p_of_type_down(:)
  character(len=500) :: tmpstring

  FLOAT :: normalizer
  R_TYPE :: wfdotp(1,1)

  PUSH_SUB(X(modelmb_sym_state_1diag))

  sym_ok_alltypes = 0

  ndimmb = modelmbparticles%ndim

  normalizer = product(mesh%spacing(1:space%dim)) !1/units_out%length**space%dim

  ! this is in case _none_ of the particles is symmetrized,
  !   then the whole itype loop is skipped
  norm = M_ONE

  tmpstring = ""

  ofst_so_far = 0
  ! for each particle type
  do itype = 1, modelmbparticles%ntype_of_particle

    ! FIXME: for multiple particle types this needs to be fixed.
    ! Also, for inequivalent spin configurations this should vary, and we get
    ! different 1 body densities, no?
    ikeeppart = modelmbparticles%particles_of_type(1,itype)

    ! if the particle is not fermionic, just cycle to next one
    ! FIXME: boson case is not treated yet
    if (modelmbparticles%bosonfermion(ikeeppart) /= 1) then ! 1 = fermion
      sym_ok_alltypes(itype) = 1
      POP_SUB(X(modelmb_sym_state_1diag))
      return
    end if

    npptype = modelmbparticles%nparticles_per_type(itype)

    SAFE_ALLOCATE(ofst(1:npptype))
    do ipart1 = 1, npptype
      ofst(ipart1) = ofst_so_far + (ipart1 - 1) * ndimmb
    end do
    ofst_so_far = ofst_so_far + npptype*ndimmb

    ! note: use of spin nomenclature is just for visualization, no real spin
    ! here.
    nspinup = npptype - nspindown_in(itype)

    call permutations_init(nspinup,perms_up)
    call permutations_init(nspindown_in(itype), perms_down)

    ! generate all Young diagrams, decorated, for this distribution of up
    ! and downs
    ! we will only use the diagram indexed iyoung_in(itype) - this is necessary for several particle types :
    ! we need 1 diagram per type.
    call young_init (young, nspinup, nspindown_in(itype))

    ! allocate pointers to pairs of up and down spins in present Young diagram
    SAFE_ALLOCATE(p_of_type_up(1:nspindown_in(itype)))
    SAFE_ALLOCATE(p_of_type_down(1:nspindown_in(itype)))
    p_of_type_up   = modelmbparticles%particles_of_type(young%young_up  (1:nspindown_in(itype), iyoung_in(itype)), itype)
    p_of_type_down = modelmbparticles%particles_of_type(young%young_down(1:nspindown_in(itype), iyoung_in(itype)), itype)

    ! symmetrize pairs of up/down spins in the YD
    call X(modelmb_sym_updown) (ndimmb, npptype, &
      ofst, nspindown_in(itype), p_of_type_up, p_of_type_down, mesh, normalizer, antisymwf)

    SAFE_DEALLOCATE_A(p_of_type_up)
    SAFE_DEALLOCATE_A(p_of_type_down)

    ! antisymmetrize up spins amongst themselves
    call X(modelmb_antisym_1spin) (nspinup,             perms_up, ndimmb, npptype, ofst, &
      young%young_up(:,iyoung_in(itype)), mesh, normalizer, antisymwf)
    ! antisymmetrize down spins amongst themselves
    call X(modelmb_antisym_1spin) (nspindown_in(itype), perms_down, ndimmb, npptype, ofst, &
      young%young_down(:,iyoung_in(itype)), mesh, normalizer, antisymwf)

    call permutations_end(perms_up)
    call permutations_end(perms_down)
    call young_end (young)

    SAFE_DEALLOCATE_A(ofst)

    if (mesh%parallel_in_domains) then
      call batch_init(antisymwfbatch, 1, 1, 1, antisymwf)
      call X(mesh_batch_dotp_self)(mesh, antisymwfbatch, wfdotp, reduce=.true.)
      norm = TOFLOAT(wfdotp(1,1))
      call antisymwfbatch%end()
    else
      norm = TOFLOAT(sum(R_CONJ(antisymwf(:,1,1))*antisymwf(:,1,1))) * normalizer
    end if

    if (abs(norm) > CNST(1.e-6)) sym_ok_alltypes(itype) = 1
  end do ! itype

  POP_SUB(X(modelmb_sym_state_1diag))
end subroutine X(modelmb_sym_state_1diag)
! ---------------------------------------------------------


!> input 1 wave function, and symmetrize wrt spin down labeled particles, according to the given young diagrams
subroutine X(modelmb_sym_updown)(ndimmb, npptype, ofst, ndown, p_of_type_up, p_of_type_down, mesh, normalizer, antisymwf)
  integer,      intent(in)    :: ndimmb
  integer,      intent(in)    :: npptype
  integer,      intent(in)    :: ndown
  integer,      intent(in)    :: ofst(1:npptype)
  integer,      intent(in)    :: p_of_type_up(ndown)
  integer,      intent(in)    :: p_of_type_down(ndown)
  class(mesh_t),intent(in)    :: mesh
  FLOAT,        intent(in)    :: normalizer
  R_TYPE,       intent(inout) :: antisymwf(:,:,:)

  integer :: idown, ipart1, ipart2, ip, ix(MAX_DIM), ixp(MAX_DIM)
  integer(i8) :: ipp
  integer(i8), allocatable :: forward_map_exchange(:)
  FLOAT :: norm
  R_TYPE :: wfdotp(1,1)
  R_TYPE, allocatable  :: antisymwf_swap(:,:,:) ! single wf term, with correct permutation of particles
  type(batch_t) :: antisymwfbatch
  logical :: debug_antisym = .false.

  PUSH_SUB(X(modelmb_sym_updown))

  SAFE_ALLOCATE(forward_map_exchange(1:mesh%np))
  SAFE_ALLOCATE(antisymwf_swap(1:mesh%np, 1, 1))

  ! first symmetrize over pairs of particles associated in the present
  ! Young diagram
  do idown = 1, ndown
    ipart1 = p_of_type_down(idown)
    ipart2 = p_of_type_up(idown)


    ! each processor needs the local map of points for send and recv
    do ip = 1, mesh%np
      ! get present position
      call mesh_local_index_to_coords(mesh, ip, ix)

      ! invert coordinates of ipart1 and ipart2
      ixp = ix
      ! permutate the particles ipart1 and its spin up partner
      ixp(ofst(ipart1)+1:ofst(ipart1)+ndimmb) = ix(ofst(ipart2)+1:ofst(ipart2)+ndimmb)
      ixp(ofst(ipart2)+1:ofst(ipart2)+ndimmb) = ix(ofst(ipart1)+1:ofst(ipart1)+ndimmb)

      ! get position of exchanged point
      ipp = mesh_global_index_from_coords(mesh, ixp)
      ASSERT (ipp <= mesh%np_global)
      forward_map_exchange(ip) = ipp
    end do ! ip

    if (mesh%parallel_in_domains) then
      antisymwf_swap = antisymwf
      ! set up batch type for global exchange operation: 1 state, the loop over MB states is outside this routine
      call batch_init(antisymwfbatch, 1, 1, 1, antisymwf_swap)
      call X(mesh_batch_exchange_points) (mesh, antisymwfbatch, forward_map=forward_map_exchange)
      call antisymwfbatch%end()
    else
      antisymwf_swap(:,1,1) = antisymwf(forward_map_exchange(:),1,1)
    end if

    antisymwf = M_HALF * (antisymwf_swap + antisymwf)

  end do ! idown

  if (debug_antisym) then
    if (mesh%parallel_in_domains) then
      call batch_init(antisymwfbatch, 1, 1, 1, antisymwf)
      call X(mesh_batch_dotp_self)(mesh, antisymwfbatch, wfdotp, reduce=.true.)
      norm = TOFLOAT(wfdotp(1,1))
      call antisymwfbatch%end()
    else
      norm = TOFLOAT(sum(R_CONJ(antisymwf(:,1,1))*antisymwf(:,1,1))) * normalizer
    end if
    write (message(1), '(a,I7,a,E20.10)') 'norm of pair-symmetrized-state with ',&
      ndown, ' spins down is ', norm
    call messages_info(1)
  end if
  SAFE_DEALLOCATE_A(forward_map_exchange)
  SAFE_DEALLOCATE_A(antisymwf_swap)

  POP_SUB(X(modelmb_sym_updown))
end subroutine X(modelmb_sym_updown)

! ---------------------------------------------------------

subroutine X(modelmb_antisym_1spin) (n1spin, perms_1spin, ndimmb, npptype, ofst, young_1spin, mesh, normalizer, antisymwf)
  integer,              intent(in)    :: n1spin
  integer,              intent(in)    :: ndimmb
  integer,              intent(in)    :: npptype
  integer,              intent(in)    :: ofst(1:npptype)
  integer,              intent(in)    :: young_1spin(1:n1spin)
  class(mesh_t),        intent(in)    :: mesh
  FLOAT,                intent(in)    :: normalizer
  type(permutations_t), intent(in)    :: perms_1spin
  R_TYPE,               intent(inout) :: antisymwf(:,:,:)

  integer :: iperm_1spin, i1spin
  integer :: ipart1, ipart2, ip, ix(MAX_DIM), ixp(MAX_DIM)
  R_TYPE :: wfdotp(1,1)
  FLOAT :: norm
  integer(i8), allocatable :: forward_map_exchange(:)
  R_TYPE, allocatable  :: antisymwf_swap(:,:,:) ! single wf term, with correct permutation of particles
  R_TYPE, allocatable  :: antisymwf_acc(:,:,:)
  type(batch_t) :: antisymwfbatch
  logical :: debug_antisym = .false.

  PUSH_SUB(X(modelmb_antisym_1spin))

  SAFE_ALLOCATE(forward_map_exchange(1:mesh%np))
  SAFE_ALLOCATE(antisymwf_swap(1:mesh%np, 1, 1))
  SAFE_ALLOCATE(antisymwf_acc(1:mesh%np, 1, 1))
  ! for each permutation of particles of this type
  !  antisymmetrize the up labeled spins, amongst themselves
  antisymwf_acc = R_TOTYPE(M_ZERO)
  do iperm_1spin = 1, perms_1spin%npermutations

    do ip = 1, mesh%np
      ! get present position
      call mesh_local_index_to_coords(mesh, ip, ix)
      ! initialize coordinates for all particles
      ixp = ix
      ! permute the particles labeled spin up
      do i1spin = 1, n1spin
        ! get image of ipart1 under permutation iperm1
        ipart1 = young_1spin(i1spin)
        ipart2 = young_1spin(perms_1spin%allpermutations(i1spin,iperm_1spin))
        ixp (ofst(ipart1)+1:ofst(ipart1)+ndimmb) = ix (ofst(ipart2)+1:ofst(ipart2)+ndimmb) ! part1 to 2
      end do
      ! get position of exchanged point
      forward_map_exchange(ip) = mesh_global_index_from_coords(mesh, ixp)
    end do ! ip

    if (mesh%parallel_in_domains) then
      antisymwf_swap=antisymwf
      ! set up batch type for global exchange operation: 1 state, the loop over MB states is outside this routine
      call batch_init (antisymwfbatch, 1, 1, 1, antisymwf_swap)
      call X(mesh_batch_exchange_points) (mesh, antisymwfbatch, forward_map=forward_map_exchange)
      call antisymwfbatch%end()
    else
      antisymwf_swap(:,1,1) = antisymwf(forward_map_exchange(:),1,1)
    end if

    antisymwf_acc = antisymwf_acc + perms_1spin%permsign(iperm_1spin)*antisymwf_swap

  end do ! iperm_1spin

  antisymwf = antisymwf_acc / TOFLOAT(perms_1spin%npermutations)

  ! the following could be removed for production
  if (debug_antisym) then
    if (mesh%parallel_in_domains) then
      call batch_init(antisymwfbatch, 1, 1, 1, antisymwf)
      call X(mesh_batch_dotp_self)(mesh, antisymwfbatch, wfdotp, reduce=.true.)
      norm = TOFLOAT(wfdotp(1,1))
      call antisymwfbatch%end()
    else
      norm = TOFLOAT(sum(R_CONJ(antisymwf(:,1,1))*antisymwf(:,1,1))) * normalizer
    end if
    write (message(1), '(a,I7,a,E20.10)') 'norm of 1spin-antisym+pairsym-state with ',&
      n1spin, ' spins of present type is ', norm
    call messages_info(1)
  end if
  SAFE_DEALLOCATE_A(forward_map_exchange)
  SAFE_DEALLOCATE_A(antisymwf_swap)
  SAFE_DEALLOCATE_A(antisymwf_acc)

  POP_SUB(X(modelmb_antisym_1spin))
end subroutine X(modelmb_antisym_1spin)

! ---------------------------------------------------------
!
!> routine to loop over projection of all states wrt fermionic or bosonic character
!
subroutine X(modelmb_sym_all_states)(space, mesh, st)
  type(space_t),          intent(in)    :: space
  class(mesh_t),          intent(in)    :: mesh
  type(states_elec_t),    intent(inout) :: st

  integer :: mm, itype
  integer :: tdrun
  integer :: ncombo
  integer, allocatable :: ndiagrams(:)
  integer, allocatable :: young_used(:)
  integer, allocatable :: young_used_save(:)
  logical :: symmetries_satisfied, impose_exch_symmetry
  R_TYPE, allocatable :: wf(:)

  PUSH_SUB(X(modelmb_sym_all_states))

  impose_exch_symmetry = .true.
  ! this is ugly, but no other way to tell if we are in td run.
  ! other option is to extract present routine and call explicitly from outside output_all. Dont wanna.
  ! TODO : find a better fix now that we are calling from scf
  tdrun = -1
  if (tdrun > 0) impose_exch_symmetry = .false.

  SAFE_ALLOCATE(wf(1:mesh%np))

  ! treat all particle types
  SAFE_ALLOCATE(ndiagrams(1:st%modelmbparticles%ntype_of_particle))
  ndiagrams = 1
  do itype = 1, st%modelmbparticles%ntype_of_particle
    call young_ndiagrams (st%modelmbparticles%nparticles_per_type(itype), ndiagrams(itype))
  end do

  ncombo = product(ndiagrams)

  SAFE_ALLOCATE(young_used(1:ncombo))
  SAFE_ALLOCATE(young_used_save(1:ncombo))
  young_used = 0
  young_used_save = 0

  do mm = 1, st%nst
    call states_elec_get_state(st, mesh, 1, mm, 1, wf)

    symmetries_satisfied = .true.
    if (impose_exch_symmetry) then
      if (mm > 1) then
        ! if eigenval is not degenerate reset young_used
        if (abs(st%eigenval(mm,1) - st%eigenval(mm-1,1)) > CNST(1.e-5)) then
          young_used_save = 0
        end if
      end if

      call X(modelmb_sym_state)(space, mesh, st%modelmbparticles, ncombo, young_used_save, &
        wf, symmetries_satisfied, .true., &
        st%mmb_nspindown(:,mm), st%mmb_iyoung(:,mm), st%mmb_proj(mm))

      young_used = 0
      call X(modelmb_sym_state)(space, mesh, st%modelmbparticles, ncombo, young_used, &
        wf, symmetries_satisfied, .true., &
        st%mmb_nspindown(:,mm), st%mmb_iyoung(:,mm), st%mmb_proj(mm))
    end if

    ! push back the projected state - this may overwrite with a bunch of 0s if we found a bosonic state...
    call states_elec_set_state(st, mesh, 1, mm, 1, wf)

  end do

  SAFE_DEALLOCATE_A(ndiagrams)
  SAFE_DEALLOCATE_A(young_used)
  SAFE_DEALLOCATE_A(young_used_save)


  SAFE_DEALLOCATE_A(wf)

  POP_SUB(X(modelmb_sym_all_states))

end subroutine X(modelmb_sym_all_states)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
