!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2009 X. Andrade
!! Copyright (C) 2020 M. Oliveira
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

#include "global.h"

module electrons_oct_m
  use accel_oct_m
  use absorbing_boundaries_oct_m
  use algorithm_oct_m
  use calc_mode_par_oct_m
  use clock_oct_m
  use debug_oct_m
  use density_oct_m
  use elf_oct_m
  use energy_calc_oct_m
  use gauge_field_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_elec_base_oct_m
  use interaction_oct_m
  use interaction_partner_oct_m
  use ion_dynamics_oct_m
  use ions_oct_m
  use kick_oct_m
  use kpoints_oct_m
  use lattice_vectors_oct_m
  use lasers_oct_m
  use lda_u_oct_m
  use loct_oct_m
  use mesh_oct_m
  use messages_oct_m
  use modelmb_particles_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use pes_oct_m
  use poisson_oct_m
  use propagator_oct_m
  use propagator_elec_oct_m
  use propagator_exp_mid_oct_m
  use propagation_ops_elec_oct_m
  use profiling_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use sort_oct_m
  use system_oct_m
  use td_oct_m
  use td_write_oct_m
  use unit_system_oct_m
  use v_ks_oct_m
  use xc_oct_m
  use xc_oep_oct_m
  use xc_interaction_oct_m

  implicit none

  private
  public ::               &
    electrons_t

  type, extends(system_t) :: electrons_t
    ! Components are public by default
    type(ions_t),     pointer    :: ions => NULL()
    type(grid_t)                 :: gr    !< the mesh
    type(states_elec_t)          :: st    !< the states
    type(v_ks_t)                 :: ks    !< the Kohn-Sham potentials
    type(output_t)               :: outp  !< the output
    type(multicomm_t)            :: mc    !< index and domain communicators
    type(hamiltonian_elec_t)     :: hm
    type(td_t)                   :: td

    type(kpoints_t) :: kpoints                   !< the k-points

    logical :: generate_epot

    type(states_elec_t)          :: st_copy  !< copy of the states

    ! At the moment this is not treated as an external potential
    class(lasers_t), pointer :: lasers => null()      !< lasers
    class(gauge_field_t), pointer :: gfield => null()      !< gauge field

    ! List with all the external partners
    ! This will become a list of interactions in the future
    type(partner_list_t) :: ext_partners

    !TODO: have a list of self interactions
    type(xc_interaction_t), pointer   :: xc_interaction => null()
  contains
    procedure :: init_interaction => electrons_init_interaction
    procedure :: init_parallelization => electrons_init_parallelization
    procedure :: initial_conditions => electrons_initial_conditions
    procedure :: do_td_operation => electrons_do_td_operation
    procedure :: is_tolerance_reached => electrons_is_tolerance_reached
    procedure :: update_quantity => electrons_update_quantity
    procedure :: update_exposed_quantity => electrons_update_exposed_quantity
    procedure :: init_interaction_as_partner => electrons_init_interaction_as_partner
    procedure :: copy_quantities_to_interaction => electrons_copy_quantities_to_interaction
    procedure :: output_start => electrons_output_start
    procedure :: output_write => electrons_output_write
    procedure :: output_finish => electrons_output_finish
    procedure :: process_is_slave  => electrons_process_is_slave
    procedure :: exec_end_of_timestep_tasks => electrons_exec_end_of_timestep_tasks
    procedure :: restart_write_data => electrons_restart_write_data
    procedure :: restart_read_data => electrons_restart_read_data
    procedure :: update_kinetic_energy => electrons_update_kinetic_energy
    final :: electrons_finalize
  end type electrons_t

  interface electrons_t
    procedure electrons_constructor
  end interface electrons_t

contains

  !----------------------------------------------------------
  function electrons_constructor(namespace, generate_epot) result(sys)
    class(electrons_t), pointer    :: sys
    type(namespace_t),  intent(in) :: namespace
    logical,  optional, intent(in) :: generate_epot

    integer :: iatom
    type(lattice_vectors_t) :: latt_inp
    type(profile_t), save :: prof

    PUSH_SUB(electrons_constructor)
    call profiling_in(prof,"ELECTRONS_CONSTRUCTOR")

    SAFE_ALLOCATE(sys)

    sys%namespace = namespace

    call messages_obsolete_variable(sys%namespace, 'SystemName')

    call space_init(sys%space, sys%namespace)
    call sys%space%write_info(sys%namespace)
    if (sys%space%has_mixed_periodicity()) then
      call messages_experimental('Support for mixed periodicity systems')
    end if

    sys%ions => ions_t(sys%namespace, latt_inp=latt_inp)

    call grid_init_stage_1(sys%gr, sys%namespace, sys%space, sys%ions%symm, latt_inp, sys%ions%natoms, sys%ions%pos)

    if (sys%space%is_periodic()) then
      call sys%ions%latt%write_info(sys%namespace)
    end if

    ! Sanity check for atomic coordinates
    do iatom = 1, sys%ions%natoms
      if (.not. sys%gr%box%contains_point(sys%ions%pos(:, iatom))) then
        if (sys%space%periodic_dim /= sys%space%dim) then
          ! FIXME: This could fail for partial periodicity systems
          ! because contains_point is too strict with atoms close to
          ! the upper boundary to the cell.
          write(message(1), '(a,i5,a)') "Atom ", iatom, " is outside the box."
          call messages_warning(1, namespace=sys%namespace)
        end if
      end if
    end do

    ! we need k-points for periodic systems
    call kpoints_init(sys%kpoints, sys%namespace, sys%gr%symm, sys%space%dim, sys%space%periodic_dim, sys%ions%latt)

    call states_elec_init(sys%st, sys%namespace, sys%space, sys%ions%val_charge(), sys%kpoints)
    call sys%st%write_info(sys%namespace)
    ! if independent particles in N dimensions are being used, need to initialize them
    !  after masses are set to 1 in grid_init_stage_1 -> derivatives_init
    call modelmb_copy_masses (sys%st%modelmbparticles, sys%gr%der%masses)
    call elf_init(sys%namespace)

    sys%generate_epot = optional_default(generate_epot, .true.)

    call profiling_out(prof)
    POP_SUB(electrons_constructor)
  end function electrons_constructor

  ! ---------------------------------------------------------
  subroutine electrons_init_interaction(this, interaction)
    class(electrons_t), target, intent(inout) :: this
    class(interaction_t),       intent(inout) :: interaction

    PUSH_SUB(electrons_init_interactions)

    select type (interaction)
    class default
      message(1) = "Trying to initialize an unsupported interaction by the electrons."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(electrons_init_interactions)
  end subroutine electrons_init_interaction

  ! ---------------------------------------------------------
  subroutine electrons_init_parallelization(this, grp)
    class(electrons_t), intent(inout) :: this
    type(mpi_grp_t),    intent(in)    :: grp

    integer(i8) :: index_range(4)
    FLOAT :: mesh_global, mesh_local, wfns
    integer :: idir
    FLOAT :: spiral_q(3), spiral_q_red(3)
    type(block_t) :: blk

    PUSH_SUB(electrons_init_parallelization)

    call mpi_grp_copy(this%grp, grp)

    ! store the ranges for these two indices (serves as initial guess
    ! for parallelization strategy)
    index_range(1) = this%gr%np_global  ! Number of points in mesh
    index_range(2) = this%st%nst             ! Number of states
    index_range(3) = this%st%d%nik           ! Number of k-points
    index_range(4) = 100000                 ! Some large number

    ! create index and domain communicators
    call multicomm_init(this%mc, this%namespace, this%grp, calc_mode_par_parallel_mask(), calc_mode_par_default_parallel_mask(), &
      mpi_world%size, index_range, (/ 5000, 1, 1, 1 /))

    call this%ions%partition(this%mc)
    call kpoints_distribute(this%st%d, this%mc)
    call states_elec_distribute_nodes(this%st, this%namespace, this%mc)


    if (parse_is_defined(this%namespace, 'TDMomentumTransfer') .or. &
      parse_is_defined(this%namespace, 'TDReducedMomentumTransfer')) then
      if (parse_block(this%namespace, 'TDMomentumTransfer', blk) == 0) then
        do idir = 1, this%space%dim
          call parse_block_float(blk, 0, idir - 1, spiral_q(idir))
        end do
      else if(parse_block(this%namespace, 'TDReducedMomentumTransfer', blk) == 0) then
        do idir = 1, this%space%dim
          call parse_block_float(blk, 0, idir - 1, spiral_q_red(idir))
        end do
        call kpoints_to_absolute(this%kpoints%latt, spiral_q_red, spiral_q) 
      end if
      call parse_block_end(blk)
      call grid_init_stage_2(this%gr, this%namespace, this%space, this%mc, spiral_q)
    else
      call grid_init_stage_2(this%gr, this%namespace, this%space, this%mc)
    end if

    if (this%st%symmetrize_density) then
      call mesh_check_symmetries(this%gr, this%gr%symm, this%ions%space%periodic_dim)
    end if

    call output_init(this%outp, this%namespace, this%space, this%st, this%st%nst, this%ks)
    call states_elec_densities_init(this%st, this%gr)
    call states_elec_exec_init(this%st, this%namespace, this%mc)

    call v_ks_init(this%ks, this%namespace, this%gr, this%st, this%ions, this%mc, this%space, this%kpoints)
    if (this%ks%theory_level == KOHN_SHAM_DFT .or. this%ks%theory_level == GENERALIZED_KOHN_SHAM_DFT) then
      this%xc_interaction => xc_interaction_t(this)
    end if

    ! Temporary place for the initialization of the lasers
    this%lasers => lasers_t(this%namespace, this%gr, this%space, this%ions%latt)
    if(this%lasers%no_lasers > 0) then
      call this%ext_partners%add(this%lasers)
      call lasers_check_symmetries(this%lasers, this%kpoints)
    else
      deallocate(this%lasers)
    end if

    ! Temporary place for the initialization of the gauge field
    this%gfield => gauge_field_t(this%namespace, this%ions%latt%rcell_volume)
    if(gauge_field_is_used(this%gfield)) then
      call this%ext_partners%add(this%gfield)
      call gauge_field_check_symmetries(this%gfield, this%kpoints)
    else
      deallocate(this%gfield)
    end if

    call hamiltonian_elec_init(this%hm, this%namespace, this%space, this%gr, this%ions, this%ext_partners, &
      this%st, this%ks%theory_level, this%ks%xc, this%mc, this%kpoints, &
      need_exchange = output_need_exchange(this%outp) .or. this%ks%oep%level /= OEP_LEVEL_NONE)

    if (this%hm%pcm%run_pcm .and. this%mc%par_strategy /= P_STRATEGY_SERIAL .and. this%mc%par_strategy /= P_STRATEGY_STATES) then
      call messages_experimental('Parallel in domain calculations with PCM')
    end if

    ! Print memory requirements
    call messages_print_stress(msg='Approximate memory requirements', namespace=this%namespace)

    mesh_global = mesh_global_memory(this%gr)
    mesh_local  = mesh_local_memory(this%gr)

    call messages_write('Mesh')
    call messages_new_line()
    call messages_write('  global  :')
    call messages_write(mesh_global, units = unit_megabytes, fmt = '(f10.1)')
    call messages_new_line()
    call messages_write('  local   :')
    call messages_write(mesh_local, units = unit_megabytes, fmt = '(f10.1)')
    call messages_new_line()
    call messages_write('  total   :')
    call messages_write(mesh_global + mesh_local, units = unit_megabytes, fmt = '(f10.1)')
    call messages_new_line()
    call messages_info()

    wfns = states_elec_wfns_memory(this%st, this%gr)
    call messages_write('States')
    call messages_new_line()
    call messages_write('  real    :')
    call messages_write(wfns, units = unit_megabytes, fmt = '(f10.1)')
    call messages_write(' (par_kpoints + par_states + par_domains)')
    call messages_new_line()
    call messages_write('  complex :')
    call messages_write(2.0_8*wfns, units = unit_megabytes, fmt = '(f10.1)')
    call messages_write(' (par_kpoints + par_states + par_domains)')
    call messages_new_line()
    call messages_info()

    call messages_print_stress(namespace=this%namespace)

    if (this%generate_epot) then
      message(1) = "Info: Generating external potential"
      call messages_info(1, namespace=this%namespace)
      call hamiltonian_elec_epot_generate(this%hm, this%namespace, this%space, this%gr, this%ions, &
        this%ext_partners, this%st)
      message(1) = "      done."
      call messages_info(1, namespace=this%namespace)
    end if

    if (this%ks%theory_level /= INDEPENDENT_PARTICLES) then
      call poisson_async_init(this%hm%psolver, this%mc)
      ! slave nodes do not call the calculation routine
      if (multicomm_is_slave(this%mc)) then
        !for the moment we only have one type of slave
        call poisson_slave_work(this%hm%psolver, this%namespace)
      end if
    end if


    POP_SUB(electrons_init_parallelization)
  end subroutine electrons_init_parallelization

  ! ---------------------------------------------------------
  subroutine electrons_initial_conditions(this)
    class(electrons_t), intent(inout) :: this

    logical :: fromScratch

    PUSH_SUB(electrons_initial_conditions)

    call td_init(this%td, this%namespace, this%space, this%gr, this%ions, this%st, this%ks, &
      this%hm, this%ext_partners, this%outp)
    ! TODO: implement restarting for electron system
    fromScratch = .true.
    call td_init_run(this%td, this%namespace, this%mc, this%gr, this%ions, this%st, this%ks, this%hm, &
      this%ext_partners, this%outp, this%space, fromScratch)
    this%td%iter = this%td%iter - 1

    POP_SUB(electrons_initial_conditions)
  end subroutine electrons_initial_conditions

  ! ---------------------------------------------------------
  subroutine electrons_do_td_operation(this, operation)
    class(electrons_t),             intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    logical :: update_energy_

    PUSH_SUB(electrons_do_td_operation)

    update_energy_ = .true.

    ! kick at t > 0 still missing!

    select case (operation%id)
    case (STORE_CURRENT_STATUS)
      ! store states at t
      call states_elec_copy(this%st_copy, this%st)

    case (EXPMID_PREDICT_DT_2)
      ! predict states at t + dt/2 from states at t
      call propagation_ops_elec_fuse_density_exp_apply(this%td%tr%te, this%namespace, &
        this%st, this%gr, this%hm, this%prop%dt*M_HALF)

    case (EXPMID_PREDICT_DT)
      ! predict states at t + dt from states at t, using H at t + dt/2
      ! restore states from t to st object
      call states_elec_end(this%st)
      call states_elec_copy(this%st, this%st_copy)
      call states_elec_end(this%st_copy)
      call propagation_ops_elec_fuse_density_exp_apply(this%td%tr%te, this%namespace, &
        this%st, this%gr, this%hm, this%prop%dt)
      this%td%iter = this%td%iter + 1

    case (UPDATE_HAMILTONIAN)
      ! get potential from the updated density
      call v_ks_calc(this%ks, this%namespace, this%space, this%hm, this%st, this%ions, this%ext_partners, &
        calc_eigenval = update_energy_, time = abs(this%prop%clock%time()), calc_energy = update_energy_)
      if (update_energy_) then
        call energy_calc_total(this%namespace, this%space, this%hm, this%gr, this%st, this%ext_partners, iunit = -1)
      end if
      ! update the occupation matrices
      call lda_u_update_occ_matrices(this%hm%lda_u, this%namespace, this%gr, &
        this%st, this%hm%hm_base, this%hm%energy)

    case (EXPMID_START)
    case (EXPMID_FINISH)
    case default
      message(1) = "Unsupported TD operation."
      write(message(2), '(A,A,A)') trim(operation%id), ": ", trim(operation%label)
      call messages_fatal(2, namespace=this%namespace)
    end select


    POP_SUB(electrons_do_td_operation)
  end subroutine electrons_do_td_operation

  ! ---------------------------------------------------------
  logical function electrons_is_tolerance_reached(this, tol) result(converged)
    class(electrons_t), intent(in) :: this
    FLOAT,              intent(in) :: tol

    PUSH_SUB(electrons_is_tolerance_reached)

    converged = .false.

    POP_SUB(electrons_is_tolerance_reached)
  end function electrons_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine electrons_update_quantity(this, iq)
    class(electrons_t),   intent(inout) :: this
    integer,              intent(in)    :: iq

    PUSH_SUB(electrons_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(electrons_update_quantity)
  end subroutine electrons_update_quantity

  ! ---------------------------------------------------------
  subroutine electrons_update_exposed_quantity(partner, iq)
    class(electrons_t), intent(inout) :: partner
    integer,            intent(in)    :: iq

    PUSH_SUB(electrons_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(electrons_update_exposed_quantity)
  end subroutine electrons_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine electrons_init_interaction_as_partner(partner, interaction)
    class(electrons_t),   intent(in)    :: partner
    class(interaction_t), intent(inout) :: interaction

    PUSH_SUB(electrons_init_interaction_as_partner)

    select type (interaction)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(electrons_init_interaction_as_partner)
  end subroutine electrons_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine electrons_copy_quantities_to_interaction(partner, interaction)
    class(electrons_t),   intent(inout) :: partner
    class(interaction_t), intent(inout) :: interaction

    PUSH_SUB(electrons_copy_quantities_to_interaction)

    select type (interaction)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(electrons_copy_quantities_to_interaction)
  end subroutine electrons_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine electrons_output_start(this)
    class(electrons_t), intent(inout) :: this

    PUSH_SUB(electrons_output_start)

    POP_SUB(electrons_output_start)
  end subroutine electrons_output_start

  ! ---------------------------------------------------------
  subroutine electrons_output_write(this)
    class(electrons_t), intent(inout) :: this

    integer :: iter, scsteps
    logical :: fromScratch
    FLOAT :: etime
    logical :: stopping

    PUSH_SUB(electrons_output_write)

    iter = this%td%iter
    scsteps = 1
    fromScratch = .true.
    stopping = .false.
    etime = loct_clock()

    call td_write_iter(this%td%write_handler, this%namespace, this%space, this%outp, this%gr, &
      this%st, this%hm, this%ions, this%ext_partners, this%hm%kick, this%ks, this%td%dt, iter)

    ! write down data
    call td_check_point(this%td, this%namespace, this%mc, this%gr, this%ions, &
      this%st, this%ks, this%hm, this%ext_partners, this%outp, this%space, iter, &
      scsteps, etime, stopping, fromScratch)

    POP_SUB(electrons_output_write)
  end subroutine electrons_output_write

  ! ---------------------------------------------------------
  subroutine electrons_output_finish(this)
    class(electrons_t), intent(inout) :: this

    PUSH_SUB(electrons_output_finish)

    POP_SUB(electrons_output_finish)
  end subroutine electrons_output_finish

  ! ---------------------------------------------------------
  logical function electrons_process_is_slave(this) result(is_slave)
    class(electrons_t), intent(in) :: this

    PUSH_SUB(electrons_process_is_slave)

    is_slave = multicomm_is_slave(this%mc)

    POP_SUB(electrons_process_is_slave)
  end function electrons_process_is_slave

  ! ---------------------------------------------------------
  subroutine electrons_exec_end_of_timestep_tasks(this)
    class(electrons_t), intent(inout) :: this
    logical :: stopping

    PUSH_SUB(electrons_exec_end_of_timestep_tasks)

    stopping = .false.

    !Apply mask absorbing boundaries
    if (this%hm%abs_boundaries%abtype == MASK_ABSORBING) call zvmask(this%gr, this%hm, this%st)

    !Photoelectron stuff
    if (this%td%pesv%calc_spm .or. this%td%pesv%calc_mask .or. this%td%pesv%calc_flux) then
      call pes_calc(this%td%pesv, this%namespace, this%space, this%gr, this%st, &
        this%td%dt, this%td%iter, this%gr%der, this%hm%kpoints, this%ext_partners, stopping)
    end if

    POP_SUB(electrons_exec_end_of_timestep_tasks)
  end subroutine electrons_exec_end_of_timestep_tasks

  ! ---------------------------------------------------------
  subroutine electrons_restart_write_data(this)
    class(electrons_t), intent(inout) :: this

    PUSH_SUB(electrons_restart_write_data)

    message(1) = "Electron system "//trim(this%namespace%get())//" cannot write restart data."
    call messages_warning(1, namespace=this%namespace)

    POP_SUB(electrons_restart_write_data)
  end subroutine electrons_restart_write_data

  ! ---------------------------------------------------------
  ! this function returns true if restart data could be read
  logical function electrons_restart_read_data(this)
    class(electrons_t), intent(inout) :: this

    PUSH_SUB(electrons_restart_read_data)

    ! restarting not yet supported
    electrons_restart_read_data = .false.

    POP_SUB(electrons_restart_read_data)
  end function electrons_restart_read_data

  !----------------------------------------------------------
  subroutine electrons_update_kinetic_energy(this)
    class(electrons_t), intent(inout) :: this

    PUSH_SUB(electrons_update_kinetic_energy)

    if (states_are_real(this%st)) then
      this%kinetic_energy = denergy_calc_electronic(this%namespace, this%hm, this%gr%der, this%st, terms = TERM_KINETIC)
    else
      this%kinetic_energy = zenergy_calc_electronic(this%namespace, this%hm, this%gr%der, this%st, terms = TERM_KINETIC)
    end if

    POP_SUB(electrons_update_kinetic_energy)

  end subroutine electrons_update_kinetic_energy


  !----------------------------------------------------------
  subroutine electrons_finalize(sys)
    type(electrons_t), intent(inout) :: sys

    type(partner_iterator_t) :: iter
    class(interaction_partner_t), pointer :: partner

    PUSH_SUB(electrons_finalize)

    if (associated(sys%prop)) then
      call td_end_run(sys%td, sys%st, sys%hm)
      call td_end(sys%td)
    end if

    if (sys%ks%theory_level /= INDEPENDENT_PARTICLES) then
      call poisson_async_end(sys%hm%psolver, sys%mc)
    end if

    call iter%start(sys%ext_partners)
    do while (iter%has_next())
      partner => iter%get_next()
      SAFE_DEALLOCATE_P(partner)
    end do
    call sys%ext_partners%empty()

    SAFE_DEALLOCATE_P(sys%xc_interaction)

    call hamiltonian_elec_end(sys%hm)

    nullify(sys%gfield)
    nullify(sys%lasers)

    call multicomm_end(sys%mc)

    call v_ks_end(sys%ks)

    call states_elec_end(sys%st)

    SAFE_DEALLOCATE_P(sys%ions)

    call kpoints_end(sys%kpoints)

    call grid_end(sys%gr)

    call system_end(sys)

    POP_SUB(electrons_finalize)
  end subroutine electrons_finalize

end module electrons_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
