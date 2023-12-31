!! Copyright (C) 2009-2011 X. Andrade, M. Oliveira
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

! ---------------------------------------------------------
!> To create an etsf file one has to do the following:
!!
!! - Calculate the dimensions and the flags with the _dims functions
!!   for all sets of values.
!! - Init the data file with the flags and the dims.
!! - Call the _write functions for all sets.
!! - Close the file.
!!
!! \note to keep things clean, new data MUST be added following this
!! scheme and using functions.
subroutine output_etsf(outp, namespace, space, dir, st, gr, kpoints, ions, iter)
  type(output_t),         intent(in) :: outp
  type(namespace_t),      intent(in) :: namespace
  type(space_t),          intent(in) :: space
  character(len=*),       intent(in) :: dir
  type(states_elec_t),    intent(in) :: st
  type(grid_t),           intent(in) :: gr
  type(kpoints_t),        intent(in) :: kpoints
  type(ions_t),           intent(in) :: ions
  integer,                intent(in) :: iter

  type(cube_t) :: dcube, zcube

#ifdef HAVE_ETSF_IO
  type(cube_function_t) :: cf
  type(fourier_shell_t) :: shell
  integer :: ncid
  logical :: lstat
  type(etsf_io_low_error)  :: error_data
  type(etsf_dims) :: geometry_dims, density_dims, wfs_dims, pw_dims
  type(etsf_groups_flags) :: geometry_flags, density_flags, wfs_flags, pw_flags
#endif

  PUSH_SUB(output_etsf)

#ifndef HAVE_ETSF_IO
  ASSERT(.false.)
#endif

  !Create a cube
  call cube_init(dcube, gr%idx%ll, namespace, space, gr%spacing, &
    gr%coord_system, fft_type=FFT_REAL, dont_optimize = .true.)
  call cube_init_cube_map(dcube, gr)
  call cube_init(zcube, gr%idx%ll, namespace, space, gr%spacing, &
    gr%coord_system, fft_type=FFT_COMPLEX, dont_optimize = .true.)
  call cube_init_cube_map(zcube, gr)

#ifdef HAVE_ETSF_IO

  ! Careful: only root processor in MPI should attempt to create or modify files.
  ! Nonetheless, routines containing MPI calls such as X(mesh_to_cube) must be called by all processors.

  ! geometry
  if (outp%what_now(OPTION__OUTPUT__GEOMETRY, iter) &
    .and. bitand(outp%how(OPTION__OUTPUT__GEOMETRY), OPTION__OUTPUTFORMAT__ETSF) /= 0) then

    if (mpi_grp_is_root(mpi_world)) then
      call output_etsf_geometry_dims(ions, gr%symm, geometry_dims, geometry_flags)

      call output_etsf_file_init(dir//"/geometry-etsf.nc", "Crystallographic_data file", &
        geometry_dims, geometry_flags, ncid, namespace)

      call output_etsf_geometry_write(ions, gr%symm, ncid, namespace, gr%box)

      call etsf_io_low_close(ncid, lstat, error_data = error_data)
      if (.not. lstat) call output_etsf_error(error_data, namespace)
    end if
  end if

  ! density
  if (outp%what_now(OPTION__OUTPUT__DENSITY, iter) &
    .and. bitand(outp%how(OPTION__OUTPUT__DENSITY), OPTION__OUTPUTFORMAT__ETSF) /= 0) then
    call dcube_function_alloc_rs(dcube, cf)

    call output_etsf_geometry_dims(ions, gr%symm, density_dims, density_flags)
    call output_etsf_density_dims(st, dcube, density_dims, density_flags)

    call output_etsf_file_init(dir//"/density-etsf.nc", "Density file", density_dims, &
      density_flags, ncid, namespace)

    call output_etsf_density_write(st, gr, dcube, cf, ncid, namespace)

    if (mpi_grp_is_root(mpi_world)) then
      call output_etsf_geometry_write(ions, gr%symm, ncid, namespace, gr%box)

      call etsf_io_low_close(ncid, lstat, error_data = error_data)
      if (.not. lstat) call output_etsf_error(error_data, namespace)
    end if

    call dcube_function_free_rs(dcube, cf)
  end if

  ! wave-functions
  if (outp%what_now(OPTION__OUTPUT__WFS, iter) &
    .and. bitand(outp%how(OPTION__OUTPUT__WFS), OPTION__OUTPUTFORMAT__ETSF) /= 0) then

    if (st%parallel_in_states) then
      call messages_not_implemented("ETSF_IO real-space wavefunctions output parallel in states", namespace=namespace)
    end if
    if (st%d%kpt%parallel) then
      call messages_not_implemented("ETSF_IO real-space wavefunctions output parallel in k", namespace=namespace)
    end if

    call dcube_function_alloc_rs(dcube, cf)

    call output_etsf_geometry_dims(ions, gr%symm, wfs_dims, wfs_flags)
    call output_etsf_kpoints_dims(kpoints, wfs_dims, wfs_flags)
    call output_etsf_electrons_dims(st, wfs_dims, wfs_flags)
    call output_etsf_wfs_rsp_dims(st, dcube, wfs_dims, wfs_flags)

    call output_etsf_file_init(dir//"/wfs-etsf.nc", "Wavefunctions file", wfs_dims, wfs_flags, ncid, &
      namespace)

    if (mpi_grp_is_root(mpi_world)) then
      call output_etsf_electrons_write(st, ncid, namespace)
      call output_etsf_geometry_write(ions, gr%symm, ncid, namespace, gr%box)
      call output_etsf_kpoints_write(kpoints, space%dim, ncid, namespace)
    end if
    call output_etsf_wfs_rsp_write(st, gr, dcube, cf, ncid, namespace)

    if (mpi_grp_is_root(mpi_world)) then
      call etsf_io_low_close(ncid, lstat, error_data = error_data)
      if (.not. lstat) call output_etsf_error(error_data, namespace)
    end if

    call dcube_function_free_rs(dcube, cf)
  end if

  ! wave-functions in fourier space
  if (outp%what_now(OPTION__OUTPUT__WFS_FOURIER, iter) &
    .and. bitand(outp%how(OPTION__OUTPUT__WFS_FOURIER), OPTION__OUTPUTFORMAT__ETSF) /= 0) then

    if (st%parallel_in_states) then
      call messages_not_implemented("ETSF_IO Fourier-space wavefunctions output parallel in states", namespace=namespace)
    end if
    if (st%d%kpt%parallel) then
      call messages_not_implemented("ETSF_IO Fourier-space wavefunctions output parallel in k", namespace=namespace)
    end if

    call zcube_function_alloc_rs(zcube, cf)
    call cube_function_alloc_fs(zcube, cf)
    call fourier_shell_init(shell, namespace, space, zcube, gr)

    call output_etsf_geometry_dims(ions, gr%symm, pw_dims, pw_flags)
    call output_etsf_kpoints_dims(kpoints, pw_dims, pw_flags)
    call output_etsf_electrons_dims(st, pw_dims, pw_flags)
    call output_etsf_basisdata_dims(pw_flags)
    call output_etsf_wfs_pw_dims(shell, pw_dims, pw_flags)

    call output_etsf_file_init(dir//"/wfs-pw-etsf.nc", "Wavefunctions file", pw_dims, pw_flags, ncid, &
      namespace)

    if (mpi_grp_is_root(mpi_world)) then
      call output_etsf_electrons_write(st, ncid, namespace)
      call output_etsf_geometry_write(ions, gr%symm, ncid, namespace, gr%box)
      call output_etsf_kpoints_write(kpoints, space%dim, ncid, namespace)
      call output_etsf_basisdata_write(gr, shell, ncid, namespace)
    end if
    call output_etsf_wfs_pw_write(st, gr, zcube, cf, shell, ncid, namespace)

    if (mpi_grp_is_root(mpi_world)) then
      call etsf_io_low_close(ncid, lstat, error_data = error_data)
      if (.not. lstat) call output_etsf_error(error_data, namespace)
    end if

    call fourier_shell_end(shell)
    call cube_function_free_fs(zcube, cf)
    call zcube_function_free_rs(zcube, cf)
  end if
#endif

  call cube_end(dcube)
  call cube_end(zcube)

  POP_SUB(output_etsf)
end subroutine output_etsf

#ifdef HAVE_ETSF_IO
! --------------------------------------------------------

subroutine output_etsf_file_init(filename, filetype, dims, flags, ncid, namespace)
  character(len=*),        intent(in)    :: filename
  character(len=*),        intent(in)    :: filetype
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags
  integer,                 intent(out)   :: ncid
  type(namespace_t),       intent(in)    :: namespace

  logical :: lstat
  type(etsf_io_low_error) :: error_data

  if (.not. mpi_grp_is_root(mpi_world)) return

  PUSH_SUB(output_etsf_file_init)

  ! Note: the presence of 'PACKAGE_STRING' here means the size of the file, and hence the
  ! result of the matches in the testsuite, will change when the version name is updated.
  call etsf_io_data_init(filename, flags, dims, filetype, "Created by " // &
    PACKAGE_STRING, lstat, error_data, overwrite = .true., k_dependent = .false.)
  if (.not. lstat) call output_etsf_error(error_data, namespace)

  call etsf_io_low_open_modify(ncid, filename, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data, namespace)

  call etsf_io_low_set_write_mode(ncid, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data, namespace)

  POP_SUB(output_etsf_file_init)
end subroutine output_etsf_file_init

! --------------------------------------------------------

subroutine output_etsf_error(error_data, namespace)
  type(etsf_io_low_error), intent(in) :: error_data
  type(namespace_t),       intent(in) :: namespace

  PUSH_SUB(output_etsf_error)

  call output_etsf_io_low_error_handle(error_data)
  message(1) = "ETSF_IO returned a fatal error. See message above."
  call messages_fatal(1, only_root_writes = .true., namespace=namespace)

  POP_SUB(output_etsf_error)
end subroutine output_etsf_error

! --------------------------------------------------------

subroutine output_etsf_geometry_dims(ions, symm, dims, flags)
  type(ions_t),            intent(in)    :: ions
  type(symmetries_t),      intent(in)    :: symm
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags

  PUSH_SUB(output_etsf_geometry_dims)

  flags%geometry = etsf_geometry_all - etsf_geometry_valence_charges - etsf_geometry_pseudo_types

  dims%number_of_symmetry_operations = symmetries_number(symm)
  dims%number_of_atom_species = ions%nspecies
  dims%number_of_atoms = ions%natoms

  POP_SUB(output_etsf_geometry_dims)
end subroutine output_etsf_geometry_dims

! --------------------------------------------------------

subroutine output_etsf_geometry_write(ions, symm, ncid, namespace, box)
  type(ions_t),           intent(in)    :: ions
  type(symmetries_t),     intent(in)    :: symm
  integer,                intent(in)    :: ncid
  type(namespace_t),      intent(in)    :: namespace
  class(box_t),           intent(in)    :: box

  type(etsf_geometry) :: geometry
  integer :: idir, isymm, ispecies, i, j
  FLOAT :: offset(1:3)
  type(etsf_io_low_error)  :: error_data
  logical :: lstat

  PUSH_SUB(output_etsf_geometry_write)

  ! Primitive vectors
  SAFE_ALLOCATE(geometry%primitive_vectors(1:3, 1:3))
  do idir = 1, ions%space%dim
    geometry%primitive_vectors(1:3, idir) = ions%latt%rlattice(1:3, idir)
  end do
  do idir = ions%space%periodic_dim + 1, ions%space%dim
    geometry%primitive_vectors(1:3, idir) = geometry%primitive_vectors(1:3, idir) * M_TWO &
      * box%bounding_box_l(idir)
  end do

  ! The symmetries
  SAFE_ALLOCATE(geometry%space_group)
  geometry%space_group = symmetries_space_group_number(symm)
  SAFE_ALLOCATE(geometry%reduced_symmetry_matrices(1:3, 1:3, 1:symmetries_number(symm)))
  SAFE_ALLOCATE(geometry%reduced_symmetry_translations(1:3, 1:symmetries_number(symm)))

  do isymm = 1, symmetries_number(symm)
    geometry%reduced_symmetry_matrices(1:3, 1:3, isymm) = symm_op_rotation_matrix_red(symm%ops(isymm))
    geometry%reduced_symmetry_translations(1:3, isymm) = symm_op_translation_vector_red(symm%ops(isymm))
  end do

  ! The species
  SAFE_ALLOCATE(geometry%atomic_numbers(1:ions%nspecies))
  SAFE_ALLOCATE(geometry%chemical_symbols(1:ions%nspecies))
  SAFE_ALLOCATE(geometry%atom_species_names(1:ions%nspecies))

  do ispecies = 1, ions%nspecies
    geometry%atomic_numbers(ispecies) = species_z(ions%species(ispecies))
    geometry%chemical_symbols(ispecies) = trim(species_label(ions%species(ispecies)))
    ! according to the specification atomic_numbers is enough, but
    ! v_sim wants atom_species_name, so we use the label as name
    geometry%atom_species_names(ispecies) = trim(species_label(ions%species(ispecies)))
  end do

  ! The atoms
  SAFE_ALLOCATE(geometry%atom_species(1:ions%natoms))

  do i = 1, ions%natoms
    do j = 1, ions%nspecies
      if (species_z(ions%atom(i)%species) == species_z(ions%species(j))) then
        geometry%atom_species(i) = j
        exit
      end if
    end do
  end do

  ! The coordinates
  SAFE_ALLOCATE(geometry%reduced_atom_positions(1:3, 1:ions%natoms))

  offset = M_ZERO
  offset(1:ions%space%dim) = -M_HALF*sum(ions%latt%rlattice, dim=2)

  do i = 1, ions%natoms
    ! this is only valid if the primitive vectors are along the x, y, and z directions.
    do idir = 1, 3
      geometry%reduced_atom_positions(idir, i) = (ions%pos(idir, i) - offset(idir))/geometry%primitive_vectors(idir, idir)
    end do
  end do

  call etsf_io_geometry_put(ncid, geometry, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data, namespace)

  ! Free the geometry container
  SAFE_DEALLOCATE_P(geometry%primitive_vectors)
  SAFE_DEALLOCATE_P(geometry%reduced_symmetry_matrices)
  SAFE_DEALLOCATE_P(geometry%reduced_symmetry_translations)
  SAFE_DEALLOCATE_P(geometry%space_group)
  SAFE_DEALLOCATE_P(geometry%reduced_atom_positions)
  SAFE_DEALLOCATE_P(geometry%atom_species)
  SAFE_DEALLOCATE_P(geometry%atomic_numbers)
  SAFE_DEALLOCATE_P(geometry%chemical_symbols)
  SAFE_DEALLOCATE_P(geometry%atom_species_names)

  POP_SUB(output_etsf_geometry_write)
end subroutine output_etsf_geometry_write

! --------------------------------------------------------

subroutine output_etsf_kpoints_dims(kpoints, dims, flags)
  type(kpoints_t),         intent(in)    :: kpoints
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags

  PUSH_SUB(output_etsf_kpoints_dims)

  flags%kpoints = etsf_kpoints_red_coord_kpt + etsf_kpoints_kpoint_weights

  dims%number_of_kpoints = kpoints%reduced%npoints

  POP_SUB(output_etsf_kpoints_dims)
end subroutine output_etsf_kpoints_dims

! --------------------------------------------------------

subroutine output_etsf_kpoints_write(kpoints, dim, ncid, namespace)
  type(kpoints_t),        intent(in)    :: kpoints
  integer,                intent(in)    :: dim
  integer,                intent(in)    :: ncid
  type(namespace_t),      intent(in)    :: namespace

  type(etsf_kpoints), target :: etsf_kpts
  integer  :: nkpoints, ikpoint
  type(etsf_io_low_error)  :: error_data
  logical :: lstat

  PUSH_SUB(output_etsf_kpoints_write)

  nkpoints = kpoints%reduced%npoints

  !Create the kpoints container
  SAFE_ALLOCATE(etsf_kpts%reduced_coordinates_of_kpoints(1:3, 1:nkpoints))
  SAFE_ALLOCATE(etsf_kpts%kpoint_weights(1:nkpoints))

  do ikpoint = 1, nkpoints
    etsf_kpts%reduced_coordinates_of_kpoints(1:3, ikpoint) = M_ZERO
    etsf_kpts%reduced_coordinates_of_kpoints(1:dim, ikpoint) = kpoints%reduced%red_point(1:dim, ikpoint)
    etsf_kpts%kpoint_weights(ikpoint) = kpoints%get_weight(ikpoint)
  end do

  call etsf_io_kpoints_put(ncid, etsf_kpts, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data, namespace)

  SAFE_DEALLOCATE_P(etsf_kpts%reduced_coordinates_of_kpoints)
  SAFE_DEALLOCATE_P(etsf_kpts%kpoint_weights)

  POP_SUB(output_etsf_kpoints_write)
end subroutine output_etsf_kpoints_write

! --------------------------------------------------------

subroutine output_etsf_electrons_dims(st, dims, flags)
  type(states_elec_t),     intent(in)    :: st
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags

  PUSH_SUB(output_etsf_electrons_dims)

  flags%electrons = etsf_electrons_eigenvalues + etsf_electrons_occupations + &
    etsf_electrons_number_of_electrons

  !Set the dimensions
  dims%number_of_spins = 1
  if (st%d%ispin == SPIN_POLARIZED) dims%number_of_spins = 2

  dims%max_number_of_states = st%nst
  dims%number_of_spinor_components = st%d%dim

  POP_SUB(output_etsf_electrons_dims)
end subroutine output_etsf_electrons_dims

! --------------------------------------------------------

subroutine output_etsf_electrons_write(st, ncid, namespace)
  type(states_elec_t), intent(in)    :: st
  integer,             intent(in)    :: ncid
  type(namespace_t),   intent(in)    :: namespace

  type(etsf_electrons) :: electrons
  type(etsf_io_low_error)  :: error_data
  logical :: lstat
  integer :: nspin, ist, ik, ispin, nkpoints

  PUSH_SUB(output_etsf_electrons_write)

  SAFE_ALLOCATE(electrons%number_of_electrons)
  electrons%number_of_electrons = int(st%qtot)

  nspin = 1
  if (st%d%ispin == SPIN_POLARIZED) nspin = 2

  nkpoints = st%d%nik/nspin

  SAFE_ALLOCATE(electrons%eigenvalues%data3D(1:st%nst, 1:nkpoints, 1:nspin))
  SAFE_ALLOCATE(electrons%occupations%data3D(1:st%nst, 1:nkpoints, 1:nspin))

  do ist = 1, st%nst
    do ik = 1, nkpoints
      do ispin = 1, nspin
        electrons%eigenvalues%data3D(ist, ik, ispin) = st%eigenval(ist, nspin*(ik-1) + ispin)
        electrons%occupations%data3D(ist, ik, ispin) = st%occ(ist, nspin*(ik-1) + ispin)
      end do
    end do
  end do

  call etsf_io_electrons_put(ncid, electrons, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data, namespace)

  SAFE_DEALLOCATE_P(electrons%number_of_electrons)
  SAFE_DEALLOCATE_P(electrons%eigenvalues%data3D)
  SAFE_DEALLOCATE_P(electrons%occupations%data3D)

  POP_SUB(output_etsf_electrons_write)
end subroutine output_etsf_electrons_write

! --------------------------------------------------------

subroutine output_etsf_density_dims(st, cube, dims, flags)
  type(states_elec_t),     intent(in)    :: st
  type(cube_t),            intent(in)    :: cube
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags

  PUSH_SUB(output_etsf_density_dims)

  flags%main = etsf_main_density

  dims%number_of_components = st%d%nspin
  dims%number_of_grid_points_vector1 = cube%rs_n_global(1)
  dims%number_of_grid_points_vector2 = cube%rs_n_global(2)
  dims%number_of_grid_points_vector3 = cube%rs_n_global(3)
  dims%real_or_complex_density = 1

  POP_SUB(output_etsf_density_dims)
end subroutine output_etsf_density_dims

! --------------------------------------------------------

subroutine output_etsf_density_write(st, mesh, cube, cf, ncid, namespace)
  type(states_elec_t),   intent(in)    :: st
  class(mesh_t),         intent(in)    :: mesh
  type(cube_t),          intent(inout) :: cube
  type(cube_function_t), intent(inout) :: cf
  integer,               intent(in)    :: ncid
  type(namespace_t),     intent(in)    :: namespace

  type(etsf_main) :: main
  type(etsf_io_low_error)  :: error_data
  logical :: lstat
  integer :: ispin, n(3)
  FLOAT, allocatable :: d(:), md(:,:)

  PUSH_SUB(output_etsf_density_write)


  n = cube%rs_n_global
  SAFE_ALLOCATE(main%density%data4D(1:n(1), 1:n(2), 1:n(3), 1:st%d%nspin))

  if (st%d%ispin /= SPINORS) then
    do ispin = 1, st%d%nspin
      call dmesh_to_cube(mesh, st%rho(:, ispin), cube, cf)
      main%density%data4D(1:n(1), 1:n(2), 1:n(3), ispin) = cf%drs(1:n(1), 1:n(2), 1:n(3))
    end do
  else
    SAFE_ALLOCATE(md(1:mesh%np, 1:3))
    SAFE_ALLOCATE(d(1:mesh%np_part))

    d = st%rho(:, 1) + st%rho(:, 2)
    call magnetic_density(mesh, st%d, st%rho, md)

    call dmesh_to_cube(mesh, d, cube, cf)
    main%density%data4D(1:n(1), 1:n(2), 1:n(3), 1) = cf%drs(1:n(1), 1:n(2), 1:n(3))
    do ispin = 1, 3
      call dmesh_to_cube(mesh, md(:, ispin), cube, cf)
      main%density%data4D(1:n(1), 1:n(2), 1:n(3), ispin + 1) = cf%drs(1:n(1), 1:n(2), 1:n(3))
    end do
    SAFE_DEALLOCATE_A(d)
    SAFE_DEALLOCATE_A(md)
  end if

  if (mpi_grp_is_root(mpi_world)) then
    call etsf_io_main_put(ncid, main, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data, namespace)
  end if

  SAFE_DEALLOCATE_P(main%density%data4D)

  POP_SUB(output_etsf_density_write)
end subroutine output_etsf_density_write


! --------------------------------------------------------

subroutine output_etsf_wfs_rsp_dims(st, cube, dims, flags)
  type(states_elec_t),     intent(in)    :: st
  type(cube_t),            intent(in)    :: cube
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags

  PUSH_SUB(output_etsf_wfs_rsp_dims)

  if (states_are_real(st)) then
    dims%real_or_complex_wavefunctions = 1
  else
    dims%real_or_complex_wavefunctions = 2
  end if

  dims%number_of_grid_points_vector1 = cube%rs_n_global(1)
  dims%number_of_grid_points_vector2 = cube%rs_n_global(2)
  dims%number_of_grid_points_vector3 = cube%rs_n_global(3)

  flags%main = etsf_main_wfs_rsp

  POP_SUB(output_etsf_wfs_rsp_dims)
end subroutine output_etsf_wfs_rsp_dims

! --------------------------------------------------------

subroutine output_etsf_wfs_rsp_write(st, mesh, cube, cf, ncid, namespace)
  type(states_elec_t),   intent(in)    :: st
  class(mesh_t),         intent(in)    :: mesh
  type(cube_t),          intent(inout) :: cube
  type(cube_function_t), intent(inout) :: cf
  integer,               intent(in)    :: ncid
  type(namespace_t),     intent(in)    :: namespace

  integer :: ist, ispin, ik, idim, nspin, zdim, nkpoints, n(3)
  type(etsf_main) :: main
  type(etsf_io_low_error)  :: error_data
  logical :: lstat
  REAL_DOUBLE, allocatable, target :: local_wfs(:,:,:,:,:,:,:)
  FLOAT, allocatable :: dpsi(:)
  CMPLX, allocatable :: zpsi(:)

  PUSH_SUB(output_etsf_wfs_rsp_write)

  nspin = 1
  if (st%d%ispin == SPIN_POLARIZED) nspin = 2

  nkpoints = st%d%nik/nspin
  if (states_are_real(st)) then
    zdim = 1
  else
    zdim = 2
  end if

  if (states_are_real(st)) then
    SAFE_ALLOCATE(dpsi(1:mesh%np))
  else
    SAFE_ALLOCATE(zpsi(1:mesh%np))
  end if

  n = cube%rs_n_global
  SAFE_ALLOCATE(local_wfs(1:zdim, 1:n(1), 1:n(2), 1:n(3), 1:st%d%dim, 1:st%nst, 1:st%d%nik))
  do ispin = 1, nspin
    do ik = 1, st%d%nik, nspin
      do ist = 1, st%nst
        do idim = 1, st%d%dim
          if (states_are_real(st)) then
            call states_elec_get_state(st, mesh, idim, ist, ik + ispin - 1, dpsi)

            call dmesh_to_cube(mesh, dpsi, cube, cf)
            local_wfs(1, 1:n(1), 1:n(2), 1:n(3), idim, ist, ik+(ispin-1)*nkpoints) = cf%drs(1:n(1), 1:n(2), 1:n(3))

          else
            call states_elec_get_state(st, mesh, idim, ist, ik + ispin - 1, zpsi)

            call dmesh_to_cube(mesh, TOFLOAT(zpsi), cube, cf)
            local_wfs(1, 1:n(1), 1:n(2), 1:n(3), idim, ist, ik+(ispin-1)*nkpoints) = cf%drs(1:n(1), 1:n(2), 1:n(3))

            call dmesh_to_cube(mesh, aimag(zpsi), cube, cf)
            local_wfs(2, 1:n(1), 1:n(2), 1:n(3), idim, ist, ik+(ispin-1)*nkpoints) = cf%drs(1:n(1), 1:n(2), 1:n(3))

          end if
        end do
      end do
    end do
  end do

  SAFE_DEALLOCATE_A(dpsi)
  SAFE_DEALLOCATE_A(zpsi)

  main%real_space_wavefunctions%data7D => local_wfs

  if (mpi_grp_is_root(mpi_world)) then
    call etsf_io_main_put(ncid, main, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data, namespace)
  end if

  SAFE_DEALLOCATE_A(local_wfs)

  POP_SUB(output_etsf_wfs_rsp_write)
end subroutine output_etsf_wfs_rsp_write

! --------------------------------------------------

subroutine output_etsf_basisdata_dims(flags)
  type(etsf_groups_flags), intent(inout) :: flags

  PUSH_SUB(output_etsf_basisdata_dims)

  flags%basisdata = etsf_basisdata_basis_set + etsf_basisdata_kin_cutoff + &
    etsf_basisdata_red_coord_pw

  POP_SUB(output_etsf_basisdata_dims)

end subroutine output_etsf_basisdata_dims

! --------------------------------------------------------

subroutine output_etsf_basisdata_write(mesh, shell, ncid, namespace)
  class(mesh_t),         intent(in)    :: mesh
  type(fourier_shell_t), intent(in)    :: shell
  integer,               intent(in)    :: ncid
  type(namespace_t),     intent(in)    :: namespace

  type(etsf_basisdata) :: basisdata
  type(etsf_io_low_error) :: error_data
  logical :: lstat
  integer :: ng

  PUSH_SUB(output_etsf_basisdata_write)

  if ((maxval(mesh%spacing(1:3)) - minval(mesh%spacing(1:3))) > CNST(1e-10)) then
    message(1) = 'Cannot generate a ETSF plane-wave wave-functions file,'
    message(2) = 'spacing is not the same for each direction.'
    call messages_fatal(2, namespace=namespace)
  end if

  SAFE_ALLOCATE(basisdata%basis_set)

  basisdata%basis_set = "plane_waves"

  SAFE_ALLOCATE(basisdata%kinetic_energy_cutoff)

  basisdata%kinetic_energy_cutoff = shell%ekin_cutoff

  ng = shell%ngvectors

  SAFE_ALLOCATE(basisdata%reduced_coordinates_of_plane_waves%data2D(1:3, 1:ng))

  basisdata%reduced_coordinates_of_plane_waves%data2D(1:3, 1:ng) = shell%red_gvec(1:3, 1:ng)

  call etsf_io_basisdata_put(ncid, basisdata, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data, namespace)

  SAFE_DEALLOCATE_P(basisdata%basis_set)
  SAFE_DEALLOCATE_P(basisdata%kinetic_energy_cutoff)
  SAFE_DEALLOCATE_P(basisdata%reduced_coordinates_of_plane_waves%data2D)

  POP_SUB(output_etsf_basisdata_write)

end subroutine output_etsf_basisdata_write

! --------------------------------------------------------

subroutine output_etsf_wfs_pw_dims(shell, dims, flags)
  type(fourier_shell_t),   intent(in)    :: shell
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags

  PUSH_SUB(output_etsf_wfs_pw_dims)

  flags%main = etsf_main_wfs_coeff

  dims%max_number_of_coefficients = shell%ngvectors
  dims%real_or_complex_coefficients = 2

  POP_SUB(output_etsf_wfs_pw_dims)

end subroutine output_etsf_wfs_pw_dims

! --------------------------------------------------------

subroutine output_etsf_wfs_pw_write(st, mesh, cube, cf, shell, ncid, namespace)
  type(states_elec_t),   intent(in)    :: st
  class(mesh_t),         intent(in)    :: mesh
  type(cube_t),          intent(inout) :: cube
  type(cube_function_t), intent(inout) :: cf
  type(fourier_shell_t), intent(in)    :: shell
  integer,               intent(in)    :: ncid
  type(namespace_t),     intent(in)    :: namespace

  type(etsf_main) :: main
  type(etsf_io_low_error)  :: error_data
  logical :: lstat
  REAL_DOUBLE, allocatable, target :: local_wfs(:, :, :, :, :, :)
  CMPLX, allocatable :: zpsi(:)
  integer :: nkpoints, nspin, zdim
  integer :: idim, ist, iq, ikpoint, ispin
  integer :: ig, ng, ix, iy, iz

  PUSH_SUB(output_etsf_wfs_pw_write)

  nspin = 1
  if (st%d%ispin == SPIN_POLARIZED) nspin = 2

  nkpoints = st%d%nik/nspin
  zdim = 2

  ng = shell%ngvectors

  !Write the wavefunctions to the file
  SAFE_ALLOCATE(local_wfs(1:zdim, 1:ng, 1:st%d%dim, 1:st%nst, 1:nkpoints, 1:nspin))

  SAFE_ALLOCATE(zpsi(1:mesh%np))

  do iq = 1, st%d%nik
    ispin = st%d%get_spin_index(iq)
    ikpoint = st%d%get_kpoint_index(iq)
    do ist = 1, st%nst
      do idim = 1, st%d%dim

        ! for the moment we treat all functions as complex
        call states_elec_get_state(st, mesh, idim, ist, iq, zpsi)
        call zmesh_to_cube(mesh, zpsi, cube, cf)
        call zcube_function_rs2fs(cube, cf)

        do ig = 1, ng
          ix = shell%coords(1, ig)
          iy = shell%coords(2, ig)
          iz = shell%coords(3, ig)

          local_wfs(1, ig, idim, ist, ikpoint, ispin) = TOFLOAT(cf%fs(ix, iy, iz))
          local_wfs(2, ig, idim, ist, ikpoint, ispin) = aimag(cf%fs(ix, iy, iz))
        end do

      end do
    end do
  end do

  SAFE_DEALLOCATE_A(zpsi)

  main%coefficients_of_wavefunctions%data6D => local_wfs

  if (mpi_grp_is_root(mpi_world)) then
    call etsf_io_main_put(ncid, main, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data, namespace)

    call etsf_io_tools_set_time_reversal_symmetry(ncid, .false., lstat, error_data)
    if (.not. lstat) call output_etsf_error(error_data, namespace)
  end if

  SAFE_DEALLOCATE_A(local_wfs)

  POP_SUB(output_etsf_wfs_pw_write)

end subroutine output_etsf_wfs_pw_write

!> DAS: copied from ETSF_IO 1.0.3 etsf_io_low_level.f90, changed to send output to standard error
!!****m* etsf_io_low_error_group/etsf_io_low_error_handle
!! NAME
!!  etsf_io_low_error_handle
!!
!! FUNCTION
!!  This method can be used to output the informations contained in an error
!!  structure. The output is done on standard output. Write your own method
!!  if custom error handling is required.
!!
!! COPYRIGHT
!!  Copyright (C) 2006, 2007 (Damien Caliste)
!!  This file is distributed under the terms of the
!!  GNU Lesser General Public License, see the COPYING file
!!  or http://www.gnu.org/copyleft/lesser.txt .
!!
!! INPUTS
!!  * error_data <type(etsf_io_low_error)>=informations about an error.
!!
!! SOURCE
subroutine output_etsf_io_low_error_handle(error_data)
  type(etsf_io_low_error), intent(in) :: error_data

  integer :: i

  if (.not. mpi_grp_is_root(mpi_world)) return

  PUSH_SUB(output_etsf_io_low_error_handle)

  ! Error handling
  write(stderr,*)
  write(stderr,*) "    ***"
  write(stderr,*) "    *** ETSF I/O ERROR"
  write(stderr,*) "    ***"
  write(stderr,*) "    *** Backtrace          : ", &
    trim(error_data%backtrace(error_data%backtraceId)), "()"
  do i = error_data%backtraceId - 1, 1, -1
    write(stderr,*) "    ***                      ", trim(error_data%backtrace(i)), "()"
  end do
  write(stderr,*) "    *** Action performed   : ", trim(error_data%access_mode_str), &
    " ", trim(error_data%target_type_str)
  if (trim(error_data%target_name) /= "") then
    write(stderr,*) "    *** Target (name)      : ", trim(error_data%target_name)
  end if
  if (error_data%target_id /= 0) then
    write(stderr,*) "    *** Target (id)        : ", error_data%target_id
  end if
  if (trim(error_data%error_message) /= "") then
    write(stderr,*) "    *** Error message      : ", trim(error_data%error_message)
  end if
  if (error_data%error_id /= nf90_noerr) then
    write(stderr,*) "    *** Error id           : ", error_data%error_id
  end if
  write(stderr,*) "    ***"
  write(stderr,*)

  POP_SUB(output_etsf_io_low_error_handle)

end subroutine output_etsf_io_low_error_handle
!!***

#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
