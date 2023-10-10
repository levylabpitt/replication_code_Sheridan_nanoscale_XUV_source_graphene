!! Copyright (C) 2019 R. Jestaedt, H. Appel, F. Bonafe, M. Oliveira, N. Tancogne-Dejean
!! Copyright (C) 2022 F. Troisi
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

subroutine X(helmholtz_decomposition_trans_field)(namespace, solver, grid, total_field, transverse_field, &
  apply_boundary, surface_correction)
  type(namespace_t),        intent(in)    :: namespace
  type(poisson_t),          intent(in)    :: solver
  type(grid_t),             intent(in)    :: grid
  R_TYPE,                   intent(inout) :: total_field(:,:)
  R_TYPE,                   intent(out)   :: transverse_field(:,:)
  logical, optional,        intent(in)    :: apply_boundary, surface_correction

  integer             :: idim
  R_TYPE, allocatable :: poisson_density(:,:), poisson_solution(:,:)
  logical             :: apply_boundary_, surface_correction_

  PUSH_SUB(X(helmholtz_decomposition_trans_field))

  SAFE_ALLOCATE(poisson_density(1:grid%np_part, 1:3))
  SAFE_ALLOCATE(poisson_solution(1:grid%np_part, 1:3))

  poisson_solution = M_ZERO
  poisson_density = M_ZERO
  ! if not specified, we apply the boundary condition we computing the curl
  apply_boundary_ = optional_default(apply_boundary, .true.)
  ! if not specified, the surface correction is not applied
  surface_correction_ = optional_default(apply_boundary, .false.)
  ASSERT(.not. surface_correction_)

  ! First of all, we have to compute the curl of the field 
  call X(derivatives_curl)(grid%der, total_field, poisson_density, set_bc = apply_boundary_)

  ! Then we solve poisson equation to compute helmholtz decomposition integral. 
  ! The curl is a vector, so we apply poisson in all directions
  do idim = 1, grid%box%dim
    call X(poisson_solve)(solver, namespace, poisson_solution(:, idim), poisson_density(:, idim))
  end do

  ! TODO: correction surface integral

  ! Add the prefactor 1/(4*pi)
  poisson_solution(:, :) = (M_ONE / (M_FOUR * M_PI)) * poisson_solution(:, :)

  ! Finally we compute the curl again to retrieve the divergence-free term in the helmholtz decomposition
  call X(derivatives_curl)(grid%der, poisson_solution, transverse_field, set_bc = .true.)

  SAFE_DEALLOCATE_A(poisson_solution)
  SAFE_DEALLOCATE_A(poisson_density)

  POP_SUB(X(helmholtz_decomposition_trans_field))
end subroutine X(helmholtz_decomposition_trans_field)

!----------------------------------------------------------
subroutine X(helmholtz_decomposition_long_field)(namespace, solver, grid, total_field, longitudinal_field, apply_boundary, &
 surface_correction)
  type(namespace_t),   intent(in)    :: namespace
  type(poisson_t),     intent(in)    :: solver
  type(grid_t),        intent(in)    :: grid
  R_TYPE,              intent(inout) :: total_field(:,:)
  R_TYPE,              intent(out)   :: longitudinal_field(:,:)
  logical, optional,   intent(in)    :: apply_boundary, surface_correction

  R_TYPE, allocatable :: poisson_density(:), poisson_solution(:)
  logical             :: apply_boundary_, surface_correction_

  PUSH_SUB(X(helmholtz_decomposition_long_field))

  SAFE_ALLOCATE(poisson_density(1:grid%np_part))
  SAFE_ALLOCATE(poisson_solution(1:grid%np_part))

  poisson_density = M_ZERO
  poisson_solution = M_ZERO
  ! if not specified, we apply the boundary condition we computing the curl
  apply_boundary_ = optional_default(apply_boundary, .true.)
  ! if not specified, the surface correction is not applied
  surface_correction_ = optional_default(apply_boundary, .false.)
  ASSERT(.not. surface_correction_)

  ! First of all, we have to compute the divergence of the field 
  call X(derivatives_div)(grid%der, total_field, poisson_density, set_bc = apply_boundary_)

  ! Then we solve poisson equation to compute helmholtz decomposition integral. 
  ! The divergence is scalar, so we apply the poisson only in one direction
  call X(poisson_solve)(solver, namespace, poisson_solution(:), poisson_density(:))
  
  ! TODO: correction surface integral

  ! Add the prefactor 1/(4*pi)
  poisson_solution(:) = - (M_ONE / (M_FOUR * M_PI)) * poisson_solution(:)

  ! Finally we compute the gradient to retrieve the curl-free term in the helmholtz decomposition
  call X(derivatives_grad)(grid%der, poisson_solution, longitudinal_field, set_bc = .true.)

  SAFE_DEALLOCATE_A(poisson_density)
  SAFE_DEALLOCATE_A(poisson_solution)

  POP_SUB(X(helmholtz_decomposition_long_field))
end subroutine X(helmholtz_decomposition_long_field)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
