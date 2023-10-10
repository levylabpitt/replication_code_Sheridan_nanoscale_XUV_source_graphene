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

#include "global.h"

!> This module is intended to contain "only mathematical" functions and procedures to compute the Helmholtz decomposition
!> of a generic field. 
!> Given a generic vector field F, it is possibile to define a curl-free component div(phi) and a divergence free component 
!> curl(psi) such that F = -div(phi) + curl(psi), where
!> phi = (1/4*pi) * (volume_integral(div(F) / |r - r1|) - surface_integral(dot_p(n, F) / |r - r1|))
!> psi = (1/4*pi) * (volume_integral(curl(F) / |r - r1|) - surface_integral(cross_p(n, F) / |r - r1|))

module helmholtz_decomposition_m
  use blas_oct_m
  use debug_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use io_oct_m
  use io_function_oct_m
  use lalg_basic_oct_m
  use loct_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use namespace_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use space_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                                       &
    helmholtz_decomposition_trans_field,          &
    helmholtz_decomposition_long_field,           &
    helmholtz_decomposition_hertzian_dipole_test, &
    helmholtz_decomposition_gaussian_test

  interface helmholtz_decomposition_trans_field
    module procedure dhelmholtz_decomposition_trans_field, zhelmholtz_decomposition_trans_field
  end interface helmholtz_decomposition_trans_field

  interface helmholtz_decomposition_long_field
    module procedure dhelmholtz_decomposition_long_field, zhelmholtz_decomposition_long_field
  end interface helmholtz_decomposition_long_field

contains

  subroutine helmholtz_decomposition_hertzian_dipole_test(namespace, poisson_solver, grid, space)
    type(namespace_t), intent(in) :: namespace
    type(poisson_t),   intent(in) :: poisson_solver
    type(grid_t),      intent(in) :: grid
    type(space_t),     intent(in) :: space

    FLOAT, allocatable :: total_field(:,:), trans_field_exact(:,:), long_field_exact(:,:)
    FLOAT, allocatable :: trans_field_computed(:,:), long_field_computed(:,:), div_rot(:)
    FLOAT              :: lambda, omega, kk, rr, pp, time
    FLOAT              :: xx(space%dim), xx_origin(space%dim), dipole(space%dim)
    integer            :: ip

    PUSH_SUB(helmholtz_decomposition_hertzian_dipole_test)
    ASSERT(space%dim == 3)

    ! Allocate and initialize the fields
    SAFE_ALLOCATE(total_field(1:grid%np_part, 1:grid%box%dim))
    SAFE_ALLOCATE(trans_field_exact(1:grid%np_part, 1:grid%box%dim))
    SAFE_ALLOCATE(long_field_exact(1:grid%np_part, 1:grid%box%dim))
    SAFE_ALLOCATE(trans_field_computed(1:grid%np_part, 1:grid%box%dim))
    SAFE_ALLOCATE(long_field_computed(1:grid%np_part, 1:grid%box%dim))
    SAFE_ALLOCATE(div_rot(1:grid%np_part))
    total_field = M_ZERO
    trans_field_computed = M_ZERO
    trans_field_exact = M_ZERO
    long_field_computed = M_ZERO
    long_field_exact = M_ZERO
    div_rot = M_ZERO
    ! Define the parameters that define the field
    lambda = CNST(1.55)
    omega = M_TWO * M_PI * (P_c / lambda)
    kk = M_TWO * M_PI / lambda
    time = M_ZERO
    pp = cos(omega * time)
    dipole = M_ZERO
    dipole(3) = M_ONE
    xx_origin = M_ZERO
    xx_origin(1) = CNST(3.11)

    ! Define the field to be tested (Hertzian dipole). 
    ! Reference: The Helmholtz Decomposition and the Coulomb Gauge, Kirk T. McDonald, Joseph Henry Laboratories, Princeton (2008)
    ! https://puhep1.princeton.edu/~kirkmcd/examples/helmholtz.pdf 
    do ip = 1, grid%np
      call mesh_r(grid, ip, rr, coords = xx, origin = xx_origin)
      if (rr > CNST(1e-6)) then
        xx = xx / rr
        ! Total field
        total_field(ip, 1) = (kk**2)*pp*(-xx(1)*xx(3))*(cos(kk*rr - omega*time)/rr)
        total_field(ip, 1) = total_field(ip, 1) + &
          pp*(M_THREE*xx(1)*xx(3) - dipole(1))*((cos(kk*rr - omega*time)/(rr**3)) + (kk*sin(kk*rr - omega*time)/(rr**2)))

        total_field(ip, 2) = (kk**2)*pp*(-xx(2)*xx(3))*(cos(kk*rr - omega*time)/rr)
        total_field(ip, 2) = total_field(ip, 2) + &
          pp*(M_THREE*xx(2)*xx(3) - dipole(2))*((cos(kk*rr - omega*time)/(rr**3)) + (kk*sin(kk*rr - omega*time)/(rr**2)))

        total_field(ip, 3) = (kk**2)*pp*((xx(1)**2) + (xx(2)**2))*(cos(kk*rr - omega*time)/rr)
        total_field(ip, 3) = total_field(ip,3) + &
          pp*(M_THREE*(xx(3)**2) - dipole(3))*((cos(kk*rr - omega*time)/(rr**3)) + (kk*sin(kk*rr - omega*time)/(rr**2)))

        ! Longitudinal field
        long_field_exact(ip, 1) = pp*(M_THREE*xx(1)*xx(3) - dipole(1))*(cos(omega*time)/(rr**3))
        long_field_exact(ip, 2) = pp*(M_THREE*xx(2)*xx(3) - dipole(2))*(cos(omega*time)/(rr**3))
        long_field_exact(ip, 3) = pp*(M_THREE*(xx(3)**2) - dipole(3))*(cos(omega*time)/(rr**3))

        ! Transverse field
        trans_field_exact(ip, 1) = (kk**2)*pp*(-xx(1)*xx(3))*(cos(kk*rr - omega*time)/rr)
        trans_field_exact(ip, 1) = trans_field_exact(ip, 1) + pp*(M_THREE*xx(1)*xx(3) - dipole(1)) * &
          (((cos(kk*rr - omega*time) - cos(omega*time))/(rr**3)) + (kk*sin(kk*rr - omega*time)/(rr**2)))

        trans_field_exact(ip, 2) = (kk**2)*pp*(-xx(2)*xx(3))*(cos(kk*rr - omega*time)/rr)
        trans_field_exact(ip, 2) = trans_field_exact(ip, 2) + pp*(M_THREE*xx(2)*xx(3) - dipole(2)) * &
          (((cos(kk*rr - omega*time) - cos(omega*time))/(rr**3)) + (kk*sin(kk*rr - omega*time)/(rr**2)))

        trans_field_exact(ip, 3) = (kk**2)*pp*((xx(1)**2) + (xx(2)**2))*(cos(kk*rr - omega*time)/rr)
        trans_field_exact(ip, 3) = trans_field_exact(ip, 3) + pp*(M_THREE*(xx(3)**2) - dipole(3)) * &
          (((cos(kk*rr - omega*time) - cos(omega*time))/(rr**3)) + (kk*sin(kk*rr - omega*time)/(rr**2)))
      end if
    end do

    ! Compute the transverse and longitudinal fields numerically
    call helmholtz_decomposition_trans_field(namespace, poisson_solver, grid, total_field, trans_field_computed)
    call helmholtz_decomposition_long_field(namespace, poisson_solver, grid, total_field, long_field_computed)

    ! Compute and output the results
    call dderivatives_div(grid%der, trans_field_computed, div_rot, set_bc = .true.)
    call check_norms_and_output_fields(grid, space, namespace, total_field, trans_field_exact, long_field_exact, &
     trans_field_computed, long_field_computed, div_rot, "hertz")

    SAFE_DEALLOCATE_A(total_field)
    SAFE_DEALLOCATE_A(trans_field_exact)
    SAFE_DEALLOCATE_A(long_field_exact)
    SAFE_DEALLOCATE_A(trans_field_computed)
    SAFE_DEALLOCATE_A(long_field_computed)
    SAFE_DEALLOCATE_A(div_rot)

    POP_SUB(helmholtz_decomposition_hertzian_dipole_test)

  end subroutine helmholtz_decomposition_hertzian_dipole_test

  subroutine helmholtz_decomposition_gaussian_test(namespace, poisson_solver, grid, space)
    type(namespace_t), intent(in) :: namespace
    type(poisson_t),   intent(in) :: poisson_solver
    type(grid_t),      intent(in) :: grid
    type(space_t),     intent(in) :: space

    FLOAT, allocatable :: total_field(:,:), trans_field_exact(:,:), long_field_exact(:,:)
    FLOAT, allocatable :: trans_field_computed(:,:), long_field_computed(:,:), div_rot(:)
    FLOAT              :: alpha, beta, a_prefactor, b_prefactor, ra, rb
    FLOAT              :: aa(space%dim), a_origin(space%dim), bb(space%dim), b_origin(space%dim)
    integer            :: ip

    PUSH_SUB(helmholtz_decomposition_gaussian_test)
    ASSERT(space%dim == 3)

    ! Allocate and initialize the fields
    SAFE_ALLOCATE(total_field(1:grid%np_part, 1:grid%box%dim))
    SAFE_ALLOCATE(trans_field_exact(1:grid%np_part, 1:grid%box%dim))
    SAFE_ALLOCATE(long_field_exact(1:grid%np_part, 1:grid%box%dim))
    SAFE_ALLOCATE(trans_field_computed(1:grid%np_part, 1:grid%box%dim))
    SAFE_ALLOCATE(long_field_computed(1:grid%np_part, 1:grid%box%dim))
    SAFE_ALLOCATE(div_rot(1:grid%np_part))
    total_field = M_ZERO
    trans_field_computed = M_ZERO
    trans_field_exact = M_ZERO
    long_field_computed = M_ZERO
    long_field_exact = M_ZERO
    div_rot = M_ZERO
    ! Define the parameters that define the gaussian
    alpha = CNST(2)
    a_prefactor = CNST(1.5)
    a_origin = M_ONE
    beta = CNST(3)
    b_prefactor = CNST(1.7)
    b_origin = -M_ONE
    
    ! Define the field to be tested 
    do ip = 1, grid%np
      call mesh_r(grid, ip, ra, coords = aa, origin = a_origin)
      long_field_exact(ip, 1) = (2 * a_prefactor / (alpha**2)) * aa(1) * exp(-(ra/alpha)**2)
      long_field_exact(ip, 2) = (2 * a_prefactor / (alpha**2)) * aa(2) * exp(-(ra/alpha)**2)
      long_field_exact(ip, 3) = (2 * a_prefactor / (alpha**2)) * aa(3) * exp(-(ra/alpha)**2)

      call mesh_r(grid, ip, rb, coords = bb, origin = b_origin)
      trans_field_exact(ip, 1) = - (2 * b_prefactor / (beta**2)) * bb(2) * exp(-(rb/beta)**2)
      trans_field_exact(ip, 2) = (2 * b_prefactor / (beta**2)) * bb(1) * exp(-(rb/beta)**2)

      total_field(ip, :) = long_field_exact(ip, :) + trans_field_exact(ip, :)
    end do

    ! Compute the transverse and longitudinal fields numerically
    call helmholtz_decomposition_trans_field(namespace, poisson_solver, grid, total_field, trans_field_computed)
    call helmholtz_decomposition_long_field(namespace, poisson_solver, grid, total_field, long_field_computed)

    ! Compute and output the results
    call dderivatives_div(grid%der, trans_field_computed, div_rot, set_bc = .true.)
    call check_norms_and_output_fields(grid, space, namespace, total_field, trans_field_exact, long_field_exact, &
     trans_field_computed, long_field_computed, div_rot, "gauss")

    SAFE_DEALLOCATE_A(total_field)
    SAFE_DEALLOCATE_A(trans_field_exact)
    SAFE_DEALLOCATE_A(long_field_exact)
    SAFE_DEALLOCATE_A(trans_field_computed)
    SAFE_DEALLOCATE_A(long_field_computed)
    SAFE_DEALLOCATE_A(div_rot)

    POP_SUB(helmholtz_decomposition_gaussian_test)

  end subroutine helmholtz_decomposition_gaussian_test

  subroutine check_norms_and_output_fields(grid, space, namespace, total_field, trans_field_exact, long_field_exact, &
   trans_field_computed, long_field_computed, div_rot, test)
    type(grid_t),      intent(in) :: grid
    type(space_t),     intent(in) :: space
    type(namespace_t), intent(in) :: namespace
    FLOAT,             intent(in) :: total_field(:,:), trans_field_exact(:,:), long_field_exact(:,:)
    FLOAT,             intent(in) :: trans_field_computed(:,:), long_field_computed(:,:), div_rot(:)
    character(len=*),  intent(in) :: test

    FLOAT              :: delta
    integer            :: ierr, iunit

    PUSH_SUB(check_norms_and_output_fields)

    iunit = io_open(test // "_helmholtz_decomposition", namespace, action='write')

    ! Check that the exact fields ere defined correctly
    delta = dmf_nrm2(grid, 3, total_field - trans_field_exact - long_field_exact)
    write(iunit, '(a,f19.13)') 'Norm of the difference between the exact fields = ', delta

    ! Check that the computed fields are equal to the exact ones (in norm)
    delta = dmf_nrm2(grid, 3, trans_field_computed - trans_field_exact)
    write(iunit, '(a,f19.13)') 'Helmholtz decomposition transverse field test (abs.) = ', delta
    delta = delta / dmf_nrm2(grid, 3, trans_field_exact)
    write(iunit, '(a,f19.13)') 'Helmholtz decomposition transverse field test (rel.) = ', delta
    delta = dmf_nrm2(grid, 3, long_field_computed - long_field_exact)
    write(iunit, '(a,f19.13)') 'Helmholtz decomposition longitudinal field test (abs.) = ', delta
    delta = delta / dmf_nrm2(grid, 3, long_field_exact)
    write(iunit, '(a,f19.13)') 'Helmholtz decomposition longitudinal field test (rel.) = ', delta
    delta = dmf_nrm2(grid, 3, total_field - trans_field_computed - long_field_computed)
    write(iunit, '(a,f19.13)') 'Self consistency of computed fields (transverse + longitudinal = total_field) (abs.) = ',delta
    delta = delta / dmf_nrm2(grid, 3, total_field)
    write(iunit, '(a,f19.13)') 'Self consistency of computed fields (transverse + longitudinal = total_field) (rel.) = ',delta

    call io_close(iunit)
    
    ! Output the values of the fields on the grid
    call io_function_output_vector (io_function_fill_how('AxisX'), ".", test // "_total_field", namespace, space, &
      grid, total_field, unit_one, ierr)
    call io_function_output_vector (io_function_fill_how('AxisX'), ".", test // "_trans_field_exact", namespace, space, &
      grid, trans_field_exact, unit_one, ierr)
    call io_function_output_vector (io_function_fill_how('AxisX'), ".", test // "_trans_field_comp", namespace, space, &
      grid, trans_field_computed, unit_one, ierr)
    call io_function_output_vector (io_function_fill_how('AxisX'), ".", test // "_long_field_exact", namespace, space, &
      grid, long_field_exact, unit_one, ierr)
    call io_function_output_vector (io_function_fill_how('AxisX'), ".", test // "_long_field_comp", namespace, space, &
      grid, long_field_computed, unit_one, ierr)

    POP_SUB(check_norms_and_output_fields)
  end subroutine check_norms_and_output_fields

#include "undef.F90"
#include "complex.F90"
#include "helmholtz_decomposition_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "helmholtz_decomposition_inc.F90"

end module helmholtz_decomposition_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
