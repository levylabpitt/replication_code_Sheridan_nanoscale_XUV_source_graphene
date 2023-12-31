!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

module fourier_shell_oct_m
  use cube_oct_m
  use debug_oct_m
  use fft_oct_m
  use global_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use space_oct_m
  use sort_oct_m

  implicit none

  private
  public ::                        &
    fourier_shell_t,               &
    fourier_shell_cutoff,          &
    fourier_shell_init,            &
    fourier_shell_end

  type fourier_shell_t
    ! Components are public by default
    integer              :: ngvectors
    FLOAT                :: ekin_cutoff
    integer, allocatable :: coords(:, :)
    integer, allocatable :: red_gvec(:, :)
  end type fourier_shell_t

contains

  FLOAT function fourier_shell_cutoff(space, cube, mesh, is_wfn, dg)
    type(space_t),   intent(in)  :: space
    type(cube_t),    intent(in)  :: cube
    class(mesh_t),   intent(in)  :: mesh
    logical,         intent(in)  :: is_wfn
    FLOAT, optional, intent(out) :: dg(:) !< (3)

    FLOAT :: dg_(3)

    PUSH_SUB(fourier_shell_cutoff)

    ! FIXME: what about anisotropic spacing?
    dg_(1:3) = M_PI/(cube%rs_n_global(1:3)/2*mesh%spacing(1:3))
    if (present(dg)) dg(1:3) = dg_(1:3)
    if (is_wfn .and. space%is_periodic()) then
      fourier_shell_cutoff = (dg_(1)*(cube%rs_n_global(1)/2-2))**2/M_TWO
    else
      fourier_shell_cutoff = (dg_(1)*(cube%rs_n_global(1)/2))**2/M_TWO
    end if

    POP_SUB(fourier_shell_cutoff)
  end function fourier_shell_cutoff

  subroutine fourier_shell_init(this, namespace, space, cube, mesh, kk)
    type(fourier_shell_t), intent(inout) :: this
    type(namespace_t),     intent(in)    :: namespace
    type(space_t),         intent(in)    :: space
    type(cube_t),          intent(in)    :: cube
    class(mesh_t),         intent(in)    :: mesh
    FLOAT, optional,       intent(in)    :: kk(:) !< (3)

    integer :: ig, ix, iy, iz, ixx(1:3), imap
    FLOAT :: dg(1:3), gvec(1:3)
    FLOAT, allocatable :: modg2(:)
    integer, allocatable :: map(:), ucoords(:, :), ured_gvec(:, :)
    integer(i8) :: number_points

    PUSH_SUB(fourier_shell_init)

    this%ekin_cutoff = fourier_shell_cutoff(space, cube, mesh, present(kk), dg = dg)

    ! make sure we do not run into integer overflow here
    number_points = cube%rs_n_global(1) * cube%rs_n_global(2)
    number_points = number_points * cube%rs_n_global(3)
    if (number_points >= HUGE(0)) then
      message(1) = "Error: too many points for the normal cube. Please try to use a distributed FFT."
      call messages_fatal(1, namespace=namespace)
    end if
    SAFE_ALLOCATE(modg2(1:product(cube%rs_n_global(1:3))))
    SAFE_ALLOCATE(ucoords(1:3, 1:product(cube%rs_n_global(1:3))))
    SAFE_ALLOCATE(ured_gvec(1:3, 1:product(cube%rs_n_global(1:3))))

    ig = 0
    ! According to the conventions of plane-wave codes, e.g. Quantum ESPRESSO,
    ! PARATEC, EPM, and BerkeleyGW, if the FFT grid is even, then neither
    ! nfft/2 nor -nfft/2 should be a valid G-vector component.
    do ix = 1, cube%rs_n_global(1)
      ixx(1) = pad_feq(ix, cube%rs_n_global(1), .true.)
      if (2 * ixx(1) == cube%rs_n_global(1)) cycle
      do iy = 1, cube%rs_n_global(2)
        ixx(2) = pad_feq(iy, cube%rs_n_global(2), .true.)
        if (2 * ixx(2) == cube%rs_n_global(2)) cycle
        do iz = 1, cube%rs_n_global(3)
          ixx(3) = pad_feq(iz, cube%rs_n_global(3), .true.)
          if (2 * ixx(3) == cube%rs_n_global(3)) cycle

          if (present(kk)) then
            gvec(1:3) = dg(1:3)*(ixx(1:3) + kk(1:3))
          else
            gvec(1:3) = dg(1:3)*ixx(1:3)
          end if

          if (sum(gvec(1:3)**2)/M_TWO <= this%ekin_cutoff + CNST(1e-10)) then
            ig = ig + 1
            ucoords(1:3, ig) = (/ ix, iy, iz /)
            ured_gvec(1:3, ig) = ixx(1:3)
            modg2(ig) = sum(gvec(1:3)**2)
          end if

        end do
      end do
    end do

    this%ngvectors = ig

    SAFE_ALLOCATE(this%coords(1:3, 1:this%ngvectors))
    SAFE_ALLOCATE(this%red_gvec(1:3, 1:this%ngvectors))
    SAFE_ALLOCATE(map(1:this%ngvectors))

    do ig = 1, this%ngvectors
      map(ig) = ig
    end do

    call sort(modg2(1:this%ngvectors), map)

    do ig = 1, this%ngvectors
      imap = map(ig)
      this%coords(1:3, ig) = ucoords(1:3, imap)
      this%red_gvec(1:3, ig) = ured_gvec(1:3, imap)
    end do

    SAFE_DEALLOCATE_A(ucoords)
    SAFE_DEALLOCATE_A(ured_gvec)
    SAFE_DEALLOCATE_A(modg2)
    SAFE_DEALLOCATE_A(map)

    POP_SUB(fourier_shell_init)
  end subroutine fourier_shell_init

  ! -----------------------------------------------------

  subroutine fourier_shell_end(this)
    type(fourier_shell_t), intent(inout) :: this

    PUSH_SUB(fourier_shell_end)

    SAFE_DEALLOCATE_A(this%coords)
    SAFE_DEALLOCATE_A(this%red_gvec)

    POP_SUB(fourier_shell_end)
  end subroutine fourier_shell_end

end module fourier_shell_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
