!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

! ------------------------------------------------------------------
! Preprocessor directives
! ------------------------------------------------------------------

#if TYPE == 2
#  define TYPE1 real(r8)
#  define TYPE2 real(r8)
#elif TYPE == 4
#  define TYPE1 complex(r8)
#  define TYPE2 real(r8)
#endif

#ifdef __GFORTRAN__
#define PASTE(x) x/**/_
#define FNAME(x) PASTE(x)TYPE
#else
#define FNAME(x) xFNAME(x, TYPE)
#define xFNAME(x,y) yFNAME(x,y)
#define yFNAME(x,y) x ## _ ## y
#endif


!> ------------------------------------------------------------------
!! BLAS level I
!! ------------------------------------------------------------------

!! ------------------------------------------------------------------
!! swap two vectors
!! ------------------------------------------------------------------

subroutine FNAME(swap_1)(n1, dx, dy)
  integer, intent(in)    :: n1
  TYPE1,   intent(inout) :: dx(:), dy(:)

  if (n1 < 1) return

  PUSH_SUB(FNAME(swap_1))

  ASSERT(ubound(dx, dim = 1) >= n1)
  ASSERT(ubound(dy, dim = 1) >= n1)
  ASSERT(not_in_openmp())

  call blas_swap(n1, dx(1), 1, dy(1), 1)

  POP_SUB(FNAME(swap_1))
end subroutine FNAME(swap_1)

subroutine FNAME(swap_2)(n1, n2, dx, dy)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(inout) :: dx(:,:), dy(:,:)

  integer :: ii

  if (n1*n2 < 1) return

  PUSH_SUB(FNAME(swap_2))

  ASSERT(ubound(dy, dim = 1) == ubound(dx, dim = 1))
  ASSERT(ubound(dx, dim = 2) >= n2)
  ASSERT(ubound(dy, dim = 2) >= n2)
  ASSERT(not_in_openmp())

  if (ubound(dx, dim = 1)  == n1) then
    call blas_swap(n1*n2, dx(1,1), 1, dy(1,1), 1)
  else
    do ii = 1, n2
      call blas_swap(n1, dx(1, ii), 1, dy(1, ii), 1)
    end do
  end if

  POP_SUB(FNAME(swap_2))
end subroutine FNAME(swap_2)

subroutine FNAME(swap_3)(n1, n2, n3, dx, dy)
  integer, intent(in)    :: n1, n2, n3
  TYPE1,   intent(inout) :: dx(:,:,:), dy(:,:,:)

  if (n1*n2*n3 < 1) return

  PUSH_SUB(FNAME(swap_3))

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dy, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) >= n3)
  ASSERT(ubound(dy, dim = 3) >= n3)
  ASSERT(not_in_openmp())

  call blas_swap(n1*n2*n3, dx(1,1,1), 1, dy(1,1,1), 1)

  POP_SUB(FNAME(swap_3))
end subroutine FNAME(swap_3)

subroutine FNAME(swap_4)(n1, n2, n3, n4, dx, dy)
  integer, intent(in)    :: n1, n2, n3, n4
  TYPE1,   intent(inout) :: dx(:,:,:,:), dy(:,:,:,:)

  if (n1*n2*n3*n4 < 1) return

  PUSH_SUB(FNAME(swap_4))

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dy, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) == n3)
  ASSERT(ubound(dy, dim = 3) == n3)
  ASSERT(ubound(dx, dim = 4) >= n4)
  ASSERT(ubound(dy, dim = 4) >= n4)
  ASSERT(not_in_openmp())

  call blas_swap(n1*n2*n3*n4, dx(1,1,1,1), 1, dy(1,1,1,1), 1)

  PUSH_SUB(FNAME(swap_4))
end subroutine FNAME(swap_4)

!> ------------------------------------------------------------------
!! scales a vector by a constant
!! ------------------------------------------------------------------

subroutine FNAME(scal_1)(n1, da, dx)
  integer, intent(in)    :: n1
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:)

  if (n1 < 1) return

  PUSH_SUB(FNAME(scal_1))

  ASSERT(ubound(dx, dim = 1) >= n1)
  ASSERT(not_in_openmp())

  call blas_scal(n1, da, dx(1), 1)

  POP_SUB(FNAME(scal_1))
end subroutine FNAME(scal_1)

subroutine FNAME(scal_2)(n1, n2, da, dx)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:, :)

  integer :: ii

  if (n1*n2 < 1) return

  PUSH_SUB(FNAME(scal_2))

  ASSERT(ubound(dx, dim = 2) >= n2)
  ASSERT(not_in_openmp())

  if (ubound(dx, dim = 1) == n1) then
    call blas_scal(n1*n2, da, dx(1,1), 1)
  else
    do ii = 1, n2
      call blas_scal(n1, da, dx(1, ii), 1)
    end do
  end if

  POP_SUB(FNAME(scal_2))
end subroutine FNAME(scal_2)

subroutine FNAME(scal_3)(n1, n2, n3, da, dx)
  integer, intent(in)    :: n1, n2, n3
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:, :, :)

  if (n1*n2*n3 < 1) return

  PUSH_SUB(FNAME(scal_3))

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) >= n3)
  ASSERT(not_in_openmp())

  call blas_scal(n1*n2*n3, da, dx(1,1,1), 1)

  POP_SUB(FNAME(scal_3))
end subroutine FNAME(scal_3)

subroutine FNAME(scal_4)(n1, n2, n3, n4, da, dx)
  integer, intent(in)    :: n1, n2, n3, n4
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:, :, :, :)

  if (n1*n2*n3*n4 < 1) return

  PUSH_SUB(FNAME(scal_4))

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) == n3)
  ASSERT(ubound(dx, dim = 4) >= n4)
  ASSERT(not_in_openmp())

  call blas_scal(n1*n2*n3*n4, da, dx(1,1,1,1), 1)

  POP_SUB(FNAME(scal_4))
end subroutine FNAME(scal_4)

#if TYPE == 4
subroutine FNAME(scal_5)(n1, da, dx)
  integer, intent(in)    :: n1
  TYPE2,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:)

  if (n1 < 1) return

  PUSH_SUB(FNAME(scal_5))

  ASSERT(ubound(dx, dim = 1) >= n1)
  ASSERT(not_in_openmp())

  call blas_scal(n1, da, dx(1))

  POP_SUB(FNAME(scal_5))
end subroutine FNAME(scal_5)

subroutine FNAME(scal_6)(n1, n2, da, dx)
  integer, intent(in)    :: n1, n2
  TYPE2,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:, :)

  integer :: ii

  if (n1*n2 < 1) return

  PUSH_SUB(FNAME(scal_6))

  ASSERT(ubound(dx, dim = 2) >= n2)
  ASSERT(not_in_openmp())

  if (ubound(dx, dim = 1) == n1) then
    call blas_scal(n1*n2, da, dx(1,1))
  else
    do ii = 1, n2
      call blas_scal(n1, da, dx(1, ii))
    end do
  end if

  POP_SUB(FNAME(scal_6))
end subroutine FNAME(scal_6)

#endif

!> ------------------------------------------------------------------
!! constant times a vector plus a vector
!! ------------------------------------------------------------------

subroutine FNAME(axpy_1)(n1, da, dx, dy)
  integer, intent(in)    :: n1
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:)
  TYPE1,   intent(inout) :: dy(:)

  if (n1 < 1) return

  PUSH_SUB(FNAME(axpy_1))

  ASSERT(ubound(dx, dim = 1) >= n1)
  ASSERT(ubound(dy, dim = 1) >= n1)
  ASSERT(not_in_openmp())

  call profiling_in(FNAME(axpy_profile), TOSTRING(FNAME(BLAS_AXPY)))

  call blas_axpy(n1, da, dx(1), 1, dy(1), 1)

#if TYPE == 2
  call profiling_count_operations(n1*2)
#else
  call profiling_count_operations(n1*8)
#endif

  call profiling_out(FNAME(axpy_profile))

  POP_SUB(FNAME(axpy_1))
end subroutine FNAME(axpy_1)

subroutine FNAME(axpy_2)(n1, n2, da, dx, dy)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:, :)
  TYPE1,   intent(inout) :: dy(:, :)

  integer :: ii

  if (n1*n2 < 1) return

  PUSH_SUB(FNAME(axpy_2))

  call profiling_in(FNAME(axpy_profile), TOSTRING(FNAME(BLAS_AXPY)))

  ASSERT(ubound(dx, dim = 1) >= n1)
  ASSERT(ubound(dy, dim = 1) >= n1)
  ASSERT(ubound(dx, dim = 2) >= n2)
  ASSERT(ubound(dy, dim = 2) >= n2)
  ASSERT(not_in_openmp())

  if (ubound(dx, dim = 1) == n1 .and. ubound(dy, dim = 1) == n1) then
    ASSERT(ubound(dy, dim = 1) == n1)
    call blas_axpy(n1*n2, da, dx(1,1), 1, dy(1,1), 1)
  else
    do ii = 1, n2
      call blas_axpy(n1, da, dx(1, ii), 1, dy(1, ii), 1)
    end do
  end if

#if TYPE == 2
  call profiling_count_operations(n1*n2*2)
#else
  call profiling_count_operations(n1*n2*8)
#endif

  call profiling_out(FNAME(axpy_profile))
  POP_SUB(FNAME(axpy_2))
end subroutine FNAME(axpy_2)

subroutine FNAME(axpy_3)(n1, n2, n3, da, dx, dy)
  integer, intent(in)    :: n1, n2, n3
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:, :, :)
  TYPE1,   intent(inout) :: dy(:, :, :)

  if (n1*n2*n3 < 1) return

  PUSH_SUB(FNAME(axpy_3))

  call profiling_in(FNAME(axpy_profile), TOSTRING(FNAME(BLAS_AXPY)))

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dy, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) >= n3)
  ASSERT(ubound(dy, dim = 3) >= n3)
  ASSERT(not_in_openmp())

  call blas_axpy(n1*n2*n3, da, dx(1,1,1), 1, dy(1,1,1), 1)

#if TYPE == 2
  call profiling_count_operations(n1*n2*n3*2)
#else
  call profiling_count_operations(n1*n2*n3*8)
#endif

  call profiling_out(FNAME(axpy_profile))
  POP_SUB(FNAME(axpy_3))
end subroutine FNAME(axpy_3)

subroutine FNAME(axpy_4)(n1, n2, n3, n4, da, dx, dy)
  integer, intent(in)    :: n1, n2, n3, n4
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:, :, :, :)
  TYPE1,   intent(inout) :: dy(:, :, :, :)

  if (n1*n2*n3*n4 < 1) return

  PUSH_SUB(FNAME(axpy_4))

  call profiling_in(FNAME(axpy_profile), TOSTRING(FNAME(BLAS_AXPY)))

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dy, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) == n3)
  ASSERT(ubound(dy, dim = 3) == n3)
  ASSERT(ubound(dx, dim = 4) >= n4)
  ASSERT(ubound(dy, dim = 4) >= n4)
  ASSERT(not_in_openmp())

  call blas_axpy(n1*n2*n3*n4, da, dx(1,1,1,1), 1, dy(1,1,1,1), 1)

#if TYPE == 2
  call profiling_count_operations(n1*n2*n3*n4*2)
#else
  call profiling_count_operations(n1*n2*n3*n4*8)
#endif

  call profiling_out(FNAME(axpy_profile))
  POP_SUB(FNAME(axpy_2))
end subroutine FNAME(axpy_4)

#if TYPE == 4
subroutine FNAME(axpy_5)(n1, da, dx, dy)
  integer, intent(in)    :: n1
  TYPE2,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:)
  TYPE1,   intent(inout) :: dy(:)

  if (n1 < 1) return

  PUSH_SUB(FNAME(axpy_5))

  ASSERT(ubound(dx, dim = 1) >= n1)
  ASSERT(ubound(dy, dim = 1) >= n1)
  ASSERT(not_in_openmp())

  call profiling_in(FNAME(axpy_profile), TOSTRING(FNAME(BLAS_AXPY)))

  call blas_axpy(n1, da, dx(1), dy(1))

  call profiling_count_operations(n1*4)

  call profiling_out(FNAME(axpy_profile))

  POP_SUB(FNAME(axpy_5))
end subroutine FNAME(axpy_5)

subroutine FNAME(axpy_6)(n1, n2, da, dx, dy)
  integer, intent(in)    :: n1
  integer, intent(in)    :: n2
  TYPE2,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:, :)
  TYPE1,   intent(inout) :: dy(:, :)

  integer :: ii

  if (n1 < 1) return

  PUSH_SUB(FNAME(axpy_6))

  call profiling_in(FNAME(axpy_profile), TOSTRING(FNAME(BLAS_AXPY)))

  ASSERT(ubound(dx, dim = 1) >= n1)
  ASSERT(ubound(dy, dim = 1) >= n1)
  ASSERT(ubound(dx, dim = 2) >= n2)
  ASSERT(ubound(dy, dim = 2) >= n2)
  ASSERT(not_in_openmp())

  if (ubound(dx, dim = 1) == n1 .and. ubound(dy, dim = 1) == n1) then
    call blas_axpy(n1*n2, da, dx(1, 1), dy(1, 1))
  else
    do ii = 1, n2
      call blas_axpy(n1, da, dx(1, ii), dy(1, ii))
    end do
  end if

  call profiling_count_operations(n1*n2*4)

  call profiling_out(FNAME(axpy_profile))

  POP_SUB(FNAME(axpy_6))
end subroutine FNAME(axpy_6)

subroutine FNAME(axpy_7)(n1, n2, n3, da, dx, dy)
  integer, intent(in)    :: n1, n2, n3
  TYPE2,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:, :, :)
  TYPE1,   intent(inout) :: dy(:, :, :)

  if (n1*n2*n3 < 1) return

  PUSH_SUB(FNAME(axpy_7))

  call profiling_in(FNAME(axpy_profile), TOSTRING(FNAME(BLAS_AXPY)))

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dy, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) >= n3)
  ASSERT(ubound(dy, dim = 3) >= n3)
  ASSERT(not_in_openmp())

  call blas_axpy(n1*n2*n3, da, dx(1,1,1), dy(1,1,1))

  call profiling_count_operations(n1*n2*n3*4)

  call profiling_out(FNAME(axpy_profile))
  POP_SUB(FNAME(axpy_7))
end subroutine FNAME(axpy_7)
#endif

!> ------------------------------------------------------------------
!! Copies a vector x, to a vector y
!! ------------------------------------------------------------------

subroutine FNAME(copy_1)(n1, dx, dy)
  integer, intent(in)    :: n1
  TYPE1,   intent(in)    :: dx(:)
  TYPE1,   intent(inout) :: dy(:)

  if (n1 < 1) return

  PUSH_SUB(FNAME(copy_1))

  ASSERT(ubound(dx, dim = 1) >= n1)
  ASSERT(ubound(dy, dim = 1) >= n1)
  ASSERT(not_in_openmp())

  call profiling_in(FNAME(copy_profile), TOSTRING(FNAME(BLAS_COPY)))

  call blas_copy(n1, dx(1), 1, dy(1), 1)

  call profiling_count_transfers(n1, dx(1))

  call profiling_out(FNAME(copy_profile))
  POP_SUB(FNAME(copy_1))
end subroutine FNAME(copy_1)

subroutine FNAME(copy_2)(n1, n2, dx, dy)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(in)    :: dx(:,:)
  TYPE1,   intent(inout) :: dy(:,:)

  integer :: ii

  if (n1*n2 < 1) return

  PUSH_SUB(FNAME(copy_2))

  call profiling_in(FNAME(copy_profile), TOSTRING(FNAME(BLAS_COPY)))

  ASSERT(ubound(dx, dim = 1) >= n1)
  ASSERT(ubound(dy, dim = 1) >= n1)
  ASSERT(ubound(dx, dim = 2) >= n2)
  ASSERT(ubound(dy, dim = 2) >= n2)
  ASSERT(not_in_openmp())

  if (ubound(dx, dim = 1) == n1 .and. ubound(dy, dim = 1) == n1) then
    call blas_copy(n1*n2, dx(1,1), 1, dy(1,1), 1)
  else
    do ii = 1, n2
      call blas_copy(n1, dx(1, ii), 1, dy(1, ii), 1)
    end do
  end if

  call profiling_out(FNAME(copy_profile))
  POP_SUB(FNAME(copy_2))
end subroutine FNAME(copy_2)

subroutine FNAME(copy_3)(n1, n2, n3, dx, dy)
  integer, intent(in)    :: n1, n2, n3
  TYPE1,   intent(in)    :: dx(:,:,:)
  TYPE1,   intent(inout) :: dy(:,:,:)

  if (n1*n2*n3 < 1) return

  PUSH_SUB(FNAME(copy_3))

  call profiling_in(FNAME(copy_profile), TOSTRING(FNAME(BLAS_COPY)))

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dy, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) >= n3)
  ASSERT(ubound(dy, dim = 3) >= n3)
  ASSERT(not_in_openmp())

  call blas_copy (n1*n2*n3, dx(1,1,1), 1, dy(1,1,1), 1)

  call profiling_out(FNAME(copy_profile))
  POP_SUB(FNAME(copy_3))
end subroutine FNAME(copy_3)

subroutine FNAME(copy_4)(n1, n2, n3, n4, dx, dy)
  integer, intent(in)    :: n1, n2, n3, n4
  TYPE1,   intent(in)    :: dx(:,:,:,:)
  TYPE1,   intent(inout) :: dy(:,:,:,:)

  if (n1*n2*n3*n4 < 1) return

  PUSH_SUB(FNAME(copy_4))

  call profiling_in(FNAME(copy_profile), TOSTRING(FNAME(BLAS_COPY)))

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dy, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) == n3)
  ASSERT(ubound(dy, dim = 3) == n3)
  ASSERT(ubound(dx, dim = 4) >= n4)
  ASSERT(ubound(dy, dim = 4) >= n4)
  ASSERT(not_in_openmp())

  call blas_copy (n1*n2*n3*n4, dx(1,1,1,1), 1, dy(1,1,1,1), 1)

  call profiling_out(FNAME(copy_profile))
  POP_SUB(FNAME(copy_4))
end subroutine FNAME(copy_4)

!> ------------------------------------------------------------------
!! Returns the euclidean norm of a vector
!! ------------------------------------------------------------------

TYPE2 function FNAME(nrm2)(n, dx) result(nrm2)
  integer, intent(in) :: n
  TYPE1,   intent(in) :: dx(:)

  PUSH_SUB(FNAME(nrm2))

  nrm2 = CNST(0.0)
  if (n < 1) then
    POP_SUB(FNAME(nrm2))
    return
  end if

  ASSERT(ubound(dx, dim = 1) >= n)
  ASSERT(not_in_openmp())

  nrm2 = blas_nrm2(n, dx(1), 1)

  POP_SUB(FNAME(nrm2))
end function FNAME(nrm2)

!> ------------------------------------------------------------------
!! BLAS level II
!! ------------------------------------------------------------------

!> ------------------------------------------------------------------
!! Matrix-vector multiplication plus vector.
!! ------------------------------------------------------------------

subroutine FNAME(symv_1)(n, alpha, a, x, beta, y)
  integer, intent(in)    :: n
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:, :)
  TYPE1,   intent(in)    :: x(:)
  TYPE1,   intent(inout) :: y(:)

  ! no push_sub, called too frequently

  ASSERT(ubound(a, dim=1) >= n)
  ASSERT(not_in_openmp())

  call profiling_in(FNAME(symv_profile), TOSTRING(FNAME(BLAS_SYMV)))
  call blas_symv('U', n, alpha, a(1, 1), lead_dim(a), x(1), 1, beta, y(1), 1)
  call profiling_out(FNAME(symv_profile))

end subroutine FNAME(symv_1)

subroutine FNAME(symv_2)(n1, n2, alpha, a, x, beta, y)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:, :, :)
  TYPE1,   intent(in)    :: x(:)
  TYPE1,   intent(inout) :: y(:, :)

  PUSH_SUB(FNAME(symv_2))

  ASSERT(ubound(a, dim=1) == n1)
  ASSERT(ubound(a, dim=2) == n2)
  ASSERT(ubound(y, dim=1) == n1)
  ASSERT(ubound(y, dim=2) >= n2)
  ASSERT(not_in_openmp())

  call profiling_in(FNAME(symv_profile), TOSTRING(FNAME(BLAS_SYMV)))
  call blas_symv('U', n1*n2, alpha, a(1, 1, 1), n1*n2, x(1), 1, beta, y(1, 1), 1)
  call profiling_out(FNAME(symv_profile))

  POP_SUB(FNAME(symv_2))
end subroutine FNAME(symv_2)

subroutine FNAME(gemv_1)(m, n, alpha, a, x, beta, y)
  integer, intent(in)    :: m, n
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:,:)
  TYPE1,   intent(in)    :: x(:)
  TYPE1,   intent(inout) :: y(:)

  PUSH_SUB(FNAME(gemv_1))

  ASSERT(ubound(a, dim=1) >= m)
  ASSERT(not_in_openmp())

  call profiling_in(FNAME(gemv_profile), TOSTRING(FNAME(BLAS_GEMV)))
  call blas_gemv('N', m, n, alpha, a(1,1), lead_dim(a), x(1), 1, beta, y(1), 1)
  call profiling_out(FNAME(gemv_profile))

  POP_SUB(FNAME(gemv_1))
end subroutine FNAME(gemv_1)

subroutine FNAME(gemv_2)(m1, m2, n, alpha, a, x, beta, y)
  integer, intent(in)    :: m1, m2, n
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:,:,:)
  TYPE1,   intent(in)    :: x(:)
  TYPE1,   intent(inout) :: y(:,:)

  PUSH_SUB(FNAME(gemv_2))

  ASSERT(ubound(a, dim=1) == m1)
  ASSERT(ubound(a, dim=2) == m2)
  ASSERT(ubound(y, dim=1) == m1)
  ASSERT(ubound(y, dim=2) >= m2)
  ASSERT(not_in_openmp())

  call profiling_in(FNAME(gemv_profile), TOSTRING(FNAME(BLAS_GEMV)))
  call blas_gemv('N', m1*m2, n, alpha, a(1,1,1), m1*m2, x(1), 1, beta, y(1,1), 1)
  call profiling_out(FNAME(gemv_profile))

  POP_SUB(FNAME(gemv_2))
end subroutine FNAME(gemv_2)


!> ------------------------------------------------------------------
!! BLAS level III
!! ------------------------------------------------------------------

!> ------------------------------------------------------------------
!! Matrix-matrix multiplication plus matrix.
!! ------------------------------------------------------------------

subroutine FNAME(gemm_1)(m, n, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m, n, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:,:)  !< a(m, k)
  TYPE1,   intent(in)    :: b(:,:)  !< b(k, n)
  TYPE1,   intent(inout) :: c(:,:)  !< c(m, n)

  ! no PUSH SUB, called too often

  ASSERT(ubound(a, dim=1) >= m)
  ASSERT(ubound(a, dim=2) >= k)
  ASSERT(ubound(b, dim=1) >= k)
  ASSERT(ubound(b, dim=2) >= n)
  ASSERT(ubound(c, dim=1) >= m)
  ASSERT(ubound(c, dim=2) >= n)
  ASSERT(not_in_openmp())

  call blas_gemm('N', 'N', m, n, k, alpha, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), beta, c(1, 1), lead_dim(c))

end subroutine FNAME(gemm_1)

subroutine FNAME(gemm_2)(m1, m2, n, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m1, m2, n, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:, :, :)  !< a(m1, m2, k)
  TYPE1,   intent(in)    :: b(:, :)     !< b(k, n)
  TYPE1,   intent(inout) :: c(:, :, :)  !< c(m1, m2, n)

  PUSH_SUB(FNAME(gemm_2))

  ASSERT(ubound(a, dim=1) == m1)
  ASSERT(ubound(a, dim=2) == m2)
  ASSERT(ubound(a, dim=3) >= k)
  ASSERT(ubound(b, dim=1) >= k)
  ASSERT(ubound(c, dim=1) == m1)
  ASSERT(ubound(c, dim=2) == m2)
  ASSERT(ubound(c, dim=3) >= n)
  ASSERT(not_in_openmp())

  call blas_gemm('N', 'N', m1*m2, n, k, alpha, a(1, 1, 1), m1*m2, &
    b(1, 1), lead_dim(b), beta, c(1, 1, 1), m1*m2)

  POP_SUB(FNAME(gemm_2))
end subroutine FNAME(gemm_2)

!> The same as above but with (Hermitian) transpose of a.
subroutine FNAME(gemm_cn_1)(m, n, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m, n, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:,:)  !< a(k, m)
  TYPE1,   intent(in)    :: b(:,:)  !< b(k, n)
  TYPE1,   intent(inout) :: c(:,:)  !< c(m, n)

  ! no PUSH_SUB, called too often

  ASSERT(ubound(a, dim=1) >= k)
  ASSERT(ubound(a, dim=2) >= m)
  ASSERT(ubound(b, dim=1) >= k)
  ASSERT(ubound(b, dim=2) >= n)
  ASSERT(ubound(c, dim=1) >= m)
  ASSERT(ubound(c, dim=2) >= n)
  ASSERT(not_in_openmp())

  call blas_gemm('C', 'N', m, n, k, alpha, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), beta, c(1, 1), lead_dim(c))

end subroutine FNAME(gemm_cn_1)

subroutine FNAME(gemm_cn_2)(m1, m2, n1, n2, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m1, m2, n1, n2, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:, :, :)  !< a(k, m2, m1)
  TYPE1,   intent(in)    :: b(:, :, :)  !< b(k, n2, n1)
  TYPE1,   intent(inout) :: c(:, :)     !< c(m1*m2, n1*n2)

  PUSH_SUB(FNAME(gemm_cn_2))

  ASSERT(ubound(a, dim=1) >= k)
  ASSERT(ubound(a, dim=2) == m2)
  ASSERT(ubound(a, dim=3) == m1)
  ASSERT(ubound(b, dim=1) >= k)
  ASSERT(ubound(b, dim=2) == n2)
  ASSERT(ubound(b, dim=3) >= n1)
  ASSERT(ubound(c, dim=1) >= m1*m2)
  ASSERT(ubound(c, dim=2) >= n1*n2)
  ASSERT(not_in_openmp())

  call blas_gemm('C', 'N', m1*m2, n1*n2, k, alpha, a(1, 1, 1), lead_dim(a), &
    b(1, 1, 1), lead_dim(b), beta, c(1, 1), lead_dim(c))

  POP_SUB(FNAME(gemm_cn_2))
end subroutine FNAME(gemm_cn_2)

!> The same as gemm but with (Hermitian) transpose of b.
subroutine FNAME(gemm_nc_1)(m, n, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m, n, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:,:)  !< a(m, k)
  TYPE1,   intent(in)    :: b(:,:)  !< b(n, k)
  TYPE1,   intent(inout) :: c(:,:)  !< c(m, n)

  ! no PUSH_SUB, called too often

  ASSERT(ubound(a, dim=1) >= m)
  ASSERT(ubound(a, dim=2) >= k)
  ASSERT(ubound(b, dim=1) >= k)
  ASSERT(ubound(b, dim=2) >= n)
  ASSERT(ubound(c, dim=1) >= m)
  ASSERT(ubound(c, dim=2) >= n)
  ASSERT(not_in_openmp())

  call blas_gemm('N', 'C', m, n, k, alpha, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), beta, c(1, 1), lead_dim(c))

end subroutine FNAME(gemm_nc_1)

subroutine FNAME(gemm_nc_2)(m1, m2, n1, n2, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m1, m2, n1, n2, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:, :, :)  !< a(k, m2, m1)
  TYPE1,   intent(in)    :: b(:, :, :)  !< b(k, n2, n1)
  TYPE1,   intent(inout) :: c(:, :)     !< c(m1*m2, n1*n2)

  PUSH_SUB(FNAME(gemm_nc_2))

  ASSERT(ubound(a, dim=1) == m1)
  ASSERT(ubound(a, dim=2) == m2)
  ASSERT(ubound(a, dim=3) >= k)
  ASSERT(ubound(b, dim=1) == n1)
  ASSERT(ubound(b, dim=2) == n2)
  ASSERT(ubound(b, dim=3) >= k)
  ASSERT(ubound(c, dim=1) >= m1*m2)
  ASSERT(ubound(c, dim=2) >= n1*n2)
  ASSERT(not_in_openmp())

  call blas_gemm('N', 'C', m1*m2, n1*n2, k, alpha, a(1, 1, 1), lead_dim(a), &
    b(1, 1, 1), lead_dim(b), beta, c(1, 1), lead_dim(c))

  POP_SUB(FNAME(gemm_nc_2))
end subroutine FNAME(gemm_nc_2)

!> The following matrix multiplications all expect upper triangular matrices for a.
!! For real matrices, a = a^T, for complex matrices a = a^H.
subroutine FNAME(symm_1)(m, n, side, alpha, a, b, beta, c)
  integer,      intent(in)    :: m, n
  character(1), intent(in)    :: side
  TYPE1,        intent(in)    :: alpha, beta, a(:, :), b(:, :)
  TYPE1,        intent(inout) :: c(:, :) !c(m, n)

  ! no push_sub, called too frequently
  !The size specified are for the matrix C
  ASSERT(ubound(c, dim=1) >= m)
  ASSERT(ubound(c, dim=2) >= n)
  ASSERT(not_in_openmp())

  select case (side)
  case ('l', 'L') ! Here we compute C := alpha*A*B + beta*C
    ASSERT(ubound(a, dim=1) >= m)
    ASSERT(ubound(b, dim=1) >= n)
  case ('r', 'R') ! Here we compute C := alpha*B*A + beta*C
    ASSERT(ubound(a, dim=1) >= n)
    ASSERT(ubound(b, dim=1) >= m)
  end select

  call blas_symm(side, 'U', m, n, alpha, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), beta, c(1, 1), lead_dim(c))

end subroutine FNAME(symm_1)

!> ------------------------------------------------------------------
!! Matrix-matrix multiplication.
!! ------------------------------------------------------------------

subroutine FNAME(trmm_1)(m, n, uplo, transa, side, alpha, a, b)
  integer,      intent(in)    :: m, n
  character(1), intent(in)    :: side, transa, uplo
  TYPE1,        intent(in)    :: alpha
  TYPE1,        intent(in)    :: a(:, :) !< a(m, m), upper triangular matrix.
  TYPE1,        intent(inout) :: b(:, :) !< b(m, n).

  ! no push_sub, called too frequently

  ASSERT(ubound(b, dim=1) >= m)
  ASSERT(ubound(b, dim=2) >= n)
  ASSERT(not_in_openmp())

  select case (side)
  case ('L', 'l')
    ASSERT(ubound(a, dim=1) >= m)
    ASSERT(ubound(a, dim=2) >= m)
  case ('R', 'r')
    ASSERT(ubound(a, dim=1) >= n)
    ASSERT(ubound(a, dim=2) >= n)
  end select

  call blas_trmm(side, uplo, transa, 'N', m, n, alpha, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b))

end subroutine FNAME(trmm_1)


! ------------------------------------------------------------------
! Clean up preprocessor directives
! ------------------------------------------------------------------

#undef ARG_LIST
#undef ARG_CALL
#undef TYPE1
#undef TYPE2
#undef FNAME
#undef xFNAME
#undef yFNAME


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
