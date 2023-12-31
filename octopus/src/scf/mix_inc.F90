!! Copyright (C) 2002-2015 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
!! Copyright (C) 2021 N. Tancogne-Dejean
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
subroutine X(mixing)(namespace, smix, vin, vout, vnew)
  type(namespace_t), intent(in)    :: namespace
  type(mix_t),       intent(inout) :: smix
  R_TYPE,            intent(in)    :: vin(:, :, :), vout(:, :, :)
  R_TYPE,            intent(out)   :: vnew(:, :, :)

  integer :: ii

  PUSH_SUB(X(mixing))

  ASSERT(associated(smix%der))

  smix%iter = smix%iter + 1

  select case (smix%scheme)
  case (OPTION__MIXINGSCHEME__LINEAR)
    call X(mixing_linear)(smix%coeff, smix%mixfield%d1, smix%mixfield%d2, smix%mixfield%d3, vin, vout, vnew)
    do ii = 1, smix%nauxmixfield
      if (smix%auxmixfield(ii)%p%func_type == TYPE_FLOAT) then
        call dmixing_linear(smix%coeff, smix%auxmixfield(ii)%p%d1, smix%auxmixfield(ii)%p%d2, smix%auxmixfield(ii)%p%d3, &
          smix%auxmixfield(ii)%p%dvin, smix%auxmixfield(ii)%p%dvout, smix%auxmixfield(ii)%p%dvnew)
      else
        call zmixing_linear(smix%coeff, smix%auxmixfield(ii)%p%d1, smix%auxmixfield(ii)%p%d2, smix%auxmixfield(ii)%p%d3, &
          smix%auxmixfield(ii)%p%zvin, smix%auxmixfield(ii)%p%zvout, smix%auxmixfield(ii)%p%zvnew)
      end if
    end do

  case (OPTION__MIXINGSCHEME__BROYDEN)
    call X(mixing_broyden)(namespace, smix, vin, vout, vnew, smix%iter)

  case (OPTION__MIXINGSCHEME__DIIS)
    call X(mixing_diis)(smix, vin, vout, vnew, smix%iter)

  case (OPTION__MIXINGSCHEME__BOWLER_GILLAN)
    call X(mixing_grpulay)(smix, vin, vout, vnew, smix%iter)

  end select

  POP_SUB(X(mixing))
end subroutine X(mixing)


! ---------------------------------------------------------
subroutine X(mixing_linear)(coeff, d1, d2, d3, vin, vout, vnew)
  FLOAT,        intent(in) :: coeff
  integer,      intent(in) :: d1, d2, d3
  R_TYPE,       intent(in) :: vin(:, :, :), vout(:, :, :)
  R_TYPE,       intent(out):: vnew(:, :, :)

  PUSH_SUB(X(mixing_linear))

  vnew(1:d1, 1:d2, 1:d3) = vin(1:d1, 1:d2, 1:d3)*(M_ONE - coeff) + coeff*vout(1:d1, 1:d2, 1:d3)

  POP_SUB(X(mixing_linear))
end subroutine X(mixing_linear)


! ---------------------------------------------------------
subroutine X(mixing_broyden)(namespace, smix, vin, vout, vnew, iter)
  type(namespace_t), intent(in)    :: namespace
  type(mix_t),       intent(inout) :: smix
  R_TYPE,            intent(in)    :: vin(:, :, :), vout(:, :, :)
  R_TYPE,            intent(out)   :: vnew(:, :, :)
  integer,           intent(in)    :: iter

  integer :: ipos, iter_used, d1, d2, d3
  R_TYPE, allocatable :: f(:, :, :)

  PUSH_SUB(X(mixing_broyden))

  d1 = smix%mixfield%d1
  d2 = smix%mixfield%d2
  d3 = smix%mixfield%d3

  SAFE_ALLOCATE(f(1:d1, 1:d2, 1:d3))

  f(1:d1, 1:d2, 1:d3) = vout(1:d1, 1:d2, 1:d3) - vin(1:d1, 1:d2, 1:d3)
  ! According to Johnson (1988), df and dv should be normalized here.
  ! However, it turns out that this blocks convergence to very high
  ! accuracies. Also, the normalization is not present in the description
  ! by Kresse & Furthmueller (1996). Thus, we do not normalize here.
  if (iter > 1) then
    ! Store df and dv from current iteration
    ipos = mod(smix%last_ipos, smix%ns) + 1
    smix%ipos = ipos

    call lalg_copy(d1, d2, d3, f(:, :, :), smix%mixfield%X(df)(:, :, :, ipos))
    call lalg_copy(d1, d2, d3, vin(:, :, :), smix%mixfield%X(dv)(:, :, :, ipos))
    call lalg_axpy(d1, d2, d3, -M_ONE, smix%mixfield%X(f_old)(:, :, :),   smix%mixfield%X(df)(:, :, :, ipos))
    call lalg_axpy(d1, d2, d3, -M_ONE, smix%mixfield%X(vin_old)(:, :, :), smix%mixfield%X(dv)(:, :, :, ipos))

    smix%last_ipos = ipos
  end if

  ! Store residual and vin for next iteration
  smix%mixfield%X(vin_old)(1:d1, 1:d2, 1:d3) = vin(1:d1, 1:d2, 1:d3)
  smix%mixfield%X(f_old)  (1:d1, 1:d2, 1:d3) = f  (1:d1, 1:d2, 1:d3)

  ! extrapolate new vector
  iter_used = min(iter - 1, smix%ns)
  call X(broyden_extrapolation)(smix, namespace, smix%coeff, d1, d2, d3, &
    vin, vnew, iter_used, f, smix%mixfield%X(df), smix%mixfield%X(dv))

  SAFE_DEALLOCATE_A(f)

  POP_SUB(X(mixing_broyden))
end subroutine X(mixing_broyden)


! ---------------------------------------------------------
subroutine X(broyden_extrapolation)(this, namespace, coeff, d1, d2, d3, vin, vnew, iter_used, f, df, dv)
  type(mix_t),       intent(inout) :: this
  type(namespace_t), intent(in)    :: namespace
  FLOAT,             intent(in)    :: coeff
  integer,           intent(in)    :: d1, d2, d3, iter_used
  R_TYPE,            intent(in)    :: vin(:, :, :), f(:, :, :)
  R_TYPE,            intent(in)    :: df(:, :, :, :), dv(:, :, :, :)
  R_TYPE,            intent(out)   :: vnew(:, :, :)

  ! The parameter ww does not influence the mixing (see description
  ! in Kresse & Furthmueller 1996). Thus we choose it to be one.
  FLOAT, parameter :: ww = M_ONE
  integer  :: i, k, l
  R_TYPE    :: gamma
  R_TYPE, allocatable :: beta(:, :), work(:)

  PUSH_SUB(X(broyden_extrapolation))

  if (iter_used == 0) then
    ! linear mixing...
    vnew(1:d1, 1:d2, 1:d3) = vin(1:d1, 1:d2, 1:d3) + coeff*f(1:d1, 1:d2, 1:d3)

    do i = 1, this%nauxmixfield
      if (this%auxmixfield(i)%p%func_type == TYPE_FLOAT) then
        call dbroyden_extrapolation_aux(this, namespace, i, coeff, iter_used)
      else
        call zbroyden_extrapolation_aux(this, namespace, i, coeff, iter_used)
      end if
    end do

    POP_SUB(X(broyden_extrapolation))
    return
  end if

  SAFE_ALLOCATE(beta(1:iter_used, 1:iter_used))
  SAFE_ALLOCATE(work(1:iter_used))

  ! compute matrix beta, Johnson eq. 13a
  call X(mixing_build_matrix)(this, df, iter_used, d2, d3, ww*ww, beta)

  ! invert matrix beta
  call lalg_inverter(iter_used, beta)

  do i = 1, iter_used
    work(i) = M_ZERO
    do l = 1, d3
      do k = 1, d2
        work(i) = work(i) + X(mix_dotp)(this, df(:, k, l, i), f(:, k, l), reduce = .false.)
      end do
    end do
    if (this%der%mesh%parallel_in_domains) call this%der%mesh%allreduce(work(i))
  end do

  ! linear mixing term
  vnew(1:d1, 1:d2, 1:d3) = vin(1:d1, 1:d2, 1:d3) + coeff*f(1:d1, 1:d2, 1:d3)

  ! other terms
  do i = 1, iter_used
    gamma = ww*sum(beta(:, i)*work(:))
    vnew(1:d1, 1:d2, 1:d3) = vnew(1:d1, 1:d2, 1:d3) - ww*gamma*(coeff*df(1:d1, 1:d2, 1:d3, i) + dv(1:d1, 1:d2, 1:d3, i))
  end do

  do i = 1, this%nauxmixfield
    if (this%auxmixfield(i)%p%func_type == TYPE_FLOAT) then
      call dbroyden_extrapolation_aux(this, namespace, i, coeff, iter_used, X(beta)=beta, X(work)=work)
    else
      call zbroyden_extrapolation_aux(this, namespace, i, coeff, iter_used, X(beta)=beta, X(work)=work)
    end if
  end do

  SAFE_DEALLOCATE_A(beta)
  SAFE_DEALLOCATE_A(work)

  POP_SUB(X(broyden_extrapolation))
end subroutine X(broyden_extrapolation)

!--------------------------------------------------------------------

subroutine X(broyden_extrapolation_aux)(this, namespace, ii, coeff, iter_used, dbeta, dwork, zbeta, zwork)
  type(mix_t),     intent(inout) :: this
  type(namespace_t),  intent(in) :: namespace
  integer,            intent(in) :: ii
  FLOAT,              intent(in) :: coeff
  integer,            intent(in) :: iter_used
  FLOAT, optional,    intent(in) :: dbeta(:, :), dwork(:)
  CMPLX, optional,    intent(in) :: zbeta(:, :), zwork(:)

  FLOAT, parameter :: ww = M_ONE
  integer :: d1,d2,d3, ipos, i
  R_TYPE  :: gamma
  R_TYPE, allocatable :: f(:, :, :)
  type(mixfield_t), pointer :: mf

  PUSH_SUB(X(broyden_extrapolation_aux))

#ifdef R_TREAL
  ! We cannot have a complex auxiliary field being mixed along with a real field
  ASSERT(this%auxmixfield(ii)%p%func_type == TYPE_FLOAT)
#endif

  mf => this%auxmixfield(ii)%p

  d1 = mf%d1
  d2 = mf%d2
  d3 = mf%d3

  SAFE_ALLOCATE(f(1:d1, 1:d2, 1:d3))

  f(1:d1, 1:d2, 1:d3) = mf%X(vout)(1:d1, 1:d2, 1:d3) - mf%X(vin)(1:d1, 1:d2, 1:d3)

  if (this%iter > 1) then
    ! Store df and dv from current iteration
    ipos = this%ipos

    call lalg_copy(d1, d2, d3, f(:, :, :), mf%X(df)(:, :, :, ipos))
    call lalg_copy(d1, d2, d3, mf%X(vin)(:, :, :), mf%X(dv)(:, :, :, ipos))
    call lalg_axpy(d1, d2, d3, -M_ONE, mf%X(f_old)(:, :, :), mf%X(df)(:, :, :, ipos))
    call lalg_axpy(d1, d2, d3, -M_ONE, mf%X(vin_old)(:, :, :), mf%X(dv)(:, :, :, ipos))
  end if

  ! Store residual and vin for next iteration
  mf%X(vin_old)(1:d1, 1:d2, 1:d3) = mf%X(vin)(1:d1, 1:d2, 1:d3)
  mf%X(f_old)  (1:d1, 1:d2, 1:d3) = f  (1:d1, 1:d2, 1:d3)

  ! linear mixing term
  mf%X(vnew)(1:d1, 1:d2, 1:d3) = mf%X(vin)(1:d1, 1:d2, 1:d3) + coeff*f(1:d1, 1:d2, 1:d3)

  if (iter_used == 0) then
    !We stop here
    POP_SUB(X(broyden_extrapolation_aux))
    return
  end if


  ! other terms
  do i = 1, iter_used
    !Here gamma might be of a different type than the main mixfield type, so we convert it to the proper type
    if (present(dbeta) .and. present(dwork)) then
      gamma = ww*sum(dbeta(:, i)*dwork(:))
    else if (present(zbeta) .and. present(zwork)) then
#ifdef R_TCOMPLX
      gamma = ww*sum(zbeta(:, i)*zwork(:))
#else
      ASSERT(.false.)
#endif
    else
      write(message(1), '(a)') 'Internal error in broyden_extrapolation_aux'
      call messages_fatal(1, namespace=namespace)
    end if
    mf%X(vnew)(1:d1, 1:d2, 1:d3) = mf%X(vnew)(1:d1, 1:d2, 1:d3) &
      - ww*gamma*(coeff*mf%X(df)(1:d1, 1:d2, 1:d3, i) + mf%X(dv)(1:d1, 1:d2, 1:d3, i))
  end do

  SAFE_DEALLOCATE_A(f)

  POP_SUB(X(broyden_extrapolation_aux))
endsubroutine X(broyden_extrapolation_aux)

!--------------------------------------------------------------------

subroutine X(mixing_diis)(this, vin, vout, vnew, iter)
  type(mix_t), intent(inout) :: this
  R_TYPE,      intent(in)    :: vin(:, :, :)
  R_TYPE,      intent(in)    :: vout(:, :, :)
  R_TYPE,      intent(out)   :: vnew(:, :, :)
  integer,     intent(in)    :: iter

  integer :: size, ii
  integer :: d1, d2, d3
  R_TYPE :: sumalpha
  R_TYPE, allocatable :: aa(:, :), alpha(:), rhs(:)

  PUSH_SUB(X(mixing_diis))

  d1 = this%mixfield%d1
  d2 = this%mixfield%d2
  d3 = this%mixfield%d3

  if (iter <= this%mixfield%d4) then

    size = iter

  else

    size = this%mixfield%d4

    do ii = 2, size
      call lalg_copy(d1, d2, d3, this%mixfield%X(dv)(:, :, :, ii), this%mixfield%X(dv)(:, :, :, ii - 1))
      call lalg_copy(d1, d2, d3, this%mixfield%X(df)(:, :, :, ii), this%mixfield%X(df)(:, :, :, ii - 1))
    end do

  end if

  call lalg_copy(d1, d2, d3, vin, this%mixfield%X(dv)(:, :, :, size))
  this%mixfield%X(df)(1:d1, 1:d2, 1:d3, size) = vout(1:d1, 1:d2, 1:d3) - vin(1:d1, 1:d2, 1:d3)

  if (iter == 1 .or. mod(iter, this%interval) /= 0) then

    vnew(1:d1, 1:d2, 1:d3) = (M_ONE - this%coeff)*vin(1:d1, 1:d2, 1:d3) &
      + this%coeff*vout(1:d1, 1:d2, 1:d3)

    POP_SUB(X(mixing_diis))
    return
  end if


  SAFE_ALLOCATE(aa(1:size + 1, 1:size + 1))
  SAFE_ALLOCATE(alpha(1:size + 1))
  SAFE_ALLOCATE(rhs(1:size + 1))

  call X(mixing_build_matrix)(this, this%mixfield%X(df), size, d2, d3, M_ONE, aa)

  aa(1:size, size + 1) = -M_ONE
  aa(size + 1, 1:size) = -M_ONE
  aa(size + 1, size + 1) = M_ZERO

  rhs(1:size) = M_ZERO
  rhs(size + 1) = -M_ONE

  call lalg_least_squares(size + 1, aa, rhs, alpha, preserve_mat=.false.)

  sumalpha = sum(alpha(1:size))
  alpha = alpha/sumalpha

  vnew(1:d1, 1:d2, 1:d3) = M_ZERO

  do ii = 1, size
    vnew(1:d1, 1:d2, 1:d3) = vnew(1:d1, 1:d2, 1:d3) &
      + alpha(ii)*(this%mixfield%X(dv)(1:d1, 1:d2, 1:d3, ii) &
      + this%residual_coeff*this%mixfield%X(df)(1:d1, 1:d2, 1:d3, ii))
  end do

  POP_SUB(X(mixing_diis))
end subroutine X(mixing_diis)

! --------------------------------------------------------
! Guaranteed-reduction Pulay
! ---------------------------------------------------------
subroutine X(mixing_grpulay)(smix, vin, vout, vnew, iter)
  type(mix_t), intent(inout) :: smix
  integer,      intent(in)   :: iter
  R_TYPE,        intent(in)   :: vin(:, :, :), vout(:, :, :)
  R_TYPE,        intent(out)  :: vnew(:, :, :)

  integer :: d1, d2, d3
  integer :: ipos, iter_used
  R_TYPE, allocatable :: f(:, :, :)

  PUSH_SUB(X(mixing_grpulay))

  d1 = smix%mixfield%d1
  d2 = smix%mixfield%d2
  d3 = smix%mixfield%d3

  SAFE_ALLOCATE(f(1:d1, 1:d2, 1:d3))
  f(1:d1, 1:d2, 1:d3) = vout(1:d1, 1:d2, 1:d3) - vin(1:d1, 1:d2, 1:d3)

  ! we only extrapolate a new vector every two iterations
  select case (mod(iter, 2))
  case (1)
    ! Store df and dv from current iteration
    if (iter > 1) then
      ipos = smix%last_ipos
      call lalg_copy(d1, d2, d3, f(:, :, :), smix%mixfield%X(df)(:, :, :, ipos))
      call lalg_copy(d1, d2, d3, vin(:, :, :), smix%mixfield%X(dv)(:, :, :, ipos))
      call lalg_axpy(d1, d2, d3, -M_ONE, smix%mixfield%X(f_old)(:, :, :), smix%mixfield%X(df)(:, :, :, ipos))
      call lalg_axpy(d1, d2, d3, -M_ONE, smix%mixfield%X(vin_old)(:, :, :), smix%mixfield%X(dv)(:, :, :, ipos))
    end if

    ! Store residual and vin for next extrapolation
    smix%mixfield%X(vin_old) = vin
    smix%mixfield%X(f_old) = f

    ! we need the output vector for vout. So let`s do vnew = vout to get that information
    vnew = vout
  case (0)
    ! Store df and dv from current iteration in arrays df and dv so that we can use them
    ! for the extrapolation. Next iterations they will be lost.
    ipos = mod(smix%last_ipos, smix%ns + 1) + 1
    call lalg_copy(d1, d2, d3, f(:, :, :), smix%mixfield%X(df)(:, :, :, ipos))
    call lalg_copy(d1, d2, d3, vin(:, :, :), smix%mixfield%X(dv)(:, :, :, ipos))
    call lalg_axpy(d1, d2, d3, -M_ONE, smix%mixfield%X(f_old)(:, :, :), smix%mixfield%X(df)(:, :, :, ipos))
    call lalg_axpy(d1, d2, d3, -M_ONE, smix%mixfield%X(vin_old)(:, :, :), smix%mixfield%X(dv)(:, :, :, ipos))

    smix%last_ipos = ipos

    ! extrapolate new vector
    iter_used = min(iter/2, smix%ns + 1)
    call X(pulay_extrapolation)(smix, d2, d3, vin, vnew, iter_used, f, &
      smix%mixfield%X(df)(1:d1, 1:d2, 1:d3, 1:iter_used), smix%mixfield%X(dv)(1:d1, 1:d2, 1:d3, 1:iter_used))

  end select

  SAFE_DEALLOCATE_A(f)

  POP_SUB(X(mixing_grpulay))
end subroutine X(mixing_grpulay)


! ---------------------------------------------------------
subroutine X(pulay_extrapolation)(this, d2, d3, vin, vnew, iter_used, f, df, dv)
  type(mix_t), intent(in) :: this
  integer, intent(in) :: d2, d3
  integer, intent(in)   :: iter_used
  R_TYPE,  intent(in)  :: vin(:, :, :), f(:, :, :), df(:, :, :, :), dv(:, :, :, :)
  R_TYPE,  intent(out) :: vnew(:, :, :)

  integer :: i, j, k, l
  R_TYPE :: alpha
  R_TYPE, allocatable :: a(:, :)

  PUSH_SUB(X(pulay_extrapolation))

  SAFE_ALLOCATE(a(1:iter_used, 1:iter_used))

  ! set matrix A
  call X(mixing_build_matrix)(this, df, iter_used, d2, d3, M_ONE, a)

  call lalg_inverter(iter_used, a)

  ! compute new vector
  vnew = vin
  do i = 1, iter_used
    alpha = M_ZERO
    do j = 1, iter_used
      do l = 1, d3
        do k = 1, d2
          alpha = alpha - a(i, j)*X(mix_dotp)(this, df(:, k, l, j), f(:, k, l), reduce = .false.)
        end do
      end do
    end do
    if (this%der%mesh%parallel_in_domains) call this%der%mesh%allreduce(alpha)
    vnew(:, :, :) = vnew(:, :, :) + alpha * dv(:, :, :, i)
  end do

  SAFE_DEALLOCATE_A(a)

  POP_SUB(X(pulay_extrapolation))
end subroutine X(pulay_extrapolation)

! --------------------------------------------------------------
! In all mixing schemes apart from the linear mixing one, we need to build the same matrix
! This is done in this routine, together with the stabilization scheme proposed by Johnson (1988)
subroutine X(mixing_build_matrix)(this, df, size, d2, d3, ww, beta)
  type(mix_t),       intent(in)  :: this
  R_TYPE,            intent(in)  :: df(:,:,:,:)
  integer,           intent(in)  :: size, d2, d3
  FLOAT,             intent(in)  :: ww
  R_TYPE,            intent(out) :: beta(:,:)

  integer :: i, j, k, l
  FLOAT :: w0

  PUSH_SUB(X(mixing_build_matrix))

  beta = M_ZERO
  do i = 1, size
    do j = i, size
      beta(i, j) = M_ZERO
      do l = 1, d3
        do k = 1, d2
          beta(i, j) = beta(i, j) + ww*X(mix_dotp)(this, df(:, k, l, j), df(:, k, l, i), reduce = .false.)
        end do
      end do
      if (i /= j) beta(j, i) = R_CONJ(beta(i, j))
    end do
  end do

  if (this%der%mesh%parallel_in_domains) call this%der%mesh%allreduce(beta)

  ! According to Johnson (1988), w0 is chosen as 0.01. Because we do not
  ! normalize the components, we need to choose w0 differently. Its purpose
  ! is to stabilize the inversion by adding a small number to the diagonal.
  ! Thus we compute w0 as 0.01 of the minimum of the values on the diagonal
  ! to improve the numerical stability. This enables convergence to very
  ! high accuracies.
  !
  ! NTD: I changed the logic to be 0.01^2 times the smallest diagonal value.
  ! Else, for a small value of 1e-10, we were getting w0 to be 1e-24
  ! For 1e-6, we were getting w0=1e-16.
  ! Overall, this was having no practical effect apart for large numbers.
  w0 = M_HUGE
  do i = 1, size
    w0 = min(w0, abs(beta(i, i)))
  end do
  w0 = w0 * CNST(1e-4)
  ! safeguard if w0 should be exactly zero or underflow
  if (w0 < CNST(1e-150)) then
    w0 = CNST(1e-4)
  end if
  do i = 1, size
    beta(i, i) = w0 + beta(i, i)
  end do

  POP_SUB(X(mixing_build_matrix))
end subroutine X(mixing_build_matrix)

! --------------------------------------------------------------

R_TYPE function X(mix_dotp)(this, xx, yy, reduce) result(dotp)
  type(mix_t),       intent(in) :: this
  R_TYPE,            intent(in) :: xx(:)
  R_TYPE,            intent(in) :: yy(:)
  logical, optional, intent(in) :: reduce

  R_TYPE, allocatable :: ff(:), precff(:)
  logical :: reduce_

  PUSH_SUB(X(mix_dotp))

  reduce_ = .true.
  if (present(reduce)) reduce_ = reduce

  ASSERT(this%mixfield%d1 == this%der%mesh%np)

  if (this%precondition) then

    SAFE_ALLOCATE(ff(1:this%der%mesh%np_part))
    SAFE_ALLOCATE(precff(1:this%der%mesh%np))

    ff(1:this%der%mesh%np) = yy(1:this%der%mesh%np)
    call X(derivatives_perform)(this%preconditioner, this%der, ff, precff)
    dotp = X(mf_dotp)(this%der%mesh, xx, precff, reduce = reduce_)

    SAFE_DEALLOCATE_A(precff)
    SAFE_DEALLOCATE_A(ff)

  else
    dotp = X(mf_dotp)(this%der%mesh, xx, yy, reduce = reduce_)
  end if

  POP_SUB(X(mix_dotp))

end function X(mix_dotp)

! --------------------------------------------------------------
subroutine X(mixfield_set_vin2)(mixfield, vin)
  type(mixfield_t), intent(inout) :: mixfield
  R_TYPE,           intent(in)    :: vin(:, :)

  PUSH_SUB(X(mixfield_set_vin2))

  mixfield%X(vin)(1:mixfield%d1, 1, 1:mixfield%d3) = vin(1:mixfield%d1, 1:mixfield%d3)

  POP_SUB(X(mixfield_set_vin2))
end subroutine X(mixfield_set_vin2)

! --------------------------------------------------------------
subroutine X(mixfield_set_vin3)(mixfield, vin)
  type(mixfield_t), intent(inout) :: mixfield
  R_TYPE,           intent(in)    :: vin(:, :, :)

  PUSH_SUB(X(mixfield_set_vin3))

  call lalg_copy(mixfield%d1, mixfield%d2, mixfield%d3, vin(:, :, :), mixfield%X(vin)(:, :, :))

  POP_SUB(X(mixfield_set_vin3))
end subroutine X(mixfield_set_vin3)

! --------------------------------------------------------------
subroutine X(mixfield_set_vout2)(mixfield, vout)
  type(mixfield_t), intent(inout) :: mixfield
  R_TYPE,           intent(in)    :: vout(:, :)

  PUSH_SUB(X(mixfield_set_vout2))

  mixfield%X(vout)(1:mixfield%d1, 1, 1:mixfield%d3) = vout(1:mixfield%d1, 1:mixfield%d3)

  POP_SUB(X(mixfield_set_vout2))
end subroutine X(mixfield_set_vout2)

! --------------------------------------------------------------
subroutine X(mixfield_set_vout3)(mixfield, vout)
  type(mixfield_t), intent(inout) :: mixfield
  R_TYPE,           intent(in)    :: vout(:, :, :)

  PUSH_SUB(X(mixfield_set_vout3))

  call lalg_copy(mixfield%d1, mixfield%d2, mixfield%d3, vout(:, :, :), mixfield%X(vout)(:, :, :))

  POP_SUB(X(mixfield_set_vout3))
end subroutine X(mixfield_set_vout3)

! --------------------------------------------------------------
subroutine X(mixfield_get_vnew)(mixfield, vnew)
  type(mixfield_t), intent(in)    :: mixfield
  R_TYPE,           intent(inout) :: vnew(:,:)

  PUSH_SUB(X(mixfield_get_vnew))

  vnew(1:mixfield%d1, 1:mixfield%d3) = mixfield%X(vnew)(1:mixfield%d1, 1, 1:mixfield%d3)

  POP_SUB(X(mixfield_get_vnew))
end subroutine X(mixfield_get_vnew)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
