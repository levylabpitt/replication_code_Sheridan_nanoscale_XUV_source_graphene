!! Copyright (C) 2020 Heiko Appel
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

module lorentz_force_oct_m
  use debug_oct_m
  use force_interaction_oct_m
  use global_oct_m
  use interaction_with_partner_oct_m
  use interaction_partner_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use quantity_oct_m

  implicit none

  private
  public ::                &
    lorentz_force_t

  !> Lorenz force between a systems of particles and an electromagnetic field.
  type, extends(force_interaction_t) :: lorentz_force_t
    private
    FLOAT, pointer :: system_charge(:) !< pointer to array storing the charges of the particles
    FLOAT, pointer, public :: system_pos(:,:) !< pointer to array storing the positions of the particles
    FLOAT, pointer :: system_vel(:,:) !< pointer to array storing the velocities of the particles

    FLOAT, allocatable, public :: partner_E_field(:,:) !< E field generated by partner at the positions of the system particles
    FLOAT, allocatable, public :: partner_B_field(:,:) !< B field generated by partner at the positions of the system particles

  contains
    procedure :: init => lorentz_force_init
    procedure :: calculate => lorentz_force_calculate
    procedure :: calculate_energy => lorentz_force_calculate_energy
    final :: lorentz_force_finalize
  end type lorentz_force_t

  interface lorentz_force_t
    module procedure lorentz_force_constructor
  end interface lorentz_force_t

contains

  ! ---------------------------------------------------------
  function lorentz_force_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(lorentz_force_t),               pointer       :: this

    PUSH_SUB(lorentz_force_constructor)

    SAFE_ALLOCATE(this)

    this%label = "lorenz_force"

    this%partner => partner

    ! The Lorentz force needs the position, velocity and charge of the system
    this%n_system_quantities = 3
    SAFE_ALLOCATE(this%system_quantities(1:this%n_system_quantities))
    this%system_quantities(1) = POSITION
    this%system_quantities(2) = VELOCITY
    this%system_quantities(3) = CHARGE
    nullify(this%system_pos)
    nullify(this%system_vel)

    ! The Lorenz force needs the E and B field of the interaction partner at the particle position
    this%n_partner_quantities = 2
    SAFE_ALLOCATE(this%partner_quantities(1:this%n_partner_quantities))
    this%partner_quantities(1) = E_FIELD
    this%partner_quantities(2) = B_FIELD
    this%intra_interaction = .false. ! This interaction does not support intra-interactions


    POP_SUB(lorentz_force_constructor)
  end function lorentz_force_constructor

  ! ---------------------------------------------------------
  subroutine lorentz_force_init(this, dim, system_np, system_quantities, system_charge, system_pos, system_vel, namespace)
    class(lorentz_force_t),               intent(inout) :: this
    integer,                              intent(in)    :: dim !< number of dimensions in space
    integer,                              intent(in)    :: system_np  !< number of particles in the system that owns this interaction
    type(quantity_t),                     intent(inout) :: system_quantities(:)
    FLOAT,                        target, intent(in)    :: system_charge(:)
    FLOAT,                        target, intent(in)    :: system_pos(:,:)
    FLOAT,                        target, intent(in)    :: system_vel(:,:)
    type(namespace_t),                    intent(in)    :: namespace

    PUSH_SUB(lorentz_force_init)

    message(1) = "Energies for Lorentz forces are not yet implemented, and currently set to 0."
    call messages_warning(1, namespace=namespace)

    this%dim = dim
    this%system_np = system_np
    this%energy = M_ZERO
    SAFE_ALLOCATE(this%force(dim, system_np))

    SAFE_ALLOCATE(this%partner_E_field(1:dim, 1:system_np))
    SAFE_ALLOCATE(this%partner_B_field(1:dim, 1:system_np))

    this%system_charge => system_charge
    this%system_pos => system_pos
    this%system_vel => system_vel

    POP_SUB(lorentz_force_init)
  end subroutine lorentz_force_init

  ! ---------------------------------------------------------
  subroutine lorentz_force_calculate(this)
    class(lorentz_force_t),             intent(inout) :: this

    integer :: ip

    PUSH_SUB(lorentz_force_calculate)

    ! Curl is defined only in 3D
    ASSERT(this%dim == 3)
    ASSERT(.not. this%intra_interaction)

    do ip = 1, this%system_np
      this%force(1, ip) = this%partner_E_field(1, ip) + &
        this%system_vel(2, ip)*this%partner_B_field(3, ip) - this%system_vel(3, ip)*this%partner_B_field(2, ip)

      this%force(2, ip) = this%partner_E_field(2, ip) + &
        this%system_vel(3, ip)*this%partner_B_field(1, ip) - this%system_vel(1, ip)*this%partner_B_field(3, ip)

      this%force(3, ip) = this%partner_E_field(3, ip) + &
        this%system_vel(1, ip)*this%partner_B_field(2, ip) - this%system_vel(2, ip)*this%partner_B_field(1, ip)

      this%force(1:3, ip) = this%force(1:3, ip)*this%system_charge(ip)
    end do

    POP_SUB(lorentz_force_calculate)
  end subroutine lorentz_force_calculate

  ! ---------------------------------------------------------
  subroutine lorentz_force_calculate_energy(this)
    class(lorentz_force_t),             intent(inout) :: this

    integer :: ip
    FLOAT   :: power

    PUSH_SUB(lorentz_force_calculate_energy)

    ! the rate at which the energy is transferred from the EM field to the particle
    ! the B field does not contribute any energy to the particles
    power = M_ZERO
    do ip = 1, this%system_np
      power = power - dot_product(this%system_vel(1:3,ip), this%force(1:3,ip))
    end do

    this%energy = M_ZERO
    !TODO: We need to implement the proper time integration of the power.
    !      However, to that end, we first need to implement another clock at which the energies are
    !      calculated.
    ! this%energy = this%energy + power*this%partner%clock%get_time_step()

    POP_SUB(lorentz_force_calculate_energy)
  end subroutine lorentz_force_calculate_energy


  ! ---------------------------------------------------------
  subroutine lorentz_force_finalize(this)
    type(lorentz_force_t), intent(inout) :: this

    PUSH_SUB(lorentz_force_finalize)

    this%force = M_ZERO
    nullify(this%system_charge)
    nullify(this%system_pos)
    nullify(this%system_vel)
    SAFE_DEALLOCATE_A(this%force)
    SAFE_DEALLOCATE_A(this%partner_E_field)
    SAFE_DEALLOCATE_A(this%partner_B_field)

    call interaction_with_partner_end(this)

    POP_SUB(lorentz_force_finalize)
  end subroutine lorentz_force_finalize

end module lorentz_force_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
