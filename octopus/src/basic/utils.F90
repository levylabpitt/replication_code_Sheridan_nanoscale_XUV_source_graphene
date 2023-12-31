!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2021 S. Ohlmann
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

!> This module is intended to contain simple general-purpose utility functions
!! and procedures.

module utils_oct_m
  use debug_oct_m
  use global_oct_m
  use io_oct_m
  use loct_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use string_oct_m

  implicit none

  private
  public ::                       &
    get_divisors,                 &
    index2axis,                   &
    index2axisBZ,                 &
    output_tensor,                &
    output_dipole,                &
    print_header,                 &
    leading_dimension_is_known,   &
    lead_dim,                     &
    get_config_opts,              &
    get_optional_libraries,       &
    make_array_larger

  interface leading_dimension_is_known
    module procedure dleading_dimension_is_known, zleading_dimension_is_known
    module procedure dleading_dimension_is_known2, zleading_dimension_is_known2
  end interface leading_dimension_is_known

  interface lead_dim
    module procedure dlead_dim, zlead_dim
    module procedure dlead_dim2, zlead_dim2
  end interface lead_dim

contains

  ! ---------------------------------------------------------
  subroutine get_divisors(nn, n_divisors, divisors)
    integer, intent(in)    :: nn
    integer, intent(inout) :: n_divisors
    integer, intent(out)   :: divisors(:)

    integer :: ii, max_d

    PUSH_SUB(get_divisors)

    ASSERT(n_divisors > 1)
    max_d = n_divisors

    n_divisors = 1
    divisors(n_divisors) = 1
    do ii = 2, nn / 2
      if (mod(nn, ii) == 0) then
        n_divisors = n_divisors + 1

        if (n_divisors > max_d - 1) then
          message(1) = "Internal error in get_divisors. Please increase n_divisors"
          call messages_fatal(1)
        end if

        divisors(n_divisors) = ii
      end if
    end do
    n_divisors = n_divisors + 1
    divisors(n_divisors) = nn

    POP_SUB(get_divisors)
  end subroutine get_divisors


  ! ---------------------------------------------------------
  character pure function index2axis(idir) result(ch)
    integer, intent(in) :: idir

    select case (idir)
    case (1)
      ch = 'x'
    case (2)
      ch = 'y'
    case (3)
      ch = 'z'
    case (4)
      ch = 'w'
    case default
      write(ch,'(i1)') idir
    end select

  end function index2axis

  ! ---------------------------------------------------------
  pure function index2axisBZ(idir) result(ch)
    integer, intent(in) :: idir
    character(len=2) :: ch

    select case (idir)
    case (1)
      ch = "kx"
    case (2)
      ch = "ky"
    case (3)
      ch = "kz"
    case (4)
      ch = "kw"
    case default
      write(ch,'(i2)') idir
    end select

  end function index2axisBZ


  ! ---------------------------------------------------------
  subroutine output_tensor(tensor, ndim, unit, write_average, iunit, namespace)
    FLOAT,                       intent(in) :: tensor(:,:)
    integer,                     intent(in) :: ndim
    type(unit_t),                intent(in) :: unit
    logical,           optional, intent(in) :: write_average
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    FLOAT :: trace
    integer :: jj, kk
    logical :: write_average_

    PUSH_SUB(output_tensor)

    write_average_ = optional_default(write_average, .true.)

    trace = M_ZERO
    message(1) = ""
    do jj = 1, ndim
      do kk = 1, ndim
        write(message(1), '(a,f20.6)') trim(message(1)), units_from_atomic(unit, tensor(jj, kk))
      end do
      trace = trace + tensor(jj, jj)
      call messages_info(1, iunit=iunit, namespace=namespace)
    end do

    if (write_average_) then
      write(message(1), '(a, f20.6)') 'Isotropic average', units_from_atomic(unit, trace/TOFLOAT(ndim))
      call messages_info(1, iunit=iunit, namespace=namespace)
    end if

    POP_SUB(output_tensor)
  end subroutine output_tensor


  ! ---------------------------------------------------------
  subroutine output_dipole(dipole, ndim, iunit, namespace)
    FLOAT,                       intent(in) :: dipole(:)
    integer,                     intent(in) :: ndim
    integer,           optional, intent(in) :: iunit
    type(namespace_t), optional, intent(in) :: namespace

    integer :: idir

    PUSH_SUB(output_dipole)

    write(message(1), '(a,a20,a17)') 'Dipole:', '[' // trim(units_abbrev(units_out%length)) // ']', &
      '[' // trim(units_abbrev(unit_debye)) // ']'
    do idir = 1, ndim
      write(message(1+idir), '(6x,3a,es14.5,3x,2es14.5)') '<', index2axis(idir), '> = ', &
        units_from_atomic(units_out%length, dipole(idir)), units_from_atomic(unit_debye, dipole(idir))
    end do
    call messages_info(1+ndim, iunit=iunit, namespace=namespace)

    POP_SUB(output_dipole)
  end subroutine output_dipole

  !> This subroutine prints the logo followed by information about
  !! the compilation and the system. It also prints the start time
  !! of the execution.
  ! ---------------------------------------------------------
  subroutine print_header()
#ifdef HAVE_FC_COMPILER_VERSION
    use iso_fortran_env
#endif

    character(len=256) :: sys_name

    ! Let us print our logo
    call io_dump_file(stdout, trim(trim(conf%share) // '/logo'))

    ! Let us print the version
    message(1) = ""
    message(2) = str_center("Running octopus", 70)
    message(3) = ""
    call messages_info(3)

    message(1) = &
      "Version                : " // trim(conf%version)
    message(2) = &
      "Commit                 : "// trim(conf%git_commit)
    message(3) = &
      "Configuration time     : "// trim(conf%config_time)
    call messages_info(3)

    message(1) = 'Configuration options  : ' // trim(get_config_opts())
    message(2) = 'Optional libraries     :'  // trim(get_optional_libraries())

    message(3) = 'Architecture           : ' + TOSTRING(OCT_ARCH)
    call messages_info(3)

    message(1) = &
      "C compiler             : "//trim(conf%cc)
    message(2) = &
      "C compiler flags       : "//trim(conf%cflags)
    message(3) = &
      "C++ compiler           : "//trim(conf%cxx)
    message(4) = &
      "C++ compiler flags     : "//trim(conf%cxxflags)
#ifdef HAVE_FC_COMPILER_VERSION
    message(5) = "Fortran compiler       : "//trim(conf%fc) //" ("//compiler_version()//")"
#else
    message(5) = "Fortran compiler       : "//trim(conf%fc)
#endif
    message(6) = &
      "Fortran compiler flags : "//trim(conf%fcflags)
    call messages_info(6)

    message(1) = ""
    call messages_info(1)

    ! Let us print where we are running
    call loct_sysname(sys_name)
    write(message(1), '(a)') str_center("The octopus is swimming in " // trim(sys_name), 70)
    message(2) = ""
    call messages_info(2)

    call mpi_world%barrier()

    call print_date("Calculation started on ")
  end subroutine print_header

  character(len=256) function get_config_opts()

    write(get_config_opts, '(a, i1)') 'maxdim', MAX_DIM
#ifdef HAVE_OPENMP
    get_config_opts = trim(get_config_opts)//' openmp'
#endif
#ifdef HAVE_MPI
    get_config_opts = trim(get_config_opts)//' mpi'
#endif
#ifdef HAVE_OPENCL
    get_config_opts = trim(get_config_opts)//' opencl'
#endif
#ifdef HAVE_CUDA
    get_config_opts = trim(get_config_opts)//' cuda'
#endif
#ifdef HAVE_M128D
    get_config_opts = trim(get_config_opts)//' sse2'
#endif
#ifdef HAVE_M256D
    get_config_opts = trim(get_config_opts)//' avx'
#endif
#ifdef HAVE_BLUE_GENE
    get_config_opts = trim(get_config_opts)//' bluegene/p'
#endif
#ifdef HAVE_BLUE_GENE_Q
    get_config_opts = trim(get_config_opts)//' bluegene/q'
#endif
#ifdef HAVE_LIBXC5
    get_config_opts = trim(get_config_opts)//' libxc5'
#endif
#ifdef HAVE_LIBXC_FXC
    get_config_opts = trim(get_config_opts)//' libxc_fxc'
#endif
#ifdef HAVE_LIBXC_KXC
    get_config_opts = trim(get_config_opts)//' libxc_kxc'
#endif

  end function get_config_opts

  character(len=256) function get_optional_libraries()

    ! keep in alphabetical order, for ease in seeing if something is listed
    get_optional_libraries = ''
#ifdef HAVE_BERKELEYGW
    get_optional_libraries = trim(get_optional_libraries)//' berkeleygw'
#endif
#ifdef HAVE_CGAL
    get_optional_libraries = trim(get_optional_libraries)//' cgal'
#endif
#ifdef HAVE_CLFFT
    get_optional_libraries = trim(get_optional_libraries)//' clamdfft'
#endif
#ifdef HAVE_CLBLAS
    get_optional_libraries = trim(get_optional_libraries)//' clblas'
#endif
#ifdef HAVE_DFTBPLUS
    get_optional_libraries = trim(get_optional_libraries)//' DFTBPlus'
#endif
#ifdef HAVE_DFTBPLUS_DEVEL
    get_optional_libraries = trim(get_optional_libraries)//' DFTBPlus_devel'
#endif
#ifdef HAVE_ELPA
    get_optional_libraries = trim(get_optional_libraries)//' ELPA'
#endif
#ifdef HAVE_ETSF_IO
    get_optional_libraries = trim(get_optional_libraries)//' etsf_io'
#endif
#ifdef HAVE_GDLIB
    get_optional_libraries = trim(get_optional_libraries)//' gdlib'
#endif
#ifdef HAVE_LIBFM
    get_optional_libraries = trim(get_optional_libraries)//' libfm'
#endif
#if (defined HAVE_LIBISF) || (defined HAVE_PSOLVER)
    get_optional_libraries = trim(get_optional_libraries)//' psolver'
#endif
#ifdef HAVE_LIBVDWXC
    get_optional_libraries = trim(get_optional_libraries)//' libvdwxc'
#endif
#ifdef HAVE_METIS
    get_optional_libraries = trim(get_optional_libraries)//' metis'
#endif
#ifdef HAVE_NETCDF
    get_optional_libraries = trim(get_optional_libraries)//' netcdf'
#endif
#ifdef HAVE_NFFT
    get_optional_libraries = trim(get_optional_libraries)//' nfft'
#endif
#ifdef HAVE_PARMETIS
    get_optional_libraries = trim(get_optional_libraries)//' parmetis'
#endif
#ifdef HAVE_PFFT
    get_optional_libraries = trim(get_optional_libraries)//' pfft'
#endif
#ifdef HAVE_PNFFT
    get_optional_libraries = trim(get_optional_libraries)//' pnfft'
#endif
#ifdef HAVE_PSPIO
    get_optional_libraries = trim(get_optional_libraries)//' pspio'
#endif
#ifdef HAVE_SCALAPACK
    get_optional_libraries = trim(get_optional_libraries)//' scalapack'
#endif
#ifdef HAVE_SPARSKIT
    get_optional_libraries = trim(get_optional_libraries)//' sparskit'
#endif
#ifdef HAVE_NLOPT
    get_optional_libraries = trim(get_optional_libraries)//' nlopt'
#endif

  end function get_optional_libraries

  ! ---------------------------------------------------------

  logical function dleading_dimension_is_known(array) result(known)
    real(r8), intent(in) :: array(:, :)

    known = .true.

#if defined(HAVE_FORTRAN_LOC) && defined(HAVE_FC_SIZEOF)
    if (ubound(array, dim = 2) > 1) then
      known = ubound(array, dim = 1) == (loc(array(1, 2)) - loc(array(1, 1)))/sizeof(array(1, 1))
    end if
#endif

  end function dleading_dimension_is_known


  ! ---------------------------------------------------------

  logical function zleading_dimension_is_known(array) result(known)
    complex(r8), intent(in) :: array(:, :)

    known = .true.

#if defined(HAVE_FORTRAN_LOC) && defined(HAVE_FC_SIZEOF)
    if (ubound(array, dim = 2) > 1) then
      known = ubound(array, dim = 1) == (loc(array(1, 2)) - loc(array(1, 1)))/sizeof(array(1, 1))
    end if
#endif

  end function zleading_dimension_is_known

  ! ---------------------------------------------------------

  logical function dleading_dimension_is_known2(array) result(known)
    real(r8), intent(in) :: array(:, :, :)

    known = .true.

#if defined(HAVE_FORTRAN_LOC) && defined(HAVE_FC_SIZEOF)
    if (ubound(array, dim = 2) > 1) then
      known = ubound(array, dim = 1) == (loc(array(1, 2, 1)) - loc(array(1, 1, 1)))/sizeof(array(1, 1, 1))
    end if
#endif

  end function dleading_dimension_is_known2


  ! ---------------------------------------------------------

  logical function zleading_dimension_is_known2(array) result(known)
    complex(r8), intent(in) :: array(:, :, :)

    known = .true.

#if defined(HAVE_FORTRAN_LOC) && defined(HAVE_FC_SIZEOF)
    if (ubound(array, dim = 2) > 1) then
      known = ubound(array, dim = 1) == (loc(array(1, 2, 1)) - loc(array(1, 1, 1)))/sizeof(array(1, 1, 1))
    end if
#endif

  end function zleading_dimension_is_known2

  ! ---------------------------------------------------------

  integer function dlead_dim(array) result(lead_dim)
    real(r8), intent(in) :: array(:, :)

    ASSERT(leading_dimension_is_known(array))

    lead_dim = ubound(array, dim = 1)
  end function dlead_dim

  ! ---------------------------------------------------------

  integer function zlead_dim(array) result(lead_dim)
    complex(r8), intent(in) :: array(:, :)

    ASSERT(leading_dimension_is_known(array))

    lead_dim = ubound(array, dim = 1)
  end function zlead_dim

  ! ---------------------------------------------------------

  integer function dlead_dim2(array) result(lead_dim)
    real(r8), intent(in) :: array(:, :, :)

    ASSERT(leading_dimension_is_known(array))

    lead_dim = ubound(array, dim = 1) * ubound(array, dim = 2)
  end function dlead_dim2

  ! ---------------------------------------------------------

  integer function zlead_dim2(array) result(lead_dim)
    complex(r8), intent(in) :: array(:, :, :)

    ASSERT(leading_dimension_is_known(array))

    lead_dim = ubound(array, dim = 1) * ubound(array, dim = 2)
  end function zlead_dim2

  subroutine make_array_larger(array, new_size)
    integer(i8), allocatable, intent(inout) :: array(:)
    integer,                  intent(in)    :: new_size

    integer(i8), allocatable :: tmp(:)
    integer :: copy_size

    PUSH_SUB(make_array_larger)

    SAFE_ALLOCATE(tmp(1:new_size))
    copy_size = min(new_size, size(array))
    tmp(1:copy_size) = array(1:copy_size)
    SAFE_DEALLOCATE_A(array)
    call move_alloc(tmp, array)

    POP_SUB(make_array_larger)
  end subroutine make_array_larger

end module utils_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
