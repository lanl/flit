!
! © 2024. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by
! Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by
! Triad National Security, LLC, and the U.S. Department of Energy/National
! Nuclear Security Administration. The Government is granted for itself and
! others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare. derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!


module libflit_readpar

    use libflit_string
    use libflit_utility
    use libflit_io
    use iso_fortran_env

    implicit none

    private

    integer, parameter :: nstring_max = 100

    ! To improve readness of an end-user code,
    ! here I still make the subroutines publicly accessible
    interface readpar_int
        module procedure :: readpar_int2
        module procedure :: readpar_int4
        module procedure :: readpar_int8
    end interface

    interface readpar_nint
        module procedure :: readnpar_int2
        module procedure :: readnpar_int4
        module procedure :: readnpar_int8
    end interface readpar_nint

    interface readpar_xint
        module procedure :: readxpar_int2
        module procedure :: readxpar_int4
        module procedure :: readxpar_int8
    end interface

    interface readpar_nfloat
        module procedure :: readnpar_float
    end interface
    interface readpar_xfloat
        module procedure :: readxpar_float
    end interface

    interface readpar_ndouble
        module procedure :: readnpar_double
    end interface
    interface readpar_xdouble
        module procedure :: readxpar_double
    end interface

    interface readpar_ncomplex
        module procedure :: readnpar_complex
    end interface
    interface readpar_xcomplex
        module procedure :: readxpar_complex
    end interface

    interface readpar_ndcomplex
        module procedure :: readnpar_dcomplex
    end interface
    interface readpar_xdcomplex
        module procedure :: readxpar_dcomplex
    end interface

    interface readpar_nlogical
        module procedure :: readnpar_logical
    end interface
    interface readpar_xlogical
        module procedure :: readxpar_logical
    end interface

    interface readpar_nstring
        module procedure :: readnpar_string
    end interface
    interface readpar_xstring
        module procedure :: readxpar_string
    end interface

    interface getpar_int
        module procedure :: getpar_int2
        module procedure :: getpar_int4
        module procedure :: getpar_int8
    end interface

    interface getpar_nint
        module procedure :: getnpar_int2
        module procedure :: getnpar_int4
        module procedure :: getnpar_int8
    end interface getpar_nint

    interface getpar_nfloat
        module procedure :: getnpar_float
    end interface

    interface getpar_ndouble
        module procedure :: getnpar_double
    end interface

    interface getpar_ncomplex
        module procedure :: getnpar_complex
    end interface

    interface getpar_ndcomplex
        module procedure :: getnpar_dcomplex
    end interface

    interface getpar_nlogical
        module procedure :: getnpar_logical
    end interface

    interface getpar_nstring
        module procedure :: getnpar_string
    end interface

    interface getpar_xint
        module procedure :: getxpar_int2
        module procedure :: getxpar_int4
        module procedure :: getxpar_int8
    end interface getpar_xint

    interface getpar_xfloat
        module procedure :: getxpar_float
    end interface

    interface getpar_xdouble
        module procedure :: getxpar_double
    end interface

    interface getpar_xcomplex
        module procedure :: getxpar_complex
    end interface

    interface getpar_xdcomplex
        module procedure :: getxpar_dcomplex
    end interface

    interface getpar_xlogical
        module procedure :: getxpar_logical
    end interface

    interface getpar_xstring
        module procedure :: getxpar_string
    end interface

    interface parsepar_int
        module procedure :: parsepar_int2
        module procedure :: parsepar_int4
        module procedure :: parsepar_int8
    end interface

    interface parsepar_nint
        module procedure :: parsenpar_int2
        module procedure :: parsenpar_int4
        module procedure :: parsenpar_int8
    end interface parsepar_nint

    interface parsepar_nfloat
        module procedure :: parsenpar_float
    end interface

    interface parsepar_ndouble
        module procedure :: parsenpar_double
    end interface

    interface parsepar_ncomplex
        module procedure :: parsenpar_complex
    end interface

    interface parsepar_ndcomplex
        module procedure :: parsenpar_dcomplex
    end interface

    interface parsepar_nlogical
        module procedure :: parsenpar_logical
    end interface

    interface parsepar_nstring
        module procedure :: parsenpar_string
    end interface

    interface parsepar_xint
        module procedure :: parsexpar_int2
        module procedure :: parsexpar_int4
        module procedure :: parsexpar_int8
    end interface parsepar_xint

    interface parsepar_xfloat
        module procedure :: parsexpar_float
    end interface

    interface parsepar_xdouble
        module procedure :: parsexpar_double
    end interface

    interface parsepar_xcomplex
        module procedure :: parsexpar_complex
    end interface

    interface parsepar_xdcomplex
        module procedure :: parsexpar_dcomplex
    end interface

    interface parsepar_xlogical
        module procedure :: parsexpar_logical
    end interface

    interface parsepar_xstring
        module procedure :: parsexpar_string
    end interface

    public :: readpar_int
    public :: readpar_nint
    public :: readpar_xint
    public :: readpar_float
    public :: readpar_nfloat
    public :: readpar_xfloat
    public :: readpar_double
    public :: readpar_ndouble
    public :: readpar_xdouble
    public :: readpar_complex
    public :: readpar_ncomplex
    public :: readpar_xcomplex
    public :: readpar_dcomplex
    public :: readpar_ndcomplex
    public :: readpar_xdcomplex
    public :: readpar_logical
    public :: readpar_nlogical
    public :: readpar_xlogical
    public :: readpar_string
    public :: readpar_nstring
    public :: readpar_xstring

    public :: getpar_int
    public :: getpar_nint
    public :: getpar_xint
    public :: getpar_float
    public :: getpar_nfloat
    public :: getpar_xfloat
    public :: getpar_double
    public :: getpar_ndouble
    public :: getpar_xdouble
    public :: getpar_complex
    public :: getpar_ncomplex
    public :: getpar_xcomplex
    public :: getpar_ndcomplex
    public :: getpar_xdcomplex
    public :: getpar_logical
    public :: getpar_nlogical
    public :: getpar_xlogical
    public :: getpar_string
    public :: getpar_nstring
    public :: getpar_xstring

    public :: parsepar_int
    public :: parsepar_nint
    public :: parsepar_xint
    public :: parsepar_float
    public :: parsepar_nfloat
    public :: parsepar_xfloat
    public :: parsepar_double
    public :: parsepar_ndouble
    public :: parsepar_xdouble
    public :: parsepar_complex
    public :: parsepar_ncomplex
    public :: parsepar_xcomplex
    public :: parsepar_ndcomplex
    public :: parsepar_xdcomplex
    public :: parsepar_logical
    public :: parsepar_nlogical
    public :: parsepar_xlogical
    public :: parsepar_string
    public :: parsepar_nstring
    public :: parsepar_xstring

    public :: checkpar
    public :: checkarg

contains

    !
    !> Check if a parameter file contains any invalid parameter
    !
    subroutine checkpar(filename, valid_parnames)

        character(len=*), intent(in) :: filename
        character(len=*), dimension(:), intent(in) :: valid_parnames

        character(len=1024) :: line
        integer :: funit, eqindex
        integer :: ioerr
        logical :: invalid

        open (newunit=funit, file=tidy(filename), status='old', action='read')
        do
            read (funit, '(a)', iostat=ioerr) line
            if (ioerr /= 0) exit
            if (tidy(line) == 'exit') exit
            if (len(tidy(line)) > 0) then
                eqindex = index(line, '=')
                if (eqindex >= 1 .and. .not. any(to_lower(tidy(line(1:eqindex - 1))) == valid_parnames)) then
                    write (error_unit, *) ' Invalid parameter: '//to_lower(tidy(line(1:eqindex - 1)))
                    invalid = .true.
                end if
            end if
        end do
        close (funit)

        if (invalid) then
            write (error_unit, *) ' Parameter file contains invalid parameters. Exit.'
            stop
        end if

    end subroutine checkpar

    !
    !> Check if a list of command-line arguments contains any invalid parameter
    !
    subroutine checkarg(valid_parnames)

        character(len=*), dimension(:), intent(in) :: valid_parnames

        integer :: i
        integer :: len_arg, eqindex
        character(len=256) :: arg
        logical :: invalid

        invalid = .false.
        do i = 1, command_argument_count()

            call get_command_argument(i, arg)
            len_arg = len_trim(arg)

            if (len_arg > 1) then
                eqindex = index(arg, '=')
                if (.not. any(tidy(arg(1:eqindex - 1)) == valid_parnames)) then
                    write(error_unit, *) ' Invalid argument: '//tidy(arg(1:eqindex - 1))
                    invalid = .true.
                end if
            end if

        end do

        if(invalid) then
            write(error_unit, *) ' Argument list contains invalid arguments. Exit.'
            stop
        end if

    end subroutine checkarg

#define T int2
#define TT integer(kind=2)
#define TTT integer(kind=2)
#define extract_TT extract_int2
#define extract_nTT extract_nint2
#define extract_xTT extract_xint2
#define nearest_interp
#include "template_readpar.f90"

#define T int4
#define TT integer(kind=4)
#define TTT integer(kind=4)
#define extract_TT extract_int4
#define extract_nTT extract_nint4
#define extract_xTT extract_xint4
#define nearest_interp
#include "template_readpar.f90"

#define T int8
#define TT integer(kind=8)
#define TTT integer(kind=8)
#define extract_TT extract_int8
#define extract_nTT extract_nint8
#define extract_xTT extract_xint8
#define nearest_interp
#include "template_readpar.f90"

#define T float
#define TT real
#define TTT real
#define extract_TT extract_float
#define extract_nTT extract_nfloat
#define extract_xTT extract_xfloat
#define linear_interp
#include "template_readpar.f90"

#define T double
#define TT double precision
#define TTT double precision
#define extract_TT extract_double
#define extract_nTT extract_ndouble
#define extract_xTT extract_xdouble
#define linear_interp
#include "template_readpar.f90"

#define T complex
#define TT complex
#define TTT complex
#define extract_TT extract_complex
#define extract_nTT extract_ncomplex
#define extract_xTT extract_xcomplex
#define linear_interp
#include "template_readpar.f90"

#define T dcomplex
#define TT double complex
#define TTT double complex
#define extract_TT extract_dcomplex
#define extract_nTT extract_ndcomplex
#define extract_xTT extract_xdcomplex
#define linear_interp
#include "template_readpar.f90"

#define T logical
#define TT logical
#define TTT logical
#define extract_TT extract_logical
#define extract_nTT extract_nlogical
#define extract_xTT extract_xlogical
#define nearest_interp
#include "template_readpar.f90"

#define T string
#define TT character(len=*)
#define TTT character(len=64)
#define extract_TT extract_string
#define extract_nTT extract_nstring
#define extract_xTT extract_xstring
#define nearest_interp
#include "template_readpar.f90"

end module libflit_readpar
