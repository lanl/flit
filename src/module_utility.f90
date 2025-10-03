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


module libflit_utility

    use iso_fortran_env
    use libflit_constants
    use libflit_string
    use libflit_array
    use, intrinsic :: ieee_arithmetic, only: ieee_is_normal

    implicit none

    private

    interface warn
        module procedure :: warn_string
        module procedure :: warn_int2
        module procedure :: warn_int4
        module procedure :: warn_float
        module procedure :: warn_double
        module procedure :: warn_cmplx
        module procedure :: warn_logical
    end interface

    interface clip
        module procedure :: clip_int
        module procedure :: clip_float
        module procedure :: clip_double
        module procedure :: clip_complex
    end interface

    interface swap
        module procedure :: swap_int
        module procedure :: swap_float
        module procedure :: swap_double
        module procedure :: swap_complex
        module procedure :: swap_dcomplex
        module procedure :: swap_logical
    end interface

    interface return_normal
        module procedure :: return_normal_int
        module procedure :: return_normal_float
        module procedure :: return_normal_double
        module procedure :: return_normal_complex
    end interface

    interface index_first_nonzero
        module procedure :: index_first_nonzero_int
        module procedure :: index_first_nonzero_float
        module procedure :: index_first_nonzero_double
        module procedure :: index_first_nonzero_logical
        module procedure :: index_first_nonzero_complex
    end interface

    interface equivalent
        module procedure :: equivalent_float
        module procedure :: equivalent_double
    end interface equivalent

    interface print_array
        module procedure :: print_array_2d_int
        module procedure :: print_array_2d_float
        module procedure :: print_array_2d_double
        module procedure :: print_array_2d_complex
        module procedure :: print_array_2d_dcomplex
    end interface print_array

    interface order_of_magnitude
        module procedure :: order_of_magnitude_int
        module procedure :: order_of_magnitude_float
        module procedure :: order_of_magnitude_double
    end interface

    interface rov
        module procedure :: rov_1d_int
        module procedure :: rov_1d_float
        module procedure :: rov_1d_double
        module procedure :: rov_2d_int
        module procedure :: rov_2d_float
        module procedure :: rov_2d_double
        module procedure :: rov_3d_int
        module procedure :: rov_3d_float
        module procedure :: rov_3d_double
        module procedure :: rov_4d_int
        module procedure :: rov_4d_float
        module procedure :: rov_4d_double
        module procedure :: rov_axis_2d_int
        module procedure :: rov_axis_2d_float
        module procedure :: rov_axis_2d_double
        module procedure :: rov_axis_3d_int
        module procedure :: rov_axis_3d_float
        module procedure :: rov_axis_3d_double
    end interface rov

    interface round
        module procedure :: round_float
        module procedure :: round_double
    end interface round

    interface nice
        module procedure :: nice_float
        module procedure :: nice_double
    end interface nice

    public :: warn
    public :: clip
    public :: swap
    public :: return_normal
    public :: index_first_nonzero
    public :: equivalent
    public :: print_array
    public :: get_hostname
    public :: get_integer_digits
    public :: order_of_magnitude
    public :: rov
    public :: round
    public :: nice

contains

    !
    !> Print array nicely
    !
    subroutine print_array_2d_int(a)

        integer, dimension(:, :), intent(in) :: a

        integer :: i

        do i = 1, size(a, 1)
            write(error_unit, *) a(i, :)
        end do

    end subroutine print_array_2d_int

    subroutine print_array_2d_float(a)

        real, dimension(:, :), intent(in) :: a

        integer :: i

        do i = 1, size(a, 1)
            write(error_unit, *) a(i, :)
        end do

    end subroutine print_array_2d_float

    subroutine print_array_2d_double(a)

        double precision, dimension(:, :), intent(in) :: a

        integer :: i

        do i = 1, size(a, 1)
            write(error_unit, *) a(i, :)
        end do

    end subroutine print_array_2d_double

    subroutine print_array_2d_complex(a)

        complex, dimension(:, :), intent(in) :: a

        integer :: i

        do i = 1, size(a, 1)
            write(error_unit, *) a(i, :)
        end do

    end subroutine print_array_2d_complex

    subroutine print_array_2d_dcomplex(a)

        double complex, dimension(:, :), intent(in) :: a

        integer :: i

        do i = 1, size(a, 1)
            write(error_unit, *) a(i, :)
        end do

    end subroutine print_array_2d_dcomplex

    elemental function return_normal_int(x) result(xn)

        integer, intent(in) :: x
        integer :: xn

        xn = x
        if (.not. ieee_is_normal(real(xn))) then
            xn = 0
        end if

    end function return_normal_int

    elemental function return_normal_float(x) result(xn)

        real, intent(in) :: x
        real :: xn

        xn = x
        if (.not. ieee_is_normal(xn)) then
            xn = 0.0
        end if

    end function return_normal_float

    elemental function return_normal_double(x) result(xn)

        double precision, intent(in) :: x
        double precision :: xn

        xn = x
        if (.not. ieee_is_normal(xn)) then
            xn = 0.0
        end if

    end function return_normal_double

    elemental function return_normal_complex(x) result(xn)

        complex, intent(in) :: x
        complex :: xn

        xn = x
        if (.not. ieee_is_normal(real(xn))) then
            xn = cmplx(0.0, imag(xn))
        end if
        if (.not. ieee_is_normal(imag(xn))) then
            xn = cmplx(real(xn), 0.0)
        end if

    end function return_normal_complex

    !
    !> Get the name of the host
    !
    function get_hostname() result(hostname)

        character(len=:), allocatable :: hostname
        character(len=128) :: hst

        call hostnm(hst)
        allocate (character(len=len_trim(hst)) :: hostname)
        hostname = trim(adjustl(hst))

    end function get_hostname

    !
    !> Show warn message
    !
    subroutine warn_string(message)
        character(len=*), intent(in) :: message
        write (error_unit, '(a)') message
    end subroutine warn_string

    subroutine warn_int2(message)
        integer(2), intent(in) :: message
        write (error_unit, '(a)') num2str(message)
    end subroutine warn_int2

    subroutine warn_int4(message)
        integer(4), intent(in) :: message
        write (error_unit, '(a)') num2str(message)
    end subroutine warn_int4

    subroutine warn_float(message)
        real, intent(in) :: message
        write (error_unit, '(a)') num2str(message)
    end subroutine warn_float

    subroutine warn_cmplx(message)
        complex, intent(in) :: message
        write (error_unit, '(a)') num2str(message)
    end subroutine warn_cmplx

    subroutine warn_double(message)
        double precision, intent(in) :: message
        write (error_unit, '(a)') num2str(message)
    end subroutine warn_double

    subroutine warn_logical(message)
        logical, intent(in) :: message
        write (error_unit, '(a)') num2str(message)
    end subroutine warn_logical

    !===========================================
    ! Index of the first non-zero element

    function index_first_nonzero_int(w) result(i)

        integer, dimension(:) :: w

        integer :: i

        do i = 1, size(w)
            if (w(i) /= 0) then
                exit
            end if
        end do

    end function index_first_nonzero_int

    function index_first_nonzero_float(w) result(i)

        real, dimension(:) :: w

        integer :: i

        do i = 1, size(w)
            if (w(i) /= 0) then
                exit
            end if
        end do

    end function index_first_nonzero_float

    function index_first_nonzero_double(w) result(i)

        double precision, dimension(:) :: w

        integer :: i

        do i = 1, size(w)
            if (w(i) /= 0) then
                exit
            end if
        end do

    end function index_first_nonzero_double

    function index_first_nonzero_logical(w) result(i)

        logical, dimension(:) :: w

        integer :: i

        do i = 1, size(w)
            if (.not. w(i)) then
                exit
            end if
        end do

    end function index_first_nonzero_logical

    function index_first_nonzero_complex(w) result(i)

        complex, dimension(:) :: w

        integer :: i

        do i = 1, size(w)
            if (real(w(i)) /= 0 .or. imag(w(i)) /= 0) then
                exit
            end if
        end do

    end function index_first_nonzero_complex

    !
    !> Find the order of magnitude of a quantity
    !
    elemental function order_of_magnitude_int(x) result(om)

        integer, intent(in) :: x
        integer :: om

        if (x == 0) then
            om = 0
        else
            om = floor(log10(abs(x*1.0)))
        end if

    end function order_of_magnitude_int

    elemental function order_of_magnitude_float(x) result(om)

        real, intent(in) :: x
        integer :: om

        if (x == 0) then
            om = 0
        else
            om = floor(log10(abs(x)))
        end if

    end function order_of_magnitude_float

    elemental function order_of_magnitude_double(x) result(om)

        double precision, intent(in) :: x
        integer :: om

        if (x == 0) then
            om = 0
        else
            om = floor(log10(abs(x)))
        end if

    end function order_of_magnitude_double

    elemental function round_float(x, base, precision) result(r)

        real, intent(in) :: x, base
        integer, intent(in), optional :: precision
        real :: r

        integer :: prec

        if (present(precision)) then
            prec = precision
        else
            prec = 2
        end if

        r = base*nint(x/base)
        r = r*10.0d0**prec
        r = nint(r)
        r = r/10.0d0**prec

    end function round_float

    elemental function round_double(x, base, precision) result(r)

        double precision, intent(in) :: x, base
        integer, intent(in), optional :: precision
        double precision :: r

        integer :: prec

        if (present(precision)) then
            prec = precision
        else
            prec = 2
        end if

        r = base*nint(x/base)
        r = r*10.0d0**prec
        r = nint(r)
        r = r/10.0d0**prec

    end function round_double

    elemental function nice_float(x, base) result(xnice)

        real, intent(in) :: x
        real, intent(in), optional :: base
        real :: xnice

        integer :: nm
        real :: b

        if (present(base)) then
            b = base
        else
            b = 0.5
        end if

        nm = order_of_magnitude(x)
        xnice = x/10.0d0**nm
        xnice = round(xnice, b)
        xnice = xnice*10.0d0**nm

    end function nice_float

    elemental function nice_double(x, base) result(xnice)

        double precision, intent(in) :: x
        double precision, intent(in), optional :: base
        double precision :: xnice

        integer :: nm
        double precision :: b

        if (present(base)) then
            b = base
        else
            b = 0.5
        end if

        nm = order_of_magnitude(x)
        xnice = x/10.0d0**nm
        xnice = round(xnice, b)
        xnice = xnice*10.0d0**nm

    end function nice_double

    !
    !> Separate an integer into digits
    !
    subroutine get_integer_digits(w, nw, wd)

        integer, intent(in) :: w, nw
        integer, dimension(:), intent(inout) :: wd

        integer :: rem, i

        rem = w
        do i = 1, nw
            wd(nw - i + 1) = rem - (rem/10)*10  ! Take advantage of integer division
            rem = rem/10
        end do

    end subroutine get_integer_digits

    elemental function clip_int(w, wmin, wmax) result(v)

        integer, intent(in) :: w, wmin, wmax
        integer :: v

        if (w < wmin) then
            v = wmin
        else if (w > wmax) then
            v = wmax
        else
            v = w
        end if

    end function clip_int

    elemental function clip_float(w, wmin, wmax) result(v)

        real, intent(in) :: w, wmin, wmax
        real :: v

        if (w < wmin) then
            v = wmin
        else if (w > wmax) then
            v = wmax
        else
            v = w
        end if

    end function clip_float

    elemental function clip_double(w, wmin, wmax) result(v)

        double precision, intent(in) :: w, wmin, wmax
        double precision :: v

        if (w < wmin) then
            v = wmin
        else if (w > wmax) then
            v = wmax
        else
            v = w
        end if

    end function clip_double

    elemental function clip_complex(w, wmin, wmax) result(v)

        complex, intent(in) :: w, wmin, wmax
        complex :: v

        if (abs(w) < abs(wmin)) then
            v = wmin
        else if (abs(w) > abs(wmax)) then
            v = wmax
        else
            v = w
        end if

    end function clip_complex

    elemental function arithop(w, op) result(wa)

        real, intent(in) :: w
        character(len=*), intent(in) :: op
        real :: wa

        wa = w

        select case (op)
            case ('abs')
                wa = abs(w)
            case ('-abs')
                wa = -abs(w)
            case ('log10')
                wa = sign(1.0, w)*log10(abs(w))
            case ('log2')
                wa = sign(1.0, w)*log(abs(w))/log(2.0)
            case ('ln')
                wa = sign(1.0, w)*log(abs(w))
            case ('sign')
                wa = sign(1.0, w)
            case ('sign*abs')
                wa = sign(1.0, w)*abs(w)
            case ('>0')
                wa = max(w, tiny(0.0))
            case ('<0')
                wa = min(w, -tiny(0.0))
            case ('>=0')
                wa = max(w, 0.0)
            case ('<=0')
                wa = min(w, 0.0)
            case ('sin')
                wa = sin(w)
            case ('cos')
                wa = cos(w)
            case ('tan')
                wa = tan(w)
            case ('asin')
                wa = asin(w)
            case ('acos')
                wa = acos(w)
            case ('atan')
                wa = atan(w)
            case ('sinh')
                wa = sinh(w)
            case ('cosh')
                wa = cosh(w)
            case ('tanh')
                wa = tanh(w)
        end select

    end function arithop

    !=========================================================
    !
    !> @brife Swap two numbers
    !
#define T int
#define TT integer
#include "template_swap.f90"

#define T float
#define TT real
#include "template_swap.f90"

#define T double
#define TT double precision
#include "template_swap.f90"

#define T complex
#define TT complex
#include "template_swap.f90"

#define T dcomplex
#define TT double complex
#include "template_swap.f90"

#define T logical
#define TT logical
#include "template_swap.f90"

    !===============================================
    ! Check equivalence
    !
    elemental function equivalent_float(a, b) result(yn)

        real, intent(in) :: a, b
        logical :: yn

        yn = .false.
        if (abs(a - b) < float_tiny) then
            yn = .true.
        end if

    end function equivalent_float

    elemental function equivalent_double(a, b) result(yn)

        double precision, intent(in) :: a, b
        logical :: yn

        yn = .false.
        if (abs(a - b) < double_tiny) then
            yn = .true.
        end if

    end function equivalent_double

    !===============================================
    ! Range of values
    !
#define T int
#define TT integer
#include "template_rov.f90"

#define T float
#define TT real
#include "template_rov.f90"

#define T double
#define TT double precision
#include "template_rov.f90"

end module libflit_utility
