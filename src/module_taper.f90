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


module libflit_taper

    use libflit_constants
    use libflit_array
    use libflit_specialfunc
    use libflit_utility
    use libflit_error

    implicit none

    private

    interface window_function
        module procedure :: window_function_float
        module procedure :: window_function_double
    end interface window_function

    interface taper
        module procedure :: taper_1d_float
        module procedure :: taper_2d_float
        module procedure :: taper_3d_float
        module procedure :: taper_1d_double
        module procedure :: taper_2d_double
        module procedure :: taper_3d_double
        module procedure :: taper_1d_complex
        module procedure :: taper_2d_complex
        module procedure :: taper_3d_complex
        module procedure :: taper_1d_dcomplex
        module procedure :: taper_2d_dcomplex
        module procedure :: taper_3d_dcomplex
    end interface taper

    public :: taper_window
    public :: taper
    public :: window_function

contains

    !
    !> Various types of taper window functions
    !
    function left_taper(len, method, alpha) result(taper)

        integer :: len
        character(len=*), optional :: method
        real, optional :: alpha
        real, dimension(:), allocatable :: taper

        integer :: n, i
        double precision :: a0, a1, a2, a3, a
        character(len=24) :: taper_method

        if (present(alpha)) then
            a = alpha
        else
            a = 0.0
        end if

        if (present(method)) then
            taper_method = method
        else
            taper_method = 'hann'
        end if

        taper = ones(len)
        n = 2*len - 1

        do i = 0, len - 1

            select case (taper_method)

                case ('step')
                    taper(i + 1) = 0.0

                case ('linear')
                    taper(i + 1) = i*1.0/(len - 1.0)

                case ('parzen')
                    if (i <= n/4.0) then
                        taper(i + 1) = 6.0*(i/(n/2.0))**2*(1.0 - i/(n/2.0))
                    else
                        taper(i + 1) = 1.0 - 2.0*(1.0 - i/(n/2.0))**3
                    end if

                case ('welch')
                    taper(i + 1) = 1.0 - ((i - (n - 1.0)/2.0)/((n - 1.0)/2.0))**2

                case ('sine')
                    taper(i + 1) = sin(const_pi*i/(n - 1.0))

                case ('power-sine')
                    taper(i + 1) = sin(const_pi*i/(n - 1.0))**a

                case ('hann')
                    taper(i + 1) = sin(const_pi*i/(n - 1.0))**2

                case ('hamming')
                    a0 = 25.0d0/46.0d0
                    taper(i + 1) = a0 - (1.0 - a0)*cos(2.0*const_pi*i/(n - 1.0))

                case ('blackman')
                    a0 = 0.42d0
                    a1 = 0.5d0
                    a2 = 0.08d0
                    taper(i + 1) = a0 - a1*cos(2*const_pi*i/(n - 1.0)) + a2*cos(4*const_pi*i/(n - 1.0))

                case ('nuttall')
                    a0 = 0.355768d0
                    a1 = 0.487396d0
                    a2 = 0.144232d0
                    a3 = 0.012604d0
                    taper(i + 1) = a0 - a1*cos(2*const_pi*i/(n - 1.0)) + a2*cos(4*const_pi*i/(n - 1.0)) &
                        - a3*cos(6*const_pi*i/(n - 1.0))

                case ('blackman-nuttall')
                    a0 = 0.3635819d0
                    a1 = 0.4891775d0
                    a2 = 0.1365995d0
                    a3 = 0.0106411d0
                    taper(i + 1) = a0 - a1*cos(2*const_pi*i/(n - 1.0)) + a2*cos(4*const_pi*i/(n - 1.0)) &
                        - a3*cos(6*const_pi*i/(n - 1.0))

                case ('blackman-harris')
                    a0 = 0.35875d0
                    a1 = 0.48829d0
                    a2 = 0.14128d0
                    a3 = 0.01168d0
                    taper(i + 1) = a0 - a1*cos(2*const_pi*i/(n - 1.0)) + a2*cos(4*const_pi*i/(n - 1.0)) &
                        - a3*cos(6*const_pi*i/(n - 1.0))

                case ('kaiser')
                    taper(i + 1) = bessel_i0(dble(const_pi*a*sqrt(1.0 - (2.0*i/(n - 1.0) - 1.0)**2)))/bessel_i0(dble(a*const_pi))

                case ('gauss')
                    taper(i + 1) = exp(-0.5*((i - (n - 1.0)/2.0)/(a*(n - 1.0)/2.0))**2)

            end select

        end do

        select case(method)
            case('step', 'linear', 'welch', 'sine', 'power-sine', 'hann', 'blackman', 'nuttall')
                taper(1) = 0.0
        end select

    end function left_taper

    ! Window function
#define T float
#define TT real
#include "template_window.f90"

#define T double
#define TT double precision
#include "template_window.f90"

    !
    !> Various types of taper window functions
    !
    function taper_window(ns, len, method, alpha) result(taper)

        integer :: ns
        integer, dimension(:), optional :: len
        character(len=*), dimension(:), optional :: method
        real, dimension(:), optional :: alpha
        real, dimension(:), allocatable :: taper

        integer :: n1, n2
        real, dimension(1:2) :: taper_alpha
        character(len=24), dimension(1:2) :: taper_method

        if (present(len)) then
            n1 = len(1)
            n2 = len(2)
        else
            n1 = min(ns, max(5, nint(0.05*ns)))
            n2 = min(ns, max(5, nint(0.05*ns)))
        end if

        call assert(all(len >= 0), ' <taper_window> Error: Taper length must >= 0.')
        call assert(ns >= n1 + n2 - 1, ' <taper_window> Error: Number of samples is too small for creating the taper.')

        if (present(alpha)) then
            taper_alpha = alpha
        else
            taper_alpha = 0.0
        end if

        if (present(method)) then
            taper_method = method
        else
            taper_method = ['hann', 'hann']
        end if

        taper = ones(ns)

        ! The following two approaches are essentially the same
        if (n1 >= 1) then
            taper(1:n1) =  window_function(linspace(0.0, 0.5, n1), taper_method(1), taper_alpha(1))
        end if
        if (n2 >= 1) then
            taper(ns - n2 + 1:ns) = window_function(linspace(0.5, 1.0, n2), taper_method(2), taper_alpha(2))
        end if

        !        if (n1 >= 1) then
        !            taper(1:n1) = left_taper(n1, taper_method(1), taper_alpha(1))
        !        end if
        !        if (n2 >= 1) then
        !            taper(ns:ns - n2 + 1:-1) = left_taper(n2, taper_method(2), taper_alpha(2))
        !        end if

    end function taper_window

    ! Tapering array
#define T float
#define TT real
#include "template_taper.f90"

#define T double
#define TT double precision
#include "template_taper.f90"

#define T complex
#define TT complex
#include "template_taper.f90"

#define T dcomplex
#define TT double complex
#include "template_taper.f90"

end module libflit_taper
