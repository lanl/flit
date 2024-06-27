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


#define PASTE(X)            X
#define PASTE2(X)           PASTE(X)_
#define CONCATHELP(X, Y)    PASTE2(X)Y
#define CONCAT(X, Y)        CONCATHELP(X, Y)

#define window_function_      CONCAT(window_function, T)

function window_function_(x0, method, alpha) result(w)

    !> Distance to zero; normalized to [0, 1]
    TT, dimension(:) :: x0
    !> Window type
    character(len=*) :: method
    !> Coefficient that has different meanings for different windows
    TT, optional :: alpha
    !> Value of window function
    TT, allocatable, dimension(:) :: w

    TT, allocatable, dimension(:) :: x
    TT :: a0, a1, a2, a3

    x = x0
    w = ones_like(x)

    select case (method)

        case ('step')
            w = 0

        case ('linear')
            where (x <= 0.5d0)
                w = 2*x
            end where
            where (x >= 0.5d0)
                w = 1 - 2*(x - 0.5d0)
            end where

        case ('parzen')
            x = x - 0.5d0
            x = x/0.5d0
            where (abs(x) <= 1.0d0/2.0d0)
                w = 1.0d0 - 6.0d0*x**2*(1.0d0 - abs(x))
            end where
            where (abs(x) > 1.0d0/2.0d0)
                w = 2*(1 - abs(x))**3
            end where

        case ('welch')
            x = x - 0.5d0
            x = x/0.5d0
            w = 1 - x**2

        case ('sine')
            w = sin(const_pi*x)

        case ('power-sine')
            w = sin(const_pi*x)**alpha

        case ('hann')
            w = sin(const_pi*x)**2

        case ('hamming')
            a0 = 25.0d0/46.0d0
            w = a0 - (1 - a0)*cos(2*const_pi*x)

        case ('blackman')
            a0 = 0.42d0
            a1 = 0.5d0
            a2 = 0.08d0
            w = a0 - a1*cos(2*const_pi*x) + a2*cos(4*const_pi*x)

        case ('nuttall')
            a0 = 0.355768d0
            a1 = 0.487396d0
            a2 = 0.144232d0
            a3 = 0.012604d0
            w = a0 - a1*cos(2*const_pi*x) + a2*cos(4*const_pi*x) - a3*cos(6*const_pi*x)

        case ('blackman-nuttall')
            a0 = 0.3635819d0
            a1 = 0.4891775d0
            a2 = 0.1365995d0
            a3 = 0.0106411d0
            w = a0 - a1*cos(2*const_pi*x) + a2*cos(4*const_pi*x) - a3*cos(6*const_pi*x)

        case ('blackman-harris')
            a0 = 0.35875d0
            a1 = 0.48829d0
            a2 = 0.14128d0
            a3 = 0.01168d0
            w = a0 - a1*cos(2*const_pi*x) + a2*cos(4*const_pi*x) - a3*cos(6*const_pi*x)

        case ('kaiser')
            x = x - 0.5d0
            x = x/0.5d0
            w = bessel_i0(const_pi*alpha*sqrt(1.0d0 - x**2))/bessel_i0(const_pi*alpha)

        case ('gauss')
            x = x - 0.5d0
            x = x/0.5d0
            w = exp(-0.5d0*(x/alpha)**2)

    end select

    where (w < 0)
        w = 0
    end where
    where (w > 1)
        w = 1
    end where

end function window_function_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef window_function_
