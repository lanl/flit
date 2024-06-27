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

!
!> The module implements LOWESS -- locally weighted scatterplot smoothing
!>        LOWESS is in fact a curve-fitting method, but could be considered equivalently as
!>        a smoothing/filtering technique.
!
module libflit_lowessfilt

    use libflit_array
    use libflit_error
    use libflit_array_operation
    use libflit_linear_algebra
    use libflit_utility
    use libflit_statistics

    implicit none

    interface lowess_filt
        module procedure :: lowess_1d_float
        module procedure :: lowess_2d_float
        module procedure :: lowess_3d_float
        module procedure :: lowess_1d_double
        module procedure :: lowess_2d_double
        module procedure :: lowess_3d_double
    end interface lowess_filt

    private
    public :: lowess_filt

contains

#define T float
#define TT real
#include "template_lowessfilt.f90"

#define T double
#define TT double precision
#include "template_lowessfilt.f90"

end module libflit_lowessfilt

