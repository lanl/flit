!
! © 2024-2026. Triad National Security, LLC. All rights reserved.
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


module libflit_gaussfilt

    use libflit_array
    use libflit_constants
    use libflit_transform
    use libflit_array_operation
    use libflit_error
    use libflit_linear_algebra
    use libflit_statistics
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32

    implicit none

    interface rotation_matrix
        module procedure :: rotation_matrix_2d_float
        module procedure :: rotation_matrix_2d_double
        module procedure :: rotation_matrix_3d_euler_float
        module procedure :: rotation_matrix_3d_euler_double
    end interface

    interface gauss_filt
        module procedure :: gauss_filt_1d_float
        module procedure :: gauss_filt_2d_float
        module procedure :: gauss_filt_3d_float
        module procedure :: gauss_filt_1d_double
        module procedure :: gauss_filt_2d_double
        module procedure :: gauss_filt_3d_double
    end interface

    private
    public :: gauss_filt

contains

#define T float
#define TT real
#define fp sp
#include "template_rotate.f90"
#undef fp

#define T double
#define TT double precision
#define fp dp
#include "template_rotate.f90"
#undef fp

#define T float
#define TT real
#include "template_gaussfilt.f90"

#define T double
#define TT double precision
#include "template_gaussfilt.f90"

end module libflit_gaussfilt
