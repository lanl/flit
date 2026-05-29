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

module libflit_dip

    use libflit_array
    use libflit_array_extension
    use libflit_array_operation
    use libflit_gaussfilt
    use libflit_utility
    use libflit_andffilt
    use libflit_random
    use libflit_transform
    use libflit_interp
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32

    implicit none

    interface local_dip
        module procedure :: gst_local_dip_2d_float
        module procedure :: gst_local_dip_3d_float
        module procedure :: gst_local_dip_2d_double
        module procedure :: gst_local_dip_3d_double
    end interface

    interface dip_filt
        module procedure :: dip_filt_2d_float
        module procedure :: dip_filt_3d_float
        module procedure :: dip_filt_2d_double
        module procedure :: dip_filt_3d_double
    end interface

    private
    public :: local_dip
    public :: dip_filt

contains

#define T float
#define TT real
#define TTT meta_array1_real
#define TF real
#define fp sp
#include "template_dipfilt.f90"

#define T double
#define TT double precision
#define TTT meta_array1_double
#define TF dble
#define fp dp
#include "template_dipfilt.f90"

end module
