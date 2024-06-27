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


module libflit_fourierfilt

    use libflit_constants
    use libflit_transform
    use libflit_array
    use libflit_error
    use libflit_taper
    use libflit_array_operation

    implicit none

    private

    ! 1D Fourier filter defined by freqs and amps
    interface fourier_filter
        module procedure :: fourier_filter_1d_float
        module procedure :: fourier_filter_1d_double
    end interface fourier_filter

    ! Dimension-separable Fourier filtering defined by freqs and amps
    interface fourier_filt
        module procedure :: fourier_filt_1d_float
        module procedure :: fourier_filt_2d_float
        module procedure :: fourier_filt_3d_float
        module procedure :: fourier_filt_1d_double
        module procedure :: fourier_filt_2d_double
        module procedure :: fourier_filt_3d_double
    end interface fourier_filt

    ! Fourier high-frequency emphasis filtering
    interface fourier_sharpen
        module procedure :: fourier_sharpen_1d_float
        module procedure :: fourier_sharpen_2d_float
        module procedure :: fourier_sharpen_3d_float
        module procedure :: fourier_sharpen_1d_double
        module procedure :: fourier_sharpen_2d_double
        module procedure :: fourier_sharpen_3d_double
    end interface fourier_sharpen

    ! Fourier low-frequency emphasis filtering
    interface fourier_smooth
        module procedure :: fourier_smooth_1d_float
        module procedure :: fourier_smooth_2d_float
        module procedure :: fourier_smooth_3d_float
        module procedure :: fourier_smooth_1d_double
        module procedure :: fourier_smooth_2d_double
        module procedure :: fourier_smooth_3d_double
    end interface fourier_smooth

    public :: fourier_filter
    public :: fourier_filt
    public :: fourier_sharpen
    public :: fourier_smooth

contains

#define T float
#define TT real
#define TTT real
#include "template_fourierfilt.f90"

#define T double
#define TT double precision
#define TTT dble
#include "template_fourierfilt.f90"

end module libflit_fourierfilt
