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


module libflit_statistics

    use libflit_array
    use libflit_sort
    use libflit_constants
    use libflit_linear_algebra
    use libflit_transform
    use libflit_error
    use iso_fortran_env
    use libflit_array_operation

    implicit none

    ! Cross- and auto-correlation
    interface xcorrd
        module procedure :: discrete_xcorr_1d_float
        module procedure :: discrete_xcorr_2d_float
        module procedure :: discrete_xcorr_3d_float
        module procedure :: discrete_xcorr_1d_double
        module procedure :: discrete_xcorr_2d_double
        module procedure :: discrete_xcorr_3d_double
    end interface xcorrd

    interface acorrd
        module procedure :: discrete_acorr_1d_float
        module procedure :: discrete_acorr_2d_float
        module procedure :: discrete_acorr_3d_float
        module procedure :: discrete_acorr_1d_double
        module procedure :: discrete_acorr_2d_double
        module procedure :: discrete_acorr_3d_double
    end interface acorrd

    interface xcorr
        module procedure :: xcorr_1d_float
        module procedure :: xcorr_2d_float
        module procedure :: xcorr_3d_float
        module procedure :: xcorr_1d_double
        module procedure :: xcorr_2d_double
        module procedure :: xcorr_3d_double
    end interface xcorr

    interface acorr
        module procedure :: acorr_1d_float
        module procedure :: acorr_2d_float
        module procedure :: acorr_3d_float
        module procedure :: acorr_1d_double
        module procedure :: acorr_2d_double
        module procedure :: acorr_3d_double
    end interface acorr

    interface lxcorr
        module procedure :: local_xcorr_1d_float
        module procedure :: local_xcorr_1d_double
    end interface lxcorr

    interface xcorr_coef
        module procedure :: xcorr_coef_1d_float
        module procedure :: xcorr_coef_2d_float
        module procedure :: xcorr_coef_3d_float
        module procedure :: xcorr_coef_1d_double
        module procedure :: xcorr_coef_2d_double
        module procedure :: xcorr_coef_3d_double
    end interface xcorr_coef

    ! Mean
    interface mean
        module procedure :: mean_1d_float
        module procedure :: mean_2d_float
        module procedure :: mean_3d_float
        module procedure :: mean_1d_double
        module procedure :: mean_2d_double
        module procedure :: mean_3d_double
        module procedure :: mean_1d_complex
        module procedure :: mean_2d_complex
        module procedure :: mean_3d_complex
        module procedure :: mean_1d_dcomplex
        module procedure :: mean_2d_dcomplex
        module procedure :: mean_3d_dcomplex
    end interface mean

    ! Median
    interface median
        module procedure :: median_1d_float
        module procedure :: median_2d_float
        module procedure :: median_3d_float
        module procedure :: median_1d_double
        module procedure :: median_2d_double
        module procedure :: median_3d_double
    end interface median

    ! Standard deviation
    interface std
        module procedure :: standard_deviation_1d_float
        module procedure :: standard_deviation_2d_float
        module procedure :: standard_deviation_3d_float
        module procedure :: standard_deviation_1d_double
        module procedure :: standard_deviation_2d_double
        module procedure :: standard_deviation_3d_double
    end interface std

    ! Covariance
    interface covar
        module procedure :: covariance_1d_float
        module procedure :: covariance_2d_float
        module procedure :: covariance_3d_float
        module procedure :: covariance_1d_double
        module procedure :: covariance_2d_double
        module procedure :: covariance_3d_double
    end interface covar

    ! Compute histogram
    interface histogram
        module procedure :: histogram_1d_float
        module procedure :: histogram_2d_float
        module procedure :: histogram_3d_float
        module procedure :: histogram_1d_double
        module procedure :: histogram_2d_double
        module procedure :: histogram_3d_double
    end interface histogram

    ! Plot histogram
    interface plot_histogram
        module procedure :: plot_histogram_1d_float
        module procedure :: plot_histogram_2d_float
        module procedure :: plot_histogram_3d_float
        module procedure :: plot_histogram_1d_double
        module procedure :: plot_histogram_2d_double
        module procedure :: plot_histogram_3d_double
    end interface plot_histogram

    ! Gaussian pdf
    interface gaussian
        module procedure :: gaussian_1d_float
        module procedure :: gaussian_2d_float
        module procedure :: gaussian_3d_float
        module procedure :: gaussian_1d_double
        module procedure :: gaussian_2d_double
        module procedure :: gaussian_3d_double
    end interface gaussian

    ! Covariance matrix
    interface covar_matrix
        module procedure :: covariance_matrix_1d_float
        module procedure :: covariance_matrix_2d_float
        module procedure :: covariance_matrix_3d_float
        module procedure :: cross_covariance_matrix_1d_float
        module procedure :: cross_covariance_matrix_2d_float
        module procedure :: cross_covariance_matrix_3d_float
        module procedure :: row_covariance_matrix_2d_float
        module procedure :: row_cross_covariance_matrix_2d_float
        module procedure :: covariance_matrix_1d_double
        module procedure :: covariance_matrix_2d_double
        module procedure :: covariance_matrix_3d_double
        module procedure :: cross_covariance_matrix_1d_double
        module procedure :: cross_covariance_matrix_2d_double
        module procedure :: cross_covariance_matrix_3d_double
        module procedure :: row_covariance_matrix_2d_double
        module procedure :: row_cross_covariance_matrix_2d_double
    end interface covar_matrix

    ! Variance
    interface var
        module procedure :: variance_1d_float
        module procedure :: variance_2d_float
        module procedure :: variance_3d_float
        module procedure :: variance_1d_double
        module procedure :: variance_2d_double
        module procedure :: variance_3d_double
    end interface var

    ! Kernels
    interface kernel_triangular
        module procedure :: kernel_triangular_float
        module procedure :: kernel_triangular_double
    end interface kernel_triangular

    interface kernel_bisquare
        module procedure :: kernel_bisquare_float
        module procedure :: kernel_bisquare_double
    end interface kernel_bisquare

    interface kernel_trisquare
        module procedure :: kernel_trisquare_float
        module procedure :: kernel_trisquare_double
    end interface kernel_trisquare

    interface kernel_tricube
        module procedure :: kernel_tricube_float
        module procedure :: kernel_tricube_double
    end interface kernel_tricube

    interface kernel_epanechnikov
        module procedure :: kernel_epanechnikov_float
        module procedure :: kernel_epanechnikov_double
    end interface kernel_epanechnikov

    interface kernel_cosine
        module procedure :: kernel_cosine_float
        module procedure :: kernel_cosine_double
    end interface kernel_cosine

    private
    public :: xcorrd, acorrd
    public :: xcorr, acorr
    public :: mean
    public :: median
    public :: var, covar, covar_matrix
    public :: std
    public :: histogram, plot_histogram
    public :: lxcorr
    public :: xcorr_coef
    public :: gaussian
    public :: kernel_triangular
    public :: kernel_bisquare
    public :: kernel_trisquare
    public :: kernel_tricube
    public :: kernel_epanechnikov
    public :: kernel_cosine

contains

#define T float
#define TT real
#define TTT cmplx
#include "template_statistics.f90"

#define T double
#define TT double precision
#define TTT dcmplx
#include "template_statistics.f90"

    ! Some can be extended to complex domain
#define T float
#define TT real
#define nTT dble
#define mTT real
#include "template_mean.f90"

#define T double
#define TT double precision
#define nTT
#define mTT
#include "template_mean.f90"

#define T complex
#define TT complex
#define nTT dcmplx
#define mTT cmplx
#include "template_mean.f90"

#define T dcomplex
#define TT double complex
#define nTT
#define mTT
#include "template_mean.f90"

end module libflit_statistics
