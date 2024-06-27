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


module libflit_linear_algebra

    use libflit_constants
    use libflit_sort
    use libflit_array
    use libflit_error
    use libflit_sort
    use libflit_random
    use libflit_utility
    use blas95
    use lapack95
    use libflit_array_operation

    implicit none

    interface trace
        module procedure :: trace_float
        module procedure :: trace_double
        module procedure :: trace_complex
        module procedure :: trace_dcomplex
    end interface trace

    interface det
        module procedure :: det_nxn_float
        module procedure :: det_nxn_double
        module procedure :: det_nxn_complex
        module procedure :: det_nxn_dcomplex
    end interface det

    interface inv
        module procedure :: inv_nxn_float
        module procedure :: inv_nxn_double
        module procedure :: inv_nxn_complex
        module procedure :: inv_nxn_dcomplex
    end interface

    interface solve
        module procedure :: solve_single_rhs_float
        module procedure :: solve_multiple_rhs_float
        module procedure :: solve_single_rhs_double
        module procedure :: solve_multiple_rhs_double
        module procedure :: solve_single_rhs_complex
        module procedure :: solve_multiple_rhs_complex
        module procedure :: solve_single_rhs_dcomplex
        module procedure :: solve_multiple_rhs_dcomplex
    end interface solve

    interface solve_band
        module procedure :: solve_band_single_rhs_float
        module procedure :: solve_band_single_rhs_double
        module procedure :: solve_band_single_rhs_complex
        module procedure :: solve_band_single_rhs_dcomplex
    end interface solve_band

    interface lsqsolve
        module procedure :: least_squares_solve_single_rhs_float
        module procedure :: least_squares_solve_multiple_rhs_float
        module procedure :: least_squares_solve_single_rhs_double
        module procedure :: least_squares_solve_multiple_rhs_double
        module procedure :: least_squares_solve_single_rhs_complex
        module procedure :: least_squares_solve_multiple_rhs_complex
        module procedure :: least_squares_solve_single_rhs_dcomplex
        module procedure :: least_squares_solve_multiple_rhs_dcomplex
    end interface lsqsolve

    interface svd
        module procedure :: svd_float
        module procedure :: svd_double
        module procedure :: svd_complex
        module procedure :: svd_dcomplex
    end interface svd

    interface eigen
        module procedure :: eigen_simple_symmetric_l_or_r_float
        module procedure :: eigen_simple_symmetric_l_or_r_double

        module procedure :: eigen_simple_l_or_r_float
        module procedure :: eigen_simple_l_or_r_double
        module procedure :: eigen_simple_l_or_r_complex
        module procedure :: eigen_simple_l_or_r_dcomplex

        module procedure :: eigen_general_l_or_r_float
        module procedure :: eigen_general_l_or_r_double
        module procedure :: eigen_general_l_or_r_complex
        module procedure :: eigen_general_l_or_r_dcomplex

        module procedure :: eigen_simple_l_and_r_float
        module procedure :: eigen_simple_l_and_r_double
        module procedure :: eigen_simple_l_and_r_complex
        module procedure :: eigen_simple_l_and_r_dcomplex

        module procedure :: eigen_general_l_and_r_float
        module procedure :: eigen_general_l_and_r_double
        module procedure :: eigen_general_l_and_r_complex
        module procedure :: eigen_general_l_and_r_dcomplex
    end interface eigen

    interface matx
        module procedure :: matx_mat_vec_float
        module procedure :: matx_mat_vec_double
        module procedure :: matx_mat_vec_complex
        module procedure :: matx_mat_vec_dcomplex
        module procedure :: matx_mat_mat_float
        module procedure :: matx_mat_mat_double
        module procedure :: matx_mat_mat_complex
        module procedure :: matx_mat_mat_dcomplex
    end interface matx

    interface diagx
        module procedure :: diagx_mat_mat_float
        module procedure :: diagx_mat_mat_double
        module procedure :: diagx_mat_mat_complex
        module procedure :: diagx_mat_mat_dcomplex
    end interface diagx

    interface xdiag
        module procedure :: xdiag_mat_mat_float
        module procedure :: xdiag_mat_mat_double
        module procedure :: xdiag_mat_mat_complex
        module procedure :: xdiag_mat_mat_dcomplex
    end interface xdiag

    interface toeplitz_matrix
        module procedure :: toeplitz_float
        module procedure :: toeplitz_double
        module procedure :: toeplitz_complex
        module procedure :: toeplitz_dcomplex
        module procedure :: circulant_float
        module procedure :: circulant_complex
        module procedure :: circulant_double
        module procedure :: circulant_dcomplex
    end interface toeplitz_matrix

    interface hankel_matrix
        module procedure :: hankel_float
        module procedure :: hankel_double
        module procedure :: hankel_complex
        module procedure :: hankel_dcomplex
    end interface hankel_matrix

    interface vandermonde_matrix
        module procedure :: vandermonde_float
        module procedure :: vandermonde_double
        module procedure :: vandermonde_complex
        module procedure :: vandermonde_dcomplex
    end interface vandermonde_matrix

    interface spectral_radius
        module procedure :: spectral_radius_float
        module procedure :: spectral_radius_double
    end interface spectral_radius

    interface eigen_symm2x2
        module procedure :: eigen_symm2x2_float
        module procedure :: eigen_symm2x2_double
    end interface eigen_symm2x2

    interface eigen_symm3x3
        module procedure :: eigen_symm3x3_float
        module procedure :: eigen_symm3x3_double
    end interface eigen_symm3x3

    private
    public :: trace
    public :: det
    public :: inv
    public :: solve
    public :: lsqsolve
    public :: svd
    public :: eigen
    public :: matx
    public :: diagx
    public :: xdiag
    public :: solve_band
    public :: spectral_radius
    public :: eigen_symm2x2
    public :: eigen_symm3x3
    public :: toeplitz_matrix
    public :: hankel_matrix
    public :: vandermonde_matrix

contains

    ! Some linear algebra functions
#define T float
#define TT real
#define TTT real
#include "template_linear_algebra.f90"

#define T double
#define TT double precision
#define TTT double precision
#include "template_linear_algebra.f90"

#define T complex
#define TT complex
#define TTT real
#include "template_linear_algebra.f90"

#define T dcomplex
#define TT double complex
#define TTT double precision
#include "template_linear_algebra.f90"

    ! Eigensystem solvers
#define T float
#define TT real
#include "template_eigen.f90"

#define T double
#define TT double precision
#include "template_eigen.f90"

    ! General eigensystem solvers
#define T float
#define TT real
#define TTT complex
#define nTT cmplx
#include "template_eigen_general.f90"

#define T double
#define TT double precision
#define TTT double complex
#define nTT dcmplx
#include "template_eigen_general.f90"

#define T complex
#define TT complex
#define TTT complex
#define nTT
#include "template_eigen_general.f90"

#define T dcomplex
#define TT double complex
#define TTT double complex
#define nTT
#include "template_eigen_general.f90"

end module libflit_linear_algebra
