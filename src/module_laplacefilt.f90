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


module libflit_laplacefilt

    use libflit_array
    use libflit_error
    use libflit_array_operation

    implicit none

    interface laplace_filt
        module procedure :: laplace_filt_1d_float
        module procedure :: laplace_filt_2d_float
        module procedure :: laplace_filt_3d_float
        module procedure :: laplace_filt_1d_double
        module procedure :: laplace_filt_2d_double
        module procedure :: laplace_filt_3d_double
        module procedure :: laplace_filt_1d_complex
        module procedure :: laplace_filt_2d_complex
        module procedure :: laplace_filt_3d_complex
        module procedure :: laplace_filt_1d_dcomplex
        module procedure :: laplace_filt_2d_dcomplex
        module procedure :: laplace_filt_3d_dcomplex
    end interface

    private
    public :: laplace_filt

contains

#define T float
#define TT real
#include "template_laplacefilt.f90"

#define T double
#define TT double precision
#include "template_laplacefilt.f90"

#define T complex
#define TT complex
#include "template_laplacefilt.f90"

#define T dcomplex
#define TT double complex
#include "template_laplacefilt.f90"

end module libflit_laplacefilt
