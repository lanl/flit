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


module libflit_unique

    use libflit_array
    use ieee_arithmetic

    implicit none

    private
    interface unique
        module procedure :: unique_1d_int
        module procedure :: unique_1d_float
        module procedure :: unique_1d_double
        module procedure :: unique_1d_complex
        module procedure :: unique_1d_dcomplex
        module procedure :: unique_1d_string
        module procedure :: unique_2d_int
        module procedure :: unique_2d_float
        module procedure :: unique_2d_double
        module procedure :: unique_2d_complex
        module procedure :: unique_2d_dcomplex
        module procedure :: unique_2d_string
    end interface unique

    public :: unique

contains

#define T int
#define TT integer
#define TTT integer
#include "template_unique.f90"

#define T float
#define TT real
#define TTT real
#include "template_unique.f90"

#define T double
#define TT double precision
#define TTT double precision
#include "template_unique.f90"

#define T complex
#define TT complex
#define TTT complex
#include "template_unique.f90"

#define T dcomplex
#define TT double complex
#define TTT double complex
#include "template_unique.f90"

#define T logical
#define TT logical
#define TTT logical
#include "template_unique.f90"

#define T string
#define TT character(len=*)
#define TTT character(len=1024)
#include "template_unique.f90"

end module libflit_unique
