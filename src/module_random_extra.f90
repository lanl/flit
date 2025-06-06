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


module libflit_random_extra

    use libflit_array
    use libflit_array_operation
    use libflit_random
    use libflit_calculus
    use libflit_interp

    implicit none

    ! The module is to generate random values with customized PDF

    interface random_pdf
        module procedure :: random_pdf_1d_float
        module procedure :: random_pdf_1d_double
        module procedure :: random_pdf_2d_float
        module procedure :: random_pdf_2d_double
        module procedure :: random_pdf_3d_float
        module procedure :: random_pdf_3d_double
    end interface

    private
    public :: random_pdf

contains

#define T float
#define TT real
#define _random_ random
#include "template_random_pdf.f90"

#define T double
#define TT double precision
#define _random_ drandom
#include "template_random_pdf.f90"

end module
