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


#define T int
#define TT integer
#include "template_meta_array.f90"

#define T real
#define TT real
#include "template_meta_array.f90"

#define T double
#define TT double precision
#include "template_meta_array.f90"

#define T complex
#define TT complex
#include "template_meta_array.f90"

#define T dcomplex
#define TT double complex
#include "template_meta_array.f90"

#define T logical
#define TT logical
#include "template_meta_array.f90"

module libflit_array_extension

    use libflit_meta_array_int
    use libflit_meta_array_real
    use libflit_meta_array_double
    use libflit_meta_array_complex
    use libflit_meta_array_dcomplex
    use libflit_meta_array_logical

    implicit none

    ! Get a meta-array from an array of meta-arrays by name; the return is a meta-array
    interface get_meta_array
        module procedure :: get_meta_array1_int
        module procedure :: get_meta_array1_real
        module procedure :: get_meta_array1_double
        module procedure :: get_meta_array1_complex
        module procedure :: get_meta_array1_dcomplex
        module procedure :: get_meta_array1_logical
        module procedure :: get_meta_array2_int
        module procedure :: get_meta_array2_real
        module procedure :: get_meta_array2_double
        module procedure :: get_meta_array2_complex
        module procedure :: get_meta_array2_dcomplex
        module procedure :: get_meta_array2_logical
        module procedure :: get_meta_array3_int
        module procedure :: get_meta_array3_real
        module procedure :: get_meta_array3_double
        module procedure :: get_meta_array3_complex
        module procedure :: get_meta_array3_dcomplex
        module procedure :: get_meta_array3_logical
        module procedure :: get_meta_array4_int
        module procedure :: get_meta_array4_real
        module procedure :: get_meta_array4_double
        module procedure :: get_meta_array4_complex
        module procedure :: get_meta_array4_dcomplex
        module procedure :: get_meta_array4_logical
    end interface

    ! Get a meta-array core from an array of meta-arrays by name; the return is the core array of a meta-array
    interface get_meta_array_core
        module procedure :: get_meta_array_core1_int
        module procedure :: get_meta_array_core1_real
        module procedure :: get_meta_array_core1_double
        module procedure :: get_meta_array_core1_complex
        module procedure :: get_meta_array_core1_dcomplex
        module procedure :: get_meta_array_core1_logical
        module procedure :: get_meta_array_core2_int
        module procedure :: get_meta_array_core2_real
        module procedure :: get_meta_array_core2_double
        module procedure :: get_meta_array_core2_complex
        module procedure :: get_meta_array_core2_dcomplex
        module procedure :: get_meta_array_core2_logical
        module procedure :: get_meta_array_core3_int
        module procedure :: get_meta_array_core3_real
        module procedure :: get_meta_array_core3_double
        module procedure :: get_meta_array_core3_complex
        module procedure :: get_meta_array_core3_dcomplex
        module procedure :: get_meta_array_core3_logical
        module procedure :: get_meta_array_core4_int
        module procedure :: get_meta_array_core4_real
        module procedure :: get_meta_array_core4_double
        module procedure :: get_meta_array_core4_complex
        module procedure :: get_meta_array_core4_dcomplex
        module procedure :: get_meta_array_core4_logical
    end interface

    ! Assign an array to some meta-array in an array of meta-arrays by name
    interface set_meta_array_core
        module procedure :: set_meta_array_core1_int
        module procedure :: set_meta_array_core1_real
        module procedure :: set_meta_array_core1_double
        module procedure :: set_meta_array_core1_complex
        module procedure :: set_meta_array_core1_dcomplex
        module procedure :: set_meta_array_core1_logical
        module procedure :: set_meta_array_core2_int
        module procedure :: set_meta_array_core2_real
        module procedure :: set_meta_array_core2_double
        module procedure :: set_meta_array_core2_complex
        module procedure :: set_meta_array_core2_dcomplex
        module procedure :: set_meta_array_core2_logical
        module procedure :: set_meta_array_core3_int
        module procedure :: set_meta_array_core3_real
        module procedure :: set_meta_array_core3_double
        module procedure :: set_meta_array_core3_complex
        module procedure :: set_meta_array_core3_dcomplex
        module procedure :: set_meta_array_core3_logical
        module procedure :: set_meta_array_core4_int
        module procedure :: set_meta_array_core4_real
        module procedure :: set_meta_array_core4_double
        module procedure :: set_meta_array_core4_complex
        module procedure :: set_meta_array_core4_dcomplex
        module procedure :: set_meta_array_core4_logical
    end interface

end module libflit_array_extension


