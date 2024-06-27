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


module libflit_sort

    implicit none

    interface sort
        module procedure :: qksort_int
        module procedure :: qksort_float
        module procedure :: qksort_double
        module procedure :: qksort_2d_int
        module procedure :: qksort_2d_float
        module procedure :: qksort_2d_double
    end interface sort

    interface sort_index
        module procedure :: qksort_index_int
        module procedure :: qksort_index_float
        module procedure :: qksort_index_double
    end interface sort_index

    interface mergesort_index
        module procedure :: merge_sort_index_int
        module procedure :: merge_sort_index_float
        module procedure :: merge_sort_index_double
    end interface mergesort_index

    private
    public :: sort
    public :: sort_index
    public :: mergesort_index

contains

#define T int
#define TT integer
#include "template_sort.f90"

#define T float
#define TT real
#include "template_sort.f90"

#define T double
#define TT double precision
#include "template_sort.f90"

end module libflit_sort
