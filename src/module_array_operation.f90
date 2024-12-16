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


module libflit_array_operation

    use libflit_error
    use libflit_array
    use libflit_constants

    implicit none

    ! Adjust the shape of an array
    interface adjust
        module procedure :: adjust_1d_int
        module procedure :: adjust_2d_int
        module procedure :: adjust_3d_int
        module procedure :: adjust_4d_int
        module procedure :: adjust_1d_float
        module procedure :: adjust_2d_float
        module procedure :: adjust_3d_float
        module procedure :: adjust_4d_float
        module procedure :: adjust_1d_double
        module procedure :: adjust_2d_double
        module procedure :: adjust_3d_double
        module procedure :: adjust_4d_double
        module procedure :: adjust_1d_complex
        module procedure :: adjust_2d_complex
        module procedure :: adjust_3d_complex
        module procedure :: adjust_4d_complex
        module procedure :: adjust_1d_dcomplex
        module procedure :: adjust_2d_dcomplex
        module procedure :: adjust_3d_dcomplex
        module procedure :: adjust_4d_dcomplex
        module procedure :: adjust_1d_logical
        module procedure :: adjust_2d_logical
        module procedure :: adjust_3d_logical
        module procedure :: adjust_4d_logical
    end interface adjust

    ! Convert a n-dimensional 1-sized array to a scalar
    interface as_scalar
        module procedure :: as_scalar_1d_int
        module procedure :: as_scalar_2d_int
        module procedure :: as_scalar_3d_int
        module procedure :: as_scalar_4d_int
        module procedure :: as_scalar_1d_float
        module procedure :: as_scalar_2d_float
        module procedure :: as_scalar_3d_float
        module procedure :: as_scalar_4d_float
        module procedure :: as_scalar_1d_double
        module procedure :: as_scalar_2d_double
        module procedure :: as_scalar_3d_double
        module procedure :: as_scalar_4d_double
        module procedure :: as_scalar_1d_complex
        module procedure :: as_scalar_2d_complex
        module procedure :: as_scalar_3d_complex
        module procedure :: as_scalar_4d_complex
        module procedure :: as_scalar_1d_dcomplex
        module procedure :: as_scalar_2d_dcomplex
        module procedure :: as_scalar_3d_dcomplex
        module procedure :: as_scalar_4d_dcomplex
        module procedure :: as_scalar_1d_logical
        module procedure :: as_scalar_2d_logical
        module procedure :: as_scalar_3d_logical
        module procedure :: as_scalar_4d_logical
    end interface as_scalar

    ! Crop an array by range in-place
    ! After cropping, the array's lower and upper bounds may change.
    ! This is in constrast to functions crop, where the lower bound of
    ! the new array is always 1.
    interface crop_array
        module procedure :: crop_array_1d_int
        module procedure :: crop_array_2d_int
        module procedure :: crop_array_3d_int
        module procedure :: crop_array_4d_int
        module procedure :: crop_array_1d_float
        module procedure :: crop_array_2d_float
        module procedure :: crop_array_3d_float
        module procedure :: crop_array_4d_float
        module procedure :: crop_array_1d_double
        module procedure :: crop_array_2d_double
        module procedure :: crop_array_3d_double
        module procedure :: crop_array_4d_double
        module procedure :: crop_array_1d_complex
        module procedure :: crop_array_2d_complex
        module procedure :: crop_array_3d_complex
        module procedure :: crop_array_4d_complex
        module procedure :: crop_array_1d_logical
        module procedure :: crop_array_2d_logical
        module procedure :: crop_array_3d_logical
        module procedure :: crop_array_4d_logical
    end interface crop_array

    ! Crop an array by range
    interface crop
        module procedure :: crop_1d_int
        module procedure :: crop_2d_int
        module procedure :: crop_3d_int
        module procedure :: crop_4d_int
        module procedure :: crop_1d_float
        module procedure :: crop_2d_float
        module procedure :: crop_3d_float
        module procedure :: crop_4d_float
        module procedure :: crop_1d_double
        module procedure :: crop_2d_double
        module procedure :: crop_3d_double
        module procedure :: crop_4d_double
        module procedure :: crop_1d_complex
        module procedure :: crop_2d_complex
        module procedure :: crop_3d_complex
        module procedure :: crop_4d_complex
        module procedure :: crop_1d_logical
        module procedure :: crop_2d_logical
        module procedure :: crop_3d_logical
        module procedure :: crop_4d_logical
    end interface crop

    ! Cross product of 1D vectors
    interface cross
        module procedure :: cross_product_1d_int
        module procedure :: cross_product_1d_float
        module procedure :: cross_product_1d_double
        module procedure :: cross_product_1d_complex
        module procedure :: cross_product_1d_dcomplex
    end interface cross

    ! Flatten an n-dimensional array to 1D
    interface flatten
        module procedure :: flatten_2d_int
        module procedure :: flatten_3d_int
        module procedure :: flatten_4d_int
        module procedure :: flatten_2d_float
        module procedure :: flatten_3d_float
        module procedure :: flatten_4d_float
        module procedure :: flatten_2d_double
        module procedure :: flatten_3d_double
        module procedure :: flatten_4d_double
        module procedure :: flatten_2d_complex
        module procedure :: flatten_3d_complex
        module procedure :: flatten_4d_complex
        module procedure :: flatten_2d_dcomplex
        module procedure :: flatten_3d_dcomplex
        module procedure :: flatten_4d_dcomplex
        module procedure :: flatten_2d_logical
        module procedure :: flatten_3d_logical
        module procedure :: flatten_4d_logical
    end interface flatten

    ! Flip array along one or several axes
    interface flip
        module procedure :: flip_1d_int
        module procedure :: flip_2d_int
        module procedure :: flip_3d_int
        module procedure :: flip_1d_float
        module procedure :: flip_2d_float
        module procedure :: flip_3d_float
        module procedure :: flip_1d_double
        module procedure :: flip_2d_double
        module procedure :: flip_3d_double
        module procedure :: flip_1d_complex
        module procedure :: flip_2d_complex
        module procedure :: flip_3d_complex
        module procedure :: flip_1d_dcomplex
        module procedure :: flip_2d_dcomplex
        module procedure :: flip_3d_dcomplex
        module procedure :: flip_1d_logical
        module procedure :: flip_2d_logical
        module procedure :: flip_3d_logical
    end interface flip

    ! If-then-else ternary operation, including scalar and array
    interface ifelse
        module procedure :: ifelse_int
        module procedure :: ifelse_float
        module procedure :: ifelse_double
        module procedure :: ifelse_complex
        module procedure :: ifelse_dcomplex
        module procedure :: ifelse_logical
        module procedure :: ifelse_1d_int
        module procedure :: ifelse_2d_int
        module procedure :: ifelse_3d_int
        module procedure :: ifelse_4d_int
        module procedure :: ifelse_1d_float
        module procedure :: ifelse_2d_float
        module procedure :: ifelse_3d_float
        module procedure :: ifelse_4d_float
        module procedure :: ifelse_1d_double
        module procedure :: ifelse_2d_double
        module procedure :: ifelse_3d_double
        module procedure :: ifelse_4d_double
        module procedure :: ifelse_1d_complex
        module procedure :: ifelse_2d_complex
        module procedure :: ifelse_3d_complex
        module procedure :: ifelse_4d_complex
        module procedure :: ifelse_1d_dcomplex
        module procedure :: ifelse_2d_dcomplex
        module procedure :: ifelse_3d_dcomplex
        module procedure :: ifelse_4d_dcomplex
        module procedure :: ifelse_1d_logical
        module procedure :: ifelse_2d_logical
        module procedure :: ifelse_3d_logical
        module procedure :: ifelse_4d_logical
        module procedure :: ifelse_string
    end interface ifelse

    ! Mask an array
    interface mask
        module procedure :: mask_1d_int
        module procedure :: mask_2d_int
        module procedure :: mask_3d_int
        module procedure :: mask_4d_int
        module procedure :: mask_1d_float
        module procedure :: mask_2d_float
        module procedure :: mask_3d_float
        module procedure :: mask_4d_float
        module procedure :: mask_1d_double
        module procedure :: mask_2d_double
        module procedure :: mask_3d_double
        module procedure :: mask_4d_double
        module procedure :: mask_1d_complex
        module procedure :: mask_2d_complex
        module procedure :: mask_3d_complex
        module procedure :: mask_4d_complex
        module procedure :: mask_1d_dcomplex
        module procedure :: mask_2d_dcomplex
        module procedure :: mask_3d_dcomplex
        module procedure :: mask_4d_dcomplex
    end interface mask

    ! Pad an array
    interface pad_array
        module procedure :: pad_array_1d_int
        module procedure :: pad_array_2d_int
        module procedure :: pad_array_3d_int
        module procedure :: pad_array_1d_float
        module procedure :: pad_array_2d_float
        module procedure :: pad_array_3d_float
        module procedure :: pad_array_1d_double
        module procedure :: pad_array_2d_double
        module procedure :: pad_array_3d_double
        module procedure :: pad_array_1d_complex
        module procedure :: pad_array_2d_complex
        module procedure :: pad_array_3d_complex
        module procedure :: pad_array_1d_dcomplex
        module procedure :: pad_array_2d_dcomplex
        module procedure :: pad_array_3d_dcomplex
        module procedure :: pad_array_1d_logical
        module procedure :: pad_array_2d_logical
        module procedure :: pad_array_3d_logical
        module procedure :: pad_array_1d_string
        module procedure :: pad_array_2d_string
        module procedure :: pad_array_3d_string
    end interface pad_array

    ! Pad an array
    interface pad
        module procedure :: pad_1d_int
        module procedure :: pad_2d_int
        module procedure :: pad_3d_int
        module procedure :: pad_1d_float
        module procedure :: pad_2d_float
        module procedure :: pad_3d_float
        module procedure :: pad_1d_double
        module procedure :: pad_2d_double
        module procedure :: pad_3d_double
        module procedure :: pad_1d_complex
        module procedure :: pad_2d_complex
        module procedure :: pad_3d_complex
        module procedure :: pad_1d_dcomplex
        module procedure :: pad_2d_dcomplex
        module procedure :: pad_3d_dcomplex
        module procedure :: pad_1d_logical
        module procedure :: pad_2d_logical
        module procedure :: pad_3d_logical
        module procedure :: pad_1d_string
        module procedure :: pad_2d_string
        module procedure :: pad_3d_string
    end interface pad

    ! Permute an array in-place
    interface permute_array
        module procedure :: permute_array_3d_int
        module procedure :: permute_array_3d_float
        module procedure :: permute_array_3d_double
        module procedure :: permute_array_3d_complex
        module procedure :: permute_array_3d_dcomplex
        module procedure :: permute_array_3d_logical
        module procedure :: permute_array_4d_int
        module procedure :: permute_array_4d_float
        module procedure :: permute_array_4d_double
        module procedure :: permute_array_4d_complex
        module procedure :: permute_array_4d_dcomplex
        module procedure :: permute_array_4d_logical
    end interface permute_array

    ! Permute an array out-of-place
    interface permute
        module procedure :: permute_3d_int
        module procedure :: permute_3d_float
        module procedure :: permute_3d_double
        module procedure :: permute_3d_complex
        module procedure :: permute_3d_dcomplex
        module procedure :: permute_3d_logical
        module procedure :: permute_4d_int
        module procedure :: permute_4d_float
        module procedure :: permute_4d_double
        module procedure :: permute_4d_complex
        module procedure :: permute_4d_dcomplex
        module procedure :: permute_4d_logical
    end interface permute

    ! Linearly scale array to a specific range
    interface rescale
        module procedure :: rescale_1d_float
        module procedure :: rescale_2d_float
        module procedure :: rescale_3d_float
        module procedure :: rescale_4d_float
        module procedure :: rescale_1d_double
        module procedure :: rescale_2d_double
        module procedure :: rescale_3d_double
        module procedure :: rescale_4d_double
    end interface rescale

    ! Rotate array CW or CCW by 90 degrees
    interface rot90
        module procedure :: rot90_1d_int
        module procedure :: rot90_2d_int
        module procedure :: rot90_3d_int
        module procedure :: rot90_1d_float
        module procedure :: rot90_2d_float
        module procedure :: rot90_3d_float
        module procedure :: rot90_1d_double
        module procedure :: rot90_2d_double
        module procedure :: rot90_3d_double
        module procedure :: rot90_1d_complex
        module procedure :: rot90_2d_complex
        module procedure :: rot90_3d_complex
        module procedure :: rot90_1d_dcomplex
        module procedure :: rot90_2d_dcomplex
        module procedure :: rot90_3d_dcomplex
        module procedure :: rot90_1d_logical
        module procedure :: rot90_2d_logical
        module procedure :: rot90_3d_logical
    end interface rot90

    ! Slice an (n - 1)-dimensional array from an n-dimensional array
    interface slice
        module procedure :: slice_1d_int
        module procedure :: slice_1d_float
        module procedure :: slice_1d_double
        module procedure :: slice_1d_complex
        module procedure :: slice_1d_dcomplex
        module procedure :: slice_1d_logical
        module procedure :: slice_2d_int
        module procedure :: slice_2d_float
        module procedure :: slice_2d_double
        module procedure :: slice_2d_complex
        module procedure :: slice_2d_dcomplex
        module procedure :: slice_2d_logical
        module procedure :: slice_3d_int
        module procedure :: slice_3d_float
        module procedure :: slice_3d_double
        module procedure :: slice_3d_complex
        module procedure :: slice_3d_dcomplex
        module procedure :: slice_3d_logical
    end interface slice

    ! Any element in another array
    interface any_in
        module procedure :: any_in_1d_int
        module procedure :: any_in_1d_float
        module procedure :: any_in_1d_double
        module procedure :: any_in_1d_complex
        module procedure :: any_in_1d_dcomplex
        module procedure :: any_in_1d_logical
        module procedure :: any_in_1d_string
    end interface any_in

    ! All elements must in another array
    interface all_in
        module procedure :: all_in_1d_int
        module procedure :: all_in_1d_float
        module procedure :: all_in_1d_double
        module procedure :: all_in_1d_complex
        module procedure :: all_in_1d_dcomplex
        module procedure :: all_in_1d_logical
        module procedure :: all_in_1d_string
    end interface all_in

    interface remove_any_in
        module procedure :: remove_any_in_1d_int
        module procedure :: remove_any_in_1d_float
        module procedure :: remove_any_in_1d_double
        module procedure :: remove_any_in_1d_complex
        module procedure :: remove_any_in_1d_dcomplex
        module procedure :: remove_any_in_1d_logical
        module procedure :: remove_any_in_1d_string
    end interface remove_any_in

    interface binarize
        module procedure :: binarize_1d_int
        module procedure :: binarize_1d_float
        module procedure :: binarize_1d_double
        module procedure :: binarize_2d_int
        module procedure :: binarize_2d_float
        module procedure :: binarize_2d_double
        module procedure :: binarize_3d_int
        module procedure :: binarize_3d_float
        module procedure :: binarize_3d_double
        module procedure :: binarize_4d_int
        module procedure :: binarize_4d_float
        module procedure :: binarize_4d_double
    end interface binarize

    ! Repeat an array
    interface tile
        module procedure :: tile_1d_int
        module procedure :: tile_1d_float
        module procedure :: tile_1d_double
        module procedure :: tile_1d_complex
        module procedure :: tile_1d_dcomplex
        module procedure :: tile_1d_logical
        module procedure :: tile_2d_int
        module procedure :: tile_2d_float
        module procedure :: tile_2d_double
        module procedure :: tile_2d_complex
        module procedure :: tile_2d_dcomplex
        module procedure :: tile_2d_logical
        module procedure :: tile_3d_int
        module procedure :: tile_3d_float
        module procedure :: tile_3d_double
        module procedure :: tile_3d_complex
        module procedure :: tile_3d_dcomplex
        module procedure :: tile_3d_logical
        module procedure :: tile_4d_int
        module procedure :: tile_4d_float
        module procedure :: tile_4d_double
        module procedure :: tile_4d_complex
        module procedure :: tile_4d_dcomplex
        module procedure :: tile_4d_logical
        ! For string, the only useful might be 1d tiling
        module procedure :: tile_1d_string
    end interface tile

    private
    public :: adjust
    public :: as_scalar
    public :: crop_array
    public :: crop
    public :: cross
    public :: flatten
    public :: flip
    public :: ifelse
    public :: mask
    public :: pad_array
    public :: pad
    public :: permute_array
    public :: permute
    public :: rescale
    public :: rot90
    public :: slice
    public :: ndigits
    public :: breakint
    public :: any_in
    public :: all_in
    public :: remove_any_in
    public :: binarize
    public :: tile

contains

    !================================================================
    ! Adjust the shape of an array
#define T int
#define TT integer
#include "template_adjust.f90"

#define T float
#define TT real
#include "template_adjust.f90"

#define T double
#define TT double precision
#include "template_adjust.f90"

#define T complex
#define TT complex
#include "template_adjust.f90"

#define T dcomplex
#define TT double complex
#include "template_adjust.f90"

#define T logical
#define TT logical
#include "template_adjust.f90"

    !================================================================
    ! As scalar
#define T int
#define TT integer
#include "template_as_scalar.f90"

#define T float
#define TT real
#include "template_as_scalar.f90"

#define T double
#define TT double precision
#include "template_as_scalar.f90"

#define T complex
#define TT complex
#include "template_as_scalar.f90"

#define T dcomplex
#define TT double complex
#include "template_as_scalar.f90"

#define T logical
#define TT logical
#include "template_as_scalar.f90"

    !================================================================
    ! Cross product
#define T int
#define TT integer
#include "template_cross_product.f90"

#define T float
#define TT real
#include "template_cross_product.f90"

#define T double
#define TT double precision
#include "template_cross_product.f90"

#define T complex
#define TT complex
#include "template_cross_product.f90"

#define T dcomplex
#define TT double complex
#include "template_cross_product.f90"

    !================================================================
    ! Crop
#define T int
#define TT integer
#include "template_crop.f90"

#define T float
#define TT real
#include "template_crop.f90"

#define T double
#define TT double precision
#include "template_crop.f90"

#define T complex
#define TT complex
#include "template_crop.f90"

#define T dcomplex
#define TT double complex
#include "template_crop.f90"

#define T logical
#define TT logical
#include "template_crop.f90"

    !================================================================
    ! Flatten array
#define T int
#define TT integer
#include "template_flatten.f90"

#define T float
#define TT real
#include "template_flatten.f90"

#define T double
#define TT double precision
#include "template_flatten.f90"

#define T complex
#define TT complex
#include "template_flatten.f90"

#define T dcomplex
#define TT double complex
#include "template_flatten.f90"

#define T logical
#define TT logical
#include "template_flatten.f90"

    !================================================================
    ! Flip array
#define T int
#define TT integer
#include "template_flip.f90"

#define T float
#define TT real
#include "template_flip.f90"

#define T double
#define TT double precision
#include "template_flip.f90"

#define T complex
#define TT complex
#include "template_flip.f90"

#define T dcomplex
#define TT double complex
#include "template_flip.f90"

#define T logical
#define TT logical
#include "template_flip.f90"

    !================================================================
    ! If-then-else ternary operation for scalar and array
#define T int
#define TT integer
#include "template_ifelse.f90"

#define T float
#define TT real
#include "template_ifelse.f90"

#define T double
#define TT double precision
#include "template_ifelse.f90"

#define T complex
#define TT complex
#include "template_ifelse.f90"

#define T dcomplex
#define TT double complex
#include "template_ifelse.f90"

#define T logical
#define TT logical
#include "template_ifelse.f90"

    ! String requires special routine
    function ifelse_string(condition, a, b) result(c)

        logical :: condition
        character(len=*) :: a, b
        character(:), allocatable :: c

        if (condition) then
            c = a
        else
            c = b
        end if

    end function ifelse_string

    !================================================================
    ! Mask array
#define T int
#define TT integer
#include "template_mask.f90"

#define T float
#define TT real
#include "template_mask.f90"

#define T double
#define TT double precision
#include "template_mask.f90"

#define T complex
#define TT complex
#include "template_mask.f90"

#define T dcomplex
#define TT double complex
#include "template_mask.f90"

#define T logical
#define TT logical
#include "template_mask.f90"

    !================================================================
    ! Pad array
#define T int
#define TT integer
#define TTT integer
#define DEFAULT_VALUE 0
#include "template_pad.f90"

#define T float
#define TT real
#define TTT real
#define DEFAULT_VALUE 0.0
#include "template_pad.f90"

#define T double
#define TT double precision
#define TTT double precision
#define DEFAULT_VALUE 0.0d0
#include "template_pad.f90"

#define T complex
#define TT complex
#define TTT complex
#define DEFAULT_VALUE cmplx(0.0, 0.0)
#include "template_pad.f90"

#define T dcomplex
#define TT double complex
#define TTT double complex
#define DEFAULT_VALUE dcmplx(0.0d0, 0.0d0)
#include "template_pad.f90"

#define T logical
#define TT logical
#define TTT logical
#define DEFAULT_VALUE .false.
#include "template_pad.f90"

#define T string
#define TT character(len=*)
#define TTT character(len=1024)
#define DEFAULT_VALUE ''
#include "template_pad.f90"

    !================================================================
    ! Set operations
#define T int
#define TT integer
#define TTT integer
#include "template_set.f90"

#define T float
#define TT real
#define TTT real
#include "template_set.f90"

#define T double
#define TT double precision
#define TTT double precision
#include "template_set.f90"

#define T complex
#define TT complex
#define TTT complex
#include "template_set.f90"

#define T dcomplex
#define TT double complex
#define TTT double complex
#include "template_set.f90"

#define T logical
#define TT logical
#define TTT logical
#include "template_set.f90"

#define T string
#define TT character(len=*)
#define TTT character(len=1024)
#include "template_set.f90"

    !================================================================
    ! Permute array
    !
    ! ................ Fortran reshape usage ................
    !
    !       reshape(w, shape=[], order=[])
    !
    ! Note that the order indicates "where do the old dims go in new array"
    ! e.g., reshape the array w(35,30,40) to w(40,35,30):
    !
    !       reshape(w, shape=[40,35,30], order=[?])
    !
    ! then order should be order=[2, 3, 1], because:
    ! 35 becomes dimension 2 in the new array
    ! 30 becomes dimension 3 in the new array
    ! 40 becomes dimension 1 in the new array
    !
    function ndigits(x) result(n)

        integer, intent(in) :: x
        integer :: n

        if (x == 0) then
            n = 1
        else
            n = ceiling(log10(abs(x) + 1.0))
        end if

    end function ndigits

    function breakint(x) result(w)

        integer, intent(in) :: x
        integer, allocatable, dimension(:) :: w

        integer :: rem, i
        integer :: nw

        nw = ndigits(x)
        allocate (w(1:nw))

        rem = x
        do i = 1, nw
            ! Take advantage of integer division
            w(nw - i + 1) = rem - (rem/10)*10
            rem = rem/10
        end do

    end function breakint

#define T int
#define TT integer
#include "template_permute.f90"

#define T float
#define TT real
#include "template_permute.f90"

#define T double
#define TT double precision
#include "template_permute.f90"

#define T complex
#define TT complex
#include "template_permute.f90"

#define T dcomplex
#define TT double complex
#include "template_permute.f90"

#define T logical
#define TT logical
#include "template_permute.f90"

    !================================================================
    ! Rescale array
#define T float
#define TT real
#include "template_rescale.f90"

#define T double
#define TT double precision
#include "template_rescale.f90"

    !================================================================
    ! Rotate array
#define T int
#define TT integer
#include "template_rot90.f90"

#define T float
#define TT real
#include "template_rot90.f90"

#define T double
#define TT double precision
#include "template_rot90.f90"

#define T complex
#define TT complex
#include "template_rot90.f90"

#define T dcomplex
#define TT double complex
#include "template_rot90.f90"

#define T logical
#define TT logical
#include "template_rot90.f90"

    !================================================================
    ! Slice array
#define T int
#define TT integer
#include "template_slice.f90"

#define T float
#define TT real
#include "template_slice.f90"

#define T double
#define TT double precision
#include "template_slice.f90"

#define T complex
#define TT complex
#include "template_slice.f90"

#define T dcomplex
#define TT double complex
#include "template_slice.f90"

#define T logical
#define TT logical
#include "template_slice.f90"

    !================================================================
    ! Repeat array
#define T int
#define TT integer
#include "template_tile.f90"

#define T float
#define TT real
#include "template_tile.f90"

#define T double
#define TT double precision
#include "template_tile.f90"

#define T complex
#define TT complex
#include "template_tile.f90"

#define T dcomplex
#define TT double complex
#include "template_tile.f90"

#define T logical
#define TT logical
#include "template_tile.f90"

    function tile_1d_string(w, n) result(wr)

        character(len=*) :: w
        integer :: n
        character(len=:), allocatable :: wr

        integer :: n1, i

        n1 = len(trim(adjustl(w)))

        allocate (character(len=n*n1) :: wr)
        do i = 1, n
            wr((i - 1)*n1 + 1:i*n1) = w
        end do

    end function tile_1d_string

    !================================================================
    ! Binarize array
#define T int
#define TT integer
#include "template_binarize.f90"

#define T float
#define TT real
#include "template_binarize.f90"

#define T double
#define TT double precision
#include "template_binarize.f90"

end  module libflit_array_operation

