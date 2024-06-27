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


module libflit_array

    use libflit_constants
    use libflit_error

    implicit none

    interface alloc_array
        module procedure :: alloc_array_1d_int
        module procedure :: alloc_array_2d_int
        module procedure :: alloc_array_3d_int
        module procedure :: alloc_array_4d_int
        module procedure :: alloc_array_1d_float
        module procedure :: alloc_array_2d_float
        module procedure :: alloc_array_3d_float
        module procedure :: alloc_array_4d_float
        module procedure :: alloc_array_1d_double
        module procedure :: alloc_array_2d_double
        module procedure :: alloc_array_3d_double
        module procedure :: alloc_array_4d_double
        module procedure :: alloc_array_1d_complex
        module procedure :: alloc_array_2d_complex
        module procedure :: alloc_array_3d_complex
        module procedure :: alloc_array_4d_complex
        module procedure :: alloc_array_1d_dcomplex
        module procedure :: alloc_array_2d_dcomplex
        module procedure :: alloc_array_3d_dcomplex
        module procedure :: alloc_array_4d_dcomplex
        module procedure :: alloc_array_1d_logical
        module procedure :: alloc_array_2d_logical
        module procedure :: alloc_array_3d_logical
        module procedure :: alloc_array_4d_logical
        module procedure :: alloc_array_1d_string
        module procedure :: alloc_array_2d_string
        module procedure :: alloc_array_3d_string
        module procedure :: alloc_array_4d_string
        module procedure :: alloc_large_array_1d_int
        module procedure :: alloc_large_array_1d_float
        module procedure :: alloc_large_array_1d_double
        module procedure :: alloc_large_array_1d_complex
        module procedure :: alloc_large_array_1d_dcomplex
        module procedure :: alloc_large_array_1d_logical
        module procedure :: alloc_large_array_1d_string
    end interface alloc_array

    interface regspace
        module procedure :: regspace_int
        module procedure :: regspace_float
        module procedure :: regspace_double
    end interface regspace

    interface linspace
        module procedure :: linspace_float
        module procedure :: linspace_double
        module procedure :: linspace_complex
        module procedure :: linspace_dcomplex
    end interface linspace

    interface logspace
        module procedure :: logspace_float
        module procedure :: logspace_double
        module procedure :: logspace_complex
        module procedure :: logspace_dcomplex
    end interface logspace

    interface const
        module procedure :: const_1d_int
        module procedure :: const_2d_int
        module procedure :: const_3d_int
        module procedure :: const_4d_int
        module procedure :: const_1d_float
        module procedure :: const_2d_float
        module procedure :: const_3d_float
        module procedure :: const_4d_float
        module procedure :: const_1d_double
        module procedure :: const_2d_double
        module procedure :: const_3d_double
        module procedure :: const_4d_double
        module procedure :: const_1d_complex
        module procedure :: const_2d_complex
        module procedure :: const_3d_complex
        module procedure :: const_4d_complex
        module procedure :: const_1d_dcomplex
        module procedure :: const_2d_dcomplex
        module procedure :: const_3d_dcomplex
        module procedure :: const_4d_dcomplex
        module procedure :: const_1d_logical
        module procedure :: const_2d_logical
        module procedure :: const_3d_logical
        module procedure :: const_4d_logical
        module procedure :: const_1d_string
        module procedure :: const_2d_string
        module procedure :: const_3d_string
        module procedure :: const_4d_string
    end interface const

    interface const_like
        module procedure :: const_like_1d_int
        module procedure :: const_like_2d_int
        module procedure :: const_like_3d_int
        module procedure :: const_like_4d_int
        module procedure :: const_like_1d_float
        module procedure :: const_like_2d_float
        module procedure :: const_like_3d_float
        module procedure :: const_like_4d_float
        module procedure :: const_like_1d_double
        module procedure :: const_like_2d_double
        module procedure :: const_like_3d_double
        module procedure :: const_like_4d_double
        module procedure :: const_like_1d_complex
        module procedure :: const_like_2d_complex
        module procedure :: const_like_3d_complex
        module procedure :: const_like_4d_complex
        module procedure :: const_like_1d_dcomplex
        module procedure :: const_like_2d_dcomplex
        module procedure :: const_like_3d_dcomplex
        module procedure :: const_like_4d_dcomplex
        module procedure :: const_like_1d_logical
        module procedure :: const_like_2d_logical
        module procedure :: const_like_3d_logical
        module procedure :: const_like_4d_logical
        module procedure :: const_like_1d_string
        module procedure :: const_like_2d_string
        module procedure :: const_like_3d_string
        module procedure :: const_like_4d_string
    end interface const_like

    interface zeros
        module procedure :: zeros_1d_float
        module procedure :: zeros_2d_float
        module procedure :: zeros_3d_float
        module procedure :: zeros_4d_float
        module procedure :: zeros_5d_float
    end interface zeros

    interface zeros_like
        module procedure :: zeros_like_1d_int
        module procedure :: zeros_like_2d_int
        module procedure :: zeros_like_3d_int
        module procedure :: zeros_like_4d_int
        module procedure :: zeros_like_1d_float
        module procedure :: zeros_like_2d_float
        module procedure :: zeros_like_3d_float
        module procedure :: zeros_like_4d_float
        module procedure :: zeros_like_1d_double
        module procedure :: zeros_like_2d_double
        module procedure :: zeros_like_3d_double
        module procedure :: zeros_like_4d_double
        module procedure :: zeros_like_1d_complex
        module procedure :: zeros_like_2d_complex
        module procedure :: zeros_like_3d_complex
        module procedure :: zeros_like_4d_complex
        module procedure :: zeros_like_1d_dcomplex
        module procedure :: zeros_like_2d_dcomplex
        module procedure :: zeros_like_3d_dcomplex
        module procedure :: zeros_like_4d_dcomplex
    end interface zeros_like

    interface ones
        module procedure :: ones_1d_float
        module procedure :: ones_2d_float
        module procedure :: ones_3d_float
        module procedure :: ones_4d_float
    end interface ones

    interface ones_like
        module procedure :: ones_like_1d_int
        module procedure :: ones_like_2d_int
        module procedure :: ones_like_3d_int
        module procedure :: ones_like_4d_int
        module procedure :: ones_like_1d_float
        module procedure :: ones_like_2d_float
        module procedure :: ones_like_3d_float
        module procedure :: ones_like_4d_float
        module procedure :: ones_like_1d_double
        module procedure :: ones_like_2d_double
        module procedure :: ones_like_3d_double
        module procedure :: ones_like_4d_double
        module procedure :: ones_like_1d_complex
        module procedure :: ones_like_2d_complex
        module procedure :: ones_like_3d_complex
        module procedure :: ones_like_4d_complex
        module procedure :: ones_like_1d_dcomplex
        module procedure :: ones_like_2d_dcomplex
        module procedure :: ones_like_3d_dcomplex
        module procedure :: ones_like_4d_dcomplex
    end interface ones_like

    interface falses
        module procedure :: falses_1d
        module procedure :: falses_2d
        module procedure :: falses_3d
        module procedure :: falses_4d
    end interface falses

    interface falses_like
        module procedure :: falses_like_1d
        module procedure :: falses_like_2d
        module procedure :: falses_like_3d
        module procedure :: falses_like_4d
    end interface falses_like

    interface trues
        module procedure :: trues_1d
        module procedure :: trues_2d
        module procedure :: trues_3d
        module procedure :: trues_4d
    end interface trues

    interface trues_like
        module procedure :: trues_like_1d
        module procedure :: trues_like_2d
        module procedure :: trues_like_3d
        module procedure :: trues_like_4d
    end interface trues_like

    interface diag
        module procedure :: get_diagonal_2d_int
        module procedure :: get_diagonal_2d_float
        module procedure :: get_diagonal_2d_double
        module procedure :: get_diagonal_2d_complex
        module procedure :: get_diagonal_2d_dcomplex
        module procedure :: get_diagonal_2d_logical
    end interface diag

    interface eye
        module procedure :: eye_2d_int
        module procedure :: eye_2d_float
        module procedure :: eye_2d_complex
        module procedure :: eye_2d_double
        module procedure :: eye_2d_dcomplex
        module procedure :: eye_2d_logical
        module procedure :: eye_2d_from_vector_int
        module procedure :: eye_2d_from_vector_float
        module procedure :: eye_2d_from_vector_complex
        module procedure :: eye_2d_from_vector_double
        module procedure :: eye_2d_from_vector_dcomplex
        module procedure :: eye_2d_from_vector_logical
    end interface eye

    private
    public :: alloc_array
    public :: const
    public :: zeros
    public :: ones
    public :: const_like
    public :: zeros_like
    public :: ones_like
    public :: falses
    public :: trues
    public :: falses_like
    public :: trues_like
    public :: eye
    public :: diag
    public :: regspace
    public :: linspace
    public :: logspace

    ! - Array allocation
    !   allocatable array can directly assign, like x = y even if x
    !   is not allocated, and x's shape and bounds will be same
    !   as thos of y
    ! - Array fast allocation using constant
    !   b = [(0., integer :: i = 1,5)]

contains

    ! Allocate array with safe deallocation
#define T int
#define TT integer
#define TTT integer
#define DEFAULT_VALUE 0
#include "template_alloc.f90"

#define T float
#define TT real
#define TTT real
#define DEFAULT_VALUE 0.0
#include "template_alloc.f90"

#define T double
#define TT double precision
#define TTT double precision
#define DEFAULT_VALUE 0.0d0
#include "template_alloc.f90"

#define T complex
#define TT complex
#define TTT complex
#define DEFAULT_VALUE cmplx(0.0, 0.0)
#include "template_alloc.f90"

#define T dcomplex
#define TT double complex
#define TTT double complex
#define DEFAULT_VALUE dcmplx(0.0d0, 0.0d0)
#include "template_alloc.f90"

#define T logical
#define TT logical
#define TTT logical
#define DEFAULT_VALUE .false.
#include "template_alloc.f90"

#define T string
#define TT character(len=*)
#define TTT character(len=1024)
#define DEFAULT_VALUE ''
#include "template_alloc.f90"

    ! Ones and zeros array
#define T int
#define TT integer
#define TTT integer
#include "template_zero_one.f90"

#define T float
#define TT real
#define TTT real
#include "template_zero_one.f90"

#define T double
#define TT double precision
#define TTT double precision
#include "template_zero_one.f90"

#define T complex
#define TT complex
#define TTT complex
#include "template_zero_one.f90"

#define T dcomplex
#define TT double complex
#define TTT double complex
#include "template_zero_one.f90"

    !===========================================================================
    ! Eye array
#define T int
#define TT integer
#define DEFAULT_VALUE 0
#include "template_eye.f90"

#define T float
#define TT real
#define DEFAULT_VALUE 0.0
#include "template_eye.f90"

#define T double
#define TT double precision
#define DEFAULT_VALUE 0.0d0
#include "template_eye.f90"

#define T complex
#define TT complex
#define DEFAULT_VALUE cmplx(0.0, 0.0)
#include "template_eye.f90"

#define T dcomplex
#define TT double complex
#define DEFAULT_VALUE dcmplx(0.0d0, 0.0d0)
#include "template_eye.f90"

#define T logical
#define TT logical
#define DEFAULT_VALUE .false.
#include "template_eye.f90"

    !===========================================================================
    ! Constant-value array
#define T int
#define TT integer
#define TTT integer
#include "template_const.f90"

#define T float
#define TT real
#define TTT real
#include "template_const.f90"

#define T double
#define TT double precision
#define TTT double precision
#include "template_const.f90"

#define T complex
#define TT complex
#define TTT complex
#include "template_const.f90"

#define T dcomplex
#define TT double complex
#define TTT double complex
#include "template_const.f90"

#define T logical
#define TT logical
#define TTT logical
#include "template_const.f90"

#define T string
#define TT character(len=*)
#define TTT character(len=1024)
#include "template_const.f90"

    !
    !> Return a linearly increasing or decreasing integer array
    !>        specified by beg, end, and interval
    !
    function regspace_int(beg, step, end) result(w)

        integer, intent(in) :: beg, step, end

        integer, allocatable, dimension(:) :: w
        integer :: n, i

        n = nint((end - beg + 0.0d0)/step) + 1
        if (beg + (n - 1)*step > end) then
            n = n - 1
        end if
        if (n == 0) then
            n = 1
        end if

        allocate (w(1:n))
        do i = 1, n
            w(i) = beg + (i - 1)*step
        end do

    end function regspace_int

    !
    !> Return a linearly increasing or decreasing float array
    !>        specified by beg, end, and interval
    !
    function regspace_float(beg, step, end) result(w)

        real, intent(in) :: beg, step, end

        real, allocatable, dimension(:) :: w
        integer :: n, i

        n = nint((end - beg)/step) + 1
        if (beg + (n - 1)*step > end) then
            n = n - 1
        end if
        if (n == 0) then
            n = 1
        end if

        allocate (w(1:n))
        do i = 1, n
            w(i) = beg + (i - 1)*step
        end do

    end function regspace_float

    !
    !> Return a linearly increasing or decreasing double array
    !>        specified by beg, end, and interval
    !
    function regspace_double(beg, step, end) result(w)

        double precision, intent(in) :: beg, step, end

        double precision, allocatable, dimension(:) :: w
        integer :: n, i

        n = nint((end - beg)/step) + 1
        if (beg + (n - 1)*step > end) then
            n = n - 1
        end if
        if (n == 0) then
            n = 1
        end if

        allocate (w(1:n))
        do i = 1, n
            w(i) = beg + (i - 1)*step
        end do

    end function regspace_double

    !
    !> Return a linearly increasing or decreasing float array
    !>        specified by beg, end, and number of output elements
    !
    function linspace_float(beg, end, nstep) result(w)

        real, intent(in) :: beg, end
        integer, intent(in) :: nstep

        real, allocatable, dimension(:) :: w
        real :: step
        integer :: i

        if (nstep == 1) then
            step = end - beg
        else
            step = (end - beg)/(nstep - 1.0)
        end if

        allocate (w(1:nstep))
        do i = 1, nstep
            w(i) = beg + (i - 1)*step
        end do

    end function linspace_float

    !
    !> Return a linearly increasing or decreasing double array
    !>        specified by beg, end, and number of output elements
    !
    function linspace_double(beg, end, nstep) result(w)

        double precision, intent(in) :: beg, end
        integer, intent(in) :: nstep

        double precision, allocatable, dimension(:) :: w
        double precision :: step
        integer :: i

        if (nstep == 1) then
            step = end - beg
        else
            step = (end - beg)/(nstep - 1.0)
        end if

        allocate (w(1:nstep))
        do i = 1, nstep
            w(i) = beg + (i - 1)*step
        end do

    end function linspace_double

    !
    !> Return a linearly increasing or decreasing complex array
    !>        specified by beg, end, and number of output elements
    !
    function linspace_complex(beg, end, nstep) result(w)

        complex, intent(in) :: beg, end
        integer, intent(in) :: nstep

        complex, allocatable, dimension(:) :: w
        complex :: step
        integer :: i

        if (nstep == 1) then
            step = end - beg
        else
            step = (end - beg)/(nstep - 1.0)
        end if

        allocate (w(1:nstep))
        do i = 1, nstep
            w(i) = beg + (i - 1)*step
        end do

    end function linspace_complex

    !
    !> Return a linearly increasing or decreasing double complex array
    !>        specified by beg, end, and number of output elements
    !
    function linspace_dcomplex(beg, end, nstep) result(w)

        double complex, intent(in) :: beg, end
        integer, intent(in) :: nstep

        double complex, allocatable, dimension(:) :: w
        double complex :: step
        integer :: i

        if (nstep == 1) then
            step = end - beg
        else
            step = (end - beg)/(nstep - 1.0)
        end if

        allocate (w(1:nstep))
        do i = 1, nstep
            w(i) = beg + (i - 1)*step
        end do

    end function linspace_dcomplex

    !
    !> Return a linearly increasing or decreasing float array
    !>        specified by beg, end, and number of output elements
    !
    function logspace_float(beg, end, nstep, base) result(w)

        real, intent(in) :: beg, end
        integer, intent(in) :: nstep
        real, intent(in), optional :: base

        real, allocatable, dimension(:) :: w
        real :: step, logbase
        integer :: i

        if (nstep == 1) then
            step = end - beg
        else
            step = (end - beg)/(nstep - 1.0)
        end if

        if (present(base)) then
            logbase = base
        else
            logbase = 10.0
        end if

        allocate (w(1:nstep))
        do i = 1, nstep
            w(i) = logbase**(beg + (i - 1)*step)
        end do

    end function logspace_float

    !
    !> Return a linearly increasing or decreasing double array
    !>        specified by beg, end, and number of output elements
    !
    function logspace_double(beg, end, nstep, base) result(w)

        double precision, intent(in) :: beg, end
        integer, intent(in) :: nstep
        double precision, intent(in), optional :: base

        double precision, allocatable, dimension(:) :: w
        double precision :: step, logbase
        integer :: i

        if (nstep == 1) then
            step = end - beg
        else
            step = (end - beg)/(nstep - 1.0)
        end if

        if (present(base)) then
            logbase = base
        else
            logbase = 10.0d0
        end if

        allocate (w(1:nstep))
        do i = 1, nstep
            w(i) = logbase**(beg + (i - 1)*step)
        end do

    end function logspace_double

    !
    !> Return a linearly increasing or decreasing complex array
    !>        specified by beg, end, and number of output elements
    !
    function logspace_complex(beg, end, nstep, base) result(w)

        complex, intent(in) :: beg, end
        integer, intent(in) :: nstep
        complex, intent(in), optional :: base

        complex, allocatable, dimension(:) :: w
        complex :: step, logbase
        integer :: i

        if (nstep == 1) then
            step = end - beg
        else
            step = (end - beg)/(nstep - 1.0)
        end if

        if (present(base)) then
            logbase = base
        else
            logbase = cmplx(10.0, 0.0)
        end if

        allocate (w(1:nstep))
        do i = 1, nstep
            w(i) = logbase**(beg + (i - 1)*step)
        end do

    end function logspace_complex

    !
    !> Return a linearly increasing or decreasing double complex array
    !>        specified by beg, end, and number of output elements
    !
    function logspace_dcomplex(beg, end, nstep, base) result(w)

        double complex, intent(in) :: beg, end
        integer, intent(in) :: nstep
        double complex, intent(in), optional :: base

        double complex, allocatable, dimension(:) :: w
        double complex :: step, logbase
        integer :: i

        if (nstep == 1) then
            step = end - beg
        else
            step = (end - beg)/(nstep - 1.0)
        end if

        if (present(base)) then
            logbase = base
        else
            logbase = dcmplx(10.0d0, 0.0d0)
        end if

        allocate (w(1:nstep))
        do i = 1, nstep
            w(i) = logbase**(beg + (i - 1)*step)
        end do

    end function logspace_dcomplex

    !============================================================
    ! Zero-valued array

    function zeros_1d_float(n1) result(w)

        integer :: n1
        real, allocatable, dimension(:) :: w

        allocate (w(1:n1))
        w = 0.0

    end function zeros_1d_float

    function zeros_2d_float(n1, n2) result(w)

        integer :: n1, n2
        real, allocatable, dimension(:, :) :: w

        allocate (w(1:n1, 1:n2))
        w = 0.0

    end function zeros_2d_float

    function zeros_3d_float(n1, n2, n3) result(w)

        integer :: n1, n2, n3
        real, allocatable, dimension(:, :, :) :: w

        allocate (w(1:n1, 1:n2, 1:n3))
        w = 0.0

    end function zeros_3d_float

    function zeros_4d_float(n1, n2, n3, n4) result(w)

        integer :: n1, n2, n3, n4
        real, allocatable, dimension(:, :, :, :) :: w

        allocate (w(1:n1, 1:n2, 1:n3, 1:n4))
        w = 0.0

    end function zeros_4d_float

    function zeros_5d_float(n1, n2, n3, n4, n5) result(w)

        integer :: n1, n2, n3, n4, n5
        real, allocatable, dimension(:, :, :, :, :) :: w

        allocate (w(1:n1, 1:n2, 1:n3, 1:n4, 1:n5))
        w = 0.0

    end function zeros_5d_float

    !============================================================
    ! 1-valued array

    function ones_1d_float(n1) result(w)

        integer :: n1
        real, allocatable, dimension(:) :: w

        allocate (w(1:n1))
        w = 1.0

    end function ones_1d_float

    function ones_2d_float(n1, n2) result(w)

        integer :: n1, n2
        real, allocatable, dimension(:, :) :: w

        allocate (w(1:n1, 1:n2))
        w = 1.0

    end function ones_2d_float

    function ones_3d_float(n1, n2, n3) result(w)

        integer :: n1, n2, n3
        real, allocatable, dimension(:, :, :) :: w

        allocate (w(1:n1, 1:n2, 1:n3))
        w = 1.0

    end function ones_3d_float

    function ones_4d_float(n1, n2, n3, n4) result(w)

        integer :: n1, n2, n3, n4
        real, allocatable, dimension(:, :, :, :) :: w

        allocate (w(1:n1, 1:n2, 1:n3, 1:n4))
        w = 1.0

    end function ones_4d_float

    !===================================================
    ! trues
    function falses_1d(n1) result(w)

        integer :: n1
        logical, allocatable, dimension(:) :: w

        allocate (w(1:n1))
        w = .false.

    end function falses_1d

    function falses_2d(n1, n2) result(w)

        integer :: n1, n2
        logical, allocatable, dimension(:, :) :: w

        allocate (w(1:n1, 1:n2))
        w = .false.

    end function falses_2d

    function falses_3d(n1, n2, n3) result(w)

        integer :: n1, n2, n3
        logical, allocatable, dimension(:, :, :) :: w

        allocate (w(1:n1, 1:n2, 1:n3))
        w = .false.

    end function falses_3d

    function falses_4d(n1, n2, n3, n4) result(w)

        integer :: n1, n2, n3, n4
        logical, allocatable, dimension(:, :, :, :) :: w

        allocate (w(1:n1, 1:n2, 1:n3, 1:n4))
        w = .false.

    end function falses_4d

    function trues_1d(n1) result(w)

        integer :: n1
        logical, allocatable, dimension(:) :: w

        allocate (w(1:n1))
        w = .true.

    end function trues_1d

    function trues_2d(n1, n2) result(w)

        integer :: n1, n2
        logical, allocatable, dimension(:, :) :: w

        allocate (w(1:n1, 1:n2))
        w = .true.

    end function trues_2d

    function trues_3d(n1, n2, n3) result(w)

        integer :: n1, n2, n3
        logical, allocatable, dimension(:, :, :) :: w

        allocate (w(1:n1, 1:n2, 1:n3))
        w = .true.

    end function trues_3d

    function trues_4d(n1, n2, n3, n4) result(w)

        integer :: n1, n2, n3, n4
        logical, allocatable, dimension(:, :, :, :) :: w

        allocate (w(1:n1, 1:n2, 1:n3, 1:n4))
        w = .true.

    end function trues_4d

    function falses_like_1d(w) result(wr)

        logical, dimension(:), intent(in) :: w
        logical, allocatable, dimension(:) :: wr

        allocate (wr(1:size(w)))
        wr = .false.

    end function falses_like_1d

    function trues_like_1d(w) result(wr)

        logical, dimension(:), intent(in) :: w
        logical, allocatable, dimension(:) :: wr

        allocate (wr(1:size(w)))
        wr = .true.

    end function trues_like_1d

    function falses_like_2d(w) result(wr)

        logical, dimension(:, :), intent(in) :: w
        logical, allocatable, dimension(:, :) :: wr

        allocate (wr(1:size(w, 1), 1:size(w, 2)))
        wr = .false.

    end function falses_like_2d

    function trues_like_2d(w) result(wr)

        logical, dimension(:, :), intent(in) :: w
        logical, allocatable, dimension(:, :) :: wr

        allocate (wr(1:size(w, 1), 1:size(w, 2)))
        wr = .true.

    end function trues_like_2d

    function falses_like_3d(w) result(wr)

        logical, dimension(:, :, :), intent(in) :: w
        logical, allocatable, dimension(:, :, :) :: wr

        allocate (wr(1:size(w, 1), 1:size(w, 2), 1:size(w, 3)))
        wr = .false.

    end function falses_like_3d

    function trues_like_3d(w) result(wr)

        logical, dimension(:, :, :), intent(in) :: w
        logical, allocatable, dimension(:, :, :) :: wr

        allocate (wr(1:size(w, 1), 1:size(w, 2), 1:size(w, 3)))
        wr = .true.

    end function trues_like_3d

    function falses_like_4d(w) result(wr)

        logical, dimension(:, :, :, :), intent(in) :: w
        logical, allocatable, dimension(:, :, :, :) :: wr

        allocate (wr(1:size(w, 1), 1:size(w, 2), 1:size(w, 3), 1:size(w, 4)))
        wr = .false.

    end function falses_like_4d

    function trues_like_4d(w) result(wr)

        logical, dimension(:, :, :, :), intent(in) :: w
        logical, allocatable, dimension(:, :, :, :) :: wr

        allocate (wr(1:size(w, 1), 1:size(w, 2), 1:size(w, 3), 1:size(w, 4)))
        wr = .true.

    end function trues_like_4d

end module libflit_array

!    real :: a(n)                ! An explicit shape array
!    real :: b(:)                ! An assumed shape array
!    real, allocatable :: c(:)   ! A deferred shape array (allocatable)
!    real, pointer :: d(:)       ! A deferred shape array (pointer)
!    real :: e(*)                ! An assumed size array
!    real :: f(..)               ! An assumed rank array/scalar
