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


module libflit_domain_decomposition

    use libflit_array
    use libflit_interp
    use libflit_error
    use libflit_calculus
    use libflit_utility
    use libflit_array_operation

    implicit none

    private

    interface divide_domain
        module procedure :: divide_domain_1
        module procedure :: divide_domain_1_weight
    end interface

    interface cut
        module procedure :: cut_simple
        module procedure :: cut_simple_large
        module procedure :: cut_pwconst
    end interface

    interface split_list
        module procedure :: split_list_int
        module procedure :: split_list_float
        module procedure :: split_list_double
        module procedure :: split_list_complex
        module procedure :: split_list_dcomplex
        module procedure :: split_list_logical
    end interface split_list

    public :: divide_domain
    public :: cut
    public :: split_list

contains

    !
    !> Divide a 1D domain
    !
    !> @param[in] n - number of grid points
    !> @param[in] nhead - number of additional points at the head
    !> @param[in] ntail - number of additional points at the tail
    !> @param[in] weight - computational complexity of the additional points
    !> @param[in] nblk - number of blocks after division
    !
    !> @param[out] blkrange - the calculated ranges of blocks
    !
    subroutine divide_domain_1(n, nhead, ntail, weight, nblk, blkrange)

        integer, intent(in) :: n, nblk, nhead, ntail
        real, intent(in) :: weight
        integer, allocatable, dimension(:, :), intent(inout) :: blkrange

        integer, dimension(1:nblk) :: blocksizes
        real, allocatable, dimension(:) :: w
        integer :: r_start, r_end, i, j
        real :: plength, d
        real :: avg, last
        integer :: lb
        integer :: bkbeg, bkend
        integer, allocatable, dimension(:, :) :: tmpbk

        if (.not. allocated(blkrange)) then
            bkbeg = 0
            bkend = nblk - 1
            allocate (blkrange(bkbeg:bkend, 1:2))
            blkrange = 0
        else
            bkbeg = lbound(blkrange, 1)
            bkend = ubound(blkrange, 1)
        end if

        if (nhead == 0 .and. ntail == 0) then

            avg = (n + 0.0d0)/(nblk + 0.0d0)
            last = 1.0
            lb = bkbeg
            do while (last <= n)
                blkrange(lb, 1) = nint(last)
                blkrange(lb, 2) = nint(last + avg) - 1
                last = last + avg
                lb = lb + 1
            end do

        else

            ! Weight
            w = ones(n + nhead + ntail + 1)
            w(1:nhead) = weight
            w(n + nhead + 1:n + nhead + ntail) = weight
            w(n + nhead + ntail + 1) = 0.0

            ! Pseudo-length of the domain
            plength = n + (nhead + ntail)*weight

            ! Average pseudo-length
            d = plength/nblk

            ! Initialization
            r_start = 1
            r_end = 1
            j = 1

            ! Solve for size
            do while (j <= nblk)
                do while ((r_end < n + nhead + ntail) &
                        .and. (abs(sum(w(r_start:r_end + 1)) - d) <= abs(sum(w(r_start:r_end)) - d)))
                    r_end = r_end + 1
                end do
                blocksizes(j) = r_end - r_start + 1
                r_start = r_end + 1
                r_end = r_end + 1
                j = j + 1
            end do

            ! Ensure the last blocks has correct size
            blocksizes(nblk) = n + nhead + ntail - sum(blocksizes(1:nblk - 1))

            ! Generate blkrange
            allocate (tmpbk(1:nblk, 1:2))
            ! Here to keep consistency with other functions in this module, the additional nhead points are not subtracted
            ! otherwise, it should be - nhead + 1
            tmpbk(1, 1) = 1
            do i = 2, nblk
                ! Same here, otherwise it should be -nhead + ...
                tmpbk(i, 1) = 1 + sum(blocksizes(1:i - 1))
            end do
            do i = 1, nblk - 1
                tmpbk(i, 2) = tmpbk(i + 1, 1) - 1
            end do
            ! Same here, otherwise it should be n + ntail
            tmpbk(nblk, 2) = n + nhead + ntail

            blkrange(bkbeg:bkend, 1) = tmpbk(1:nblk, 1)
            blkrange(bkbeg:bkend, 2) = tmpbk(1:nblk, 2)

        end if

    end subroutine divide_domain_1

    !
    !> Divide a 1D domain where elements may have different weights
    !> into blocks with approximately equal sum
    !
    subroutine divide_domain_1_weight(weight, nblk, blkrange)

        real, dimension(:), intent(in) :: weight
        integer, intent(in) :: nblk
        integer, allocatable, dimension(:, :), intent(inout) :: blkrange

        integer, dimension(1:nblk) :: blocksizes
        integer :: r_start, r_end, i, j, bkbeg, bkend
        real :: d
        integer, allocatable, dimension(:, :) :: tmpbk
        integer :: n

        if (.not. allocated(blkrange)) then
            bkbeg = 0
            bkend = nblk - 1
            allocate (blkrange(bkbeg:bkend, 1:2))
            blkrange = 0
        else
            bkbeg = lbound(blkrange, 1)
            bkend = ubound(blkrange, 1)
        end if

        ! Length of the array
        n = size(weight)

        ! Average pseudo-length
        d = sum(weight)/nblk

        ! Initialization
        r_start = 1
        r_end = 1
        j = 1

        ! Solve for size
        do while (j <= nblk)
            do while (r_end < n &
                    .and. abs(sum(weight(r_start:r_end + 1)) - d) <= abs(sum(weight(r_start:r_end)) - d))
                r_end = r_end + 1
            end do
            blocksizes(j) = r_end - r_start + 1
            r_start = r_end + 1
            r_end = r_end + 1
            j = j + 1
        end do

        ! Ensure the last blocks has correct size
        blocksizes(nblk) = n - sum(blocksizes(1:nblk - 1))

        ! Generate blkrange
        allocate (tmpbk(1:nblk, 1:2))
        tmpbk(1, 1) = 1
        do i = 2, nblk
            tmpbk(i, 1) = 1 + sum(blocksizes(1:i - 1))
        end do
        do i = 1, nblk - 1
            tmpbk(i, 2) = tmpbk(i + 1, 1) - 1
        end do
        tmpbk(nblk, 2) = n

        blkrange(bkbeg:bkend, 1) = tmpbk(1:nblk, 1)
        blkrange(bkbeg:bkend, 2) = tmpbk(1:nblk, 2)

    end subroutine divide_domain_1_weight

    !
    !> Cut a 1D domain into blocks
    !
    subroutine cut_simple(nbeg, nend, nblk, blkrange)

        ! arguments
        integer, intent(in) :: nbeg, nend, nblk
        integer, allocatable, dimension(:, :), intent(inout) :: blkrange

        ! local variables
        real :: avg, last
        integer :: lb
        integer :: bkbeg, bkend

        if (.not. allocated(blkrange)) then
            bkbeg = 1
            bkend = nblk
            allocate (blkrange(bkbeg:bkend, 1:2))
            blkrange(:, 1) = nbeg - 1
            blkrange(:, 2) = nbeg - 2
        else
            bkbeg = lbound(blkrange, 1)
            bkend = ubound(blkrange, 1)
            blkrange(:, 1) = nbeg - 1
            blkrange(:, 2) = nbeg - 2
        end if

        ! Cut into blocks with min block size = 1
        avg = max(1.0, (nend - nbeg + 1.0d0)/(nblk + 0.0d0))
        last = nbeg
        lb = bkbeg
        do while (last <= nend)
            blkrange(lb, 1) = nint(last)
            blkrange(lb, 2) = nint(last + avg) - 1
            last = last + avg
            lb = lb + 1
        end do

    end subroutine cut_simple

    !
    !> Cut a large 1D domain into blocks
    !
    subroutine cut_simple_large(nbeg, nend, nblk, blkrange)

        ! arguments
        integer(kind=8), intent(in) :: nbeg, nend, nblk
        integer(kind=8), allocatable, dimension(:, :), intent(inout) :: blkrange

        ! local variables
        real :: avg, last
        integer(kind=8) :: lb
        integer(kind=8) :: bkbeg, bkend

        if (.not. allocated(blkrange)) then
            bkbeg = 1
            bkend = nblk
            allocate (blkrange(bkbeg:bkend, 1:2))
            blkrange(:, 1) = nbeg - 1
            blkrange(:, 2) = nbeg - 2
        else
            bkbeg = lbound(blkrange, 1)
            bkend = ubound(blkrange, 1)
            blkrange(:, 1) = nbeg - 1
            blkrange(:, 2) = nbeg - 2
        end if

        ! Cut into blocks with min block size = 1
        avg = max(1.0, (nend - nbeg + 1.0d0)/(nblk + 0.0d0))
        last = nbeg
        lb = bkbeg
        do while (last <= nend)
            blkrange(lb, 1) = nint(last)
            blkrange(lb, 2) = nint(last + avg) - 1
            last = last + avg
            lb = lb + 1
        end do

    end subroutine cut_simple_large

    !
    !> Cut a 1D piecewise domain into blocks based on piece interfaces
    !
    subroutine cut_pwconst(pwconst, nblk, blkrange)

        real, dimension(:), intent(in) :: pwconst
        integer, intent(in) :: nblk
        integer, allocatable, dimension(:, :), intent(inout) :: blkrange

        integer :: i, l, nw
        real, allocatable, dimension(:) :: w

        call alloc_array(blkrange, [1, nblk, 1, 2])

        ! Find piecewise interfaces using finite difference
        nw = size(pwconst)
        w = deriv(pwconst, method='forward')

        l = index_first_nonzero(w)
        blkrange(1, 1) = 1
        blkrange(1, 2) = l
        w(l) = 0

        do i = 2, nblk
            l = min(index_first_nonzero(w), nw)
            blkrange(i, 1) = blkrange(i - 1, 2) + 1
            blkrange(i, 2) = l
            w(l) = 0
        end do

    end subroutine cut_pwconst

    ! Split list into chuncks
#define T int
#define TT integer
#define TTT integer
#include "template_split_list.f90"

#define T float
#define TT real
#include "template_split_list.f90"

#define T double
#define TT double precision
#include "template_split_list.f90"

#define T complex
#define TT complex
#include "template_split_list.f90"

#define T dcomplex
#define TT double complex
#include "template_split_list.f90"

#define T logical
#define TT logical
#include "template_split_list.f90"

end module libflit_domain_decomposition
