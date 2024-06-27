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

#define PASTE(X)            X
#define PASTE2(X)           PASTE(X)_
#define CONCATHELP(X, Y)    PASTE2(X)Y
#define CONCAT(X, Y)        CONCATHELP(X, Y)

#define merge_sort_index_     CONCAT(merge_sort_index, T)
#define quick_sort_     CONCAT(quick_sort, T)
#define qksort_     CONCAT(qksort, T)
#define split_      CONCAT(split, T)
#define quick_sort_index_       CONCAT(quick_sort_index, T)
#define qksort_index_       CONCAT(qksort_index, T)
#define qksort_2d_      CONCAT(qksort_2d, T)

subroutine merge_sort_index_(r, d)

    TT, intent(in), dimension(:) :: r
    integer, intent(out), dimension(size(r)) :: d

    integer, dimension(size(r)) :: il

    integer :: stepsize
    integer :: i, j, n, left, k, ksize

    n = size(r)

    do i = 1, n
        d(i) = i
    end do

    if (n == 1) then
        return
    end if

    stepsize = 1
    do while (stepsize < n)
        do left = 1, n - stepsize, stepsize*2
            i = left
            j = left + stepsize
            ksize = min(stepsize*2, n - left + 1)
            k = 1

            do while (i < left + stepsize .and. j < left + ksize)
                if (r(d(i)) > r(d(j))) then
                    il(k) = d(i)
                    i = i + 1
                    k = k + 1
                else
                    il(k) = d(j)
                    j = j + 1
                    k = k + 1
                endif
            enddo

            if (i < left + stepsize) then
                ! fill up remaining from left
                il(k:ksize) = d(i:left + stepsize - 1)
            else
                ! fill up remaining from right
                il(k:ksize) = d(j:left + ksize - 1)
            endif
            d(left:left + ksize - 1) = il(1:ksize)
        end do
        stepsize = stepsize*2

    end do

end subroutine merge_sort_index_

!
!> Recursive quick-sort routine sorting array into ascending or descending order
!
!> Reformatted and modified from
!>      https://gist.github.com/t-nissie/479f0f16966925fa29ea#file-quicksort-f
!> Original author:
!>      Takeshi Nishimatsu
!> License:
!>      GPLv3
!
recursive subroutine quick_sort_(a, order)

    TT, dimension(:), intent(inout) :: a
    integer, intent(in), optional :: order

    TT :: x, t
    integer :: first, last
    integer :: i, j

    first = 1
    last = size(a, 1)

    x = a((first + last)/2)
    i = first
    j = last
    do
        do while (a(i) < x)
            i = i + 1
        end do
        do while (x < a(j))
            j = j - 1
        end do
        if (i >= j) exit
        t = a(i)
        a(i) = a(j)
        a(j) = t

        i = i + 1
        j = j - 1
    end do

    if (first < i - 1) then
        call quick_sort_(a(first:i - 1))
    end if
    if (j + 1 < last) then
        call quick_sort_(a(j + 1:last))
    end if

    if (present(order)) then
        if (order == -1) then
            a(1:size(a, 1)) = a(size(a, 1):1:-1)
        end if
    end if

end subroutine quick_sort_

!
!> Quick sort array into ascending or descending order
!
function qksort_(a, order) result(b)

    TT, dimension(:) :: a
    integer, intent(in), optional :: order

    TT, allocatable, dimension(:) :: b

    allocate (b(1:size(a)), source=a)

    if (present(order)) then
        call quick_sort_(b, order)
    else
        call quick_sort_(b)
    end if

end function qksort_

!
! Split a list into two sublists, using the first element
! as a pivot, and return the position of the element about which the
! list was divided. Local variables used are:
! left      : position of the first element
! right     : position of the last element
! pivot     : pivot element
! swap      : used to swap elements
!
! accepts:  array a and positions low and high of the first and last elements
! returns:  array a (modified) with elements in ascending order
!
subroutine split_(a, low, high, mid, indices)

    TT, dimension(:), intent(inout) :: a
    integer, intent(in) :: low, high
    integer, intent(out) :: mid
    integer, dimension(:), intent(inout) :: indices

    integer ::   left, right
    TT ::  pivot, swap
    integer :: ipivot, iswap

    left = low
    right = high
    pivot = a(low)
    ipivot = indices(low)

    ! repeat the following while left and right haven't met
    do
        if (left >= right) exit

        ! scan right to left to find element < pivot
        do
            if (left >= right .or. a(right) < pivot) exit
            right = right - 1
        end do

        ! scan left to right to find element > pivot
        do
            if (a(left) > pivot) exit
            left = left + 1
        end do

        ! if left and right haven't met, exchange the as
        if (left < right) then
            swap = a(left) ! exchange the array as
            a(left) = a(right)
            a(right) = swap

            iswap = indices(left) ! exchange the indices as
            indices(left) = indices(right)
            indices(right) = iswap
        end if

    end do

    ! switch element in split position with pivot
    a(low) = a(right) ! switch array elems
    a(right) = pivot
    mid = right

    indices(low) = indices(right) ! switch array elems
    indices(right) = ipivot

end subroutine split_

recursive subroutine quick_sort_index_(a, first, last, indices)

    TT, dimension(:), intent(inout) :: a ! array of values
    integer, intent(in)   :: first, last
    integer, dimension(:), intent(inout) :: indices

    integer :: mid

    if (first < last) then ! if list size >= 2
        call split_(a, first, last, mid, indices) ! split it
        call quick_sort_index_(a, first, mid - 1, indices) ! sort left  half
        call quick_sort_index_(a, mid + 1, last, indices) ! sort right half
    end if

end subroutine quick_sort_index_

subroutine qksort_index_(a, indices, order)

    TT, dimension(:), intent(inout) :: a
    integer, allocatable, dimension(:), intent(inout) :: indices
    integer, optional, intent(in) :: order

    integer :: i

    indices = [(i, i=1, size(a))]
    call quick_sort_index_(a, 1, size(a), indices)
    if (present(order)) then
        if (order == -1) then
            a = a(size(a):1:-1)
            indices = indices(size(indices):1:-1)
        end if
    end if

end subroutine qksort_index_

!
!> Sort a float matrix based on one of the columns
!
function qksort_2d_(a, col, order) result(asort)

    TT, dimension(:, :) :: a
    integer, optional :: col, order
    TT, allocatable :: asort(:, :)

    TT, allocatable, dimension(:) :: b
    integer, allocatable, dimension(:) :: indices
    integer :: which_col, n, i, sort_order

    ! Size of array
    n = size(a, 1)

    ! Sort based on which column
    if (present(col)) then
        which_col = col
    else
        which_col = 1
    end if
    which_col = max(1, min(size(a, 2), which_col))

    ! Sort based on which order
    if (present(order)) then
        sort_order = order
    else
        sort_order = 1
    end if

    ! Sort
    allocate (b(1:n), source=a(:, which_col))
    allocate (indices(1:n))
    call qksort_index_(b, indices, sort_order)

    ! Output
    allocate (asort(1:n, 1:size(a, 2)))
    do i = 1, n
        asort(i, :) = a(indices(i), :)
    end do

    deallocate (b, indices)

end function qksort_2d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef merge_sort_index_
#undef quick_sort_
#undef qksort_
#undef split_
#undef quick_sort_index_
#undef qksort_index_
#undef qksort_2d_
