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

#define unique_1d_     CONCAT(unique_1d, T)
#define unique_2d_     CONCAT(unique_2d, T)

function unique_1d_(w) result(wr)

    TT, dimension(:) :: w
    TTT, allocatable, dimension(:) :: wr

    integer :: n, i, l
    logical, allocatable, dimension(:) :: q

    n = size(w)
    call alloc_array(wr, [1, n])
    q = trues(n)

    wr(1) = w(1)
    q(1) = .false.

    l = 2
    ! Check every element in the source array
    do i = 1, n
        ! If the element has not been checked to be unique
        if (q(i)) then
            ! To be a unique element, it must not be equivalent with any confirmed unique element
            if (.not. any(w(i) == wr(1:l - 1))) then
                wr(l) = w(i)
                q(i) = .false.
                l = l + 1
            end if
        end if
    end do

    wr = wr(1:l - 1)

end function unique_1d_

function unique_2d_(w, cols) result(wr)

    TT, dimension(:, :) :: w
    integer, optional :: cols(:)
    TTT, allocatable, dimension(:, :) :: wr

    integer :: n1, n2, i, j, k, l
    logical, allocatable, dimension(:) :: q, duplicate
    integer, allocatable, dimension(:) :: cl

    n1 = size(w, 1)
    n2 = size(w, 2)
    call alloc_array(wr, [1, n1, 1, n2])
    q = trues(n1)

    if (present(cols)) then
        cl = cols
    else
        cl = regspace(1, 1, n2)
    end if
    duplicate = falses(size(cl))

    wr(1, :) = w(1, :)
    l = 2
    ! Check every row in the source array
    do i = 1, n1

        ! If the row has not been checked to be unique
        if (q(i)) then

            ! Check duplication with all previous unique rows
            check_duplicate: do j = 1, l - 1

                ! For each row, assume the row is unique
                duplicate = .false.

                ! Check duplication of each column
                do k = 1, size(cl)
                    if (w(i, cl(k)) == wr(j, cl(k))) then
                        duplicate(k) = .true.
                    end if
                end do

                ! If all elements are the same with some previous row, then this
                ! row is a duplicate row and there is no need to further check
                if (all(duplicate)) then
                    exit check_duplicate
                end if

            end do check_duplicate

            ! If a unique row, then add to unique list
            if (any(.not. duplicate)) then
                wr(l, :) = w(i, :)
                q(i) = .false.
                l = l + 1
            end if

        end if

    end do

    wr = wr(1:l - 1, :)

end function unique_2d_

#undef T
#undef TT
#undef TTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef unique_1d_
#undef unique_2d_
