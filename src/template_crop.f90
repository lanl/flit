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

#define crop_1d_      CONCAT(crop_1d, T)
#define crop_2d_      CONCAT(crop_2d, T)
#define crop_3d_      CONCAT(crop_3d, T)
#define crop_4d_      CONCAT(crop_4d, T)
#define crop_array_1d_      CONCAT(crop_array_1d, T)
#define crop_array_2d_      CONCAT(crop_array_2d, T)
#define crop_array_3d_      CONCAT(crop_array_3d, T)
#define crop_array_4d_      CONCAT(crop_array_4d, T)

subroutine crop_array_1d_(w, range)

    TT, allocatable, dimension(:), intent(inout) :: w
    integer, dimension(1:2), intent(in) :: range

    call assert(size(range) == 2 &
        .and. range(1) >= lbound(w, 1) .and. range(2) <= ubound(w, 1), &
        ' <crop_array_1d> Error: range specification is wrong. ')

    call alloc_array(w, range, &
        source=w(range(1):range(2)))

end subroutine crop_array_1d_

subroutine crop_array_2d_(w, range)

    TT, allocatable, dimension(:, :), intent(inout) :: w
    integer, dimension(1:4), intent(in) :: range

    call assert(size(range) == 4 &
        .and. range(1) >= lbound(w, 1) .and. range(2) <= ubound(w, 1) &
        .and. range(3) >= lbound(w, 2) .and. range(4) <= ubound(w, 2), &
        ' <crop_array_2d> Error: range specification is wrong. ')

    call alloc_array(w, range, &
        source=w(range(1):range(2), &
        range(3):range(4)))

end subroutine crop_array_2d_

subroutine crop_array_3d_(w, range)

    TT, allocatable, dimension(:, :, :), intent(inout) :: w
    integer, dimension(1:6), intent(in) :: range

    call assert(size(range) == 6 &
        .and. range(1) >= lbound(w, 1) .and. range(2) <= ubound(w, 1) &
        .and. range(3) >= lbound(w, 2) .and. range(4) <= ubound(w, 2) &
        .and. range(5) >= lbound(w, 3) .and. range(6) <= ubound(w, 3), &
        ' <crop_array_3d> Error: range specification is wrong. ')

    call alloc_array(w, range, &
        source=w(range(1):range(2), &
        range(3):range(4), &
        range(5):range(6)))

end subroutine crop_array_3d_

subroutine crop_array_4d_(w, range)

    TT, allocatable, dimension(:, :, :, :), intent(inout) :: w
    integer, dimension(1:8), intent(in) :: range

    call assert(size(range) == 8 &
        .and. range(1) >= lbound(w, 1) .and. range(2) <= ubound(w, 1) &
        .and. range(3) >= lbound(w, 2) .and. range(4) <= ubound(w, 2) &
        .and. range(5) >= lbound(w, 3) .and. range(6) <= ubound(w, 3) &
        .and. range(7) >= lbound(w, 4) .and. range(8) <= ubound(w, 4), &
        ' <crop_array_4d> Error: range specification is wrong. ')

    call alloc_array(w, range, &
        source=w(range(1):range(2), &
        range(3):range(4), &
        range(5):range(6), &
        range(7):range(8)))

end subroutine crop_array_4d_

function crop_1d_(w, range) result(wt)

    TT, dimension(:), intent(in) :: w
    integer, dimension(1:2), intent(in) :: range
    TT, allocatable, dimension(:) :: wt

    call assert(size(range) == 2 &
        .and. range(1) >= lbound(w, 1) .and. range(2) <= ubound(w, 1), &
        ' <crop_array_1d> Error: range specification is wrong. ')

    call alloc_array(wt, [1, range(2) - range(1) + 1], &
        source=w(range(1):range(2)))

end function crop_1d_

function crop_2d_(w, range) result(wt)

    TT, dimension(:, :), intent(in) :: w
    integer, dimension(1:4), intent(in) :: range
    TT, allocatable, dimension(:, :) :: wt

    call assert(size(range) == 4 &
        .and. range(1) >= lbound(w, 1) .and. range(2) <= ubound(w, 1) &
        .and. range(3) >= lbound(w, 2) .and. range(4) <= ubound(w, 2), &
        ' <crop_array_2d> Error: range specification is wrong. ')

    call alloc_array(wt, [ &
        1, range(2) - range(1) + 1, &
        1, range(4) - range(3) + 1], &
        source=w( &
        range(1):range(2), &
        range(3):range(4)))

end function crop_2d_

function crop_3d_(w, range) result(wt)

    TT, dimension(:, :, :), intent(in) :: w
    integer, dimension(1:6), intent(in) :: range
    TT, allocatable, dimension(:, :, :) :: wt

    call assert(size(range) == 6 &
        .and. range(1) >= lbound(w, 1) .and. range(2) <= ubound(w, 1) &
        .and. range(3) >= lbound(w, 2) .and. range(4) <= ubound(w, 2) &
        .and. range(5) >= lbound(w, 3) .and. range(6) <= ubound(w, 3), &
        ' <crop_array_3d> Error: range specification is wrong. ')

    call alloc_array(wt, [ &
        1, range(2) - range(1) + 1, &
        1, range(4) - range(3) + 1, &
        1, range(6) - range(5) + 1], &
        source=w( &
        range(1):range(2), &
        range(3):range(4), &
        range(5):range(6)))

end function crop_3d_

function crop_4d_(w, range) result(wt)

    TT, dimension(:, :, :, :), intent(in) :: w
    integer, dimension(1:8), intent(in) :: range
    TT, allocatable, dimension(:, :, :, :) :: wt

    call assert(size(range) == 8 &
        .and. range(1) >= lbound(w, 1) .and. range(2) <= ubound(w, 1) &
        .and. range(3) >= lbound(w, 2) .and. range(4) <= ubound(w, 2) &
        .and. range(5) >= lbound(w, 3) .and. range(6) <= ubound(w, 3) &
        .and. range(7) >= lbound(w, 4) .and. range(8) <= ubound(w, 4), &
        ' <crop_array_4d> Error: range specification is wrong. ')

    call alloc_array(wt, [ &
        1, range(2) - range(1) + 1, &
        1, range(4) - range(3) + 1, &
        1, range(6) - range(5) + 1, &
        1, range(8) - range(7) + 1], &
        source=w( &
        range(1):range(2), &
        range(3):range(4), &
        range(5):range(6), &
        range(7):range(8)))

end function crop_4d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef crop_1d_
#undef crop_2d_
#undef crop_3d_
#undef crop_4d_
#undef crop_array_1d_
#undef crop_array_2d_
#undef crop_array_3d_
#undef crop_array_4d_
