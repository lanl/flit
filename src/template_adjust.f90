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

#define adjust_1d_     CONCAT(adjust_1d, T)
#define adjust_2d_     CONCAT(adjust_2d, T)
#define adjust_3d_     CONCAT(adjust_3d, T)
#define adjust_4d_     CONCAT(adjust_4d, T)

function adjust_1d_(w, n) result(wp)

    TT, dimension(:), intent(in) :: w
    integer, intent(in) :: n
    TT, allocatable, dimension(:) :: wp

    integer :: nw

    nw = size(w)
    allocate (wp(1:n))
    wp = 0
    wp(1:min(nw, n)) = w(1:min(nw, n))

end function adjust_1d_

function adjust_2d_(w, n) result(wp)

    TT, dimension(:, :), intent(in) :: w
    integer, dimension(:), intent(in) :: n
    TT, allocatable, dimension(:, :) :: wp

    integer :: nw(1:2)

    nw = min(shape(w), n)
    allocate (wp(1:n(1), 1:n(2)))
    wp = 0
    wp(1:nw(1), 1:nw(2)) = w(1:nw(1), 1:nw(2))

end function adjust_2d_

function adjust_3d_(w, n) result(wp)

    TT, dimension(:, :, :), intent(in) :: w
    integer, dimension(:), intent(in) :: n
    TT, allocatable, dimension(:, :, :) :: wp

    integer :: nw(1:3)

    nw = min(shape(w), n)
    allocate (wp(1:n(1), 1:n(2), 1:n(3)))
    wp = 0
    wp(1:nw(1), 1:nw(2), 1:nw(3)) = w(1:nw(1), 1:nw(2), 1:nw(3))

end function adjust_3d_

function adjust_4d_(w, n) result(wp)

    TT, dimension(:, :, :, :), intent(in) :: w
    integer, dimension(:), intent(in) :: n
    TT, allocatable, dimension(:, :, :, :) :: wp

    integer :: nw(1:4)

    nw = min(shape(w), n)
    allocate (wp(1:n(1), 1:n(2), 1:n(3), 1:n(4)))
    wp = 0
    wp(1:nw(1), 1:nw(2), 1:nw(3), 1:nw(4)) = w(1:nw(1), 1:nw(2), 1:nw(3), 1:nw(4))

end function adjust_4d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef adjust_1d_
#undef adjust_2d_
#undef adjust_3d_
#undef adjust_4d_
