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

#define eigen_simple_l_or_r_     CONCAT(eigen_simple_l_or_r, T)
#define eigen_general_l_or_r_     CONCAT(eigen_general_l_or_r, T)
#define eigen_simple_l_and_r_     CONCAT(eigen_simple_l_and_r, T)
#define eigen_general_l_and_r_     CONCAT(eigen_general_l_and_r, T)

! ==============================================================
! A X = lambada X
! Returns eigenvalues and left or right eigenvectors
!
subroutine eigen_simple_l_or_r_(a, w, v, lr)

    TT, dimension(:, :), intent(in) :: a
    TTT, allocatable, dimension(:), intent(inout) :: w
    TTT, allocatable, dimension(:, :), intent(inout), optional :: v
    character(len=*), intent(in), optional :: lr

    character(len=12) :: left_or_right
    integer :: n
    TTT, allocatable :: ac(:, :), vl(:, :), vr(:, :)

    call assert(size(a, 1) /= size(a, 2), &
        ' <eigen_general_l_or_r> Error: a must be square. ')

    if (present(lr)) then
        left_or_right = lr
    else
        left_or_right = 'right'
    end if

    n = size(a, 1)
    ac = nTT(a)
    w = zeros(n)
    vl = zeros(n, n)
    vr = zeros(n, n)

    call geev(ac, w, vl, vr)

    if (present(v)) then
        select case (left_or_right)
            case ('right')
                v = vr
            case ('left')
                v = vl
        end select
    end if

end subroutine eigen_simple_l_or_r_

! ==============================================================
! A X = lambada B X
! Returns eigenvalues and left or right eigenvectors
!
subroutine eigen_general_l_or_r_(a, b, w, v, lr)

    TT, dimension(:, :), intent(in) :: a, b
    TTT, allocatable, dimension(:), intent(inout) :: w
    TTT, allocatable, dimension(:, :), intent(inout), optional :: v
    character(len=*), intent(in), optional :: lr

    character(len=12) :: left_or_right
    integer :: n
    TTT, allocatable :: ac(:, :), bc(:, :), ww(:), vl(:, :), vr(:, :)

    call assert(size(a, 1) /= size(a, 2) .and. size(b, 1) /= size(b, 2), &
        ' <eigen_general_l_or_r> Error: both a and b must be square. ')

    if (present(lr)) then
        left_or_right = lr
    else
        left_or_right = 'right'
    end if

    n = size(a, 1)
    ac = nTT(a)
    bc = nTT(b)
    w = zeros(n)
    ww = zeros(n)
    vl = zeros(n, n)
    vr = zeros(n, n)

    call ggev(ac, bc, w, ww, vl, vr)
    where (abs(ww) /= 0)
        w = w/ww
    end where

    if (present(v)) then
        select case (left_or_right)
            case ('right')
                v = vr
            case ('left')
                v = vl
        end select
    end if

end subroutine eigen_general_l_or_r_

! ==============================================================
! A X = lambda X
! Returns eigenvalues and both left and right eigenvectors
!
subroutine eigen_simple_l_and_r_(a, w, vl, vr)

    TT, dimension(:, :), intent(in) :: a
    TTT, allocatable, dimension(:), intent(inout) :: w
    TTT, allocatable, dimension(:, :), intent(inout) :: vl, vr

    integer :: n
    TTT, allocatable, dimension(:, :) :: ac

    call assert(size(a, 1) == size(a, 2), &
        ' <eigen_simple_l_and_r> Error: A must be square.')

    n = size(a, 1)
    ac = nTT(a)
    w = zeros(n)
    vl = zeros(n, n)
    vr = zeros(n, n)

    call geev(ac, w, vl, vr)

end subroutine eigen_simple_l_and_r_

! ==============================================================
! A X = lambda B X
! Returns eigenvalues and both left and right eigenvectors
!
subroutine eigen_general_l_and_r_(a, b, w, vl, vr)

    TT, dimension(:, :), intent(in) :: a, b
    TTT, allocatable, dimension(:), intent(inout) :: w
    TTT, allocatable, dimension(:, :), intent(inout) :: vl, vr

    integer :: n
    TTT, allocatable, dimension(:, :) :: ac, bc
    TTT, allocatable, dimension(:) :: cc

    call assert(size(a, 1) == size(a, 2) .and. size(b, 1) == size(b, 2), &
        ' <eigen_general_l_and_r> Error: Both A and B must be square.')

    n = size(a, 1)
    ac = nTT(a)
    bc = nTT(b)
    cc = zeros(n)
    w = zeros(n)
    vl = zeros(n, n)
    vr = zeros(n, n)

    call ggev(ac, bc, w, cc, vl, vr)

    where (abs(cc) /= 0)
        w = w/cc
    end where

end subroutine eigen_general_l_and_r_

#undef T
#undef TT
#undef TTT
#undef nTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef eigen_simple_l_or_r_
#undef eigen_general_l_or_r_
#undef eigen_simple_l_and_r_
#undef eigen_general_l_and_r_

