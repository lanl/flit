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

#define rov_1d_     CONCAT(rov_1d, T)
#define rov_2d_     CONCAT(rov_2d, T)
#define rov_3d_     CONCAT(rov_3d, T)
#define rov_4d_     CONCAT(rov_4d, T)
#define rov_axis_2d_        CONCAT(rov_axis_2d, T)
#define rov_axis_3d_        CONCAT(rov_axis_3d, T)

function rov_1d_(w) result(r)

    TT, dimension(:), intent(in) :: w
    TT :: r

    r = maxval(w) - minval(w)

end function rov_1d_

function rov_2d_(w) result(r)

    TT, dimension(:, :), intent(in) :: w
    TT :: r

    r = maxval(w) - minval(w)

end function rov_2d_

function rov_3d_(w) result(r)

    TT, dimension(:, :, :), intent(in) :: w
    TT :: r

    r = maxval(w) - minval(w)

end function rov_3d_

function rov_4d_(w) result(r)

    TT, dimension(:, :, :, :), intent(in) :: w
    TT :: r

    r = maxval(w) - minval(w)

end function rov_4d_

function rov_axis_2d_(w, axis) result(r)

    TT, dimension(:, :), intent(in) :: w
    integer, intent(in) :: axis
    TT, allocatable, dimension(:) :: r

    integer :: n1, n2
    integer :: i

    n1 = size(w, 1)
    n2 = size(w, 2)

    select case(axis)
        case(1)
            r = const(n2, 0)
            do i = 1, n2
                r(i) = rov_1d_(w(:, i))
            end do
        case(2)
            r = const(n1, 0)
            do i = 1, n1
                r(i) = rov_1d_(w(i, :))
            end do
    end select

end function rov_axis_2d_

function rov_axis_3d_(w, axis) result(r)

    TT, dimension(:, :, :), intent(in) :: w
    integer, intent(in) :: axis
    TT, allocatable, dimension(:, :) :: r

    integer :: n1, n2, n3
    integer :: i, j

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)

    select case(axis)
        case(1)
            r = const(n2, n3, 0)
            do j = 1, n3
                do i = 1, n2
                    r(i, j) = rov_1d_(w(:, i, j))
                end do
            end do
        case(2)
            r = const(n1, n3, 0)
            do j = 1, n3
                do i = 1, n1
                    r(i, j) = rov_1d_(w(i, :, j))
                end do
            end do
        case(3)
            r = const(n1, n2, 0)
            do j = 1, n2
                do i = 1, n1
                    r(i, j) = rov_1d_(w(i, j, :))
                end do
            end do
    end select

end function rov_axis_3d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef rov_1d_
#undef rov_2d_
#undef rov_3d_
#undef rov_4d_
#undef rov_axis_2d_
#undef rov_axis_3d_

