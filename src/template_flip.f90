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

#define flip_1d_      CONCAT(flip_1d, T)
#define flip_2d_      CONCAT(flip_2d, T)
#define flip_3d_      CONCAT(flip_3d, T)

function flip_1d_(w) result(wf)

    TT, dimension(:), intent(in) :: w

    TT, allocatable, dimension(:) :: wf
    integer :: n1

    n1 = size(w, 1)
    allocate (wf(1:n1))

    wf = w(n1:1:-1)

end function flip_1d_

function flip_2d_(w, axis) result(wf)

    TT, dimension(:, :), intent(in) :: w
    integer, dimension(:), intent(in) :: axis

    TT, allocatable, dimension(:, :) :: wf
    integer :: n1, n2, nf
    integer, allocatable, dimension(:, :) :: f
    integer :: i, r(1:3)

    n1 = size(w, 1)
    n2 = size(w, 2)
    nf = size(axis)

    allocate (wf(1:n1, 1:n2))
    allocate (f(1:2, 1:3))
    f(1, :) = [1, n1, 1]
    f(2, :) = [1, n2, 1]

    do i = 1, size(axis)
        select case (axis(i))
            case (1)
                r = f(1, :)
                f(1, :) = [r(2), r(1), -r(3)]
            case (2)
                r = f(2, :)
                f(2, :) = [r(2), r(1), -r(3)]
        end select
    end do

    wf = w(f(1, 1):f(1, 2):f(1, 3), f(2, 1):f(2, 2):f(2, 3))

end function flip_2d_

function flip_3d_(w, axis) result(wf)

    TT, dimension(:, :, :), intent(in) :: w
    integer, dimension(:), intent(in) :: axis

    TT, allocatable, dimension(:, :, :) :: wf
    integer :: n1, n2, n3, nf
    integer, allocatable, dimension(:, :) :: f
    integer :: i, r(1:3)

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)
    nf = size(axis)

    allocate (wf(1:n1, 1:n2, 1:n3))
    allocate (f(1:3, 1:3))
    f(1, :) = [1, n1, 1]
    f(2, :) = [1, n2, 1]
    f(3, :) = [1, n3, 1]

    do i = 1, size(axis)
        select case (axis(i))
            case (1)
                r = f(1, :)
                f(1, :) = [r(2), r(1), -r(3)]
            case (2)
                r = f(2, :)
                f(2, :) = [r(2), r(1), -r(3)]
            case (3)
                r = f(3, :)
                f(3, :) = [r(2), r(1), -r(3)]
        end select
    end do

    wf = w(f(1, 1):f(1, 2):f(1, 3), &
        f(2, 1):f(2, 2):f(2, 3), &
        f(3, 1):f(3, 2):f(3, 3))

end function flip_3d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef flip_1d_
#undef flip_2d_
#undef flip_3d_
