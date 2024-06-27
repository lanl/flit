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

#define rot90_1d_      CONCAT(rot90_1d, T)
#define rot90_2d_      CONCAT(rot90_2d, T)
#define rot90_3d_      CONCAT(rot90_3d, T)

function rot90_1d_(w) result(wr)

    TT, dimension(:) :: w
    TT, allocatable, dimension(:, :) :: wr

    integer :: n

    n = size(w)

    allocate(wr(n, 1))
    wr(:, 1) = w(:)

end function rot90_1d_

function rot90_2d_(w, cw) result(wr)

    TT, dimension(:, :), intent(in) :: w
    integer, intent(in), optional :: cw

    TT, allocatable, dimension(:, :) :: wr

    integer :: n1, n2, i, j
    integer :: rotate_direction

    if (present(cw)) then
        rotate_direction = cw
    else
        rotate_direction = 1
    end if
    rotate_direction = mod(rotate_direction, 4)

    n1 = size(w, 1)
    n2 = size(w, 2)

    ! Allocate memory
    if (abs(rotate_direction) == 1 .or. abs(rotate_direction) == 3) then
        allocate (wr(1:n2, 1:n1))
    else
        allocate (wr(1:n1, 1:n2))
    end if

    ! If no rotation eventually
    if (rotate_direction == 0) then
        wr = w
        return
    end if

    ! Otherwise
    select case (rotate_direction)
        case (1, -3)
            ! Clockwise by 90 degree or equivalently counter-clockwise by 270 degree
            do j = 1, n2
                do i = 1, n1
                    wr(j, n1 - i + 1) = w(i, j)
                end do
            end do

        case (2, -2)
            ! Clockwise by 180 degree or equivalently counter-clockwise by 180 degree
            do j = 1, n2
                do i = 1, n1
                    wr(n1 - i + 1, n2 - j + 1) = w(i, j)
                end do
            end do

        case (3, -1)
            ! Clockwise by 270 degree or equivalently counter-clockwise by 90 degree
            do j = 1, n2
                do i = 1, n1
                    wr(n2 - j + 1, i) = w(i, j)
                end do
            end do

    end select

end function rot90_2d_

function rot90_3d_(w, cw, dim) result(wr)

    TT, dimension(:, :, :), intent(in) :: w
    integer, intent(in), optional :: cw, dim

    TT, allocatable, dimension(:, :, :) :: wr

    integer :: n1, n2, n3, i, j, k
    integer :: m1, m2, m3
    integer :: rotate_direction, rotate_axis

    if (present(cw)) then
        rotate_direction = cw
    else
        rotate_direction = 1
    end if
    rotate_direction = mod(rotate_direction, 4)

    if (present(dim)) then
        call assert(dim == 1 .or. dim == 2 .or. dim == 3, ' dim /= 1, 2, 3')
        rotate_axis = dim
    else
        rotate_axis = 1
    end if

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)

    ! Allocate memory
    if (abs(rotate_direction) == 1 .or. abs(rotate_direction) == 3) then
        select case (dim)
            case (1)
                m1 = n1
                m2 = n3
                m3 = n2
            case (2)
                m1 = n3
                m2 = n2
                m3 = n1
            case (3)
                m1 = n2
                m2 = n1
                m3 = n3
        end select
    else
        m1 = n1
        m2 = n2
        m3 = n3
    end if

    allocate (wr(1:m1, 1:m2, 1:m3))

    select case (dim)
        case (1)
            !$omp parallel do private(i)
            do i = 1, m1
                wr(i, :, :) = rot90_2d_(w(i, :, :), cw=rotate_direction)
            end do
            !$omp end parallel do
        case (2)
            !$omp parallel do private(j)
            do j = 1, m2
                wr(:, j, :) = rot90_2d_(w(:, j, :), cw=rotate_direction)
            end do
            !$omp end parallel do
        case (3)
            !$omp parallel do private(k)
            do k = 1, m3
                wr(:, :, k) = rot90_2d_(w(:, :, k), cw=rotate_direction)
            end do
            !$omp end parallel do
    end select

end function rot90_3d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef rot90_1d_
#undef rot90_2d_
#undef rot90_3d_
