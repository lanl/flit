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

#define pad_array_1d_     CONCAT(pad_array_1d, T)
#define pad_array_2d_     CONCAT(pad_array_2d, T)
#define pad_array_3d_     CONCAT(pad_array_3d, T)
!#define pad_array_4d_     CONCAT(pad_array_4d, T)
!
#define pad_1d_     CONCAT(pad_1d, T)
#define pad_2d_     CONCAT(pad_2d, T)
#define pad_3d_     CONCAT(pad_3d, T)
!#define pad_4d_     CONCAT(pad_4d, T)

subroutine pad_array_1d_(w, pad, method, const)

    ! arguments
    TT, allocatable, dimension(:), intent(inout) :: w
    integer, dimension(1:2), intent(in) :: pad
    character(len=*), dimension(1:2), intent(in), optional :: method
    TT, intent(in), optional :: const

    ! local variables
    integer :: i
    integer :: n1beg, n1end
    TTT, allocatable, dimension(:) :: wp
    integer :: l1, u1
    character(len=16), dimension(1:2) :: pad_method
    TTT :: pad_const

    if (present(const)) then
        pad_const = const
    else
        pad_const = DEFAULT_VALUE
    end if

    ! bounds
    n1beg = lbound(w, 1)
    n1end = ubound(w, 1)

    l1 = pad(1)
    u1 = pad(2)

    ! new array
    allocate(wp(n1beg - l1:n1end + u1))
    do i = n1beg, n1end
        wp(i) = w(i)
    end do

    ! padding method
    if (present(method)) then
        pad_method = method
    else
        pad_method = ['edge', 'edge']
    end if

    do i = 1, size(pad_method)
        call assert(any(pad_method(i) == ['edge', 'symm', 'const']), &
            ' <pad_1d> Error: Padding method must be one of edge, symm, const.')
    end do

    ! axis 1 lower boundary
    select case (pad_method(1))
        case ('edge')
            do i = 1, l1
                wp(n1beg - i) = wp(n1beg)
            end do
        case ('symm')
            do i = 1, l1
                wp(n1beg - i) = wp(n1beg + i - 1)
            end do
        case ('const')
            do i = 1, l1
                wp(n1beg - i) = pad_const
            end do
    end select

    ! axis 1 upper boundary
    select case (pad_method(2))
        case ('edge')
            do i = 1, u1
                wp(n1end + i) = wp(n1end)
            end do
        case ('symm')
            do i = 1, u1
                wp(n1end + i) = wp(n1end - i + 1)
            end do
        case ('const')
            do i = 1, u1
                wp(n1end + i) = pad_const
            end do
    end select

    ! reallocate input array
    w = wp

end subroutine pad_array_1d_

subroutine pad_array_2d_(w, pad, method, const)

    ! arguments
    TT, allocatable, dimension(:, :), intent(inout) :: w
    integer, dimension(1:4), intent(in) :: pad
    character(len=*), dimension(1:4), intent(in), optional :: method
    TT, intent(in), optional :: const

    ! local variables
    integer :: i, j
    integer :: n1beg, n1end, n2beg, n2end
    TTT, allocatable, dimension(:, :) :: wp
    integer :: l1, u1, l2, u2
    character(len=16), dimension(1:4) :: pad_method
    TTT :: pad_const

    if (present(const)) then
        pad_const = const
    else
        pad_const = DEFAULT_VALUE
    end if

    ! bounds
    n1beg = lbound(w, 1)
    n1end = ubound(w, 1)
    n2beg = lbound(w, 2)
    n2end = ubound(w, 2)

    l1 = pad(1)
    u1 = pad(2)
    l2 = pad(3)
    u2 = pad(4)

    ! new array
    allocate(wp(n1beg - l1:n1end + u1, n2beg - l2:n2end + u2))
    !$omp parallel do private(i, j) collapse(2)
    do j = n2beg, n2end
        do i = n1beg, n1end
            wp(i, j) = w(i, j)
        end do
    end do
    !$omp end parallel do

    ! padding method
    if (present(method)) then
        pad_method = method
    else
        pad_method = ['edge', 'edge', 'edge', 'edge']
    end if

    do i = 1, size(pad_method)
        call assert(any(pad_method(i) == ['edge', 'symm', 'const']), &
            ' <pad_2d> Error: Padding method must be one of edge, symm, const.')
    end do

    ! axis 1 lower boundary
    select case (pad_method(1))
        case ('edge')
            !$omp parallel do private(i, j) collapse(2)
            do j = n2beg - l2, n2end + u2
                do i = 1, l1
                    wp(n1beg - i, j) = wp(n1beg, j)
                end do
            end do
            !$omp end parallel do
        case ('symm')
            !$omp parallel do private(i, j) collapse(2)
            do j = n2beg - l2, n2end + u2
                do i = 1, l1
                    wp(n1beg - i, j) = wp(n1beg + i - 1, j)
                end do
            end do
            !$omp end parallel do
        case ('const')
            !$omp parallel do private(i, j) collapse(2)
            do j = n2beg - l2, n2end + u2
                do i = 1, l1
                    wp(n1beg - i, j) = pad_const
                end do
            end do
            !$omp end parallel do
    end select

    ! axis 1 upper boundary
    select case (pad_method(2))
        case ('edge')
            !$omp parallel do private(i, j) collapse(2)
            do j = n2beg - l2, n2end + u2
                do i = 1, u1
                    wp(n1end + i, j) = wp(n1end, j)
                end do
            end do
            !$omp end parallel do
        case ('symm')
            !$omp parallel do private(i, j) collapse(2)
            do j = n2beg - l2, n2end + u2
                do i = 1, u1
                    wp(n1end + i, j) = wp(n1end - i + 1, j)
                end do
            end do
            !$omp end parallel do
        case ('const')
            !$omp parallel do private(i, j) collapse(2)
            do j = n2beg - l2, n2end + u2
                do i = 1, u1
                    wp(n1end + i, j) = pad_const
                end do
            end do
            !$omp end parallel do
    end select

    ! axis 2 lower boundary
    select case (pad_method(3))
        case ('edge')
            !$omp parallel do private(i, j) collapse(2)
            do j = 1, l2
                do i = n1beg - l1, n1end + u1
                    wp(i, n2beg - j) = wp(i, n2beg)
                end do
            end do
            !$omp end parallel do
        case ('symm')
            !$omp parallel do private(i, j) collapse(2)
            do j = 1, l2
                do i = n1beg - l1, n1end + u1
                    wp(i, n2beg - j) = wp(i, n2beg + j - 1)
                end do
            end do
            !$omp end parallel do
        case ('const')
            !$omp parallel do private(i, j) collapse(2)
            do j = 1, l2
                do i = n1beg - l1, n1end + u1
                    wp(i, n2beg - j) = pad_const
                end do
            end do
            !$omp end parallel do
    end select

    ! axis 2 upper boundary
    select case (pad_method(4))
        case ('edge')
            !$omp parallel do private(i, j) collapse(2)
            do j = 1, u2
                do i = n1beg - l1, n1end + u1
                    wp(i, n2end + j) = wp(i, n2end)
                end do
            end do
            !$omp end parallel do
        case ('symm')
            !$omp parallel do private(i, j) collapse(2)
            do j = 1, u2
                do i = n1beg - l1, n1end + u1
                    wp(i, n2end + j) = wp(i, n2end - j + 1)
                end do
            end do
            !$omp end parallel do
        case ('const')
            !$omp parallel do private(i, j) collapse(2)
            do j = 1, u2
                do i = n1beg - l1, n1end + u1
                    wp(i, n2end + j) = pad_const
                end do
            end do
            !$omp end parallel do
    end select

    ! reallocate input array
    w = wp

end subroutine pad_array_2d_

subroutine pad_array_3d_(w, pad, method, const)

    ! arguments
    TT, allocatable, dimension(:, :, :), intent(inout) :: w
    integer, dimension(1:6), intent(in) :: pad
    character(len=*), dimension(1:6), intent(in), optional :: method
    TT, intent(in), optional :: const

    ! local variables
    integer :: i, j, k
    integer :: n1beg, n1end, n2beg, n2end, n3beg, n3end
    integer :: l1, u1, l2, u2, l3, u3
    TTT, allocatable, dimension(:, :, :) :: wp
    character(len=16), dimension(1:6) :: pad_method
    TTT :: pad_const

    if (present(const)) then
        pad_const = const
    else
        pad_const = DEFAULT_VALUE
    end if

    ! bounds
    n1beg = lbound(w, 1)
    n1end = ubound(w, 1)
    n2beg = lbound(w, 2)
    n2end = ubound(w, 2)
    n3beg = lbound(w, 3)
    n3end = ubound(w, 3)

    l1 = pad(1)
    u1 = pad(2)
    l2 = pad(3)
    u2 = pad(4)
    l3 = pad(5)
    u3 = pad(6)

    ! padding method
    if (present(method)) then
        pad_method = method
    else
        pad_method = ['edge', 'edge', 'edge', 'edge', 'edge', 'edge']
    end if

    do i = 1, size(pad_method)
        call assert(any(pad_method(i) == ['edge', 'symm', 'const']), &
            ' <pad_3d> Error: Padding method must be one of edge, symm, const.')
    end do

    ! new array
    allocate(wp(n1beg - l1:n1end + u1, n2beg - l2:n2end + u2, n3beg - l3:n3end + u3))
    !$omp parallel do private(i, j, k) collapse(3)
    do k = n3beg, n3end
        do j = n2beg, n2end
            do i = n1beg, n1end
                wp(i, j, k) = w(i, j, k)
            end do
        end do
    end do
    !$omp end parallel do

    ! Left boundary
    select case (pad_method(1))
        case ('edge')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = n3beg - l3, n3end + u3
                do j = n2beg - l2, n2end + u2
                    do i = 1, l1
                        wp(n1beg - i, j, k) = wp(n1beg, j, k)
                    end do
                end do
            end do
            !$omp end parallel do
        case ('symm')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = n3beg - l3, n3end + u3
                do j = n2beg - l2, n2end + u2
                    do i = 1, l1
                        wp(n1beg - i, j, k) = wp(n1beg + i - 1, j, k)
                    end do
                end do
            end do
            !$omp end parallel do
        case ('const')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = n3beg - l3, n3end + u3
                do j = n2beg - l2, n2end + u2
                    do i = 1, l1
                        wp(n1beg - i, j, k) = pad_const
                    end do
                end do
            end do
            !$omp end parallel do
    end select

    ! Right boundary
    select case (pad_method(2))
        case ('edge')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = n3beg - l3, n3end + u3
                do j = n2beg - l2, n2end + u2
                    do i = 1, u1
                        wp(n1end + i, j, k) = wp(n1end, j, k)
                    end do
                end do
            end do
            !$omp end parallel do
        case ('symm')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = n3beg - l3, n3end + u3
                do j = n2beg - l2, n2end + u2
                    do i = 1, u1
                        wp(n1end + i, j, k) = wp(n1end - i + 1, j, k)
                    end do
                end do
            end do
            !$omp end parallel do
        case ('const')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = n3beg - l3, n3end + u3
                do j = n2beg - l2, n2end + u2
                    do i = 1, u1
                        wp(n1end + i, j, k) = pad_const
                    end do
                end do
            end do
            !$omp end parallel do
    end select

    ! Front boundary
    select case (pad_method(3))
        case ('edge')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = n3beg - l3, n3end + u3
                do j = 1, l2
                    do i = n1beg - l1, n1end + u1
                        wp(i, n2beg - j, k) = wp(i, n2beg, k)
                    end do
                end do
            end do
            !$omp end parallel do
        case ('symm')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = n3beg - l3, n3end + u3
                do j = 1, l2
                    do i = n1beg - l1, n1end + u1
                        wp(i, n2beg - j, k) = wp(i, n2beg + j - 1, k)
                    end do
                end do
            end do
            !$omp end parallel do
        case ('const')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = n3beg - l3, n3end + u3
                do j = 1, l2
                    do i = n1beg - l1, n1end + u1
                        wp(i, n2beg - j, k) = pad_const
                    end do
                end do
            end do
            !$omp end parallel do
    end select

    ! Back boundary
    select case (pad_method(4))
        case ('edge')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = n3beg - l3, n3end + u3
                do j = 1, u2
                    do i = n1beg - l1, n1end + u1
                        wp(i, n2end + j, k) = wp(i, n2end, k)
                    end do
                end do
            end do
            !$omp end parallel do
        case ('symm')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = n3beg - l3, n3end + u3
                do j = 1, u2
                    do i = n1beg - l1, n1end + u1
                        wp(i, n2end + j, k) = wp(i, n2end - j + 1, k)
                    end do
                end do
            end do
            !$omp end parallel do
        case ('const')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = n3beg - l3, n3end + u3
                do j = 1, u2
                    do i = n1beg - l1, n1end + u1
                        wp(i, n2end + j, k) = pad_const
                    end do
                end do
            end do
            !$omp end parallel do
    end select

    ! Top boundary
    select case (pad_method(5))
        case ('edge')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = 1, l3
                do j = n2beg - l2, n2end + u2
                    do i = n1beg - l1, n1end + u1
                        wp(i, j, n3beg - k) = wp(i, j, n3beg)
                    end do
                end do
            end do
            !$omp end parallel do
        case ('symm')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = 1, l3
                do j = n2beg - l2, n2end + u2
                    do i = n1beg - l1, n1end + u1
                        wp(i, j, n3beg - k) = wp(i, j, n3beg + k - 1)
                    end do
                end do
            end do
            !$omp end parallel do
        case ('const')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = 1, l3
                do j = n2beg - l2, n2end + u2
                    do i = n1beg - l1, n1end + u1
                        wp(i, j, n3beg - k) = pad_const
                    end do
                end do
            end do
            !$omp end parallel do
    end select

    ! Bottom boundary
    select case (pad_method(6))
        case ('edge')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = 1, u3
                do j = n2beg - l2, n2end + u2
                    do i = n1beg - l1, n1end + u1
                        wp(i, j, n3end + k) = wp(i, j, n3end)
                    end do
                end do
            end do
            !$omp end parallel do
        case ('symm')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = 1, u3
                do j = n2beg - l2, n2end + u2
                    do i = n1beg - l1, n1end + u1
                        wp(i, j, n3end + k) = wp(i, j, n3end - k + 1)
                    end do
                end do
            end do
            !$omp end parallel do
        case ('const')
            !$omp parallel do private(i, j, k) collapse(3)
            do k = 1, u3
                do j = n2beg - l2, n2end + u2
                    do i = n1beg - l1, n1end + u1
                        wp(i, j, n3end + k) = pad_const
                    end do
                end do
            end do
            !$omp end parallel do
    end select

    ! reallocate array
    w = wp

end subroutine pad_array_3d_

function pad_1d_(w, pad, method, const) result(wp)

    ! arguments
    TT, dimension(:) :: w
    integer, dimension(1:2) :: pad
    character(len=*), dimension(1:2), optional :: method
    TT, optional :: const
    TTT, allocatable, dimension(:) :: wp

    ! local variables
    character(len=16), dimension(1:2) :: pad_method
    TTT :: pad_const

    if (present(const)) then
        pad_const = const
    else
        pad_const = DEFAULT_VALUE
    end if

    ! padding method
    if (present(method)) then
        pad_method = method
    else
        pad_method = ['edge', 'edge']
    end if

    wp = w
    call pad_array_1d_(wp, pad, pad_method, pad_const)

end function pad_1d_

function pad_2d_(w, pad, method, const) result(wp)

    ! arguments
    TT, dimension(:, :) :: w
    integer, dimension(1:4) :: pad
    character(len=*), dimension(1:4), optional :: method
    TT, optional :: const
    TTT, allocatable, dimension(:, :) :: wp

    ! local variables
    character(len=16), dimension(1:4) :: pad_method
    TTT :: pad_const

    if (present(const)) then
        pad_const = const
    else
        pad_const = DEFAULT_VALUE
    end if

    ! padding method
    if (present(method)) then
        pad_method = method
    else
        pad_method = ['edge', 'edge', 'edge', 'edge']
    end if

    wp = w
    call pad_array_2d_(wp, pad, pad_method, pad_const)

end function pad_2d_

function pad_3d_(w, pad, method, const) result(wp)

    ! arguments
    TT, dimension(:, :, :) :: w
    integer, dimension(1:6) :: pad
    character(len=*), dimension(1:6), optional :: method
    TT, optional :: const
    TTT, allocatable, dimension(:, :, :) :: wp

    ! local variables
    character(len=16), dimension(1:6) :: pad_method
    TTT :: pad_const

    if (present(const)) then
        pad_const = const
    else
        pad_const = DEFAULT_VALUE
    end if

    ! padding method
    if (present(method)) then
        pad_method = method
    else
        pad_method = ['edge', 'edge', 'edge', 'edge', 'edge', 'edge']
    end if

    wp = w
    call pad_array_3d_(wp, pad, pad_method, pad_const)

end function pad_3d_

#undef T
#undef TT
#undef TTT
#undef DEFAULT_VALUE

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef pad_array_1d_
#undef pad_array_2d_
#undef pad_array_3d_

#undef pad_1d_
#undef pad_2d_
#undef pad_3d_

