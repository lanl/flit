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

#define laplace_filt_1d_      CONCAT(laplace_filt_1d, T)
#define laplace_filt_2d_      CONCAT(laplace_filt_2d, T)
#define laplace_filt_3d_      CONCAT(laplace_filt_3d, T)

!
!> 1D discrete Laplacian filter
!
function laplace_filt_1d_(w, type) result(wt)

    ! Arguments
    TT, dimension(:), intent(inout) :: w
    integer, intent(in), optional :: type

    ! Local variables
    integer :: i, n1
    TT, allocatable, dimension(:) :: wt, ww
    integer :: filter_type

    if (present(type)) then
        filter_type = type
    else
        filter_type = 1
    end if

    ! dimension
    n1 = size(w)

    ! Allocate memory
    wt = zeros(n1)
    ww = w
    call pad_array(ww, [1, 1])

    ! Finte-difference Laplacian operator
    select case (filter_type)

        case (1)
            !$omp parallel do private(i)
            do i = 1, n1
                wt(i) = ww(i - 1) + ww(i + 1) - 2.0*ww(i)
            end do
            !$omp end parallel do

        case (2)
            !$omp parallel do private(i)
            do i = 1, n1
                wt(i) = sum(ww(i - 1:i + 1)) - 5.0*ww(i)
            end do
            !$omp end parallel do

    end select

end function laplace_filt_1d_

!
!> 2D discrete Laplacian filter
!
function laplace_filt_2d_(w, type) result(wt)

    ! Arguments
    TT, dimension(:, :), intent(inout) :: w
    integer, intent(in), optional :: type

    ! Local variables
    integer :: i, j, n1, n2
    TT, allocatable, dimension(:, :) :: wt, ww
    integer :: filter_type

    if (present(type)) then
        filter_type = type
    else
        filter_type = 1
    end if

    ! dimension
    n1 = size(w, 1)
    n2 = size(w, 2)

    ! Allocate memory
    wt = zeros(n1, n2)
    ww = w
    call pad_array(ww, [1, 1, 1, 1])

    ! Finte-difference Laplacian operator
    select case (filter_type)

        case (1)
            !$omp parallel do private(i, j)
            do j = 1, n2
                do i = 1, n1
                    wt(i, j) = ww(i - 1, j) + ww(i + 1, j) &
                        + ww(i, j - 1) + ww(i, j + 1) &
                        - 4.0*ww(i, j)
                end do
            end do
            !$omp end parallel do

        case (2)
            !$omp parallel do private(i, j)
            do j = 1, n2
                do i = 1, n1
                    wt(i, j) = sum(ww(i - 1:i + 1, j - 1:j + 1)) &
                        - 9.0*ww(i, j)
                end do
            end do
            !$omp end parallel do

    end select

    deallocate (ww)

end function laplace_filt_2d_

!
!> 3D discrete Laplacian filter
!
function laplace_filt_3d_(w, type) result(wt)

    ! Arguments
    TT, dimension(:, :, :), intent(inout) :: w
    integer, intent(in), optional :: type

    ! Local variables
    integer :: i, j, k, n1, n2, n3
    TT, allocatable, dimension(:, :, :) :: wt, ww
    integer :: filter_type

    if (present(type)) then
        filter_type = type
    else
        filter_type = 1
    end if

    ! dimension of arrays
    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)

    ! Allocate memory
    wt = zeros(n1, n2, n3)
    ww = w
    call pad_array(ww, [1, 1, 1, 1, 1, 1])

    ! Finte-difference Laplacian operator
    select case (filter_type)

        case (1)
            !$omp parallel do private(i, j, k)
            do k = 1, n3
                do j = 1, n2
                    do i = 1, n1
                        wt(i, j, k) = ww(i - 1, j, k) + ww(i + 1, j, k) &
                            + ww(i, j - 1, k) + ww(i, j + 1, k) &
                            + ww(i, j, k - 1) + ww(i, j, k + 1) &
                            -6.0*ww(i, j, k)
                    end do
                end do
            end do
            !$omp end parallel do

        case (2)
            !$omp parallel do private(i, j, k)
            do k = 1, n3
                do j = 1, n2
                    do i = 1, n1
                        wt(i, j, k) = sum(ww(i - 1:i + 1, j - 1:j + 1, k - 1:k + 1)) &
                            -13.0*ww(i, j, k)
                    end do
                end do
            end do
            !$omp end parallel do

    end select

end function laplace_filt_3d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef laplace_filt_1d_
#undef laplace_filt_2d_
#undef laplace_filt_3d_
