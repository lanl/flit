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

#define balance_filt_1d_      CONCAT(balance_filt_1d, T)
#define balance_filt_2d_      CONCAT(balance_filt_2d, T)
#define balance_filt_3d_      CONCAT(balance_filt_3d, T)

!
!> 1D moving balance by dividing energy
!
function balance_filt_1d_(w, radius, eps) result(wp)

    TT, dimension(:) :: w
    integer :: radius
    real, optional :: eps
    TT, allocatable, dimension(:) :: wp

    integer :: n1, i, od
    TT, allocatable, dimension(:) :: ww, wr
    real :: epslevel, scalar

    if (present(eps)) then
        epslevel = eps
    else
        epslevel = 0.01
    end if

    od = order_of_magnitude(maxval(abs(w)))

    n1 = size(w)
    ww = zeros(n1)
    wp = (w/10.0d0**od)**2
    if (maxval(abs(w)) == 0) then
        wp = 0.0
        return
    end if
    call pad_array(wp, [radius, radius], ['symm', 'symm'])
    call alloc_array(wr, [-radius, radius])

    do i = 1, n1
        wr = wp(i - radius:i + radius)
        scalar = epslevel*maxval(abs(wr))
        if (scalar /= 0) then
            ww(i) = sum(wr) + scalar
        else
            ww(i) = float_tiny*(radius + 1)*2
        end if
    end do
    ww = sqrt(ww/(2*radius + 1))

    wp = TTT((w/10.0d0**od)/ww)

end function balance_filt_1d_

!
!> 2D moving balance by dividing energy
!
function balance_filt_2d_(w, radius, eps) result(wp)

    TT, dimension(:, :) :: w
    integer, dimension(:) :: radius
    real, optional :: eps
    TT, allocatable, dimension(:, :) :: wp

    integer :: n1, n2, i, j, od
    TT, allocatable, dimension(:, :) :: ww, wr
    real :: epslevel, scalar

    if (present(eps)) then
        epslevel = eps
    else
        epslevel = 0.01
    end if

    od = order_of_magnitude(maxval(abs(w)))

    n1 = size(w, 1)
    n2 = size(w, 2)
    ww = zeros(n1, n2)
    wp = (w/10.0d0**od)**2
    if (maxval(abs(w)) == 0) then
        wp = 0.0
        return
    end if
    call pad_array(wp, [radius(1), radius(1), radius(2), radius(2)], &
        ['symm', 'symm', 'symm', 'symm'])
    call alloc_array(wr, [-radius(1), radius(1), -radius(2), radius(2)])

    !$omp parallel do private(i, j, wr, scalar) collapse(2)
    do j = 1, n2
        do i = 1, n1
            wr = wp(i - radius(1):i + radius(1), j - radius(2):j + radius(2))
            scalar = epslevel*maxval(abs(wr))
            if (scalar /= 0) then
                ww(i, j) = sum(wr) + scalar
            else
                ww(i, j) = float_tiny*product(radius + 1)*4
            end if
        end do
    end do
    !$omp end parallel do
    ww = sqrt(ww/(2*radius(1) + 1)/(2*radius(2) + 1))

    wp = TTT((w/10.0d0**od)/ww)

end function balance_filt_2d_

!
!> 3D moving balance by dividing energy
!
function balance_filt_3d_(w, radius, eps) result(wp)

    TT, dimension(:, :, :) :: w
    integer, dimension(:) :: radius
    real, optional :: eps
    TT, allocatable, dimension(:, :, :) :: wp

    integer :: n1, n2, n3, i, j, k, od
    TT, allocatable, dimension(:, :, :) :: ww, wr
    real :: epslevel, scalar

    if (present(eps)) then
        epslevel = eps
    else
        epslevel = 0.01
    end if

    od = order_of_magnitude(maxval(abs(w)))

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)
    ww = zeros(n1, n2, n3)
    wp = (w/10.0d0**od)**2
    if (maxval(abs(w)) == 0) then
        wp = 0.0
        return
    end if
    call pad_array(wp, [radius(1), radius(1), radius(2), radius(2), radius(3), radius(3)], &
        ['symm', 'symm', 'symm', 'symm', 'symm', 'symm'])
    call alloc_array(wr, [-radius(1), radius(1), -radius(2), radius(2), -radius(3), radius(3)])

    !$omp parallel do private(i, j, k, wr, scalar) collapse(3)
    do k = 1, n3
        do j = 1, n2
            do i = 1, n1
                wr = wp(i - radius(1):i + radius(1), &
                    j - radius(2):j + radius(2), &
                    k - radius(3):k + radius(3))
                scalar = epslevel*maxval(abs(wr))
                if (scalar /= 0) then
                    ww(i, j, k) = sum(wr) + scalar
                else
                    ww(i, j, k) = float_tiny*product(radius + 1)*8
                end if
            end do
        end do
    end do
    !$omp end parallel do
    ww = sqrt(ww/(2*radius(1) + 1)/(2*radius(2) + 1)/(2*radius(3) + 1))

    wp = TTT((w/10.0d0**od)/ww)

end function balance_filt_3d_

#undef T
#undef TT
#undef TTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef balance_filt_1d_
#undef balance_filt_2d_
#undef balance_filt_3d_
