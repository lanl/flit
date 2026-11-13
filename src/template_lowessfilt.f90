!
! © 2024-2026. Triad National Security, LLC. All rights reserved.
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

#define lowess_1d_      CONCAT(lowess_1d, T)
#define lowess_2d_      CONCAT(lowess_2d, T)
#define lowess_3d_      CONCAT(lowess_3d, T)

recursive function lowess_1d_(x, y, xx, order, window, kernel, robust) result(yy)

    TT, dimension(:) :: x, y, xx
    integer, optional :: order
    TT, optional :: window
    character(len=*), optional :: kernel
    logical, optional :: robust
    TT, allocatable, dimension(:) :: yy

    TT, allocatable, dimension(:, :) :: b, btwb, btw
    TT, allocatable, dimension(:) :: w, beta
    integer :: n, nn, i, j
    integer :: lowess_order
    TT :: lowess_window
    TT, allocatable, dimension(:) :: rw, y_est, residual
    logical :: lowess_robust
    character(len=64) :: lowess_kernel
    TT :: rx

    if (present(order)) then
        lowess_order = order
    else
        lowess_order = 1
    end if

    if (present(window)) then
        lowess_window = window
    else
        lowess_window = 0.1
    end if
    call assert(lowess_window > 0, ' <lowess_1d> Error: window must > 0.')

    if (present(kernel)) then
        lowess_kernel = kernel
    else
        lowess_kernel = 'tri-cube'
    end if

    if (present(robust)) then
        lowess_robust = robust
    else
        lowess_robust = .false.
    end if

    n = size(x)
    nn = size(xx)
    rx = lowess_window*rov(x)

    b = zeros(n, lowess_order + 1)
    b(:, 1) = 1.0

    w = zeros(n)
    yy = zeros(nn)
    if (lowess_robust) then
        rw = zeros(n)
        y_est = lowess_1d_(x, y, x, window=lowess_window, kernel='epanechnikov', robust=.false.)
        residual = y_est - y
        residual = residual/(6*median(abs(return_normal(residual))))
        rw = return_normal(kernel_bisquare(residual))
    end if
    do i = 1, nn

        ! Polynomial basis centered at the output point and scaled by the window,
        ! which keeps the normal equations well conditioned
        do j = 1, lowess_order
            b(:, j + 1) = ((x - xx(i))/rx)**j
        end do

        w = abs(xx(i) - x)/rx
        select case(lowess_kernel)
            case('tri-cube')
                w = kernel_tricube(w)
            case('bi-square')
                w = kernel_bisquare(w)
            case('epanechnikov')
                w = kernel_epanechnikov(w)
        end select
        if (lowess_robust) then
            w = w*rw
        end if

        btw = xdiag(transpose(b), w)
        btwb = matx(btw, b)

        beta = matx(solve(btwb, btw), y)

        ! All basis functions except the constant one are zero at the output point
        yy(i) = beta(1)

    end do

end function lowess_1d_

recursive function lowess_2d_(x, y, z, xx, yy, order, window, kernel, robust) result(zz)

    TT, dimension(:) :: x, y, z, xx, yy
    integer, optional :: order
    TT, dimension(1:2), optional :: window
    character(len=*), optional :: kernel
    logical, optional :: robust
    TT, allocatable, dimension(:) :: zz

    TT, allocatable, dimension(:, :) :: b, btwb, btw
    TT, allocatable, dimension(:) :: w, beta
    integer :: n, nn, i, j
    integer :: lowess_order
    TT, dimension(1:2) :: lowess_window
    logical :: lowess_robust
    character(len=64) :: lowess_kernel
    TT, allocatable, dimension(:) :: rw, z_est, residual
    TT :: rx, ry

    if (present(order)) then
        lowess_order = order
    else
        lowess_order = 1
    end if

    if (present(window)) then
        lowess_window = window
    else
        lowess_window = [0.1, 0.1]
    end if
    call assert(all(lowess_window > 0), ' <lowess_2d> Error: window must > 0.')

    if (present(kernel)) then
        lowess_kernel = kernel
    else
        lowess_kernel = 'tri-cube'
    end if

    if (present(robust)) then
        lowess_robust = robust
    else
        lowess_robust = .false.
    end if

    n = size(x)
    nn = size(xx)
    rx = lowess_window(1)*rov(x)
    ry = lowess_window(2)*rov(y)

    b = zeros(n, 2*lowess_order + 1)
    b(:, 1) = 1.0

    w = zeros(n)
    zz = zeros(nn)

    if (lowess_robust) then
        rw = zeros(n)
        z_est = lowess_2d_(x, y, z, x, y, lowess_order, lowess_window, kernel='epanechnikov', robust=.false.)
        residual = z_est - z
        residual = residual/(6*median(abs(return_normal(residual))))
        rw = return_normal(kernel_bisquare(residual))
    end if

    do i = 1, nn

        ! Polynomial basis centered at the output point and scaled by the window,
        ! which keeps the normal equations well conditioned
        do j = 1, lowess_order
            b(:, 1 + 2*(j - 1) + 1) = ((x - xx(i))/rx)**j
            b(:, 1 + 2*j) = ((y - yy(i))/ry)**j
        end do

        w = sqrt(((xx(i) - x)/rx)**2 + ((yy(i) - y)/ry)**2)
        select case(lowess_kernel)
            case('tri-cube')
                w = kernel_tricube(w)
            case('bi-square')
                w = kernel_bisquare(w)
            case('epanechnikov')
                w = kernel_epanechnikov(w)
        end select
        if (lowess_robust) then
            w = w*rw
        end if

        btw = xdiag(transpose(b), w)
        btwb = matx(btw, b)

        beta = matx(solve(btwb, btw), z)

        ! All basis functions except the constant one are zero at the output point
        zz(i) = beta(1)

    end do

end function lowess_2d_

recursive function lowess_3d_(x, y, z, v, xx, yy, zz, order, window, kernel, robust) result(vv)

    TT, dimension(:) :: x, y, z, v, xx, yy, zz
    integer, optional :: order
    TT, dimension(1:3), optional :: window
    character(len=*), optional :: kernel
    logical, optional :: robust
    TT, allocatable, dimension(:) :: vv

    TT, allocatable, dimension(:, :) :: b, btwb, btw
    TT, allocatable, dimension(:) :: w, beta
    integer :: n, nn, i, j
    integer :: lowess_order
    TT, dimension(1:3) :: lowess_window
    logical :: lowess_robust
    character(len=64) :: lowess_kernel
    TT, allocatable, dimension(:) :: rw, v_est, residual
    TT :: rx, ry, rz

    if (present(order)) then
        lowess_order = order
    else
        lowess_order = 1
    end if

    if (present(window)) then
        lowess_window = window
    else
        lowess_window = [0.1, 0.1, 0.1]
    end if
    call assert(all(lowess_window > 0), ' <lowess_3d> Error: window must > 0.')

    if (present(kernel)) then
        lowess_kernel = kernel
    else
        lowess_kernel = 'tri-cube'
    end if

    if (present(robust)) then
        lowess_robust = robust
    else
        lowess_robust = .false.
    end if

    n = size(x)
    nn = size(xx)
    rx = lowess_window(1)*rov(x)
    ry = lowess_window(2)*rov(y)
    rz = lowess_window(3)*rov(z)

    b = zeros(n, 3*lowess_order + 1)
    b(:, 1) = 1.0

    w = zeros(n)
    vv = zeros(nn)

    if (lowess_robust) then
        rw = zeros(n)
        v_est = lowess_3d_(x, y, z, v, x, y, z, lowess_order, lowess_window, kernel='epanechnikov', robust=.false.)
        residual = v_est - v
        residual = residual/(6*median(abs(return_normal(residual))))
        rw = return_normal(kernel_bisquare(residual))
    end if

    do i = 1, nn

        ! Polynomial basis centered at the output point and scaled by the window,
        ! which keeps the normal equations well conditioned
        do j = 1, lowess_order
            b(:, 1 + 3*(j - 1) + 1) = ((x - xx(i))/rx)**j
            b(:, 1 + 3*(j - 1) + 2) = ((y - yy(i))/ry)**j
            b(:, 1 + 3*j) = ((z - zz(i))/rz)**j
        end do

        w = sqrt(((xx(i) - x)/rx)**2 + ((yy(i) - y)/ry)**2 + ((zz(i) - z)/rz)**2)
        select case(lowess_kernel)
            case('tri-cube')
                w = kernel_tricube(w)
            case('bi-square')
                w = kernel_bisquare(w)
            case('epanechnikov')
                w = kernel_epanechnikov(w)
        end select
        if (lowess_robust) then
            w = w*rw
        end if

        btw = xdiag(transpose(b), w)
        btwb = matx(btw, b)

        beta = matx(solve(btwb, btw), v)

        ! All basis functions except the constant one are zero at the output point
        vv(i) = beta(1)

    end do

end function lowess_3d_

#undef T
#undef TT
#undef TTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef lowess_1d_
#undef lowess_2d_
#undef lowess_3d_
