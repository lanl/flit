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

#define interp_nearest_1d_     CONCAT(interp_nearest_1d, T)
#define interp_nearest_2d_     CONCAT(interp_nearest_2d, T)
#define interp_nearest_3d_     CONCAT(interp_nearest_3d, T)
#define interp_linear_1d_     CONCAT(interp_linear_1d, T)
#define interp_cubic_spline_1d_     CONCAT(interp_cubic_spline_1d, T)
#define interp_biharmonic_1d_     CONCAT(interp_biharmonic_1d, T)
#define interp_biharmonic_2d_     CONCAT(interp_biharmonic_2d, T)
#define interp_biharmonic_3d_     CONCAT(interp_biharmonic_3d, T)
#define interp_sinc_1d_     CONCAT(interp_sinc_1d, T)
#define meshgrid_      CONCAT(meshgrid, T)
#define meshgrid_1d_      CONCAT(meshgrid_1d, T)
#define meshgrid_2d_      CONCAT(meshgrid_2d, T)
#define meshgrid_3d_      CONCAT(meshgrid_3d, T)
#define resample_1d_      CONCAT(resample_1d, T)
#define resample_2d_      CONCAT(resample_2d, T)
#define resample_3d_      CONCAT(resample_3d, T)
#define interp_to_1d_      CONCAT(interp_to_1d, T)
#define interp_to_2d_      CONCAT(interp_to_2d, T)
#define interp_to_3d_      CONCAT(interp_to_3d, T)
#define interp_like_1d_      CONCAT(interp_like_1d, T)
#define interp_like_2d_      CONCAT(interp_like_2d, T)
#define interp_like_3d_      CONCAT(interp_like_3d, T)
#define reg_to_reg_interp_1d_     CONCAT(reg_to_reg_interp_1d, T)
#define reg_to_reg_interp_2d_     CONCAT(reg_to_reg_interp_2d, T)
#define reg_to_reg_interp_3d_     CONCAT(reg_to_reg_interp_3d, T)
#define irreg_to_irreg_interp_1d_     CONCAT(irreg_to_irreg_interp_1d, T)
#define irreg_to_irreg_interp_2d_     CONCAT(irreg_to_irreg_interp_2d, T)
#define irreg_to_irreg_interp_3d_     CONCAT(irreg_to_irreg_interp_3d, T)
#define irreg_to_reg_interp_1d_     CONCAT(irreg_to_reg_interp_1d, T)
#define irreg_to_reg_interp_2d_     CONCAT(irreg_to_reg_interp_2d, T)
#define irreg_to_reg_interp_3d_     CONCAT(irreg_to_reg_interp_3d, T)
#define inpaint_1d_     CONCAT(inpaint_1d, T)
#define inpaint_2d_     CONCAT(inpaint_2d, T)
#define inpaint_3d_     CONCAT(inpaint_3d, T)

#define point_interp_linear_1d_     CONCAT(point_interp_linear_1d, T)
#define point_interp_linear_2d_     CONCAT(point_interp_linear_2d, T)
#define point_interp_linear_3d_     CONCAT(point_interp_linear_3d, T)
#define point_interp_barycentric_2d_     CONCAT(point_interp_barycentric_2d, T)
#define point_interp_barycentric_3d_     CONCAT(point_interp_barycentric_3d, T)

! external
#define interp_cspline_1d_     CONCAT(interp_cspline_1d, T)
#define interp_pchip_1d_     CONCAT(interp_pchip_1d, T)
#define interp_quintic_1d_     CONCAT(interp_quintic_1d, T)
#define interp_mba_1d_     CONCAT(interp_mba_1d, T)
#define interp_mba_2d_     CONCAT(interp_mba_2d, T)
#define interp_mba_3d_     CONCAT(interp_mba_3d, T)

!
!> Nearest neighbour interpolation for 1D data
!
subroutine interp_nearest_1d_(n, x, y, nn, xx, yy)

    integer, intent(in) :: n, nn
    TT, dimension(:), intent(in) :: x, xx
    TTT, dimension(:), intent(in) :: y
    TTT, dimension(:), intent(out) :: yy

    integer :: i

    call assert(size(x) == size(y), &
        ' <interp_nearest_1d> Error: size(x) /= size(y). ')

    i = n

    ! A brute-foce nearest neighbour implementation; could be slow for large data
    do i = 1, nn
        yy(i) = y(minloc(abs(xx(i) - x), dim=1))
    end do

end subroutine interp_nearest_1d_

subroutine interp_nearest_2d_(n, x, y, z, nn, xx, yy, zz)

    integer, intent(in) :: n, nn
    TT, dimension(:), intent(in) :: x, y, xx, yy
    TTT, dimension(:), intent(in) :: z
    TTT, dimension(:), intent(out) :: zz

    integer :: i

    call assert(size(x) == size(y) .and. size(y) == size(z), &
        ' <interp_nearest_3d> Error: Sizes of x, y, z are inconsistent. ')

    i = n

    ! A brute-foce nearest neighbour implementation; could be slow for large data
    !$omp parallel do private(i)
    do i = 1, nn
        zz(i) = z(minloc((xx(i) - x)**2 + (yy(i) - y)**2, dim=1))
    end do
    !$omp end parallel do

end subroutine interp_nearest_2d_

subroutine interp_nearest_3d_(n, x, y, z, v, nn, xx, yy, zz, vv)

    integer, intent(in) :: n, nn
    TT, dimension(:), intent(in) :: x, y, z, xx, yy, zz
    TTT, dimension(:), intent(in) :: v
    TTT, dimension(:), intent(out) :: vv

    integer :: i

    call assert(size(x) == size(y) .and. size(y) == size(z) .and. size(z) == size(v), &
        ' <interp_nearest_3d> Error: Sizes of x, y, z, v are inconsistent. ')

    i = n

    ! A brute-foce nearest neighbour implementation; could be slow for large data
    !$omp parallel do private(i)
    do i = 1, nn
        vv(i) = v(minloc((xx(i) - x)**2 + (yy(i) - y)**2 + (zz(i) - z)**2, dim=1))
    end do
    !$omp end parallel do

end subroutine interp_nearest_3d_

!
!> Linear interpolation for 1D data
!
subroutine interp_linear_1d_(n, x, y, nn, xx, yy)

    integer, intent(in) :: n, nn
    TT, dimension(:), intent(in) :: x, xx
    TTT, dimension(:), intent(in) :: y
    TTT, dimension(:), intent(out) :: yy

    integer :: i, k
    TT :: t

    call assert(size(x) == size(y), ' <interp_linear_1d> Error: size(x) /= size(y). ')

    yy = 0.0d0
    if (n == 1) then
        yy(1:nn) = y(1)
        return
    end if

    do i = 1, nn

        if (xx(i) <= x(1)) then

            t = (xx(i) - x(1))/(x(2) - x(1))
            yy(i) = (1.0d0 - t)*y(1) + t*y(2)

        else if (x(n) <= xx(i)) then

            t = (xx(i) - x(n - 1))/(x(n) - x(n - 1))
            yy(i) = (1.0 - t)*y(n - 1) + t*y(n)

        else

            do k = 2, n

                if (xx(i) >= x(k - 1) .and. xx(i) < x(k)) then

                    t = (xx(i) - x(k - 1))/(x(k) - x(k - 1))
                    yy(i) = (1.0d0 - t)*y(k - 1) + t*y(k)
                    exit

                end if

            end do

        end if

    end do

end subroutine interp_linear_1d_

!
!> Cubic spline interpolation for 1D data; contains three methods:
!>      - c2: cubic spline, a c2-continuous spline
!>      - hermite: Hermite spline, where each piece is a third-degree polynomial
!>          specified in Hermite form, that is, by its values and
!>          first derivatives at the end points of the corresponding domain interval.
!>      - monotonic: cubic monotonic Hermite spline, a.k.a. PCHIP,
!>          is a variant of cubic interpolation
!>          that preserves monotonicity of the data set being interpolated.
!>          Here I implemented the algorithm described in
!>          https://en.wikipedia.org/wiki/Monotone_cubic_interpolation,
!>          with modifications (slope selection and boundary points handling).
!
subroutine interp_cubic_spline_1d_(n, x, y, nn, xx, yy, method)

    integer, intent(in) :: n, nn
    TT, dimension(:), intent(in) :: x, xx
    TTT, dimension(:), intent(in) :: y
    character(len=*), intent(in) :: method
    TTT, dimension(:), intent(out) :: yy

    TTT, allocatable, dimension(:) :: a, b, c, d, h, r
    TTT, allocatable, dimension(:, :) :: m
    integer :: i, j
    TT :: dist
    TTT :: f1, f2

    TTT, allocatable, dimension(:) :: delta, mm, alpha, beta, tau
    TTT :: h00, h01, h10, h11
    TTT :: z(1:3)

    call assert(size(x) == size(y), ' <interp_cubic_spline_1d> Error: size(x) /= size(y). ')
    call assert(size(x) >= 3, ' <interp_cubic_spline_1d> Error: size(x) must >= 3. ')

    a = y
    b = zeros(n)
    d = zeros(n)
    h = zeros(n)
    do i = 1, n - 1
        h(i) = x(i + 1) - x(i)
    end do
    h(n) = h(n - 1)

    select case (method)

        case ('c2')
            ! Cubic spline

            m = zeros(3, n - 2)
            r = zeros(n - 2)
            ! Super-diagonal elements
            do i = 2, n - 2
                m(1, i) = h(i)/3.0
            end do
            ! Diagonal elements
            do i = 2, n - 1
                m(2, i - 1) = (h(i - 1) + h(i))*2.0/3.0
            end do
            ! Sub-diagonal elements
            do i = 2, n - 2
                m(3, i - 1) = h(i - 1)/3.0
            end do
            ! RHS
            do i = 2, n - 1
                r(i - 1) = (y(i + 1) - y(i))/h(i) - (y(i) - y(i - 1))/h(i - 1)
            end do

            c = [nTTT(0.0/2.0), solve_band(m, 1, 1, r), nTTT(0.0/2.0)]
            do i = 1, n - 1
                d(i) = (c(i + 1) - c(i))/(3.0*h(i))
                b(i) = (y(i + 1) - y(i))/h(i) - (2*c(i) + c(i + 1))*h(i)/3.0
            end do
            d(n) = 0.0
            b(n) = 3*d(n - 1)*h(n - 1)**2 + 2*c(n - 1)*h(n - 1) + b(n - 1)

            do i = 1, nn
                if (xx(i) < x(1)) then
                    dist = xx(i) - x(1)
                    yy(i) = a(1) + b(1)*dist + c(1)*dist**2

                else if (xx(i) >= x(n)) then
                    dist = xx(i) - x(n)
                    yy(i) = a(n) + b(n)*dist + c(n)*dist**2

                else
                    do j = 1, n - 1
                        if (xx(i) >= x(j) .and. xx(i) < x(j + 1)) then
                            dist = xx(i) - x(j)
                            yy(i) = a(j) + b(j)*dist + c(j)*dist**2 + d(j)*dist**3
                            exit
                        end if
                    end do
                end if
            end do

        case ('hermite')
            ! Cubic Hermite spline

            c = zeros(n)
            do i = 2, n - 1
                f1 = (y(i) - y(i - 1))/h(i - 1)
                f2 = (y(i + 1) - y(i))/h(i)
                b(i) = h(i - 1)/(h(i - 1) + h(i))*f2 + h(i)/(h(i) + h(i - 1))*f1
            end do
            b(1) = 0.5*(-b(2) + 3*(y(2) - y(1))/h(1))
            b(n) = 0.5*(-b(n - 1) + 3*(y(n) - y(n - 1))/h(n - 1))

            do i = 1, n - 1
                c(i) = -(2*b(i) + b(i + 1))/h(i) + 3*(a(i + 1) - a(i))/h(i)**2
                d(i) = -2*c(i)/(3*h(i)) + (b(i + 1) - b(i))/(3*h(i)**2)
            end do

            do i = 1, nn
                if (xx(i) < x(1)) then
                    dist = xx(i) - x(1)
                    yy(i) = a(1) + b(1)*dist + c(1)*dist**2

                else if (xx(i) >= x(n)) then
                    dist = xx(i) - x(n)
                    yy(i) = a(n) + b(n)*dist + c(n)*dist**2

                else
                    do j = 1, n - 1
                        if (xx(i) >= x(j) .and. xx(i) < x(j + 1)) then
                            dist = xx(i) - x(j)
                            yy(i) = a(j) + b(j)*dist + c(j)*dist**2 + d(j)*dist**3
                            exit
                        end if
                    end do
                end if
            end do

        case ('monotonic')
            ! PCHIP spline

            ! For handling boundary points
            c = zeros(n)
            do i = 2, n - 1
                f1 = (y(i) - y(i - 1))/h(i - 1)
                f2 = (y(i + 1) - y(i))/h(i)
                b(i) = h(i - 1)/(h(i - 1) + h(i))*f2 + h(i)/(h(i) + h(i - 1))*f1
            end do
            b(1) = 0.5*(-b(2) + 3*(y(2) - y(1))/h(1))
            b(n) = 0.5*(-b(n - 1) + 3*(y(n) - y(n - 1))/h(n - 1))

            do i = n - 1, n - 1
                c(i) = -(2*b(i) + b(i + 1))/h(i) + 3*(a(i + 1) - a(i))/h(i)**2
            end do
            c(n) = c(n - 1)

            ! Monotonicity-perserving cubic splines
            delta = zeros(n)
            do i = 1, n - 1
                delta(i) = (y(i + 1) - y(i))/(x(i + 1) - x(i))
            end do
            delta(n) = delta(n - 1)
            mm = zeros(n)
            do i = 2, n - 1
                ! Choosing the min slope in L1 sense
                z = [0.5*(delta(i - 1) + delta(i)), delta(i - 1), delta(i)]
                mm(i) = z(minloc(abs(z), dim=1))
                if (delta(i - 1)*delta(i) < 0) then
                    mm(i) = 0
                end if
            end do
            mm(1) = delta(1)
            mm(n) = delta(n - 1)
            do i = 1, n - 1
                if (delta(i) == 0) then
                    mm(i:i + 1) = 0
                end if
            end do

            alpha = zeros(n)
            beta = zeros(n)
            do i = 1, n - 1
                if (mm(i) /= 0) then
                    alpha(i) = mm(i)/delta(i)
                    beta(i) = mm(i + 1)/delta(i)
                end if
            end do
            do i = 1, n - 1
                if (alpha(i) < 0) then
                    mm(i) = 0
                end if
                if (beta(i) < 0) then
                    mm(i + 1) = 0
                end if
            end do

            tau = zeros(n)
            do i = 1, n - 1
                if (alpha(i)**2 + beta(i)**2 > 9) then
                    tau(i) = 9/sqrt(alpha(i)**2 + beta(i)**2)
                end if
            end do
            do i = 1, n - 1
                if (alpha(i)**2 + beta(i)**2 > 9 .and. mm(i) /= 0) then
                    mm(i) = tau(i)*alpha(i)*delta(i)
                    mm(i + 1) = tau(i)*beta(i)*delta(i)
                end if
            end do

            yy = 0
            do i = 1, nn
                if (xx(i) < x(1)) then
                    dist = xx(i) - x(1)
                    yy(i) = a(1) + b(1)*dist + c(1)*dist**2
                else if (xx(i) >= x(n)) then
                    dist = xx(i) - x(n)
                    yy(i) = a(n) + b(n)*dist + c(n)*dist**2
                else
                    do j = 1, n - 1
                        if (xx(i) >= x(j) .and. xx(i) < x(j + 1)) then
                            dist = (xx(i) - x(j))/h(j)
                            h00 = (1 + 2*dist)*(1 - dist)**2
                            h10 = dist*(1 - dist)**2
                            h01 = dist**2*(3 - 2*dist)
                            h11 = dist**2*(dist - 1)
                            yy(i) = h00*y(j) + h10*h(j)*mm(j) + h01*y(j + 1) + h11*h(j)*mm(j + 1)
                            exit
                        end if
                    end do
                end if
            end do

    end select

    !    yy = 0
    !    select case (order)
    !        case (0)
    !            ! The function itself
    !            do i = 1, nn
    !                if (xx(i) < x(1)) then
    !                    dist = xx(i) - x(1)
    !                    yy(i) = a(1) + b(1)*dist + c(1)*dist**2
    !
    !                else if (xx(i) >= x(n)) then
    !                    dist = xx(i) - x(n)
    !                    yy(i) = a(n) + b(n)*dist + c(n)*dist**2
    !                else
    !                    do j = 1, n - 1
    !                        if (xx(i) >= x(j) .and. xx(i) < x(j + 1)) then
    !                            dist = xx(i) - x(j)
    !                            yy(i) = a(j) + b(j)*dist + c(j)*dist**2 + d(j)*dist**3
    !                            exit
    !                        end if
    !                    end do
    !                end if
    !            end do
    !        case (1)
    !            ! First-order derivative
    !            do i = 1, nn
    !                if (xx(i) < x(1)) then
    !                    dist = xx(i) - x(1)
    !                    yy(i) = b(1) + 2*c(1)*dist
    !                else if (xx(i) >= x(n)) then
    !                    dist = xx(i) - x(n)
    !                    yy(i) = b(n) + 2*c(n)*dist
    !                else
    !                    do j = 1, n - 1
    !                        if (xx(i) >= x(j) .and. xx(i) < x(j + 1)) then
    !                            dist = xx(i) - x(j)
    !                            yy(i) = b(j) + 2*c(j)*dist + 3*d(j)*dist**2
    !                            exit
    !                        end if
    !                    end do
    !                end if
    !            end do
    !        case (2)
    !            ! Second-order derivative
    !            do i = 1, nn
    !                if (xx(i) < x(1)) then
    !                    yy(i) = 2*c(1)
    !                else if (xx(i) >= x(n)) then
    !                    yy(i) = 2*c(n)
    !                else
    !                    do j = 1, n - 1
    !                        if (xx(i) >= x(j) .and. xx(i) < x(j + 1)) then
    !                            dist = xx(i) - x(j)
    !                            yy(i) = 2*c(j) + 6*d(j)*dist
    !                            exit
    !                        end if
    !                    end do
    !                end if
    !            end do
    !    end select

end subroutine interp_cubic_spline_1d_

!
!> Biharmonic interpolation; could be very slow and also memory-intensive
!> because needs to solve a large linear system;
!> https://doi.org/10.1029/GL014i002p00139
!
subroutine interp_biharmonic_1d_(n_, x, y, nn, xx, yy)

    integer, intent(in) :: n_, nn
    TT, dimension(:), intent(in) :: x, xx
    TTT, dimension(:), intent(in) :: y
    TTT, dimension(:), intent(out) :: yy

    integer :: i, j, n
    TTT, allocatable, dimension(:, :) :: dist
    TTT, allocatable, dimension(:) :: w
    TT, allocatable, dimension(:) :: xd
    TTT, allocatable, dimension(:) :: yd
    TTT, allocatable, dimension(:, :) :: xyd

    call assert(size(x) == size(y), ' <interp_biharmonic_1d> Error: size(x) /= size(y). ')

    ! First find unique pairs of original scatter points
    ! Otherwise the linear system is not solvable
    xyd = unique(reshape([nTTT(x), y], [n_, 2]), cols=[1])
    n = size(xyd, 1)
    xd = xyd(:, 1)
    yd = xyd(:, 2)

    ! Compute distance and Green's function
    dist = zeros(n, n)
    do j = 1, n
        do i = 1, n
            if (i /= j) then
                dist(i, j) = abs(xd(i) - xd(j))
                dist(i, j) = dist(i, j)**2*(log(dist(i, j)) - 1.0)
            else
                dist(i, j) = 0.0
            end if
        end do
    end do

    ! Compute weights
    w = solve(dist, yd)

    ! Comptue interpolation values
    dist = zeros(nn, n)
    do j = 1, n
        do i = 1, nn
            dist(i, j) = abs(xx(i) - xd(j))
            if (abs(dist(i, j)) /= 0) then
                dist(i, j) = dist(i, j)**2*(log(dist(i, j)) - 1.0)
            end if
        end do
    end do

    yy = matx(dist, w)

end subroutine interp_biharmonic_1d_

subroutine interp_biharmonic_2d_(n_, x, y, z, nn, xx, yy, zz)

    integer, intent(in) :: n_, nn
    TT, dimension(:), intent(in) :: x, y, xx, yy
    TTT, dimension(:), intent(in) :: z
    TTT, dimension(:), intent(out) :: zz

    integer :: i, j, n
    TTT, allocatable, dimension(:, :) :: dist
    TTT, allocatable, dimension(:) :: w
    TT, allocatable, dimension(:) :: xd, yd
    TTT, allocatable, dimension(:) :: zd
    TTT, allocatable, dimension(:, :) :: xyzd

    call assert(size(x) == size(y) .and. size(y) == size(z), &
        ' <interp_biharmonic_2d> Error: Sizes of x, y, and z are inconsistent. ')

    ! First find unique pairs of original scatter points
    ! Otherwise the linear system is not solvable
    xyzd = unique(reshape([nTTT(x), nTTT(y), z], [n_, 3]), cols=[1, 2])
    n = size(xyzd, 1)
    xd = xyzd(:, 1)
    yd = xyzd(:, 2)
    zd = xyzd(:, 3)

    dist = zeros(n, n)

    ! Compute distance and Green's function
    !$omp parallel do private(i, j)
    do j = 1, n
        do i = 1, n
            dist(i, j) = sqrt((xd(i) - xd(j))**2 + (yd(i) - yd(j))**2)
            if (i /= j) then
                dist(i, j) = dist(i, j)**2*(log(dist(i, j)) - 1.0)
            else
                dist(i, j) = 0.0
            end if
        end do
    end do
    !$omp end parallel do

    ! Compute weights
    w = solve(dist, zd)

    ! Comptue interpolation values
    dist = zeros(nn, n)
    !$omp parallel do private(i, j)
    do j = 1, n
        do i = 1, nn
            dist(i, j) = sqrt((xx(i) - xd(j))**2 + (yy(i) - yd(j))**2)
            if (abs(dist(i, j)) /= 0) then
                dist(i, j) = dist(i, j)**2*(log(dist(i, j)) - 1.0)
            end if
        end do
    end do
    !$omp end parallel do

    zz = matx(dist, w)

end subroutine interp_biharmonic_2d_

subroutine interp_biharmonic_3d_(n_, x, y, z, v, nn, xx, yy, zz, vv)

    integer, intent(in) :: n_, nn
    TT, dimension(:), intent(in) :: x, y, z, xx, yy, zz
    TTT, dimension(:), intent(in) :: v
    TTT, dimension(:), intent(out) :: vv

    integer :: i, j, n
    TTT, allocatable, dimension(:, :) :: dist
    TTT, allocatable, dimension(:) :: w
    TT, allocatable, dimension(:) :: xd, yd, zd
    TTT, allocatable, dimension(:) :: vd
    TTT, allocatable, dimension(:, :) :: xyzvd

    call assert(size(x) == size(y) .and. size(y) == size(z) .and. size(z) == size(v), &
        ' <interp_biharmonic_3d> Error: Sizes of x, y, z, and v are inconsistent. ')

    ! First find unique pairs of original scatter points
    ! Otherwise the linear system is not solvable
    xyzvd = unique(reshape([nTTT(x), nTTT(y), nTTT(z), v], [n_, 4]), cols=[1, 2, 3])
    n = size(xyzvd, 1)
    xd = xyzvd(:, 1)
    yd = xyzvd(:, 2)
    zd = xyzvd(:, 3)
    vd = xyzvd(:, 4)

    ! Compute distance and Green's function
    dist = zeros(n, n)
    !$omp parallel do private(i, j)
    do j = 1, n
        do i = 1, n
            dist(i, j) = sqrt((xd(i) - xd(j))**2 + (yd(i) - yd(j))**2 &
                + (zd(i) - zd(j))**2)
            if (i /= j) then
                dist(i, j) = dist(i, j)**2*(log(dist(i, j)) - 1.0)
            else
                dist(i, j) = 0.0
            end if
        end do
    end do
    !$omp end parallel do

    ! Compute weights
    w = solve(dist, vd)

    ! Comptue interpolation values
    dist = zeros(nn, n)
    !$omp parallel do private(i, j)
    do j = 1, n
        do i = 1, nn
            dist(i, j) = sqrt((xx(i) - xd(j))**2 + (yy(i) - yd(j))**2 &
                + (zz(i) - zd(j))**2)
            if (abs(dist(i, j)) /= 0) then
                dist(i, j) = dist(i, j)**2*(log(dist(i, j)) - 1.0)
            end if
        end do
    end do
    !$omp end parallel do

    vv = matx(dist, w)

end subroutine interp_biharmonic_3d_

!
!> Windowed sinc interpolation for 1D data
!
subroutine interp_sinc_1d_(n, x, y, nn, xx, yy)

    integer, intent(in) :: n, nn
    TT, dimension(:), intent(in) :: x, xx
    TTT, dimension(:), intent(in) :: y
    TTT, dimension(:), intent(out) :: yy

    integer :: i, h, nk1, nk2
    TT :: d, dist
    TTT :: ya, yb
    TT, allocatable, dimension(:) :: ksinc, xd
    TTT, allocatable, dimension(:) :: yd
    integer :: nkw
    TTT :: b0

    nkw = 8
    b0 = 3.0d0

    xd = x
    yd = y
    d = x(2) - x(1)
    if (xx(1) < x(1)) then
        nk1 = nkw + ceiling(abs(xx(1) - x(1))/d)
    else
        nk1 = nkw
    end if
    if (xx(nn) > x(n)) then
        nk2 = nkw + ceiling(abs(xx(nn) - x(n))/d)
    else
        nk2 = nkw
    end if
    call pad_array(xd, [nk1, nk2])
    call pad_array(yd, [nk1, nk2])
    ya = y(2) - y(1)
    yb = y(n) - y(n - 1)
    do i = 1, nk1
        xd(1 - i) = xd(1) - i*d
        yd(1 - i) = yd(1) - i*ya
    end do
    do i = 1, nk2
        xd(n + i) = xd(n) + i*d
        yd(n + i) = yd(n) + i*yb
    end do

    do i = 1, nn
        h = nint((xx(i) - xd(1))/d) + 1
        dist = (xx(i) - xd(h))/d
        ksinc = dist - regspace(-nkw, 1, nkw)
        ksinc = sinc(ksinc*const_pi)*window_function(ksinc/nkw + 0.5d0, method='kaiser', alpha=b0 + 0.0d0)
        yy(i) = sum(ksinc*yd(h - nkw:h + nkw))
    end do

end subroutine interp_sinc_1d_

!==================================================================================================
!
!> Interpolate regularly sampled data to regularly sampled data
!
function reg_to_reg_interp_1d_(f, n1, d1, o1, nn1, dd1, oo1, method) result(ff)

    TTT, dimension(:) :: f
    integer :: n1, nn1
    TT :: d1, dd1, o1, oo1
    character(len=*), optional :: method
    TTT, allocatable, dimension(:) :: ff

    TT, allocatable, dimension(:) :: x, xx
    character(len=24) :: interp_method
    integer :: i

    call assert(d1 > 0, ' <reg_to_reg_interp_1d> Error: The original interval must be > 0')
    call assert(dd1 > 0, ' <reg_to_reg_interp_1d> Error: The resampling interval must be > 0')

    if (present(method)) then
        interp_method = method
    else
        interp_method = 'linear'
    end if

    allocate (ff(1:nn1))

    if (n1 == nn1 .and. d1 == dd1 .and. o1 == oo1) then
        ! Trivial case
        ff = f

    else
        ! Do interpolation

        ! Point coordinates before and after interpolation
        allocate (x(1:n1))
        allocate (xx(1:nn1))
        do i = 1, n1
            x(i) = o1 + (i - 1)*d1
        end do
        do i = 1, nn1
            xx(i) = oo1 + (i - 1)*dd1
        end do

        select case (interp_method)
            case ('nearest')
                call interp_nearest_1d_(n1, x, f, nn1, xx, ff)
            case ('linear')
                call interp_linear_1d_(n1, x, f, nn1, xx, ff)
            case ('sinc')
                call interp_sinc_1d_(n1, x, f, nn1, xx, ff)
            case ('cubic')
                call interp_cspline_1d_(n1, x, f, nn1, xx, ff)
            case ('pchip')
                call interp_pchip_1d_(n1, x, f, nn1, xx, ff)
            case ('quintic')
                call interp_quintic_1d_(n1, x, f, nn1, xx, ff)
            case ('mba')
                call interp_mba_1d_(n1, x, f, nn1, xx, ff)
            case ('biharmonic')
                call interp_biharmonic_1d_(n1, x, f, nn1, xx, ff)
            case ('cubic_spline')
                call interp_cubic_spline_1d_(n1, x, f, nn1, xx, ff, method='c2')
            case ('hermite_spline')
                call interp_cubic_spline_1d_(n1, x, f, nn1, xx, ff, method='hermite')
            case ('monotonic_spline')
                call interp_cubic_spline_1d_(n1, x, f, nn1, xx, ff, method='monotonic')
            case default
                call interp_linear_1d_(n1, x, f, nn1, xx, ff)
        end select

    end if

end function reg_to_reg_interp_1d_

function reg_to_reg_interp_2d_(f, n, d, o, nn, dd, oo, method) result(ff)

    TTT, dimension(:, :) :: f
    integer, dimension(:) :: n, nn
    TT, dimension(1:2) :: d, dd, o, oo
    character(len=*), dimension(1:2), optional :: method
    TTT, allocatable, dimension(:, :) :: ff

    character(len=24), dimension(1:2) :: interp_method
    integer :: i, j
    TTT, allocatable, dimension(:, :) :: w

    call assert(all(d > 0), ' <reg_to_reg_interp_2d_> Error: All original intervals must be > 0')
    call assert(all(dd > 0), ' <reg_to_reg_interp_2d_> Error: All resampling intervals must be > 0')

    if (present(method)) then
        interp_method = method
    else
        interp_method = ['linear', 'linear']
    end if

    allocate (w(1:nn(1), 1:n(2)))
    allocate (ff(1:nn(1), 1:nn(2)))

    ! Interpolate along the 1st dimension
    !$omp parallel do private(j)
    do j = 1, n(2)
        w(:, j) = reg_to_reg_interp_1d_(f(:, j), &
            n(1), d(1), o(1), nn(1), dd(1), oo(1), interp_method(1))
    end do
    !$omp end parallel do

    ! Interpolate along the 2nd dimension
    !$omp parallel do private(i)
    do i = 1, nn(1)
        ff(i, :) = reg_to_reg_interp_1d_(w(i, :), &
            n(2), d(2), o(2), nn(2), dd(2), oo(2), interp_method(2))
    end do
    !$omp end parallel do

end function reg_to_reg_interp_2d_

function reg_to_reg_interp_3d_(f, n, d, o, nn, dd, oo, method) result(ff)

    TTT, dimension(:, :, :) :: f
    integer, dimension(:) :: n, nn
    TT, dimension(1:3) :: d, dd, o, oo
    character(len=*), dimension(:), optional :: method
    TTT, allocatable, dimension(:, :, :) :: ff

    character(len=24), dimension(1:3) :: interp_method
    integer :: i, j, k
    TTT, allocatable, dimension(:, :, :) :: w1, w2

    call assert(all(d > 0), ' <reg_to_reg_interp_3d_> Error: All original intervals must be > 0')
    call assert(all(dd > 0), ' <reg_to_reg_interp_3d_> Error: All resampling intervals must be > 0')

    if (present(method)) then
        interp_method = method
    else
        interp_method = ['linear', 'linear', 'linear']
    end if

    allocate (w1(1:nn(1), 1:n(2), 1:n(3)))
    allocate (w2(1:nn(1), 1:nn(2), 1:n(3)))
    allocate (ff(1:nn(1), 1:nn(2), 1:nn(3)))

    ! Interpolate along the 1st dimension
    !$omp parallel do private(j, k) collapse(2)
    do k = 1, n(3)
        do j = 1, n(2)
            w1(:, j, k) = reg_to_reg_interp_1d_(f(:, j, k), &
                n(1), d(1), o(1), nn(1), dd(1), oo(1), interp_method(1))
        end do
    end do
    !$omp end parallel do

    ! Interpolate along the 2nd dimension
    !$omp parallel do private(i, k) collapse(2)
    do k = 1, n(3)
        do i = 1, nn(1)
            w2(i, :, k) = reg_to_reg_interp_1d_(w1(i, :, k), &
                n(2), d(2), o(2), nn(2), dd(2), oo(2), interp_method(2))
        end do
    end do
    !$omp end parallel do

    ! Interpolate along the 3rd dimension
    !$omp parallel do private(i, j) collapse(2)
    do j = 1, nn(2)
        do i = 1, nn(1)
            ff(i, j, :) = reg_to_reg_interp_1d_(w2(i, j, :), &
                n(3), d(3), o(3), nn(3), dd(3), oo(3), interp_method(3))
        end do
    end do
    !$omp end parallel do

end function reg_to_reg_interp_3d_

!
!> Interpolate a regularly sampled data using ratios; the new dimensions will be old dimensions / ratios
!
function resample_1d_(f, r, method) result(ff)

    TT, dimension(:), intent(in) :: f
    real, intent(in) :: r
    character(len=*), intent(in), optional :: method
    TT, allocatable, dimension(:) :: ff

    character(len=24) :: interp_method

    if (present(method)) then
        interp_method = method
    else
        interp_method = 'linear'
    end if

    call assert(r > 0, ' <resample_1d> Error: Resampling ratio must > 0.')
    ff = reg_to_reg_interp_1d_(f, size(f), nTT(1.0), nTT(0.0), &
        nint(r*(size(f) - 1) + 1), nTT(1.0/r), nTT(0.0), interp_method)

end function resample_1d_

function resample_2d_(f, r, method) result(ff)

    TT, dimension(:, :), intent(in) :: f
    real, dimension(:), intent(in) :: r
    character(len=*), dimension(:), intent(in), optional :: method
    TT, allocatable, dimension(:, :) :: ff

    character(len=24), dimension(1:2) :: interp_method

    if (present(method)) then
        interp_method = method
    else
        interp_method = ['linear', 'linear']
    end if

    call assert(all(r > 0), ' <resample_2d> Error: Resampling ratio must > 0.')
    ff = reg_to_reg_interp_2d_(f, shape(f), nTT([1.0, 1.0]), nTT([0.0, 0.0]), &
        nint(r*(shape(f) - 1) + 1), nTT(1.0/r), nTT([0.0, 0.0]), interp_method)

end function resample_2d_

function resample_3d_(f, r, method) result(ff)

    TT, dimension(:, :, :), intent(in) :: f
    real, dimension(:), intent(in) :: r
    character(len=*), dimension(:), intent(in), optional :: method
    TT, allocatable, dimension(:, :, :) :: ff

    character(len=24), dimension(1:3) :: interp_method

    if (present(method)) then
        interp_method = method
    else
        interp_method = ['linear', 'linear', 'linear']
    end if

    call assert(all(r > 0), ' <resample_3d> Error: Resampling ratio must > 0.')
    ff = reg_to_reg_interp_3d_(f, shape(f), nTT([1.0, 1.0, 1.0]), nTT([0.0, 0.0, 0.0]), &
        nint(r*(shape(f) - 1) + 1), nTT(1.0/r), nTT([0.0, 0.0, 0.0]), interp_method)

end function resample_3d_

!
!> Interpolate a regularly sampled data like another array
!
function interp_like_1d_(f, ff, method) result(g)

    TT, dimension(:), intent(in) :: f, ff
    character(len=*), intent(in), optional :: method
    TT, allocatable, dimension(:) :: g

    character(len=24) :: interp_method
    integer :: n, nn

    if (present(method)) then
        interp_method = method
    else
        interp_method = 'linear'
    end if

    n = size(f)
    nn = size(ff)

    g = reg_to_reg_interp_1d_(f, n, nTT(1.0), nTT(0.0), &
        nn, nTT((n - 1.0)/(nn - 1.0)), nTT(0.0), interp_method)

end function interp_like_1d_

function interp_like_2d_(f, ff, method) result(g)

    TT, dimension(:, :), intent(in) :: f, ff
    character(len=*), dimension(:), intent(in), optional :: method
    TT, allocatable, dimension(:, :) :: g

    character(len=24), dimension(1:2) :: interp_method
    integer, dimension(1:2) :: n, nn

    if (present(method)) then
        interp_method = method
    else
        interp_method = ['linear', 'linear']
    end if

    n = shape(f)
    nn = shape(ff)

    g = reg_to_reg_interp_2d_(f, n, nTT([1.0, 1.0]), nTT([0.0, 0.0]), &
        nn, nTT((n - 1.0)/(nn - 1.0)), nTT([0.0, 0.0]), interp_method)

end function interp_like_2d_

function interp_like_3d_(f, ff, method) result(g)

    TT, dimension(:, :, :), intent(in) :: f, ff
    character(len=*), dimension(:), intent(in), optional :: method
    TT, allocatable, dimension(:, :, :) :: g

    character(len=24), dimension(1:3) :: interp_method
    integer, dimension(1:3) :: n, nn

    if (present(method)) then
        interp_method = method
    else
        interp_method = ['linear', 'linear', 'linear']
    end if

    n = shape(f)
    nn = shape(ff)

    g = reg_to_reg_interp_3d_(f, n, nTT([1.0, 1.0, 1.0]), nTT([0.0, 0.0, 0.0]), &
        nn, nTT((n - 1.0)/(nn - 1.0)), nTT([0.0, 0.0, 0.0]), interp_method)

end function interp_like_3d_

!
!> Interpolate from regularly sampled data based on target dimensions
!
function interp_to_1d_(f, nn, method) result(g)

    TT, dimension(:), intent(in) :: f
    integer, intent(in) :: nn
    character(len=*), intent(in), optional :: method
    TT, allocatable, dimension(:) :: g

    character(len=24) :: interp_method
    integer :: n

    if (present(method)) then
        interp_method = method
    else
        interp_method = 'linear'
    end if

    n = size(f)

    g = reg_to_reg_interp_1d_(f, n, nTT(1.0), nTT(0.0), &
        nn, nTT((n - 1.0)/(nn - 1.0)), nTT(0.0), interp_method)

end function interp_to_1d_

function interp_to_2d_(f, nn, method) result(g)

    TT, dimension(:, :), intent(in) :: f
    integer, dimension(1:2), intent(in) :: nn
    character(len=*), dimension(:), intent(in), optional :: method
    TT, allocatable, dimension(:, :) :: g

    character(len=24), dimension(1:2) :: interp_method
    integer, dimension(1:2) :: n

    if (present(method)) then
        interp_method = method
    else
        interp_method = ['linear', 'linear']
    end if

    n = shape(f)

    g = reg_to_reg_interp_2d_(f, n, nTT([1.0, 1.0]), nTT([0.0, 0.0]), &
        nn, nTT((n - 1.0)/(nn - 1.0)), nTT([0.0, 0.0]), interp_method)

end function interp_to_2d_

function interp_to_3d_(f, nn, method) result(g)

    TT, dimension(:, :, :), intent(in) :: f
    integer, dimension(1:3), intent(in) :: nn
    character(len=*), dimension(:), intent(in), optional :: method
    TT, allocatable, dimension(:, :, :) :: g

    character(len=24), dimension(1:3) :: interp_method
    integer, dimension(1:3) :: n

    if (present(method)) then
        interp_method = method
    else
        interp_method = ['linear', 'linear', 'linear']
    end if

    n = shape(f)

    g = reg_to_reg_interp_3d_(f, n, nTT([1.0, 1.0, 1.0]), nTT([0.0, 0.0, 0.0]), &
        nn, nTT((n - 1.0)/(nn - 1.0)), nTT([0.0, 0.0, 0.0]), interp_method)

end function interp_to_3d_

!
!> Interpolate from irregularly sampled data to irregularly sampled data
!
function irreg_to_irreg_interp_1d_(x, f, xx, method) result(ff)

    TT, dimension(:), intent(in) :: x, xx
    TTT, dimension(:), intent(in) :: f
    TTT, allocatable, dimension(:) :: ff
    character(len=*), optional :: method

    character(len=24) :: interp_method
    integer :: n, nn

    if (present(method)) then
        interp_method = method
    else
        interp_method = 'mba'
    end if

    n = size(x)
    nn = size(xx)
    allocate (ff(1:nn))

    select case (interp_method)
        case ('nearest')
            call interp_nearest_1d_(n, x, f, nn, xx, ff)
        case ('linear')
            call interp_linear_1d_(n, x, f, nn, xx, ff)
        case ('cubic')
            call interp_cspline_1d_(n, x, f, nn, xx, ff)
        case ('pchip')
            call interp_pchip_1d_(n, x, f, nn, xx, ff)
        case ('quintic')
            call interp_quintic_1d_(n, x, f, nn, xx, ff)
        case ('mba')
            call interp_mba_1d_(n, x, f, nn, xx, ff)
        case ('biharmonic')
            call interp_biharmonic_1d_(n, x, f, nn, xx, ff)
        case ('cubic_spline')
            call interp_cubic_spline_1d_(n, x, f, nn, xx, ff, method='c2')
        case ('hermite_spline')
            call interp_cubic_spline_1d_(n, x, f, nn, xx, ff, method='hermite')
        case ('monotonic_spline')
            call interp_cubic_spline_1d_(n, x, f, nn, xx, ff, method='monotonic')
        case default
            call interp_linear_1d_(n, x, f, nn, xx, ff)
    end select

end function irreg_to_irreg_interp_1d_

function irreg_to_irreg_interp_2d_(x, y, f, xx, yy, method) result(ff)

    TT, dimension(:), intent(in) :: x, y, xx, yy
    TTT, dimension(:), intent(in) :: f
    TTT, allocatable, dimension(:) :: ff
    character(len=*), optional :: method

    character(len=24) :: interp_method
    integer :: n, nn

    if (present(method)) then
        interp_method = method
    else
        interp_method = 'mba'
    end if

    n = size(x)
    nn = size(xx)
    allocate (ff(1:nn))

    select case (interp_method)
        case ('nearest')
            call interp_nearest_2d_(n, x, y, f, nn, xx, yy, ff)
        case ('mba')
            call interp_mba_2d_(n, x, y, f, nn, xx, yy, ff)
        case ('biharmonic')
            call interp_biharmonic_2d_(n, x, y, f, nn, xx, yy, ff)
    end select

end function irreg_to_irreg_interp_2d_

function irreg_to_irreg_interp_3d_(x, y, z, f, xx, yy, zz, method) result(ff)

    TT, dimension(:) :: x, y, z, xx, yy, zz
    TTT, dimension(:) :: f
    TTT, allocatable, dimension(:) :: ff
    character(len=*), optional :: method

    character(len=24) :: interp_method
    integer :: n, nn

    if (present(method)) then
        interp_method = method
    else
        interp_method = 'mba'
    end if

    n = size(x)
    nn = size(xx)
    allocate (ff(1:nn))

    select case (interp_method)
        case ('nearest')
            call interp_nearest_3d_(n, x, y, z, f, nn, xx, yy, zz, ff)
        case ('mba')
            call interp_mba_3d_(n, x, y, z, f, nn, xx, yy, zz, ff)
        case ('biharmonic')
            call interp_biharmonic_3d_(n, x, y, z, f, nn, xx, yy, zz, ff)
    end select

end function irreg_to_irreg_interp_3d_

!
!> Generate mesh grids
!
function meshgrid_(n, d, o, dim) result(g)

    integer, intent(in) :: n(:)
    TT, intent(in) :: o(:), d(:)
    integer, intent(in) :: dim

    TT, allocatable, dimension(:) :: g
    integer :: i, j, k, l, n1, n2, n3
    TT :: o1, o2, o3, d1, d2, d3

    select case (size(n))

        case (1)
            g = zeros(n(1))
            do i = 1, n(1)
                g(i) = o(1) + (i - 1)*d(1)
            end do

        case (2)
            n1 = n(1)
            n2 = n(2)
            o1 = o(1)
            o2 = o(2)
            d1 = d(1)
            d2 = d(2)
            allocate (g(1:n1*n2))
            select case (dim)
                case (1)
                    l = 1
                    do j = 1, n2
                        do i = 1, n1
                            g(l) = o1 + (i - 1)*d1
                            l = l + 1
                        end do
                    end do
                case (2)
                    l = 1
                    do j = 1, n2
                        do i = 1, n1
                            g(l) = o2 + (j - 1)*d2
                            l = l + 1
                        end do
                    end do
            end select

        case (3)
            n1 = n(1)
            n2 = n(2)
            n3 = n(3)
            o1 = o(1)
            o2 = o(2)
            o3 = o(3)
            d1 = d(1)
            d2 = d(2)
            d3 = d(3)
            allocate (g(1:n1*n2*n3))
            select case (dim)
                case (1)
                    l = 1
                    do k = 1, n3
                        do j = 1, n2
                            do i = 1, n1
                                g(l) = o1 + (i - 1)*d1
                                l = l + 1
                            end do
                        end do
                    end do
                case (2)
                    l = 1
                    do k = 1, n3
                        do j = 1, n2
                            do i = 1, n1
                                g(l) = o2 + (j - 1)*d2
                                l = l + 1
                            end do
                        end do
                    end do
                case (3)
                    l = 1
                    do k = 1, n3
                        do j = 1, n2
                            do i = 1, n1
                                g(l) = o3 + (k - 1)*d3
                                l = l + 1
                            end do
                        end do
                    end do
            end select

    end select

end function meshgrid_

!
!> Generate mesh grids from vectors; the vectors may be irregularly sampled
!
function meshgrid_1d_(v) result(g)

    TT, dimension(:) :: v
    TT, allocatable, dimension(:) :: g

    g = v

end function meshgrid_1d_

function meshgrid_2d_(v1, v2) result(g)

    TT, dimension(:) :: v1, v2
    TT, allocatable, dimension(:, :) :: g

    integer :: n1, n2, l, i, j

    n1 = size(v1)
    n2 = size(v2)

    g = zeros(n1*n2, 2)

    l = 1
    do j = 1, n2
        do i = 1, n1
            g(l, :) = [v1(i), v2(j)]
            l = l + 1
        end do
    end do

end function meshgrid_2d_

function meshgrid_3d_(v1, v2, v3) result(g)

    TT, dimension(:) :: v1, v2, v3
    TT, allocatable, dimension(:, :) :: g

    integer :: n1, n2, n3, l, i, j, k

    n1 = size(v1)
    n2 = size(v2)
    n3 = size(v3)

    g = zeros(n1*n2*n3, 3)

    l = 1
    do k = 1, n3
        do j = 1, n2
            do i = 1, n1
                g(l, :) = [v1(i), v2(j), v3(k)]
                l = l + 1
            end do
        end do
    end do

end function meshgrid_3d_

!
!> Interpolate from irregularly sampled data to regularly sampled data
!
function irreg_to_reg_interp_1d_(x, f, n, d, o, method) result(ff)

    TT, dimension(:) :: x
    TTT, dimension(:) :: f
    integer :: n
    TT :: d, o
    TTT, allocatable, dimension(:) :: ff
    character(len=*), optional :: method

    character(len=24) :: interp_method
    integer :: l
    TT, allocatable, dimension(:) :: xx

    if (present(method)) then
        interp_method = method
    else
        interp_method = 'mba'
    end if

    l = size(x)
    call assert(n > 0 .and. d > 0, ' <irreg_to_reg_interp_1d> Error: n and d must > 0. ')
    xx = regspace(o, d, o + (n - 1)*d)
    ff = zeros(n)

    select case (interp_method)
        case ('nearest')
            call interp_nearest_1d_(l, x, f, n, xx, ff)
        case ('linear')
            call interp_linear_1d_(l, x, f, n, xx, ff)
        case ('cubic')
            call interp_cspline_1d_(l, x, f, n, xx, ff)
        case ('pchip')
            call interp_pchip_1d_(l, x, f, n, xx, ff)
        case ('quintic')
            call interp_quintic_1d_(l, x, f, n, xx, ff)
        case ('mba')
            call interp_mba_1d_(l, x, f, n, xx, ff)
        case ('biharmonic')
            call interp_biharmonic_1d_(l, x, f, n, xx, ff)
        case ('cubic_spline')
            call interp_cubic_spline_1d_(l, x, f, n, xx, ff, method='c2')
        case ('hermite_spline')
            call interp_cubic_spline_1d_(l, x, f, n, xx, ff, method='hermite')
        case ('monotonic_spline')
            call interp_cubic_spline_1d_(l, x, f, n, xx, ff, method='monotonic')
        case default
            call interp_linear_1d_(l, x, f, n, xx, ff)
    end select

end function irreg_to_reg_interp_1d_

function irreg_to_reg_interp_2d_(x, y, f, n, d, o, method) result(ff)

    TT, dimension(:) :: x, y
    TTT, dimension(:) :: f
    integer, dimension(1:2) :: n
    TT, dimension(1:2) :: d, o
    TTT, allocatable, dimension(:, :) :: ff
    character(len=*), optional :: method

    character(len=24) :: interp_method
    integer :: nn, l
    TT, allocatable, dimension(:) :: xx, yy
    TTT, allocatable, dimension(:) :: tf

    if (present(method)) then
        interp_method = method
    else
        interp_method = 'mba'
    end if

    l = size(x)
    call assert(all(n > 0) .and. all(d > 0), ' <irreg_to_reg_interp_2d> Error: All n and d must > 0. ')
    xx = meshgrid(n=n, d=d, o=o, dim=1)
    yy = meshgrid(n=n, d=d, o=o, dim=2)
    nn = product(n)
    tf = zeros(nn)

    select case (interp_method)
        case ('nearest')
            call interp_nearest_2d_(l, x, y, f, nn, xx, yy, tf)
        case ('mba')
            call interp_mba_2d_(l, x, y, f, nn, xx, yy, tf)
        case ('biharmonic')
            call interp_biharmonic_2d_(l, x, y, f, nn, xx, yy, tf)
    end select

    ff = reshape(tf, n)

end function irreg_to_reg_interp_2d_

function irreg_to_reg_interp_3d_(x, y, z, f, n, d, o, method) result(ff)

    TT, dimension(:) :: x, y, z
    TTT, dimension(:) :: f
    integer, dimension(1:3) :: n
    TT, dimension(1:3) :: d, o
    TTT, allocatable, dimension(:, :, :) :: ff
    character(len=*), optional :: method

    character(len=24) :: interp_method
    integer :: nn, l
    TT, allocatable, dimension(:) :: xx, yy, zz
    TTT, allocatable, dimension(:) :: tf

    if (present(method)) then
        interp_method = method
    else
        interp_method = 'mba'
    end if

    l = size(x)
    call assert(all(n > 0) .and. all(d > 0), ' <irreg_to_reg_interp_3d> Error: All n and d must > 0. ')
    xx = meshgrid(n=n, d=d, o=o, dim=1)
    yy = meshgrid(n=n, d=d, o=o, dim=2)
    zz = meshgrid(n=n, d=d, o=o, dim=3)
    nn = product(n)
    tf = zeros(nn)

    select case (interp_method)
        case ('nearest')
            call interp_nearest_3d_(l, x, y, z, f, nn, xx, yy, zz, tf)
        case ('mba')
            call interp_mba_3d_(l, x, y, z, f, nn, xx, yy, zz, tf)
        case ('biharmonic')
            call interp_biharmonic_3d_(l, x, y, z, f, nn, xx, yy, zz, tf)
    end select

    ff = reshape(tf, n)

end function irreg_to_reg_interp_3d_

!
!> Remove NaN values from regularly sampled data
!
function inpaint_1d_(w, method) result(ww)

    TTT, dimension(:), intent(in) :: w
    character(len=*), optional :: method
    TTT, allocatable, dimension(:) :: ww

    character(len=24) :: interp_method
    integer :: n1, n, nn
    integer :: i, l, h
    TT, allocatable, dimension(:) :: x1, xx1
    TTT, allocatable, dimension(:) :: f, ff

    if (present(method)) then
        interp_method = method
    else
        interp_method = 'mba'
    end if

    n1 = size(w)

    ! Find NaN pixels
    n = count(.not. isnan(w))
    nn = n1

    x1 = zeros(n)
    f = zeros(n)
    xx1 = zeros(nn)
    ff = zeros(nn)

    l = 1
    h = 1
    do i = 1, n1
        if (.not. isnan(w(i))) then
            x1(l) = i - 1.0
            f(l) = w(i)
            l = l + 1
        end if
        xx1(h) = i - 1.0
        h = h + 1
    end do

    select case (interp_method)
        case ('nearest')
            call interp_nearest_1d_(n, x1, f, nn, xx1, ff)
        case ('linear')
            call interp_linear_1d_(n, x1, f, nn, xx1, ff)
        case ('sinc')
            call interp_sinc_1d_(n, x1, f, nn, xx1, ff)
        case ('cubic')
            call interp_cspline_1d_(n, x1, f, nn, xx1, ff)
        case ('pchip')
            call interp_pchip_1d_(n, x1, f, nn, xx1, ff)
        case ('quintic')
            call interp_quintic_1d_(n, x1, f, nn, xx1, ff)
        case ('mba')
            call interp_mba_1d_(n, x1, f, nn, xx1, ff)
        case ('biharmonic')
            call interp_biharmonic_1d_(n, x1, f, nn, xx1, ff)
        case ('cubic_spline')
            call interp_cubic_spline_1d_(n, x1, f, nn, xx1, ff, method='c2')
        case ('hermite_spline')
            call interp_cubic_spline_1d_(n, x1, f, nn, xx1, ff, method='hermite')
        case ('monotonic_spline')
            call interp_cubic_spline_1d_(n, x1, f, nn, xx1, ff, method='monotonic')
        case default
            call interp_linear_1d_(n, x1, f, nn, xx1, ff)
    end select

    ww = ff

end function inpaint_1d_

function inpaint_2d_(w, method) result(ww)

    TTT, dimension(:, :), intent(in) :: w
    character(len=*), optional :: method
    TTT, allocatable, dimension(:, :) :: ww

    character(len=24) :: interp_method
    integer :: n1, n2, n, nn
    integer :: i, j, l, h
    TT, allocatable, dimension(:) :: x1, xx1, x2, xx2
    TTT, allocatable, dimension(:) :: f, ff

    if (present(method)) then
        interp_method = method
    else
        interp_method = 'mba'
    end if

    n1 = size(w, 1)
    n2 = size(w, 2)

    ! Find NaN pixels
    n = count(.not. isnan(w))
    nn = size(w)

    x1 = zeros(n)
    x2 = zeros(n)
    f = zeros(n)
    xx1 = zeros(nn)
    xx2 = zeros(nn)
    ff = zeros(nn)

    l = 1
    h = 1
    do j = 1, n2
        do i = 1, n1
            if (.not. isnan(w(i, j))) then
                x1(l) = i - 1.0
                x2(l) = j - 1.0
                f(l) = w(i, j)
                l = l + 1
            end if
            xx1(h) = i - 1.0
            xx2(h) = j - 1.0
            h = h + 1
        end do
    end do

    select case (interp_method)
        case ('nearest')
            call interp_nearest_2d_(n, x1, x2, f, nn, xx1, xx2, ff)
        case ('mba')
            call interp_mba_2d_(n, x1, x2, f, nn, xx1, xx2, ff)
        case ('biharmonic')
            call interp_biharmonic_2d_(n, x1, x2, f, nn, xx1, xx2, ff)
    end select

    ww = reshape(ff, [n1, n2])

end function inpaint_2d_

function inpaint_3d_(w, method) result(ww)

    TTT, dimension(:, :, :), intent(in) :: w
    character(len=*), optional :: method
    TTT, allocatable, dimension(:, :, :) :: ww

    character(len=24) :: interp_method
    integer :: n1, n2, n3, n, nn
    integer :: i, j, k, l, h
    TT, allocatable, dimension(:) :: x1, xx1, x2, xx2, x3, xx3
    TTT, allocatable, dimension(:) :: f, ff

    if (present(method)) then
        interp_method = method
    else
        interp_method = 'mba'
    end if

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)

    ! Find NaN pixels
    n = count(.not. isnan(w))
    nn = size(w)

    x1 = zeros(n)
    x2 = zeros(n)
    x3 = zeros(n)
    f = zeros(n)
    xx1 = zeros(nn)
    xx2 = zeros(nn)
    xx3 = zeros(nn)
    ff = zeros(nn)

    l = 1
    h = 1
    do k = 1, n3
        do j = 1, n2
            do i = 1, n1
                if (.not. isnan(w(i, j, k))) then
                    x1(l) = i - 1.0
                    x2(l) = j - 1.0
                    x3(l) = k - 1.0
                    f(l) = w(i, j, k)
                    l = l + 1
                end if
                xx1(h) = i - 1.0
                xx2(h) = j - 1.0
                xx3(h) = k - 1.0
                h = h + 1
            end do
        end do
    end do

    select case (interp_method)
        case ('nearest')
            call interp_nearest_3d_(n, x1, x2, x3, f, nn, xx1, xx2, xx3, ff)
        case ('mba')
            call interp_mba_3d_(n, x1, x2, x3, f, nn, xx1, xx2, xx3, ff)
        case ('biharmonic')
            call interp_biharmonic_3d_(n, x1, x2, x3, f, nn, xx1, xx2, xx3, ff)
    end select

    ww = reshape(ff, [n1, n2, n3])

end function inpaint_3d_

function point_interp_linear_1d_(x, f, xx) result(ff)

    TT, dimension(:) :: x
    TT, dimension(:) :: f
    TT :: xx, ff

    TT, dimension(1:2) :: w
    integer :: i

    do i = 1, 2
        w(i) = abs(x(3 - i) - xx)/(x(2) - x(1))
    end do

    ff = sum(w*f)

end function point_interp_linear_1d_

function point_interp_linear_2d_(x, y, f, xx, yy) result(ff)

    TT, dimension(:) :: x, y
    TT, dimension(:, :) :: f
    TT :: xx, yy, ff

    TT, dimension(1:2, 1:2) :: w
    integer :: i, j

    do i = 1, 2
        do j = 1, 2
            w(i, j) = abs((x(3 - i) - xx)*(y(3 - j) - yy))/((x(2) - x(1))*(y(2) - y(1)))
        end do
    end do

    ff = sum(w*f)

end function point_interp_linear_2d_

function point_interp_linear_3d_(x, y, z, f, xx, yy, zz) result(ff)

    TT, dimension(:) :: x, y, z
    TT, dimension(:, :, :) :: f
    TT :: xx, yy, zz, ff

    TT, dimension(1:2, 1:2, 1:2) :: w
    integer :: i, j, k

    do i = 1, 2
        do j = 1, 2
            do k = 1, 2
                w(i, j, k) = abs((x(3 - i) - xx)*(y(3 - j) - yy)*(z(3 - k) - zz)) &
                    /((x(2) - x(1))*(y(2) - y(1))*(z(2) - z(1)))
            end do
        end do
    end do

    ff = sum(w*f)

end function point_interp_linear_3d_

function point_interp_barycentric_2d_(v1, v2, v3, f, p) result(ff)

    TT, dimension(:) :: v1, v2, v3, f, p
    TT :: ff

    TT, dimension(1:3) :: w
    TT :: a123

    a123 = area(v1, v2, v3)

    if (abs(a123) < 1e-6) then
        print *, ' <point_interp_barycentric_2d> Error: The points form a degenerate triangle! '
        stop
    end if

    ! Note the circulant order
    w(1) = area(p, v2, v3)/a123
    w(2) = area(v1, p, v3)/a123
    w(3) = area(v1, v2, p)/a123

    if (any(w < 0)) then
        !        print *, ' <point_interp_barycentric_2d> Warning: The point is outside of the triangle! '
        ff = 0
        return
    end if

    ff = sum(w*f)

contains

    function area(v1, v2, v3) result(a)

        TT, dimension(:) :: v1, v2, v3
        TT :: a

        a = det(transpose(reshape([v2 - v1, v3 - v1], [2, 2])))

    end function area

end function point_interp_barycentric_2d_

function point_interp_barycentric_3d_(v1, v2, v3, v4, f, p) result(ff)

    TT, dimension(:) :: v1, v2, v3, v4, f, p
    TT :: ff

    TT, dimension(1:4) :: w
    TT :: vol1234

    ! Calculate the determinant for the tetrahedron volume
    vol1234 = vol(v1, v2, v3, v4)

    if (abs(vol1234) < 1e-6) then
        print *, ' <point_interp_barycentric_3d> Error: The points form a degenerate tetrahedron. '
        stop
    end if

    ! Calculate the determinants for sub-tetrahedra
    ! Note the circulant order
    w(1) = vol(p, v2, v3, v4)/vol1234
    w(2) = vol(v1, p, v3, v4)/vol1234
    w(3) = vol(v1, v2, p, v4)/vol1234
    w(4) = vol(v1, v2, v3, p)/vol1234

    if (any(w < 0)) then
        !        print *, ' <point_interp_barycentric_3d> Warning: The point is outside of the tetrahedron! '
        ff = 0
        return
    end if

    ff = sum(w*f)

contains

    function vol(v1, v2, v3, v4) result(v)

        TT, dimension(:) :: v1, v2, v3, v4
        TT :: v

        v = det(transpose(reshape([v2 - v1, v3 - v1, v4 - v1], [3, 3])))

    end function vol

end function point_interp_barycentric_3d_

#undef T
#undef TT
#undef TTT
#undef nTT
#undef nTTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef interp_nearest_1d_
#undef interp_nearest_2d_
#undef interp_nearest_3d_
#undef interp_linear_1d_
#undef interp_cubic_spline_1d_
#undef interp_biharmonic_1d_
#undef interp_biharmonic_2d_
#undef interp_biharmonic_3d_
#undef interp_sinc_1d_
#undef meshgrid_
#undef meshgrid_1d_
#undef meshgrid_2d_
#undef meshgrid_3d_
#undef resample_1d_
#undef resample_2d_
#undef resample_3d_
#undef interp_to_1d_
#undef interp_to_2d_
#undef interp_to_3d_
#undef interp_like_1d_
#undef interp_like_2d_
#undef interp_like_3d_
#undef reg_to_reg_interp_1d_
#undef reg_to_reg_interp_2d_
#undef reg_to_reg_interp_3d_
#undef irreg_to_irreg_interp_1d_
#undef irreg_to_irreg_interp_2d_
#undef irreg_to_irreg_interp_3d_
#undef irreg_to_reg_interp_1d_
#undef irreg_to_reg_interp_2d_
#undef irreg_to_reg_interp_3d_
#undef inpaint_1d_
#undef inpaint_2d_
#undef inpaint_3d_

#undef point_interp_linear_1d_
#undef point_interp_linear_2d_
#undef point_interp_linear_3d_
#undef point_interp_barycentric_2d_
#undef point_interp_barycentric_3d_

! external
#undef interp_cspline_1d_
#undef interp_pchip_1d_
#undef interp_quintic_1d_
#undef interp_mba_1d_
#undef interp_mba_2d_
#undef interp_mba_3d_
