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

#define gaussian_curvature_      CONCAT(gaussian_curvature, T)
#define spherical_to_cartesian_      CONCAT(spherical_to_cartesian, T)
#define cartesian_to_spherical_      CONCAT(cartesian_to_spherical, T)
#define convexhull_2d_      CONCAT(convexhull_2d, T)
#define convexhull_3d_      CONCAT(convexhull_3d, T)
#define distance_between_line_segments_      CONCAT(distance_between_line_segments, T)
#define point_distance_to_line_      CONCAT(point_distance_to_line, T)
#define point_signed_distance_to_line_      CONCAT(point_signed_distance_to_line, T)
#define point_distance_to_line_segement_      CONCAT(point_distance_to_line_segement, T)
#define polygon_contains_point_      CONCAT(polygon_contains_point, T)
#define fit_curve_      CONCAT(fit_curve, T)
#define fit_surface_      CONCAT(fit_surface, T)

!
!> Compute Gaussian curvate for a 2D image
!
function gaussian_curvature_(w) result(c)

    TT, dimension(:, :), intent(in) :: w
    TT, allocatable, dimension(:, :) :: c

    TT, allocatable, dimension(:, :) :: hx, hy, hxx, hyy, hxy

    hx = deriv(w, dim=1, order=1, method='center')
    hy = deriv(w, dim=2, order=1, method='center')
    hxx = deriv(w, dim=1, order=2, method='center')
    hyy = deriv(w, dim=2, order=2, method='center')
    hxy = deriv(hx, dim=2, order=1, method='center')
    c = (hxx*hyy - hxy**2)/(1 + hx**2 + hy**2)**2

end function

!
!> Spherical coordinates to Cartesian coordinates
!
function spherical_to_cartesian_(r, polar, azimuth) result(p)

    TT :: r, polar, azimuth
    TT, allocatable, dimension(:) :: p

    p = zeros(3)

    p(1) = r*sin(polar)*cos(azimuth)
    p(2) = r*sin(polar)*sin(azimuth)
    p(3) = r*cos(polar)

end function

!
!> Cartesian coordinates to spherical coordinates
!
function cartesian_to_spherical_(x, y, z) result(s)

    TT :: x, y, z
    TT, allocatable, dimension(:) :: s

    TT :: r, polar, azimuth

    r = norm2([x, y, z])
    polar = acos(z/r)
    azimuth = sign(1.0, y)*acos(x/r)

    s = [r, polar, azimuth]

end function

!
!> Get indicies of convex hull points in 2D
!
function convexhull_2d_(x, y) result(f)

    use convex_hull_2d_mod

    TT, dimension(:) :: x, y
    integer, allocatable, dimension(:) :: f

    integer :: n
    double precision, allocatable, dimension(:, :) :: points

    call assert(size(x) == size(y), ' <convexhull_2d> Error: size(x) must = size(y)')

    n = size(x)
    f = zeros(n)

    points = zeros(2, n)
    points(1, :) = x
    points(2, :) = y

    call convex_hull_2d(points, f)

end function

!
!> Get indicies of convex hull points in 3D
!
function convexhull_3d_(x, y, z) result(f)

    use convex_hull_3d_mod

    TT, dimension(:) :: x, y, z
    integer, allocatable, dimension(:, :) :: f

    integer :: n
    type(face3d_t), allocatable :: faces(:)
    double precision, allocatable, dimension(:, :) :: points
    integer :: i, stat

    call assert(size(x) == size(y) .and. size(x) == size(z), &
        ' <convexhull_3d> Error: size(x) must = size(y) and = size(z)')

    n = size(x)

    points = zeros(3, n)
    points(1, :) = x
    points(2, :) = y
    points(3, :) = z

    call convex_hull_3d(points, faces, stat=stat)

    f = zeros(3, size(faces))
    do i = 1, size(faces)
        f(:, i) = faces(i)%v
    end do

end function

!
!  Compute shortest distance between two line segments
!
subroutine distance_between_line_segments_(point1s, point1e, point2s, point2e, dist, p)

    TT, dimension(:), intent(in) :: point1s, point1e, point2s, point2e
    TT, intent(out) :: dist
    TT, allocatable, dimension(:, :), intent(out), optional :: p

    TT, allocatable, dimension(:) :: p1, p2, p12
    TT :: d1, d2, s1, s2, r, den
    TT :: u, t, uf

    p1 = point1e - point1s
    p2 = point2e - point2s
    p12 = point2s - point1s
    D1 = sum(p1*p1)
    D2 = sum(p2*p2)
    S1 = sum(p1*p12)
    S2 = sum(p2*p12)
    R = sum(p1*p2)
    den = D1*D2 - R**2

    if (D1 == 0 .or. D2 == 0) then  ! if one of the segments is a point

        if (D1 /= 0) then ! if line1 is a segment and line2 is a point
            u = 0
            t = S1/D1
            t = clip(t, real(0.0, fp), real(1.0, fp))

        else if (D2 /= 0) then  ! if line2 is a segment and line 1 is a point

            t = 0
            u = -S2/D2
            u = clip(u, real(0.0, fp), real(1.0, fp))

        else                !    % both segments are points
            t = 0
            u = 0
        end if

    else if (den == 0) then    !      % if lines are parallel
        t = 0
        u = -S2/D2

        uf = clip(u, real(0.0, fp), real(1.0, fp))

        if (uf /= u) then
            t = (uf*R + S1)/D1
            t = clip(t, real(0.0, fp), real(1.0, fp))

            u = uf
        end if
    else      !                  % general case

        t = (S1*D2 - S2*R)/den

        t = clip(t, real(0.0, fp), real(1.0, fp))

        u = (t*R - S2)/D2
        uf = clip(u, real(0.0, fp), real(1.0, fp))

        if (uf /= u) then
            t = (uf*R + S1)/D1
            t = clip(t, real(0.0, fp), real(1.0, fp))

            u = uf
        end if
    end if

    dist = norm2(p1*t - p2*u - p12)

    if (present(p)) then
        p = zeros(2, size(point1s))
        p(1, :) = point1s + p1*t
        p(2, :) = point2s + p2*u
    end if

end subroutine

!
!  Compute shortest distance between a line and a point
!
!  Author:
!       John Burkardt, 2005
!       K.G., 2024
!
function point_distance_to_line_(p1, p2, p) result(dist)

    TT, dimension(:) :: p1, p2, p
    TT :: dist

    TT :: bot, dot, t
    TT, allocatable, dimension(:) :: pn

    call assert(size(p1) == size(p2) .and. size(p1) == size(p), &
        ' <point_distance_to_line> Error: size(p1) must = size(p2) = size(p)')

    if (all(p1 == p2)) then
        ! Line is a point

        pn = p1

    else
        !
        !  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
        !
        !  (P-P1) dot (P2-P1) / Norm(P-P1)^2 = normalized coordinate T
        !  of the projection of (P-P1) onto (P2-P1).
        !
        dot = sum((p - p1)*(p2 - p1))
        bot = sum((p2 - p1)**2)
        t = dot/bot
        pn = p1 + t*(p2 - p1)

    end if

    dist = norm2(p - pn)

end function

!
!  Compute signed shortest distance between a line and a point
!
!  Author:
!       John Burkardt, 2005
!       K.G., 2024
!
function point_signed_distance_to_line_(p1, p2, p) result(dist_signed)

    TT, dimension(1:2) :: p1, p2, p
    TT :: dist_signed

    TT :: a, b, c

    if (all(p1 == p2)) then

        dist_signed = norm2(p1 - p1)

    else
        !
        !  Convert the explicit line to the implicit form A * P(1) + B * P(2) + C = 0.
        !  This makes the computation of the signed distance to (X,Y) easy.
        !
        a = p2(2) - p1(2)
        b = p1(1) - p2(1)
        c = p2(1)*p1(2) - p1(1)*p2(2)

        dist_signed = (a*p(1) + b*p(2) + c)/sqrt(a**2 + b**2)

    end if

end function

!
!> Compute point distance to a line segement and the nearest point on the line segment
!>
!> Original author:
!>      John Burkardt, 2005
!
!  Parameters:
!
!    Input, real (kind = real_kind ) P1(2), P2(2), the endpoints of the line segment.
!
!    Input, real (kind = real_kind ) P(2), the point whose nearest neighbor
!    on the line segment is to be determined.
!
!    Output, real (kind = real_kind ) PN(2), the point on the line segment which is
!    nearest the point P.
!
!    Output, real (kind = real_kind ) DIST, the distance from the point to the
!    nearest point on the line segment.
!
!    Output, real (kind = real_kind ) T, the relative position of the point PN
!    to the points P1 and P2.
!
subroutine point_distance_to_line_segement_(p1, p2, p, pn, dist, t)

    TT, dimension(:), intent(in) :: p, p1, p2
    TT, dimension(:), intent(out) :: pn
    TT, intent(out) :: t, dist

    !        integer, parameter :: ndim = 2
    TT :: bot

    if (all(p1 == p2)) then
        ! Line is a point

        t = 0.0d0

    else

        bot = sum((p2 - p1)**2)

        t = sum((p - p1)*(p2 - p1))/bot

        t = max(t, 0.0d0)
        t = min(t, 1.0d0)

    end if

    pn = p1 + t*(p2 - p1)

    dist = norm2(p - pn)

end subroutine

!
!> Check if a point is in polygon
!
!> @note Given a polygonal line connecting the vertices (x(i), y(i)) (i = 1,...,n)
!> taken in this order.  it is assumed that the polygonal path is a loop,
!> where (x(n), y(n)) = (x(1), y(1)) or there is an arc from (x(n), y(n)) to
!> (x(1), y(1)). The polygon may cross itself any number of times.
!>
!>    (x0, y0) is an arbitrary point
!>    l = -1   if (x0,y0) is outside the polygonal path
!>    l =  0   if (x0,y0) lies on the polygonal path
!>    l =  1   if (x0,y0) is inside the polygonal path
!>    m = 0 if (x0,y0) is on or outside the path.  if (x0,y0) is inside the
!>    path then m is the winding number of the path around the point (x0,y0).
!>
!> @author fortran 66 version by A.H. Morris
!> @author converted to elf90 compatibility by Alan Miller, 15 February 1997
!> @author Further refined by K.G.
!
function polygon_contains_point_(polygon, p) result(l)

    TT, dimension(:, :) :: polygon
    TT, dimension(:) :: p

    TT, allocatable, dimension(:) :: x, y
    TT :: x0, y0
    integer :: l, n, m, i, n0
    TT :: angle, eps, pi, pi2, sum, theta, theta1, thetai, tol, u, v

    x = polygon(:, 1)
    y = polygon(:, 2)
    x0 = p(1)
    y0 = p(2)

    ! Small number
    eps = epsilon(1.0)

    ! start checking
    n = size(x)
    n0 = n
    if (x(1) == x(n) .and. y(1) == y(n)) then
        n0 = n - 1
    end if
    pi = const_pi
    pi2 = 2.0*const_pi
    tol = 4.0*eps*const_pi
    l = -1
    m = 0

    u = x(1) - x0
    v = y(1) - y0
    if (u == 0.0 .and. v == 0.0) then
        l = 0
        return
    end if
    if (n0 < 2) then
        return
    end if
    theta1 = atan2(v, u)

    sum = 0.0
    theta = theta1
    do i = 2, n0
        u = x(i) - x0
        v = y(i) - y0
        if (u == 0.0 .and. v == 0.0) then
            l = 0
            return
        end if
        thetai = atan2(v, u)

        angle = abs(thetai - theta)
        if (abs(angle - pi) < tol) then
            l = 0
            return
        end if
        if (angle > pi) then
            angle = angle - pi2
        end if
        if (theta > thetai) then
            angle = -angle
        end if
        sum = sum + angle
        theta = thetai
    end do

    angle = abs(theta1 - theta)
    if (abs(angle - pi) < tol) then
        l = 0
        return
    end if
    if (angle > pi) then
        angle = angle - pi2
    end if
    if (theta > theta1) then
        angle = -angle
    end if
    sum = sum + angle

    ! sum = 2*const_pi*m where m is the winding number
    m = abs(sum)/pi2 + 0.2
    if (m == 0) then
        return
    end if
    l = 1
    if (sum < 0.0) then
        m = -m
    end if

end function

!
!> Find a best-fit, continuous, smooth curve for a 2D point cloud using polynomial or LOWESS
!> The return are (1) the points projects onto the fitted curve, and
!> (2) a 2D pixel image where the fit curve is 1 while the
!> background is 0. To ensure sufficient smoothness, the sampling interval d
!> should be properly small.
!
subroutine fit_curve_(p1, p2, p1f, p2f, m, n, d, o, method, smooth, order, range)

    TT, dimension(:), intent(in) :: p1, p2
    integer, dimension(1:2), intent(in), optional :: n
    TT, dimension(1:2), intent(in), optional :: d, o
    character(len=*), intent(in), optional :: method
    TT, intent(in), optional :: smooth
    integer, intent(in), optional :: order
    character(len=*), intent(in), optional :: range
    TT, allocatable, dimension(:), intent(out) :: p1f, p2f
    TT, allocatable, dimension(:, :), intent(out), optional :: m

    integer :: i, np
    type(polyfit1) :: ff
    TT, allocatable :: p(:, :), pr(:, :), xx(:), yy(:)
    TT, allocatable :: eigenval(:), eigenvec(:, :), mmt(:, :)
    TT, allocatable :: r(:, :)
    TT :: xmin, xmax
    TT :: c1, c2
    integer :: n1, n2
    TT :: d1, d2
    TT :: o1, o2
    integer :: ind, jnd
    character(len=32) :: fit_method
    TT :: fit_smooth
    integer :: fit_order
    TT :: xextra

    if (present(method)) then
        fit_method = method
    else
        fit_method = 'polynomial'
    end if
    if (present(smooth)) then
        fit_smooth = smooth
    else
        fit_smooth = 0.5
    end if
    if (present(order)) then
        fit_order = order
    else
        fit_order = 2
    end if

    ! Find center
    np = size(p1)
    p = zeros(np, 2)
    c1 = mean(p1)
    c2 = mean(p2)

    ! Compute principal direction
    p(:, 1) = p1 - c1
    p(:, 2) = p2 - c2
    mmt = matx(transpose(p), p)
    call eigen(mmt, eigenval, eigenvec)

    ! Rotate data to so that the (p1, p2) -> (x2, x1)
    ! -- just a choice, could be (x1, x2) but the fitting should be
    ! changed accordingly
    !
    ! The change of basis problem is
    !       R eigenvec = axis
    ! therefore,
    !       R eigenvec eigenvec^-1 = axis eigenvec^-1
    ! Because eigenvec matrix could be singular, here I use lsqsolve
    ! to avoid singularity.
    r = matx(flip(eye(2, real(1.0, fp)), axis=[2]), lsqsolve(eigenvec, eye(2, real(1.0, fp))))
    pr = zeros(np, 2)
    !$omp parallel do private(i)
    do i = 1, np
        pr(i, :) = matx(r, p(i, :))
    end do
    !$omp end parallel do

    ! The points on the fitted surface
    xx = pr(:, 2)
    yy = pr(:, 1)
    np = size(xx)
    select case (fit_method)
        case ('polynomial')
            ff%order = fit_order
            call ff%build(TTT(xx), TTT(yy), verbose=.false.)
            do i = 1, np
                yy(i) = ff%eval(TTT(xx(i)))
            end do
        case ('lowess')
            yy = lowess_filt(xx, yy, xx, fit_order, fit_smooth)
    end select

    ! Rotate back
    p(:, 1) = yy
    p(:, 2) = xx
    !$omp parallel do private(i)
    do i = 1, np
        p(i, :) = matx(transpose(r), p(i, :)) + [c1, c2]
    end do
    !$omp end parallel do

    p1f = p(:, 1)
    p2f = p(:, 2)

    ! If output a regular-grid array
    if (present(m)) then

        call assert(present(n) .and. present(d) .and. present(o), ' <fit_curve> Error: n, d, o must be set. ')

        n1 = n(1)
        n2 = n(2)
        d1 = d(1)
        d2 = d(2)
        o1 = o(1)
        o2 = o(2)

        ! Fit by polynomial or LOWESS
        xextra = 10*d2
        if (present(range)) then
            if (range == 'strict') then
                xextra = 0.0
            end if
        end if
        xmin = minval(pr(:, 2)) - xextra
        xmax = maxval(pr(:, 2)) + xextra
        xx = regspace(xmin, 0.5*d2, xmax)
        yy = zeros_like(xx)
        np = size(xx)
        select case (fit_method)
            case ('polynomial')
                do i = 1, np
                    yy(i) = ff%eval(TTT(xx(i)))
                end do
            case ('lowess')
                yy = lowess_filt(pr(:, 2), pr(:, 1), xx, fit_order, fit_smooth)
        end select

        ! Rotate back
        pr = zeros(np, 2)
        pr(:, 1) = yy
        pr(:, 2) = xx
        !$omp parallel do private(i)
        do i = 1, np
            pr(i, :) = matx(transpose(r), pr(i, :)) + [c1, c2]
        end do
        !$omp end parallel do

        m = zeros(n1, n2)
        !$omp parallel do private(i, ind, jnd)
        do i = 1, np
            ind = nint((pr(i, 1) - o1)/d1) + 1
            jnd = nint((pr(i, 2) - o2)/d2) + 1
            if (ind >= 1 .and. ind <= n1 .and. jnd >= 1 .and. jnd <= n2) then
                m(ind, jnd) = 1.0
            end if
        end do
        !$omp end parallel do

    end if

end subroutine

!
!> Find a best-fit, continuous, smooth curve for a 3D point cloud using polynomial or LOWESS
!> The returns are (1) the points projected on the fitted surface, and
!> (2) a 3D pixel image where the fit curve is 1 while the
!> background is 0. To ensure sufficient smoothness, the sampling interval d
!> should be properly small.
!
subroutine fit_surface_(p1, p2, p3, p1f, p2f, p3f, m, n, d, o, method, smooth, order, range)

    TT, dimension(:), intent(in) :: p1, p2, p3
    integer, dimension(1:3), intent(in), optional :: n
    TT, dimension(1:3), intent(in), optional :: d, o
    character(len=*), intent(in), optional :: method
    TT, intent(in), optional :: smooth
    integer, intent(in), optional :: order
    character(len=*), intent(in), optional :: range
    TT, allocatable, dimension(:), intent(out) :: p1f, p2f, p3f
    TT, allocatable, dimension(:, :, :), intent(out), optional :: m

    integer :: i, np, m3, m2
    type(polyfit2) :: ff
    TT, allocatable :: p(:, :), pr(:, :), xx(:), yy(:), zz(:)
    TT, allocatable :: eigenval(:), eigenvec(:, :), mmt(:, :)
    TT, allocatable :: r(:, :)
    TT :: xmin, xmax, ymin, ymax
    TT :: c1, c2, c3
    integer :: n1, n2, n3
    TT :: d1, d2, d3
    TT :: o1, o2, o3
    integer :: ind, jnd, knd
    integer :: l
    integer, allocatable, dimension(:) :: qindex
    TT, allocatable, dimension(:, :) :: pcx
    character(len=32) :: fit_method
    TT :: fit_smooth
    integer :: fit_order

    if (present(method)) then
        fit_method = method
    else
        fit_method = 'polynomial'
    end if
    if (present(smooth)) then
        fit_smooth = smooth
    else
        fit_smooth = 0.5
    end if
    if (present(order)) then
        fit_order = order
    else
        fit_order = 2
    end if

    ! Find center
    np = size(p1)
    p = zeros(np, 3)
    c1 = mean(p1)
    c2 = mean(p2)
    c3 = mean(p3)

    ! Compute principal direction
    p(:, 1) = p1 - c1
    p(:, 2) = p2 - c2
    p(:, 3) = p3 - c3
    mmt = matx(transpose(p), p)
    call eigen(mmt, eigenval, eigenvec)

    ! Rotate data to so that the (p1, p2, p3) -> (x3, x2, x1)
    ! -- just a choice, could be (x1, x2, x3) but the fitting should be
    ! changed accordingly
    r = matx(flip(eye(3, real(1.0, fp)), axis=[2]), lsqsolve(eigenvec, eye(3, real(1.0, fp))))
    pr = zeros(np, 3)
    !$omp parallel do private(i)
    do i = 1, np
        pr(i, :) = matx(r, p(i, :))
    end do
    !$omp end parallel do

    ! The points on the fitted surface
    xx = pr(:, 3)
    yy = pr(:, 2)
    zz = pr(:, 1)
    select case (fit_method)
        case ('polynomial')
            ff%order = fit_order
            call ff%build(TTT(xx), TTT(yy), TTT(zz), verbose=.false.)
            do i = 1, np
                zz(i) = ff%eval(TTT(xx(i)), TTT(yy(i)))
            end do
        case ('lowess')
            zz = lowess_filt(xx, yy, zz, xx, yy, fit_order, [fit_smooth, fit_smooth])
    end select

    ! Rotate back
    p(:, 1) = zz
    p(:, 2) = yy
    p(:, 3) = xx
    !$omp parallel do private(i)
    do i = 1, np
        p(i, :) = matx(transpose(r), p(i, :)) + [c1, c2, c3]
    end do
    !$omp end parallel do

    p1f = p(:, 1)
    p2f = p(:, 2)
    p3f = p(:, 3)

    ! If output a regular-grid array
    if (present(m)) then

        n1 = n(1)
        n2 = n(2)
        n3 = n(3)
        d1 = d(1)
        d2 = d(2)
        d3 = d(3)
        o1 = o(1)
        o2 = o(2)
        o3 = o(3)

        call assert(present(n) .and. present(d) .and. present(o), ' <fit_curve> Error: n, d, o must be set. ')

        ! Fit surface by polynomial or LOWESS
        xmin = minval(pr(:, 3)) - 10*d3
        xmax = maxval(pr(:, 3)) + 10*d3
        ymin = minval(pr(:, 2)) - 10*d2
        ymax = maxval(pr(:, 2)) + 10*d2
        xx = regspace(xmin, 0.5*d3, xmax)
        yy = regspace(ymin, 0.5*d2, ymax)
        m3 = size(xx)
        m2 = size(yy)
        xx = meshgrid([m2, m3], [0.5*d2, 0.5*d3], [ymin, xmin], dim=2)
        yy = meshgrid([m2, m3], [0.5*d2, 0.5*d3], [ymin, xmin], dim=1)
        zz = zeros_like(xx)
        np = size(xx)

        ! If range = strict then use convex hull to get the exact range of the points projected
        ! on the principal plan (x1, x2); the resulting plane may be irregular
        if (present(range)) then
            if (range == 'strict') then

                qindex = convexhull(pr(:, 3), pr(:, 2))
                pcx = pr(qindex, [3, 2])

                qindex = zeros(np)
                l = 1
                do i = 1, np
                    if (polygon_contains_point(pcx, [xx(i), yy(i)]) /= -1) then
                        qindex(l) = i
                        l = l + 1
                    end if
                end do

                qindex = qindex(1:l - 1)

                xx = xx(qindex)
                yy = yy(qindex)
                zz = zz(qindex)
                np = size(xx)

            end if
        end if

        select case (fit_method)
            case ('polynomial')
                do i = 1, np
                    zz(i) = ff%eval(TTT(xx(i)), TTT(yy(i)))
                end do
            case ('lowess')
                zz = lowess_filt(pr(:, 3), pr(:, 2), pr(:, 1), xx, yy, fit_order, [fit_smooth, fit_smooth])
        end select

        ! Rotate back
        pr = zeros(np, 3)
        pr(:, 1) = zz
        pr(:, 2) = yy
        pr(:, 3) = xx
        !$omp parallel do private(i)
        do i = 1, np
            pr(i, :) = matx(transpose(r), pr(i, :)) + [c1, c2, c3]
        end do
        !$omp end parallel do

        m = zeros(n1, n2, n3)
        !$omp parallel do private(i, ind, jnd, knd)
        do i = 1, np
            ind = nint((pr(i, 1) - o1)/d1) + 1
            jnd = nint((pr(i, 2) - o2)/d2) + 1
            knd = nint((pr(i, 3) - o3)/d3) + 1
            if (ind >= 1 .and. ind <= n1 .and. jnd >= 1 .and. jnd <= n2 .and. knd >= 1 .and. knd <= n3) then
                m(ind, jnd, knd) = 1.0
            end if
        end do
        !$omp end parallel do

    end if

end subroutine

#undef T
#undef TT
#undef TTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef gaussian_curvature_
#undef spherical_to_cartesian_
#undef cartesian_to_spherical_
#undef convexhull_2d_
#undef convexhull_3d_
#undef distance_between_line_segments_
#undef point_distance_to_line_
#undef point_signed_distance_to_line_
#undef point_distance_to_line_segement_
#undef polygon_contains_point_
#undef fit_curve_
#undef fit_surface_
