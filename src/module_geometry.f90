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


module libflit_geometry

    use libflit_linear_algebra
    use libflit_constants
    use libflit_utility
    use libflit_array
    use libflit_calculus
    use libflit_error
    use libflit_fit
    use libflit_lowessfilt
    use libflit_statistics
    use libflit_interp
    use libflit_array_operation

    implicit none

    ! For using external convex hull codes
    interface

        subroutine convhull_2d(n, x, y, nf, f) bind(c, name='convhull_2d')

            use iso_c_binding, only: c_int

            integer(kind=c_int), value :: n
            real, dimension(*), intent(in) :: x, y
            integer(kind=c_int), intent(out) :: nf
            integer(kind=c_int), dimension(*), intent(out) :: f

        end subroutine convhull_2d

        subroutine convhull_3d(n, x, y, z, nf, f) bind(c, name='convhull_3d')

            use iso_c_binding, only: c_int

            integer(kind=c_int), value :: n
            real, dimension(*), intent(in) :: x, y, z
            integer(kind=c_int), intent(out) :: nf
            integer(kind=c_int), dimension(*), intent(out) :: f

        end subroutine convhull_3d

    end interface

    interface rotation_matrix
        module procedure :: rotation_matrix_along_arbitrary
        module procedure :: rotation_matrix_along_axis
        module procedure :: rotation_matrix_2d
        module procedure :: rotation_matrix_3d_spherical
        module procedure :: rotation_matrix_3d_euler
        module procedure :: rotation_matrix_from_to
    end interface rotation_matrix

    interface rotate_point
        module procedure :: rotate_point_2d
        module procedure :: rotate_points_2d
        module procedure :: rotate_point_3d
        module procedure :: rotate_points_3d
    end interface rotate_point

    interface convexhull
        module procedure :: convexhull_2d
        module procedure :: convexhull_3d
    end interface convexhull

    private
    public :: convexhull
    public :: polygon_contains_point
    public :: rotate_point
    public :: rotation_matrix
    public :: distance_between_line_segments
    public :: point_distance_to_line_segement
    public :: point_distance_to_line
    public :: point_signed_distance_to_line
    public :: gaussian_curvature
    public :: spherical_to_cartesian
    public :: cartesian_to_spherical
    public :: fit_curve
    public :: fit_surface

contains

    !
    !> Spherical coordinates to Cartesian coordinates
    !
    function spherical_to_cartesian(r, polar, azimuth) result(p)

        real :: r, polar, azimuth
        real, allocatable, dimension(:) :: p

        p = zeros(3)

        p(1) = r*sin(polar)*cos(azimuth)
        p(2) = r*sin(polar)*sin(azimuth)
        p(3) = r*cos(polar)

    end function spherical_to_cartesian

    !
    !> Cartesian coordinates to spherical coordinates
    !
    function cartesian_to_spherical(x, y, z) result(s)

        real :: x, y, z
        real, allocatable, dimension(:) :: s

        real :: r, polar, azimuth

        r = norm2([x, y, z])
        polar = acos(z/r)
        azimuth = sign(1.0, y)*acos(x/r)

        s = [r, polar, azimuth]

    end function cartesian_to_spherical

    !
    !> Get indicies of convex hull points in 2D
    !
    function convexhull_2d(x, y) result(f)

        real, dimension(:) :: x, y
        integer, allocatable, dimension(:) :: f

        integer :: n, nf

        call assert(size(x) == size(y), ' <convexhull_2d> Error: size(x) must = size(y)')

        n = size(x)
        f = zeros(n)

        call convhull_2d(n, x, y, nf, f)

        ! Original output indicies are 0-based
        f = f(1:nf) + 1

    end function convexhull_2d

    !
    !> Get indicies of convex hull points in 3D
    !
    function convexhull_3d(x, y, z) result(f)

        real, dimension(:) :: x, y, z
        integer, allocatable, dimension(:, :) :: f

        integer :: n, nf
        integer, allocatable, dimension(:) ::  ft

        call assert(size(x) == size(y) .and. size(x) == size(z), &
            ' <convexhull_3d> Error: size(x) must = size(y) and = size(z)')

        n = size(x)
        ft = zeros(3*n)

        call convhull_3d(n, x, y, z, nf, ft)

        ! Original output indicies are 0-based
        f = transpose(reshape(ft(1:3*nf), [3, nf])) + 1

    end function convexhull_3d

    subroutine distance_between_line_segments(point1s, point1e, point2s, point2e, dist, p)

        real, dimension(:), intent(in) :: point1s, point1e, point2s, point2e
        real, intent(out) :: dist
        real, allocatable, dimension(:, :), intent(out), optional :: p

        real, allocatable, dimension(:) :: p1, p2, p12
        real :: d1, d2, s1, s2, r, den
        real :: u, t, uf

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
                t = clip(t, 0.0, 1.0)

            else if (D2 /= 0) then  ! if line2 is a segment and line 1 is a point

                t = 0
                u = -S2/D2
                u = clip(u, 0.0, 1.0)

            else                !    % both segments are points
                t = 0
                u = 0
            end if

        else if (den == 0) then    !      % if lines are parallel
            t = 0
            u = -S2/D2

            uf = clip(u, 0.0, 1.0)

            if (uf /= u) then
                t = (uf*R + S1)/D1
                t = clip(t, 0.0, 1.0)

                u = uf
            end if
        else      !                  % general case

            t = (S1*D2 - S2*R)/den

            t = clip(t, 0.0, 1.0)

            u = (t*R - S2)/D2
            uf = clip(u, 0.0, 1.0)

            if (uf /= u) then
                t = (uf*R + S1)/D1
                t = clip(t, 0.0, 1.0)

                u = uf
            end if
        end if

        dist = norm2(p1*t - p2*u - p12)

        if (present(p)) then
            p = zeros(2, size(point1s))
            p(1, :) = point1s + p1*t
            p(2, :) = point2s + p2*u
        end if

    end subroutine distance_between_line_segments

    function gaussian_curvature(w) result(c)

        real, dimension(:, :), intent(in) :: w
        real, allocatable, dimension(:, :) :: c

        real, allocatable, dimension(:, :) :: hx, hy, hxx, hyy, hxy

        hx = deriv(w, dim=1, order=1, method='center')
        hy = deriv(w, dim=2, order=1, method='center')
        hxx = deriv(w, dim=1, order=2, method='center')
        hyy = deriv(w, dim=2, order=2, method='center')
        hxy = deriv(hx, dim=2, order=1, method='center')
        c = (hxx*hyy - hxy**2)/(1 + hx**2 + hy**2)**2

    end function gaussian_curvature

    !
    !> 3D rotation matrix along an arbitrary direction
    !
    pure function rotation_matrix_along_arbitrary(theta, direction) result(r)

        !> Counter-clockwise rotation angle in radian
        real, intent(in) :: theta
        !> 1D array of size 3 to represent the rotation axis
        real, dimension(1:3), intent(in) :: direction
        !> 3x3 rotation matrix
        real, allocatable, dimension(:, :) :: r

        real :: s, c, t, rn(1:3), ux, uy, uz

        ! Normalized directional vector
        rn = direction/norm2(direction)
        ux = rn(1)
        uy = rn(2)
        uz = rn(3)

        ! Replacement
        c = cos(theta)
        s = sin(theta)
        t = 1.0 - c

        ! Rotation matrix
        allocate (r(1:3, 1:3))
        r(1, :) = [t*ux**2 + c, t*ux*uy - s*uz, t*ux*uz + s*uy]
        r(2, :) = [t*ux*uy + s*uz, t*uy**2 + c, t*uy*uz - s*ux]
        r(3, :) = [t*ux*uz - s*uy, t*uy*uz + s*uz, t*uz**2 + c]

    end function rotation_matrix_along_arbitrary

    !
    !> 3D rotation matrix along an axis
    !
    !> @param[in] theta Counter-clockwise rotation angle in radian
    !> @param[in] axis character, = x, y or z to indicate the rotation axis
    !> @return r 3x3 rotation matrix
    !
    pure function rotation_matrix_along_axis(theta, axis) result(r)

        real, intent(in) :: theta
        character(len=*), intent(in) :: axis
        real, allocatable, dimension(:, :) :: r

        allocate (r(1:3, 1:3))
        select case (axis)

            case ('x')
                r(1, :) = [1.0, 0.0, 0.0]
                r(2, :) = [0.0, cos(theta), -sin(theta)]
                r(3, :) = [0.0, sin(theta), cos(theta)]

            case ('y')
                r(1, :) = [cos(theta), 0.0, sin(theta)]
                r(2, :) = [0.0, 1.0, 0.0]
                r(3, :) = [-sin(theta), 0.0, cos(theta)]

            case ('z')
                r(1, :) = [cos(theta), -sin(theta), 0.0]
                r(2, :) = [sin(theta), cos(theta), 0.0]
                r(3, :) = [0.0, 0.0, 1.0]

            case default
                r(1, :) = [1.0, 0.0, 0.0]
                r(2, :) = [0.0, 1.0, 0.0]
                r(3, :) = [0.0, 0.0, 1.0]

        end select

    end function rotation_matrix_along_axis

    function rotation_matrix_from_to(v1, v2) result(r)

        real, dimension(:) :: v1, v2
        real, allocatable, dimension(:, :) :: r

        integer :: n
        real, allocatable :: a(:, :), eigenval(:), eigenvecs(:, :)

        n = size(v1)
        a = matx(reshape(v1/norm2(v1), [n, 1]), reshape(v2/norm2(v2), [1, n]))
        call eigen(a + transpose(a), eigenval, eigenvecs)

        r = eye(n, 1.0) - 2*matx(eigenvecs(:, n:n), transpose(eigenvecs(:, n:n)))

    end function rotation_matrix_from_to

    !
    !> Rotation matrix in 2D, just for convenience
    !
    pure function rotation_matrix_2d(theta) result(r)

        real, intent(in) :: theta
        real, allocatable, dimension(:, :) :: r

        allocate (r(1:2, 1:2))
        r(1, :) = [cos(theta), -sin(theta)]
        r(2, :) = [sin(theta), cos(theta)]

    end function rotation_matrix_2d

    pure function rotation_matrix_3d_spherical(theta, phi) result(r)

        real, intent(in) :: theta, phi
        real, allocatable, dimension(:, :) :: r

        allocate (r(1:3, 1:3))
        !        r(1, :) = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
        !        r(2, :) = [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)]
        !        r(3, :) = [-sin(phi), cos(phi), 0.0]
        r = matmul(rotation_matrix_along_axis(phi, 'z'), rotation_matrix_along_axis(theta, 'y'))

    end function rotation_matrix_3d_spherical

    function rotation_matrix_3d_euler(angle, order) result(r)

        real, dimension(1:3) :: angle
        character(len=*), intent(in) :: order
        real, allocatable, dimension(:, :) :: r

        real, dimension(1:3, 1:3) :: rx, ry, rz

        rx = rotation_matrix_along_axis(angle(1), 'x')
        ry = rotation_matrix_along_axis(angle(2), 'y')
        rz = rotation_matrix_along_axis(angle(3), 'z')

        select case (order)
            case ('xyz')
                r = matx(rz, matx(ry, rx))
            case ('xzy')
                r = matx(ry, matx(rz, rx))
            case ('yxz')
                r = matx(rz, matx(rx, ry))
            case ('yzx')
                r = matx(rx, matx(rz, ry))
            case ('zxy')
                r = matx(ry, matx(rx, rz))
            case ('zyx')
                r = matx(rx, matx(ry, rz))
        end select

    end function rotation_matrix_3d_euler

    !
    !> Counter-clockwise rotation by angle in 2D
    !
    function rotate_point_2d(p, angle, origin) result(pr)

        real :: p(1:2), angle, origin(1:2)
        real, allocatable, dimension(:) :: pr

        pr = matx(rotation_matrix_along_axis(angle, 'z'), p - origin) + origin

    end function rotate_point_2d

    function rotate_points_2d(p, angle, origin) result(pr)

        real, dimension(:, :) :: p
        real :: angle, origin(1:2)
        real, allocatable, dimension(:, :) :: pr

        integer :: i, n
        real :: r(1:2, 1:2)

        r = rotation_matrix_along_axis(angle, 'z')

        n = size(p, 1)
        pr = p
        !$omp parallel do private(i)
        do i = 1, n
            pr(i, :) = matx(r, p(i, :) - origin) + origin
        end do
        !$omp end parallel do

    end function rotate_points_2d

    !
    !> Counter-clockwise rotation by angle(x, y, z) in 3D
    !
    function rotate_point_3d(p, angle, origin, order) result(pr)

        real, dimension(1:3) :: p, angle, origin
        character(len=*), intent(in) :: order
        real, allocatable, dimension(:) :: pr

        real, dimension(1:3, 1:3) :: rx, ry, rz

        rx = rotation_matrix_along_axis(angle(1), 'x')
        ry = rotation_matrix_along_axis(angle(2), 'y')
        rz = rotation_matrix_along_axis(angle(3), 'z')

        select case (order)
            case ('xyz')
                pr = matx(rz, matx(ry, matx(rx, p - origin))) + origin
            case ('xzy')
                pr = matx(ry, matx(rz, matx(rx, p - origin))) + origin
            case ('yxz')
                pr = matx(rz, matx(rx, matx(ry, p - origin))) + origin
            case ('yzx')
                pr = matx(rx, matx(rz, matx(ry, p - origin))) + origin
            case ('zxy')
                pr = matx(ry, matx(rx, matx(rz, p - origin))) + origin
            case ('zyx')
                pr = matx(rx, matx(ry, matx(rz, p - origin))) + origin
        end select

    end function rotate_point_3d

    function rotate_points_3d(p, angle, origin, order) result(pr)

        real, dimension(:, :) :: p
        real, dimension(1:3) :: angle, origin
        character(len=*), intent(in) :: order
        real, allocatable, dimension(:, :) :: pr

        real, dimension(1:3, 1:3) :: rx, ry, rz
        integer :: i, n

        rx = rotation_matrix_along_axis(angle(1), 'x')
        ry = rotation_matrix_along_axis(angle(2), 'y')
        rz = rotation_matrix_along_axis(angle(3), 'z')

        n = size(p, 1)
        pr = p

        select case (order)
            case ('xyz')
                !$omp parallel do private(i)
                do i = 1, n
                    pr(i, :) = matx(rz, matx(ry, matx(rx, p(i, :) - origin))) + origin
                end do
                !$omp end parallel do
            case ('xzy')
                !$omp parallel do private(i)
                do i = 1, n
                    pr(i, :) = matx(ry, matx(rz, matx(rx, p(i, :) - origin))) + origin
                end do
                !$omp end parallel do
            case ('yxz')
                !$omp parallel do private(i)
                do i = 1, n
                    pr(i, :) = matx(rz, matx(rx, matx(ry, p(i, :) - origin))) + origin
                end do
                !$omp end parallel do
            case ('yzx')
                !$omp parallel do private(i)
                do i = 1, n
                    pr(i, :) = matx(rx, matx(rz, matx(ry, p(i, :) - origin))) + origin
                end do
                !$omp end parallel do
            case ('zxy')
                !$omp parallel do private(i)
                do i = 1, n
                    pr(i, :) = matx(ry, matx(rx, matx(rz, p(i, :) - origin))) + origin
                end do
                !$omp end parallel do
            case ('zyx')
                !$omp parallel do private(i)
                do i = 1, n
                    pr(i, :) = matx(rx, matx(ry, matx(rz, p(i, :) - origin))) + origin
                end do
                !$omp end parallel do
        end select

    end function rotate_points_3d

    !
    !  Compute shortest distance between a line and a point
    !
    !  Author:
    !       John Burkardt, 2005
    !       K.G., 2024
    !
    function point_distance_to_line(p1, p2, p) result(dist)

        real, dimension(:) :: p1, p2, p
        real :: dist

        real :: bot, dot, t
        real, allocatable, dimension(:) :: pn

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

    end function point_distance_to_line

    !
    !  Compute signed shortest distance between a line and a point
    !
    !  Author:
    !       John Burkardt, 2005
    !       K.G., 2024
    !
    function point_signed_distance_to_line(p1, p2, p) result(dist_signed)

        real, dimension(1:2) :: p1, p2, p
        real :: dist_signed

        real :: a, b, c

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

    end function point_signed_distance_to_line

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
    subroutine point_distance_to_line_segement(p1, p2, p, pn, dist, t)

        real, dimension(:), intent(in) :: p, p1, p2
        real, dimension(:), intent(out) :: pn
        real, intent(out) :: t, dist

        !        integer, parameter :: ndim = 2
        real :: bot

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

    end subroutine point_distance_to_line_segement

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
    function polygon_contains_point(polygon, p) result(l)

        real, dimension(:, :) :: polygon
        real, dimension(:) :: p

        real, allocatable, dimension(:) :: x, y
        real :: x0, y0
        integer :: l, n, m, i, n0
        real :: angle, eps, pi, pi2, sum, theta, theta1, thetai, tol, u, v

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

    end function polygon_contains_point

    !
    !> Rotate the principal axis of an n-dimensional array p
    !> to target vector, where axis = 1, 2, ..., n
    !> indicating the largest to smallest eigenvalues
    !
    function rotate_principal_axis(p, axis, target, recenter) result(pr)

        real, dimension(:, :) :: p
        integer :: axis
        real, dimension(:) :: target
        logical, optional :: recenter

        real, allocatable :: mmt(:, :), eigenval(:), eigenvec(:, :), rr(:, :), pr(:, :)
        real, allocatable :: center(:)
        integer :: n, i

        n = size(p, 1)

        ! Shift data to (0, 0, ...)
        pr = p
        center = zeros(size(p, 2))
        !$omp parallel do private(i)
        do i = 1, size(p, 2)
            center(i) = mean(p(:, i))
            pr(:, i) = pr(:, i) - center(i)
        end do
        !$omp end parallel do

        ! Compute principal axes
        mmt = matx(transpose(pr), pr)
        call eigen(mmt, eigenval, eigenvec)

        ! Get rotation matrix for a principal axis to target
        rr = rotation_matrix(eigenvec(:, axis), target)

        ! Rotate
        !$omp parallel do private(i)
        do i = 1, n
            pr(i, :) = matx(rr, pr(i, :))
        end do
        !$omp end parallel do

        ! Shift back to original origin
        if (present(recenter)) then
            if (recenter) then
                !$omp parallel do private(i)
                do i = 1, size(p, 2)
                    pr(:, i) = pr(:, i) + center(i)
                end do
                !$omp end parallel do
            end if
        end if

    end function rotate_principal_axis

    !
    !> Find a best-fit, continuous, smooth curve for a 2D point cloud using polynomial or LOWESS
    !> The return are (1) the points projects onto the fitted curve, and
    !> (2) a 2D pixel image where the fit curve is 1 while the
    !> background is 0. To ensure sufficient smoothness, the sampling interval d
    !> should be properly small.
    !
    subroutine fit_curve(p1, p2, p1f, p2f, m, n, d, o, method, smooth, order, range)

        real, dimension(:), intent(in) :: p1, p2
        integer, dimension(1:2), intent(in), optional :: n
        real, dimension(1:2), intent(in), optional :: d, o
        character(len=*), intent(in), optional :: method
        real, intent(in), optional :: smooth
        integer, intent(in), optional :: order
        character(len=*), intent(in), optional :: range
        real, allocatable, dimension(:), intent(out) :: p1f, p2f
        real, allocatable, dimension(:, :), intent(out), optional :: m

        integer :: i, np
        type(polyfit1) :: ff
        real, allocatable :: p(:, :), pr(:, :), xx(:), yy(:)
        real, allocatable :: eigenval(:), eigenvec(:, :), mmt(:, :)
        real, allocatable :: r(:, :)
        real :: xmin, xmax
        real :: c1, c2
        integer :: n1, n2
        real :: d1, d2
        real :: o1, o2
        integer :: ind, jnd
        character(len=32) :: fit_method
        real :: fit_smooth
        integer :: fit_order
        real :: xextra

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
        r = matx(flip(eye(2, 1.0), axis=[2]), lsqsolve(eigenvec, eye(2, 1.0)))
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
                call ff%build(xx, yy, verbose=.false.)
                do i = 1, np
                    yy(i) = ff%eval(xx(i))
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
                        yy(i) = ff%eval(xx(i))
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

    end subroutine fit_curve

    !
    !> Find a best-fit, continuous, smooth curve for a 3D point cloud using polynomial or LOWESS
    !> The returns are (1) the points projected on the fitted surface, and
    !> (2) a 3D pixel image where the fit curve is 1 while the
    !> background is 0. To ensure sufficient smoothness, the sampling interval d
    !> should be properly small.
    !
    subroutine fit_surface(p1, p2, p3, p1f, p2f, p3f, m, n, d, o, method, smooth, order, range)

        real, dimension(:), intent(in) :: p1, p2, p3
        integer, dimension(1:3), intent(in), optional :: n
        real, dimension(1:3), intent(in), optional :: d, o
        character(len=*), intent(in), optional :: method
        real, intent(in), optional :: smooth
        integer, intent(in), optional :: order
        character(len=*), intent(in), optional :: range
        real, allocatable, dimension(:), intent(out) :: p1f, p2f, p3f
        real, allocatable, dimension(:, :, :), intent(out), optional :: m

        integer :: i, np, m3, m2
        type(polyfit2) :: ff
        real, allocatable :: p(:, :), pr(:, :), xx(:), yy(:), zz(:)
        real, allocatable :: eigenval(:), eigenvec(:, :), mmt(:, :)
        real, allocatable :: r(:, :)
        real :: xmin, xmax, ymin, ymax
        real :: c1, c2, c3
        integer :: n1, n2, n3
        real :: d1, d2, d3
        real :: o1, o2, o3
        integer :: ind, jnd, knd
        integer :: l
        integer, allocatable, dimension(:) :: qindex
        real, allocatable, dimension(:, :) :: pcx
        character(len=32) :: fit_method
        real :: fit_smooth
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
        r = matx(flip(eye(3, 1.0), axis=[2]), lsqsolve(eigenvec, eye(3, 1.0)))
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
                call ff%build(xx, yy, zz, verbose=.false.)
                do i = 1, np
                    zz(i) = ff%eval(xx(i), yy(i))
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
                        zz(i) = ff%eval(xx(i), yy(i))
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

    end subroutine fit_surface

end module libflit_geometry
