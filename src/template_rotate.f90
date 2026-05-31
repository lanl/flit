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

#define rotation_matrix_along_arbitrary_      CONCAT(rotation_matrix_along_arbitrary, T)
#define rotation_matrix_along_axis_      CONCAT(rotation_matrix_along_axis, T)
#define rotation_matrix_from_to_      CONCAT(rotation_matrix_from_to, T)
#define rotation_matrix_2d_      CONCAT(rotation_matrix_2d, T)
#define rotation_matrix_3d_spherical_      CONCAT(rotation_matrix_3d_spherical, T)
#define rotation_matrix_3d_euler_      CONCAT(rotation_matrix_3d_euler, T)
#define rotation_matrix_along_axis_      CONCAT(rotation_matrix_along_axis, T)
#define rotate_point_2d_      CONCAT(rotate_point_2d, T)
#define rotate_points_2d_      CONCAT(rotate_points_2d, T)
#define rotate_point_3d_      CONCAT(rotate_point_3d, T)
#define rotate_points_3d_      CONCAT(rotate_points_3d, T)
#define rotate_principal_axis_      CONCAT(rotate_principal_axis, T)

!
!> 3D rotation matrix along an arbitrary direction
!
pure function rotation_matrix_along_arbitrary_(theta, direction) result(r)

    !> Counter-clockwise rotation angle in radian
    TT, intent(in) :: theta
    !> 1D array of size 3 to represent the rotation axis
    TT, dimension(1:3), intent(in) :: direction
    !> 3x3 rotation matrix
    TT, allocatable, dimension(:, :) :: r

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

end function

!
!> 3D rotation matrix along an axis
!
!> @param[in] theta Counter-clockwise rotation angle in radian
!> @param[in] axis character, = x, y or z to indicate the rotation axis
!> @return r 3x3 rotation matrix
!
pure function rotation_matrix_along_axis_(theta, axis) result(r)

    TT, intent(in) :: theta
    character(len=*), intent(in) :: axis
    TT, allocatable, dimension(:, :) :: r

    allocate (r(1:3, 1:3))
    select case (axis)

        case ('x')
            r(1, :) = [1.0, 0.0, 0.0]
            r(2, :) = [real(0.0, fp), cos(theta), -sin(theta)]
            r(3, :) = [real(0.0, fp), sin(theta), cos(theta)]

        case ('y')
            r(1, :) = [cos(theta), real(0.0, fp), sin(theta)]
            r(2, :) = [0.0, 1.0, 0.0]
            r(3, :) = [-sin(theta), real(0.0, fp), cos(theta)]

        case ('z')
            r(1, :) = [cos(theta), -sin(theta), real(0.0, fp)]
            r(2, :) = [sin(theta), cos(theta), real(0.0, fp)]
            r(3, :) = [0.0, 0.0, 1.0]

        case default
            r(1, :) = [1.0, 0.0, 0.0]
            r(2, :) = [0.0, 1.0, 0.0]
            r(3, :) = [0.0, 0.0, 1.0]

    end select

end function

!
!> Rotation matrix corresponding to the rotation from v1 to v2
!
function rotation_matrix_from_to_(v1, v2) result(r)

    TT, dimension(:) :: v1, v2
    TT, allocatable, dimension(:, :) :: r

    integer :: n
    TT, allocatable :: a(:, :), eigenval(:), eigenvecs(:, :)

    n = size(v1)
    a = matx(reshape(v1/norm2(v1), [n, 1]), reshape(v2/norm2(v2), [1, n]))
    call eigen(a + transpose(a), eigenval, eigenvecs)

    r = eye(n, 1.0) - 2*matx(eigenvecs(:, n:n), transpose(eigenvecs(:, n:n)))

end function

!
!> Rotation matrix in 2D, just for convenience
!
pure function rotation_matrix_2d_(theta) result(r)

    TT, intent(in) :: theta
    TT, allocatable, dimension(:, :) :: r

    allocate (r(1:2, 1:2))
    r(1, :) = [cos(theta), -sin(theta)]
    r(2, :) = [sin(theta), cos(theta)]

end function

!
!> Rotation matrix in 3D with polar and azimuth angles only
!
pure function rotation_matrix_3d_spherical_(theta, phi) result(r)

    TT, intent(in) :: theta, phi
    TT, allocatable, dimension(:, :) :: r

    allocate (r(1:3, 1:3))
    !        r(1, :) = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
    !        r(2, :) = [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)]
    !        r(3, :) = [-sin(phi), cos(phi), 0.0]
    r = matmul(rotation_matrix_along_axis_(phi, 'z'), rotation_matrix_along_axis_(theta, 'y'))

end function

!
!> Rotation matrix in 3D with full angles and order
!
pure function rotation_matrix_3d_euler_(angle, order) result(r)

    TT, dimension(1:3), intent(in) :: angle
    character(len=*), intent(in) :: order
    TT, allocatable, dimension(:, :) :: r

    TT, dimension(1:3, 1:3) :: rx, ry, rz

    rx = rotation_matrix_along_axis_(angle(1), 'x')
    ry = rotation_matrix_along_axis_(angle(2), 'y')
    rz = rotation_matrix_along_axis_(angle(3), 'z')

    select case (order)
        case ('xyz')
            r = matmul(rz, matmul(ry, rx))
        case ('xzy')
            r = matmul(ry, matmul(rz, rx))
        case ('yxz')
            r = matmul(rz, matmul(rx, ry))
        case ('yzx')
            r = matmul(rx, matmul(rz, ry))
        case ('zxy')
            r = matmul(ry, matmul(rx, rz))
        case ('zyx')
            r = matmul(rx, matmul(ry, rz))
    end select

end function

!
!> Counter-clockwise rotation of a point by angle in 2D
!
function rotate_point_2d_(p, angle, origin) result(pr)

    TT, intent(in) :: p(1:2), angle, origin(1:2)
    TT, allocatable, dimension(:) :: pr

    pr = matx(rotation_matrix_along_axis_(angle, 'z'), p - origin) + origin

end function

!
!> Counter-clockwise rotation of points by angle in 2D
!
function rotate_points_2d_(p, angle, origin) result(pr)

    TT, dimension(:, :), intent(in) :: p
    TT, intent(in) :: angle, origin(1:2)
    TT, allocatable, dimension(:, :) :: pr

    integer :: i, n
    TT :: r(1:2, 1:2)

    r = rotation_matrix_along_axis_(angle, 'z')

    n = size(p, 1)
    pr = p
    !$omp parallel do private(i)
    do i = 1, n
        pr(i, :) = matx(r, p(i, :) - origin) + origin
    end do
    !$omp end parallel do

end function

!
!> Counter-clockwise rotation of a point by angle(x, y, z) in 3D
!
function rotate_point_3d_(p, angle, origin, order) result(pr)

    TT, dimension(1:3), intent(in) :: p, angle, origin
    character(len=*), intent(in) :: order
    TT, allocatable, dimension(:) :: pr

    TT, dimension(1:3, 1:3) :: rx, ry, rz

    rx = rotation_matrix_along_axis_(angle(1), 'x')
    ry = rotation_matrix_along_axis_(angle(2), 'y')
    rz = rotation_matrix_along_axis_(angle(3), 'z')

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

end function

!
!> Counter-clockwise rotation of points by angle(x, y, z) in 3D
!
function rotate_points_3d_(p, angle, origin, order) result(pr)

    TT, dimension(:, :), intent(in) :: p
    TT, dimension(1:3), intent(in) :: angle, origin
    character(len=*), intent(in) :: order
    TT, allocatable, dimension(:, :) :: pr

    TT, dimension(1:3, 1:3) :: rx, ry, rz
    integer :: i, n

    rx = rotation_matrix_along_axis_(angle(1), 'x')
    ry = rotation_matrix_along_axis_(angle(2), 'y')
    rz = rotation_matrix_along_axis_(angle(3), 'z')

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

end function

!
!> Rotate the principal axis of an n-dimensional array p
!> to target vector, where axis = 1, 2, ..., n
!> indicating the largest to smallest eigenvalues
!
function rotate_principal_axis_(p, axis, target, recenter) result(pr)

    TT, dimension(:, :), intent(in) :: p
    integer, intent(in) :: axis
    TT, dimension(:), intent(in) :: target
    logical, intent(in), optional :: recenter

    TT, allocatable :: mmt(:, :), eigenval(:), eigenvec(:, :), rr(:, :), pr(:, :)
    TT, allocatable :: center(:)
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
    rr = rotation_matrix_from_to_(eigenvec(:, axis), target)

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

end function

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef rotation_matrix_along_arbitrary_
#undef rotation_matrix_along_axis_
#undef rotation_matrix_from_to_
#undef rotation_matrix_2d_
#undef rotation_matrix_3d_spherical_
#undef rotation_matrix_3d_euler_
#undef rotation_matrix_along_axis_
#undef rotate_point_2d_
#undef rotate_points_2d_
#undef rotate_point_3d_
#undef rotate_points_3d_
#undef rotate_principal_axis_


