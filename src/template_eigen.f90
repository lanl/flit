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

#define spectral_radius_     CONCAT(spectral_radius, T)
#define eigen_symm2x2_     CONCAT(eigen_symm2x2, T)
#define eigen_symm3x3_     CONCAT(eigen_symm3x3, T)
#define eigen_simple_symmetric_l_or_r_     CONCAT(eigen_simple_symmetric_l_or_r, T)

subroutine spectral_radius_(a, v, b, niter)

    TT, dimension(:, :), intent(in) :: a
    TT, intent(inout) :: v
    TT, allocatable, dimension(:), intent(inout) :: b
    integer, intent(in), optional :: niter

    TT, allocatable, dimension(:) :: b_k, b_k1
    integer :: n, i, spectral_radius_niter

    n = size(a, 2)
    b_k = random(n)
    b_k1 = zeros(n)

    if (present(niter)) then
        spectral_radius_niter = niter
    else
        spectral_radius_niter = 10
    end if

    i = 1
    do while (i <= spectral_radius_niter)

        b_k1 = matx(a, b_k)
        b_k = b_k1/norm2(b_k1)
        i = i + 1

    end do

    v = sum(b_k1*b_k)/sum(b_k**2)
    b = b_k

end subroutine spectral_radius_

!
!> Eigen solver for 2x2 symmetric matrix
!
subroutine eigen_symm2x2_(a, d, v, order)

    TT, dimension(:, :), intent(in) :: a
    TT, allocatable, dimension(:), intent(inout) :: d
    TT, allocatable, dimension(:, :), intent(inout) :: v
    integer, intent(in), optional :: order

    TT, parameter :: flt_epsilon = 1.0e-6
    TT :: a00, a01, a11
    TT :: v00, v01, v10, v11
    TT :: c, r, s, t, u, vpr, vqr, edt
    TT, allocatable :: vt(:)
    TT :: vtiny
    integer :: eigenval_sort

    if (present(order)) then
        eigenval_sort = order
    else
        eigenval_sort = -1
    end if

    ! Initialize
    d = zeros(2)
    v = zeros(2, 2)
    vt = zeros(2)

    ! Copy matrix to local variables.
    a00 = a(1, 1)
    a01 = a(1, 2)
    a11 = a(2, 2)

    ! Initial eigenvectors.
    v00 = 1.0d0
    v01 = 0.0d0
    v10 = 0.0d0
    v11 = 1.0d0

    ! If off-diagonal element is non-zero, zero it with a Jacobi rotation.
    if (a01 /= 0.0d0) then

        ! avoid overflow in r*r below
        vtiny = 0.1d0*sqrt(flt_epsilon)

        u = a11 - a00
        if (abs(a01) < vtiny*abs(u)) then
            t = a01/u
        else
            r = 0.5d0*u/a01
            if (r >= 0.0d0) then
                t = 1.0d0/(r + sqrt(1.0d0 + r**2))
            else
                t = 1.0d0/(r - sqrt(1.0d0 + r**2))
            end if
        end if

        c = 1.0d0/sqrt(1.0d0 + t**2)
        s = t*c
        u = s/(1.0d0 + c)
        r = t*a01
        a00 = a00 - r
        a11 = a11 + r
        vpr = v00
        vqr = v10
        v00 = vpr - s*(vqr + vpr*u)
        v10 = vqr + s*(vpr - vqr*u)
        vpr = v01
        vqr = v11
        v01 = vpr - s*(vqr + vpr*u)
        v11 = vqr + s*(vpr - vqr*u)

    end if

    ! Copy eigenvalues and eigenvectors to output arrays.
    d(1) = a00
    d(2) = a11
    v(1, 1) = v00
    v(1, 2) = v01
    v(2, 1) = v10
    v(2, 2) = v11

    select case (eigenval_sort)
        case (1)
            ! Sort eigenvalues (and eigenvectors) in descending order.
            if (d(1) > d(2)) then

                edt = d(2)
                d(2) = d(1)
                d(1) = edt
                vt(1) = v(2, 1)
                vt(2) = v(2, 2)
                v(2, 1) = v(1, 1)
                v(2, 2) = v(1, 2)
                v(1, 1) = vt(1)
                v(1, 2) = vt(2)

            end if
        case (-1)
            ! Sort eigenvalues (and eigenvectors) in descending order.
            if (d(1) < d(2)) then

                edt = d(2)
                d(2) = d(1)
                d(1) = edt
                vt(1) = v(2, 1)
                vt(2) = v(2, 2)
                v(2, 1) = v(1, 1)
                v(2, 2) = v(1, 2)
                v(1, 1) = vt(1)
                v(1, 2) = vt(2)

            end if
    end select

    ! Transpose to make each colum an eigenvector
    v = transpose(v)

end subroutine eigen_symm2x2_

!
!> Eigen solver for 3x3 symmetric matrix
!
!> Efficient numerical diagonalization of Hermitian 3x3 matrices, Int. J. Mod. Phys. C 19 (2008) 523-548
!> https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/
!
subroutine eigen_symm3x3_(ar, w, q, order)

    TT, dimension(:, :), intent(in) :: ar
    TT, allocatable, dimension(:), intent(inout) :: w
    TT, allocatable, dimension(:, :), intent(inout) :: q
    integer, optional :: order

    TT :: sd, so, s, c, t, g, h, z, theta, thresh
    integer :: i, x, y, r, j
    TT, allocatable :: wd(:), a(:, :), qd(:, :)
    integer :: n, eigenval_sort

    if (present(order)) then
        eigenval_sort = order
    else
        eigenval_sort = -1
    end if

    ! Initialize
    w = zeros(3)
    q = zeros(3, 3)
    wd = zeros(3)
    qd = zeros(3, 3)

    n = 3
    a = ar

    ! initialize q to the identitity matrix
    ! --- this loop can be omitted if only the eigenvalues are desired ---
    q = 0.0d0
    do x = 1, n
        q(x, x) = 1.0d0
    end do

    ! initialize w to diag(a)
    do x = 1, n
        w(x) = a(x, x)
    end do

    ! calculate sqr(tr(a))
    sd = sum(abs(w))**2

    ! main iteration loop
    do i = 1, 50

        ! test for convergence
        so = 0.0d0
        do x = 1, n
            do y = x + 1, n
                so = so + abs(a(x, y))
            end do
        end do
        if (so == 0) then
            exit
        end if

        if (i < 4) then
            thresh = 0.2d0*so/n**2
        else
            thresh = 0.0d0
        end if

        ! do sweep
        do x = 1, n
            do y = x + 1, n
                g = 100.0d0*(abs(a(x, y)))
                if ((i > 4) .and. (abs(w(x)) + g == abs(w(x))) .and. (abs(w(y)) + g == abs(w(y)))) then
                    a(x, y) = 0.0d0
                else if (abs(a(x, y)) > thresh) then
                    ! calculate jacobi transformation
                    h = w(y) - w(x)
                    if (abs(h) + g == abs(h)) then
                        t = a(x, y)/h
                    else
                        theta = 0.5d0*h/a(x, y)
                        if (theta < 0.0d0) then
                            t = -1.0d0/(sqrt(1.0d0 + theta**2) - theta)
                        else
                            t = 1.0d0/(sqrt(1.0d0 + theta**2) + theta)
                        end if
                    end if

                    c = 1.0d0/sqrt(1.0d0 + t**2)
                    s = t*c
                    z = t*a(x, y)

                    ! apply jacobi transformation
                    a(x, y) = 0.0d0
                    w(x) = w(x) - z
                    w(y) = w(y) + z
                    do r = 1, x - 1
                        t = a(r, x)
                        a(r, x) = c*t - s*a(r, y)
                        a(r, y) = s*t + c*a(r, y)
                    end do
                    do r = x + 1, y - 1
                        t = a(x, r)
                        a(x, r) = c*t - s*a(r, y)
                        a(r, y) = s*t + c*a(r, y)
                    end do
                    do r = y + 1, n
                        t = a(x, r)
                        a(x, r) = c*t - s*a(y, r)
                        a(y, r) = s*t + c*a(y, r)
                    end do
                    ! update eigenvectors
                    ! --- this loop can be omitted if only the eigenvalues are desired ---
                    do r = 1, n
                        t = q(r, x)
                        q(r, x) = c*t - s*q(r, y)
                        q(r, y) = s*t + c*q(r, y)
                    end do
                end if
            end do
        end do

    end do

    ! Sort eigenvalues (and eigenvectors) to descending order
    wd = sort(w, eigenval_sort)
    do i = 1, 3
        do j = 1, 3
            if (wd(i) == w(j)) then
                qd(:, i) = q(:, j)
                exit
            end if
        end do
    end do

    ! Copy to output
    w = wd
    q = qd

end subroutine eigen_symm3x3_

!
! A X = lambada X with A being real, symmetric
! Returns eigenvalues and left or right eigenvectors
!
subroutine eigen_simple_symmetric_l_or_r_(a, w, v, lr, order)

    TT, dimension(:, :), intent(in) :: a
    TT, allocatable, dimension(:), intent(inout) :: w
    TT, allocatable, dimension(:, :), intent(inout), optional :: v
    character(len=*), intent(in), optional :: lr
    integer, intent(in), optional :: order

    character(len=12) :: left_or_right
    integer :: n
    TT, allocatable :: aa(:, :), wr(:), wi(:), vl(:, :), vr(:, :)
    integer, allocatable, dimension(:) :: windex
    integer :: eigenval_sort

    call assert(size(a, 1) == size(a, 2), 'Error: A must be square. ')
    call assert(all(a == transpose(a)), 'Error: A must be symmetric. ')

    if (present(lr)) then
        left_or_right = lr
    else
        left_or_right = 'right'
    end if

    if (present(order)) then
        eigenval_sort = order
    else
        eigenval_sort = -1
    end if

    n = size(a, 1)
    aa = a
    wr = zeros(n)
    wi = zeros(n)
    vl = zeros(n, n)
    vr = zeros(n, n)
    v = zeros(n, n)

    call geev(aa, wr, wi, vl, vr)

    w = wr

    ! Sort the eigenvalues descendingly
    allocate (windex(1:n))
    call sort_index(w, windex, eigenval_sort)

    if (present(v)) then
        select case (left_or_right)
            case ('right')
                v = vr
            case ('left')
                v = vl
        end select

        ! Arrange the eigenvectors accordingly
        v = v(:, windex)
    end if

end subroutine eigen_simple_symmetric_l_or_r_

#undef T
#undef TT
#undef TTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef spectral_radius_
#undef eigen_symm2x2_
#undef eigen_symm3x3_
#undef eigen_simple_symmetric_l_or_r_
