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

#define trace_     CONCAT(trace, T)
#define circulant_     CONCAT(circulant, T)
#define toeplitz_     CONCAT(toeplitz, T)
#define hankel_     CONCAT(hankel, T)
#define vandermonde_     CONCAT(vandermonde, T)
#define svd_     CONCAT(svd, T)
#define det_nxn_     CONCAT(det_nxn, T)
#define inv_nxn_     CONCAT(inv_nxn, T)
#define solve_single_rhs_     CONCAT(solve_single_rhs, T)
#define solve_multiple_rhs_     CONCAT(solve_multiple_rhs, T)
#define least_squares_solve_single_rhs_     CONCAT(least_squares_solve_single_rhs, T)
#define least_squares_solve_multiple_rhs_     CONCAT(least_squares_solve_multiple_rhs, T)
#define solve_band_single_rhs_     CONCAT(solve_band_single_rhs, T)
#define matx_mat_vec_     CONCAT(matx_mat_vec, T)
#define matx_mat_mat_     CONCAT(matx_mat_mat, T)
#define diagx_mat_mat_     CONCAT(diagx_mat_mat, T)
#define xdiag_mat_mat_     CONCAT(xdiag_mat_mat, T)
!!!#define matx_mat_vec_   CONCAT(CONCAT(matx_mat_vec_, T), T)

function trace_(m) result(tr)

    TT, dimension(:, :), intent(in) :: m
    TT :: tr

    integer :: i

    tr = 0.0
    do i = 1, min(size(m, 1), size(m, 2))
        tr = tr + m(i, i)
    end do

end function trace_

function circulant_(c, nrow) result(t)

    TT, dimension(:) :: c
    integer, intent(in), optional :: nrow
    TT, allocatable, dimension(:, :) :: t

    integer :: n, i
    integer :: n2

    n = size(c)
    if (present(nrow)) then
        n2 = nrow
    else
        n2 = n
    end if

    call alloc_array(t, [1, n, 1, n2])

    t(:, 1) = c(:)
    do i = 2, n2
        t(:, i) = cshift(t(:, i - 1), -1)
    end do

end function circulant_

function toeplitz_(c, r) result(t)

    TT, dimension(:) :: c, r
    TT, allocatable, dimension(:, :) :: t

    integer :: n1, n2, i, j

    call assert(c(1) == r(1), ' Error: c(1) /= r(1) ')

    n1 = size(c)
    n2 = size(r)

    call alloc_array(t, [1, n1, 1, n2])

    do i = 1, n1
        t(i, i:n2) = r(1:n2 - i + 1)
    end do
    do j = 1, n2
        t(j:n1, j) = c(1:n1 - j + 1)
    end do

end function toeplitz_

function hankel_(x) result(h)

    TT, dimension(:), intent(in) :: x
    TT, allocatable, dimension(:, :) :: h

    integer :: nx, m, n, i

    nx = size(x)
    if (mod(nx, 2) == 0) then
        n = nx/2
        m = nx/2 + 1
    else
        n = (nx + 1)/2
        m = n
    end if

    call alloc_array(h, [1, n, 1, n])

    do i = 1, n
        h(:, i) = x(i:i + m - 1)
    end do

end function hankel_

function vandermonde_(x, n) result(h)

    TT, dimension(:) :: x
    integer :: n
    TT, allocatable, dimension(:, :) :: h

    integer :: i

    call alloc_array(h, [1, size(x), 1, n])

    do i = 1, n
        h(:, i) = x**(i - 1.0)
    end do

end function vandermonde_


function det_nxn_(w) result(wdet)

    TT, dimension(:, :), intent(in) :: w

    real :: wdet, sgn
    integer :: n1, n2, i
    integer, allocatable, dimension(:) :: ipiv
    TT, allocatable, dimension(:, :) :: wt

    n1 = size(w, 1)
    n2 = size(w, 2)

    if (n1 == 2 .and. n2 == 2) then
        ! 2x2 matrix

        wdet = w(1, 1)*w(2, 2) - w(1, 2)*w(2, 1)

    else if (n1 == 3 .and. n2 == 3) then
        ! 3x3 matrix

        wdet = &
            -w(1, 3)*w(2, 2)*w(3, 1) &
            + w(1, 2)*w(2, 3)*w(3, 1) &
            + w(1, 3)*w(2, 1)*w(3, 2) &
            - w(1, 1)*w(2, 3)*w(3, 2) &
            - w(1, 2)*w(2, 1)*w(3, 3) &
            + w(1, 1)*w(2, 2)*w(3, 3)

    else
        ! nxn matrix

        allocate (wt(1:n1, 1:n2), source=w)

        allocate (ipiv(1:n1))
        ipiv = 0
        call getrf(wt, ipiv)

        wdet = 1.0
        do i = 1, n1
            wdet = wdet*wt(i, i)
        end do

        sgn = 1.0
        do i = 1, n1
            if (ipiv(i) /= i) then
                sgn = -sgn
            end if
        end do

        wdet = sgn*wdet

    end if

end function det_nxn_

function inv_nxn_(w) result(winv)

    TT, dimension(:, :), intent(in) :: w

    integer :: n1, n2
    TT, allocatable, dimension(:, :) :: winv
    integer, allocatable, dimension(:) :: ipiv
    TT :: wdet

    n1 = size(w, 1)
    n2 = size(w, 2)

    if (n1 == 2 .and. n2 == 2) then
        ! 2x2 matrix

        allocate (winv(1:n1, 1:n2))

        wdet = det_nxn_(w)
        winv(1, 1) = w(2, 2)
        winv(1, 2) = -w(1, 2)
        winv(2, 1) = -w(2, 1)
        winv(2, 2) = w(1, 1)
        winv = winv/wdet

    else if (n1 == 3 .and. n2 == 3) then
        ! 3x3 matrix

        allocate (winv(1:n1, 1:n2))

        wdet = det_nxn_(w)
        winv(1, 1) = (-w(2, 3)*w(3, 2) + w(2, 2)*w(3, 3))/wdet
        winv(1, 2) = (w(1, 3)*w(3, 2) - w(1, 2)*w(3, 3))/wdet
        winv(1, 3) = (-w(1, 3)*w(2, 2) + w(1, 2)*w(2, 3))/wdet
        winv(2, 1) = (w(2, 3)*w(3, 1) - w(2, 1)*w(3, 3))/wdet
        winv(2, 2) = (-w(1, 3)*w(3, 1) + w(1, 1)*w(3, 3))/wdet
        winv(2, 3) = (w(1, 3)*w(2, 1) - w(1, 1)*w(2, 3))/wdet
        winv(3, 1) = (-w(2, 2)*w(3, 1) + w(2, 1)*w(3, 2))/wdet
        winv(3, 2) = (w(1, 2)*w(3, 1) - w(1, 1)*w(3, 2))/wdet
        winv(3, 3) = (-w(1, 2)*w(2, 1) + w(1, 1)*w(2, 2))/wdet

    else
        ! nxn matrix

        allocate (winv(1:n1, 1:n2), source=w)
        allocate (ipiv(1:n1))
        ipiv = 0

        call getrf(winv, ipiv)
        call getri(winv, ipiv)

    end if

end function inv_nxn_

subroutine svd_(a, s, u, vt)

    TT, dimension(:, :), intent(in) :: a
    TTT, allocatable, dimension(:), intent(out) :: s
    TT, allocatable, dimension(:, :), intent(out) :: u, vt

    integer :: m, n
    TT, allocatable, dimension(:, :) :: w

    m = size(a, 1)
    n = size(a, 2)

    call alloc_array(w, [1, m, 1, n], source=a)
    call alloc_array(s, [1, min(m, n)])
    call alloc_array(u, [1, m, 1, m])
    call alloc_array(vt, [1, n, 1, n])

    call gesvd(w, s, u, vt)

end subroutine svd_

function solve_single_rhs_(a, b) result(x)

    TT, dimension(:, :) :: a
    TT, dimension(:) :: b
    TT, allocatable, dimension(:) :: x

    TT, allocatable, dimension(:, :) :: ac

    ac = a
    x = b
    call gesv(ac, x)

end function solve_single_rhs_

function solve_multiple_rhs_(a, b) result(x)

    TT, dimension(:, :) :: a
    TT, dimension(:, :) :: b
    TT, allocatable, dimension(:, :) :: x

    TT, allocatable, dimension(:, :) :: ac

    ac = a
    x = b

    call gesv(ac, x)

end function solve_multiple_rhs_

function least_squares_solve_single_rhs_(a, b) result(x)

    TT, dimension(:, :) :: a
    TT, dimension(:) :: b
    TT, allocatable, dimension(:) :: x, xt

    integer :: m, n
    TT, allocatable, dimension(:, :) :: ac
    m = size(a, 1)
    n = size(a, 2)

    ac = a

    if (size(b) /= m) then
        print *, ' Error: size(a) and size(b) mismatch. '
        stop
    end if

    xt = zeros(max(m, n))
    xt(1:m) = b
    call gels(ac, xt)

    x = xt(1:n)

end function least_squares_solve_single_rhs_

function least_squares_solve_multiple_rhs_(a, b) result(x)

    TT, dimension(:, :) :: a
    TT, dimension(:, :) :: b
    TT, allocatable, dimension(:, :) :: x, xt

    integer :: m, n, nrhs
    TT, allocatable, dimension(:, :) :: ac
    m = size(a, 1)
    n = size(a, 2)
    nrhs = size(b, 2)
    !        allocate (ac(1:m, 1:n), source=a)
    ac = a

    if (size(b, 1) /= m) then
        print *, ' Error: size(a) and size(b) mismatch. '
        stop
    end if

    xt = zeros(max(m, n), nrhs)
    xt(1:m, 1:nrhs) = b
    call gels(ac, xt)

    x = xt(1:n, 1:nrhs)

end function least_squares_solve_multiple_rhs_

function solve_band_single_rhs_(a, nu, nl, b) result(x)

    TT, dimension(:, :), intent(in) :: a
    integer, intent(in) :: nu, nl
    TT, dimension(:), intent(in) :: b
    TT, allocatable, dimension(:) :: x

    TT, allocatable, dimension(:, :) :: ac
    integer :: m, n, i

    ! band storage: https://www.netlib.org/lapack/lug/node124.html
    ! https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-fortran/top/lapack-routines/lapack-linear-equation-routines/lapack-linear-equation-driver-routines/gbsv.html
    ! http://www.netlib.org/lapack/explore-html/dc/db2/group__real_g_bsolve_ga3656935309a19ed624052103572a4a47.html#ga3656935309a19ed624052103572a4a47

    m = size(a, 1)
    call assert(m == nu + nl + 1, &
        '<solve_band_single_rhs> Error: Dimension of a is inconsistent.')
    n = size(a, 2)
    do i = 1, nu
        call assert(all(a(nu + 1 - i, 1:i) == 0), &
            '<solve_band_single_rhs> Error: Upper half is not consistent.')
    end do
    do i = 1, nl
        call assert(all(a(nu + 1 + i, n - i + 1:n) == 0), &
            '<solve_band_single_rhs> Error: Lower half is not consistent.')
    end do

    ! dim1 of input ac = m + nl, where the first nl rows are zeros, and nl + 1: rows are a
    ac = zeros(m + nl, n)
    ac(nl + 1:m + nl, :) = a
    x = b

    call gbsv(ac, x, nl)

end function solve_band_single_rhs_

function matx_mat_vec_(a, x) result(y)

    TT, dimension(:, :) :: a
    TT, dimension(:) :: x

    TT, allocatable, dimension(:) :: y
    integer :: n

    n = size(a, 1)
    call alloc_array(y, [1, n])

    call gemv(a, x, y)

end function matx_mat_vec_

function matx_mat_mat_(a, x) result(y)

    TT, dimension(:, :) :: a
    TT, dimension(:, :) :: x

    TT, allocatable, dimension(:, :) :: y
    integer :: n1, n2

    n1 = size(a, 1)
    n2 = size(x, 2)
    call alloc_array(y, [1, n1, 1, n2])

    call gemm(a, x, y)

end function matx_mat_mat_

function diagx_mat_mat_(a, x) result(y)

    TT, dimension(:) :: a
    TT, dimension(:, :) :: x

    TT, allocatable, dimension(:, :) :: y
    integer :: n1, n2, i

    n1 = size(a)
    n2 = size(x, 2)
    call alloc_array(y, [1, n1, 1, n2])

    do i = 1, n1
        y(i, :) = x(i, :)*a(i)
    end do

end function diagx_mat_mat_

function xdiag_mat_mat_(x, a) result(y)

    TT, dimension(:, :) :: x
    TT, dimension(:) :: a

    TT, allocatable, dimension(:, :) :: y
    integer :: n1, n2, i

    n1 = size(x, 1)
    n2 = size(a)
    call alloc_array(y, [1, n1, 1, n2])

    do i = 1, n2
        y(:, i) = x(:, i)*a(i)
    end do

end function xdiag_mat_mat_

#undef T
#undef TT
#undef TTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef trace_
#undef circulant_
#undef toeplitz_
#undef hankel_
#undef vandermonde_
#undef svd_
#undef det_nxn_
#undef inv_nxn_
#undef solve_single_rhs_
#undef solve_multiple_rhs_
#undef least_squares_solve_single_rhs_
#undef least_squares_solve_multiple_rhs_
#undef solve_band_single_rhs_
#undef matx_mat_vec_
#undef matx_mat_mat_
#undef diagx_mat_mat_
#undef xdiag_mat_mat_
