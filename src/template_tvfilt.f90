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

#define tv_iso_filt_1d_      CONCAT(tv_iso_filt_1d, T)
#define tv_iso_filt_2d_      CONCAT(tv_iso_filt_2d, T)
#define tv_iso_filt_3d_      CONCAT(tv_iso_filt_3d, T)
#define sparse_tv_filt_1d_      CONCAT(sparse_tv_filt_1d, T)
#define sparse_tv_filt_2d_      CONCAT(sparse_tv_filt_2d, T)
#define sparse_tv_filt_3d_      CONCAT(sparse_tv_filt_3d, T)
#define tgpv_filt_1d_      CONCAT(tgpv_filt_1d, T)
#define tgpv_filt_2d_      CONCAT(tgpv_filt_2d, T)
#define tgpv_filt_3d_      CONCAT(tgpv_filt_3d, T)
#define tgpv_filt_2d_mpi_      CONCAT(tgpv_filt_2d_mpi, T)
#define tgpv_filt_3d_mpi_      CONCAT(tgpv_filt_3d_mpi, T)
#define soft_shrinkage_1d_      CONCAT(soft_shrinkage_1d, T)
#define soft_shrinkage_2d_      CONCAT(soft_shrinkage_2d, T)
#define soft_shrinkage_3d_      CONCAT(soft_shrinkage_3d, T)

!
!> 1D isotropic TV denoising
!
function tv_iso_filt_1d_(f, mu, niter, verbose) result(u)

    TT, dimension(:), intent(in) :: f
    TT, intent(in) :: mu
    integer, intent(in) :: niter
    logical, optional :: verbose
    TT, allocatable, dimension(:) :: u

    integer :: i, iter
    TT :: tmpx, sumu
    TT, allocatable, dimension(:) :: d1, d1t
    TT, allocatable, dimension(:) :: pu
    integer :: n1
    integer :: i1, i2
    integer :: ii1, ii2
    TT :: lambda, sk
    logical :: tv_verbose

    n1 = size(f)

    if (present(verbose)) then
        tv_verbose = verbose
    else
        tv_verbose = .true.
    end if

    ! memory
    d1 = zeros(n1)
    d1t = zeros(n1)

    ! read noisy image
    u = f

    ! lambdas for u and w subproblems
    lambda = 2.0*mu

    ! TGpV iteration
    do iter = 1, niter

        pu = u

        ! minimization u
        do i = 1, n1

            if (i == 1) then
                i1 = 0
                ii1 = i
            else
                i1 = 1
                ii1 = i - 1
            end if

            if (i == n1) then
                i2 = 0
                ii2 = i
            else
                i2 = 1
                ii2 = i + 1
            end if

            sumu = lambda*( &
                i2*u(ii2) + i1*u(ii1)) &
                + lambda*( &
                (i1*d1(ii1) - i2*d1(i)) &
                - (i1*d1t(ii1) - i2*d1t(i)) &
                ) &
                + mu*f(i)

            u(i) = sumu/(mu + (2.0 - (1 - i1) - (1 - i2))*lambda)

        end do

        ! update multiplietgpv_filt_2d_rs
        ! minimization d
        !$omp parallel do private(i, i1, i2, tmpx, sk)
        do i = 1, n1

            if (i == n1) then
                i1 = i - 1
                i2 = i
            else
                i1 = i
                i2 = i + 1
            end if
            tmpx = (u(i2) - u(i1)) + d1t(i)

            sk = abs(tmpx)

            d1(i) = max(sk - 1.0/lambda, 0.0)*tmpx/(sk + float_tiny)

        end do
        !$omp end parallel do

        ! d
        !$omp parallel do private(i, i1, i2)
        do i = 1, n1

            if (i == n1) then
                i1 = i - 1
                i2 = i
            else
                i1 = i
                i2 = i + 1
            end if
            d1t(i) = d1t(i) + (u(i2) - u(i1) - d1(i))

        end do
        !$omp end parallel do

        ! progress
        if (tv_verbose .and. (mod(iter, max(nint(niter/10.0), 1)) == 0 .or. iter == 1)) then
            call warn(' >> Isotropic TV iteration '//tidy(num2str(iter, '(i)')) &
                //' of '//tidy(num2str(niter, '(i)')//' relative l2-norm difference = ' &
                //tidy(num2str(norm2(u - pu), '(es)'))))
        end if

    end do

end function tv_iso_filt_1d_

!
!> 2D isotropic TV denoising
!
function tv_iso_filt_2d_(f, mu, niter, verbose) result(u)

    TT, dimension(:, :), intent(in) :: f
    TT, intent(in) :: mu
    integer, intent(in) :: niter
    logical, optional :: verbose
    TT, allocatable, dimension(:, :) :: u

    integer :: i, j, iter
    TT :: tmpx, tmpz, sumu
    TT, allocatable, dimension(:, :) :: d1, d2, d1t, d2t
    TT, allocatable, dimension(:, :) :: pu
    integer :: n1, n2
    integer :: i1, i2, j1, j2
    integer :: ii1, ii2, jj1, jj2
    TT :: lambda, sk
    logical :: tv_verbose

    n1 = size(f, 1)
    n2 = size(f, 2)

    if (present(verbose)) then
        tv_verbose = verbose
    else
        tv_verbose = .true.
    end if

    ! memory
    d1 = zeros(n1, n2)
    d2 = zeros(n1, n2)
    d1t = zeros(n1, n2)
    d2t = zeros(n1, n2)

    ! read noisy image
    u = f
    pu = zeros_like(u)

    ! lambdas for u and w subproblems
    lambda = 2.0*mu

    ! TGpV iteration
    do iter = 1, niter

        pu = u

        ! minimization u
        do j = 1, n2
            do i = 1, n1

                if (i == 1) then
                    i1 = 0
                    ii1 = i
                else
                    i1 = 1
                    ii1 = i - 1
                end if

                if (i == n1) then
                    i2 = 0
                    ii2 = i
                else
                    i2 = 1
                    ii2 = i + 1
                end if

                if (j == 1) then
                    j1 = 0
                    jj1 = j
                else
                    j1 = 1
                    jj1 = j - 1
                end if
                if (j == n2) then
                    j2 = 0
                    jj2 = j
                else
                    j2 = 1
                    jj2 = j + 1
                end if

                sumu = lambda*( &
                    i2*u(ii2, j) + i1*u(ii1, j) &
                    + j2*u(i, jj2) + j1*u(i, jj1)) &
                    + lambda*( &
                    (i1*d1(ii1, j) - i2*d1(i, j)) &
                    + (j1*d2(i, jj1) - j2*d2(i, j)) &
                    - (i1*d1t(ii1, j) - i2*d1t(i, j)) &
                    - (j1*d2t(i, jj1) - j2*d2t(i, j)) &
                    ) &
                    + mu*f(i, j)

                u(i, j) = sumu/(mu + (4.0 - (1 - i1) - (1 - i2) - (1 - j1) - (1 - j2))*lambda)

            end do
        end do

        ! update multiplietgpv_filt_2d_rs
        ! minimization d
        !$omp parallel do private(i, j, i1, i2, j1, j2, tmpx, tmpz, sk)
        do j = 1, n2
            do i = 1, n1

                if (i == n1) then
                    i1 = i - 1
                    i2 = i
                else
                    i1 = i
                    i2 = i + 1
                end if
                tmpx = (u(i2, j) - u(i1, j)) + d1t(i, j)

                if (j == n2) then
                    j1 = j - 1
                    j2 = j
                else
                    j1 = j
                    j2 = j + 1
                end if
                tmpz = (u(i, j2) - u(i, j1)) + d2t(i, j)

                sk = sqrt(tmpx**2 + tmpz**2)

                d1(i, j) = max(sk - 1.0/lambda, 0.0)*tmpx/(sk + float_tiny)
                d2(i, j) = max(sk - 1.0/lambda, 0.0)*tmpz/(sk + float_tiny)

            end do
        end do
        !$omp end parallel do

        ! d
        !$omp parallel do private(i, j, i1, i2, j1, j2)
        do j = 1, n2
            do i = 1, n1

                if (i == n1) then
                    i1 = i - 1
                    i2 = i
                else
                    i1 = i
                    i2 = i + 1
                end if
                d1t(i, j) = d1t(i, j) + (u(i2, j) - u(i1, j) - d1(i, j))

                if (j == n2) then
                    j1 = j - 1
                    j2 = j
                else
                    j1 = j
                    j2 = j + 1
                end if
                d2t(i, j) = d2t(i, j) + (u(i, j2) - u(i, j1) - d2(i, j))

            end do
        end do
        !$omp end parallel do

        ! progress
        if (tv_verbose .and. (mod(iter, max(nint(niter/10.0), 1)) == 0 .or. iter == 1)) then
            call warn(' >> Isotropic TV iteration '//tidy(num2str(iter, '(i)')) &
                //' of '//tidy(num2str(niter, '(i)')//' relative l2-norm difference = ' &
                //tidy(num2str(norm2(u - pu), '(es)'))))
        end if

    end do

end function tv_iso_filt_2d_

!
!> 3D isotropic TV denoising
!
function tv_iso_filt_3d_(f, mu, niter, verbose) result(u)

    TT, dimension(:, :, :), intent(in) :: f
    TT, intent(in) :: mu
    integer, intent(in) :: niter
    logical, optional :: verbose
    TT, allocatable, dimension(:, :, :) :: u

    integer :: i, j, k, iter
    TT :: tmpx, tmpy, tmpz, sumu
    TT, allocatable, dimension(:, :, :) :: d1, d2, d3, d1t, d2t, d3t
    TT, allocatable, dimension(:, :, :) :: pu
    integer :: n1, n2, n3
    integer :: i1, i2, j1, j2, k1, k2
    integer :: ii1, ii2, jj1, jj2, kk1, kk2
    TT :: lambda, sk
    logical :: tv_verbose

    n1 = size(f, 1)
    n2 = size(f, 2)
    n3 = size(f, 3)

    if (present(verbose)) then
        tv_verbose = verbose
    else
        tv_verbose = .true.
    end if

    ! memory
    d1 = zeros(n1, n2, n3)
    d2 = zeros(n1, n2, n3)
    d3 = zeros(n1, n2, n3)
    d1t = zeros(n1, n2, n3)
    d2t = zeros(n1, n2, n3)
    d3t = zeros(n1, n2, n3)

    ! read noisy image
    u = f
    pu = zeros_like(u)

    ! lambdas for u and w subproblems
    lambda = 2.0*mu

    ! TGpV iteration
    do iter = 1, niter

        pu = u

        ! minimization u
        do k = 1, n3
            do j = 1, n2
                do i = 1, n1

                    if (i == 1) then
                        i1 = 0
                        ii1 = i
                    else
                        i1 = 1
                        ii1 = i - 1
                    end if

                    if (i == n1) then
                        i2 = 0
                        ii2 = i
                    else
                        i2 = 1
                        ii2 = i + 1
                    end if

                    if (j == 1) then
                        j1 = 0
                        jj1 = j
                    else
                        j1 = 1
                        jj1 = j - 1
                    end if
                    if (j == n2) then
                        j2 = 0
                        jj2 = j
                    else
                        j2 = 1
                        jj2 = j + 1
                    end if

                    if (k == 1) then
                        k1 = 0
                        kk1 = k
                    else
                        k1 = 1
                        kk1 = k - 1
                    end if
                    if (k == n3) then
                        k2 = 0
                        kk2 = k
                    else
                        k2 = 1
                        kk2 = k + 1
                    end if

                    sumu = lambda*( &
                        +i2*u(ii2, j, k) + i1*u(ii1, j, k) &
                        + j2*u(i, jj2, k) + j1*u(i, jj1, k) &
                        + k2*u(i, j, kk2) + k1*u(i, j, kk1)) &
                        + lambda*( &
                        +(i1*d1(ii1, j, k) - i2*d1(i, j, k)) &
                        + (j1*d2(i, jj1, k) - j2*d2(i, j, k)) &
                        + (k1*d3(i, j, kk1) - k2*d3(i, j, k)) &
                        - (i1*d1t(ii1, j, k) - i2*d1t(i, j, k)) &
                        - (j1*d2t(i, jj1, k) - j2*d2t(i, j, k)) &
                        - (k1*d3t(i, j, kk1) - k2*d3t(i, j, k))) &
                        + mu*f(i, j, k)

                    u(i, j, k) = sumu/(mu + (6.0 - &
                        (1 - i1) - (1 - i2) - (1 - j1) - (1 - j2) - (1 - k1) - (1 - k2))*lambda)

                end do
            end do
        end do

        ! update multiplietgpv_filt_2d_rs
        ! minimization d
        !$omp parallel do private(i, j, k, i1, i2, j1, j2, k1, k2, tmpx, tmpy, tmpz, sk)
        do k = 1, n3
            do j = 1, n2
                do i = 1, n1

                    if (i == n1) then
                        i1 = i - 1
                        i2 = i
                    else
                        i1 = i
                        i2 = i + 1
                    end if
                    tmpx = (u(i2, j, k) - u(i1, j, k)) + d1t(i, j, k)

                    if (j == n2) then
                        j1 = j - 1
                        j2 = j
                    else
                        j1 = j
                        j2 = j + 1
                    end if
                    tmpy = (u(i, j2, k) - u(i, j1, k)) + d2t(i, j, k)

                    if (k == n3) then
                        k1 = k - 1
                        k2 = k
                    else
                        k1 = k
                        k2 = k + 1
                    end if
                    tmpz = (u(i, j, k2) - u(i, j, k1)) + d3t(i, j, k)

                    sk = sqrt(tmpx**2 + tmpy**2 + tmpz**2)

                    d1(i, j, k) = max(sk - 1.0/lambda, 0.0)*tmpx/(sk + float_tiny)
                    d2(i, j, k) = max(sk - 1.0/lambda, 0.0)*tmpy/(sk + float_tiny)
                    d3(i, j, k) = max(sk - 1.0/lambda, 0.0)*tmpz/(sk + float_tiny)

                end do
            end do
        end do
        !$omp end parallel do

        ! d
        !$omp parallel do private(i, j, k, i1, i2, j1, j2, k1, k2)
        do k = 1, n3
            do j = 1, n2
                do i = 1, n1

                    if (i == n1) then
                        i1 = i - 1
                        i2 = i
                    else
                        i1 = i
                        i2 = i + 1
                    end if
                    d1t(i, j, k) = d1t(i, j, k) + (u(i2, j, k) - u(i1, j, k) - d1(i, j, k))

                    if (j == n2) then
                        j1 = j - 1
                        j2 = j
                    else
                        j1 = j
                        j2 = j + 1
                    end if
                    d2t(i, j, k) = d2t(i, j, k) + (u(i, j2, k) - u(i, j1, k) - d2(i, j, k))

                    if (k == n3) then
                        k1 = k - 1
                        k2 = k
                    else
                        k1 = k
                        k2 = k + 1
                    end if
                    d3t(i, j, k) = d3t(i, j, k) + (u(i, j, k2) - u(i, j, k1) - d3(i, j, k))

                end do
            end do
        end do
        !$omp end parallel do

        ! progress
        if (tv_verbose .and. (mod(iter, max(nint(niter/10.0), 1)) == 0 .or. iter == 1)) then
            call warn(' >> Isotropic TV iteration '//tidy(num2str(iter, '(i)')) &
                //' of '//tidy(num2str(niter, '(i)')//' relative l2-norm difference = ' &
                //tidy(num2str(norm2(u - pu), '(es)'))))
        end if

    end do

end function tv_iso_filt_3d_

!
!> 1D TGpV denoising
!
function tgpv_filt_1d_(f, mu, alpha0, alpha1, niter, p, verbose) result(u)

    TT, dimension(:), intent(in) :: f
    TT, intent(in) :: mu, alpha0, alpha1, p
    integer, intent(in) :: niter
    logical, optional :: verbose
    TT, allocatable, dimension(:) :: u

    integer :: i, iter
    TT :: tmp, sumu, sumw1
    TT, allocatable, dimension(:) :: d1, d1t
    TT, allocatable, dimension(:) :: s11
    TT, allocatable, dimension(:) :: s11t
    TT, allocatable, dimension(:) :: w1
    TT, allocatable, dimension(:) :: pu
    integer :: n1
    integer :: inner, ninner
    integer :: i1, i2
    integer :: ii1, ii2
    TT :: lambda0, lambda1
    logical :: tv_verbose

    n1 = size(f)
    ninner = 2

    if (present(verbose)) then
        tv_verbose = verbose
    else
        tv_verbose = .true.
    end if

    ! memory
    d1 = zeros(n1)
    d1t = zeros(n1)
    s11 = zeros(n1)
    s11t = zeros(n1)
    w1 = zeros(n1)
    pu = zeros(n1)

    ! read noisy image
    u = f

    ! lambdas for u and w subproblems
    lambda0 = 2.0*mu
    lambda1 = 2.0*mu*(alpha1/alpha0)

    ! TGpV iteration
    do iter = 1, niter

        pu = u

        do inner = 1, ninner

            ! minimization u
            do i = 1, n1

                if (i == 1) then
                    i1 = 0
                    ii1 = i
                else
                    i1 = 1
                    ii1 = i - 1
                end if

                if (i == n1) then
                    i2 = 0
                    ii2 = i
                else
                    i2 = 1
                    ii2 = i + 1
                end if

                sumu = lambda0*( &
                    +i2*u(ii2) + i1*u(ii1)) &
                    + lambda0*( &
                    +(i1*d1(ii1) - i2*d1(i)) &
                    - (i1*d1t(ii1) - i2*d1t(i)) &
                    + (i1*w1(ii1) - i2*w1(i)) &
                    ) &
                    + mu*f(i)

                u(i) = sumu/(mu + (2.0 - (1 - i1) - (1 - i2))*lambda0)

            end do

            if (alpha1 /= 0) then

                ! minimization of w
                do i = 1, n1

                    if (i == 1) then
                        i1 = 0
                        ii1 = i
                    else
                        i1 = 1
                        ii1 = i - 1
                    end if

                    if (i == n1) then
                        i2 = 0
                        ii2 = i
                    else
                        i2 = 1
                        ii2 = i + 1
                    end if

                    ! w1
                    sumw1 = &
                        +lambda1*(i2*w1(ii2) + i1*w1(ii1)) &
                        - lambda0*(d1(i) - d1t(i) - (i2*u(ii2) - u(i))) &
                        + lambda1*(i1*s11(ii1) - i2*s11(i) - (i1*s11t(ii1) - i2*s11t(i)))

                    w1(i) = sumw1/(lambda0 + (2.0 - (1 - i1) - (1 - i2))*lambda1)

                end do

            end if

            ! update multipliers
            ! minimization d
            !$omp parallel do private(i, i1, i2, tmp)
            do i = 1, n1

                if (i == n1) then
                    i1 = i - 1
                    i2 = i
                else
                    i1 = i
                    i2 = i + 1
                end if

                tmp = (u(i2) - u(i1)) - w1(i) + d1t(i)
                d1(i) = max(1.0 - (lambda0/alpha0*abs(tmp))**(p - 2), 0.0)*tmp

            end do
            !$omp end parallel do

            if (alpha1 /= 0) then

                !$omp parallel do private(i, i1, i2, tmp)
                do i = 1, n1

                    if (i == n1) then
                        i1 = i - 1
                        i2 = i
                    else
                        i1 = i
                        i2 = i + 1
                    end if

                    tmp = w1(i2) - w1(i1) + s11t(i)
                    s11(i) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp

                end do
                !$omp end parallel do

            end if

        end do

        ! d
        !$omp parallel do private(i, i1, i2)
        do i = 1, n1

            if (i == n1) then
                i1 = i - 1
                i2 = i
            else
                i1 = i
                i2 = i + 1
            end if

            d1t(i) = d1t(i) + (u(i2) - u(i1) - w1(i) - d1(i))

        end do
        !$omp end parallel do

        if (alpha1 /= 0) then

            !$omp parallel do private(i, i1, i2)
            do i = 1, n1

                if (i == n1) then
                    i1 = i - 1
                    i2 = i
                else
                    i1 = i
                    i2 = i + 1
                end if

                s11t(i) = s11t(i) + (w1(i2) - w1(i1)) - s11(i)

            end do
            !$omp end parallel do

        end if

        ! progress
        if (tv_verbose .and. (mod(iter, max(nint(niter/10.0), 1)) == 0 .or. iter == 1)) then
            call warn(' >> TGpV iteration '//tidy(num2str(iter, '(i)')) &
                //' of '//tidy(num2str(niter, '(i)')//' relative l2-norm difference = ' &
                //tidy(num2str(norm2(u - pu), '(es)'))))
        end if

    end do

end function tgpv_filt_1d_

!
!> 2D TGpV denoising
!
function tgpv_filt_2d_(f, mu, alpha0, alpha1, niter, p, verbose) result(u)

    TT, dimension(:, :), intent(in) :: f
    TT, intent(in) :: mu, alpha0, alpha1, p
    integer, intent(in) :: niter
    logical, optional :: verbose
    TT, allocatable, dimension(:, :) :: u

    integer :: i, j, iter
    TT :: tmp, sumu, sumw1, sumw2
    TT, allocatable, dimension(:, :) :: d1, d2, d1t, d2t
    TT, allocatable, dimension(:, :) :: s11, s12, s22
    TT, allocatable, dimension(:, :) :: s11t, s12t, s22t
    TT, allocatable, dimension(:, :) :: w1, w2
    TT, allocatable, dimension(:, :) :: tmpd1w2
    TT, allocatable, dimension(:, :) :: tmpd2w1
    TT, allocatable, dimension(:, :) :: pu
    integer :: n1, n2
    integer :: inner, ninner
    integer :: i1, i2, j1, j2
    integer :: ii1, ii2, jj1, jj2
    TT :: lambda0, lambda1
    logical :: tv_verbose

    n1 = size(f, 1)
    n2 = size(f, 2)
    ninner = 2

    if (present(verbose)) then
        tv_verbose = verbose
    else
        tv_verbose = .true.
    end if

    ! memory
    d1 = zeros(n1, n2)
    d2 = zeros(n1, n2)
    d1t = zeros(n1, n2)
    d2t = zeros(n1, n2)
    s11 = zeros(n1, n2)
    s12 = zeros(n1, n2)
    s22 = zeros(n1, n2)
    s11t = zeros(n1, n2)
    s12t = zeros(n1, n2)
    s22t = zeros(n1, n2)
    w1 = zeros(n1, n2)
    w2 = zeros(n1, n2)
    tmpd1w2 = zeros(n1, n2)
    tmpd2w1 = zeros(n1, n2)
    pu = zeros(n1, n2)

    ! read noisy image
    u = f

    ! lambdas for u and w subproblems
    lambda0 = 2.0*mu
    lambda1 = 2.0*mu*(alpha1/alpha0)

    ! TGpV iteration
    do iter = 1, niter

        pu = u

        do inner = 1, ninner

            ! minimization u
            do j = 1, n2
                do i = 1, n1

                    if (i == 1) then
                        i1 = 0
                        ii1 = i
                    else
                        i1 = 1
                        ii1 = i - 1
                    end if

                    if (i == n1) then
                        i2 = 0
                        ii2 = i
                    else
                        i2 = 1
                        ii2 = i + 1
                    end if

                    if (j == 1) then
                        j1 = 0
                        jj1 = j
                    else
                        j1 = 1
                        jj1 = j - 1
                    end if
                    if (j == n2) then
                        j2 = 0
                        jj2 = j
                    else
                        j2 = 1
                        jj2 = j + 1
                    end if

                    sumu = lambda0*( &
                        +i2*u(ii2, j) + i1*u(ii1, j) &
                        + j2*u(i, jj2) + j1*u(i, jj1)) &
                        + lambda0*( &
                        +(i1*d1(ii1, j) - i2*d1(i, j)) &
                        + (j1*d2(i, jj1) - j2*d2(i, j)) &
                        - (i1*d1t(ii1, j) - i2*d1t(i, j)) &
                        - (j1*d2t(i, jj1) - j2*d2t(i, j)) &
                        + (i1*w1(ii1, j) - i2*w1(i, j)) &
                        + (j1*w2(i, jj1) - j2*w2(i, j)) &
                        ) &
                        + mu*f(i, j)

                    u(i, j) = sumu/(mu + (4.0 - (1 - i1) - (1 - i2) - (1 - j1) - (1 - j2))*lambda0)

                end do
            end do

            if (alpha1 /= 0) then

                ! minimization w
                !$omp parallel do private(i, j, i1, i2, j1, j2)
                do j = 1, n2
                    do i = 1, n1

                        if (i == n1) then
                            i1 = i - 1
                            i2 = i
                        else
                            i1 = i
                            i2 = i + 1
                        end if

                        if (j == n2) then
                            j1 = j - 1
                            j2 = j
                        else
                            j1 = j
                            j2 = j + 1
                        end if

                        tmpd1w2(i, j) = 0.5*(w2(i2, j) - w2(i1, j))
                        tmpd2w1(i, j) = 0.5*(w1(i, j2) - w1(i, j1))

                    end do
                end do
                !$omp end parallel do

                ! minimization of w
                do j = 1, n2
                    do i = 1, n1

                        if (i == 1) then
                            i1 = 0
                            ii1 = i
                        else
                            i1 = 1
                            ii1 = i - 1
                        end if

                        if (i == n1) then
                            i2 = 0
                            ii2 = i
                        else
                            i2 = 1
                            ii2 = i + 1
                        end if

                        if (j == 1) then
                            j1 = 0
                            jj1 = j
                        else
                            j1 = 1
                            jj1 = j - 1
                        end if
                        if (j == n2) then
                            j2 = 0
                            jj2 = j
                        else
                            j2 = 1
                            jj2 = j + 1
                        end if

                        ! w1
                        sumw1 = &
                            +lambda1*(i2*w1(ii2, j) + i1*w1(ii1, j)) &
                            + 0.5*lambda1*(j2*w1(i, jj2) + j1*w1(i, jj1)) &
                            - lambda0*(d1(i, j) - d1t(i, j) - (i2*u(ii2, j) - u(i, j))) &
                            + lambda1*(i1*s11(ii1, j) - i2*s11(i, j) - (i1*s11t(ii1, j) - i2*s11t(i, j))) &
                            + lambda1*(j1*s12(i, jj1) - j2*s12(i, j) - (j1*s12t(i, jj1) - j2*s12t(i, j)) - (j1*tmpd1w2(i, jj1) - j2*tmpd1w2(i, j)))

                        w1(i, j) = sumw1/(lambda0 + (2.0 - (1 - i1) - (1 - i2))*lambda1 + (2.0 - (1 - j1) - (1 - j2))*0.5*lambda1)

                        ! w2
                        sumw2 = &
                            +0.5*lambda1*(i2*w2(ii2, j) + i1*w2(ii1, j)) &
                            + lambda1*(j2*w2(i, jj2) + j1*w2(i, jj1)) &
                            - lambda0*(d2(i, j) - d2t(i, j) - (j2*u(i, jj2) - u(i, j))) &
                            + lambda1*(i1*s12(ii1, j) - i2*s12(i, j) - (i1*s12t(ii1, j) - i2*s12t(i, j)) - (i1*tmpd2w1(ii1, j) - i2*tmpd2w1(i, j))) &
                            + lambda1*(j1*s22(i, jj1) - j2*s22(i, j) - (j1*s22t(i, jj1) - j2*s22t(i, j)))

                        w2(i, j) = sumw2/(lambda0 + (2.0 - (1 - i1) - (1 - i2))*0.5*lambda1 + (2.0 - (1 - j1) - (1 - j2))*lambda1)

                    end do
                end do

            end if

            ! update multipliers
            ! minimization d
            !$omp parallel do private(i, j, i1, i2, j1, j2, tmp)
            do j = 1, n2
                do i = 1, n1

                    if (i == n1) then
                        i1 = i - 1
                        i2 = i
                    else
                        i1 = i
                        i2 = i + 1
                    end if
                    tmp = (u(i2, j) - u(i1, j)) - w1(i, j) + d1t(i, j)
                    d1(i, j) = max(1.0 - (lambda0/alpha0*abs(tmp))**(p - 2), 0.0)*tmp

                    if (j == n2) then
                        j1 = j - 1
                        j2 = j
                    else
                        j1 = j
                        j2 = j + 1
                    end if
                    tmp = (u(i, j2) - u(i, j1)) - w2(i, j) + d2t(i, j)
                    d2(i, j) = max(1.0 - (lambda0/alpha0*abs(tmp))**(p - 2), 0.0)*tmp

                end do
            end do
            !$omp end parallel do

            if (alpha1 /= 0) then

                !$omp parallel do private(i, j, i1, i2, j1, j2, tmp)
                do j = 1, n2
                    do i = 1, n1

                        if (i == n1) then
                            i1 = i - 1
                            i2 = i
                        else
                            i1 = i
                            i2 = i + 1
                        end if
                        if (j == n2) then
                            j1 = j - 1
                            j2 = j
                        else
                            j1 = j
                            j2 = j + 1
                        end if

                        tmp = w1(i2, j) - w1(i1, j) + s11t(i, j)
                        s11(i, j) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp
                        tmp = w2(i, j2) - w2(i, j1) + s22t(i, j)
                        s22(i, j) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp
                        tmp = 0.5*(w2(i2, j) - w2(i1, j) + w1(i, j2) - w1(i, j1)) + s12t(i, j)
                        s12(i, j) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp

                    end do
                end do
                !$omp end parallel do

            end if

        end do

        ! d
        !$omp parallel do private(i, j, i1, i2, j1, j2)
        do j = 1, n2
            do i = 1, n1

                if (i == n1) then
                    i1 = i - 1
                    i2 = i
                else
                    i1 = i
                    i2 = i + 1
                end if
                d1t(i, j) = d1t(i, j) + (u(i2, j) - u(i1, j) - w1(i, j) - d1(i, j))

                if (j == n2) then
                    j1 = j - 1
                    j2 = j
                else
                    j1 = j
                    j2 = j + 1
                end if
                d2t(i, j) = d2t(i, j) + (u(i, j2) - u(i, j1) - w2(i, j) - d2(i, j))

            end do
        end do
        !$omp end parallel do

        if (alpha1 /= 0) then

            !$omp parallel do private(i, j, i1, i2, j1, j2, tmp)
            do j = 1, n2
                do i = 1, n1

                    if (i == n1) then
                        i1 = i - 1
                        i2 = i
                    else
                        i1 = i
                        i2 = i + 1
                    end if
                    if (j == n2) then
                        j1 = j - 1
                        j2 = j
                    else
                        j1 = j
                        j2 = j + 1
                    end if

                    s11t(i, j) = s11t(i, j) + (w1(i2, j) - w1(i1, j)) - s11(i, j)
                    s22t(i, j) = s22t(i, j) + (w2(i, j2) - w2(i, j1)) - s22(i, j)
                    s12t(i, j) = s12t(i, j) + 0.5*(w1(i, j2) - w1(i, j1) + w2(i2, j) - w2(i1, j)) - s12(i, j)

                end do
            end do
            !$omp end parallel do

        end if

        ! progress
        if (tv_verbose .and. (mod(iter, max(nint(niter/10.0), 1)) == 0 .or. iter == 1)) then
            call warn(' >> TGpV iteration '//tidy(num2str(iter, '(i)')) &
                //' of '//tidy(num2str(niter, '(i)')//' relative l2-norm difference = ' &
                //tidy(num2str(norm2(u - pu), '(es)'))))
        end if

    end do

end function tgpv_filt_2d_

!
!> 3D TGpV denoising
!
function tgpv_filt_3d_(f, mu, alpha0, alpha1, niter, p) result(u)

    TT, dimension(:, :, :), intent(in) :: f
    TT, intent(in) :: mu, alpha0, alpha1, p
    integer, intent(in) :: niter
    TT, allocatable, dimension(:, :, :) :: u

    integer :: i, j, k
    TT :: tmp, sumu, sumw1, sumw2, sumw3
    TT, allocatable, dimension(:, :, :) ::  d1, d2, d3, d1t, d2t, d3t
    TT, allocatable, dimension(:, :, :) :: s11, s12, s13, s22, s23, s33
    TT, allocatable, dimension(:, :, :) :: s11t, s12t, s13t, s22t, s23t, s33t
    TT, allocatable, dimension(:, :, :) :: w1, w2, w3
    TT, allocatable, dimension(:, :, :) :: tmpd1w2, tmpd1w3
    TT, allocatable, dimension(:, :, :) :: tmpd2w1, tmpd2w3
    TT, allocatable, dimension(:, :, :) :: tmpd3w1, tmpd3w2
    TT, allocatable, dimension(:, :, :) :: pu
    integer :: n1, n2, n3
    integer :: iter, inner, ninner
    integer :: i1, i2, j1, j2, k1, k2
    integer :: ii1, ii2, jj1, jj2, kk1, kk2
    TT :: lambda0, lambda1
    TT :: bsum

    ! number of inner iteration
    ninner = 2

    ! lambdas for u and w subproblems
    lambda0 = 2.0*mu
    lambda1 = 2.0*mu*(alpha1/alpha0)

    n1 = size(f, 1)
    n2 = size(f, 2)
    n3 = size(f, 3)

    ! memory
    d1 = zeros(n1, n2, n3)
    d2 = zeros(n1, n2, n3)
    d3 = zeros(n1, n2, n3)
    d1t = zeros(n1, n2, n3)
    d2t = zeros(n1, n2, n3)
    d3t = zeros(n1, n2, n3)
    s11 = zeros(n1, n2, n3)
    s12 = zeros(n1, n2, n3)
    s13 = zeros(n1, n2, n3)
    s22 = zeros(n1, n2, n3)
    s23 = zeros(n1, n2, n3)
    s33 = zeros(n1, n2, n3)
    s11t = zeros(n1, n2, n3)
    s12t = zeros(n1, n2, n3)
    s13t = zeros(n1, n2, n3)
    s22t = zeros(n1, n2, n3)
    s23t = zeros(n1, n2, n3)
    s33t = zeros(n1, n2, n3)
    w1 = zeros(n1, n2, n3)
    w2 = zeros(n1, n2, n3)
    w3 = zeros(n1, n2, n3)
    u = zeros(n1, n2, n3)
    tmpd1w2 = zeros(n1, n2, n3)
    tmpd1w3 = zeros(n1, n2, n3)
    tmpd2w1 = zeros(n1, n2, n3)
    tmpd2w3 = zeros(n1, n2, n3)
    tmpd3w1 = zeros(n1, n2, n3)
    tmpd3w2 = zeros(n1, n2, n3)
    pu = zeros(n1, n2, n3)

    ! noisy image
    u = f

    ! TGpV iteration
    do iter = 1, niter

        pu = u

        do inner = 1, ninner

            ! Minimization of u
            do k = 1, n3
                do j = 1, n2
                    do i = 1, n1

                        if (i == 1) then
                            i1 = 0
                            ii1 = i
                        else
                            i1 = 1
                            ii1 = i - 1
                        end if

                        if (i == n1) then
                            i2 = 0
                            ii2 = i
                        else
                            i2 = 1
                            ii2 = i + 1
                        end if

                        if (j == 1) then
                            j1 = 0
                            jj1 = j
                        else
                            j1 = 1
                            jj1 = j - 1
                        end if
                        if (j == n2) then
                            j2 = 0
                            jj2 = j
                        else
                            j2 = 1
                            jj2 = j + 1
                        end if

                        if (k == 1) then
                            k1 = 0
                            kk1 = k
                        else
                            k1 = 1
                            kk1 = k - 1
                        end if
                        if (k == n3) then
                            k2 = 0
                            kk2 = k
                        else
                            k2 = 1
                            kk2 = k + 1
                        end if

                        sumu = lambda0*( &
                            +i2*u(ii2, j, k) + i1*u(ii1, j, k) &
                            + j2*u(i, jj2, k) + j1*u(i, jj1, k) &
                            + k2*u(i, j, kk2) + k1*u(i, j, kk1)) &
                            + lambda0*( &
                            +(i1*d1(ii1, j, k) - i2*d1(i, j, k)) &
                            + (j1*d2(i, jj1, k) - j2*d2(i, j, k)) &
                            + (k1*d3(i, j, kk1) - k2*d3(i, j, k)) &
                            - (i1*d1t(ii1, j, k) - i2*d1t(i, j, k)) &
                            - (j1*d2t(i, jj1, k) - j2*d2t(i, j, k)) &
                            - (k1*d3t(i, j, kk1) - k2*d3t(i, j, k)) &
                            + (i1*w1(ii1, j, k) - i2*w1(i, j, k)) &
                            + (j1*w2(i, jj1, k) - j2*w2(i, j, k)) &
                            + (k1*w3(i, j, kk1) - k2*w3(i, j, k)) &
                            ) &
                            + mu*f(i, j, k)

                        u(i, j, k) = sumu/(mu + (6.0 - (1 - i1) - (1 - i2) - (1 - j1) - (1 - j2) - (1 - k1) - (1 - k2))*lambda0)

                    end do
                end do
            end do

            if (alpha1 /= 0) then

                ! minimization w
                do k = 1, n3
                    do j = 1, n2
                        do i = 1, n1

                            if (i == n1) then
                                i1 = i - 1
                                i2 = i
                            else
                                i1 = i
                                i2 = i + 1
                            end if

                            if (j == n2) then
                                j1 = j - 1
                                j2 = j
                            else
                                j1 = j
                                j2 = j + 1
                            end if

                            if (k == n3) then
                                k1 = k - 1
                                k2 = k
                            else
                                k1 = k
                                k2 = k + 1
                            end if

                            tmpd1w2(i, j, k) = 0.5*(w2(i2, j, k) - w2(i1, j, k))
                            tmpd1w3(i, j, k) = 0.5*(w3(i2, j, k) - w3(i1, j, k))
                            tmpd2w1(i, j, k) = 0.5*(w1(i, j2, k) - w1(i, j1, k))
                            tmpd2w3(i, j, k) = 0.5*(w3(i, j2, k) - w3(i, j1, k))
                            tmpd3w1(i, j, k) = 0.5*(w1(i, j, k2) - w1(i, j, k1))
                            tmpd3w2(i, j, k) = 0.5*(w2(i, j, k2) - w2(i, j, k1))

                        end do
                    end do
                end do

                ! minimization of w
                do k = 1, n3
                    do j = 1, n2
                        do i = 1, n1

                            if (i == 1) then
                                i1 = 0
                                ii1 = i
                            else
                                i1 = 1
                                ii1 = i - 1
                            end if

                            if (i == n1) then
                                i2 = 0
                                ii2 = i
                            else
                                i2 = 1
                                ii2 = i + 1
                            end if

                            if (j == 1) then
                                j1 = 0
                                jj1 = j
                            else
                                j1 = 1
                                jj1 = j - 1
                            end if
                            if (j == n2) then
                                j2 = 0
                                jj2 = j
                            else
                                j2 = 1
                                jj2 = j + 1
                            end if

                            if (k == 1) then
                                k1 = 0
                                kk1 = k
                            else
                                k1 = 1
                                kk1 = k - 1
                            end if
                            if (k == n3) then
                                k2 = 0
                                kk2 = k
                            else
                                k2 = 1
                                kk2 = k + 1
                            end if

                            ! w1
                            sumw1 = &
                                +lambda1*(i2*w1(ii2, j, k) + i1*w1(ii1, j, k)) &
                                + 0.5*lambda1*(j2*w1(i, jj2, k) + j1*w1(i, jj1, k)) &
                                + 0.5*lambda1*(k2*w1(i, j, kk2) + k1*w1(i, j, kk1)) &
                                - lambda0*(d1(i, j, k) - d1t(i, j, k) - (i2*u(ii2, j, k) - u(i, j, k))) &
                                + lambda1*(i1*s11(ii1, j, k) - i2*s11(i, j, k) - (i1*s11t(ii1, j, k) - i2*s11t(i, j, k))) &
                                + lambda1*(j1*s12(i, jj1, k) - j2*s12(i, j, k) - (j1*s12t(i, jj1, k) - j2*s12t(i, j, k)) - (j1*tmpd1w2(i, jj1, k) - j2*tmpd1w2(i, j, k))) &
                                + lambda1*(k1*s13(i, j, kk1) - k2*s13(i, j, k) - (k1*s13t(i, j, kk1) - k2*s13t(i, j, k)) - (k1*tmpd1w3(i, j, kk1) - k2*tmpd1w3(i, j, k)))

                            w1(i, j, k) = sumw1/(lambda0 + (2.0 - (1-i1) - (1-i2))*lambda1 + (2.0 - (1-j1) - (1-j2))*0.5*lambda1 + (2.0 - (1-k1) - (1-k2))*0.5*lambda1)

                            ! w2
                            sumw2 = &
                                +0.5*lambda1*(i2*w2(ii2, j, k) + i1*w2(ii1, j, k)) &
                                + lambda1*(j2*w2(i, jj2, k) + j1*w2(i, jj1, k)) &
                                + 0.5*lambda1*(k2*w2(i, j, kk2) + k1*w2(i, j, kk1)) &
                                - lambda0*(d2(i, j, k) - d2t(i, j, k) - (j2*u(i, jj2, k) - u(i, j, k))) &
                                + lambda1*(i1*s12(ii1, j, k) - i2*s12(i, j, k) - (i1*s12t(ii1, j, k) - i2*s12t(i, j, k)) - (i1*tmpd2w1(ii1, j, k) - i2*tmpd2w1(i, j, k))) &
                                + lambda1*(j1*s22(i, jj1, k) - j2*s22(i, j, k) - (j1*s22t(i, jj1, k) - j2*s22t(i, j, k))) &
                                + lambda1*(k1*s23(i, j, kk1) - k2*s23(i, j, k) - (k1*s23t(i, j, kk1) - k2*s23t(i, j, k)) - (k1*tmpd2w3(i, j, kk1) - k2*tmpd2w3(i, j, k)))

                            w2(i, j, k) = sumw2/(lambda0 + (2.0 - (1-i1) - (1-i2))*0.5*lambda1 + (2.0 - (1-j1) - (1-j2))*lambda1 + (2.0 - (1-k1) - (1-k2))*0.5*lambda1)

                            ! w3
                            sumw3 = &
                                +0.5*lambda1*(i2*w3(ii2, j, k) + i1*w3(ii1, j, k)) &
                                + 0.5*lambda1*(j2*w3(i, jj2, k) + j1*w3(i, jj1, k)) &
                                + lambda1*(k2*w3(i, j, kk2) + k1*w3(i, j, kk1)) &
                                - lambda0*(d3(i, j, k) - d3t(i, j, k) - (k2*u(i, j, kk2) - u(i, j, k))) &
                                + lambda1*(i1*s13(ii1, j, k) - i2*s13(i, j, k) - (i1*s13t(ii1, j, k) - i2*s13t(i, j, k)) - (i1*tmpd3w1(ii1, j, k) - i2*tmpd3w1(i, j, k))) &
                                + lambda1*(j1*s23(i, jj1, k) - j2*s23(i, j, k) - (j1*s23t(i, jj1, k) - j2*s23t(i, j, k)) - (j1*tmpd3w2(i, jj1, k) - j2*tmpd3w2(i, j, k))) &
                                + lambda1*(k1*s33(i, j, kk1) - k2*s33(i, j, k) - (k1*s33t(i, j, kk1) - k2*s33t(i, j, k)))

                            w3(i, j, k) = sumw3/(lambda0 + (2.0 - (1-i1) - (1-i2))*0.5*lambda1 + (2.0 - (1-j1) - (1-j2))*0.5*lambda1 + (2.0 - (1-k1) - (1-k2))*lambda1)

                        end do
                    end do
                end do

            end if

            ! update multipliers
            ! minimization d
            do k = 1, n3
                do j = 1, n2
                    do i = 1, n1

                        if (i == n1) then
                            i1 = i - 1
                            i2 = i
                        else
                            i1 = i
                            i2 = i + 1
                        end if
                        tmp = (u(i2, j, k) - u(i1, j, k)) - w1(i, j, k) + d1t(i, j, k)
                        d1(i, j, k) = max(1.0 - (lambda0/alpha0*abs(tmp))**(p - 2), 0.0)*tmp

                        if (j == n2) then
                            j1 = j - 1
                            j2 = j
                        else
                            j1 = j
                            j2 = j + 1
                        end if
                        tmp = (u(i, j2, k) - u(i, j1, k)) - w2(i, j, k) + d2t(i, j, k)
                        d2(i, j, k) = max(1.0 - (lambda0/alpha0*abs(tmp))**(p - 2), 0.0)*tmp

                        if (k == n3) then
                            k1 = k - 1
                            k2 = k
                        else
                            k1 = k
                            k2 = k + 1
                        end if
                        tmp = (u(i, j, k2) - u(i, j, k1)) - w3(i, j, k) + d3t(i, j, k)
                        d3(i, j, k) = max(1.0 - (lambda0/alpha0*abs(tmp))**(p - 2), 0.0)*tmp

                    end do
                end do
            end do

            if (alpha1 /= 0) then

                do k = 1, n3
                    do j = 1, n2
                        do i = 1, n1

                            if (i == n1) then
                                i1 = i - 1
                                i2 = i
                            else
                                i1 = i
                                i2 = i + 1
                            end if
                            if (j == n2) then
                                j1 = j - 1
                                j2 = j
                            else
                                j1 = j
                                j2 = j + 1
                            end if
                            if (k == n3) then
                                k1 = k - 1
                                k2 = k
                            else
                                k1 = k
                                k2 = k + 1
                            end if

                            tmp = w1(i2, j, k) - w1(i1, j, k) + s11t(i, j, k)
                            s11(i, j, k) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp
                            tmp = w2(i, j2, k) - w2(i, j1, k) + s22t(i, j, k)
                            s22(i, j, k) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp
                            tmp = w3(i, j, k2) - w3(i, j, k1) + s33t(i, j, k)
                            s33(i, j, k) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp
                            tmp = 0.5*(w2(i2, j, k) - w2(i1, j, k) + w1(i, j2, k) - w1(i, j1, k)) + s12t(i, j, k)
                            s12(i, j, k) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp
                            tmp = 0.5*(w3(i2, j, k) - w3(i1, j, k) + w1(i, j, k2) - w1(i, j, k1)) + s13t(i, j, k)
                            s13(i, j, k) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp
                            tmp = 0.5*(w3(i, j2, k) - w3(i, j1, k) + w2(i, j, k2) - w2(i, j, k1)) + s23t(i, j, k)
                            s23(i, j, k) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp

                        end do
                    end do
                end do

            end if

        end do

        ! d
        do k = 1, n3
            do j = 1, n2
                do i = 1, n1

                    if (i == n1) then
                        i1 = i - 1
                        i2 = i
                    else
                        i1 = i
                        i2 = i + 1
                    end if
                    d1t(i, j, k) = d1t(i, j, k) + (u(i2, j, k) - u(i1, j, k) - w1(i, j, k) - d1(i, j, k))

                    if (j == n2) then
                        j1 = j - 1
                        j2 = j
                    else
                        j1 = j
                        j2 = j + 1
                    end if
                    d2t(i, j, k) = d2t(i, j, k) + (u(i, j2, k) - u(i, j1, k) - w2(i, j, k) - d2(i, j, k))

                    if (k == n3) then
                        k1 = k - 1
                        k2 = k
                    else
                        k1 = k
                        k2 = k + 1
                    end if
                    d3t(i, j, k) = d3t(i, j, k) + (u(i, j, k2) - u(i, j, k1) - w3(i, j, k) - d3(i, j, k))

                end do
            end do
        end do

        if (alpha1 /= 0) then

            do k = 1, n3
                do j = 1, n2
                    do i = 1, n1

                        if (i == n1) then
                            i1 = i - 1
                            i2 = i
                        else
                            i1 = i
                            i2 = i + 1
                        end if
                        if (j == n2) then
                            j1 = j - 1
                            j2 = j
                        else
                            j1 = j
                            j2 = j + 1
                        end if
                        if (k == n3) then
                            k1 = k - 1
                            k2 = k
                        else
                            k1 = k
                            k2 = k + 1
                        end if

                        s11t(i, j, k) = s11t(i, j, k) + (w1(i2, j, k) - w1(i1, j, k)) - s11(i, j, k)
                        s22t(i, j, k) = s22t(i, j, k) + (w2(i, j2, k) - w2(i, j1, k)) - s22(i, j, k)
                        s33t(i, j, k) = s33t(i, j, k) + (w3(i, j, k2) - w3(i, j, k1)) - s33(i, j, k)
                        s12t(i, j, k) = s12t(i, j, k) + 0.5*(w1(i, j2, k) - w1(i, j1, k) + w2(i2, j, k) - w2(i1, j, k)) - s12(i, j, k)
                        s13t(i, j, k) = s13t(i, j, k) + 0.5*(w1(i, j, k2) - w1(i, j, k1) + w3(i2, j, k) - w3(i1, j, k)) - s13(i, j, k)
                        s23t(i, j, k) = s23t(i, j, k) + 0.5*(w2(i, j, k2) - w2(i, j, k1) + w3(i, j2, k) - w3(i, j1, k)) - s23(i, j, k)

                    end do
                end do
            end do

        end if

        ! progress
        if (mod(iter, max(nint(niter/10.0), 1)) == 0 .or. iter == 1) then
            bsum = sum((u - pu)**2)
            call warn(' >> TGpV iteration '//num2str(iter, '(i)') &
                //' of '//num2str(niter, '(i)')//' relative l2-norm difference = '//num2str(sqrt(bsum), '(es)'))
        end if

    end do

end function tgpv_filt_3d_

!
!> 2D TGpV denoising MPI version
!
function tgpv_filt_2d_mpi_(f0, mu, alpha0, alpha1, niter, p) result(u0)

    TT, dimension(:, :), intent(in) :: f0
    TT, intent(in) :: mu, alpha0, alpha1, p
    integer, intent(in) :: niter
    TT, allocatable, dimension(:, :) :: u0

    integer :: i, j
    TT :: tmp, sumu, sumw1, sumw2
    TT, allocatable, dimension(:, :) ::  d1, d2, d1t, d2t
    TT, allocatable, dimension(:, :) :: s11, s12, s22
    TT, allocatable, dimension(:, :) :: s11t, s12t, s22t
    TT, allocatable, dimension(:, :) :: w1, w2, u, f
    TT, allocatable, dimension(:, :) :: tmpd1w2
    TT, allocatable, dimension(:, :) :: tmpd2w1
    TT, allocatable, dimension(:, :) :: pu
    integer :: n1, n2
    integer :: iter, inner, ninner
    integer :: i1, i2, j1, j2
    integer :: ii1, ii2, jj1, jj2
    TT :: lambda0, lambda1
    integer :: n1beg, n1end, n2beg, n2end
    TT :: bsum
    logical :: nonempty

    ! number of inner iteration
    ninner = 2

    ! lambdas for u and w subproblems
    lambda0 = 2.0*mu
    lambda1 = 2.0*mu*(alpha1/alpha0)

    n1 = size(f0, 1)
    n2 = size(f0, 2)
    call domain_decomp_regular(n1, n2, n1beg, n1end, n2beg, n2end)

    if (rankid <= rank1*rank2 - 1) then
        nonempty = .true.
    else
        nonempty = .false.
    end if

    ! memory
    call alloc_array(d1, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(d2, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(d1t, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(d2t, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(s11, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(s12, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(s22, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(s11t, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(s12t, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(s22t, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(w1, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(w2, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(u, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(f, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(tmpd1w2, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(tmpd2w1, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(pu, [n1beg, n1end, n2beg, n2end], pad=1)

    ! noisy image
    f(n1beg:n1end, n2beg:n2end) = f0(n1beg:n1end, n2beg:n2end)
    u(n1beg:n1end, n2beg:n2end) = f0(n1beg:n1end, n2beg:n2end)

    ! TGpV iteration
    do iter = 1, niter

        pu = u

        do inner = 1, ninner

            if (nonempty) then
                call commute_array(u, 1)
                call commute_array(d1, 1)
                call commute_array(d2, 1)
                call commute_array(d1t, 1)
                call commute_array(d2t, 1)
                call commute_array(w1, 1)
                call commute_array(w2, 1)
            end if

            ! minimization u
            do j = n2beg, n2end
                do i = n1beg, n1end

                    if (i == 1) then
                        i1 = 0
                        ii1 = i
                    else
                        i1 = 1
                        ii1 = i - 1
                    end if

                    if (i == n1) then
                        i2 = 0
                        ii2 = i
                    else
                        i2 = 1
                        ii2 = i + 1
                    end if

                    if (j == 1) then
                        j1 = 0
                        jj1 = j
                    else
                        j1 = 1
                        jj1 = j - 1
                    end if
                    if (j == n2) then
                        j2 = 0
                        jj2 = j
                    else
                        j2 = 1
                        jj2 = j + 1
                    end if

                    sumu = lambda0*( &
                        +i2*u(ii2, j) + i1*u(ii1, j) &
                        + j2*u(i, jj2) + j1*u(i, jj1)) &
                        + lambda0*( &
                        +(i1*d1(ii1, j) - i2*d1(i, j)) &
                        + (j1*d2(i, jj1) - j2*d2(i, j)) &
                        - (i1*d1t(ii1, j) - i2*d1t(i, j)) &
                        - (j1*d2t(i, jj1) - j2*d2t(i, j)) &
                        + (i1*w1(ii1, j) - i2*w1(i, j)) &
                        + (j1*w2(i, jj1) - j2*w2(i, j)) &
                        ) &
                        + mu*f(i, j)

                    u(i, j) = sumu/(mu + (4.0 - (1 - i1) - (1 - i2) - (1 - j1) - (1 - j2))*lambda0)

                end do
            end do

            if (alpha1 /= 0) then

                ! minimization w
                do j = n2beg, n2end
                    do i = n1beg, n1end

                        if (i == n1) then
                            i1 = i - 1
                            i2 = i
                        else
                            i1 = i
                            i2 = i + 1
                        end if

                        if (j == n2) then
                            j1 = j - 1
                            j2 = j
                        else
                            j1 = j
                            j2 = j + 1
                        end if

                        tmpd1w2(i, j) = 0.5*(w2(i2, j) - w2(i1, j))
                        tmpd2w1(i, j) = 0.5*(w1(i, j2) - w1(i, j1))

                    end do
                end do

                if (nonempty) then
                    call commute_array(u, 1)
                    call commute_array(tmpd1w2, 1)
                    call commute_array(tmpd2w1, 1)
                    call commute_array(s11, 1)
                    call commute_array(s12, 1)
                    call commute_array(s22, 1)
                    call commute_array(s11t, 1)
                    call commute_array(s12t, 1)
                    call commute_array(s22t, 1)
                end if

                ! minimization of w
                do j = n2beg, n2end
                    do i = n1beg, n1end

                        if (i == 1) then
                            i1 = 0
                            ii1 = i
                        else
                            i1 = 1
                            ii1 = i - 1
                        end if

                        if (i == n1) then
                            i2 = 0
                            ii2 = i
                        else
                            i2 = 1
                            ii2 = i + 1
                        end if

                        if (j == 1) then
                            j1 = 0
                            jj1 = j
                        else
                            j1 = 1
                            jj1 = j - 1
                        end if
                        if (j == n2) then
                            j2 = 0
                            jj2 = j
                        else
                            j2 = 1
                            jj2 = j + 1
                        end if

                        ! w1
                        sumw1 = &
                            +lambda1*(i2*w1(ii2, j) + i1*w1(ii1, j)) &
                            + 0.5*lambda1*(j2*w1(i, jj2) + j1*w1(i, jj1)) &
                            - lambda0*(d1(i, j) - d1t(i, j) - (i2*u(ii2, j) - u(i, j))) &
                            + lambda1*(i1*s11(ii1, j) - i2*s11(i, j) - (i1*s11t(ii1, j) - i2*s11t(i, j))) &
                            + lambda1*(j1*s12(i, jj1) - j2*s12(i, j) - (j1*s12t(i, jj1) - j2*s12t(i, j)) - (j1*tmpd1w2(i, jj1) - j2*tmpd1w2(i, j)))

                        w1(i, j) = sumw1/(lambda0 + (2.0 - (1 - i1) - (1 - i2))*lambda1 + (2.0 - (1 - j1) - (1 - j2))*0.5*lambda1)

                        ! w2
                        sumw2 = &
                            +0.5*lambda1*(i2*w2(ii2, j) + i1*w2(ii1, j)) &
                            + lambda1*(j2*w2(i, jj2) + j1*w2(i, jj1)) &
                            - lambda0*(d2(i, j) - d2t(i, j) - (j2*u(i, jj2) - u(i, j))) &
                            + lambda1*(i1*s12(ii1, j) - i2*s12(i, j) - (i1*s12t(ii1, j) - i2*s12t(i, j)) - (i1*tmpd2w1(ii1, j) - i2*tmpd2w1(i, j))) &
                            + lambda1*(j1*s22(i, jj1) - j2*s22(i, j) - (j1*s22t(i, jj1) - j2*s22t(i, j)))

                        w2(i, j) = sumw2/(lambda0 + (2.0 - (1 - i1) - (1 - i2))*0.5*lambda1 + (2.0 - (1 - j1) - (1 - j2))*lambda1)

                    end do
                end do

            end if

            ! update multipliers
            ! minimization d
            do j = n2beg, n2end
                do i = n1beg, n1end

                    if (i == n1) then
                        i1 = i - 1
                        i2 = i
                    else
                        i1 = i
                        i2 = i + 1
                    end if
                    tmp = (u(i2, j) - u(i1, j)) - w1(i, j) + d1t(i, j)
                    d1(i, j) = max(1.0 - (lambda0/alpha0*abs(tmp))**(p - 2), 0.0)*tmp

                    if (j == n2) then
                        j1 = j - 1
                        j2 = j
                    else
                        j1 = j
                        j2 = j + 1
                    end if
                    tmp = (u(i, j2) - u(i, j1)) - w2(i, j) + d2t(i, j)
                    d2(i, j) = max(1.0 - (lambda0/alpha0*abs(tmp))**(p - 2), 0.0)*tmp

                end do
            end do

            if (alpha1 /= 0) then

                if (nonempty) then
                    call commute_array(w1, 1)
                    call commute_array(w2, 1)
                end if

                do j = n2beg, n2end
                    do i = n1beg, n1end

                        if (i == n1) then
                            i1 = i - 1
                            i2 = i
                        else
                            i1 = i
                            i2 = i + 1
                        end if
                        if (j == n2) then
                            j1 = j - 1
                            j2 = j
                        else
                            j1 = j
                            j2 = j + 1
                        end if

                        tmp = w1(i2, j) - w1(i1, j) + s11t(i, j)
                        s11(i, j) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp
                        tmp = w2(i, j2) - w2(i, j1) + s22t(i, j)
                        s22(i, j) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp
                        tmp = 0.5*(w2(i2, j) - w2(i1, j) + w1(i, j2) - w1(i, j1)) + s12t(i, j)
                        s12(i, j) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp

                    end do
                end do

            end if

            call mpi_barrier(mpi_comm_world, mpi_ierr)

        end do

        ! d
        do j = n2beg, n2end
            do i = n1beg, n1end

                if (i == n1) then
                    i1 = i - 1
                    i2 = i
                else
                    i1 = i
                    i2 = i + 1
                end if
                d1t(i, j) = d1t(i, j) + (u(i2, j) - u(i1, j) - w1(i, j) - d1(i, j))

                if (j == n2) then
                    j1 = j - 1
                    j2 = j
                else
                    j1 = j
                    j2 = j + 1
                end if
                d2t(i, j) = d2t(i, j) + (u(i, j2) - u(i, j1) - w2(i, j) - d2(i, j))

            end do
        end do

        if (alpha1 /= 0) then

            do j = n2beg, n2end
                do i = n1beg, n1end

                    if (i == n1) then
                        i1 = i - 1
                        i2 = i
                    else
                        i1 = i
                        i2 = i + 1
                    end if
                    if (j == n2) then
                        j1 = j - 1
                        j2 = j
                    else
                        j1 = j
                        j2 = j + 1
                    end if

                    s11t(i, j) = s11t(i, j) + (w1(i2, j) - w1(i1, j)) - s11(i, j)
                    s22t(i, j) = s22t(i, j) + (w2(i, j2) - w2(i, j1)) - s22(i, j)
                    s12t(i, j) = s12t(i, j) + 0.5*(w1(i, j2) - w1(i, j1) + w2(i2, j) - w2(i1, j)) - s12(i, j)

                end do
            end do

        end if

        call mpi_barrier(mpi_comm_world, mpi_ierr)

        ! progress
        if (mod(iter, max(nint(niter/10.0), 1)) == 0 .or. iter == 1) then
            if (nonempty) then
                bsum = sum((u - pu)**2)
            else
                bsum = 0.0
            end if
            call allreduce(bsum)
            if (rankid == 0) then
                call warn(' >> TGpV iteration '//num2str(iter, '(i)') &
                    //' of '//num2str(niter, '(i)')//' relative l2-norm difference = '//num2str(sqrt(bsum), '(es)'))
            end if
        end if

    end do

    call mpi_barrier(mpi_comm_world, mpi_ierr)

    call alloc_array(u0, [1, n1, 1, n2])
    if (nonempty) then
        u0(n1beg:n1end, n2beg:n2end) = u(n1beg:n1end, n2beg:n2end)
    end if
    call allreduce_array(u0)

    call mpi_barrier(mpi_comm_world, mpi_ierr)

end function tgpv_filt_2d_mpi_

!
!> 3D TGpV denoising MPI version
!
function tgpv_filt_3d_mpi_(f0, mu, alpha0, alpha1, niter, p) result(u0)

    TT, dimension(:, :, :), intent(in) :: f0
    TT, intent(in) :: mu, alpha0, alpha1, p
    integer, intent(in) :: niter
    TT, allocatable, dimension(:, :, :) :: u0

    integer :: i, j, k
    TT :: tmp, sumu, sumw1, sumw2, sumw3
    TT, allocatable, dimension(:, :, :) ::  d1, d2, d3, d1t, d2t, d3t
    TT, allocatable, dimension(:, :, :) :: s11, s12, s13, s22, s23, s33
    TT, allocatable, dimension(:, :, :) :: s11t, s12t, s13t, s22t, s23t, s33t
    TT, allocatable, dimension(:, :, :) :: w1, w2, w3, u, f
    TT, allocatable, dimension(:, :, :) :: tmpd1w2, tmpd1w3
    TT, allocatable, dimension(:, :, :) :: tmpd2w1, tmpd2w3
    TT, allocatable, dimension(:, :, :) :: tmpd3w1, tmpd3w2
    TT, allocatable, dimension(:, :, :) :: pu
    integer :: n1, n2, n3
    integer :: iter, inner, ninner
    integer :: i1, i2, j1, j2, k1, k2
    integer :: ii1, ii2, jj1, jj2, kk1, kk2
    TT :: lambda0, lambda1
    integer :: n1beg, n1end, n2beg, n2end, n3beg, n3end
    TT :: bsum
    logical :: nonempty

    ! number of inner iteration
    ninner = 2

    ! lambdas for u and w subproblems
    lambda0 = 2.0*mu
    lambda1 = 2.0*mu*(alpha1/alpha0)

    n1 = size(f0, 1)
    n2 = size(f0, 2)
    n3 = size(f0, 3)
    call domain_decomp_regular(n1, n2, n3, n1beg, n1end, n2beg, n2end, n3beg, n3end)

    if (rankid <= rank1*rank2*rank3 - 1) then
        nonempty = .true.
    else
        nonempty = .false.
    end if

    ! memory
    call alloc_array(d1, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(d2, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(d3, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(d1t, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(d2t, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(d3t, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s11, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s12, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s13, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s22, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s23, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s33, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s11t, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s12t, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s13t, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s22t, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s23t, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s33t, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(w1, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(w2, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(w3, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(u, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(f, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(tmpd1w2, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(tmpd1w3, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(tmpd2w1, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(tmpd2w3, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(tmpd3w1, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(tmpd3w2, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(pu, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)

    ! noisy image
    f(n1beg:n1end, n2beg:n2end, n3beg:n3end) = f0(n1beg:n1end, n2beg:n2end, n3beg:n3end)
    u(n1beg:n1end, n2beg:n2end, n3beg:n3end) = f0(n1beg:n1end, n2beg:n2end, n3beg:n3end)

    ! TGpV iteration
    do iter = 1, niter

        pu = u

        do inner = 1, ninner

            if (nonempty) then
                call commute_array(u, 1)
                call commute_array(d1, 1)
                call commute_array(d2, 1)
                call commute_array(d3, 1)
                call commute_array(d1t, 1)
                call commute_array(d2t, 1)
                call commute_array(d3t, 1)
                call commute_array(w1, 1)
                call commute_array(w2, 1)
                call commute_array(w3, 1)
            end if

            ! minimization u
            do k = n3beg, n3end
                do j = n2beg, n2end
                    do i = n1beg, n1end

                        if (i == 1) then
                            i1 = 0
                            ii1 = i
                        else
                            i1 = 1
                            ii1 = i - 1
                        end if

                        if (i == n1) then
                            i2 = 0
                            ii2 = i
                        else
                            i2 = 1
                            ii2 = i + 1
                        end if

                        if (j == 1) then
                            j1 = 0
                            jj1 = j
                        else
                            j1 = 1
                            jj1 = j - 1
                        end if
                        if (j == n2) then
                            j2 = 0
                            jj2 = j
                        else
                            j2 = 1
                            jj2 = j + 1
                        end if

                        if (k == 1) then
                            k1 = 0
                            kk1 = k
                        else
                            k1 = 1
                            kk1 = k - 1
                        end if
                        if (k == n3) then
                            k2 = 0
                            kk2 = k
                        else
                            k2 = 1
                            kk2 = k + 1
                        end if

                        sumu = lambda0*( &
                            +i2*u(ii2, j, k) + i1*u(ii1, j, k) &
                            + j2*u(i, jj2, k) + j1*u(i, jj1, k) &
                            + k2*u(i, j, kk2) + k1*u(i, j, kk1)) &
                            + lambda0*( &
                            +(i1*d1(ii1, j, k) - i2*d1(i, j, k)) &
                            + (j1*d2(i, jj1, k) - j2*d2(i, j, k)) &
                            + (k1*d3(i, j, kk1) - k2*d3(i, j, k)) &
                            - (i1*d1t(ii1, j, k) - i2*d1t(i, j, k)) &
                            - (j1*d2t(i, jj1, k) - j2*d2t(i, j, k)) &
                            - (k1*d3t(i, j, kk1) - k2*d3t(i, j, k)) &
                            + (i1*w1(ii1, j, k) - i2*w1(i, j, k)) &
                            + (j1*w2(i, jj1, k) - j2*w2(i, j, k)) &
                            + (k1*w3(i, j, kk1) - k2*w3(i, j, k)) &
                            ) &
                            + mu*f(i, j, k)

                        u(i, j, k) = sumu/(mu + (6.0 - (1 - i1) - (1 - i2) - (1 - j1) - (1 - j2) - (1 - k1) - (1 - k2))*lambda0)

                    end do
                end do
            end do

            if (alpha1 /= 0) then

                ! minimization w
                do k = n3beg, n3end
                    do j = n2beg, n2end
                        do i = n1beg, n1end

                            if (i == n1) then
                                i1 = i - 1
                                i2 = i
                            else
                                i1 = i
                                i2 = i + 1
                            end if

                            if (j == n2) then
                                j1 = j - 1
                                j2 = j
                            else
                                j1 = j
                                j2 = j + 1
                            end if

                            if (k == n3) then
                                k1 = k - 1
                                k2 = k
                            else
                                k1 = k
                                k2 = k + 1
                            end if

                            tmpd1w2(i, j, k) = 0.5*(w2(i2, j, k) - w2(i1, j, k))
                            tmpd1w3(i, j, k) = 0.5*(w3(i2, j, k) - w3(i1, j, k))
                            tmpd2w1(i, j, k) = 0.5*(w1(i, j2, k) - w1(i, j1, k))
                            tmpd2w3(i, j, k) = 0.5*(w3(i, j2, k) - w3(i, j1, k))
                            tmpd3w1(i, j, k) = 0.5*(w1(i, j, k2) - w1(i, j, k1))
                            tmpd3w2(i, j, k) = 0.5*(w2(i, j, k2) - w2(i, j, k1))

                        end do
                    end do
                end do

                if (nonempty) then
                    call commute_array(u, 1)
                    call commute_array(tmpd1w2, 1)
                    call commute_array(tmpd1w3, 1)
                    call commute_array(tmpd2w1, 1)
                    call commute_array(tmpd2w3, 1)
                    call commute_array(tmpd3w1, 1)
                    call commute_array(tmpd3w2, 1)
                    call commute_array(s11, 1)
                    call commute_array(s12, 1)
                    call commute_array(s13, 1)
                    call commute_array(s22, 1)
                    call commute_array(s23, 1)
                    call commute_array(s33, 1)
                    call commute_array(s11t, 1)
                    call commute_array(s12t, 1)
                    call commute_array(s13t, 1)
                    call commute_array(s22t, 1)
                    call commute_array(s23t, 1)
                    call commute_array(s33t, 1)
                end if

                ! minimization of w
                do k = n3beg, n3end
                    do j = n2beg, n2end
                        do i = n1beg, n1end

                            if (i == 1) then
                                i1 = 0
                                ii1 = i
                            else
                                i1 = 1
                                ii1 = i - 1
                            end if

                            if (i == n1) then
                                i2 = 0
                                ii2 = i
                            else
                                i2 = 1
                                ii2 = i + 1
                            end if

                            if (j == 1) then
                                j1 = 0
                                jj1 = j
                            else
                                j1 = 1
                                jj1 = j - 1
                            end if
                            if (j == n2) then
                                j2 = 0
                                jj2 = j
                            else
                                j2 = 1
                                jj2 = j + 1
                            end if

                            if (k == 1) then
                                k1 = 0
                                kk1 = k
                            else
                                k1 = 1
                                kk1 = k - 1
                            end if
                            if (k == n3) then
                                k2 = 0
                                kk2 = k
                            else
                                k2 = 1
                                kk2 = k + 1
                            end if

                            ! w1
                            sumw1 = &
                                +lambda1*(i2*w1(ii2, j, k) + i1*w1(ii1, j, k)) &
                                + 0.5*lambda1*(j2*w1(i, jj2, k) + j1*w1(i, jj1, k)) &
                                + 0.5*lambda1*(k2*w1(i, j, kk2) + k1*w1(i, j, kk1)) &
                                - lambda0*(d1(i, j, k) - d1t(i, j, k) - (i2*u(ii2, j, k) - u(i, j, k))) &
                                + lambda1*(i1*s11(ii1, j, k) - i2*s11(i, j, k) - (i1*s11t(ii1, j, k) - i2*s11t(i, j, k))) &
                                + lambda1*(j1*s12(i, jj1, k) - j2*s12(i, j, k) - (j1*s12t(i, jj1, k) - j2*s12t(i, j, k)) - (j1*tmpd1w2(i, jj1, k) - j2*tmpd1w2(i, j, k))) &
                                + lambda1*(k1*s13(i, j, kk1) - k2*s13(i, j, k) - (k1*s13t(i, j, kk1) - k2*s13t(i, j, k)) - (k1*tmpd1w3(i, j, kk1) - k2*tmpd1w3(i, j, k)))

                            w1(i, j, k) = sumw1/(lambda0 + (2.0 - (1-i1) - (1-i2))*lambda1 + (2.0 - (1-j1) - (1-j2))*0.5*lambda1 + (2.0 - (1-k1) - (1-k2))*0.5*lambda1)

                            ! w2
                            sumw2 = &
                                +0.5*lambda1*(i2*w2(ii2, j, k) + i1*w2(ii1, j, k)) &
                                + lambda1*(j2*w2(i, jj2, k) + j1*w2(i, jj1, k)) &
                                + 0.5*lambda1*(k2*w2(i, j, kk2) + k1*w2(i, j, kk1)) &
                                - lambda0*(d2(i, j, k) - d2t(i, j, k) - (j2*u(i, jj2, k) - u(i, j, k))) &
                                + lambda1*(i1*s12(ii1, j, k) - i2*s12(i, j, k) - (i1*s12t(ii1, j, k) - i2*s12t(i, j, k)) - (i1*tmpd2w1(ii1, j, k) - i2*tmpd2w1(i, j, k))) &
                                + lambda1*(j1*s22(i, jj1, k) - j2*s22(i, j, k) - (j1*s22t(i, jj1, k) - j2*s22t(i, j, k))) &
                                + lambda1*(k1*s23(i, j, kk1) - k2*s23(i, j, k) - (k1*s23t(i, j, kk1) - k2*s23t(i, j, k)) - (k1*tmpd2w3(i, j, kk1) - k2*tmpd2w3(i, j, k)))

                            w2(i, j, k) = sumw2/(lambda0 + (2.0 - (1-i1) - (1-i2))*0.5*lambda1 + (2.0 - (1-j1) - (1-j2))*lambda1 + (2.0 - (1-k1) - (1-k2))*0.5*lambda1)

                            ! w3
                            sumw3 = &
                                +0.5*lambda1*(i2*w3(ii2, j, k) + i1*w3(ii1, j, k)) &
                                + 0.5*lambda1*(j2*w3(i, jj2, k) + j1*w3(i, jj1, k)) &
                                + lambda1*(k2*w3(i, j, kk2) + k1*w3(i, j, kk1)) &
                                - lambda0*(d3(i, j, k) - d3t(i, j, k) - (k2*u(i, j, kk2) - u(i, j, k))) &
                                + lambda1*(i1*s13(ii1, j, k) - i2*s13(i, j, k) - (i1*s13t(ii1, j, k) - i2*s13t(i, j, k)) - (i1*tmpd3w1(ii1, j, k) - i2*tmpd3w1(i, j, k))) &
                                + lambda1*(j1*s23(i, jj1, k) - j2*s23(i, j, k) - (j1*s23t(i, jj1, k) - j2*s23t(i, j, k)) - (j1*tmpd3w2(i, jj1, k) - j2*tmpd3w2(i, j, k))) &
                                + lambda1*(k1*s33(i, j, kk1) - k2*s33(i, j, k) - (k1*s33t(i, j, kk1) - k2*s33t(i, j, k)))

                            w3(i, j, k) = sumw3/(lambda0 + (2.0 - (1-i1) - (1-i2))*0.5*lambda1 + (2.0 - (1-j1) - (1-j2))*0.5*lambda1 + (2.0 - (1-k1) - (1-k2))*lambda1)

                        end do
                    end do
                end do

            end if

            ! update multipliers
            ! minimization d
            do k = n3beg, n3end
                do j = n2beg, n2end
                    do i = n1beg, n1end

                        if (i == n1) then
                            i1 = i - 1
                            i2 = i
                        else
                            i1 = i
                            i2 = i + 1
                        end if
                        tmp = (u(i2, j, k) - u(i1, j, k)) - w1(i, j, k) + d1t(i, j, k)
                        d1(i, j, k) = max(1.0 - (lambda0/alpha0*abs(tmp))**(p - 2), 0.0)*tmp

                        if (j == n2) then
                            j1 = j - 1
                            j2 = j
                        else
                            j1 = j
                            j2 = j + 1
                        end if
                        tmp = (u(i, j2, k) - u(i, j1, k)) - w2(i, j, k) + d2t(i, j, k)
                        d2(i, j, k) = max(1.0 - (lambda0/alpha0*abs(tmp))**(p - 2), 0.0)*tmp

                        if (k == n3) then
                            k1 = k - 1
                            k2 = k
                        else
                            k1 = k
                            k2 = k + 1
                        end if
                        tmp = (u(i, j, k2) - u(i, j, k1)) - w3(i, j, k) + d3t(i, j, k)
                        d3(i, j, k) = max(1.0 - (lambda0/alpha0*abs(tmp))**(p - 2), 0.0)*tmp

                    end do
                end do
            end do

            if (alpha1 /= 0) then

                if (nonempty) then
                    call commute_array(w1, 1)
                    call commute_array(w2, 1)
                    call commute_array(w3, 1)
                end if

                do k = n3beg, n3end
                    do j = n2beg, n2end
                        do i = n1beg, n1end

                            if (i == n1) then
                                i1 = i - 1
                                i2 = i
                            else
                                i1 = i
                                i2 = i + 1
                            end if
                            if (j == n2) then
                                j1 = j - 1
                                j2 = j
                            else
                                j1 = j
                                j2 = j + 1
                            end if
                            if (k == n3) then
                                k1 = k - 1
                                k2 = k
                            else
                                k1 = k
                                k2 = k + 1
                            end if

                            tmp = w1(i2, j, k) - w1(i1, j, k) + s11t(i, j, k)
                            s11(i, j, k) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp
                            tmp = w2(i, j2, k) - w2(i, j1, k) + s22t(i, j, k)
                            s22(i, j, k) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp
                            tmp = w3(i, j, k2) - w3(i, j, k1) + s33t(i, j, k)
                            s33(i, j, k) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp
                            tmp = 0.5*(w2(i2, j, k) - w2(i1, j, k) + w1(i, j2, k) - w1(i, j1, k)) + s12t(i, j, k)
                            s12(i, j, k) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp
                            tmp = 0.5*(w3(i2, j, k) - w3(i1, j, k) + w1(i, j, k2) - w1(i, j, k1)) + s13t(i, j, k)
                            s13(i, j, k) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp
                            tmp = 0.5*(w3(i, j2, k) - w3(i, j1, k) + w2(i, j, k2) - w2(i, j, k1)) + s23t(i, j, k)
                            s23(i, j, k) = max(1.0 - (lambda1/alpha1*abs(tmp))**(p - 2), 0.0)*tmp

                        end do
                    end do
                end do

            end if

            call mpi_barrier(mpi_comm_world, mpi_ierr)

        end do

        ! d
        do k = n3beg, n3end
            do j = n2beg, n2end
                do i = n1beg, n1end

                    if (i == n1) then
                        i1 = i - 1
                        i2 = i
                    else
                        i1 = i
                        i2 = i + 1
                    end if
                    d1t(i, j, k) = d1t(i, j, k) + (u(i2, j, k) - u(i1, j, k) - w1(i, j, k) - d1(i, j, k))

                    if (j == n2) then
                        j1 = j - 1
                        j2 = j
                    else
                        j1 = j
                        j2 = j + 1
                    end if
                    d2t(i, j, k) = d2t(i, j, k) + (u(i, j2, k) - u(i, j1, k) - w2(i, j, k) - d2(i, j, k))

                    if (k == n3) then
                        k1 = k - 1
                        k2 = k
                    else
                        k1 = k
                        k2 = k + 1
                    end if
                    d3t(i, j, k) = d3t(i, j, k) + (u(i, j, k2) - u(i, j, k1) - w3(i, j, k) - d3(i, j, k))

                end do
            end do
        end do

        if (alpha1 /= 0) then

            do k = n3beg, n3end
                do j = n2beg, n2end
                    do i = n1beg, n1end

                        if (i == n1) then
                            i1 = i - 1
                            i2 = i
                        else
                            i1 = i
                            i2 = i + 1
                        end if
                        if (j == n2) then
                            j1 = j - 1
                            j2 = j
                        else
                            j1 = j
                            j2 = j + 1
                        end if
                        if (k == n3) then
                            k1 = k - 1
                            k2 = k
                        else
                            k1 = k
                            k2 = k + 1
                        end if

                        s11t(i, j, k) = s11t(i, j, k) + (w1(i2, j, k) - w1(i1, j, k)) - s11(i, j, k)
                        s22t(i, j, k) = s22t(i, j, k) + (w2(i, j2, k) - w2(i, j1, k)) - s22(i, j, k)
                        s33t(i, j, k) = s33t(i, j, k) + (w3(i, j, k2) - w3(i, j, k1)) - s33(i, j, k)
                        s12t(i, j, k) = s12t(i, j, k) + 0.5*(w1(i, j2, k) - w1(i, j1, k) + w2(i2, j, k) - w2(i1, j, k)) - s12(i, j, k)
                        s13t(i, j, k) = s13t(i, j, k) + 0.5*(w1(i, j, k2) - w1(i, j, k1) + w3(i2, j, k) - w3(i1, j, k)) - s13(i, j, k)
                        s23t(i, j, k) = s23t(i, j, k) + 0.5*(w2(i, j, k2) - w2(i, j, k1) + w3(i, j2, k) - w3(i, j1, k)) - s23(i, j, k)

                    end do
                end do
            end do

        end if

        call mpi_barrier(mpi_comm_world, mpi_ierr)

        ! progress
        if (mod(iter, max(nint(niter/10.0), 1)) == 0 .or. iter == 1) then
            if (nonempty) then
                bsum = sum((u - pu)**2)
            else
                bsum = 0.0
            end if
            call allreduce(bsum)
            if (rankid == 0) then
                call warn(' >> TGpV iteration '//num2str(iter, '(i)') &
                    //' of '//num2str(niter, '(i)')//' relative l2-norm difference = '//num2str(sqrt(bsum), '(es)'))
            end if
        end if

    end do

    call mpi_barrier(mpi_comm_world, mpi_ierr)

    call alloc_array(u0, [1, n1, 1, n2, 1, n3])
    if (nonempty) then
        u0(n1beg:n1end, n2beg:n2end, n3beg:n3end) = u(n1beg:n1end, n2beg:n2end, n3beg:n3end)
    end if
    call allreduce_array(u0)

    call mpi_barrier(mpi_comm_world, mpi_ierr)

end function tgpv_filt_3d_mpi_

!
!> 2D TV + sparsity denoising
!
function sparse_tv_filt_1d_(f, mu, l_sparse, l_tv, niter, verbose) result(u)

    TT, dimension(:) :: f
    TT :: mu, l_sparse, l_tv
    integer :: niter
    logical, optional :: verbose
    TT, allocatable, dimension(:) :: u

    integer :: i, iter
    TT :: tmp, tmpx, sk, sumu
    TT, allocatable, dimension(:) ::  d1, d1t, s, st, pu
    integer :: n1
    TT :: alpha0, alpha1, p
    integer :: inner, ninner
    integer :: i1, i2
    integer :: ii1, ii2
    TT :: lambda0, lambda1
    logical :: tv_verbose

    ! read parameter
    n1 = size(f)
    alpha0 = l_tv
    alpha1 = l_sparse
    call assert(alpha0 <= 1, ' <sparse_tv_filt_2d_> Error: l_tv must >= 0 and <= 1')
    call assert(alpha1 <= 1, ' <sparse_tv_filt_2d_> Error: l_sparse must >= 0 and <= 1')
    p = 1.0
    ninner = 2

    if (present(verbose)) then
        tv_verbose = verbose
    else
        tv_verbose = .true.
    end if

    ! memory
    d1 = zeros(n1)
    d1t = zeros(n1)
    s = zeros(n1)
    st = zeros(n1)
    u = zeros(n1)
    pu = zeros(n1)

    ! read noisy image
    u = f

    ! lambdas for u and w subproblems
    lambda0 = 2*mu*alpha0
    lambda1 = 2*mu*alpha1

    ! TGpV iteration
    do iter = 1, niter

        pu = u

        do inner = 1, ninner

            ! minimization u
            do i = 1, n1

                if (i == 1) then
                    i1 = 0
                    ii1 = i
                else
                    i1 = 1
                    ii1 = i - 1
                end if

                if (i == n1) then
                    i2 = 0
                    ii2 = i
                else
                    i2 = 1
                    ii2 = i + 1
                end if

                sumu = lambda0*( &
                    +i2*u(ii2) + i1*u(ii1)) &
                    + lambda1*( &
                    +s(i) - st(i)) &
                    + lambda0*( &
                    +(i1*d1(ii1) - i2*d1(i)) &
                    - (i1*d1t(ii1) - i2*d1t(i)) &
                    ) &
                    + mu*f(i)

                u(i) = sumu/(mu + lambda1 + (2.0 - (1 - i1) - (1 - i2))*lambda0)

            end do

            ! update multipliers
            !$omp parallel do private(i, i1, i2, tmpx, sk)
            do i = 1, n1

                if (i == n1) then
                    i1 = i - 1
                    i2 = i
                else
                    i1 = i
                    i2 = i + 1
                end if
                tmpx = (u(i2) - u(i1)) + d1t(i)

                sk = abs(tmpx)

                d1(i) = max(sk - 1.0/(lambda0/alpha0), 0.0)*tmpx/(sk + float_tiny)

            end do
            !$omp end parallel do

            !$omp parallel do private(i, tmp, sk)
            do i = 1, n1

                tmp = u(i) + st(i)
                sk = abs(tmp)
                s(i) = max(sk - 1.0/(lambda1/alpha1), 0.0)*tmp/(sk + float_tiny)

            end do
            !$omp end parallel do

        end do

        ! update split-Bregman variables
        !$omp parallel do private(i, i1, i2)
        do i = 1, n1

            if (i == n1) then
                i1 = i - 1
                i2 = i
            else
                i1 = i
                i2 = i + 1
            end if
            d1t(i) = d1t(i) + (u(i2) - u(i1) - d1(i))

        end do
        !$omp end parallel do

        !$omp parallel do private(i)
        do i = 1, n1

            st(i) = st(i) + u(i) - s(i)

        end do
        !$omp end parallel do

        ! progress
        if (tv_verbose .and. (mod(iter, max(nint(niter/10.0), 1)) == 0 .or. iter == 1)) then
            call warn(' >> Sparse TV iteration '//tidy(num2str(iter, '(i)')) &
                //' of '//tidy(num2str(niter, '(i)')//' relative l2-norm difference = ' &
                //tidy(num2str(norm2(u - pu), '(es)'))))
        end if

    end do

end function sparse_tv_filt_1d_

!
!> 2D TV + sparsity denoising
!
function sparse_tv_filt_2d_(f, mu, l_sparse, l_tv, niter, verbose) result(u)

    TT, dimension(:, :) :: f
    TT :: mu, l_sparse, l_tv
    integer :: niter
    logical, optional :: verbose
    TT, allocatable, dimension(:, :) :: u

    integer :: i, j, iter
    TT :: tmp, tmpx, tmpz, sk, sumu
    TT, allocatable, dimension(:, :) ::  d1, d2, d1t, d2t, s, st, pu
    integer :: n1, n2
    TT :: alpha0, alpha1, p
    integer :: inner, ninner
    integer :: i1, i2, j1, j2
    integer :: ii1, ii2, jj1, jj2
    TT :: lambda0, lambda1
    logical :: tv_verbose

    ! read parameter
    n1 = size(f, 1)
    n2 = size(f, 2)
    alpha0 = l_tv
    alpha1 = l_sparse
    call assert(alpha0 <= 1, ' <sparse_tv_filt_2d_> Error: l_tv must >= 0 and <= 1')
    call assert(alpha1 <= 1, ' <sparse_tv_filt_2d_> Error: l_sparse must >= 0 and <= 1')
    p = 1.0
    ninner = 2

    if (present(verbose)) then
        tv_verbose = verbose
    else
        tv_verbose = .true.
    end if

    ! memory
    d1 = zeros(n1, n2)
    d2 = zeros(n1, n2)
    d1t = zeros(n1, n2)
    d2t = zeros(n1, n2)
    s = zeros(n1, n2)
    st = zeros(n1, n2)
    u = zeros(n1, n2)
    pu = zeros(n1, n2)

    ! read noisy image
    u = f

    ! lambdas for u and w subproblems
    lambda0 = 2*mu*alpha0
    lambda1 = 2*mu*alpha1

    ! TGpV iteration
    do iter = 1, niter

        pu = u

        do inner = 1, ninner

            ! minimization u
            do j = 1, n2
                do i = 1, n1

                    if (i == 1) then
                        i1 = 0
                        ii1 = i
                    else
                        i1 = 1
                        ii1 = i - 1
                    end if

                    if (i == n1) then
                        i2 = 0
                        ii2 = i
                    else
                        i2 = 1
                        ii2 = i + 1
                    end if

                    if (j == 1) then
                        j1 = 0
                        jj1 = j
                    else
                        j1 = 1
                        jj1 = j - 1
                    end if
                    if (j == n2) then
                        j2 = 0
                        jj2 = j
                    else
                        j2 = 1
                        jj2 = j + 1
                    end if

                    sumu = lambda0*( &
                        +i2*u(ii2, j) + i1*u(ii1, j) &
                        + j2*u(i, jj2) + j1*u(i, jj1)) &
                        + lambda1*( &
                        +s(i, j) - st(i, j)) &
                        + lambda0*( &
                        +(i1*d1(ii1, j) - i2*d1(i, j)) &
                        + (j1*d2(i, jj1) - j2*d2(i, j)) &
                        - (i1*d1t(ii1, j) - i2*d1t(i, j)) &
                        - (j1*d2t(i, jj1) - j2*d2t(i, j)) &
                        ) &
                        + mu*f(i, j)

                    u(i, j) = sumu/(mu + lambda1 + (4.0 - (1 - i1) - (1 - i2) - (1 - j1) - (1 - j2))*lambda0)

                end do
            end do

            ! update multipliers
            !$omp parallel do private(i, j, i1, i2, j1, j2, tmpx, tmpz, sk)
            do j = 1, n2
                do i = 1, n1

                    if (i == n1) then
                        i1 = i - 1
                        i2 = i
                    else
                        i1 = i
                        i2 = i + 1
                    end if
                    tmpx = (u(i2, j) - u(i1, j)) + d1t(i, j)

                    if (j == n2) then
                        j1 = j - 1
                        j2 = j
                    else
                        j1 = j
                        j2 = j + 1
                    end if
                    tmpz = (u(i, j2) - u(i, j1)) + d2t(i, j)

                    sk = sqrt(tmpx**2 + tmpz**2)

                    d1(i, j) = max(sk - 1.0/(lambda0/alpha0), 0.0)*tmpx/(sk + float_tiny)
                    d2(i, j) = max(sk - 1.0/(lambda0/alpha0), 0.0)*tmpz/(sk + float_tiny)

                end do
            end do
            !$omp end parallel do

            !$omp parallel do private(i, j, tmp, sk)
            do j = 1, n2
                do i = 1, n1

                    tmp = u(i, j) + st(i, j)
                    sk = abs(tmp)
                    s(i, j) = max(sk - 1.0/(lambda1/alpha1), 0.0)*tmp/(sk + float_tiny)

                end do
            end do
            !$omp end parallel do

        end do

        ! update split-Bregman variables
        !$omp parallel do private(i, j, i1, i2, j1, j2)
        do j = 1, n2
            do i = 1, n1

                if (i == n1) then
                    i1 = i - 1
                    i2 = i
                else
                    i1 = i
                    i2 = i + 1
                end if
                d1t(i, j) = d1t(i, j) + (u(i2, j) - u(i1, j) - d1(i, j))

                if (j == n2) then
                    j1 = j - 1
                    j2 = j
                else
                    j1 = j
                    j2 = j + 1
                end if
                d2t(i, j) = d2t(i, j) + (u(i, j2) - u(i, j1) - d2(i, j))

            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j)
        do j = 1, n2
            do i = 1, n1

                st(i, j) = st(i, j) + u(i, j) - s(i, j)

            end do
        end do
        !$omp end parallel do

        ! progress
        if (tv_verbose .and. (mod(iter, max(nint(niter/10.0), 1)) == 0 .or. iter == 1)) then
            call warn(' >> Sparse TV iteration '//tidy(num2str(iter, '(i)')) &
                //' of '//tidy(num2str(niter, '(i)')//' relative l2-norm difference = ' &
                //tidy(num2str(norm2(u - pu), '(es)'))))
        end if

    end do

end function sparse_tv_filt_2d_

!
!> 3D TV + sparsity denoising
!
function sparse_tv_filt_3d_(f, mu, l_sparse, l_tv, niter, verbose) result(u)

    TT, dimension(:, :, :) :: f
    TT :: mu, l_sparse, l_tv
    integer :: niter
    logical, optional :: verbose
    TT, allocatable, dimension(:, :, :) :: u

    integer :: i, j, k, iter
    TT :: tmp, tmpx, tmpy, tmpz, sk, sumu
    TT, allocatable, dimension(:, :, :) ::  d1, d2, d3, d1t, d2t, d3t, s, st, pu
    integer :: n1, n2, n3
    TT :: alpha0, alpha1, p
    integer :: inner, ninner
    integer :: i1, i2, j1, j2, k1, k2
    integer :: ii1, ii2, jj1, jj2, kk1, kk2
    TT :: lambda0, lambda1
    logical :: tv_verbose

    ! read parameter
    n1 = size(f, 1)
    n2 = size(f, 2)
    n3 = size(f, 3)
    alpha0 = l_tv
    alpha1 = l_sparse
    call assert(alpha0 <= 1, ' <sparse_tv_filt_3d_> Error: l_tv must >= 0 and <= 1')
    call assert(alpha1 <= 1, ' <sparse_tv_filt_3d_> Error: l_sparse must >= 0 and <= 1')
    p = 1.0
    ninner = 2

    if (present(verbose)) then
        tv_verbose = verbose
    else
        tv_verbose = .true.
    end if

    ! memory
    d1 = zeros(n1, n2, n3)
    d2 = zeros(n1, n2, n3)
    d3 = zeros(n1, n2, n3)
    d1t = zeros(n1, n2, n3)
    d2t = zeros(n1, n2, n3)
    d3t = zeros(n1, n2, n3)
    s = zeros(n1, n2, n3)
    st = zeros(n1, n2, n3)
    u = zeros(n1, n2, n3)
    pu = zeros(n1, n2, n3)

    ! read noisy image
    u = f

    ! lambdas for u and w subproblems
    lambda0 = 2*mu*alpha0
    lambda1 = 2*mu*alpha1

    ! TGpV iteration
    do iter = 1, niter

        pu = u

        do inner = 1, ninner

            ! minimization u
            do k = 1, n3
                do j = 1, n2
                    do i = 1, n1

                        if (i == 1) then
                            i1 = 0
                            ii1 = i
                        else
                            i1 = 1
                            ii1 = i - 1
                        end if

                        if (i == n1) then
                            i2 = 0
                            ii2 = i
                        else
                            i2 = 1
                            ii2 = i + 1
                        end if

                        if (j == 1) then
                            j1 = 0
                            jj1 = j
                        else
                            j1 = 1
                            jj1 = j - 1
                        end if
                        if (j == n2) then
                            j2 = 0
                            jj2 = j
                        else
                            j2 = 1
                            jj2 = j + 1
                        end if

                        if (k == 1) then
                            k1 = 0
                            kk1 = k
                        else
                            k1 = 1
                            kk1 = k - 1
                        end if
                        if (k == n3) then
                            k2 = 0
                            kk2 = k
                        else
                            k2 = 1
                            kk2 = k + 1
                        end if

                        sumu = lambda0*( &
                            +i2*u(ii2, j, k) + i1*u(ii1, j, k) &
                            + j2*u(i, jj2, k) + j1*u(i, jj1, k) &
                            + k2*u(i, j, kk2) + k1*u(i, j, kk1)) &
                            + lambda1*( &
                            +s(i, j, k) - st(i, j, k)) &
                            + lambda0*( &
                            +(i1*d1(ii1, j, k) - i2*d1(i, j, k)) &
                            + (j1*d2(i, jj1, k) - j2*d2(i, j, k)) &
                            + (k1*d3(i, j, kk1) - k2*d3(i, j, k)) &
                            - (i1*d1t(ii1, j, k) - i2*d1t(i, j, k)) &
                            - (j1*d2t(i, jj1, k) - j2*d2t(i, j, k)) &
                            - (k1*d3t(i, j, kk1) - k2*d3t(i, j, k))) &
                            + mu*f(i, j, k)

                        u(i, j, k) = sumu/(mu + lambda1 + (6.0 - &
                            (1 - i1) - (1 - i2) - (1 - j1) - (1 - j2) - (1 - k1) - (1 - k2))*lambda0)

                    end do
                end do
            end do

            ! update multipliers
            !$omp parallel do private(i, j, k, i1, i2, j1, j2, k1, k2, tmpx, tmpy, tmpz, sk)
            do k = 1, n3
                do j = 1, n2
                    do i = 1, n1

                        if (i == n1) then
                            i1 = i - 1
                            i2 = i
                        else
                            i1 = i
                            i2 = i + 1
                        end if
                        tmpx = (u(i2, j, k) - u(i1, j, k)) + d1t(i, j, k)

                        if (j == n2) then
                            j1 = j - 1
                            j2 = j
                        else
                            j1 = j
                            j2 = j + 1
                        end if
                        tmpy = (u(i, j2, k) - u(i, j1, k)) + d2t(i, j, k)

                        if (k == n3) then
                            k1 = k - 1
                            k2 = k
                        else
                            k1 = k
                            k2 = k + 1
                        end if
                        tmpz = (u(i, j, k2) - u(i, j, k1)) + d3t(i, j, k)

                        sk = sqrt(tmpx**2 + tmpy**2 + tmpz**2)

                        d1(i, j, k) = max(sk - 1.0/(lambda0/alpha0), 0.0)*tmpx/(sk + float_tiny)
                        d2(i, j, k) = max(sk - 1.0/(lambda0/alpha0), 0.0)*tmpy/(sk + float_tiny)
                        d3(i, j, k) = max(sk - 1.0/(lambda0/alpha0), 0.0)*tmpz/(sk + float_tiny)

                    end do
                end do
            end do
            !$omp end parallel do

            !$omp parallel do private(i, j, k, tmp, sk)
            do k = 1, n3
                do j = 1, n2
                    do i = 1, n1

                        tmp = u(i, j, k) + st(i, j, k)
                        sk = abs(tmp)
                        s(i, j, k) = max(sk - 1.0/(lambda1/alpha1), 0.0)*tmp/(sk + float_tiny)

                    end do
                end do
            end do
            !$omp end parallel do

        end do

        ! update split-Bregman variables
        !$omp parallel do private(i, j, k, i1, i2, j1, j2, k1, k2)
        do k = 1, n3
            do j = 1, n2
                do i = 1, n1

                    if (i == n1) then
                        i1 = i - 1
                        i2 = i
                    else
                        i1 = i
                        i2 = i + 1
                    end if
                    d1t(i, j, k) = d1t(i, j, k) + (u(i2, j, k) - u(i1, j, k) - d1(i, j, k))

                    if (j == n2) then
                        j1 = j - 1
                        j2 = j
                    else
                        j1 = j
                        j2 = j + 1
                    end if
                    d2t(i, j, k) = d2t(i, j, k) + (u(i, j2, k) - u(i, j1, k) - d2(i, j, k))

                    if (k == n3) then
                        k1 = k - 1
                        k2 = k
                    else
                        k1 = k
                        k2 = k + 1
                    end if
                    d3t(i, j, k) = d3t(i, j, k) + (u(i, j, k2) - u(i, j, k1) - d3(i, j, k))

                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k)
        do k = 1, n3
            do j = 1, n2
                do i = 1, n1

                    st(i, j, k) = st(i, j, k) + u(i, j, k) - s(i, j, k)

                end do
            end do
        end do
        !$omp end parallel do

        ! progress
        if (tv_verbose .and. (mod(iter, max(nint(niter/10.0), 1)) == 0 .or. iter == 1)) then
            call warn(' >> Sparse TV iteration '//tidy(num2str(iter, '(i)')) &
                //' of '//tidy(num2str(niter, '(i)')//' relative l2-norm difference = ' &
                //tidy(num2str(norm2(u - pu), '(es)'))))
        end if

    end do

end function sparse_tv_filt_3d_

!
!> 1D soft thresholding
!
function soft_shrinkage_1d_(x, b) result(y)

    TT, dimension(:), intent(in) :: x
    TT, intent(in) :: b
    TT, allocatable, dimension(:) :: y

    y = x*0
    where (abs(x) >= b)
        y = x*(1.0 - b/abs(x))
    end where

end function soft_shrinkage_1d_

!
!> 2D soft thresholding
!
function soft_shrinkage_2d_(x, b) result(y)

    TT, dimension(:, :), intent(in) :: x
    TT, intent(in) :: b
    TT, allocatable, dimension(:, :) :: y

    y = x*0
    where (abs(x) >= b)
        y = x*(1.0 - b/abs(x))
    end where

end function soft_shrinkage_2d_

!
!> 3D soft thresholding
!
function soft_shrinkage_3d_(x, b) result(y)

    TT, dimension(:, :, :), intent(in) :: x
    TT, intent(in) :: b
    TT, allocatable, dimension(:, :, :) :: y

    y = x*0
    where (abs(x) >= b)
        y = x*(1.0 - b/abs(x))
    end where

end function soft_shrinkage_3d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef tv_iso_filt_1d_
#undef tv_iso_filt_2d_
#undef tv_iso_filt_3d_
#undef sparse_tv_filt_1d_
#undef sparse_tv_filt_2d_
#undef sparse_tv_filt_3d_
#undef tgpv_filt_1d_
#undef tgpv_filt_2d_
#undef tgpv_filt_3d_
#undef tgpv_filt_2d_mpi_
#undef tgpv_filt_3d_mpi_
#undef soft_shrinkage_1d_
#undef soft_shrinkage_2d_
#undef soft_shrinkage_3d_
