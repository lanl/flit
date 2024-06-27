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

#define compute_deriche_coefs_      CONCAT(compute_deriche_coefs, T)
#define apply_deriche_filter_1d_      CONCAT(apply_deriche_filter_1d, T)
#define deriche_gaussian_filt_1d_      CONCAT(deriche_gaussian_filt_1d, T)
#define deriche_gaussian_filt_2d_      CONCAT(deriche_gaussian_filt_2d, T)
#define deriche_gaussian_filt_3d_      CONCAT(deriche_gaussian_filt_3d, T)
#define dct_gaussian_filt_1d_      CONCAT(dct_gaussian_filt_1d, T)
#define dct_gaussian_filt_2d_      CONCAT(dct_gaussian_filt_2d, T)
#define dct_gaussian_filt_3d_      CONCAT(dct_gaussian_filt_3d, T)
#define dft_gaussian_filt_1d_      CONCAT(dft_gaussian_filt_1d, T)
#define dft_gaussian_filt_2d_      CONCAT(dft_gaussian_filt_2d, T)
#define dft_gaussian_filt_3d_      CONCAT(dft_gaussian_filt_3d, T)
#define athead_      CONCAT(athead, T)
#define attail_      CONCAT(attail, T)
#define compute_rsgf_coef_      CONCAT(compute_rsgf_coef, T)
#define rsgf_1d_      CONCAT(rsgf_1d, T)
#define spectral_gaussian_filt_1d_      CONCAT(spectral_gaussian_filt_1d, T)
#define spectral_gaussian_filt_2d_      CONCAT(spectral_gaussian_filt_2d, T)
#define spectral_gaussian_filt_3d_      CONCAT(spectral_gaussian_filt_3d, T)
#define conv_gaussian_filt_1d_      CONCAT(conv_gaussian_filt_1d, T)
#define conv_gaussian_filt_2d_      CONCAT(conv_gaussian_filt_2d, T)
#define conv_gaussian_filt_3d_      CONCAT(conv_gaussian_filt_3d, T)
#define gauss_filt_1d_      CONCAT(gauss_filt_1d, T)
#define gauss_filt_2d_      CONCAT(gauss_filt_2d, T)
#define gauss_filt_3d_      CONCAT(gauss_filt_3d, T)

!==========================================================================
! ITK's Deriche implementation, which has slightly exaggerated boundaries
! in the output, but better outputs for constant, linear and quatratic
! signals (but who will try to smooth these?...)

!
!> Compute Deriche coefficients
!
subroutine compute_deriche_coefs_(sigma, order, coefd, coefn, coefm, coefbn, coefbm)

    TT, intent(in) :: sigma
    integer, intent(in) :: order
    TT, dimension(:), intent(inout) :: coefd, coefn, coefm, coefbn, coefbm

    TT :: sin1, sin2, cos1, cos2
    TT :: exp1, exp2
    TT, dimension(0:2) :: a1, a2
    TT, dimension(0:2) :: b1, b2
    TT :: w1, w2
    TT :: l1, l2
    TT, dimension(1:4) :: cd
    TT, dimension(0:2, 0:3) :: cn

    TT :: sn(0:2), sd
    TT :: dn(0:2), dd
    TT :: en(0:2), ed
    TT :: nc0a, nc1a, nc2a, nc2b

    if (sigma > 0) then

        a1 = [+1.3530, -0.6724, -1.3563]
        a2 = [-0.3531, +0.6724, +0.3446]
        b1 = [+1.8151, -3.4327, +5.2318]
        b2 = [+0.0902, +0.6100, -2.2355]
        w1 = +0.6681
        w2 = +2.0787
        l1 = -1.3932
        l2 = -1.3732

        sin1 = sin(w1/sigma)
        sin2 = sin(w2/sigma)
        cos1 = cos(w1/sigma)
        cos2 = cos(w2/sigma)
        exp1 = exp(l1/sigma)
        exp2 = exp(l2/sigma)

        ! compute the coefficients
        cd(1) = -2.0*exp2*cos2 - 2.0*exp1*cos1
        cd(2) = 4.0*exp1*exp2*cos1*cos2 + exp1**2 + exp2**2
        cd(3) = -2.0*exp1*exp2**2*cos1 - 2.0*exp1**2*exp2*cos2
        cd(4) = exp1**2*exp2**2

        cn(:, 0) = a1 + a2
        cn(:, 1) = exp2*(b2*sin2 - (a2 + 2.0*a1)*cos2) &
            + exp1*(b1*sin1 - (a1 + 2.0*a2)*cos1)
        cn(:, 2) = 2.0*exp1*exp2 &
            *((a1 + a2)*cos1*cos2 - b1*sin1*cos2 - b2*sin2*cos1) &
            + a2*exp1**2 + a1*exp2**2
        cn(:, 3) = exp1**2*exp2*(b2*sin2 - a2*cos2) &
            + exp1*exp2**2*(b1*sin1 - a1*cos1)

        ! normalize the coefficients
        sn(:) = sum(cn, dim=2)
        sd = 1.0 + sum(cd)
        dn(:) = cn(:, 1) + 2.0*cn(:, 2) + 3.0*cn(:, 3)
        dd = cd(1) + 2.0*cd(2) + 3.0*cd(3) + 4.0*cd(4)
        en(:) = cn(:, 1) + 4.0*cn(:, 2) + 9.0*cn(:, 3)
        ed = cd(1) + 4.0*cd(2) + 9.0*cd(3) + 16.0*cd(4)

        select case (order)
            case (0)
                nc0a = 2.0*sn(0)/sd - cn(0, 0)
                coefn = cn(0, :)/nc0a
            case (1)
                nc1a = 2.0*(sn(1)*dd - dn(1)*sd)/sd**2
                coefn = cn(1, :)/nc1a
            case (2)
                nc2b = -(2.0*sn(2) - sd*cn(2, 0))/(2.0*sn(0) - sd*cn(0, 0))
                nc2a = 1.0/sd**3 &
                    *((en(2) + nc2b*en(0))*sd**2 &
                    - ed*(sn(2) + nc2b*sn(0))*sd &
                    - 2.0*(dn(2) + nc2b*dn(2))*dd*sd &
                    + 2.0*dd**2*(sn(2) + nc2b*sn(2)))
                coefn = (cn(2, :) + nc2b*cn(0, :))/nc2a
        end select

        coefd = cd

        ! coefficient m for the anticausal part
        select case (order)
            case (0, 2)
                coefm(1:3) = coefn(2:4) - coefd(1:3)*coefn(1)
                coefm(4) = -coefd(4)*coefn(1)
            case (1)
                coefm(1:3) = -coefn(2:4) + coefd(1:3)*coefn(1)
                coefm(4) = coefd(4)*coefn(1)
        end select

        ! boundary coefficients
        coefbn = coefd*sum(coefn)/sd
        coefbm = coefd*sum(coefm)/sd

    end if

end subroutine compute_deriche_coefs_

!
!> Apply 1D Deriche filter
!
subroutine apply_deriche_filter_1d_(w, coefd, coefn, coefm, coefbn, coefbm)

    TT, dimension(:), intent(inout) :: w
    TT, dimension(:), intent(in) :: coefd, coefn, coefm, coefbn, coefbm

    integer :: i, n
    TT, allocatable, dimension(:) :: ww, wl, wr

    ! length of signal
    n = size(w)

    ! copy input array
    allocate (ww(1:n), source=w)
    allocate (wl(1:n))
    allocate (wr(1:n))
    wl = 0.0
    wr = 0.0

    ! causal part
    wl(1) = ww(1)*coefn(1) + ww(1)*coefn(2) + ww(1)*coefn(3) + ww(1)*coefn(4)
    wl(2) = ww(2)*coefn(1) + ww(1)*coefn(2) + ww(1)*coefn(3) + ww(1)*coefn(4)
    wl(3) = ww(3)*coefn(1) + ww(2)*coefn(2) + ww(1)*coefn(3) + ww(1)*coefn(4)
    wl(4) = ww(4)*coefn(1) + ww(3)*coefn(2) + ww(2)*coefn(3) + ww(1)*coefn(4)

    wl(1) = wl(1) - (ww(1)*coefbn(1) + ww(1)*coefbn(2) + ww(1)*coefbn(3) + ww(1)*coefbn(4))
    wl(2) = wl(2) - (wl(1)*coefd(1) + ww(1)*coefbn(2) + ww(1)*coefbn(3) + ww(1)*coefbn(4))
    wl(3) = wl(3) - (wl(2)*coefd(1) + wl(1)*coefd(2) + ww(1)*coefbn(3) + ww(1)*coefbn(4))
    wl(4) = wl(4) - (wl(3)*coefd(1) + wl(2)*coefd(2) + wl(1)*coefd(3) + ww(1)*coefbn(4))

    do i = 5, n
        ! the original expression, could be a little bit of slower than the next
        !       wl(i) = sum(coefn*ww(i:i-3:-1)) - sum(coefd*wl(i-1:i-4:-1))
        ! reversed version, could be a little faster
        wl(i) = sum(coefn(4:1:-1)*ww(i - 3:i)) - sum(coefd(4:1:-1)*wl(i - 4:i - 1))
    end do

    ! anti-causal part
    wr(n) = ww(n)*coefm(1) + ww(n)*coefm(2) + ww(n)*coefm(3) + ww(n)*coefm(4)
    wr(n - 1) = ww(n)*coefm(1) + ww(n)*coefm(2) + ww(n)*coefm(3) + ww(n)*coefm(4)
    wr(n - 2) = ww(n - 1)*coefm(1) + ww(n)*coefm(2) + ww(n)*coefm(3) + ww(n)*coefm(4)
    wr(n - 3) = ww(n - 2)*coefm(1) + ww(n - 1)*coefm(2) + ww(n)*coefm(3) + ww(n)*coefm(4)

    wr(n) = wr(n) - (ww(n)*coefbm(1) + ww(n)*coefbm(2) + ww(n)*coefbm(3) + ww(n)*coefbm(4))
    wr(n - 1) = wr(n - 1) - (wr(n)*coefd(1) + ww(n)*coefbm(2) + ww(n)*coefbm(3) + ww(n)*coefbm(4))
    wr(n - 2) = wr(n - 2) - (wr(n - 1)*coefd(1) + wr(n)*coefd(2) + ww(n)*coefbm(3) + ww(n)*coefbm(4))
    wr(n - 3) = wr(n - 3) - (wr(n - 2)*coefd(1) + wr(n - 1)*coefd(2) + wr(n)*coefd(3) + ww(n)*coefbm(4))

    do i = n - 4, 1, -1
        wr(i) = sum(coefm*ww(i + 1:i + 4)) - sum(coefd*wr(i + 1:i + 4))
    end do

    ! summation
    w = wl + wr

end subroutine apply_deriche_filter_1d_

!
!> 1D Gaussian filter based on ITK's Deriche implementation
!
subroutine deriche_gaussian_filt_1d_(w, sigma, order)

    TT, dimension(:), intent(inout) :: w
    TT, intent(in) :: sigma
    integer, intent(in) :: order

    TT, dimension(1:4) :: coefd, coefn, coefm, coefbn, coefbm

    ! compute Deriche coefficients
    call compute_deriche_coefs_(sigma, order, coefd, coefn, coefm, coefbn, coefbm)

    ! apply Deriche filtering
    call apply_deriche_filter_1d_(w, coefd, coefn, coefm, coefbn, coefbm)

end subroutine deriche_gaussian_filt_1d_

!
!> 2D Gaussian filter based on ITK's Deriche implementation
!
subroutine deriche_gaussian_filt_2d_(w, sigma, order)

    TT, dimension(:, :), intent(inout) :: w
    TT, dimension(:), intent(in) :: sigma
    integer, dimension(:), intent(in) :: order

    TT, dimension(1:4) :: coefd1, coefn1, coefm1, coefbn1, coefbm1
    TT, dimension(1:4) :: coefd2, coefn2, coefm2, coefbn2, coefbm2
    integer :: i, j, n1, n2

    ! compute Deriche coefficients
    call compute_deriche_coefs_(sigma(1), order(1), coefd1, coefn1, coefm1, coefbn1, coefbm1)
    call compute_deriche_coefs_(sigma(2), order(2), coefd2, coefn2, coefm2, coefbn2, coefbm2)

    n1 = size(w, 1)
    n2 = size(w, 2)

    ! apply Deriche filtering
    if (n2 > 1) then
        !$omp parallel do private(i)
        do i = 1, n1
            call apply_deriche_filter_1d_(w(i, :), coefd2, coefn2, coefm2, coefbn2, coefbm2)
        end do
        !$omp end parallel do
    end if

    if (n1 > 1) then
        !$omp parallel do private(j)
        do j = 1, n2
            call apply_deriche_filter_1d_(w(:, j), coefd1, coefn1, coefm1, coefbn1, coefbm1)
        end do
        !$omp end parallel do
    end if

end subroutine deriche_gaussian_filt_2d_

!
!> 3D Gaussian filter based on ITK's Deriche implementation
!
subroutine deriche_gaussian_filt_3d_(w, sigma, order)

    TT, dimension(:, :, :), intent(inout) :: w
    TT, dimension(:), intent(in) :: sigma
    integer, dimension(:), intent(in) :: order

    TT, dimension(1:4) :: coefd1, coefn1, coefm1, coefbn1, coefbm1
    TT, dimension(1:4) :: coefd2, coefn2, coefm2, coefbn2, coefbm2
    TT, dimension(1:4) :: coefd3, coefn3, coefm3, coefbn3, coefbm3
    integer :: i, j, k, n1, n2, n3

    ! compute Deriche coefficients
    call compute_deriche_coefs_(sigma(1), order(1), coefd1, coefn1, coefm1, coefbn1, coefbm1)
    call compute_deriche_coefs_(sigma(2), order(2), coefd2, coefn2, coefm2, coefbn2, coefbm2)
    call compute_deriche_coefs_(sigma(3), order(3), coefd3, coefn3, coefm3, coefbn3, coefbm3)

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)

    ! apply Deriche filtering
    if (n3 > 1) then
        !$omp parallel do private(i,j)
        do j = 1, n2
            do i = 1, n1
                call apply_deriche_filter_1d_(w(i, j, :), coefd3, coefn3, coefm3, coefbn3, coefbm3)
            end do
        end do
        !$omp end parallel do
    end if

    if (n2 > 1) then
        !$omp parallel do private(i,k)
        do k = 1, n3
            do i = 1, n1
                call apply_deriche_filter_1d_(w(i, :, k), coefd2, coefn2, coefm2, coefbn2, coefbm2)
            end do
        end do
        !$omp end parallel do
    end if

    if (n1 > 1) then
        !$omp parallel do private(j,k)
        do k = 1, n3
            do j = 1, n2
                call apply_deriche_filter_1d_(w(:, j, k), coefd1, coefn1, coefm1, coefbn1, coefbm1)
            end do
        end do
        !$omp end parallel do
    end if

end subroutine deriche_gaussian_filt_3d_

!==========================================================================
! Gaussian blur (i.e., the 0th order derivative) implementation using
! Discrete Cosine Transform, which outputs exact Gaussian smoothing results
! but has the highest computational cost. DCT based Gaussian filtering
! has no problem of boundary handling.

!
!> 1D Gaussian filtering by DCT
!
subroutine dct_gaussian_filt_1d_(w, sigma)

    TT, dimension(:), intent(inout) :: w
    TT, intent(in) :: sigma

    integer :: i, n
    TT :: cf

    n = size(w)
    if (n == 1) then
        return
    end if
    cf = -2*const_pi**2*sigma**2/(2.0*n)**2

    call cosine_transform(w)
    do i = 1, n
        w(i) = w(i)*exp(cf*(i - 1.0)**2)
    end do
    call inverse_cosine_transform(w)

end subroutine dct_gaussian_filt_1d_

!
!> 2D Gaussian filtering by DCT
!
subroutine dct_gaussian_filt_2d_(w, sigma)

    TT, dimension(:, :), intent(inout) :: w
    TT, dimension(:), intent(in) :: sigma

    integer :: i, j, n1, n2
    TT :: cf1, cf2

    n1 = size(w, 1)
    n2 = size(w, 2)
    cf1 = -2*const_pi**2*sigma(1)**2/(2.0*n1)**2
    cf2 = -2*const_pi**2*sigma(2)**2/(2.0*n2)**2

    call cosine_transform(w)
    !$omp parallel do private(i, j)
    do j = 1, n2
        do i = 1, n1
            w(i, j) = w(i, j) &
                *exp(cf1*(i - 1.0)**2) &
                *exp(cf2*(j - 1.0)**2)
        end do
    end do
    !$omp end parallel do
    call inverse_cosine_transform(w)

    !
    ! An alternative approach is doing DCT Gaussian filtering
    ! for each dimension strides
    ! Based on the separability of Gaussian filtering,
    ! this approach is equivalent to the 2D DCT approach
    !
    !  !$omp parallel do private(i)
    !  do i = 1,size(w,1)
    !      call dct_gaussian_filt_1d_(w(i,:), sigma(2))
    !  end do
    !  !$omp end parallel do
    !
    !  !$omp parallel do private(j)
    !  do j = 1,size(w,2)
    !      call dct_gaussian_filt_1d_(w(:,j), sigma(1))
    !  end do
    !  !$omp end parallel do
    !

end subroutine dct_gaussian_filt_2d_

!
!> 3D Gaussian filtering by DCT
!
subroutine dct_gaussian_filt_3d_(w, sigma)

    TT, dimension(:, :, :), intent(inout) :: w
    TT, dimension(:), intent(in) :: sigma

    integer :: i, j, k, n1, n2, n3
    TT :: cf1, cf2, cf3

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)
    cf1 = -2*const_pi**2*sigma(1)**2/(2.0*n1)**2
    cf2 = -2*const_pi**2*sigma(2)**2/(2.0*n2)**2
    cf3 = -2*const_pi**2*sigma(3)**2/(2.0*n3)**2

    call cosine_transform(w)
    !$omp parallel do private(i, j, k)
    do k = 1, n3
        do j = 1, n2
            do i = 1, n1
                w(i, j, k) = w(i, j, k) &
                    *exp(cf1*(i - 1.0)**2) &
                    *exp(cf2*(j - 1.0)**2) &
                    *exp(cf3*(k - 1.0)**2)
            end do
        end do
    end do
    !$omp end parallel do
    call inverse_cosine_transform(w)

end subroutine dct_gaussian_filt_3d_

!==========================================================================
! Gaussian filtering by DFT
!

!
!> 1D Gaussian filtering by DFT
!
subroutine dft_gaussian_filt_1d_(w, sigma)

    TT, dimension(:), intent(inout) :: w
    TT, intent(in) :: sigma

    integer :: i, n
    TTT, allocatable, dimension(:) :: wt

    if (sigma > 0) then

        n = size(w)
        allocate (wt(1:2*n))
        wt(1:2*n) = cmplx([w, w(n:1:-1)])

        call fourier_transform(wt)
        do i = 1, n
            wt(i) = wt(i)*exp(-2*const_pi**2*sigma**2*((i - 1.0)/(2.0*n))**2)
        end do
        do i = n + 1, 2*n
            wt(i) = wt(i)*exp(-2*const_pi**2*sigma**2*((i - 1.0)/(2.0*n) - 1.0)**2)
        end do
        call inverse_fourier_transform(wt)

        w = TTTT(wt(1:n))

    end if

end subroutine dft_gaussian_filt_1d_

!
!> 2D Gaussian filtering by DFT
!
subroutine dft_gaussian_filt_2d_(w, sigma)

    TT, dimension(:, :), intent(inout) :: w
    TT, dimension(:), intent(in) :: sigma

    integer :: i, j, n1, n2

    n1 = size(w, 1)
    n2 = size(w, 2)

    if (sigma(2) > 0) then
        !$omp parallel do private(i)
        do i = 1, n1
            call dft_gaussian_filt_1d_(w(i, :), sigma(2))
        end do
        !$omp end parallel do
    end if

    if (sigma(1) > 0) then
        !$omp parallel do private(j)
        do j = 1, n2
            call dft_gaussian_filt_1d_(w(:, j), sigma(1))
        end do
        !$omp end parallel do
    end if

end subroutine dft_gaussian_filt_2d_

!
!> 3D Gaussian filtering by DFT
!
subroutine dft_gaussian_filt_3d_(w, sigma)

    TT, dimension(:, :, :), intent(inout) :: w
    TT, dimension(:), intent(in) :: sigma

    integer :: n1, n2, n3, i, j, k

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)

    if (sigma(3) > 0) then
        !$omp parallel do private(i, j)
        do j = 1, n2
            do i = 1, n1
                call dft_gaussian_filt_1d_(w(i, j, :), sigma(3))
            end do
        end do
        !$omp end parallel do
    end if

    if (sigma(2) > 0) then
        !$omp parallel do private(i, k)
        do k = 1, n3
            do i = 1, n1
                call dft_gaussian_filt_1d_(w(i, :, k), sigma(2))
            end do
        end do
        !$omp end parallel do
    end if

    if (sigma(1) > 0) then
        !$omp parallel do private(j, k)
        do k = 1, n3
            do j = 1, n2
                call dft_gaussian_filt_1d_(w(:, j, k), sigma(1))
            end do
        end do
        !$omp end parallel do
    end if

end subroutine dft_gaussian_filt_3d_

!=========================================================================
! Recursive spectral Gaussian filter, has the overall best performance:
! speed best + accuracy just a slightly worse than DFT/DCT-based
!
! See Sugimoto and Kamata, 2013, Fast Gaussian filter with second-order shift property of DCT-5
! https://ieeexplore.ieee.org/document/6738106/
!
function athead_(x) result(xx)

    integer :: x
    integer :: xx

    xx = abs(x)

end function athead_

function attail_(x, n) result(xx)

    integer :: x, n
    integer :: xx

    xx = n - 1 - abs(n - 1 - x)

end function attail_

subroutine compute_rsgf_coef_(sigma, accuracy, radius, spectrum, table)

    TT, intent(in) :: sigma
    integer, intent(in) :: accuracy
    integer, intent(inout) :: radius
    TT, allocatable, dimension(:), intent(inout) :: spectrum, table

    TT :: phi
    integer :: k, u

    ! estimate radius
    select case (accuracy)
        case (2)
            if (sigma < 4.0) then
                radius = nint(3.0*sigma - 0.2 + 0.5)
            else
                radius = nint(3.0*sigma + 0.5)
            end if
        case (3)
            if (sigma < 4.0) then
                radius = nint(3.3333*sigma - 0.3333 + 0.5)
            else
                radius = nint(3.4113*sigma - 0.6452 + 0.5)
            end if
    end select

    ! phase
    phi = 2.0*const_pi/(radius + 1 + radius)

    ! generate spectrum
    allocate (spectrum(1:accuracy))
    do k = 1, accuracy
        spectrum(k) = 2.0*exp(-0.5*sigma**2*phi**2*k**2)
    end do

    ! build look-up table
    allocate (table(1:accuracy*(1 + radius)))
    do u = 0, radius
        do k = 1, accuracy
            table(accuracy*u + k) = cos(k*phi*u)*spectrum(k)
        end do
    end do

end subroutine compute_rsgf_coef_

subroutine rsgf_1d_(w, k, r, spectrum, table)

    TT, dimension(:), intent(inout) :: w
    integer, intent(in) :: k, r ! accuracy and radius
    TT, dimension(0:), intent(in) :: spectrum, table

    integer :: n
    TT :: norm
    TT :: cf11, cf12, cfR1, cfR2
    TT :: sum, a1, a2, b1, b2
    TT :: dA, dB, delta
    TT :: sumA, sumB
    TT, allocatable, dimension(:) :: buf
    integer :: x, u

    n = size(w)
    if (n == 1) then
        return
    end if

    norm = 1.0/(r + 1 + r)
    cf11 = table(k*1 + 0)*2.0/spectrum(0)
    cf12 = table(k*1 + 1)*2.0/spectrum(1)
    cfR1 = table(k*r + 0)
    cfR2 = table(k*r + 1)

    allocate (buf(0:n - 1), source=w(1:n))

    sum = buf(0)
    a1 = buf(0)*table(0)
    b1 = buf(1)*table(0)
    a2 = buf(0)*table(1)
    b2 = buf(1)*table(1)

    do u = 1, r
        sumA = buf(athead_(0 - u)) + buf(0 + u)
        sumB = buf(athead_(1 - u)) + buf(1 + u)
        sum = sum + sumA
        a1 = a1 + sumA*table(k*u + 0)
        a2 = a2 + sumA*table(k*u + 1)
        b1 = b1 + sumB*table(k*u + 0)
        b2 = b2 + sumB*table(k*u + 1)
    end do

    ! the first pixel
    w(0 + 1) = norm*(sum + a1 + a2)
    dA = buf(attail_(0 + r + 1, n)) - buf(athead_(0 - r))
    sum = sum + dA

    ! the other pixels
    x = 1
    do !four-length ring buffers

        w(x + 1) = norm*(sum + b1 + b2)
        dB = buf(attail_(x + r + 1, n)) - buf(athead_(x - r))
        delta = dA - dB
        sum = sum + dB
        a1 = a1 - cf11*b1 + cfR1*delta
        a2 = a2 - cf12*b2 + cfR2*delta
        x = x + 1
        if (x >= n) then
            exit
        end if

        w(x + 1) = norm*(sum - a1 - a2)
        dA = buf(attail_(x + r + 1, n)) - buf(athead_(x - r))
        delta = dB - dA
        sum = sum + dA
        b1 = b1 + cf11*a1 + cfR1*delta
        b2 = b2 + cf12*a2 + cfR2*delta
        x = x + 1
        if (x >= n) then
            exit
        end if

        w(x + 1) = norm*(sum - b1 - b2)
        dB = buf(attail_(x + r + 1, n)) - buf(athead_(x - r))
        delta = dA - dB
        sum = sum + dB
        a1 = a1 - cf11*b1 - cfR1*delta
        a2 = a2 - cf12*b2 - cfR2*delta
        x = x + 1
        if (x >= n) then
            exit
        end if

        w(x + 1) = norm*(sum + a1 + a2)
        dA = buf(attail_(x + r + 1, n)) - buf(athead_(x - r))
        delta = dB - dA
        sum = sum + dA
        b1 = b1 + cf11*a1 - cfR1*delta
        b2 = b2 + cf12*a2 - cfR2*delta
        x = x + 1
        if (x >= n) then
            exit
        end if

    end do

end subroutine rsgf_1d_

subroutine spectral_gaussian_filt_1d_(w, sigma)

    TT, dimension(:), intent(inout) :: w
    TT, intent(in) :: sigma

    integer :: radius
    TT, allocatable, dimension(:) :: spectrum, table

    if (sigma > 0) then
        call compute_rsgf_coef_(sigma, 2, radius, spectrum, table)
        call rsgf_1d_(w, 2, radius, spectrum, table)
    end if

end subroutine spectral_gaussian_filt_1d_

subroutine spectral_gaussian_filt_2d_(w, sigma)

    TT, dimension(:, :), intent(inout) :: w
    TT, dimension(:), intent(in) :: sigma

    integer :: radius1, radius2
    TT, allocatable, dimension(:) :: spectrum1, table1
    TT, allocatable, dimension(:) :: spectrum2, table2
    integer :: i, j, n1, n2

    n1 = size(w, 1)
    n2 = size(w, 2)

    ! along 2
    if (sigma(2) > 0) then
        call compute_rsgf_coef_(sigma(2), 2, radius2, spectrum2, table2)
        !$omp parallel do private(i)
        do i = 1, n1
            call rsgf_1d_(w(i, :), 2, radius2, spectrum2, table2)
        end do
        !$omp end parallel do
    end if

    if (sigma(1) > 0) then
        call compute_rsgf_coef_(sigma(1), 2, radius1, spectrum1, table1)
        !$omp parallel do private(j)
        do j = 1, n2
            call rsgf_1d_(w(:, j), 2, radius1, spectrum1, table1)
        end do
        !$omp end parallel do
    end if

end subroutine spectral_gaussian_filt_2d_

subroutine spectral_gaussian_filt_3d_(w, sigma)

    TT, dimension(:, :, :), intent(inout) :: w
    TT, dimension(:), intent(in) :: sigma

    integer :: radius1, radius2, radius3
    TT, allocatable, dimension(:) :: spectrum1, table1
    TT, allocatable, dimension(:) :: spectrum2, table2
    TT, allocatable, dimension(:) :: spectrum3, table3
    integer :: i, j, k, n1, n2, n3

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)

    if (sigma(3) > 0) then
        call compute_rsgf_coef_(sigma(3), 2, radius3, spectrum3, table3)
        !$omp parallel do private(i, j)
        do j = 1, n2
            do i = 1, n1
                call rsgf_1d_(w(i, j, :), 2, radius3, spectrum3, table3)
            end do
        end do
        !$omp end parallel do
    end if

    if (sigma(2) > 0) then
        call compute_rsgf_coef_(sigma(2), 2, radius2, spectrum2, table2)
        !$omp parallel do private(i, k)
        do k = 1, n3
            do i = 1, n1
                call rsgf_1d_(w(i, :, k), 2, radius2, spectrum2, table2)
            end do
        end do
        !$omp end parallel do
    end if

    if (sigma(1) > 0) then
        call compute_rsgf_coef_(sigma(1), 2, radius1, spectrum1, table1)
        !$omp parallel do private(j, k)
        do k = 1, n3
            do j = 1, n2
                call rsgf_1d_(w(:, j, k), 2, radius1, spectrum1, table1)
            end do
        end do
        !$omp end parallel do
    end if

end subroutine spectral_gaussian_filt_3d_

!==========================================================================
! Gaussian blur by convolution

!
!> 1D Gaussian filtering by convolution
!
function conv_gaussian_filt_1d_(w, sigma) result(v)

    TT, dimension(:) :: w
    TT, intent(in) :: sigma
    TT, allocatable, dimension(:) :: v

    integer :: i, n1
    TT, allocatable, dimension(:) :: g
    integer :: p1

    call assert(sigma > 0, " <conv_gaussian_filt_1d> Error: sigma must > 0. ")

    n1 = size(w)
    p1 = nint(n1/4.0)

    g = zeros(n1)
    call pad_array(g, [p1, p1])
    do i = -p1 + 1, n1 + p1
        g(i) = exp(-0.5*(i - 0.5*(n1 + 2.0))**2/sigma**2)
    end do
    g = g/sum(g)

    v = conv(pad(w, [p1, p1], ['symm', 'symm']), g, method='same')
    v = v(p1 + 1:p1 + n1)

end function conv_gaussian_filt_1d_

function conv_gaussian_filt_2d_(w, sigma) result(v)

    TT, dimension(:, :) :: w
    TT, dimension(2), intent(in) :: sigma
    TT, allocatable, dimension(:, :) :: v

    integer :: i, j, n1, n2
    TT, allocatable, dimension(:, :) :: g
    integer :: p1, p2

    call assert(all(sigma > 0), " <conv_gaussian_filt_2d> Error: all sigma's must > 0. ")

    n1 = size(w, 1)
    n2 = size(w, 2)
    p1 = nint(n1/4.0)
    p2 = nint(n2/4.0)

    g = zeros(n1, n2)
    call pad_array(g, [p1, p1, p2, p2])
    !$omp parallel do private(i, j)
    do j = -p2 + 1, n2 + p2
        do i = -p1 + 1, n1 + p1
            g(i, j) = exp(-0.5*((i - 0.5*(n1 + 2.0))**2/sigma(1)**2 &
                + (j - 0.5*(n2 + 2.0))**2/sigma(2)**2))
        end do
    end do
    !$omp end parallel do
    g = g/sum(g)

    v = conv(pad(w, [p1, p1, p2, p2], ['symm', 'symm', 'symm', 'symm']), g, method='same')
    v = v(p1 + 1:p1 + n1, p2 + 1:p2 + n2)

end function conv_gaussian_filt_2d_

function conv_gaussian_filt_3d_(w, sigma) result(v)

    TT, dimension(:, :, :) :: w
    TT, dimension(3), intent(in) :: sigma
    TT, allocatable, dimension(:, :, :) :: v

    integer :: i, j, k, n1, n2, n3
    TT, allocatable, dimension(:, :, :) :: g
    integer :: p1, p2, p3

    call assert(all(sigma > 0), " <conv_gaussian_filt_3d> Error: all sigma's must > 0. ")

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)
    p1 = nint(n1/4.0)
    p2 = nint(n2/4.0)
    p3 = nint(n3/4.0)

    g = zeros(n1, n2, n3)
    call pad_array(g, [p1, p1, p2, p2, p3, p3])
    !$omp parallel do private(i, j, k)
    do k = -p3 + 1, n3 + p3
        do j = -p2 + 1, n2 + p2
            do i = -p1 + 1, n1 + p1
                g(i, j, k) = exp(-0.5*((i - 0.5*(n1 + 2.0))**2/sigma(1)**2 &
                    + (j - 0.5*(n2 + 2.0))**2/sigma(2)**2 &
                    + (k - 0.5*(n3 + 2.0))**2/sigma(3)**2))
            end do
        end do
    end do
    !$omp end parallel do
    g = g/sum(g)

    v = conv(pad(w, [p1, p1, p2, p2, p3, p3], ['symm', 'symm', 'symm', 'symm', 'symm', 'symm']), g, method='same')
    v = v(p1 + 1:p1 + n1, p2 + 1:p2 + n2, p3 + 1:p3 + n3)

end function conv_gaussian_filt_3d_

!==========================================================================
! Interfaces
function gauss_filt_1d_(w, sigma, order, method) result(wr)

    TT, dimension(:), intent(in) :: w
    TT, intent(in) :: sigma
    integer, intent(in), optional :: order
    character(len=*), intent(in), optional :: method
    TT, allocatable, dimension(:) :: wr

    integer :: gauss_filt_order
    character(len=12) :: gauss_filt_method

    if (present(method)) then
        gauss_filt_method = method
    else
        gauss_filt_method = 'rsf'
    end if

    if (present(order)) then
        gauss_filt_order = order
    else
        gauss_filt_order = 0
    end if

    wr = w

    select case (gauss_filt_method)
        case ('deriche')
            call deriche_gaussian_filt_1d_(wr, sigma, gauss_filt_order)
        case ('dct')
            call dct_gaussian_filt_1d_(wr, sigma)
        case ('dft')
            call dft_gaussian_filt_1d_(wr, sigma)
        case ('rsf')
            call spectral_gaussian_filt_1d_(wr, sigma)
        case ('conv')
            wr = conv_gaussian_filt_1d_(wr, sigma)
    end select

end function gauss_filt_1d_

function gauss_filt_2d_(w, sigma, order, method) result(wr)

    TT, dimension(:, :), intent(in) :: w
    TT, dimension(1:2), intent(in) :: sigma
    integer, dimension(1:2), intent(in), optional :: order
    character(len=*), intent(in), optional :: method
    TT, allocatable, dimension(:, :) :: wr

    integer, dimension(1:2) :: gauss_filt_order
    character(len=12) :: gauss_filt_method

    if (present(method)) then
        gauss_filt_method = method
    else
        gauss_filt_method = 'rsf'
    end if

    if (present(order)) then
        gauss_filt_order = order
    else
        gauss_filt_order = [0, 0]
    end if

    wr = w

    select case (gauss_filt_method)
        case ('deriche')
            call deriche_gaussian_filt_2d_(wr, sigma, gauss_filt_order)
        case ('dct')
            call dct_gaussian_filt_2d_(wr, sigma)
        case ('dft')
            call dft_gaussian_filt_2d_(wr, sigma)
        case ('rsf')
            call spectral_gaussian_filt_2d_(wr, sigma)
        case ('conv')
            wr = conv_gaussian_filt_2d_(wr, sigma)
    end select

end function gauss_filt_2d_

function gauss_filt_3d_(w, sigma, order, method) result(wr)

    TT, dimension(:, :, :), intent(in) :: w
    TT, dimension(1:3), intent(in) :: sigma
    integer, dimension(1:3), intent(in), optional :: order
    character(len=*), intent(in), optional :: method
    TT, allocatable, dimension(:, :, :) :: wr

    integer, dimension(1:3) :: gauss_filt_order
    character(len=12) :: gauss_filt_method

    if (present(method)) then
        gauss_filt_method = method
    else
        gauss_filt_method = 'rsf'
    end if

    if (present(order)) then
        gauss_filt_order = order
    else
        gauss_filt_order = [0, 0, 0]
    end if

    wr = w

    select case (gauss_filt_method)
        case ('deriche')
            call deriche_gaussian_filt_3d_(wr, sigma, gauss_filt_order)
        case ('dct')
            call dct_gaussian_filt_3d_(wr, sigma)
        case ('dft')
            call dft_gaussian_filt_3d_(wr, sigma)
        case ('rsf')
            call spectral_gaussian_filt_3d_(wr, sigma)
        case ('conv')
            wr = conv_gaussian_filt_3d_(wr, sigma)
    end select

end function gauss_filt_3d_

#undef T
#undef TT
#undef TTT
#undef TTTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef compute_deriche_coefs_
#undef apply_deriche_filter_1d_
#undef deriche_gaussian_filt_1d_
#undef deriche_gaussian_filt_2d_
#undef deriche_gaussian_filt_3d_
#undef dct_gaussian_filt_1d_
#undef dct_gaussian_filt_2d_
#undef dct_gaussian_filt_3d_
#undef dft_gaussian_filt_1d_
#undef dft_gaussian_filt_2d_
#undef dft_gaussian_filt_3d_
#undef athead_
#undef attail_
#undef compute_rsgf_coef_
#undef rsgf_1d_
#undef spectral_gaussian_filt_1d_
#undef spectral_gaussian_filt_2d_
#undef spectral_gaussian_filt_3d_
#undef conv_gaussian_filt_1d_
#undef conv_gaussian_filt_2d_
#undef conv_gaussian_filt_3d_
#undef gauss_filt_1d_
#undef gauss_filt_2d_
#undef gauss_filt_3d_
