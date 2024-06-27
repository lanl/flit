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

#define kernel_triangular_      CONCAT(kernel_triangular, T)
#define kernel_epanechnikov_      CONCAT(kernel_epanechnikov, T)
#define kernel_bisquare_      CONCAT(kernel_bisquare, T)
#define kernel_trisquare_      CONCAT(kernel_trisquare, T)
#define kernel_tricube_      CONCAT(kernel_tricube, T)
#define kernel_cosine_      CONCAT(kernel_cosine, T)
#define covariance_matrix_1d_      CONCAT(covariance_matrix_1d, T)
#define covariance_matrix_2d_      CONCAT(covariance_matrix_2d, T)
#define covariance_matrix_3d_      CONCAT(covariance_matrix_3d, T)
#define row_covariance_matrix_2d_      CONCAT(row_covariance_matrix_2d, T)
#define cross_covariance_matrix_1d_      CONCAT(cross_covariance_matrix_1d, T)
#define cross_covariance_matrix_2d_      CONCAT(cross_covariance_matrix_2d, T)
#define cross_covariance_matrix_3d_      CONCAT(cross_covariance_matrix_3d, T)
#define row_cross_covariance_matrix_2d_      CONCAT(row_cross_covariance_matrix_2d, T)
#define variance_1d_      CONCAT(variance_1d, T)
#define variance_2d_      CONCAT(variance_2d, T)
#define variance_3d_      CONCAT(variance_3d, T)
#define gaussian_1d_      CONCAT(gaussian_1d, T)
#define gaussian_2d_      CONCAT(gaussian_2d, T)
#define gaussian_3d_      CONCAT(gaussian_3d, T)
#define covariance_1d_      CONCAT(covariance_1d, T)
#define covariance_2d_      CONCAT(covariance_2d, T)
#define covariance_3d_      CONCAT(covariance_3d, T)
#define standard_deviation_1d_      CONCAT(standard_deviation_1d, T)
#define standard_deviation_2d_      CONCAT(standard_deviation_2d, T)
#define standard_deviation_3d_      CONCAT(standard_deviation_3d, T)
#define histogram_1d_      CONCAT(histogram_1d, T)
#define histogram_2d_      CONCAT(histogram_2d, T)
#define histogram_3d_      CONCAT(histogram_3d, T)
#define plot_histogram_1d_      CONCAT(plot_histogram_1d, T)
#define plot_histogram_2d_      CONCAT(plot_histogram_2d, T)
#define plot_histogram_3d_      CONCAT(plot_histogram_3d, T)
#define median_1d_      CONCAT(median_1d, T)
#define median_2d_      CONCAT(median_2d, T)
#define median_3d_      CONCAT(median_3d, T)
#define discrete_xcorr_1d_      CONCAT(discrete_xcorr_1d, T)
#define discrete_xcorr_2d_      CONCAT(discrete_xcorr_2d, T)
#define discrete_xcorr_3d_      CONCAT(discrete_xcorr_3d, T)
#define xcorr_1d_      CONCAT(xcorr_1d, T)
#define xcorr_2d_      CONCAT(xcorr_2d, T)
#define xcorr_3d_      CONCAT(xcorr_3d, T)
#define discrete_acorr_1d_      CONCAT(discrete_acorr_1d, T)
#define discrete_acorr_2d_      CONCAT(discrete_acorr_2d, T)
#define discrete_acorr_3d_      CONCAT(discrete_acorr_3d, T)
#define acorr_1d_      CONCAT(acorr_1d, T)
#define acorr_2d_      CONCAT(acorr_2d, T)
#define acorr_3d_      CONCAT(acorr_3d, T)
#define local_xcorr_1d_      CONCAT(local_xcorr_1d, T)
#define xcorr_coef_1d_      CONCAT(xcorr_coef_1d, T)
#define xcorr_coef_2d_      CONCAT(xcorr_coef_2d, T)
#define xcorr_coef_3d_      CONCAT(xcorr_coef_3d, T)

elemental function kernel_triangular_(u) result(w)

    TT, intent(in) :: u
    TT :: w

    w = abs(u)
    if (w < 1) then
        w = 1 - w
    else
        w = 0
    end if

end function kernel_triangular_

elemental function kernel_epanechnikov_(u) result(w)

    TT, intent(in) :: u
    TT :: w

    w = abs(u)
    if (w < 1) then
        w = 0.75*(1 - w**2)
    else
        w = 0
    end if

end function kernel_epanechnikov_

elemental function kernel_bisquare_(u) result(w)

    TT, intent(in) :: u
    TT :: w

    w = abs(u)
    if (w < 1) then
        w = (1 - w**2)**2
    else
        w = 0
    end if

end function kernel_bisquare_

elemental function kernel_trisquare_(u) result(w)

    TT, intent(in) :: u
    TT :: w

    w = abs(u)
    if (w < 1) then
        w = (1 - w**2)**3
    else
        w = 0
    end if

end function kernel_trisquare_

elemental function kernel_tricube_(u) result(w)

    TT, intent(in) :: u
    TT :: w

    w = abs(u)
    if (w < 1) then
        w = (1 - w**3)**3
    else
        w = 0
    end if

end function kernel_tricube_

elemental function kernel_cosine_(u) result(w)

    TT, intent(in) :: u
    TT :: w

    w = abs(u)
    if (w < 1) then
        w = const_pi/4.0*cos(const_pi/2.0*w)
    else
        w = 0
    end if

end function kernel_cosine_

!
! @brief Variance-covariance matrix of a vector x
!
function covariance_matrix_1d_(x) result(c)

    TT, dimension(:), intent(in) :: x
    TT, allocatable, dimension(:, :) :: c

    c = ones(2, 2)*var(x)

end function covariance_matrix_1d_

function covariance_matrix_2d_(x) result(c)

    TT, dimension(:, :), intent(in) :: x
    TT, allocatable, dimension(:, :) :: c

    c = ones(2, 2)*var(x)

end function covariance_matrix_2d_

function covariance_matrix_3d_(x) result(c)

    TT, dimension(:, :, :), intent(in) :: x
    TT, allocatable, dimension(:, :) :: c

    c = ones(2, 2)*var(x)

end function covariance_matrix_3d_

!
! @brief Variance-covariance matrix of n sets variants
!        Each column of x is an observation/realization
!
function row_covariance_matrix_2d_(x, rowvar) result(c)

    TT, dimension(:, :), intent(in) :: x
    logical, intent(in) :: rowvar
    TT, allocatable, dimension(:, :) :: c

    integer :: n, i, j

    if (rowvar) then
        n = size(x, 1)
        c = zeros(n, n)
        do j = 1, n
            do i = 1, n
                if (j >= i) then
                    c(i, j) = covar(x(i, :), x(j, :))
                end if
            end do
        end do
    else
        n = size(x, 2)
        c = zeros(n, n)
        do j = 1, n
            do i = 1, n
                if (j >= i) then
                    c(i, j) = covar(x(:, i), x(:, j))
                end if
            end do
        end do
    end if

    do j = 1, n
        do i = 1, n
            if (j < i) then
                c(i, j) = c(j, i)
            end if
        end do
    end do

end function row_covariance_matrix_2d_

function cross_covariance_matrix_1d_(x, y) result(c)

    TT, dimension(:), intent(in) :: x, y
    TT, allocatable, dimension(:, :) :: c

    c = zeros(2, 2)
    c(1, :) = [var(x), covar(x, y)]
    c(2, :) = [covar(y, x), var(y)]

end function cross_covariance_matrix_1d_

function cross_covariance_matrix_2d_(x, y) result(c)

    TT, dimension(:, :), intent(in) :: x, y
    TT, allocatable, dimension(:, :) :: c

    c = zeros(2, 2)
    c(1, :) = [var(x), covar(x, y)]
    c(2, :) = [covar(y, x), var(y)]

end function cross_covariance_matrix_2d_

function cross_covariance_matrix_3d_(x, y) result(c)

    TT, dimension(:, :, :), intent(in) :: x, y
    TT, allocatable, dimension(:, :) :: c

    c = zeros(2, 2)
    c(1, :) = [var(x), covar(x, y)]
    c(2, :) = [covar(y, x), var(y)]

end function cross_covariance_matrix_3d_

function row_cross_covariance_matrix_2d_(x, y, rowvar) result(c)

    TT, dimension(:, :), intent(in) :: x, y
    logical, intent(in) :: rowvar
    TT, allocatable, dimension(:, :) :: c

    integer :: n, i, j

    call assert(size(x, 1) == size(y, 1) .and. size(x, 2) == size(y, 2), &
        ' <cross_covariance_matrix_2d_> Error: x and y must have same shape.')

    if (rowvar) then
        n = size(x, 1)
        c = zeros(n, n)
        do j = 1, n
            do i = 1, n
                if (j >= i) then
                    c(i, j) = covar(x(i, :), y(j, :))
                end if
            end do
        end do
    else
        n = size(x, 2)
        c = zeros(n, n)
        do j = 1, n
            do i = 1, n
                if (j >= i) then
                    c(i, j) = covar(x(:, i), y(:, j))
                end if
            end do
        end do
    end if

    do j = 1, n
        do i = 1, n
            if (j < i) then
                c(i, j) = c(j, i)
            end if
        end do
    end do

end function row_cross_covariance_matrix_2d_

!
!> Compute the variance of a 1D array
!
function variance_1d_(w) result(s)

    TT, dimension(:), intent(in) :: w
    TT :: s

    s = sum((w - mean(w))**2)
    if (size(w) > 1) then
        s = s/(size(w) - 1)
    end if

end function variance_1d_

!
!> Compute the variance of a 2D array
!
function variance_2d_(w) result(s)

    TT, dimension(:, :), intent(in) :: w
    TT :: s

    s = sum((w - mean(w))**2)
    if (size(w) > 1) then
        s = s/(size(w) - 1)
    end if

end function variance_2d_

!
!> Compute the variance of a 3D array
!
function variance_3d_(w) result(s)

    TT, dimension(:, :, :), intent(in) :: w
    TT :: s

    s = sum((w - mean(w))**2)
    if (size(w) > 1) then
        s = s/(size(w) - 1)
    end if

end function variance_3d_

function gaussian_1d_(mu, sigma, x) result(f)

    TT, intent(in) :: mu, sigma
    TT, dimension(:), intent(in) :: x
    TT, allocatable, dimension(:) :: f

    f = exp(-0.5*((x - mu)/sigma)**2)/(sigma*sqrt(2*const_pi))

end function gaussian_1d_

function gaussian_2d_(mu, sigma, x, y) result(f)

    TT, dimension(1:2), intent(in) :: mu, sigma
    TT, dimension(:), intent(in) :: x, y
    TT, allocatable, dimension(:, :) :: f

    TT, allocatable, dimension(:) :: f1, f2
    integer :: i, j

    f1 = exp(-0.5*((x - mu(1))/sigma(1))**2)/(sigma(1)*sqrt(2*const_pi))
    f2 = exp(-0.5*((y - mu(2))/sigma(2))**2)/(sigma(2)*sqrt(2*const_pi))

    f = zeros(size(f1), size(f2))
    !$omp parallel do private(i, j)
    do j = 1, size(f2)
        do i = 1, size(f1)
            f(i, j) = f1(i)*f2(j)
        end do
    end do
    !$omp end parallel do

end function gaussian_2d_

function gaussian_3d_(mu, sigma, x, y, z) result(f)

    TT, dimension(1:3), intent(in) :: mu, sigma
    TT, dimension(:), intent(in) :: x, y, z
    TT, allocatable, dimension(:, :, :) :: f

    TT, allocatable, dimension(:) :: f1, f2, f3
    integer :: i, j, k

    f1 = exp(-0.5*((x - mu(1))/sigma(1))**2)/(sigma(1)*sqrt(2*const_pi))
    f2 = exp(-0.5*((y - mu(2))/sigma(2))**2)/(sigma(2)*sqrt(2*const_pi))
    f3 = exp(-0.5*((z - mu(3))/sigma(3))**2)/(sigma(3)*sqrt(2*const_pi))

    f = zeros(size(f1), size(f2), size(f3))
    !$omp parallel do private(i, j, k)
    do k = 1, size(f3)
        do j = 1, size(f2)
            do i = 1, size(f1)
                f(i, j, k) = f1(i)*f2(j)*f3(k)
            end do
        end do
    end do
    !$omp end parallel do

end function gaussian_3d_


function covariance_1d_(x, y) result(c)

    TT, dimension(:), intent(in) :: x, y
    TT :: c

    TT :: mx, my
    integer :: nx, ny

    nx = size(x)
    ny = size(y)

    call assert(nx == ny, ' <covariance_1d_> Error: size(x) /= size(y)')

    mx = mean(x)
    my = mean(y)
    c = sum((x - mx)*(y - my))
    if (nx > 1) then
        c = c/(nx - 1.0)
    end if

end function covariance_1d_

function covariance_2d_(x, y) result(c)

    TT, dimension(:, :), intent(in) :: x, y
    TT :: c

    TT :: mx, my
    integer :: nx1, ny1, nx2, ny2

    nx1 = size(x, 1)
    ny1 = size(y, 1)
    call assert(nx1 == ny1, ' <covariance_2d_> Error: size(x, 1) /= size(y, 1)')
    nx2 = size(x, 2)
    ny2 = size(y, 2)
    call assert(nx2 == ny2, ' <covariance_2d_> Error: size(x, 2) /= size(y, 2)')

    mx = mean(x)
    my = mean(y)
    c = sum((x - mx)*(y - my))
    if (nx1*nx2 > 1) then
        c = c/(nx1*nx2 - 1.0)
    end if

end function covariance_2d_

function covariance_3d_(x, y) result(c)

    TT, dimension(:, :, :), intent(in) :: x, y
    TT :: c

    TT :: mx, my
    integer :: nx1, ny1, nx2, ny2, nx3, ny3

    nx1 = size(x, 1)
    ny1 = size(y, 1)
    call assert(nx1 == ny1, ' <covariance_3d_> Error: size(x, 1) /= size(y, 1)')
    nx2 = size(x, 2)
    ny2 = size(y, 2)
    call assert(nx2 == ny2, ' <covariance_3d_> Error: size(x, 2) /= size(y, 2)')
    nx3 = size(x, 3)
    ny3 = size(y, 3)
    call assert(nx3 == ny3, ' <covariance_3d_> Error: size(x, 3) /= size(y, 3)')

    mx = mean(x)
    my = mean(y)
    c = sum((x - mx)*(y - my))
    if (nx1*nx2*nx3 > 1) then
        c = c/(nx1*nx2*nx3 - 1.0)
    end if

end function covariance_3d_

!
!> Compute the unbiased standard deviation of a 1D array
!
function standard_deviation_1d_(w) result(s)

    TT, dimension(:), intent(in) :: w
    TT :: s

    s = sum((w - mean(w))**2)
    if (size(w) > 1) then
        s = sqrt(s/(size(w) - 1))
    end if

end function standard_deviation_1d_

!
!> Compute the unbiased standard deviation of a 2D array
!
function standard_deviation_2d_(w) result(s)

    TT, dimension(:, :), intent(in) :: w
    TT :: s

    s = sum((w - mean(w))**2)
    if (size(w) > 1) then
        s = sqrt(s/(size(w) - 1))
    end if

end function standard_deviation_2d_

!
!> Compute the unbiased standard deviation of a 3D array
!
function standard_deviation_3d_(w) result(s)

    TT, dimension(:, :, :), intent(in) :: w
    TT :: s

    s = sum((w - mean(w))**2)
    if (size(w) > 1) then
        s = sqrt(s/(size(w) - 1))
    end if

end function standard_deviation_3d_

!
!> Compute the histogram of a 1D array
!
subroutine histogram_1d_(w, hist, valmin, valmax, binsize)

    TT, dimension(:), intent(in) :: w
    TT, allocatable, dimension(:, :), intent(inout) :: hist
    TT, intent(in), optional :: valmin, valmax, binsize

    TT :: histmin, histmax, histbin
    integer :: nbin, i

    if (present(valmin)) then
        histmin = valmin
    else
        histmin = minval(w)
    end if

    if (present(valmax)) then
        histmax = valmax
    else
        histmax = maxval(w)
    end if

    if (present(binsize)) then
        histbin = binsize
        nbin = ceiling((histmax - histmin)/histbin)
    else
        nbin = 9
        histbin = (histmax - histmin)/nbin
    end if

    if (histbin == 0) then
        nbin = 1
    end if

    ! allocate hist
    call alloc_array(hist, [1, nbin, 1, 4])

    ! compute hist
    !$omp parallel do private(i)
    do i = 1, nbin

        if (i == 1) then
            hist(i, 1) = histmin
        else
            hist(i, 1) = histmin + (i - 1)*histbin
        end if

        if (i == nbin) then
            hist(i, 2) = histmax
        else
            hist(i, 2) = histmin + i*histbin
        end if

        if (i < nbin) then
            hist(i, 3) = count(w >= hist(i, 1) .and. w < hist(i, 2), kind=8)
        else
            hist(i, 3) = count(w >= hist(i, 1) .and. w <= hist(i, 2), kind=8)
        end if

    end do
    !$omp end parallel do

    ! normalize
    hist(:, 4) = hist(:, 3)/size(w, kind=8)

end subroutine histogram_1d_

!
!> Compute the histogram of a 2D array
!
subroutine histogram_2d_(w, hist, valmin, valmax, binsize)

    TT, dimension(:, :), intent(in) :: w
    TT, allocatable, dimension(:, :), intent(inout) :: hist
    TT, intent(in), optional :: valmin, valmax, binsize

    TT :: histmin, histmax, histbin
    integer :: nbin, i

    if (present(valmin)) then
        histmin = valmin
    else
        histmin = minval(w)
    end if

    if (present(valmax)) then
        histmax = valmax
    else
        histmax = maxval(w)
    end if

    if (present(binsize)) then
        histbin = binsize
        nbin = ceiling((histmax - histmin)/histbin)
    else
        nbin = 9
        histbin = (histmax - histmin)/nbin
    end if

    if (histbin == 0) then
        nbin = 1
    end if

    ! allocate hist
    call alloc_array(hist, [1, nbin, 1, 4])

    ! compute hist
    !$omp parallel do private(i)
    do i = 1, nbin

        if (i == 1) then
            hist(i, 1) = histmin
        else
            hist(i, 1) = histmin + (i - 1)*histbin
        end if

        if (i == nbin) then
            hist(i, 2) = histmax
        else
            hist(i, 2) = histmin + i*histbin
        end if

        if (i < nbin) then
            hist(i, 3) = count(w >= hist(i, 1) .and. w < hist(i, 2))
        else
            hist(i, 3) = count(w >= hist(i, 1) .and. w <= hist(i, 2))
        end if

    end do
    !$omp end parallel do

    ! normalize
    hist(:, 4) = hist(:, 3)/size(w, 1)/size(w, 2)

end subroutine histogram_2d_

!
!> Compute the histogram of a 3D array
!
subroutine histogram_3d_(w, hist, valmin, valmax, binsize)

    TT, dimension(:, :, :), intent(in) :: w
    TT, allocatable, dimension(:, :), intent(inout) :: hist
    TT, intent(in), optional :: valmin, valmax, binsize

    TT :: histmin, histmax, histbin
    integer :: nbin, i

    if (present(valmin)) then
        histmin = valmin
    else
        histmin = minval(w)
    end if

    if (present(valmax)) then
        histmax = valmax
    else
        histmax = maxval(w)
    end if

    if (present(binsize)) then
        histbin = binsize
        nbin = ceiling((histmax - histmin)/histbin)
    else
        nbin = 9
        histbin = (histmax - histmin)/nbin
    end if

    if (histbin == 0) then
        nbin = 1
    end if

    ! allocate hist
    call alloc_array(hist, [1, nbin, 1, 4])

    ! compute hist
    !$omp parallel do private(i)
    do i = 1, nbin

        if (i == 1) then
            hist(i, 1) = histmin
        else
            hist(i, 1) = histmin + (i - 1)*histbin
        end if

        if (i == nbin) then
            hist(i, 2) = histmax
        else
            hist(i, 2) = histmin + i*histbin
        end if

        if (i < nbin) then
            hist(i, 3) = count(w >= hist(i, 1) .and. w < hist(i, 2))
        else
            hist(i, 3) = count(w >= hist(i, 1) .and. w <= hist(i, 2))
        end if

    end do
    !$omp end parallel do

    ! normalize
    hist(:, 4) = hist(:, 3)/size(w, 1)/size(w, 2)/size(w, 3)

end subroutine histogram_3d_

!
!> Plot the histogram of a 1D array
!
subroutine plot_histogram_1d_(w, valmin, valmax, binsize, label)

    TT, dimension(:), intent(in) :: w
    TT, intent(in), optional :: valmin, valmax, binsize
    character(len=*), intent(in), optional :: label

    TT :: histmin, histmax, histbin
    integer :: nbin, i, j
    TT, allocatable, dimension(:, :) :: hist
    character(len=20) :: bar

    if (present(valmin)) then
        histmin = valmin
    else
        histmin = minval(w)
    end if

    if (present(valmax)) then
        histmax = valmax
    else
        histmax = maxval(w)
    end if

    if (present(binsize)) then
        histbin = binsize
        nbin = ceiling((histmax - histmin)/histbin)
    else
        nbin = 9
        histbin = (histmax - histmin)/nbin
    end if

    if (histbin == 0) then
        nbin = 1
    end if

    ! compute history
    call histogram_1d_(w, hist, histmin, histmax, histbin)
    nbin = size(hist, 1)

    ! compute hist
    if (present(label)) then
        write (error_unit, '(x,a)') label
    end if
    do i = 1, nbin
        bar = ''
        do j = 1, nint(hist(i, 4)/maxval(hist(:, 4))*len(bar))
            bar(j:j) = '*'
        end do
        write (error_unit, '(x,es,x,a,x,es,x,x,a,x,i12,x,es)') &
            hist(i, 1), '~', hist(i, 2), bar, int(hist(i, 3), kind=8), hist(i, 4)
    end do

end subroutine plot_histogram_1d_

!
!> Plot the histogram of a 2D array
!
subroutine plot_histogram_2d_(w, valmin, valmax, binsize, label)

    TT, dimension(:, :), intent(in) :: w
    TT, intent(in), optional :: valmin, valmax, binsize
    character(len=*), intent(in), optional :: label

    TT :: histmin, histmax, histbin
    integer :: nbin, i, j
    TT, allocatable, dimension(:, :) :: hist
    character(len=20) :: bar

    if (present(valmin)) then
        histmin = valmin
    else
        histmin = minval(w)
    end if

    if (present(valmax)) then
        histmax = valmax
    else
        histmax = maxval(w)
    end if

    if (present(binsize)) then
        histbin = binsize
        nbin = ceiling((histmax - histmin)/histbin)
    else
        nbin = 9
        histbin = (histmax - histmin)/nbin
    end if

    if (histbin == 0) then
        nbin = 1
    end if

    ! compute history
    call histogram_2d_(w, hist, histmin, histmax, histbin)
    nbin = size(hist, 1)

    ! compute hist
    if (present(label)) then
        write (error_unit, '(x,a)') label
    end if
    do i = 1, nbin
        bar = ''
        do j = 1, nint(hist(i, 4)/maxval(hist(:, 4))*len(bar))
            bar(j:j) = '*'
        end do
        write (error_unit, '(x,es,x,a,x,es,x,x,a,x,i12,x,es)') &
            hist(i, 1), '~', hist(i, 2), bar, int(hist(i, 3), kind=8), hist(i, 4)
    end do

end subroutine plot_histogram_2d_

!
!> Plot the histogram of a 3D array
!
subroutine plot_histogram_3d_(w, valmin, valmax, binsize, label)

    TT, dimension(:, :, :), intent(in) :: w
    TT, intent(in), optional :: valmin, valmax, binsize
    character(len=*), intent(in), optional :: label

    TT :: histmin, histmax, histbin
    integer :: nbin, i, j
    TT, allocatable, dimension(:, :) :: hist
    character(len=20) :: bar

    if (present(valmin)) then
        histmin = valmin
    else
        histmin = minval(w)
    end if

    if (present(valmax)) then
        histmax = valmax
    else
        histmax = maxval(w)
    end if

    if (present(binsize)) then
        histbin = binsize
        nbin = ceiling((histmax - histmin)/histbin)
    else
        nbin = 9
        histbin = (histmax - histmin)/nbin
    end if

    if (histbin == 0) then
        nbin = 1
    end if

    ! compute history
    call histogram_3d_(w, hist, histmin, histmax, histbin)
    nbin = size(hist, 1)

    ! compute hist
    if (present(label)) then
        write (error_unit, '(a)') label
    end if
    do i = 1, nbin
        bar = ''
        do j = 1, nint(hist(i, 4)/maxval(hist(:, 4))*len(bar))
            bar(j:j) = '*'
        end do
        write (error_unit, '(x,es,x,a,x,es,x,x,a,x,i12,x,es)') &
            hist(i, 1), '~', hist(i, 2), bar, int(hist(i, 3), kind=8), hist(i, 4)
    end do

end subroutine plot_histogram_3d_

function median_1d_(w) result(m)

    TT, dimension(:), intent(in) :: w
    TT :: m

    TT, allocatable, dimension(:) :: wt
    integer(kind=8) :: n

    n = size(w, kind=8)
    !        wt = w !call alloc_array(wt, [1, n], source=w)

    wt = sort(w) !wt)
    if (mod(n, 2) == 0) then
        m = 0.5*(wt(n/2) + wt(n/2 + 1))
    else
        m = wt((n + 1)/2)
    end if

end function median_1d_

function median_2d_(w) result(m)

    TT, dimension(:, :), intent(in) :: w
    TT :: m

    TT, allocatable, dimension(:) :: wt
    integer(kind=8) :: n

    n = size(w, kind=8) !size(w, 1)*size(w, 2)
    !        wt = flatten(w) !call alloc_array(wt, [1, n], source=reshape(w, [n]))

    wt = sort(flatten(w))
    if (mod(n, 2) == 0) then
        m = 0.5*(wt(n/2) + wt(n/2 + 1))
    else
        m = wt((n + 1)/2)
    end if

end function median_2d_

function median_3d_(w) result(m)

    TT, dimension(:, :, :), intent(in) :: w
    TT :: m

    TT, allocatable, dimension(:) :: wt
    integer(kind=8) :: n

    n = size(w, kind=8) ! size(w, 1)*size(w, 2)*size(w, 3)
    !        call alloc_array(wt, [1, n], source=reshape(w, [n]))

    wt = sort(flatten(w)) !wt)
    if (mod(n, 2) == 0) then
        m = 0.5*(wt(n/2) + wt(n/2 + 1))
    else
        m = wt((n + 1)/2)
    end if

end function median_3d_

!
!> 1D discrete crosscorrelation
!
subroutine discrete_xcorr_1d_(x, y, z, kx, ky, kz)

    TT, dimension(:), intent(in) :: x, y
    TT, dimension(:), intent(inout) :: z
    integer, intent(in) :: kx, ky, kz

    integer :: lx, ly, lz
    integer :: kxr, i
    TT, allocatable, dimension(:) :: xr

    lx = size(x, 1)
    ly = size(y, 1)
    lz = size(z, 1)

    allocate (xr(1:lx))
    do i = 1, lx
        xr(i) = x(lx - i + 1)
    end do

    kxr = 1 - (kx - 1) - lx

    call convd(xr, y, z, &
        kxr, &
        ky, &
        kz)

    deallocate (xr)

end subroutine discrete_xcorr_1d_

!
!> 2D discrete crosscorrelation
!
subroutine discrete_xcorr_2d_(x, y, z, &
        kx1, kx2, &
        ky1, ky2, &
        kz1, kz2)

    TT, dimension(:, :), intent(in) :: x, y
    TT, dimension(:, :), intent(inout) :: z
    integer, intent(in) :: kx1, ky1, kz1
    integer, intent(in) :: kx2, ky2, kz2

    integer :: lx1, ly1, lz1
    integer :: lx2, ly2, lz2
    TT, allocatable, dimension(:, :) :: xr
    integer :: i, j, kx1r, kx2r

    lx1 = size(x, 1)
    lx2 = size(x, 2)
    ly1 = size(y, 1)
    ly2 = size(y, 2)
    lz1 = size(z, 1)
    lz2 = size(z, 2)

    allocate (xr(1:lx1, 1:lx2))
    do j = 1, lx2
        do i = 1, lx1
            xr(i, j) = x(lx1 - i + 1, lx2 - j + 1)
        end do
    end do

    kx1r = 1 - (kx1 - 1) - lx1
    kx2r = 1 - (kx2 - 1) - lx2

    call convd(xr, y, z, &
        kx1r, kx2r, &
        ky1, ky2, &
        kz1, kz2)

    deallocate (xr)

end subroutine discrete_xcorr_2d_

!
!> 3D discrete crosscorrelation
!
subroutine discrete_xcorr_3d_(x, y, z, &
        kx1, kx2, kx3, &
        ky1, ky2, ky3, &
        kz1, kz2, kz3)

    integer, intent(in) :: kx1, ky1, kz1
    integer, intent(in) :: kx2, ky2, kz2
    integer, intent(in) :: kx3, ky3, kz3
    TT, dimension(:, :, :), intent(in) :: x, y
    TT, dimension(:, :, :), intent(inout) :: z

    integer :: lx1, ly1, lz1
    integer :: lx2, ly2, lz2
    integer :: lx3, ly3, lz3
    TT, allocatable, dimension(:, :, :) :: xr
    integer :: i, j, k, kx1r, kx2r, kx3r

    lx1 = size(x, 1)
    lx2 = size(x, 2)
    lx3 = size(x, 3)
    ly1 = size(y, 1)
    ly2 = size(y, 2)
    ly3 = size(y, 3)
    lz1 = size(z, 1)
    lz2 = size(z, 2)
    lz3 = size(z, 3)

    allocate (xr(1:lx1, 1:lx2, 1:lx3))
    do k = 1, lx3
        do j = 1, lx2
            do i = 1, lx1
                xr(i, j, k) = x(lx1 - i + 1, lx2 - j + 1, lx3 - k + 1)
            end do
        end do
    end do

    kx1r = 1 - (kx1 - 1) - lx1
    kx2r = 1 - (kx2 - 1) - lx2
    kx3r = 1 - (kx3 - 1) - lx3

    call convd(xr, y, z, &
        kx1r, kx2r, kx3r, &
        ky1, ky2, ky3, &
        kz1, kz2, kz3)

    deallocate (xr)

end subroutine discrete_xcorr_3d_

!
!> 1D FFT-based cross-correlation
!
function xcorr_1d_(x, y, maxlag) result(z)

    TT, dimension(:), intent(in) :: x, y
    integer, intent(in), optional :: maxlag

    integer :: nx, ny, n, nlag
    TT, allocatable, dimension(:) :: z

    ! Dimensions
    nx = size(x)
    ny = size(y)
    if (present(maxlag)) then
        nlag = min(maxlag, max(nx, ny))
    else
        nlag = max(nx, ny) - 1
    end if

    ! Pad to next power 235
    n = next_power_235(max(nx, ny) + nlag)

    ! Do cross-correlation using FFT
    ! What FFT does is circular cross-correlation, therefore the
    ! result must be shifted and selected
    !
    ! Cross-correlation has different definitions in different literature
    ! The definition adopted here is
    !   C_xy(m) = sum [x(n + m) * y^*(n)]
    ! Some literature defines
    !   C_xy(m) = sum [x^*(n) * y(n + m)]
    ! which gives reversed result compared with the first definition
    ! The first definition is adopted in MATLAB and numpy, and is therefore
    ! adopted here as well.
    !
    allocate (z(1:n))
    z = fftshift(ifft( &
        fft(pad(x, [0, n - nx], ['const', 'const']))* &
        conjg(fft(pad(y, [0, n - ny], ['const', 'const']))), real=.true.))
    call alloc_array(z, [-nlag, nlag], &
        source=z(nint((n + 1)/2.0) - nlag:nint((n + 1)/2.0) + nlag))

end function xcorr_1d_

!
!> 2D FFT-based cross-correlation
!
function xcorr_2d_(x, y, maxlag) result(z)

    TT, dimension(:, :), intent(in) :: x, y
    integer, dimension(:), intent(in), optional :: maxlag

    integer :: nx1, ny1, nx2, ny2, n1, n2, nlag1, nlag2
    TT, allocatable, dimension(:, :) :: z

    ! Dimensions
    nx1 = size(x, 1)
    nx2 = size(x, 2)
    ny1 = size(y, 1)
    ny2 = size(y, 2)
    if (present(maxlag)) then
        nlag1 = min(maxlag(1), max(nx1, ny1))
        nlag2 = min(maxlag(2), max(nx2, ny2))
    else
        nlag1 = max(nx1, ny1) - 1
        nlag2 = max(nx2, ny2) - 1
    end if

    ! Pad to next power 235
    n1 = next_power_235(max(nx1, ny1) + nlag1)
    n2 = next_power_235(max(nx2, ny2) + nlag2)

    ! Do cross-correlation using FFT
    ! What FFT does is circular cross-correlation, therefore the
    ! result must be shifted and selected
    allocate (z(1:n1, 1:n2))
    z = fftshift(ifft( &
        fft(pad(x, [0, n1 - nx1, 0, n2 - nx2], ['const', 'const', 'const', 'const']))* &
        conjg(fft(pad(y, [0, n1 - ny1, 0, n2 - ny2], ['const', 'const', 'const', 'const']))), real=.true.))
    call alloc_array(z, [-nlag1, nlag1, -nlag2, nlag2], &
        source=z(nint((n1 + 1)/2.0) - nlag1:nint((n1 + 1)/2.0) + nlag1, &
        nint((n2 + 1)/2.0) - nlag2:nint((n2 + 1)/2.0) + nlag2))

end function xcorr_2d_

!
!> 3D FFT-based cross-correlation
!
function xcorr_3d_(x, y, maxlag) result(z)

    TT, dimension(:, :, :), intent(in) :: x, y
    integer, dimension(:), intent(in), optional :: maxlag

    integer :: nx1, ny1, nx2, ny2, nx3, ny3, n1, n2, n3, nlag1, nlag2, nlag3
    TT, allocatable, dimension(:, :, :) :: z

    ! Dimensions
    nx1 = size(x, 1)
    nx2 = size(x, 2)
    nx3 = size(x, 3)
    ny1 = size(y, 1)
    ny2 = size(y, 2)
    ny3 = size(y, 3)
    if (present(maxlag)) then
        nlag1 = min(maxlag(1), max(nx1, ny1))
        nlag2 = min(maxlag(2), max(nx2, ny2))
        nlag3 = min(maxlag(3), max(nx3, ny3))
    else
        nlag1 = max(nx1, ny1) - 1
        nlag2 = max(nx2, ny2) - 1
        nlag3 = max(nx3, ny3) - 1
    end if

    ! Pad to next power 235
    n1 = next_power_235(max(nx1, ny1) + nlag1)
    n2 = next_power_235(max(nx2, ny2) + nlag2)
    n3 = next_power_235(max(nx3, ny3) + nlag3)

    ! Do cross-correlation using FFT
    ! What FFT does is circular cross-correlation, therefore the
    ! result must be shifted and selected
    allocate (z(1:n1, 1:n2, 1:n3))
    z = fftshift(ifft( &
        fft(pad(x, [0, n1 - nx1, 0, n2 - nx2, 0, n3 - nx3], &
        ['const', 'const', 'const', 'const', 'const', 'const']))* &
        conjg(fft(pad(y, [0, n1 - ny1, 0, n2 - ny2, 0, n3 - ny3], &
        ['const', 'const', 'const', 'const', 'const', 'const']))), real=.true.))
    call alloc_array(z, [-nlag1, nlag1, -nlag2, nlag2, -nlag3, nlag3], &
        source=z( &
        nint((n1 + 1)/2.0) - nlag1:nint((n1 + 1)/2.0) + nlag1, &
        nint((n2 + 1)/2.0) - nlag2:nint((n2 + 1)/2.0) + nlag2, &
        nint((n3 + 1)/2.0) - nlag3:nint((n3 + 1)/2.0) + nlag3))

end function xcorr_3d_

!
!> 1D discrete autocorrelation
!
subroutine discrete_acorr_1d_(x, z, kx, kz)

    TT, dimension(:), intent(in) :: x
    TT, dimension(:), intent(inout) :: z
    integer, intent(in) :: kx, kz

    integer :: lx, lz
    integer :: kxr, i
    TT, allocatable, dimension(:) :: xr

    lx = size(x, 1)
    lz = size(z, 1)

    allocate (xr(1:lx))
    do i = 1, lx
        xr(i) = x(lx - i + 1)
    end do

    kxr = 1 - (kx - 1) - lx

    call convd(xr, x, z, &
        kxr, &
        kx, &
        kz)

    deallocate (xr)

end subroutine discrete_acorr_1d_

!
!> 2D discrete autocorrelation
!
subroutine discrete_acorr_2d_(x, z, &
        kx1, kx2, &
        kz1, kz2)

    TT, dimension(:, :), intent(in) :: x
    TT, dimension(:, :), intent(inout) :: z
    integer, intent(in) :: kx1, kz1
    integer, intent(in) :: kx2, kz2

    integer :: lx1, lz1
    integer :: lx2, lz2
    TT, allocatable, dimension(:, :) :: xr
    integer :: i, j, kx1r, kx2r

    lx1 = size(x, 1)
    lx2 = size(x, 2)
    lz1 = size(z, 1)
    lz2 = size(z, 2)

    allocate (xr(1:lx1, 1:lx2))
    do j = 1, lx2
        do i = 1, lx1
            xr(i, j) = x(lx1 - i + 1, lx2 - j + 1)
        end do
    end do

    kx1r = 1 - (kx1 - 1) - lx1
    kx2r = 1 - (kx2 - 1) - lx2

    call convd(xr, x, z, &
        kx1r, kx2r, &
        kx1, kx2, &
        kz1, kz2)

    deallocate (xr)

end subroutine discrete_acorr_2d_

!
!> 3D discrete autocorrelation
!
subroutine discrete_acorr_3d_(x, z, &
        kx1, kx2, kx3, &
        kz1, kz2, kz3)

    integer, intent(in) :: kx1, kz1
    integer, intent(in) :: kx2, kz2
    integer, intent(in) :: kx3, kz3
    TT, dimension(:, :, :), intent(in) :: x
    TT, dimension(:, :, :), intent(inout) :: z

    integer :: lx1, lz1
    integer :: lx2, lz2
    integer :: lx3, lz3
    TT, allocatable, dimension(:, :, :) :: xr
    integer :: i, j, k, kx1r, kx2r, kx3r

    lx1 = size(x, 1)
    lx2 = size(x, 2)
    lx3 = size(x, 3)
    lz1 = size(z, 1)
    lz2 = size(z, 2)
    lz3 = size(z, 3)

    allocate (xr(1:lx1, 1:lx2, 1:lx3))
    do k = 1, lx3
        do j = 1, lx2
            do i = 1, lx1
                xr(i, j, k) = x(lx1 - i + 1, lx2 - j + 1, lx3 - k + 1)
            end do
        end do
    end do

    kx1r = 1 - (kx1 - 1) - lx1
    kx2r = 1 - (kx2 - 1) - lx2
    kx3r = 1 - (kx3 - 1) - lx3

    call convd(xr, x, z, &
        kx1r, kx2r, kx3r, &
        kx1, kx2, kx3, &
        kz1, kz2, kz3)

    deallocate (xr)

end subroutine discrete_acorr_3d_

!
!> 1D FFT-based auto-correlation
!
function acorr_1d_(x, maxlag) result(z)

    TT, dimension(:), intent(in) :: x
    integer, intent(in), optional :: maxlag

    integer :: nx, n, nlag
    TT, allocatable, dimension(:) :: z

    ! Dimensions
    nx = size(x)
    if (present(maxlag)) then
        nlag = min(maxlag, nx)
    else
        nlag = nx - 1
    end if

    ! Pad to next power 235
    n = next_power_235(nx + nlag)

    ! Do cross-correlation using FFT
    ! What FFT does is circular cross-correlation, therefore the
    ! result must be shifted and selected
    allocate (z(1:n))
    z = fftshift(ifft( &
        TTT(fft(pad(x, [0, n - nx], ['const', 'const']))**2), real=.true.))
    call alloc_array(z, [-nlag, nlag], &
        source=z(nint((n + 1)/2.0) - nlag:nint((n + 1)/2.0) + nlag))

end function acorr_1d_

!
!> 2D FFT-based auto-correlation
!
function acorr_2d_(x, maxlag) result(z)

    TT, dimension(:, :), intent(in) :: x
    integer, dimension(:), intent(in), optional :: maxlag

    integer :: nx1, nx2, n1, n2, nlag1, nlag2
    TT, allocatable, dimension(:, :) :: z

    ! Dimensions
    nx1 = size(x, 1)
    nx2 = size(x, 2)
    if (present(maxlag)) then
        nlag1 = min(maxlag(1), nx1)
        nlag2 = min(maxlag(2), nx2)
    else
        nlag1 = nx1 - 1
        nlag2 = nx2 - 1
    end if

    ! Pad to next power 235
    n1 = next_power_235(nx1 + nlag1)
    n2 = next_power_235(nx2 + nlag2)

    ! Do cross-correlation using FFT
    ! What FFT does is circular cross-correlation, therefore the
    ! result must be shifted and selected
    allocate (z(1:n1, 1:n2))
    z = fftshift(ifft( &
        TTT(fft(pad(x, [0, n1 - nx1, 0, n2 - nx2], &
        ['const', 'const', 'const', 'const']))**2), real=.true.))
    call alloc_array(z, [-nlag1, nlag1, -nlag2, nlag2], &
        source=z(nint((n1 + 1)/2.0) - nlag1:nint((n1 + 1)/2.0) + nlag1, &
        nint((n2 + 1)/2.0) - nlag2:nint((n2 + 1)/2.0) + nlag2))

end function acorr_2d_

!
!> 3D FFT-based auto-correlation
!
function acorr_3d_(x, maxlag) result(z)

    TT, dimension(:, :, :), intent(in) :: x
    integer, dimension(:), intent(in), optional :: maxlag

    integer :: nx1, nx2, nx3, n1, n2, n3, nlag1, nlag2, nlag3
    TT, allocatable, dimension(:, :, :) :: z

    ! Dimensions
    nx1 = size(x, 1)
    nx2 = size(x, 2)
    nx3 = size(x, 3)
    if (present(maxlag)) then
        nlag1 = min(maxlag(1), nx1)
        nlag2 = min(maxlag(2), nx2)
        nlag3 = min(maxlag(3), nx3)
    else
        nlag1 = nx1 - 1
        nlag2 = nx2 - 1
        nlag3 = nx3 - 1
    end if

    ! Pad to next power 235
    n1 = next_power_235(nx1 + nlag1)
    n2 = next_power_235(nx2 + nlag2)
    n3 = next_power_235(nx3 + nlag3)

    ! Do cross-correlation using FFT
    ! What FFT does is circular cross-correlation, therefore the
    ! result must be shifted and selected
    allocate (z(1:n1, 1:n2, 1:n3))
    z = fftshift(ifft( &
        TTT(fft(pad(x, [0, n1 - nx1, 0, n2 - nx2, 0, n3 - nx3], &
        ['const', 'const', 'const', 'const', 'const', 'const']))**2), real=.true.))
    call alloc_array(z, [-nlag1, nlag1, -nlag2, nlag2, -nlag3, nlag3], &
        source=z( &
        nint((n1 + 1)/2.0) - nlag1:nint((n1 + 1)/2.0) + nlag1, &
        nint((n2 + 1)/2.0) - nlag2:nint((n2 + 1)/2.0) + nlag2, &
        nint((n3 + 1)/2.0) - nlag3:nint((n3 + 1)/2.0) + nlag3))

end function acorr_3d_

!
!> 1D discrete cross-correlation
!
subroutine local_xcorr_1d_(tr1, tr2, wt, ccoef, tshift)

    TT, dimension(:), intent(in) :: tr1, tr2
    TT, dimension(:), intent(inout) :: ccoef, tshift
    integer, intent(in) :: wt

    integer :: nt, i, ns
    TT, allocatable, dimension(:) :: ptr1, ptr2, pw1, pw2, pc
    TT :: norm

    ! number of samples
    nt = size(tr1, 1)

    if (maxval(abs(tr1)) == 0 .or. maxval(abs(tr2)) == 0) then
        ccoef = 1.0
        tshift = 0.0
    else
        ns = floor(wt/3.0)
        call alloc_array(ptr1, [-wt + 1, nt + wt])
        call alloc_array(ptr2, [-wt + 1, nt + wt])
        ptr1(1:nt) = tr1
        ptr2(1:nt) = tr2
        call alloc_array(pw1, [-wt, wt])
        call alloc_array(pw1, [-wt, wt])
        call alloc_array(pc, [-ns, ns])
        do i = 1, nt
            pw1 = ptr1(i - wt:i + wt)
            pw2 = ptr2(i - wt:i + wt)
            call discrete_xcorr_1d_(pw1, pw2, pc, 1, 1, -ns)
            norm = norm2(pw1)*norm2(pw2)
            if (norm == 0) then
                ccoef(i) = 0.0
                tshift(i) = 0.0
            else
                ccoef(i) = maxval(pc)/norm
                tshift(i) = maxloc(pc, 1) - ns - 1.0
            end if
        end do

    end if

end subroutine local_xcorr_1d_

!
!> Zero-lag cross-correlation coefficient
!
function xcorr_coef_1d_(a, b) result(r)

    TT, dimension(:), intent(in)  :: a
    TT, dimension(:), intent(in)  :: b
    TT :: r

    TT :: ma, mb

    call assert(size(a) == size(b), 'Error: size(a) /= size(b)')

    ma = mean(a)
    mb = mean(b)

    r = sum((a - ma)*(b - mb))/sqrt(sum((a - ma)**2))/sqrt(sum((b - mb)**2))

end function xcorr_coef_1d_

function xcorr_coef_2d_(a, b) result(r)

    TT, dimension(:, :), intent(in)  :: a
    TT, dimension(:, :), intent(in)  :: b
    TT :: r

    TT :: ma, mb

    call assert(size(a, 1) == size(b, 1), 'Error: size(a, 1) /= size(b, 1)')
    call assert(size(a, 2) == size(b, 2), 'Error: size(a, 2) /= size(b, 2)')

    ma = mean(a)
    mb = mean(b)

    r = sum((a - ma)*(b - mb))/sqrt(sum((a - ma)**2))/sqrt(sum((b - mb)**2))

end function xcorr_coef_2d_

function xcorr_coef_3d_(a, b) result(r)

    TT, dimension(:, :, :), intent(in)  :: a
    TT, dimension(:, :, :), intent(in)  :: b
    TT :: r

    TT :: ma, mb

    call assert(size(a, 1) == size(b, 1), 'Error: size(a, 1) /= size(b, 1)')
    call assert(size(a, 2) == size(b, 2), 'Error: size(a, 2) /= size(b, 2)')
    call assert(size(a, 3) == size(b, 3), 'Error: size(a, 3) /= size(b, 3)')

    ma = mean(a)
    mb = mean(b)

    r = sum((a - ma)*(b - mb))/sqrt(sum((a - ma)**2))/sqrt(sum((b - mb)**2))

end function xcorr_coef_3d_

#undef T
#undef TT
#undef TTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef kernel_triangular_
#undef kernel_epanechnikov_
#undef kernel_bisquare_
#undef kernel_trisquare_
#undef kernel_tricube_
#undef kernel_cosine_
#undef covariance_matrix_1d_
#undef covariance_matrix_2d_
#undef covariance_matrix_3d_
#undef row_covariance_matrix_2d_
#undef cross_covariance_matrix_1d_
#undef cross_covariance_matrix_2d_
#undef cross_covariance_matrix_3d_
#undef row_cross_covariance_matrix_2d_
#undef variance_1d_
#undef variance_2d_
#undef variance_3d_
#undef gaussian_1d_
#undef gaussian_2d_
#undef gaussian_3d_
#undef covariance_1d_
#undef covariance_2d_
#undef covariance_3d_
#undef standard_deviation_1d_
#undef standard_deviation_2d_
#undef standard_deviation_3d_
#undef histogram_1d_
#undef histogram_2d_
#undef histogram_3d_
#undef plot_histogram_1d_
#undef plot_histogram_2d_
#undef plot_histogram_3d_
#undef median_1d_
#undef median_2d_
#undef median_3d_
#undef discrete_xcorr_1d_
#undef discrete_xcorr_2d_
#undef discrete_xcorr_3d_
#undef xcorr_1d_
#undef xcorr_2d_
#undef xcorr_3d_
#undef discrete_acorr_1d_
#undef discrete_acorr_2d_
#undef discrete_acorr_3d_
#undef acorr_1d_
#undef acorr_2d_
#undef acorr_3d_
#undef local_xcorr_1d_
#undef xcorr_coef_1d_
#undef xcorr_coef_2d_
#undef xcorr_coef_3d_
