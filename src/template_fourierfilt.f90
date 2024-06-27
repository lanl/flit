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

#define fourier_filter_1d_      CONCAT(fourier_filter_1d, T)
#define fourier_filt_1d_      CONCAT(fourier_filt_1d, T)
#define fourier_filt_2d_      CONCAT(fourier_filt_2d, T)
#define fourier_filt_3d_      CONCAT(fourier_filt_3d, T)
#define fourier_sharpen_1d_      CONCAT(fourier_sharpen_1d, T)
#define fourier_sharpen_2d_      CONCAT(fourier_sharpen_2d, T)
#define fourier_sharpen_3d_      CONCAT(fourier_sharpen_3d, T)
#define fourier_smooth_1d_      CONCAT(fourier_smooth_1d, T)
#define fourier_smooth_2d_      CONCAT(fourier_smooth_2d, T)
#define fourier_smooth_3d_      CONCAT(fourier_smooth_3d, T)

!
!> Define a frequency-domain filter using freqs and amps
!
function fourier_filter_1d_(f, a, nt, dt, method, alpha) result(w)

    TT, dimension(:), intent(in) :: f, a
    integer, intent(in) :: nt
    TT, intent(in) :: dt
    character(len=*), optional, intent(in) :: method
    TT, optional, intent(in) :: alpha
    TT, allocatable, dimension(:) :: w

    integer :: n, i, nf, f1, f2
    TT :: df, a1, a2
    integer, allocatable, dimension(:) :: findex
    character(len=24) :: window_method
    real :: window_alpha
    TT, allocatable, dimension(:) :: window

    call assert(size(f) == size(a), ' <fourier_filter_1d_> Error: size(f) /= size(a). ')
    call assert(all(f >= 0), ' <fourier_filter_1d_> Error: f must >= 0. ')
    call assert(all(a >= 0), ' <fourier_filter_1d_> Error: a must > 0. ')

    n = size(f)
    nf = nint((nt + 1.0)/2.0)
    allocate (w(1:nt))
    df = 1.0/dt/nt
    findex = min(max(nint(f/df) + 1, 1), nf)

    if (present(method)) then
        window_method = method
    else
        window_method = 'hann'
    end if

    if (present(alpha)) then
        window_alpha = alpha
    else
        window_alpha = 4.0
    end if

    w(1:findex(1)) = a(1)

    do i = 1, n - 1
        f1 = findex(i)
        f2 = findex(i + 1)
        a1 = a(i)
        a2 = a(i + 1)
        if (a2 > a1) then
            call alloc_array(window, [f1, f2], &
                source=TTT(taper_window(f2 - f1 + 1, [f2 - f1 + 1, 0], &
                [window_method, ''], [window_alpha, 0.0])))
            w(f1:f2) = rescale(window, [a1, a2])
        else
            call alloc_array(window, [f1, f2], &
                source=TTT(taper_window(f2 - f1 + 1, [0, f2 - f1 + 1], &
                ['', window_method], [0.0, window_alpha])))
            w(f1:f2) = rescale(window, [a2, a1])
        end if
    end do

    w(f2 + 1:nf) = a(n)
    if (mod(nt, 2) == 0) then
        w(nf:nt) = w(nf:2:-1)
    else
        w(nf + 1:nt) = w(nf:2:-1)
    end if

end function fourier_filter_1d_

!
!> Zero-phase, frequency-domain, dimension-separable filtering
!
function fourier_filt_1d_(w, dt, freqs, amps, method, alpha) result(wt)

    TT, dimension(:), intent(in) :: w
    TT, intent(in) :: dt
    TT, dimension(:), intent(in) :: freqs, amps
    character(len=*), optional :: method
    TT, optional :: alpha
    TT, allocatable, dimension(:) :: wt

    integer :: nt, nnt
    character(len=24) :: window_method
    TT :: window_alpha

    if (present(method)) then
        window_method = method
    else
        window_method = 'hann'
    end if

    if (present(alpha)) then
        window_alpha = alpha
    else
        window_alpha = 4.0
    end if

    nt = size(w)
    nnt = next_power_235(nint(1.5*nt))
    wt = ifft(fft([w, TTT(zeros(nnt - nt))]) &
        *fourier_filter(freqs, amps, nnt, dt, window_method, window_alpha), &
        real=.true.)
    wt = wt(1:nt)

end function fourier_filt_1d_

function fourier_filt_2d_(w, d1, freqs1, amps1, d2, freqs2, amps2, &
        method, alpha) result(wt)

    TT, dimension(:, :), intent(inout) :: w
    TT, intent(in) :: d1, d2
    TT, dimension(:), intent(in) :: freqs1, amps1, freqs2, amps2
    character(len=*), optional :: method
    TT, optional :: alpha
    TT, allocatable, dimension(:, :) :: wt

    integer :: i, j
    integer :: n1, n2
    character(len=24) :: window_method
    TT :: window_alpha

    if (present(method)) then
        window_method = method
    else
        window_method = 'hann'
    end if

    if (present(alpha)) then
        window_alpha = alpha
    else
        window_alpha = 4.0
    end if

    n1 = size(w, 1)
    n2 = size(w, 2)
    allocate (wt(1:n1, 1:n2))

    !$omp parallel do private(i, j)
    do j = 1, n2
        wt(:, j) = fourier_filt_1d_(w(:, j), d1, freqs1, amps1, &
            window_method, window_alpha)
    end do
    !$omp end parallel do

    !$omp parallel do private(i, j)
    do i = 1, n1
        wt(i, :) = fourier_filt_1d_(wt(i, :), d2, freqs2, amps2, &
            window_method, window_alpha)
    end do
    !$omp end parallel do

end function fourier_filt_2d_

function fourier_filt_3d_(w, d1, freqs1, amps1, &
        d2, freqs2, amps2, d3, freqs3, amps3, method, alpha) result(wt)

    TT, dimension(:, :, :), intent(inout) :: w
    TT, intent(in) :: d1, d2, d3
    TT, dimension(:), intent(in) :: freqs1, amps1, freqs2, amps2, freqs3, amps3
    character(len=*), optional :: method
    TT, optional :: alpha
    TT, allocatable, dimension(:, :, :) :: wt

    integer :: i, j, k
    integer :: n1, n2, n3
    character(len=24) :: window_method
    TT :: window_alpha

    if (present(method)) then
        window_method = method
    else
        window_method = 'hann'
    end if

    if (present(alpha)) then
        window_alpha = alpha
    else
        window_alpha = 4.0
    end if

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)
    allocate (wt(1:n1, 1:n2, 1:n3))

    !$omp parallel do private(i, j, k) collapse(2)
    do k = 1, n3
        do j = 1, n2
            wt(:, j, k) = fourier_filt_1d_(w(:, j, k), d1, freqs1, amps1, &
                window_method, window_alpha)
        end do
    end do
    !$omp end parallel do

    !$omp parallel do private(i, j, k) collapse(2)
    do k = 1, n3
        do i = 1, n1
            wt(i, :, k) = fourier_filt_1d_(wt(i, :, k), d2, freqs2, amps2, &
                window_method, window_alpha)
        end do
    end do
    !$omp end parallel do

    !$omp parallel do private(i, j, k) collapse(2)
    do j = 1, n2
        do i = 1, n1
            wt(i, j, :) = fourier_filt_1d_(wt(i, j, :), d3, freqs3, amps3, &
                window_method, window_alpha)
        end do
    end do
    !$omp end parallel do

end function fourier_filt_3d_

!
!> High-frequency emphasis filtering (I name it sharpen)
!
function fourier_sharpen_1d_(f, k1, k2, sigma) result(g)

    TT, dimension(:), intent(in) :: f
    TT, intent(in) :: k1, k2, sigma
    TT, allocatable, dimension(:) :: g

    TT, allocatable, dimension(:) :: h
    integer :: i
    integer :: n1

    n1 = size(f, 1)

    h = zeros(n1)
    !$omp parallel do private(i)
    do i = 1, n1
        h(i) = 1.0 - exp(-(i - n1/2.0)**2/(2*sigma**2))
    end do
    !$omp end parallel do

    g = ifft(fftshift((k1 + k2*h)*fftshift(fft(f))), real=.true.)

end function fourier_sharpen_1d_

function fourier_sharpen_2d_(f, k1, k2, sigma) result(g)

    TT, dimension(:, :), intent(in) :: f
    TT, intent(in) :: k1, k2
    TT, dimension(1:2), intent(in) :: sigma
    TT, allocatable, dimension(:, :) :: g

    TT, allocatable, dimension(:, :) :: h
    integer :: i, j
    integer :: n1, n2

    n1 = size(f, 1)
    n2 = size(f, 2)

    h = zeros(n1, n2)
    !$omp parallel do private(i, j)
    do j = 1, n2
        do i = 1, n1
            h(i, j) = 1.0 - exp(-(i - n1/2.0)**2/(2*sigma(1)**2) - (j - n2/2.0)**2/(2*sigma(2)**2))
        end do
    end do
    !$omp end parallel do

    g = ifft(fftshift((k1 + k2*h)*fftshift(fft(f))), real=.true.)

end function fourier_sharpen_2d_

function fourier_sharpen_3d_(f, k1, k2, sigma) result(g)

    TT, dimension(:, :, :), intent(in) :: f
    TT, intent(in) :: k1, k2
    TT, dimension(1:3), intent(in) :: sigma
    TT, allocatable, dimension(:, :, :) :: g

    TT, allocatable, dimension(:, :, :) :: h
    integer :: i, j, k
    integer :: n1, n2, n3

    n1 = size(f, 1)
    n2 = size(f, 2)
    n3 = size(f, 3)

    h = zeros(n1, n2, n3)
    !$omp parallel do private(i, j, k)
    do k = 1, n3
        do j = 1, n2
            do i = 1, n1
                h(i, j, k) = 1.0 - exp(-(i - n1/2.0)**2/(2*sigma(1)**2) - (j - n2/2.0)**2/(2*sigma(2)**2) &
                    - (k - n3/2.0)**2/(2*sigma(3)**2))
            end do
        end do
    end do
    !$omp end parallel do

    g = ifft(fftshift((k1 + k2*h)*fftshift(fft(f))), real=.true.)

end function fourier_sharpen_3d_

!
!> Low-frequency emphasis filtering (I name it smooth)
!
function fourier_smooth_1d_(f, sigma) result(g)

    TT, dimension(:), intent(in) :: f
    TT, intent(in) :: sigma
    TT, allocatable, dimension(:) :: g

    TT, allocatable, dimension(:) :: h
    integer :: i
    integer :: n1

    n1 = size(f)

    h = zeros(n1)
    !$omp parallel do private(i)
    do i = 1, n1
        h(i) = exp(-(i - n1/2.0)**2/(2*sigma**2))
    end do
    !$omp end parallel do

    g = ifft(fftshift(h*fftshift(fft(f))), real=.true.)

end function fourier_smooth_1d_

function fourier_smooth_2d_(f, sigma) result(g)

    TT, dimension(:, :), intent(in) :: f
    TT, dimension(1:2), intent(in) :: sigma
    TT, allocatable, dimension(:, :) :: g

    TT, allocatable, dimension(:, :) :: h
    integer :: i, j
    integer :: n1, n2

    n1 = size(f, 1)
    n2 = size(f, 2)

    h = zeros(n1, n2)
    !$omp parallel do private(i, j)
    do j = 1, n2
        do i = 1, n1
            h(i, j) = exp(-(i - n1/2.0)**2/(2*sigma(1)**2) - (j - n2/2.0)**2/(2*sigma(2)**2))
        end do
    end do
    !$omp end parallel do

    g = ifft(fftshift(h*fftshift(fft(f))), real=.true.)

end function fourier_smooth_2d_

function fourier_smooth_3d_(f, sigma) result(g)

    TT, dimension(:, :, :), intent(in) :: f
    TT, dimension(1:3), intent(in) :: sigma
    TT, allocatable, dimension(:, :, :) :: g

    TT, allocatable, dimension(:, :, :) :: h
    integer :: i, j, k
    integer :: n1, n2, n3

    n1 = size(f, 1)
    n2 = size(f, 2)
    n3 = size(f, 3)

    h = zeros(n1, n2, n3)
    !$omp parallel do private(i, j, k)
    do k = 1, n3
        do j = 1, n2
            do i = 1, n1
                h(i, j, k) = exp(-(i - n1/2.0)**2/(2*sigma(1)**2) - (j - n2/2.0)**2/(2*sigma(2)**2) &
                    - (k - n3/2.0)**2/(2*sigma(3)**2))
            end do
        end do
    end do
    !$omp end parallel do

    g = ifft(fftshift(h*fftshift(fft(f))), real=.true.)

end function fourier_smooth_3d_

#undef T
#undef TT
#undef TTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef fourier_filter_1d_
#undef fourier_filt_1d_
#undef fourier_filt_2d_
#undef fourier_filt_3d_
#undef fourier_sharpen_1d_
#undef fourier_sharpen_2d_
#undef fourier_sharpen_3d_
#undef fourier_smooth_1d_
#undef fourier_smooth_2d_
#undef fourier_smooth_3d_

