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

module libflit_spectrum

    use libflit_array
    use libflit_constants
    use libflit_transform
    use libflit_utility
    use libflit_fourierfilt
    use libflit_taper
    use libflit_statistics
    use libflit_calculus
    use libflit_interp
    use iso_fortran_env
    use libflit_string
    use libflit_io
    use libflit_date_time
    use libflit_array_operation

    use, intrinsic :: iso_c_binding

contains

    !
    !> Compute the instaneous phase of a 1D signal based on the method developed by
    !>
    !> E. Poggiagliolmi, A. Vesnaver, 2014, 
    !> Instantaneous phase and frequency derived without user-defined parameters,
    !> Geophysical Journal International, doi: 10.1093/gji/ggu352
    !
    function instant_phase(w) result(wr)

        real, dimension(:), intent(in) :: w

        complex, allocatable, dimension(:) :: v
        real, allocatable, dimension(:) :: wr, env
        integer :: nt

        nt = size(w)

        call alloc_array(v, [1, nt])
        call alloc_array(env, [1, nt])

        v = cmplx(w, hilbert(w))
        env = abs(v)
        v = v/env
        where (env < 1.0e-3*maxval(env))
            v = 0
        end where

        wr = real(integ(conjg(v)*deriv(v)*(-const_i)))

    end function

    !
    !> Generalized s-transform using adaptive Gaussian window
    !> lambda = 1, p = 1 --> short-time Fourier transform
    !
    function gst(w, lambda, p, dt, fmin, fmax, df) result(stran)

        real, dimension(:), intent(in) :: w
        real, intent(in) :: lambda, p
        real, intent(in), optional :: dt, fmin, fmax, df
        complex, allocatable, dimension(:, :) :: stran

        complex, allocatable, dimension(:) :: ww, e
        integer :: nt, nf, iw, i, ifmin, ifmax, nf0
        real :: gst_dt, gst_fmin, gst_fmax, gst_df
        integer :: dif
        real :: default_df

        nt = size(w)
        nf0 = nint(nt/2.0)

        if (present(dt)) then
            gst_dt = dt
        else
            gst_dt = 1.0
        end if

        if (present(fmin)) then
            gst_fmin = fmin
        else
            gst_fmin = 0.0
        end if

        if (present(fmax)) then
            gst_fmax = fmax
        else
            gst_fmax = 0.5/gst_dt
        end if

        default_df = 1.0/gst_dt/nt

        if (present(df)) then
            gst_df = df
        else
            gst_df = default_df
        end if

        ! Determine the integer location of min, max frequencies in the Fourier domain
        ifmin = floor(gst_fmin/default_df) + 1
        ifmax = clip(ceiling(gst_fmax/default_df) + 1, ifmin, int(nt/2.0))
        dif = clip(floor(gst_df/default_df), 1, ifmax - ifmin)
        nf = clip(nint((ifmax - ifmin + 0.0)/dif) + 1, 1, int(nt/2.0))

        stran = zeros(nt, nf)
        ww = zeros(2*nt)
        e = zeros(nt)

        ! Fourier representation of the signal
        ww(1:nt) = fft(w)
        ww(nt + 1:2*nt) = ww(1:nt)

        ! For each frequency, multiply corresponding generalized Gaussian window
        ! 	and do inverse Fourier transform
        ! 	See equation (5) of Liu et al., 2019, 
        ! 	Self-Adaptive Generalized S-Transform and Its Application in Seismic Time-Frequency Analysis, 
        !   IEEE-TGRS, 10.1109/TGRS.2019.2916792
        !$omp parallel do private(iw, i, e)
        do iw = ifmin, ifmax, dif
            do i = 1, nt
                e(i) = exp(-2*const_pi**2*((i - nf0))**2/(lambda*(iw - 1)**p)**2)
                ! An alternative windows function is
                ! e(i) = exp(-2*const_pi**2*((i - nf0))**2/(f0 - a*log(2*f0/((iw - 1)*gst_df + float_tiny) - 1.0))**2)
                ! where a is a hyperparameter, and f0 = 0.5/gst_dt
            end do
            ! fortran's cshift is a reverse version of matlab's circshift
            e = cshift(e, -(nt - nf0 + 1))
            stran(:, (iw - ifmin)/dif + 1) = ifft(ww(iw:iw + nt - 1)*e)
        end do
        !$omp end parallel do

        ! f = 0 is the mean of the signal
        if (ifmin == 1) then
            stran(:, 1) = cmplx(mean(w), 0.0)
        end if

        write (error_unit, *) ' Actual fmin, fmax, nf, df = ', &
            num2str((ifmin - 1)*default_df, '(es)'), ' ', &
            num2str((ifmax - 1)*default_df, '(es)'), ' ', &
            num2str(nf), ' ', num2str(dif*default_df, '(es)')

    end function

    !
    !> Inverse generalized s-transform
    !
    function igst(stran, dt, fmin, fmax, df) result(w)

        complex, dimension(:, :), intent(in) :: stran
        real, intent(in), optional :: dt, fmin, fmax
        real, intent(in), optional :: df
        real, allocatable, dimension(:) :: w

        integer :: nt, nf, ifmin, ifmax
        real :: gst_dt, gst_fmin, gst_fmax, gst_df
        real :: default_df
        integer :: dif
        complex, dimension(:), allocatable :: ww

        nt = size(stran, 1)

        if (present(dt)) then
            gst_dt = dt
        else
            gst_dt = 1.0
        end if

        if (present(fmin)) then
            gst_fmin = fmin
        else
            gst_fmin = 0.0
        end if

        if (present(fmax)) then
            gst_fmax = fmax
        else
            gst_fmax = 0.5/gst_dt
        end if

        default_df = 1.0/gst_dt/nt

        if (present(df)) then
            gst_df = df
        else
            gst_df = default_df
        end if

        ! Determine the integer location of min, max frequencies in the Fourier domain
        ifmin = floor(gst_fmin/default_df) + 1
        ifmax = clip(ceiling(gst_fmax/default_df) + 1, ifmin, int(nt/2.0))
        dif = clip(floor(gst_df/default_df), 1, ifmax - ifmin)
        nf = clip(nint((ifmax - ifmin + 0.0)/dif) + 1, 1, int(nt/2.0))

        if (size(stran, 2) /= nf) then
            write (error_unit, *) ' Error: size(stran, 2) must = '//num2str(nf)
        end if

        ww = zeros(nt)
        if (gst_df == default_df) then
            ! If the given gst has the default frequency sampling interval, then directly sum over t
            ww(ifmin:ifmax) = 2.0*sum(stran, dim=1)
        else
            ! Otherwise, must first upsample to the default frequency sampling interval
            ! from dif*default_df to 1*default_df, and then sum over t
            !$omp parallel do private(i) reduction(+: w)
            do i = 1, nt
                ww(ifmin:ifmax) = ww(ifmin:ifmax) + 2*cmplx( &
                    interp(real(stran(i, :)), nf, dif*1.0, 0.0, ifmax - ifmin + 1, 1.0, 0.0, 'sinc'), &
                    interp(imag(stran(i, :)), nf, dif*1.0, 0.0, ifmax - ifmin + 1, 1.0, 0.0, 'sinc'))
            end do
            !$omp end parallel do
        end if

        w = ifft(ww, real=.true.)

    end function

    !
    !> Gabor transform
    !>
    !> The parameter sigma is like the fraction of the entire time series length (fraction of 1)
    !
    function gabor(w, sigma, dt, fmin, fmax, df) result(stran)

        real, dimension(:), intent(in) :: w
        real, intent(in) :: sigma
        real, intent(in), optional :: dt, fmin, fmax, df
        complex, allocatable, dimension(:, :) :: stran

        complex, allocatable, dimension(:) :: ww, e
        integer :: nt, nf, iw, i, ifmin, ifmax, nf0
        real :: gabor_dt, gabor_fmin, gabor_fmax, gabor_df
        integer :: dif
        real :: default_df

        nt = size(w)
        nf0 = nint(nt/2.0)

        if (present(dt)) then
            gabor_dt = dt
        else
            gabor_dt = 1.0
        end if

        if (present(fmin)) then
            gabor_fmin = fmin
        else
            gabor_fmin = 0.0
        end if

        if (present(fmax)) then
            gabor_fmax = fmax
        else
            gabor_fmax = 0.5/gabor_dt
        end if

        default_df = 1.0/gabor_dt/nt

        if (present(df)) then
            gabor_df = df
        else
            gabor_df = default_df
        end if

        ! Determine the integer location of min, max frequencies in the Fourier domain
        ifmin = floor(gabor_fmin/default_df) + 1
        ifmax = clip(ceiling(gabor_fmax/default_df) + 1, ifmin, int(nt/2.0))
        dif = clip(floor(gabor_df/default_df), 1, ifmax - ifmin)
        nf = clip(nint((ifmax - ifmin + 0.0)/dif) + 1, 1, int(nt/2.0))

        stran = zeros(nt, nf)
        ww = zeros(2*nt)
        e = zeros(nt)

        ! Fourier representation of the signal
        ww(1:nt) = fft(w)
        ww(nt + 1:2*nt) = ww(1:nt)

        ! For each frequency, multiply corresponding generalized Gaussian window
        ! and do inverse Fourier transform
        !$omp parallel do private(iw, i, e)
        do iw = ifmin, ifmax, dif
            ! This is equivalent to GST by setting lambda = 1, f = sigma^2, p = -1/2
            ! This leads to admissibility condition: \int_\infty kernel = 1
            do i = 1, nt
                ! e(i) = exp(-2*const_pi**2*((i - nf0)*gabor_dt)**2/(1.0*(sigma**2)**(-0.5d0))**2)
                e(i) = exp(-2*const_pi**2*((i - nf0))**2*sigma**2)
            end do
            ! fortran's cshift is a reverse version of matlab's circshift
            e = cshift(e, -(nt - nf0 + 1))
            stran(:, (iw - ifmin)/dif + 1) = ifft(ww(iw:iw + nt - 1)*e)
        end do
        !$omp end parallel do

        ! f = 0 is the mean of the signal
        if (ifmin == 1) then
            stran(:, 1) = cmplx(mean(w), 0.0)
        end if

        write (error_unit, *) ' Actual fmin, fmax, nf, df = ', &
            num2str((ifmin - 1)*default_df, '(es)'), ' ', &
            num2str((ifmax - 1)*default_df, '(es)'), ' ', &
            num2str(nf), ' ', num2str(dif*default_df, '(es)')

    end function

    !
    !> Inverse Gabor transform
    !
    function igabor(stran, dt, fmin, fmax, df) result(w)

        complex, dimension(:, :), intent(in) :: stran
        real, intent(in), optional :: dt, fmin, fmax
        real, intent(in), optional :: df
        real, allocatable, dimension(:) :: w

        integer :: nt, nf, ifmin, ifmax
        real :: gabor_dt, gabor_fmin, gabor_fmax, gabor_df
        real :: default_df
        integer :: dif
        complex, dimension(:), allocatable :: ww

        nt = size(stran, 1)

        if (present(dt)) then
            gabor_dt = dt
        else
            gabor_dt = 1.0
        end if

        if (present(fmin)) then
            gabor_fmin = fmin
        else
            gabor_fmin = 0.0
        end if

        if (present(fmax)) then
            gabor_fmax = fmax
        else
            gabor_fmax = 0.5/gabor_dt
        end if

        default_df = 1.0/gabor_dt/nt

        if (present(df)) then
            gabor_df = df
        else
            gabor_df = default_df
        end if

        ! Determine the integer location of min, max frequencies in the Fourier domain
        ifmin = floor(gabor_fmin/default_df) + 1
        ifmax = clip(ceiling(gabor_fmax/default_df) + 1, ifmin, int(nt/2.0))
        dif = clip(floor(gabor_df/default_df), 1, ifmax - ifmin)
        nf = clip(nint((ifmax - ifmin + 0.0)/dif) + 1, 1, int(nt/2.0))

        if (size(stran, 2) /= nf) then
            write (error_unit, *) ' Error: size(stran, 2) must = '//num2str(nf)
        end if

        ww = zeros(nt)
        if (gabor_df == default_df) then
            ! If the given gabor has the default frequency sampling interval, then directly sum over t
            ww(ifmin:ifmax) = 2.0*sum(stran, dim=1)
        else
            ! Otherwise, must first upsample to the default frequency sampling interval
            ! from dif*default_df to 1*default_df, and then sum over t
            !$omp parallel do private(i) reduction(+: w)
            do i = 1, nt
                ww(ifmin:ifmax) = ww(ifmin:ifmax) + 2*cmplx( &
                    interp(real(stran(i, :)), nf, dif*1.0, 0.0, ifmax - ifmin + 1, 1.0, 0.0, 'sinc'), &
                    interp(imag(stran(i, :)), nf, dif*1.0, 0.0, ifmax - ifmin + 1, 1.0, 0.0, 'sinc'))
            end do
            !$omp end parallel do
        end if

        w = ifft(ww, real=.true.)

    end function

    !
    !> Generalized s-transform
    !>
    !> A simple version that assume dt = 1 and the other parameters are default values
    !> derived based on the input sigmal length and dt = 1
    !
    function gst_simple(w, lambda, p) result(stran)

        real, dimension(:), intent(in) :: w
        real, intent(in) :: lambda, p
        complex, allocatable, dimension(:, :) :: stran

        complex, allocatable, dimension(:) :: ww, e
        real, allocatable, dimension(:) :: ff
        integer :: nt, iw, i

        nt = size(w)
        nf0 = nint(nt/2.0)
        stran = zeros(nt, nt)
        ww = zeros(2*nt)
        e = zeros(nt)

        ! Fourier representation of the signal
        ww(1:nt) = fft(w)
        ww(nt + 1:2*nt) = ww(1:nt)

        ff = fft_omega(nt)/(2*const_pi)*nt
        !$omp parallel do private(iw, i, e)
        do iw = 1, nt
            do i = 1, nt
                e(i) = exp(-2*const_pi**2*((i - nf0))**2/(lambda*abs(ff(iw))**p)**2)
            end do
            ! fortran's cshift is a reverse version of matlab's circshift
            e = cshift(e, -(nt - nf0 + 1))
            stran(:, iw) = ifft(ww(iw:iw + nt - 1)*e)
        end do
        !$omp end parallel do

        ! f = 0 is the mean of the signal
        stran(:, 1) = cmplx(mean(w), 0.0)

    end function

    !
    !> Inverse generalized s-transform
    !>
    !> A simple version that assume dt = 1 and the other parameters are default values
    !> derived based on the input sigmal length and dt = 1
    !
    function igst_simple(stran) result(w)

        complex, allocatable, dimension(:, :), intent(in) :: stran
        real, allocatable, dimension(:) :: w

        w = ifft(sum(stran, dim=1), real=.true.)

    end function

    !
    !> Gabor transform
    !>
    !> A simple version that assume dt = 1 and the other parameters are default values
    !> derived based on the input sigmal length and dt = 1
    !
    function gabor_simple(w, sigma) result(stran)

        real, dimension(:), intent(in) :: w
        real, intent(in) :: sigma
        complex, allocatable, dimension(:, :) :: stran

        complex, allocatable, dimension(:) :: ww, e
        integer :: nt, iw, i

        nt = size(w)
        nf0 = nint(nt/2.0)
        stran = zeros(nt, nt)
        ww = zeros(2*nt)
        e = zeros(nt)

        ! Fourier representation of the signal
        ww(1:nt) = fft(w)
        ww(nt + 1:2*nt) = ww(1:nt)

        !$omp parallel do private(iw, i, e)
        do iw = 1, nt
            do i = 1, nt
                e(i) = exp(-2*const_pi**2*((i - nf0))**2*sigma**2)
            end do
            e = cshift(e, -(nt - nf0 + 1))
            stran(:, iw) = ifft(ww(iw:iw + nt - 1)*e)
        end do
        !$omp end parallel do

    end function

    !
    !> Inverse Gabor transform
    !>
    !> A simple version that assume dt = 1 and the other parameters are default values
    !> derived based on the input sigmal length and dt = 1
    !
    function igabor2(stran) result(w)

        complex, allocatable, dimension(:, :), intent(in) :: stran
        real, allocatable, dimension(:) :: w

        w = ifft(sum(stran, dim=1), real=.true.)

    end function 

    !
    !> Compute time-frequency spectrum of a 1D signal using GST or Gabor transform
    !>
    !> In contrast to GST/Gabor, this only returns a time-frequency spectrum
    !> (a real-valued 2D array) rather than complex-valued 2D array
    !
    function tfspec(w, dt, fmin, fmax, df, window, overlap, lambda, p, sigma, method) result(spec)

        real, dimension(:), intent(in) :: w
        real, intent(in), optional :: dt, fmin, fmax, df, window, overlap, lambda, p, sigma
        character(len=*), intent(in), optional :: method
        real, allocatable, dimension(:, :) :: spec

        complex, allocatable, dimension(:, :) :: s
        integer :: nt, window_nt, overlap_nt, nf
        integer(kind=8) :: it
        real :: default_df
        real :: spec_fmin, spec_fmax, spec_df, spec_window, spec_overlap, spec_lambda, spec_p, spec_sigma
        real, allocatable, dimension(:) :: ww
        character(len=24) :: spec_method

        nt = size(w)

        if (present(dt)) then
            spec_dt = dt
        else
            spec_dt = 1.0
        end if

        if (present(fmin)) then
            spec_fmin = fmin
        else
            spec_fmin = 0.0
        end if

        if (present(fmax)) then
            spec_fmax = fmax
        else
            spec_fmax = 0.5/spec_dt
        end if

        default_df = 1.0/spec_dt/nt

        if (present(df)) then
            spec_df = df
        else
            spec_df = default_df
        end if

        if (present(window)) then
            spec_window = window
        else
            spec_window = (nt - 1.0)*spec_dt
        end if

        if (present(overlap)) then
            spec_overlap = overlap
        else
            spec_overlap = 0.5*spec_window
        end if

        if (present(lambda)) then
            spec_lambda = lambda
        else
            spec_lambda = 1.0
        end if

        if (present(p)) then
            spec_p = p
        else
            spec_p = 1.0
        end if

        if (present(sigma)) then
            spec_sigma = sigma
        else
            spec_sigma = 1.0
        end if

        if (present(method)) then
            spec_method = tidy(method)
        else
            spec_method = 'gst'
        end if

        window_nt = next_power_235(nint(spec_window/spec_dt) + 1)
        overlap_nt = nint(spec_overlap/spec_dt) + 1

        it = 1
        do while (it <= nt)

            ! Select signal within a window
            ww = w(it:min(it + window_nt - 1, nt))

            ! Pad the windowed signal if necessary (e.g., at the tail of the original signal)
            if (size(ww) < window_nt) then
                ww = pad(ww, [0, window_nt - size(ww)])
            end if

            ! GST
            select case (spec_method)
                case ('gst')
                    s = gst(ww, spec_lambda, spec_p, spec_dt, spec_fmin, spec_fmax, spec_df)
                case ('gabor')
                    s = gabor(ww, spec_sigma, spec_dt, spec_fmin, spec_fmax, spec_df)
            end select

            ! Tapering
            if (it == 1) then
                s = taper(s, [0, overlap_nt, 0, 0])
                nf = size(s, 2)
                spec = zeros(nt, nf)
            else
                s = taper(s, [overlap_nt, overlap_nt, 0, 0])
            end if

            ! Add to the entire spectrum
            if (it + window_nt - 1 > nt) then
                spec(it:nt, :) = spec(it:nt, :) + abs(s(:nt - it + 1, :))
            else
                spec(it:it + window_nt - 1, :) = spec(it:it + window_nt - 1, :) + abs(s)
            end if

            ! Roll forward
            write (error_unit, *) date_time_compact()//' >> Sample = '//num2str(it)
            it = it + window_nt - overlap_nt

        end do

    end function

end module libflit_spectrum
