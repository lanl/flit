
module iir_design_butterworth_mod

    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    integer, parameter, public :: dp = real64
    real(dp), parameter :: pi = acos(-1.0_dp)

    private
    public :: design_butterworth_sos
    public :: print_sos

contains

    !=======================================================================
    ! Main public routine
    !=======================================================================

    subroutine design_butterworth_sos(kind, order, fs, freqmin, freqmax, sos)
        character(len=*), intent(in) :: kind
        integer,          intent(in) :: order
        real(dp),         intent(in) :: fs
        real(dp),         intent(in) :: freqmin
        real(dp),         intent(in) :: freqmax
        real(dp), allocatable, intent(out) :: sos(:, :)

        complex(dp), allocatable :: p0(:)
        complex(dp), allocatable :: za(:), pa(:)
        complex(dp), allocatable :: zd(:), pd(:)

        character(len=:), allocatable :: k
        real(dp) :: w1, w2, wc, f1, f2, center, width
        integer :: n

        if (order < 1) error stop "design_butterworth_sos: order must be >= 1."
        if (fs <= 0.0_dp) error stop "design_butterworth_sos: fs must be positive."

        n = order
        k = trim(adjustl(kind))

        call butterworth_prototype_poles(n, p0)

        select case (k)

            case ("lowpass", "lp")
                if (freqmax <= 0.0_dp .or. freqmax >= 0.5_dp * fs) then
                    error stop "lowpass: freqmax must be between 0 and Nyquist."
                end if

                wc = prewarp(freqmax, fs)
                call analog_lp_to_lp(p0, wc, za, pa)

            case ("highpass", "hp")
                if (freqmin <= 0.0_dp .or. freqmin >= 0.5_dp * fs) then
                    error stop "highpass: freqmin must be between 0 and Nyquist."
                end if

                wc = prewarp(freqmin, fs)
                call analog_lp_to_hp(p0, wc, za, pa)

            case ("bandpass", "bp")
                if (freqmin <= 0.0_dp .or. freqmax <= freqmin .or. freqmax >= 0.5_dp * fs) then
                    error stop "bandpass: require 0 < freqmin < freqmax < Nyquist."
                end if

                w1 = prewarp(freqmin, fs)
                w2 = prewarp(freqmax, fs)
                call analog_lp_to_bp(p0, w1, w2, za, pa)

            case ("bandstop", "bs", "bandreject", "br")
                if (freqmin <= 0.0_dp .or. freqmax <= freqmin .or. freqmax >= 0.5_dp * fs) then
                    error stop "bandstop: require 0 < freqmin < freqmax < Nyquist."
                end if

                w1 = prewarp(freqmin, fs)
                w2 = prewarp(freqmax, fs)
                call analog_lp_to_bs(p0, w1, w2, za, pa)

            case ("freq_remove", "spike_remove", "notch")
                center = freqmin
                width  = freqmax

                if (center <= 0.0_dp .or. width <= 0.0_dp) then
                    error stop "freq_remove: center frequency and width must be positive."
                end if

                f1 = center - 0.5_dp * width
                f2 = center + 0.5_dp * width

                if (f1 <= 0.0_dp .or. f2 >= 0.5_dp * fs) then
                    error stop "freq_remove: require 0 < center-width/2 and center+width/2 < Nyquist."
                end if

                w1 = prewarp(f1, fs)
                w2 = prewarp(f2, fs)
                call analog_lp_to_bs(p0, w1, w2, za, pa)

            case default
                error stop "Unknown filter kind."
        end select

        call bilinear_zpk(za, pa, fs, zd, pd)
        call zpk_to_sos(zd, pd, sos)
        call normalize_sos_gain(k, fs, freqmin, freqmax, sos)

    end subroutine design_butterworth_sos


    !=======================================================================
    ! Butterworth normalized analog low-pass prototype
    !=======================================================================

    subroutine butterworth_prototype_poles(n, p)
        integer, intent(in) :: n
        complex(dp), allocatable, intent(out) :: p(:)

        integer :: k
        real(dp) :: theta

        allocate(p(n))

        do k = 1, n
            theta = pi * real(2*k + n - 1, dp) / real(2*n, dp)

            p(k) = cmplx(cos(theta), sin(theta), kind=dp)

            ! Numerical safety: keep stable left-half-plane pole.
            if (real(p(k), dp) > 0.0_dp) p(k) = -p(k)
        end do

    end subroutine butterworth_prototype_poles


    !=======================================================================
    ! Frequency prewarping for bilinear transform
    !=======================================================================

    pure function prewarp(f, fs) result(w)
        real(dp), intent(in) :: f
        real(dp), intent(in) :: fs
        real(dp) :: w

        w = 2.0_dp * fs * tan(pi * f / fs)

    end function prewarp


    !=======================================================================
    ! Analog low-pass to low-pass
    !=======================================================================

    subroutine analog_lp_to_lp(p0, wc, z, p)
        complex(dp), intent(in) :: p0(:)
        real(dp),    intent(in) :: wc
        complex(dp), allocatable, intent(out) :: z(:), p(:)

        integer :: n

        n = size(p0)

        allocate(z(0))
        allocate(p(n))

        p = wc * p0

    end subroutine analog_lp_to_lp


    !=======================================================================
    ! Analog low-pass to high-pass
    !
    ! Transformation:
    !     s_lp = wc / s
    !
    ! Poles:
    !     p_hp = wc / p_lp
    !
    ! Zeros:
    !     z = 0, repeated n times
    !=======================================================================

    subroutine analog_lp_to_hp(p0, wc, z, p)
        complex(dp), intent(in) :: p0(:)
        real(dp),    intent(in) :: wc
        complex(dp), allocatable, intent(out) :: z(:), p(:)

        integer :: n, i

        n = size(p0)

        allocate(z(n))
        allocate(p(n))

        do i = 1, n
            p(i) = wc / p0(i)
            z(i) = cmplx(0.0_dp, 0.0_dp, kind=dp)
        end do

    end subroutine analog_lp_to_hp


    !=======================================================================
    ! Analog low-pass to band-pass
    !
    ! Transformation:
    !     s_lp = (s^2 + w0^2) / (B s)
    !
    ! For each prototype pole p0:
    !     s^2 - B*p0*s + w0^2 = 0
    !
    ! Zeros:
    !     z = 0, repeated n times
    !=======================================================================

    subroutine analog_lp_to_bp(p0, w1, w2, z, p)
        complex(dp), intent(in) :: p0(:)
        real(dp),    intent(in) :: w1, w2
        complex(dp), allocatable, intent(out) :: z(:), p(:)

        integer :: n, i
        real(dp) :: b, w0
        complex(dp) :: disc

        n  = size(p0)
        b  = w2 - w1
        w0 = sqrt(w1 * w2)

        allocate(z(n))
        allocate(p(2*n))

        do i = 1, n
            z(i) = cmplx(0.0_dp, 0.0_dp, kind=dp)

            disc = sqrt((b * p0(i))**2 - 4.0_dp * w0**2)

            p(2*i-1) = 0.5_dp * (b * p0(i) + disc)
            p(2*i  ) = 0.5_dp * (b * p0(i) - disc)
        end do

    end subroutine analog_lp_to_bp


    !=======================================================================
    ! Analog low-pass to band-stop
    !
    ! Transformation:
    !     s_lp = B s / (s^2 + w0^2)
    !
    ! For each prototype pole p0:
    !     p0*s^2 - B*s + p0*w0^2 = 0
    !
    ! Zeros:
    !     z = +i*w0 and -i*w0, each repeated n times
    !=======================================================================

    subroutine analog_lp_to_bs(p0, w1, w2, z, p)
        complex(dp), intent(in) :: p0(:)
        real(dp),    intent(in) :: w1, w2
        complex(dp), allocatable, intent(out) :: z(:), p(:)

        integer :: n, i
        real(dp) :: b, w0
        complex(dp) :: disc
        complex(dp) :: iw0

        n  = size(p0)
        b  = w2 - w1
        w0 = sqrt(w1 * w2)

        iw0 = cmplx(0.0_dp, w0, kind=dp)

        allocate(z(2*n))
        allocate(p(2*n))

        do i = 1, n
            z(2*i-1) =  iw0
            z(2*i  ) = -iw0

            disc = sqrt(b**2 - 4.0_dp * (p0(i)**2) * w0**2)

            p(2*i-1) = (b + disc) / (2.0_dp * p0(i))
            p(2*i  ) = (b - disc) / (2.0_dp * p0(i))
        end do

    end subroutine analog_lp_to_bs


    !=======================================================================
    ! Bilinear transform from analog z/p to digital z/p
    !
    !     z_d = (2 fs + z_a) / (2 fs - z_a)
    !     p_d = (2 fs + p_a) / (2 fs - p_a)
    !
    ! If there are fewer zeros than poles, extra zeros are added at z = -1.
    !=======================================================================

    subroutine bilinear_zpk(za, pa, fs, zd, pd)
        complex(dp), intent(in) :: za(:)
        complex(dp), intent(in) :: pa(:)
        real(dp),    intent(in) :: fs
        complex(dp), allocatable, intent(out) :: zd(:), pd(:)

        integer :: nz, np, i
        complex(dp) :: twofs

        nz = size(za)
        np = size(pa)

        if (np < nz) error stop "bilinear_zpk: more zeros than poles."

        twofs = cmplx(2.0_dp * fs, 0.0_dp, kind=dp)

        allocate(zd(np))
        allocate(pd(np))

        do i = 1, nz
            zd(i) = (twofs + za(i)) / (twofs - za(i))
        end do

        do i = nz + 1, np
            zd(i) = cmplx(-1.0_dp, 0.0_dp, kind=dp)
        end do

        do i = 1, np
            pd(i) = (twofs + pa(i)) / (twofs - pa(i))
        end do

    end subroutine bilinear_zpk


    !=======================================================================
    ! Convert digital zeros/poles to SOS.
    !
    ! sos(:, section) = [b0, b1, b2, a0, a1, a2]
    !
    ! This routine pairs complex-conjugate roots into real second-order
    ! sections. Real leftover roots become first-order sections represented as:
    !
    !     b0 + b1 z^-1
    !     a0 + a1 z^-1
    !
    ! with b2 = a2 = 0.
    !=======================================================================

    subroutine zpk_to_sos(z, p, sos)
        complex(dp), intent(in) :: z(:)
        complex(dp), intent(in) :: p(:)
        real(dp), allocatable, intent(out) :: sos(:, :)

        logical, allocatable :: used_z(:), used_p(:)
        integer :: nroot, nsec
        integer :: sec
        integer :: ip1, ip2, iz1, iz2
        integer :: p_order, z_order
        complex(dp) :: r1, r2
        real(dp) :: b0, b1, b2
        real(dp) :: a0, a1, a2

        nroot = size(p)

        if (size(z) /= nroot) error stop "zpk_to_sos: zero/pole count mismatch."

        nsec = (nroot + 1) / 2

        allocate(sos(6, nsec))
        sos = 0.0_dp

        allocate(used_z(nroot))
        allocate(used_p(nroot))
        used_z = .false.
        used_p = .false.

        do sec = 1, nsec

            call pick_root_pair(p, used_p, ip1, ip2, p_order)
            call pick_root_pair(z, used_z, iz1, iz2, z_order)

            ! Numerator
            if (z_order == 2) then
                r1 = z(iz1)
                r2 = z(iz2)

                b0 = 1.0_dp
                b1 = -real(r1 + r2, dp)
                b2 =  real(r1 * r2, dp)

            else if (z_order == 1) then
                r1 = z(iz1)

                b0 = 1.0_dp
                b1 = -real(r1, dp)
                b2 = 0.0_dp

            else
                error stop "zpk_to_sos: invalid zero order."
            end if

            ! Denominator
            if (p_order == 2) then
                r1 = p(ip1)
                r2 = p(ip2)

                a0 = 1.0_dp
                a1 = -real(r1 + r2, dp)
                a2 =  real(r1 * r2, dp)

            else if (p_order == 1) then
                r1 = p(ip1)

                a0 = 1.0_dp
                a1 = -real(r1, dp)
                a2 = 0.0_dp

            else
                error stop "zpk_to_sos: invalid pole order."
            end if

            sos(:, sec) = [b0, b1, b2, a0, a1, a2]

        end do

    end subroutine zpk_to_sos


    !=======================================================================
    ! Pick a root pair.
    !
    ! Preference:
    !   1. complex root with its conjugate
    !   2. two real roots
    !   3. one leftover real root
    !=======================================================================

    subroutine pick_root_pair(r, used, i1, i2, root_order)
        complex(dp), intent(in) :: r(:)
        logical, intent(inout) :: used(:)
        integer, intent(out) :: i1, i2
        integer, intent(out) :: root_order

        integer :: n, i, j
        real(dp), parameter :: tol = 1.0e-10_dp

        n = size(r)

        i1 = 0
        i2 = 0
        root_order = 0

        ! First pick an unused complex root and its conjugate.
        do i = 1, n
            if (.not. used(i)) then
                if (abs(aimag(r(i))) > tol) then
                    i1 = i

                    do j = 1, n
                        if (j /= i .and. .not. used(j)) then
                            if (abs(r(j) - conjg(r(i))) < 1.0e-8_dp) then
                                i2 = j
                                used(i1) = .true.
                                used(i2) = .true.
                                root_order = 2
                                return
                            end if
                        end if
                    end do

                    error stop "pick_root_pair: complex root has no conjugate."
                end if
            end if
        end do

        ! Otherwise pick two real roots if possible.
        do i = 1, n
            if (.not. used(i)) then
                if (abs(aimag(r(i))) <= tol) then
                    i1 = i

                    do j = i + 1, n
                        if (.not. used(j)) then
                            if (abs(aimag(r(j))) <= tol) then
                                i2 = j
                                used(i1) = .true.
                                used(i2) = .true.
                                root_order = 2
                                return
                            end if
                        end if
                    end do

                    ! One real root left.
                    used(i1) = .true.
                    i2 = 0
                    root_order = 1
                    return
                end if
            end if
        end do

        error stop "pick_root_pair: no unused roots."

    end subroutine pick_root_pair


    !=======================================================================
    ! Normalize digital gain.
    !
    ! We scale the numerator of the first SOS section so that the passband
    ! gain is approximately one at a representative frequency.
    !
    ! lowpass:     f_ref = 0
    ! highpass:    f_ref = Nyquist
    ! bandpass:    f_ref = sqrt(f1*f2)
    ! bandstop:    f_ref = 0
    ! freq_remove: f_ref = 0
    !=======================================================================

    subroutine normalize_sos_gain(kind, fs, freqmin, freqmax, sos)
        character(len=*), intent(in) :: kind
        real(dp), intent(in) :: fs
        real(dp), intent(in) :: freqmin, freqmax
        real(dp), intent(inout) :: sos(:, :)

        character(len=:), allocatable :: k
        real(dp) :: fref, omega, mag, scale

        k = trim(adjustl(kind))

        select case (k)

            case ("lowpass", "lp")
                fref = 0.0_dp

            case ("highpass", "hp")
                fref = 0.5_dp * fs

            case ("bandpass", "bp")
                fref = sqrt(freqmin * freqmax)

            case ("bandstop", "bs", "bandreject", "br")
                fref = 0.0_dp

            case ("freq_remove", "spike_remove", "notch")
                fref = 0.0_dp

            case default
                fref = 0.0_dp

        end select

        omega = 2.0_dp * pi * fref / fs

        mag = sos_freq_response_mag(sos, omega)

        if (mag <= 0.0_dp) then
            error stop "normalize_sos_gain: zero gain at reference frequency."
        end if

        scale = 1.0_dp / mag

        sos(1, 1) = scale * sos(1, 1)
        sos(2, 1) = scale * sos(2, 1)
        sos(3, 1) = scale * sos(3, 1)

    end subroutine normalize_sos_gain


    !=======================================================================
    ! Evaluate magnitude response of SOS at digital angular frequency omega.
    !
    ! H(z), z = exp(i omega)
    !
    ! SOS convention:
    !     b0 + b1 z^-1 + b2 z^-2
    !     -------------------------
    !     a0 + a1 z^-1 + a2 z^-2
    !=======================================================================

    function sos_freq_response_mag(sos, omega) result(mag)
        real(dp), intent(in) :: sos(:, :)
        real(dp), intent(in) :: omega
        real(dp) :: mag

        integer :: s, nsec
        complex(dp) :: z1, z2
        complex(dp) :: h, num, den

        nsec = size(sos, 2)

        z1 = exp(cmplx(0.0_dp, -omega, kind=dp))
        z2 = exp(cmplx(0.0_dp, -2.0_dp * omega, kind=dp))

        h = cmplx(1.0_dp, 0.0_dp, kind=dp)

        do s = 1, nsec
            num = cmplx(sos(1, s), 0.0_dp, kind=dp) &
                + cmplx(sos(2, s), 0.0_dp, kind=dp) * z1 &
                + cmplx(sos(3, s), 0.0_dp, kind=dp) * z2

            den = cmplx(sos(4, s), 0.0_dp, kind=dp) &
                + cmplx(sos(5, s), 0.0_dp, kind=dp) * z1 &
                + cmplx(sos(6, s), 0.0_dp, kind=dp) * z2

            h = h * num / den
        end do

        mag = abs(h)

    end function sos_freq_response_mag


    !=======================================================================
    ! Utility: print SOS coefficients
    !=======================================================================

    subroutine print_sos(sos)
        real(dp), intent(in) :: sos(:, :)
        integer :: i

        do i = 1, size(sos, 2)
            write(*,'(i4,2x,6(es24.16,1x))') i, sos(:, i)
        end do

    end subroutine print_sos

end module iir_design_butterworth_mod
