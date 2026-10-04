
!
! Test TV and TGpV filtering against ground-truth images, and record the
! convergence curves (error versus number of iterations).
!
!   1. Constant image: a constant is an exact fixed point of TV and TGpV
!      denoising, so the filtered image must equal the input to rounding,
!      at every iteration count and for any p.
!
!   2. Noisy image with a known ground truth (linear trends plus steps):
!      the filtered image must end up much closer to the ground truth
!      than the noisy input is.
!
! The curves are printed and written to test_tvfilt_convergence.txt;
! plot_tvfilt_convergence.py plots them.
!
program test

    use libflit

    implicit none

    integer, parameter :: nit = 6
    integer, dimension(nit), parameter :: its = [1, 3, 10, 30, 100, 300]
    double precision, parameter :: alpha0 = 1.0d0, alpha1 = 2.0d0
    double precision, parameter :: c = 2000.0d0, mu_const = 100.0d0/c, mu_noisy = 5.0d0
    double precision, parameter :: tol_const = 1.0d-12, tol_noisy = 0.6d0
    character(len=8), dimension(3), parameter :: fname = [character(len=8) :: 'tv', 'tgpv', 'tgpv']
    double precision, dimension(3), parameter :: fp = [0.0d0, 1.0d0, 0.5d0]

    double precision, allocatable, dimension(:) :: c1, t1, g1
    double precision, allocatable, dimension(:, :) :: c2, t2, g2
    double precision, allocatable, dimension(:, :, :) :: c3, t3, g3
    double precision, dimension(nit, 3, 3, 2) :: err
    real, allocatable, dimension(:) :: r1
    real, allocatable, dimension(:, :) :: r2
    real, allocatable, dimension(:, :, :) :: r3
    double precision :: x1, x2, x3
    integer :: n1, n2, n3, i, j, k, ifilt, nd, fam, npass, nfail, fid
    logical :: ok
    character(len=32) :: case_name

    ! Constant images
    c1 = zeros(101) + c
    c2 = zeros(80, 60) + c
    c3 = zeros(30, 24, 20) + c

    ! 1D ground truth and noisy signal
    n1 = 401
    r1 = random(n1, dist='normal', seed=5)
    allocate (t1(n1))
    do i = 1, n1
        x1 = (i - 1.0d0)/(n1 - 1.0d0)
        t1(i) = 10.0d0 + 2.0d0*x1
        if (x1 > 0.35d0) then
            t1(i) = t1(i) + 1.5d0
        end if
        if (x1 > 0.7d0) then
            t1(i) = t1(i) - 1.0d0
        end if
    end do
    g1 = t1 + 0.15d0*r1

    ! 2D ground truth and noisy image
    n1 = 48
    n2 = 40
    r2 = random(n1, n2, dist='normal', seed=7)
    allocate (t2(n1, n2))
    do j = 1, n2
        do i = 1, n1
            x1 = (i - 1.0d0)/(n1 - 1.0d0)
            x2 = (j - 1.0d0)/(n2 - 1.0d0)
            t2(i, j) = 10.0d0 + 2.0d0*x1 + 1.0d0*x2
            if (x1 > 0.35d0 + 0.3d0*x2) then
                t2(i, j) = t2(i, j) + 1.5d0
            end if
        end do
    end do
    g2 = t2 + 0.15d0*r2

    ! 3D ground truth and noisy image
    n1 = 30
    n2 = 24
    n3 = 20
    r3 = random(n1, n2, n3, dist='normal', seed=9)
    allocate (t3(n1, n2, n3))
    do k = 1, n3
        do j = 1, n2
            do i = 1, n1
                x1 = (i - 1.0d0)/(n1 - 1.0d0)
                x2 = (j - 1.0d0)/(n2 - 1.0d0)
                x3 = (k - 1.0d0)/(n3 - 1.0d0)
                t3(i, j, k) = 10.0d0 + 2.0d0*x1 + 1.0d0*x2 - 1.0d0*x3
                if (x1 > 0.3d0 + 0.3d0*x2 + 0.2d0*x3) then
                    t3(i, j, k) = t3(i, j, k) + 1.5d0
                end if
            end do
        end do
    end do
    g3 = t3 + 0.15d0*r3

    ! Convergence curves; family 1 = constant image, family 2 = noisy image
    do fam = 1, 2
        do nd = 1, 3
            do ifilt = 1, 3
                do k = 1, nit
                    err(k, ifilt, nd, fam) = error_after(fam, ifilt, nd, its(k))
                end do
            end do
        end do
    end do

    ! Report
    npass = 0
    nfail = 0
    open (newunit=fid, file='test_tvfilt_convergence.txt', status='replace')
    write (fid, '(a)') '# FLIT test_tvfilt convergence curves'
    write (fid, '(a)') '# columns: image dimension filter p, then the error after each iteration count'
    write (fid, '(a, *(1x, i0))') '# iterations', its

    do fam = 1, 2
        if (fam == 1) then
            print '(/, a)', ' Constant image c = 2000, mu = 100/c, alpha0 = 1, alpha1 = 2: max|u - c|/c'
        else
            print '(/, a)', ' Noisy image, noise std = 0.15, mu = 5, alpha0 = 1, alpha1 = 2: rms(u - truth)/rms(g - truth)'
        end if
        print '(1x, a24, a6, *(i10))', 'Iterations', '', its
        do nd = 1, 3
            do ifilt = 1, 3
                if (fam == 1) then
                    ok = maxval(err(:, ifilt, nd, fam)) <= tol_const
                else
                    ok = err(nit, ifilt, nd, fam) <= tol_noisy
                end if
                if (ok) then
                    npass = npass + 1
                else
                    nfail = nfail + 1
                end if
                write (case_name, '(a, a, i0, a)') trim(fname(ifilt))//'_filt', ' ', nd, 'D'
                if (ifilt > 1) then
                    case_name = trim(case_name)//', p = '//merge('1  ', '0.5', ifilt == 2)
                end if
                print '(1x, a24, a6, *(es10.2))', case_name, merge('  PASS', '  FAIL', ok), err(:, ifilt, nd, fam)
                write (fid, '(a, 1x, i0, 1x, a, 1x, f3.1, *(1x, es14.6))') &
                    merge('constant', 'noisy   ', fam == 1), nd, trim(fname(ifilt))//'_filt', fp(ifilt), err(:, ifilt, nd, fam)
            end do
        end do
    end do
    close (fid)

    print '(/, a, es8.1, a, f4.2, a)', ' Pass: constant image max|u - c|/c <= ', tol_const, &
        ' at every count; noisy image error ratio <= ', tol_noisy, ' at the last count'
    print '(1x, i0, a, i0, a)', npass, ' passed, ', nfail, ' failed'
    if (nfail > 0) then
        error stop ' <test_tvfilt> Error: some checks failed'
    end if

contains

    !
    !> Error against the ground truth after niter iterations
    !
    function error_after(fam, ifilt, nd, niter) result(e)

        integer, intent(in) :: fam, ifilt, nd, niter
        double precision :: e

        select case (nd)
            case (1)
                if (fam == 1) then
                    e = maxval(abs(filt1(c1, mu_const, ifilt, niter) - c))/c
                else
                    e = rms1(filt1(g1, mu_noisy, ifilt, niter) - t1)/rms1(g1 - t1)
                end if
            case (2)
                if (fam == 1) then
                    e = maxval(abs(filt2(c2, mu_const, ifilt, niter) - c))/c
                else
                    e = rms1(pack(filt2(g2, mu_noisy, ifilt, niter) - t2, .true.))/rms1(pack(g2 - t2, .true.))
                end if
            case (3)
                if (fam == 1) then
                    e = maxval(abs(filt3(c3, mu_const, ifilt, niter) - c))/c
                else
                    e = rms1(pack(filt3(g3, mu_noisy, ifilt, niter) - t3, .true.))/rms1(pack(g3 - t3, .true.))
                end if
        end select

    end function error_after

    function filt1(f, mu, ifilt, niter) result(u)

        double precision, dimension(:), intent(in) :: f
        double precision, intent(in) :: mu
        integer, intent(in) :: ifilt, niter
        double precision, allocatable, dimension(:) :: u

        if (ifilt == 1) then
            u = tv_filt(f, mu, niter, verbose=.false.)
        else
            u = tgpv_filt(f, mu, alpha0, alpha1, niter, fp(ifilt), verbose=.false.)
        end if

    end function filt1

    function filt2(f, mu, ifilt, niter) result(u)

        double precision, dimension(:, :), intent(in) :: f
        double precision, intent(in) :: mu
        integer, intent(in) :: ifilt, niter
        double precision, allocatable, dimension(:, :) :: u

        if (ifilt == 1) then
            u = tv_filt(f, mu, niter, verbose=.false.)
        else
            u = tgpv_filt(f, mu, alpha0, alpha1, niter, fp(ifilt), verbose=.false.)
        end if

    end function filt2

    function filt3(f, mu, ifilt, niter) result(u)

        double precision, dimension(:, :, :), intent(in) :: f
        double precision, intent(in) :: mu
        integer, intent(in) :: ifilt, niter
        double precision, allocatable, dimension(:, :, :) :: u

        if (ifilt == 1) then
            u = tv_filt(f, mu, niter, verbose=.false.)
        else
            ! The 3D tgpv_filt has no verbose argument and reports progress on stderr
            u = tgpv_filt(f, mu, alpha0, alpha1, niter, fp(ifilt))
        end if

    end function filt3

    function rms1(w) result(r)

        double precision, dimension(:), intent(in) :: w
        double precision :: r

        r = sqrt(sum(w**2)/size(w))

    end function rms1

end program test
