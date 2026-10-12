!
! Gaussian filter gauss_filt (default method 'rsf') on short arrays and for small to huge sigma, in 1D, 2D
! and 3D, in single and double precision. The inputs and outputs go to gaussfilt_out/; check_gaussfilt.py
! compares them with a direct convolution of the same kernels and with two exact conditions (the identity
! for sigma < 1/15, and constants preserved).
!
!   case_<k>_in.bin    input (float32)
!   case_<k>_out.bin   gauss_filt of the input (float32)
!   case_<k>_outd.bin  gauss_filt of the input converted to double precision (float64)
!   cases.txt          k, ndim, n1, n2, n3, sigma1, sigma2, sigma3, signal
!
program test

    use libflit

    implicit none

    integer, dimension(8) :: ns = [1, 2, 3, 4, 5, 7, 10, 16]
    real, dimension(23) :: sigmas = [0.001, 0.01, 0.05, 0.0666, 0.0667, 0.1, 0.25, 0.3832, 0.5, 0.75, 0.999, &
        1.0, 1.5, 2.5, 3.0, 5.0, 10.0, 50.0, 100.0, 1.0e3, 1.0e4, 1.0e5, 1.0e6]
    character(len=8), dimension(5) :: signals = [character(len=8) :: 'random', 'spike1', 'spikemid', 'ramp', 'const']
    character(len=8), dimension(3) :: kinds = [character(len=8) :: 'random', 'const', 'spike1']
    integer, dimension(2, 8) :: shape2
    real, dimension(2, 10) :: sigma2
    integer, dimension(3, 8) :: shape3
    real, dimension(3, 8) :: sigma3
    real, allocatable, dimension(:) :: w1
    real, allocatable, dimension(:, :) :: w2
    real, allocatable, dimension(:, :, :) :: w3
    character(len=1024) :: dir, fk
    integer :: k, a, b, c, funit

    dir = 'gaussfilt_out'
    call make_directory(dir)
    open (newunit=funit, file=tidy(dir)//'/cases.txt', status='replace')
    k = 0

    ! 1D: every length, sigma and signal
    do a = 1, size(ns)
        do b = 1, size(sigmas)
            do c = 1, size(signals)
                k = k + 1
                fk = tidy(dir)//'/case_'//num2str(k)
                w1 = signal_1d(ns(a), signals(c), 1000*a + c)
                call output_array(w1, tidy(fk)//'_in.bin')
                call output_array(gauss_filt(w1, sigmas(b)), tidy(fk)//'_out.bin')
                call output_array(gauss_filt(dble(w1), dble(sigmas(b))), tidy(fk)//'_outd.bin')
                write (funit, '(i6, i3, 3i6, 3es16.8, 1x, a)') k, 1, ns(a), 1, 1, sigmas(b), 0.0, 0.0, trim(signals(c))
            end do
        end do
    end do

    ! 2D: small arrays, with small, large and huge sigma mixed over the axes (a zero sigma leaves that axis
    ! unfiltered); random, constant and corner-spike inputs
    shape2 = reshape([1, 7, 7, 1, 2, 3, 3, 20, 5, 8, 7, 12, 10, 10, 16, 4], [2, 8])
    sigma2 = reshape([0.05, 1.0e4, 0.25, 3.0, 1.0e3, 0.5, 10.0, 10.0, 10.0, 3.0, 20.0, 1.0, 1.0e5, 1.0e5, &
        0.0, 2.5, 0.999, 1.0, 1.0e6, 0.01], [2, 10])
    do a = 1, size(shape2, 2)
        do b = 1, size(sigma2, 2)
            do c = 1, 3
                k = k + 1
                fk = tidy(dir)//'/case_'//num2str(k)
                select case (c)
                    case (1)
                        w2 = random(shape2(1, a), shape2(2, a), range=[0.0, 1.0], seed=20000 + 100*a + b)
                    case (2)
                        w2 = zeros(shape2(1, a), shape2(2, a)) + 1.0
                    case (3)
                        w2 = zeros(shape2(1, a), shape2(2, a))
                        w2(1, 1) = 1.0
                end select
                call output_array(w2, tidy(fk)//'_in.bin')
                call output_array(gauss_filt(w2, sigma2(:, b)), tidy(fk)//'_out.bin')
                call output_array(gauss_filt(dble(w2), dble(sigma2(:, b))), tidy(fk)//'_outd.bin')
                write (funit, '(i6, i3, 3i6, 3es16.8, 1x, a)') k, 2, shape2(:, a), 1, sigma2(:, b), 0.0, trim(kinds(c))
            end do
        end do
    end do

    ! 3D: the same for small volumes
    shape3 = reshape([1, 1, 5, 1, 6, 1, 2, 3, 4, 2, 3, 30, 4, 4, 4, 4, 5, 6, 6, 8, 10, 10, 10, 10], [3, 8])
    sigma3 = reshape([0.05, 0.25, 0.5, 1.0e3, 10.0, 0.01, 0.999, 1.0, 1.0e4, 1.0e5, 1.0e5, 1.0e5, &
        3.0, 0.0, 50.0, 1.0e6, 0.3832, 2.5, 10.0, 10.0, 10.0, 5.0, 20.0, 2.0], [3, 8])
    do a = 1, size(shape3, 2)
        do b = 1, size(sigma3, 2)
            do c = 1, 3
                k = k + 1
                fk = tidy(dir)//'/case_'//num2str(k)
                select case (c)
                    case (1)
                        w3 = random(shape3(1, a), shape3(2, a), shape3(3, a), range=[0.0, 1.0], seed=30000 + 100*a + b)
                    case (2)
                        w3 = zeros(shape3(1, a), shape3(2, a), shape3(3, a)) + 1.0
                    case (3)
                        w3 = zeros(shape3(1, a), shape3(2, a), shape3(3, a))
                        w3(1, 1, 1) = 1.0
                end select
                call output_array(w3, tidy(fk)//'_in.bin')
                call output_array(gauss_filt(w3, sigma3(:, b)), tidy(fk)//'_out.bin')
                call output_array(gauss_filt(dble(w3), dble(sigma3(:, b))), tidy(fk)//'_outd.bin')
                write (funit, '(i6, i3, 3i6, 3es16.8, 1x, a)') k, 3, shape3(:, a), sigma3(:, b), trim(kinds(c))
            end do
        end do
    end do

    close (funit)
    print *, k, ' cases written to '//tidy(dir)

contains

    function signal_1d(n, kind, seed) result(w)

        integer, intent(in) :: n, seed
        character(len=*), intent(in) :: kind
        real, allocatable, dimension(:) :: w

        integer :: i

        w = zeros(n)
        select case (kind)
            case ('random')
                w = random(n, range=[0.0, 1.0], seed=seed)
            case ('spike1')
                w(1) = 1.0
            case ('spikemid')
                w((n + 1)/2) = 1.0
            case ('ramp')
                do i = 1, n
                    w(i) = i - 1.0
                end do
            case ('const')
                w = 1.0
        end select

    end function signal_1d

end program test
