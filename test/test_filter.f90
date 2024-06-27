
program test

    use libflit

    call mpistart

    if (nrank == 1) then

        block

            real, allocatable, dimension(:) :: w, v
            integer :: n1

            n1 = 301

            print *, ' 1D '

            w = random(n1, dist='normal')
            w = w - mean(w)
            w = w/std(w)
            w = w/maxval(w)
            call output_array(w, 'w.bin')

            v = gauss_filt(w, 5.0)
            call output_array(v, 'w.bin', append=.true.)

            v = laplace_filt(w)
            call output_array(v, 'w.bin', append=.true.)

            v = median_filt(w, 5)
            call output_array(v, 'w.bin', append=.true.)

            v = mean_filt(w, 2)
            call output_array(v, 'w.bin', append=.true.)

            v = lowess_filt(regspace(0.0, 1.0, n1 - 1.0), w, regspace(0.0, 1.0, n1 - 1.0), window=0.25, order=2)
            call output_array(v, 'w.bin', append=.true.)

            v = fourier_filt(w, dt=0.001, freqs=[0.0, 5.0, 30.0, 40.0], amps=[0.0, 1.0, 1.0, 0.0])
            call output_array(v, 'w.bin', append=.true.)

            v = tv_filt(w, mu=1.0, niter=100)
            call output_array(v, 'w.bin', append=.true.)

            v = balance_filt(w, 5)
            call output_array(v, 'w.bin', append=.true.)

        end block

        block

            real, allocatable, dimension(:, :) :: w, v
            integer :: n1, n2
            type(andf_param) :: p

            n1 = 201
            n2 = 101

            print *, ' 2D '

            w = random(n1, n2, dist='normal')
            w = w - mean(w)
            w = w/std(w)
            w = w/maxval(w)
            call output_array(w, 'w2.bin')

            v = gauss_filt(w, [5.0, 2.0])
            call output_array(v, 'w2.bin', append=.true.)

            v = laplace_filt(w)
            call output_array(v, 'w2.bin', append=.true.)

            v = median_filt(w, [5, 4])
            call output_array(v, 'w2.bin', append=.true.)

            v = mean_filt(w, [5, 7])
            call output_array(v, 'w2.bin', append=.true.)

            ! ! slow
            !            v = reshape(lowess_filt(meshgrid([n1, n2], [1.0, 1.0], [0.0, 0.0], dim=1), &
                !                meshgrid([n1, n2], [1.0, 1.0], [0.0, 0.0], dim=2), &
                !                flatten(w), &
                !                meshgrid([n1, n2], [1.0, 1.0], [0.0, 0.0], dim=1), &
                !                meshgrid([n1, n2], [1.0, 1.0], [0.0, 0.0], dim=2), &
                !                window=[0.25, 0.25], order=2), [n1, n2])
            !            call output_array(v, 'w2.bin', append=.true.)

            v = fourier_filt(w, d1=0.001, freqs1=[0.0, 5.0, 30.0, 40.0], amps1=[0.0, 1.0, 1.0, 0.0], &
                d2=0.001, freqs2=[0.0, 5.0, 30.0, 40.0], amps2=[0.0, 1.0, 1.0, 0.0])
            call output_array(v, 'w2.bin', append=.true.)

            v = tv_filt(w, mu=1.0, niter=100)
            call output_array(v, 'w2.bin', append=.true.)

            v = tgpv_filt(w, mu=1.0, alpha0=1.0, alpha1=0.5, p=0.5,  niter=100)
            call output_array(v, 'w2.bin', append=.true.)

            v = balance_filt(w, [5, 5])
            call output_array(v, 'w2.bin', append=.true.)

            p%smooth1 = 7
            p%smooth2 = 2
            p%sigma = 5
            p%niter = 10
            v = andf_filt(w, p)
            call output_array(v, 'w2.bin', append=.true.)

        end block

        block

            real, allocatable, dimension(:, :, :) :: w, v
            integer :: n1, n2
            type(andf_param) :: p

            n1 = 201
            n2 = 101
            n3 = 51

            print *, ' 3D '

            w = random(n1, n2, n3, dist='normal')
            w = w - mean(w)
            w = w/std(w)
            w = w/maxval(w)
            call output_array(w, 'w3.bin')

            v = gauss_filt(w, [5.0, 2.0, 3.0])
            call output_array(v, 'w3.bin', append=.true.)

            v = laplace_filt(w)
            call output_array(v, 'w3.bin', append=.true.)

            v = median_filt(w, [5, 4, 2])
            call output_array(v, 'w3.bin', append=.true.)

            v = mean_filt(w, [5, 7, 2])
            call output_array(v, 'w3.bin', append=.true.)

            v = fourier_filt(w, d1=0.001, freqs1=[0.0, 5.0, 30.0, 40.0], amps1=[0.0, 1.0, 1.0, 0.0], &
                d2=0.001, freqs2=[0.0, 5.0, 30.0, 40.0], amps2=[0.0, 1.0, 1.0, 0.0], &
                d3=0.001, freqs3=[0.0, 5.0, 30.0, 40.0], amps3=[0.0, 1.0, 1.0, 0.0])
            call output_array(v, 'w3.bin', append=.true.)

            v = tv_filt(w, mu=10.0, niter=100)
            call output_array(v, 'w3.bin', append=.true.)

            v = balance_filt(w, [5, 5, 5])
            call output_array(v, 'w3.bin', append=.true.)

            p%smooth1 = 7
            p%smooth2 = 2
            p%smooth3 = 2
            p%sigma = 5
            p%niter = 10
            v = andf_filt(w, p)
            call output_array(v, 'w3.bin', append=.true.)

        end block

    else

        block

            real, allocatable, dimension(:, :, :) :: w, v
            integer :: n1, n2
            type(andf_param) :: p

            n1 = 201
            n2 = 101
            n3 = 51

            print *, ' 3D '

            w = zeros(n1, n2, n3)

            if (rankid == 0) then
                w = random(n1, n2, n3, dist='normal')
                w = w - mean(w)
                w = w/std(w)
                w = median_filt(w)
                w = w/maxval(w)
            end if
            call mpibarrier
            call bcast_array(w)

            call getpar_int('rank1', rank1, 4)
            call getpar_int('rank2', rank2, 3)
            call getpar_int('rank3', rank3, 2)
            v = tgpv_filt_mpi(w, mu=10.0, alpha0=1.0, alpha1=0.5, p=1.0,  niter=100)
            if (rankid == 0) then
                call output_array(v, 'w3.bin', append=.true.)
            end if

            p%smooth1 = 7
            p%smooth2 = 2
            p%smooth3 = 2
            p%sigma = 5
            p%niter = 10
            v = andf_filt_mpi(w, p)
            if (rankid == 0) then
                call output_array(v, 'w3.bin', append=.true.)
            end if

        end block

    end if

    call mpibarrier
    call mpiend

end program test
