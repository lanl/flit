
program test

    use libflit

    block

        real, allocatable, dimension(:) :: w, v
        integer :: n1

        n1 = 100

        print *, ' 1D '

        w = random(n1, dist='uniform')

        print *

        call plot_histogram(w)

        v = random_permute(w)

        print *
        print *, random_order(10)

        print *
        do i = 1, 10
            print *, rand(seed=1111), rand(), random(1, seed=2222), random(1)
        end do

        print *

        do i = 1, 30
            print *, random_string(i)
        end do

        print *
        call plot_histogram(random(1000, dist='exponential', lambda=1.0))
        
        print *
        call plot_histogram(random(1000, dist='exponential', lambda=0.1))

    end block

    block

        real, allocatable, dimension(:, :) :: w, v

        print *, irandom(3, 4)

        print *
        print *, random(3, 4, dist='normal')

        print *
        print *, random(1, 2, 3, dist='exponential')

        print *
        call plot_histogram(random(300, 400, dist='exponential'))

        print *
        call plot_histogram(random(300, 400, dist='uniform'))

        print *
        call plot_histogram(random(20, 30, 40, dist='normal'))

        w = random(50, 100)
        v = random(50, 100, dist='normal')
        print *, mean(w)
        print *, median(w)
        print *, std(w)
        print *, minval(w), maxval(w)
        print *, mean(v)
        print *, median(v)
        print *, std(v)
        print *, minval(v), maxval(v)

        ! The output dimension is 99x199
        v = xcorr(w, v)
        call output_array(v, './v.bin')

        v = acorr(w)
        call output_array(v, './v.bin', append=.true.)

        v = random(4, 4, dist='normal')
        v = conv(w, v, 'same')
        call output_array(v, 'vv.bin')

    end block

end program test
