
program test

    use libflit

    block

        real, allocatable, dimension(:) :: w, ww
        character(len=32), allocatable, dimension(:) :: m
        integer :: o, d
        integer :: n, nn, i
        real :: start, finish

        n = 303
        o = 2
        d = 3

        w = gauss_filt(random(n), 4.0)
        w = rescale(w, [0.0, 1.0])
        call output_array(w, 'w.bin')
        w = w(o + 1:n:d + 1)
        nn = size(w)

        m = ['nearest', 'linear', 'cubic', 'pchip', 'sinc', 'quintic', 'mba', 'biharmonic', 'cubic_spline', 'hermite_spline', 'monotonic_spline']

        do i = 1, size(m)

            call cpu_time(start)

            do j = 1, 1
                ww = interp(w, nn, d + 1.0, o*1.0, n, 1.0, 0.0, method=m(i))
            end do

            call cpu_time(finish)

            print *, i, m(i), finish - start

            call output_array(ww, './w_'//tidy(m(i))//'.bin')
            ! This needs github.com/lanl/pymplot
            call execute_command_line('x_showgraph -in=w.bin,w_'//tidy(m(i))//'.bin -n1='//num2str(n)//','//num2str(n) &
                //' -linecolor=b,r -x2beg=-0.1 -x2end=1.1 -size2=3 -size1=10 -linestyle=solid -linewidth=2,2 -out=interp_' &
                //tidy(m(i))//'.png &')

        end do

    end block

    block

        real, allocatable, dimension(:, :, :) :: w, ww
        character(len=24) :: m = 'linear'

        w = gauss_filt(random(200, 300, 100), [4.0, 4.0, 4.0])
        w = w - mean(w)
        call output_array(w, './w3.bin')
        ww = w

        w = interp_to(w, nint(shape(ww)/2.5), [m, m, m])
        w = interp_to(w, shape(ww), [m, m, m])
        call output_array(w, './w3.bin', append=.true.)
        call output_array(w - ww, './w3.bin', append=.true.)

    end block

    block

        real, allocatable, dimension(:, :) :: w, ww
        character(len=24) :: m = 'sinc'

        w = gauss_filt(random(300, 400), [4.0, 4.0])
        w = w - mean(w)
        call output_array(w, './w2.bin')
        ww = w

        w = interp_to(w, [71, 81], [m, m])
        w = interp_to(w, [300, 400], [m, m])
        call output_array(w, './w2.bin', append=.true.)
        call output_array(w - ww, './w2.bin', append=.true.)

    end block

    stop


end program test
