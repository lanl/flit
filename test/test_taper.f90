
program test

    use libflit

    real, allocatable, dimension(:) :: w
    character(len=32), allocatable, dimension(:) :: m

    m = ['step', 'kaiser', 'linear', 'parzen', 'welch', 'sine', 'power-sine', &
        'hann', 'hamming', 'blackman', 'nuttall', 'blackman-nuttall', 'blackman-harris', 'gauss']

    open (3, file='w.bin', form='unformatted', access='stream', status='replace')
    close (3)
    do i = 1, size(m)

        if (m(i) == 'kaiser') alpha = 4.14
        if (m(i) == 'power-sine') alpha = 2.0
        if (m(i) == 'gauss') alpha = 0.4

        w = ones(100)
        w = taper(w, len=[20, 30], method=[m(i), m(i)], alpha=[alpha, alpha])
        call output_array(w, './w.bin', append=.true.)

        print *, i, m(i), minval(w), maxval(w)

    end do

end program test
