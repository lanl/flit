
program test

    use libflit

    block

        real, allocatable, dimension(:) :: a
        double precision, allocatable, dimension(:) :: b
        complex, allocatable, dimension(:) :: c
        double complex, allocatable, dimension(:) :: d

        a = random(10)
        b = random(10)
        c = cmplx(random(10), random(10))
        d = dcmplx(random(10), random(10))

        print *, fft(a)
        print *, a - ifft(fft(a), real=.true.)

        print *

        print *, fft(b)
        print *, fft(b, n=next_power_235(size(b)))

        print *

        print *, fft(c)
        print *, c - ifft(fft(c))

        print *

        print *, fft(d)
        print *, ifft(d, real=.true.)

        print *

        print *, hilbert(a)

        print *

        print *, envelope(a)

    end block

    block

        real, allocatable, dimension(:, :) :: a
        double precision, allocatable, dimension(:, :) :: b
        complex, allocatable, dimension(:, :) :: c
        double complex, allocatable, dimension(:, :) :: d

        a = random(5, 5)
        b = random(5, 5)
        c = cmplx(random(5, 5), random(5, 5))
        d = dcmplx(random(5, 5), random(5, 5))

        print *, fft(a)
        print *, a - ifft(fft(a), real=.true.)

        print *

        print *, fft(b)
        print *, fft(b, n=[next_power_235(size(b, 1)), next_power_235(size(b, 2))])

        print *

        print *, fft(c)
        print *, c - ifft(fft(c))

        print *

        print *, fft(d)
        print *, ifft(d, real=.true.)

        print *

        print *, dct(a)

    end block

end program test
