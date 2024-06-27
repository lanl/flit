
program test

    use libflit

    block

        real, allocatable, dimension(:) :: w, v
        integer :: n1

        n1 = 3

        print *, ' 1D '

        w = zeros(n1)
        print *, w

        w = ones(n1)
        print *, w

        w = const(n1, 2.0)
        print *, w

        w = regspace(1.0, 1.0, 5.0)
        print *, w

        w = linspace(2.0, 5.0, 2*n1)
        print *, w

        v = zeros_like(w)
        v = ones_like(w)

    end block

    block

        real, allocatable, dimension(:, :) :: w, v
        integer :: n1, n2
        integer :: i

        n1 = 3
        n2 = 4

        print *, ' 2D '

        w = zeros(n1, n2)
        print *, w

        v = ones_like(w)
        print *, v

        print *, ' flip '
        do i = 1, n1
            w(i, :) = linspace(2.0, 5.0, n2)
        end do
        w = flip(w, [2])
        print *, w

        print *, ' pad '
        w = pad(w, [1, 1, 2, 2], ['edge', 'symm', 'edge', 'const'])
        print *, w

        print *, ' rot90 '
        w = rot90(w)
        print *, w

        print *, ' rescale '
        w = rescale(w, [0.0, 1.0])
        print *, w

        print *, ' binarize '
        w = binarize(w, 0.3, [0.0, 1.0])
        print *, w

        print *, ' tile '
        w = w(1:n1, 1:n2)
        w = tile(w, [1, 1])
        print *, w

    end block

    block

        real, allocatable, dimension(:, :, :) :: w, v
        integer :: n1, n2, n3
        integer :: i, j

        n1 = 3
        n2 = 4
        n3 = 5

        print *, ' 3D '

        w = zeros(n1, n2, n3)
        print *, w

        v = ones_like(w)
        print *, v

        print *, ' flip '
        do i = 1, n1
            do j = 1, n2
                w(i, j, :) = linspace(1.0, 7.0, n3)
            end do
        end do
        w = flip(w, [2, 1, 3])
        print *, w

        print *, ' pad '
        w = pad(w, [1, 1, 2, 2, 0, 2], ['edge', 'symm', 'edge', 'const', 'symm', 'edge'])
        print *, w

        print *, ' rot90 '
        w = rot90(w, dim=3)
        print *, w

        print *, ' rescale '
        w = rescale(w, [0.0, 1.0])
        print *, w

        print *, ' binarize '
        w = binarize(w, 0.3, [0.0, 1.0])
        print *, w

        print *, ' tile '
        w = w(1:n1, 1:n2, 1:n3)
        w = tile(w, [1, 2, 3])
        print *, shape(w)
        print *, w

        print *, ' permute'
        print *, shape(w)
        w = permute(w, order=231)
        print *, shape(w)
        print *, w

    end block

end program test
