
program test

    use libflit

    implicit none

    real, allocatable, dimension(:) :: w

    w = [0, 4, 16, 36, 64, 100]

    print *, integ(w)
    
    print *, cumsum(w)

end program test
