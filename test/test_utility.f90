
program test

    use libflit

    print *, round(1.2222, 0.25)
    print *, nice(1.2222, 0.1)
    print *, nice(1.2222e-5, 0.25)

end program test
