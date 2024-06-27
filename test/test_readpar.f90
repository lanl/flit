
program test

    use libflit

    block

        character(len=1024) :: file_parameter
        real :: par
        real, allocatable, dimension(:) :: pars

        file_parameter = './tmp.txt'

        print *, 'readpar -- float'

        call readpar_float(file_parameter, 'par1', par, 0.0)
        print *, par

        call readpar_xfloat(file_parameter, 'par2', par, 0.0, 3.0)
        print *, par

        call readpar_nfloat(file_parameter, 'par3', pars, [1.0, 2.0, 3.0])
        print *, pars

    end block

    block

        character(len=1024) :: file_parameter
        integer :: par
        integer, allocatable, dimension(:) :: pars

        file_parameter = './tmp.txt'

        print *, 'readpar -- integer '

        call readpar_int(file_parameter, 'ipar1', par, 0)
        print *, par

        call readpar_xint(file_parameter, 'ipar2', par, 0, 3.0)
        print *, par

        call readpar_nint(file_parameter, 'ipar3', pars, [1, 2, 3])
        print *, pars

    end block

    block

        character(len=1024) :: file_parameter
        complex :: par
        complex, allocatable, dimension(:) :: pars

        file_parameter = './tmp.txt'

        print *, 'readpar -- complex'

        call readpar_complex(file_parameter, 'zcpar1', par, cmplx(0.0, 0.0))
        print *, par

        call readpar_xcomplex(file_parameter, 'zcpar2', par, cmplx(0.0, 0.0), 4.0)
        print *, par

        call readpar_ncomplex(file_parameter, 'zcpar3', pars, cmplx([1.0, 2.0, 3.0]))
        print *, pars

    end block

    block

        real :: par
        real, allocatable, dimension(:) :: pars

        print *, 'getpar -- float '

        call getpar_float('par1', par, 0.0)
        print *, par

        call getpar_xfloat('par2', par, 0.0, 1.0)
        print *, par

        call getpar_nfloat('par3', pars, [1.0, 2.0, 3.0])
        print *, pars

    end block

    block

        complex :: par
        complex, allocatable, dimension(:) :: pars

        print *, 'getpar -- complex'

        call getpar_complex('cpar1', par, cmplx(0.0, 0.0))
        print *, par

        call getpar_xcomplex('cpar2', par, cmplx(0.0, 0.0), 4.0)
        print *, par

        call getpar_ncomplex('cpar3', pars, cmplx([1.0, 2.0, 3.0]))
        print *, pars

    end block

end program test
