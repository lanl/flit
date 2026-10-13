
program test

    use libflit
    
    ! The out-of-range integer test below runs this program again with overflow=<name>,
    ! where reading the parameter must stop with an error
    block

        character(len=8) :: name
        integer(kind=2) :: par2
        integer :: par4
        integer(kind=8) :: par8

        call getpar_string('overflow', name, '')
        select case (name)
            case ('opar1')
                call readpar_int('./tmp.txt', name, par2, 0_2)
            case ('opar2', 'opar3', 'opar4')
                call readpar_int('./tmp.txt', name, par4, 0)
            case ('opar5', 'opar6')
                call readpar_int('./tmp.txt', name, par8, 0_8)
        end select
        if (name /= '') then
            stop
        end if

    end block

    block

        character(len=1024) :: file_parameter
        character(len=24) :: par
        character(len=48), allocatable, dimension(:) :: pars
        character(len=3) :: x

        file_parameter = './tmp.txt'

        print *, repeat('-', 60)
        print *, 'readpar -- string'
        print *, repeat('-', 60)

        call readpar_string(file_parameter, 'spar1', par, '')
        print *, par, merge('PASS', 'FAIL', par == 'model.bin')

        call readpar_xstring(file_parameter, 'spar2', par, '', 3.0)
        print *, par, merge('PASS', 'FAIL', par == 'vp.bin')

        call readpar_nstring(file_parameter, 'spar3', pars, [''])
        print *, pars(1), merge('PASS', 'FAIL', pars(1) == 'vp')
        print *, pars(2), merge('PASS', 'FAIL', pars(2) == 'vs')
        print *, pars(3), merge('PASS', 'FAIL', pars(3) == 'rho')
        
        call extract_nstring('vp, vs, rho', ',', pars)
        print *, pars(1), merge('PASS', 'FAIL', pars(1) == 'vp')
        print *, pars(2), merge('PASS', 'FAIL', pars(2) == 'vs')
        print *, pars(3), merge('PASS', 'FAIL', pars(3) == 'rho')
        
        print *, any('vs' == pars), merge('PASS', 'FAIL', any('vs' == pars))
        
        par = 'vs'
        print *, any(par == pars), merge('PASS', 'FAIL', any(par == pars))
        
        x = 'vs'
        print *, any(x == pars), merge('PASS', 'FAIL', any(x == pars))

    end block

    block

        character(len=1024) :: file_parameter
        real :: par
        real, allocatable, dimension(:) :: pars

        file_parameter = './tmp.txt'

        print *, repeat('-', 60)
        print *, 'readpar -- float'
        print *, repeat('-', 60)

        call readpar_float(file_parameter, 'par1', par, 0.0)
        print *, par, merge('PASS', 'FAIL', abs(par - 0.2222) < 1.0e-6)

        call readpar_xfloat(file_parameter, 'par2', par, 0.0, 3.0)
        print *, par, merge('PASS', 'FAIL', abs(par - 0.3) < 1.0e-6)

        call readpar_nfloat(file_parameter, 'par3', pars, [1.0, 2.0, 3.0])
        print *, pars, merge('PASS', 'FAIL', all(abs(pars - [1.2, 2.3, 3.4]) < 1.0e-6))

    end block

    block

        character(len=1024) :: file_parameter
        integer :: par
        integer, allocatable, dimension(:) :: pars

        file_parameter = './tmp.txt'

        print *, repeat('-', 60)
        print *, 'readpar -- integer '
        print *, repeat('-', 60)

        call readpar_int(file_parameter, 'ipar1', par, 0)
        print *, par, merge('PASS', 'FAIL', par == 1234)

        call readpar_xint(file_parameter, 'ipar2', par, -1, 3.0)
        print *, par, merge('PASS', 'FAIL', par == 0)

        call readpar_nint(file_parameter, 'ipar3', pars, [1, 2, 3])
        print *, pars, merge('PASS', 'FAIL', all(pars == [2, 3, 4, 5, 6]))

    end block

    block

        character(len=1024) :: file_parameter
        complex :: par
        complex, allocatable, dimension(:) :: pars

        file_parameter = './tmp.txt'

        print *, repeat('-', 60)
        print *, 'readpar -- complex'
        print *, repeat('-', 60)

        call readpar_complex(file_parameter, 'zcpar1', par, cmplx(0.0, 0.0))
        print *, par, merge('PASS', 'FAIL', abs(par - cmplx(1.0, 2.0)) < 1.0e-6)

        call readpar_xcomplex(file_parameter, 'zcpar2', par, cmplx(0.0, 0.0), 4.0)
        print *, par, merge('PASS', 'FAIL', abs(par - cmplx(2.0, 3.0)/7.0) < 1.0e-6)

        call readpar_ncomplex(file_parameter, 'zcpar3', pars, cmplx([1.0, 2.0, 3.0]))
        print *, pars, merge('PASS', 'FAIL', all(abs(pars - [cmplx(1.2, 0.0), cmplx(2.3, 1.0), cmplx(3.4, -4.0)]) < 1.0e-6))

    end block

    block

        real :: par
        real, allocatable, dimension(:) :: pars
        logical :: given

        print *, repeat('-', 60)
        print *, 'getpar -- float '
        print *, repeat('-', 60)

        ! With the command-line arguments in test.sh, getpar must return their values;
        ! without arguments, it must return the default values
        given = command_argument_count() > 0
        if (.not. given) then
            print *, 'No command-line arguments, so getpar must return the default values'
        end if

        call getpar_float('par1', par, 0.0)
        print *, par, merge('PASS', 'FAIL', abs(par - merge(0.12345, 0.0, given)) < 1.0e-6)

        call getpar_xfloat('par2', par, 0.0, 1.0)
        print *, par, merge('PASS', 'FAIL', abs(par - merge(0.2, 0.0, given)) < 1.0e-6)

        call getpar_nfloat('par3', pars, [1.0, 2.0, 3.0])
        print *, pars, merge('PASS', 'FAIL', all(abs(pars - merge([1.2, 2.3, 3.4], [1.0, 2.0, 3.0], given)) < 1.0e-6))

    end block

    block

        complex :: par
        complex, allocatable, dimension(:) :: pars
        logical :: given

        print *, repeat('-', 60)
        print *, 'getpar -- complex'
        print *, repeat('-', 60)

        ! With the command-line arguments in test.sh, getpar must return their values;
        ! without arguments, it must return the default values
        given = command_argument_count() > 0

        call getpar_complex('cpar1', par, cmplx(0.0, 0.0))
        print *, par, merge('PASS', 'FAIL', abs(par - merge(cmplx(0.12345, 0.1), cmplx(0.0, 0.0), given)) < 1.0e-6)

        call getpar_xcomplex('cpar2', par, cmplx(0.0, 0.0), 4.0)
        print *, par, merge('PASS', 'FAIL', abs(par - merge(cmplx(2.0, 3.0)/7.0, cmplx(0.0, 0.0), given)) < 1.0e-6)

        call getpar_ncomplex('cpar3', pars, cmplx([1.0, 2.0, 3.0]))
        print *, pars, merge('PASS', 'FAIL', &
            all(abs(pars - merge([cmplx(1.2, 0.0), cmplx(2.3, 1.0), cmplx(3.4, -4.0)], cmplx([1.0, 2.0, 3.0]), given)) < 1.0e-6))

    end block

    block

        character(len=1024) :: file_parameter
        integer :: par
        integer, allocatable, dimension(:) :: pars
        integer(kind=8) :: par8
        integer(kind=8), allocatable, dimension(:) :: pars8

        file_parameter = './tmp.txt'

        print *, repeat('-', 60)
        print *, 'readpar -- large integer'
        print *, repeat('-', 60)

        call readpar_int(file_parameter, 'ipar4', par, 0)
        print *, par, merge('PASS', 'FAIL', par == 16777217)

        call readpar_nint(file_parameter, 'ipar5', pars, [0])
        print *, pars, merge('PASS', 'FAIL', all(pars == [123456789, 2000000001]))

        call readpar_int(file_parameter, 'ipar6', par8, 0_8)
        print *, par8, merge('PASS', 'FAIL', par8 == 3000000000_8)

        call readpar_nint(file_parameter, 'ipar7', pars, [0])
        print *, pars, merge('PASS', 'FAIL', all(pars == [huge(par), -huge(par) - 1]))

        call readpar_nint(file_parameter, 'ipar8', pars8, [0_8])
        print *, pars8, merge('PASS', 'FAIL', all(pars8 == [9007199254740993_8, huge(par8), -huge(par8) - 1_8]))

        call readpar_nint(file_parameter, 'ipar9', pars, [0])
        print *, pars, merge('PASS', 'FAIL', all(pars == [huge(par), -huge(par) - 1, 1000, 2]))

    end block

    block

        character(len=1024) :: exe
        character(len=8), dimension(6) :: name = [character(len=8) :: 'opar1', 'opar2', 'opar3', 'opar4', 'opar5', 'opar6']
        integer :: k, status

        print *, repeat('-', 60)
        print *, 'readpar -- out-of-range integer'
        print *, repeat('-', 60)

        ! Reading each value must stop with an error, so each read runs in a child process
        call get_command_argument(0, exe)
        do k = 1, size(name)
            call execute_command_line(trim(exe)//' overflow='//trim(name(k))//' > /dev/null 2>&1', exitstat=status)
            print *, trim(name(k)), status, merge('PASS', 'FAIL', status /= 0)
        end do

    end block

    block

        character(len=1024) :: file_parameter
        real :: par
        double precision :: dpar
        double precision, allocatable, dimension(:) :: dpars

        file_parameter = './tmp.txt'

        print *, repeat('-', 60)
        print *, 'readpar -- high-precision float'
        print *, repeat('-', 60)

        ! Each value must be the nearest float or double to the decimal string,
        ! the same as the compiler gives for the literal in the code
        call readpar_float(file_parameter, 'fpar1', par, 0.0)
        print '(1x, es16.8, 2x, a)', par, merge('PASS', 'FAIL', par == 1.0000001)

        call readpar_float(file_parameter, 'fpar2', par, 0.0)
        print '(1x, es16.8, 2x, a)', par, merge('PASS', 'FAIL', par == 0.9999999)

        call readpar_float(file_parameter, 'fpar3', par, 0.0)
        print '(1x, es16.8, 2x, a)', par, merge('PASS', 'FAIL', par == 1.0)

        call readpar_float(file_parameter, 'fpar4', par, 0.0)
        print '(1x, es16.8, 2x, a)', par, merge('PASS', 'FAIL', par == 1.00000006)

        call readpar_float(file_parameter, 'fpar5', par, 0.0)
        print '(1x, es16.8, 2x, a)', par, merge('PASS', 'FAIL', par == 3.0e9)

        call readpar_double(file_parameter, 'dpar1', dpar, 0.0d0)
        print '(1x, es24.16, 2x, a)', dpar, merge('PASS', 'FAIL', dpar == 1.0000000000000002d0)

        call readpar_double(file_parameter, 'dpar2', dpar, 0.0d0)
        print '(1x, es24.16, 2x, a)', dpar, merge('PASS', 'FAIL', dpar == 0.9999999999999999d0)

        call readpar_double(file_parameter, 'dpar3', dpar, 0.0d0)
        print '(1x, es24.16, 2x, a)', dpar, merge('PASS', 'FAIL', dpar == 1.0000001d0)

        call readpar_double(file_parameter, 'dpar4', dpar, 0.0d0)
        print '(1x, es24.16, 2x, a)', dpar, merge('PASS', 'FAIL', dpar == 1.2345678901234567890123d-10)

        call readpar_ndouble(file_parameter, 'dpar5', dpars, [0.0d0])
        print '(1x, 3es24.16, 2x, a)', dpars, &
            merge('PASS', 'FAIL', all(dpars == [1.0000000000000002d0, 0.9999999999999999d0, 1.2345678901234567890d-10]))

    end block

    block

        character(len=1024) :: file_parameter
        real, allocatable, dimension(:) :: pars

        file_parameter = './tmp.txt'

        print *, repeat('-', 60)
        print *, 'readpar -- default value'
        print *, repeat('-', 60)

        ! par5 is not in the file, so the allocated pars takes the default value
        call readpar_nfloat(file_parameter, 'par3', pars, [0.0])
        call readpar_nfloat(file_parameter, 'par5', pars, [9.0, 9.0], required=.false.)
        if (allocated(pars)) then
            print *, pars, merge('PASS', 'FAIL', all(pars == 9.0))
        else
            print *, 'pars is not allocated: FAIL'
        end if

    end block

    block

        character(len=1024) :: file_parameter
        real, allocatable, dimension(:) :: pars

        file_parameter = './tmp.txt'

        print *, repeat('-', 60)
        print *, 'readpar -- separator'
        print *, repeat('-', 60)

        call readpar_nfloat(file_parameter, 'par4', pars, [0.0], separator=' ')
        print *, pars, merge('PASS', 'FAIL', all(pars == [1.0, 2.0, 3.0]))

    end block

    block

        character(len=32), dimension(3) :: source
        real :: par
        real, allocatable, dimension(:) :: pars

        source = [character(len=32) :: 'par1 = 0.2222', 'par2 = 0:0, 10:1, 11~20:2', 'par3 = 1.2, 2.3, 3.4']

        print *, repeat('-', 60)
        print *, 'parsepar -- float'
        print *, repeat('-', 60)

        call parsepar_float(source, 'par1', par, 0.0)
        print *, par, merge('PASS', 'FAIL', abs(par - 0.2222) < 1.0e-6)

        call parsepar_xfloat(source, 'par2', par, 0.0, 3.0)
        print *, par, merge('PASS', 'FAIL', abs(par - 0.3) < 1.0e-6)

        call parsepar_nfloat(source, 'par3', pars, [1.0, 2.0, 3.0])
        print *, pars, merge('PASS', 'FAIL', all(abs(pars - [1.2, 2.3, 3.4]) < 1.0e-6))

    end block

    block

        character(len=1024) :: file_parameter
        character(len=8), allocatable, dimension(:) :: valid

        file_parameter = './tmp.txt'

        print *, repeat('-', 60)
        print *, 'checkpar'
        print *, repeat('-', 60)

        valid = [character(len=8) :: 'spar1', 'spar2', 'spar3', &
            'ipar1', 'ipar2', 'ipar3', 'ipar4', 'ipar5', 'ipar6', 'ipar7', 'ipar8', 'ipar9', &
            'par1', 'par2', 'par3', 'par4', 'fpar1', 'fpar2', 'fpar3', 'fpar4', 'fpar5', &
            'dpar1', 'dpar2', 'dpar3', 'dpar4', 'dpar5', 'opar1', 'opar2', 'opar3', 'opar4', 'opar5', 'opar6', &
            'zcpar1', 'zcpar2', 'zcpar3', 'zcpar1r']

        call checkpar(file_parameter, valid)
        print *, 'All parameters are valid: PASS'

        ! checkpar stops the program at an invalid parameter, so this test comes last
        print *, 'Without zcpar1r in the list, checkpar must report zcpar1r and stop:'
        call checkpar(file_parameter, valid(1:size(valid) - 1))
        print *, 'FAIL'

    end block

end program test
