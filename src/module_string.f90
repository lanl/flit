!
! © 2024. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by
! Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by
! Triad National Security, LLC, and the U.S. Department of Energy/National
! Nuclear Security Administration. The Government is granted for itself and
! others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare. derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!

module libflit_string

    use iso_fortran_env

    implicit none

    private :: int1_to_string
    private :: int2_to_string
    private :: int4_to_string
    private :: int8_to_string
    private :: float_to_string
    private :: double_to_string
    private :: complex_to_string
    private :: dcomplex_to_string
    private :: bool_to_string

    private :: int2array_to_string
    private :: int4array_to_string
    private :: floatarray_to_string
    private :: doublearray_to_string
    private :: complexarray_to_string
    private :: dcomplexarray_to_string
    private :: stringarray_to_string
    private :: boolarray_to_string

    interface num2str
        module procedure :: int1_to_string
        module procedure :: int2_to_string
        module procedure :: int4_to_string
        module procedure :: int8_to_string
        module procedure :: float_to_string
        module procedure :: double_to_string
        module procedure :: complex_to_string
        module procedure :: dcomplex_to_string
        module procedure :: bool_to_string
    end interface num2str

    interface array_to_string
        module procedure :: int2array_to_string
        module procedure :: int4array_to_string
        module procedure :: floatarray_to_string
        module procedure :: doublearray_to_string
        module procedure :: complexarray_to_string
        module procedure :: dcomplexarray_to_string
        module procedure :: stringarray_to_string
        module procedure :: boolarray_to_string
    end interface array_to_string

    interface extract_nint
        module procedure :: extract_nint2
        module procedure :: extract_nint4
        module procedure :: extract_nint8
    end interface extract_nint

    interface extract_xint
        module procedure :: extract_xint2
        module procedure :: extract_xint4
        module procedure :: extract_xint8
    end interface extract_xint

contains

    !
    !> Convert logical to string
    !
    function bool_to_string(bool, format) result(str)

        logical, intent(in) :: bool
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: boolstr
        character(len=24) :: numformat

        if (present(format)) then
            numformat = format
        else
            numformat = 'simple'
        end if

        select case (numformat)
            case ('full')
                if (bool) then
                    boolstr = '.true.'
                else
                    boolstr = '.false.'
                end if
            case ('simple')
                if (bool) then
                    boolstr = 'y'
                else
                    boolstr = 'n'
                end if
        end select

        allocate (character(len=len_trim(boolstr)) :: str)
        str = trim(adjustl(boolstr))

    end function bool_to_string

    !
    !> Convert integer(kind=2) to string
    !
    function int1_to_string(num, format) result(str)

        integer(1), intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat, numstr

        if (present(format)) then
            numformat = format
        else
            numformat = '(i)'
        end if

        write (numstr, numformat) num

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function int1_to_string

    !
    !> Convert integer(kind=2) to string
    !
    function int2_to_string(num, format) result(str)

        integer(kind=2), intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat, numstr

        if (present(format)) then
            numformat = format
        else
            numformat = '(i)'
        end if

        write (numstr, numformat) num

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function int2_to_string

    !
    !> Convert integer(kind=4) to string
    !
    function int4_to_string(num, format) result(str)

        integer(kind=4), intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat, numstr

        if (present(format)) then
            numformat = format
        else
            numformat = '(i)'
        end if

        write (numstr, numformat) num

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function int4_to_string

    !
    !> Convert integer(kind=8) to string
    !
    function int8_to_string(num, format) result(str)

        integer(kind=8), intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat, numstr

        if (present(format)) then
            numformat = format
        else
            numformat = '(i)'
        end if

        write (numstr, numformat) num

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function int8_to_string

    !
    !> Convert real to string
    !
    function float_to_string(num, format) result(str)

        real, intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat, numstr

        if (present(format)) then
            numformat = format
        else
            numformat = '(f)'
        end if

        write (numstr, numformat) num

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function float_to_string

    !
    !> Convert double to string
    !
    function double_to_string(num, format) result(str)

        double precision, intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat, numstr

        if (present(format)) then
            numformat = format
        else
            numformat = '(es)'
        end if

        write (numstr, numformat) num

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function double_to_string

    !
    !> Convert complex to string
    !
    function complex_to_string(num, format) result(str)

        complex, intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat, numstr

        if (present(format)) then
            numformat = format
        else
            numformat = '(a)'
        end if

        write (numstr, numformat) float_to_string(real(num))//'+'//float_to_string(imag(num))//'i'

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function complex_to_string

    !
    !> Convert complex to string
    !
    function dcomplex_to_string(num, format) result(str)

        double complex, intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat, numstr

        if (present(format)) then
            numformat = format
        else
            numformat = '(a)'
        end if

        write (numstr, numformat) double_to_string(real(num))//'+'//double_to_string(imag(num))//'i'

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function dcomplex_to_string

    !
    !> Convert integer(kind=2) array to string, separated by coma (,)
    !
    function int2array_to_string(num, format) result(str)

        integer(kind=2), dimension(:), intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat
        character(len=size(num)*12) :: numstr
        integer :: np, i

        np = size(num)

        if (present(format)) then
            numformat = format
        else
            numformat = '(i)'
        end if

        numstr = ''
        do i = 1, np
            numstr = tidy(numstr)//int2_to_string(num(i), numformat)
            if (i < np) then
                numstr = tidy(numstr)//','
            end if
        end do

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function int2array_to_string

    !
    !> Convert integer(kind=4) array to string, separated by coma (,)
    !
    function int4array_to_string(num, format) result(str)

        integer(kind=4), dimension(:), intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat
        character(len=size(num)*12) :: numstr
        integer :: np, i

        np = size(num)

        if (present(format)) then
            numformat = format
        else
            numformat = '(i)'
        end if

        numstr = ''
        do i = 1, np
            numstr = tidy(numstr)//int4_to_string(num(i), numformat)
            if (i < np) then
                numstr = tidy(numstr)//','
            end if
        end do

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function int4array_to_string

    !
    !> Convert float array to string, separated by coma (,)
    !
    function floatarray_to_string(num, format) result(str)

        real, dimension(:), intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat
        character(len=size(num)*24) :: numstr
        integer :: np, i

        np = size(num)

        if (present(format)) then
            numformat = format
        else
            numformat = '(f)'
        end if

        numstr = ''
        do i = 1, np
            numstr = tidy(numstr)//float_to_string(num(i), numformat)
            if (i < np) then
                numstr = tidy(numstr)//','
            end if
        end do

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function floatarray_to_string

    !
    !> Convert double precision array to string, separated by coma (,)
    !
    function doublearray_to_string(num, format) result(str)

        double precision, dimension(:), intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat
        character(len=size(num)*24) :: numstr
        integer :: np, i

        np = size(num)

        if (present(format)) then
            numformat = format
        else
            numformat = '(es)'
        end if

        numstr = ''
        do i = 1, np
            numstr = tidy(numstr)//double_to_string(num(i), numformat)
            if (i < np) then
                numstr = tidy(numstr)//','
            end if
        end do

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function doublearray_to_string

    !
    !> Convert complex array to string, separated by coma (,)
    !
    function complexarray_to_string(num, format) result(str)

        complex, dimension(:), intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat
        character(len=size(num)*24) :: numstr
        integer :: np, i

        np = size(num)

        if (present(format)) then
            numformat = format
        else
            numformat = '(a)'
        end if

        numstr = ''
        do i = 1, np
            numstr = tidy(numstr)//complex_to_string(num(i), numformat)
            if (i < np) then
                numstr = tidy(numstr)//','
            end if
        end do

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function complexarray_to_string

    !
    !> Convert double complex array to string, separated by coma (,)
    !
    function dcomplexarray_to_string(num, format) result(str)

        double complex, dimension(:), intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat
        character(len=size(num)*24) :: numstr
        integer :: np, i

        np = size(num)

        if (present(format)) then
            numformat = format
        else
            numformat = '(a)'
        end if

        numstr = ''
        do i = 1, np
            numstr = tidy(numstr)//dcomplex_to_string(num(i), numformat)
            if (i < np) then
                numstr = tidy(numstr)//','
            end if
        end do

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function dcomplexarray_to_string

    !
    !> Convert string array to string, separated by coma (,)
    !
    function stringarray_to_string(num) result(str)

        character(len=*), dimension(:), intent(in) :: num
        character(len=:), allocatable :: str

        character(len=size(num)*24) :: numstr
        integer :: np, i

        np = size(num)

        numstr = ''
        do i = 1, np
            numstr = tidy(numstr)//tidy(num(i))
            if (i < np) then
                numstr = tidy(numstr)//','
            end if
        end do

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function stringarray_to_string

    !
    !> Convert logical array to string, separated by coma (,)
    !
    function boolarray_to_string(num, format) result(str)

        logical, dimension(:), intent(in) :: num
        character(len=*), intent(in), optional :: format
        character(len=:), allocatable :: str

        character(len=64) :: numformat
        character(len=size(num)*24) :: numstr
        integer :: np, i

        np = size(num)

        if (present(format)) then
            numformat = format
        else
            numformat = 'simple'
        end if

        numstr = ''
        do i = 1, np
            numstr = tidy(numstr)//bool_to_string(num(i), numformat)
            if (i < np) then
                numstr = tidy(numstr)//','
            end if
        end do

        allocate (character(len=len_trim(numstr)) :: str)
        str = trim(adjustl(numstr))

    end function boolarray_to_string

    !
    !> Remove leading spaces and tabs, and tailing spaces, from a string
    !
    function tidy(str) result(w)

        character(len=*), intent(in) :: str

        character(len=:), allocatable :: t, w
        integer :: i

        ! allocate memory
        allocate (character(len=len_trim(adjustl(str))) :: t)
        t = trim(adjustl(str))

        ! first remove the leading tabs, achar(9) is for tab
        i = index(t, achar(9))
        do while (i == 1)
            t = t(2:)
            i = index(t, achar(9))
        end do

        ! then extract the tidy string
        allocate (character(len=len_trim(adjustl(t))) :: w)
        w = trim(adjustl(t))

    end function tidy

    !
    !> Count number (and optionally positions) of a substring in a string
    !
    subroutine count_substring(s1, s2, c, s2_pos)

        character(len=*), intent(in) :: s1, s2
        integer, intent(out) :: c
        integer, allocatable, dimension(:), intent(inout), optional :: s2_pos

        integer :: p, posn

        c = 0
        if (len(s2) == 0) then
            return
        end if

        p = 1
        do
            posn = index(s1(p:), s2)
            if (posn == 0) then
                exit
            end if
            c = c + 1
            p = p + posn + len(s2)
        end do

        if (present(s2_pos)) then

            if (allocated(s2_pos)) then
                deallocate (s2_pos)
            end if
            allocate (s2_pos(1:c))

            c = 0
            p = 1
            do
                posn = index(s1(p:), s2)
                if (posn == 0) then
                    exit
                end if
                c = c + 1
                p = p + posn + len(s2)
                s2_pos(c) = p - len(s2) - 1
            end do

        end if

    end subroutine count_substring

    !
    !> Count number of words in string
    !
    function count_words(str)

        character(len=*), intent(in) :: str
        integer :: count_words

        call count_substring(tidy(str), ' ', count_words)
        count_words = count_words + 1

    end function count_words

    !
    !> Convert string from lowercase to uppercase
    !
    !>        http://www.star.le.ac.uk/~cgp/fortran.html
    !
    function to_upper(strIn) result(strOut)

        character(len=*), intent(in) :: strIn
        character(len=len(strIn)) :: strOut
        integer :: i, j

        do i = 1, len(strIn)
            j = iachar(strIn(i:i))
            if (j >= iachar("a") .and. j <= iachar("z")) then
                strOut(i:i) = achar(iachar(strIn(i:i)) - 32)
            else
                strOut(i:i) = strIn(i:i)
            end if
        end do

    end function to_upper

    !
    !> Convert string from uppercase to lower case
    !
    !>        http://www.star.le.ac.uk/~cgp/fortran.html
    !
    function to_lower(strIn) result(strOut)

        character(len=*), intent(in) :: strIn
        character(len=len(strIn)) :: strOut
        integer :: i, j

        do i = 1, len(strIn)
            j = iachar(strIn(i:i))
            if (j >= iachar("A") .and. j <= iachar("Z")) then
                strOut(i:i) = achar(iachar(strIn(i:i)) + 32)
            else
                strOut(i:i) = strIn(i:i)
            end if
        end do

    end function to_lower

    !
    !> Extract a single-precision number from string
    !
    function extract_float(str_) result(f)

        character(len=*), intent(in) :: str_
        real :: f

        character(len=:), allocatable :: str
        double precision :: var
        integer :: var_int
        integer :: pos_dot, pos_les, pos_ues
        character(len=256) :: var_char

        str = trim(adjustl(str_))

        pos_dot = index(str, '.')
        pos_les = index(str, 'e')
        pos_ues = index(str, 'E')

        if (pos_dot == 0 .and. pos_les == 0 .and. pos_ues == 0) then
            ! If the input stream is an integer
            read (str, '(i)') var_int
            var = real(var_int)
        else if (pos_dot == 0 .and. pos_les /= 0 .and. pos_ues == 0) then
            ! If the input stream is in the form of xxxxexxxx,
            ! then add a .0 to the string to convert it to xxxx.0exxxx
            var_char = str(1:pos_les - 1)//'.0'//str(pos_les:len_trim(str))
            read (var_char, '(g)') var
        else if (pos_dot == 0 .and. pos_les == 0 .and. pos_ues /= 0) then
            ! If the input stream is in the form of xxxxExxxx,
            ! then add a .0 to the string to convert it to xxxx.0Exxxx
            var_char = str(1:pos_ues - 1)//'.0'//str(pos_ues:len_trim(str))
            read (var_char, '(g)') var
        else
            read (str, '(g)') var
        end if

        f = real(var)

    end function extract_float

    !
    !> Extract a double-precision number from string
    !
    function extract_double(str_) result(f)

        character(len=*), intent(in) :: str_
        double precision :: f

        character(len=:), allocatable :: str
        double precision :: var
        integer :: var_int
        integer :: pos_dot, pos_les, pos_ues
        character(len=256) :: var_char

        str = trim(adjustl(str_))

        pos_dot = index(str, '.')
        pos_les = max(index(str, 'd'), index(str, 'e'))
        pos_ues = max(index(str, 'D'), index(str, 'E'))

        if (pos_dot == 0 .and. pos_les == 0 .and. pos_ues == 0) then
            ! If the input stream is an integer
            read (str, '(i)') var_int
            var = dble(var_int)
        else if (pos_dot == 0 .and. pos_les /= 0 .and. pos_ues == 0) then
            ! If the input stream is in the form of xxxxexxxx,
            ! then add a .0 to the string to convert it to xxxx.0exxxx
            var_char = str(1:pos_les - 1)//'.0'//str(pos_les:len_trim(str))
            read (var_char, '(g)') var
        else if (pos_dot == 0 .and. pos_les == 0 .and. pos_ues /= 0) then
            ! If the input stream is in the form of xxxxExxxx,
            ! then add a .0 to the string to convert it to xxxx.0Exxxx
            var_char = str(1:pos_ues - 1)//'.0'//str(pos_ues:len_trim(str))
            read (var_char, '(g)') var
        else
            read (str, '(g)') var
        end if

        f = var

    end function extract_double

    !
    !> Extract a 2-byte integer number from string
    !
    function extract_int2(str_)

        character(len=*), intent(in) :: str_
        integer(kind=2) :: extract_int2

        character(len=:), allocatable :: str

        str = trim(adjustl(str_))

        extract_int2 = int(extract_float(str), kind=2)

    end function extract_int2

    !
    !> Extract a 4-byte integer number from string
    !
    function extract_int4(str_)

        character(len=*), intent(in) :: str_
        integer(kind=4) :: extract_int4

        character(len=:), allocatable :: str

        str = trim(adjustl(str_))

        extract_int4 = int(extract_float(str), kind=4)

    end function extract_int4

    !
    !> Extract a 8-byte integer number from string
    !
    function extract_int8(str_)

        character(len=*), intent(in) :: str_
        integer(kind=8) :: extract_int8

        character(len=:), allocatable :: str

        str = trim(adjustl(str_))

        extract_int8 = int(extract_float(str), kind=8)

    end function extract_int8

    !
    !> Extract a complex number from string
    !
    function extract_complex(str_)

        character(len=*), intent(in) :: str_
        complex :: extract_complex

        character(len=:), allocatable :: str
        complex :: var
        real :: var_real = 0.0
        real :: var_imag = 0.0
        integer :: pos_i, pos_sign(1:2)
        integer, allocatable, dimension(:) :: pos_plus, pos_minus
        integer :: np, nm
        character(len=256) :: var_real_char, var_imag_char

        str = trim(adjustl(str_))

        pos_i = index(str, 'i')
        pos_sign = 0
        call count_substring(str, '+', np, pos_plus)
        call count_substring(str, '-', nm, pos_minus)
        pos_sign(1) = min(minval(pos_plus), minval(pos_minus))
        pos_sign(2) = max(maxval(pos_plus), maxval(pos_minus))
        where (pos_sign < 0 .or. pos_sign > len_trim(str)) pos_sign = 0

        if (pos_i == 0) then
            ! if pure real number
            var_real = extract_float(str)
            var = cmplx(var_real, 0.0)
        else

            if (pos_i == len_trim(str)) then
                ! imaginary part is at tail

                var_real_char = str(1:pos_sign(2) - 1)
                var_imag_char = str(max(pos_sign(2), 1):len_trim(str) - 1)

                if (var_real_char == var_imag_char) then
                    var_real_char = '0.0'
                end if

                select case (var_real_char)
                    case ('')
                        var_real = 0.0
                    case ('+')
                        var_real = 1.0
                    case ('-')
                        var_real = -1.0
                    case default
                        var_real = extract_float(var_real_char)
                end select

                select case (var_imag_char)
                    case ('', '+')
                        var_imag = 1.0
                    case ('-')
                        var_imag = -1.0
                    case default
                        var_imag = extract_float(var_imag_char)
                end select

            else
                ! if imaginary part is at head

                var_real_char = str(pos_i + 1:len_trim(str))
                var_imag_char = str(1:pos_i - 1)

                select case (var_real_char)
                    case ('')
                        var_real = 0.0
                    case ('+')
                        var_real = 1.0
                    case ('-')
                        var_real = -1.0
                    case default
                        var_real = extract_float(var_real_char)
                end select

                select case (var_imag_char)
                    case ('', '+')
                        var_imag = 1.0
                    case ('-')
                        var_imag = -1.0
                    case default
                        var_imag = extract_float(var_imag_char)
                end select
            end if

            var = cmplx(var_real, var_imag)
        end if

        extract_complex = var

    end function extract_complex

    !
    !> Extract a complex number from string
    !
    function extract_dcomplex(str_)

        character(len=*), intent(in) :: str_
        double complex :: extract_dcomplex

        character(len=:), allocatable :: str
        double complex :: var
        double precision :: var_real = 0.0
        double precision :: var_imag = 0.0
        integer :: pos_i, pos_sign(1:2)
        integer, allocatable, dimension(:) :: pos_plus, pos_minus
        integer :: np, nm
        character(len=256) :: var_real_char, var_imag_char

        str = trim(adjustl(str_))

        pos_i = index(str, 'i')
        pos_sign = 0
        call count_substring(str, '+', np, pos_plus)
        call count_substring(str, '-', nm, pos_minus)
        pos_sign(1) = min(minval(pos_plus), minval(pos_minus))
        pos_sign(2) = max(maxval(pos_plus), maxval(pos_minus))
        where (pos_sign < 0 .or. pos_sign > len_trim(str)) pos_sign = 0

        if (pos_i == 0) then
            ! if pure real number
            var_real = extract_double(str)
            var = dcmplx(var_real, 0.0d0)
        else

            if (pos_i == len_trim(str)) then
                ! imaginary part is at tail

                var_real_char = str(1:pos_sign(2) - 1)
                var_imag_char = str(max(pos_sign(2), 1):len_trim(str) - 1)

                if (var_real_char == var_imag_char) then
                    var_real_char = '0.0'
                end if

                select case (var_real_char)
                    case ('')
                        var_real = 0.0d0
                    case ('+')
                        var_real = 1.0d0
                    case ('-')
                        var_real = -1.0d0
                    case default
                        var_real = extract_double(var_real_char)
                end select

                select case (var_imag_char)
                    case ('', '+')
                        var_imag = 1.0d0
                    case ('-')
                        var_imag = -1.0d0
                    case default
                        var_imag = extract_double(var_imag_char)
                end select

            else
                ! if imaginary part is at head

                var_real_char = str(pos_i + 1:len_trim(str))
                var_imag_char = str(1:pos_i - 1)

                select case (var_real_char)
                    case ('')
                        var_real = 0.0d0
                    case ('+')
                        var_real = 1.0d0
                    case ('-')
                        var_real = -1.0d0
                    case default
                        var_real = extract_double(var_real_char)
                end select

                select case (var_imag_char)
                    case ('', '+')
                        var_imag = 1.0d0
                    case ('-')
                        var_imag = -1.0d0
                    case default
                        var_imag = extract_double(var_imag_char)
                end select
            end if

            var = dcmplx(var_real, var_imag)
        end if

        extract_dcomplex = var

    end function extract_dcomplex

    !
    !> Extract a logical from string
    !
    function extract_logical(str_)

        character(len=*), intent(in) :: str_
        logical :: extract_logical

        character(len=:), allocatable :: str

        str = trim(adjustl(str_))

        if (tidy(str) == 'yes' &
                .or. tidy(str) == 'YES' &
                .or. tidy(str) == 'y' &
                .or. tidy(str) == 'Y' &
                .or. tidy(str) == 'true' &
                .or. tidy(str) == 'TRUE' &
                .or. tidy(str) == '.true.' &
                .or. tidy(str) == '.TRUE.' &
                .or. tidy(str) == 't' &
                .or. tidy(str) == 'T' &
                .or. tidy(str) == '1') then
            extract_logical = .true.
        else
            extract_logical = .false.
        end if

    end function extract_logical

    !
    !> Extract a string from string, which is to handle cases containing ', ", of mixed of them
    !
    function extract_string(str)

        character(len=*), intent(in) :: str
        character(len=len(trim(adjustl(str)))) :: tmpstring
        character(len=len(trim(adjustl(str)))) :: extract_string

        tmpstring = tidy(str)
        if ((index(tmpstring, '"') == 1 &
                .and. index(tmpstring, '"', back=.true.) == len(trim(adjustl(str)))) &
                .or. &
                (index(tmpstring, "'") == 1 &
                .and. index(tmpstring, "'", back=.true.) == len(trim(adjustl(str))))) then
            extract_string = tidy(tmpstring(2:len_trim(tmpstring) - 1))
        else if (index(tmpstring, "'") == 0 .and. index(tmpstring, '"') == 0) then
            extract_string = tidy(tmpstring)
        else
            extract_string = tmpstring
        end if

    end function extract_string

    !
    !> Extract multiple integers of kind 2 from string
    !
    subroutine extract_nint2(str, separator, var)

        character(len=*), intent(in) :: str, separator
        integer(kind=2), allocatable, dimension(:), intent(inout) :: var

        integer :: i, n
        character(len=24), allocatable, dimension(:) :: vars

        vars = split_string(str, separator)
        n = size(vars)

        if (allocated(var)) then
            deallocate(var)
        end if
        allocate (var(1:n))
        do i = 1, n
            var(i) = int(extract_float(vars(i)), kind=2)
        end do

    end subroutine extract_nint2

    !
    !> Extract multiple integers of kind 4 from string
    !
    subroutine extract_nint4(str, separator, var)

        character(len=*), intent(in) :: str, separator
        integer(kind=4), allocatable, dimension(:), intent(inout) :: var

        integer :: i, n
        character(len=24), allocatable, dimension(:) :: vars

        vars = split_string(str, separator)
        n = size(vars)

        if (allocated(var)) then
            deallocate(var)
        end if
        allocate (var(1:n))
        do i = 1, n
            var(i) = int(extract_float(vars(i)), kind=4)
        end do

    end subroutine extract_nint4

    !
    !> Extract multiple integers of kind 8 from string
    !
    subroutine extract_nint8(str, separator, var)

        character(len=*), intent(in) :: str, separator
        integer(kind=8), allocatable, dimension(:), intent(inout) :: var

        integer :: i, n
        character(len=24), allocatable, dimension(:) :: vars

        vars = split_string(str, separator)
        n = size(vars)

        if (allocated(var)) then
            deallocate(var)
        end if
        allocate (var(1:n))
        do i = 1, n
            var(i) = int(extract_float(vars(i)), kind=8)
        end do

    end subroutine extract_nint8

    !
    !> Extract multiple single-precision float numbers from string
    !
    subroutine extract_nfloat(str, separator, var)

        character(len=*), intent(in) :: str, separator
        real, allocatable, dimension(:), intent(inout) :: var

        integer :: i, n
        character(len=24), allocatable, dimension(:) :: vars

        vars = split_string(str, separator)
        n = size(vars)

        if (allocated(var)) then
            deallocate(var)
        end if
        allocate (var(1:n))
        do i = 1, n
            var(i) = extract_float(vars(i))
        end do

    end subroutine extract_nfloat

    !
    !> Extract multiple double-precision float numbers from string
    !
    subroutine extract_ndouble(str, separator, var)

        character(len=*), intent(in) :: str, separator
        double precision, allocatable, dimension(:), intent(inout) :: var

        integer :: i, n
        character(len=24), allocatable, dimension(:) :: vars

        vars = split_string(str, separator)
        n = size(vars)

        if (allocated(var)) then
            deallocate(var)
        end if
        allocate (var(1:n))
        do i = 1, n
            var(i) = extract_double(vars(i))
        end do

    end subroutine extract_ndouble

    !
    !> Extract multiple single-precision complex numbers from string
    !
    subroutine extract_ncomplex(str, separator, var)

        character(len=*), intent(in) :: str, separator
        complex, allocatable, dimension(:), intent(inout) :: var

        integer :: i, n
        character(len=48), allocatable, dimension(:) :: vars

        vars = split_string(str, separator)
        n = size(vars)

        if (allocated(var)) then
            deallocate(var)
        end if
        allocate (var(1:n))
        do i = 1, n
            var(i) = extract_complex(vars(i))
        end do

    end subroutine extract_ncomplex

    !
    !> Extract multiple double-precision complex numbers from string
    !
    subroutine extract_ndcomplex(str, separator, var)

        character(len=*), intent(in) :: str, separator
        double complex, allocatable, dimension(:), intent(inout) :: var

        integer :: i, n
        character(len=48), allocatable, dimension(:) :: vars

        vars = split_string(str, separator)
        n = size(vars)

        if (allocated(var)) then
            deallocate(var)
        end if
        allocate (var(1:n))
        do i = 1, n
            var(i) = extract_dcomplex(vars(i))
        end do

    end subroutine extract_ndcomplex

    !
    !> Extract multiple strings from string
    !
    subroutine extract_nstring(str, separator, var)

        character(len=*), intent(in) :: str, separator
        character(len=*), allocatable, dimension(:), intent(inout) :: var

        var = split_string(str, separator)

    end subroutine extract_nstring

    !
    !> Extract multiple logical values from string
    !
    subroutine extract_nlogical(str, separator, var)

        character(len=*), intent(in) :: str, separator
        logical, allocatable, dimension(:), intent(inout) :: var

        integer :: i, n
        character(len=8), allocatable, dimension(:) :: vars

        vars = split_string(str, separator)
        n = size(vars)

        if (allocated(var)) then
            deallocate(var)
        end if
        allocate (var(1:n))
        do i = 1, n
            var(i) = extract_logical(vars(i))
        end do

    end subroutine extract_nlogical

    !
    !> Extract a~b:c from a string (where a, b are floats and c is int2)
    !
    subroutine extract_xint2(str, rbeg, rend, var)

        character(len=*), intent(in) :: str
        real, intent(inout) :: rbeg, rend
        integer(kind=2), intent(inout) :: var

        integer :: pos_tilde, pos_smcol

        pos_tilde = index(str, '~')
        pos_smcol = index(str, ':')

        if (pos_tilde /= 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_tilde - 1))
            rend = extract_float(str(pos_tilde + 1:pos_smcol - 1))
        else if (pos_tilde == 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_smcol - 1))
            rend = rbeg
        else
            rbeg = 0
            rend = 0
        end if

        var = extract_int2(str(pos_smcol + 1:))

    end subroutine extract_xint2

    !
    !> Extract a~b:c from a string (where a, b are floats and c is int4)
    !
    subroutine extract_xint4(str, rbeg, rend, var)

        character(len=*), intent(in) :: str
        real, intent(inout) :: rbeg, rend
        integer(kind=4), intent(inout) :: var

        integer :: pos_tilde, pos_smcol

        pos_tilde = index(str, '~')
        pos_smcol = index(str, ':')

        if (pos_tilde /= 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_tilde - 1))
            rend = extract_float(str(pos_tilde + 1:pos_smcol - 1))
        else if (pos_tilde == 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_smcol - 1))
            rend = rbeg
        else
            rbeg = 0
            rend = 0
        end if

        var = extract_int4(str(pos_smcol + 1:))

    end subroutine extract_xint4

    !
    !> Extract a~b:c from a string (where a, b are floats and c is int8)
    !
    subroutine extract_xint8(str, rbeg, rend, var)

        character(len=*), intent(in) :: str
        real, intent(inout) :: rbeg, rend
        integer(kind=8), intent(inout) :: var

        integer :: pos_tilde, pos_smcol

        pos_tilde = index(str, '~')
        pos_smcol = index(str, ':')

        if (pos_tilde /= 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_tilde - 1))
            rend = extract_float(str(pos_tilde + 1:pos_smcol - 1))
        else if (pos_tilde == 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_smcol - 1))
            rend = rbeg
        else
            rbeg = 0
            rend = 0
        end if

        var = extract_int8(str(pos_smcol + 1:))

    end subroutine extract_xint8

    !
    !> Extract a~b:c from a string (where a, b are floats and c is float)
    !
    subroutine extract_xfloat(str, rbeg, rend, var)

        character(len=*), intent(in) :: str
        real, intent(inout) :: rbeg, rend, var

        integer :: pos_tilde, pos_smcol

        pos_tilde = index(str, '~')
        pos_smcol = index(str, ':')

        if (pos_tilde /= 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_tilde - 1))
            rend = extract_float(str(pos_tilde + 1:pos_smcol - 1))
        else if (pos_tilde == 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_smcol - 1))
            rend = rbeg
        else
            rbeg = 0
            rend = 0
        end if

        var = extract_float(str(pos_smcol + 1:))

    end subroutine extract_xfloat

    !
    !> Extract a~b:c from a string (where a, b are floats and c is double)
    !
    subroutine extract_xdouble(str, rbeg, rend, var)

        character(len=*), intent(in) :: str
        real, intent(inout) :: rbeg, rend
        double precision, intent(inout) :: var

        integer :: pos_tilde, pos_smcol

        pos_tilde = index(str, '~')
        pos_smcol = index(str, ':')

        if (pos_tilde /= 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_tilde - 1))
            rend = extract_float(str(pos_tilde + 1:pos_smcol - 1))
        else if (pos_tilde == 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_smcol - 1))
            rend = rbeg
        else
            rbeg = 0
            rend = 0
        end if

        var = extract_double(str(pos_smcol + 1:))

    end subroutine extract_xdouble

    !
    !> Extract a~b:c from a string (where a, b are floats and c is single-precision complex)
    !
    subroutine extract_xcomplex(str, rbeg, rend, var)

        character(len=*), intent(in) :: str
        real, intent(inout) :: rbeg, rend
        complex, intent(inout) :: var

        integer :: pos_tilde, pos_smcol

        pos_tilde = index(str, '~')
        pos_smcol = index(str, ':')

        if (pos_tilde /= 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_tilde - 1))
            rend = extract_float(str(pos_tilde + 1:pos_smcol - 1))
        else if (pos_tilde == 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_smcol - 1))
            rend = rbeg
        else
            rbeg = 0
            rend = 0
        end if

        var = extract_complex(str(pos_smcol + 1:))

    end subroutine extract_xcomplex

    !
    !> Extract a~b:c from a string (where a, b are floats and c is double-precision complex)
    !
    subroutine extract_xdcomplex(str, rbeg, rend, var)

        character(len=*), intent(in) :: str
        real, intent(inout) :: rbeg, rend
        double complex, intent(inout) :: var

        integer :: pos_tilde, pos_smcol

        pos_tilde = index(str, '~')
        pos_smcol = index(str, ':')

        if (pos_tilde /= 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_tilde - 1))
            rend = extract_float(str(pos_tilde + 1:pos_smcol - 1))
        else if (pos_tilde == 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_smcol - 1))
            rend = rbeg
        else
            rbeg = 0
            rend = 0
        end if

        var = extract_dcomplex(str(pos_smcol + 1:))

    end subroutine extract_xdcomplex

    !
    !> Extract a~b:c from a string (where a, b are floats and c is string)
    !
    subroutine extract_xstring(str, rbeg, rend, var)

        character(len=*), intent(in) :: str
        real, intent(inout) :: rbeg, rend
        character(len=*), intent(inout) :: var

        integer :: pos_tilde, pos_smcol

        pos_tilde = index(str, '~')
        pos_smcol = index(str, ':')

        if (pos_tilde /= 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_tilde - 1))
            rend = extract_float(str(pos_tilde + 1:pos_smcol - 1))
        else if (pos_tilde == 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_smcol - 1))
            rend = rbeg
        else
            rbeg = 0
            rend = 0
        end if

        var = extract_string(str(pos_smcol + 1:))

    end subroutine extract_xstring

    !
    !> Extract a~b:c from a string (where a, b are floats and c is logical)
    !
    subroutine extract_xlogical(str, rbeg, rend, var)

        character(len=*), intent(in) :: str
        real, intent(inout) :: rbeg, rend
        logical, intent(inout) :: var

        integer :: pos_tilde, pos_smcol

        pos_tilde = index(str, '~')
        pos_smcol = index(str, ':')

        if (pos_tilde /= 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_tilde - 1))
            rend = extract_float(str(pos_tilde + 1:pos_smcol - 1))
        else if (pos_tilde == 0 .and. pos_smcol /= 0) then
            rbeg = extract_float(str(1:pos_smcol - 1))
            rend = rbeg
        else
            rbeg = 0
            rend = 0
        end if

        var = extract_logical(str(pos_smcol + 1:))

    end subroutine extract_xlogical

    !
    !> Remove substring tailing a substring
    !
    function remove_string_after(str, strbeg) result(w)

        character(len=*), intent(in) :: str
        character(len=*), dimension(:), intent(in) :: strbeg

        character(len=:), allocatable :: t, w
        integer :: i

        ! allocate memory
        allocate (character(len=len_trim(adjustl(str))) :: t)
        t = trim(adjustl(str))

        do i = 1, size(strbeg)
            ! Iterate through all given prefixes
            if (t(1:len_trim(adjustl(strbeg(i)))) == trim(adjustl(strbeg(i)))) then
                ! If the beginning of str is equal to the given prefix, then
                ! assign to result string and exit
                allocate (character(len=len_trim(adjustl(strbeg(i)))) :: w)
                w = trim(adjustl(strbeg(i)))
                return
            end if
        end do

        ! Otherwise, return the input str
        allocate (character(len=len_trim(adjustl(str))) :: w)
        w = trim(adjustl(str))

    end function remove_string_after

    !
    !> Get substring tailing a substring
    !
    function get_string_after(str, strbeg) result(w)

        character(len=*), intent(in) :: str
        character(len=*), dimension(:), intent(in) :: strbeg

        character(len=:), allocatable :: t, w
        integer :: i

        ! allocate memory
        allocate (character(len=len_trim(adjustl(str))) :: t)
        t = trim(adjustl(str))

        do i = 1, size(strbeg)
            ! Iterate through all given prefixes
            if (t(1:len_trim(adjustl(strbeg(i)))) == trim(adjustl(strbeg(i)))) then
                ! If the beginning of str is equal to the given prefix, then
                ! assign to result rest of the str and exit
                allocate (character(len=len(t) - len_trim(adjustl(strbeg(i)))) :: w)
                w = t(len_trim(adjustl(strbeg(i))) + 1:len(t))
                return
            end if
        end do

        ! Otherwise, return an empty string
        allocate (character(len=1) :: w)
        w = ''

    end function get_string_after

    !
    !> Get substring based on index of characters
    !
    function substring(str, imin, imax) result(w)

        character(len=*), intent(in) :: str

        character(len=:), allocatable :: t, w
        integer :: imin, imax

        allocate (character(len=len_trim(adjustl(str))) :: t)
        t = trim(adjustl(str))

        if (imax > len_trim(t)) then
            allocate (character(len=len_trim(t) - imin + 1) :: w)
            w = t(imin:)
        else
            allocate (character(len=imax - imin + 1) :: w)
            w = t(imin:imax)
        end if

    end function substring

    !
    !> Reverse a string
    !
    function reverse_string(str) result(s)

        character(len=*) :: str
        character(len=len(str)) :: s

        integer :: i, n

        n = len(str)
        do i = 1, n
            s(n - i + 1:n - i + 1) = str(i:i)
        end do

    end function reverse_string

    !
    !> Pad a string with spaces and place the string in the center
    !
    function center_substring(s, l) result(sr)

        character(len=*) :: s
        integer :: l
        character(len=:), allocatable :: sr

        integer :: i, n, nz

        n = len_trim(adjustl(s))
        if (n > l) then
            write (error_unit, *) ' <center_substring> Error: len(substring) must <= l. '
            stop
        end if
        nz = nint((l - n)/2.0)

        allocate (character(len=l) :: sr)
        do i = 1, l
            sr(i:i) = ' '
        end do

        sr(nz + 1:nz + n) = trim(adjustl(s))

    end function center_substring

    !
    !> Split a string with a separator
    !> Inside the string, the separator may be consecutively repeated,
    !> and in this case, consecutive separators will be treated as one.
    !
    function split_string(str, sep) result(parts)

        character(len=*), intent(in) :: str
        character(len=*), intent(in) :: sep
        character(len=:), allocatable, dimension(:) :: parts

        integer :: start, pos, sep_len, str_len
        integer :: word_count, i
        character(len=:), allocatable :: temp

        !---------------------------------------
        ! Step 0: Initialize
        !---------------------------------------
        str_len = len(str)
        sep_len = len(sep)
        word_count = 0
        start = 1

        !---------------------------------------
        ! Step 1: Count words (skip empty ones)
        !---------------------------------------
        do while (start <= str_len)
            pos = index(str(start:), sep)

            if (pos == 0) then
                temp = str(start:)
                if (len_trim(adjustl(temp)) > 0) then
                    word_count = word_count + 1
                end if
                exit
            else
                if (pos > 1) then
                    temp = str(start:start + pos - 2)
                    if (len_trim(adjustl(temp)) > 0) then
                        word_count = word_count + 1
                    end if
                end if
                start = start + pos - 1 + sep_len
                ! Skip consecutive separators
                do while (start <= str_len .and. str(start:start + sep_len - 1) == sep)
                    start = start + sep_len
                end do
            end if
        end do

        !---------------------------------------
        ! Step 2: Allocate array and fill it
        !---------------------------------------
        allocate (character(len=len_trim(adjustl(str))) :: parts(word_count))
        start = 1
        i = 0

        do while (start <= str_len)
            pos = index(str(start:), sep)

            if (pos == 0) then
                temp = str(start:)
                if (len_trim(temp) > 0) then
                    i = i + 1
                    parts(i) = trim(adjustl(temp))
                end if
                exit
            else
                if (pos > 1) then
                    temp = str(start:start + pos - 2)
                    if (len_trim(temp) > 0) then
                        i = i + 1
                        parts(i) = trim(adjustl(temp))
                    end if
                end if
                start = start + pos - 1 + sep_len
                do while (start <= str_len .and. str(start:start + sep_len - 1) == sep)
                    start = start + sep_len
                end do
            end if
        end do

    end function split_string

end module libflit_string


! History
!
! The following cannot properly process consecutive separators such as
! ,,, ::: or consecutive spaces; therefore move them to legacy.
!
!    !
!    !> Extract multiple integers of kind 2 from string
!    !
!    subroutine extract_nint2(str, separator, var)
!
!        character(len=*), intent(in) :: str, separator
!        integer(kind=2), allocatable, dimension(:), intent(inout) :: var
!
!        integer :: nsep, nvar, i
!        integer, allocatable, dimension(:) :: pos_sep
!        integer :: ibeg, iend
!
!        call count_substring(str, separator, nsep, pos_sep)
!        nvar = nsep + 1
!        if (allocated(var)) then
!            deallocate (var)
!        end if
!        allocate (var(1:nvar))
!
!        if (nvar == 1) then
!            var(1) = int(extract_float(str), kind=2)
!            return
!        end if
!
!        do i = 1, nvar
!
!            if (i == 1) then
!                ibeg = 1
!                iend = pos_sep(1) - 1
!            else if (i == nvar) then
!                ibeg = pos_sep(nsep) + 1
!                iend = len_trim(str)
!            else
!                ibeg = pos_sep(i - 1) + 1
!                iend = pos_sep(i) - 1
!            end if
!            var(i) = int(extract_float(str(ibeg:iend)), kind=2)
!
!        end do
!
!    end subroutine extract_nint2
!
!    !
!    !> Extract multiple integers of kind 4 from string
!    !
!    subroutine extract_nint4(str, separator, var)
!
!        character(len=*), intent(in) :: str, separator
!        integer(kind=4), allocatable, dimension(:), intent(inout) :: var
!
!        integer :: nsep, nvar, i
!        integer, allocatable, dimension(:) :: pos_sep
!        integer :: ibeg, iend
!
!        call count_substring(str, separator, nsep, pos_sep)
!        nvar = nsep + 1
!        if (allocated(var)) then
!            deallocate (var)
!        end if
!        allocate (var(1:nvar))
!
!        if (nvar == 1) then
!            var(1) = int(extract_float(str), kind=4)
!            return
!        end if
!
!        do i = 1, nvar
!
!            if (i == 1) then
!                ibeg = 1
!                iend = pos_sep(1) - 1
!            else if (i == nvar) then
!                ibeg = pos_sep(nsep) + 1
!                iend = len_trim(str)
!            else
!                ibeg = pos_sep(i - 1) + 1
!                iend = pos_sep(i) - 1
!            end if
!            var(i) = int(extract_float(str(ibeg:iend)), kind=4)
!
!        end do
!
!    end subroutine extract_nint4
!
!    subroutine extract_nint8(str, separator, var)
!
!        character(len=*), intent(in) :: str, separator
!        integer(kind=8), allocatable, dimension(:), intent(inout) :: var
!
!        integer :: nsep, nvar, i
!        integer, allocatable, dimension(:) :: pos_sep
!        integer :: ibeg, iend
!
!        call count_substring(str, separator, nsep, pos_sep)
!        nvar = nsep + 1
!        if (allocated(var)) then
!            deallocate (var)
!        end if
!        allocate (var(1:nvar))
!
!        if (nvar == 1) then
!            var(1) = int(extract_float(str), kind=8)
!            return
!        end if
!
!        ! !$omp parallel do private(i,ibeg,iend)
!        do i = 1, nvar
!
!            if (i == 1) then
!                ibeg = 1
!                iend = pos_sep(1) - 1
!            else if (i == nvar) then
!                ibeg = pos_sep(nsep) + 1
!                iend = len_trim(str)
!            else
!                ibeg = pos_sep(i - 1) + 1
!                iend = pos_sep(i) - 1
!            end if
!            var(i) = int(extract_float(str(ibeg:iend)), kind=8)
!
!        end do
!        ! !$omp end parallel do
!
!    end subroutine extract_nint8
!
!    !
!    !> Extract multiple real numbers from string
!    !
!    subroutine extract_nfloat(str, separator, var)
!
!        character(len=*), intent(in) :: str, separator
!        real, allocatable, dimension(:), intent(inout) :: var
!
!        integer :: nsep, nvar, i
!        integer, allocatable, dimension(:) :: pos_sep
!        integer :: ibeg, iend
!
!        call count_substring(str, separator, nsep, pos_sep)
!        nvar = nsep + 1
!        if (allocated(var)) then
!            deallocate (var)
!        end if
!        allocate (var(1:nvar))
!
!        if (nvar == 1) then
!            var(1) = extract_float(str)
!            return
!        end if
!
!        ! !$omp parallel do private(i,ibeg,iend)
!        do i = 1, nvar
!
!            if (i == 1) then
!                ibeg = 1
!                iend = pos_sep(1) - 1
!            else if (i == nvar) then
!                ibeg = pos_sep(nsep) + 1
!                iend = len_trim(str)
!            else
!                ibeg = pos_sep(i - 1) + 1
!                iend = pos_sep(i) - 1
!            end if
!            var(i) = extract_float(str(ibeg:iend))
!
!        end do
!        ! !$omp end parallel do
!
!    end subroutine extract_nfloat
!
!    subroutine extract_ndouble(str, separator, var)
!
!        character(len=*), intent(in) :: str, separator
!        double precision, allocatable, dimension(:), intent(inout) :: var
!
!        integer :: nsep, nvar, i
!        integer, allocatable, dimension(:) :: pos_sep
!        integer :: ibeg, iend
!
!        call count_substring(str, separator, nsep, pos_sep)
!        nvar = nsep + 1
!        if (allocated(var)) then
!            deallocate (var)
!        end if
!        allocate (var(1:nvar))
!
!        if (nvar == 1) then
!            var(1) = extract_double(str)
!            return
!        end if
!
!        ! !$omp parallel do private(i,ibeg,iend)
!        do i = 1, nvar
!
!            if (i == 1) then
!                ibeg = 1
!                iend = pos_sep(1) - 1
!            else if (i == nvar) then
!                ibeg = pos_sep(nsep) + 1
!                iend = len_trim(str)
!            else
!                ibeg = pos_sep(i - 1) + 1
!                iend = pos_sep(i) - 1
!            end if
!            var(i) = extract_double(str(ibeg:iend))
!
!        end do
!        ! !$omp end parallel do
!
!    end subroutine extract_ndouble
!
!    !
!    !> Extract multiple complex numbers from string
!    !
!    subroutine extract_ncomplex(str, separator, var)
!
!        character(len=*), intent(in) :: str, separator
!        complex, allocatable, dimension(:), intent(inout) :: var
!
!        integer :: nsep, nvar, i
!        integer, allocatable, dimension(:) :: pos_sep
!        integer :: ibeg, iend
!
!        call count_substring(str, separator, nsep, pos_sep)
!        nvar = nsep + 1
!        if (allocated(var)) then
!            deallocate (var)
!        end if
!        allocate (var(1:nvar))
!
!        if (nvar == 1) then
!            var(1) = extract_complex(str)
!            return
!        end if
!
!        do i = 1, nvar
!
!            if (i == 1) then
!                ibeg = 1
!                iend = pos_sep(1) - 1
!            else if (i == nvar) then
!                ibeg = pos_sep(nsep) + 1
!                iend = len_trim(str)
!            else
!                ibeg = pos_sep(i - 1) + 1
!                iend = pos_sep(i) - 1
!            end if
!            var(i) = extract_complex(str(ibeg:iend))
!
!        end do
!
!    end subroutine extract_ncomplex
!
!    !
!    !> Extract multiple complex numbers from string
!    !
!    subroutine extract_ndcomplex(str, separator, var)
!
!        character(len=*), intent(in) :: str, separator
!        double complex, allocatable, dimension(:), intent(inout) :: var
!
!        integer :: nsep, nvar, i
!        integer, allocatable, dimension(:) :: pos_sep
!        integer :: ibeg, iend
!
!        call count_substring(str, separator, nsep, pos_sep)
!        nvar = nsep + 1
!        if (allocated(var)) then
!            deallocate (var)
!        end if
!        allocate (var(1:nvar))
!
!        if (nvar == 1) then
!            var(1) = extract_dcomplex(str)
!            return
!        end if
!
!        ! !$omp parallel do private(i,ibeg,iend)
!        do i = 1, nvar
!
!            if (i == 1) then
!                ibeg = 1
!                iend = pos_sep(1) - 1
!            else if (i == nvar) then
!                ibeg = pos_sep(nsep) + 1
!                iend = len_trim(str)
!            else
!                ibeg = pos_sep(i - 1) + 1
!                iend = pos_sep(i) - 1
!            end if
!            var(i) = extract_dcomplex(str(ibeg:iend))
!
!        end do
!        ! !$omp end parallel do
!
!    end subroutine extract_ndcomplex
!
!    !
!    !> Extract multiple strings from string
!    !
!    subroutine extract_nstring(str, separator, var)
!
!        character(len=*), intent(in) :: str, separator
!        character(len=*), allocatable, dimension(:), intent(inout) :: var
!
!        integer :: nsep, nvar, i
!        integer, allocatable, dimension(:) :: pos_sep
!        integer :: ibeg, iend
!
!        call count_substring(str, separator, nsep, pos_sep)
!        nvar = nsep + 1
!        if (allocated(var)) then
!            deallocate (var)
!        end if
!        allocate (var(1:nvar))
!
!        if (nvar == 1) then
!            var(1) = extract_string(str)
!            return
!        end if
!
!        do i = 1, nvar
!
!            if (i == 1) then
!                ibeg = 1
!                iend = pos_sep(1) - 1
!            else if (i == nvar) then
!                ibeg = pos_sep(nsep) + 1
!                iend = len_trim(str)
!            else
!                ibeg = pos_sep(i - 1) + 1
!                iend = pos_sep(i) - 1
!            end if
!            var(i) = extract_string(str(ibeg:iend))
!
!        end do
!
!    end subroutine extract_nstring
!
!    !
!    !> Extract multiple strings from string
!    !
!    subroutine extract_nlogical(str, separator, var)
!
!        character(len=*), intent(in) :: str, separator
!        logical, allocatable, dimension(:), intent(inout) :: var
!
!        integer :: nsep, nvar, i
!        integer, allocatable, dimension(:) :: pos_sep
!        integer :: ibeg, iend
!
!        call count_substring(str, separator, nsep, pos_sep)
!        nvar = nsep + 1
!        if (allocated(var)) then
!            deallocate (var)
!        end if
!        allocate (var(1:nvar))
!
!        if (nvar == 1) then
!            var(1) = extract_logical(str)
!            return
!        end if
!
!        ! !$omp parallel do private(i,ibeg,iend)
!        do i = 1, nvar
!
!            if (i == 1) then
!                ibeg = 1
!                iend = pos_sep(1) - 1
!            else if (i == nvar) then
!                ibeg = pos_sep(nsep) + 1
!                iend = len_trim(str)
!            else
!                ibeg = pos_sep(i - 1) + 1
!                iend = pos_sep(i) - 1
!            end if
!            var(i) = extract_logical(str(ibeg:iend))
!
!        end do
!        ! !$omp end parallel do
!
!    end subroutine extract_nlogical
