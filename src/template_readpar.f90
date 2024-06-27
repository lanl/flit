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


#define PASTE(X)            X
#define PASTE2(X)           PASTE(X)_
#define CONCATHELP(X, Y)    PASTE2(X)Y
#define CONCAT(X, Y)        CONCATHELP(X, Y)

#define readpar_      CONCAT(readpar, T)
#define readnpar_      CONCAT(readnpar, T)
#define readxpar_      CONCAT(readxpar, T)
#define getpar_      CONCAT(getpar, T)
#define getnpar_      CONCAT(getnpar, T)
#define getxpar_      CONCAT(getxpar, T)
#define parsepar_      CONCAT(parsepar, T)
#define parsenpar_      CONCAT(parsenpar, T)
#define parsexpar_      CONCAT(parsexpar, T)

!
!> Read a parameter from file by name, e.g.,
!> parname = 10
!
subroutine readpar_(filename, parname, par, default_value, required)

    character(len=*), intent(in) :: filename, parname
    TT, intent(inout) :: par
    TT, intent(in) :: default_value
    logical, intent(in), optional :: required

    character(len=1024) :: line
    integer :: npar, funit, eqindex
    integer :: ioerr

    npar = 0
    open (newunit=funit, file=tidy(filename), status='old', action='read')
    do
        read (funit, '(a)', iostat=ioerr) line
        if (ioerr /= 0) exit
        if (tidy(line) == 'exit') exit
        if (len(tidy(line)) > 0) then
            eqindex = index(line, '=')
            if (to_lower(tidy(line(1:eqindex - 1))) == to_lower(tidy(parname))) then
                par = extract_TT(tidy(line(eqindex + 1:)))
                npar = npar + 1
            end if
        end if
    end do
    close (funit)

    if (npar == 0) then
        ! if not exists then set to default
        if (present(required)) then
            if (required) then
                ! if required then stop
                call warn(' Error: Parameter '//tidy(parname)//' is required but not set. ')
                stop
            else
                ! if not required then set to default
                par = default_value
            end if
        else
            ! if not required then set to default value
            par = default_value
        end if
    end if

end subroutine readpar_

!
!> Read a parameter from file by name, e.g.,
!> parname = 10, 20, 30, ...
!
subroutine readnpar_(filename, parname, par, default_value, required, separator)

    character(len=*), intent(in) :: filename, parname
    TT, allocatable, dimension(:), intent(inout) :: par
    TT, dimension(:), intent(in) :: default_value
    logical, intent(in), optional :: required
    character(len=*), intent(in), optional :: separator

    character(len=1024) :: line
    integer :: npar, funit, eqindex
    character(len=10) :: sepr
    integer :: ioerr

    if (present(separator)) then
        sepr = separator
    else
        sepr = ','
    end if

    npar = 0
    open (newunit=funit, file=tidy(filename), status='old', action='read')
    do
        read (funit, '(a)', iostat=ioerr) line
        if (ioerr /= 0) exit
        if (tidy(line) == 'exit') exit
        if (len(tidy(line)) > 0) then
            eqindex = index(line, '=')
            if (to_lower(tidy(line(1:eqindex - 1))) == to_lower(tidy(parname))) then
                call extract_nTT(tidy(line(eqindex + 1:)), trim(adjustl(sepr)), par)
                npar = npar + 1
            end if
        end if
    end do
    close (funit)

    if (npar == 0) then
        ! if not exists then set to default
        if (present(required)) then
            if (required) then
                ! if required then stop
                call warn(' Error: Parameter '//tidy(parname)//' is required but not set. ')
                stop
            else
                ! if not required then set to default
                if (allocated(par)) then
                    deallocate (par)
                else
                    par = default_value
                end if
            end if
        else
            ! if not required then set to default value
            if (allocated(par)) then
                deallocate (par)
            end if
            par = default_value
        end if
    end if

end subroutine readnpar_

!
!> Read a parameter from file using the format
!> parname = range1:value1, range2:value2, ...
!> e.g., parname = 1~10:0.0, 11~20:1.0, 21~30:2.0
!> where the variable ranges from 1~10, the parameter value is 10, ...
!
subroutine readxpar_(filename, parname, par, default_value, var, required)

    character(len=*), intent(in) :: filename, parname
    TT, intent(inout) :: par
    TT, intent(in) :: default_value
    real, intent(in) :: var
    logical, intent(in), optional :: required

    character(len=64), allocatable, dimension(:) :: tmpar
    real, allocatable, dimension(:) :: rbeg, rend
    TTT, allocatable, dimension(:) :: rpar
    character(len=1024) :: line
    integer :: npar, funit, eqindex, i, l, nr
    integer :: ioerr

    ! this initialization is necessary
    ! otherwise, segmentation fault
    allocate (tmpar(1:nstring_max))

    npar = 0
    open (newunit=funit, file=tidy(filename), status='old', action='read')
    do
        read (funit, '(a)', iostat=ioerr) line
        if (ioerr /= 0) exit
        if (tidy(line) == 'exit') exit
        if (len(tidy(line)) > 0) then
            eqindex = index(line, '=')
            if (to_lower(tidy(line(1:eqindex - 1))) == to_lower(tidy(parname))) then
                call extract_nstring(tidy(line(eqindex + 1:)), ',', tmpar)
                npar = npar + 1
            end if
        end if
    end do
    close (funit)

    if (npar >= 1) then

        nr = size(tmpar)

        ! extracting range and values
        allocate (rbeg(1:nr))
        allocate (rend(1:nr))
        allocate (rpar(1:nr))

        do i = 1, nr
            call extract_xTT(tidy(tmpar(i)), rbeg(i), rend(i), rpar(i))
        end do

        ! find target value correponding to var
        if (var < rbeg(1)) then
            ! if var < min of given range
            par = rpar(1)

        else if (var > rend(nr)) then
            ! if var > max of given range
            par = rpar(nr)

        else
            ! if var between min and max of given range
            ! and in some given interval
            l = 0
            do i = 1, nr
                if (var >= rbeg(i) .and. var <= rend(i)) then
                    par = rpar(i)
                    l = l + 1
                end if
            end do

            ! if var is located somewhere between some non-defined interval [rend, rbeg_next]
            ! then linearly interpolate to obtain par value
            if (l == 0) then

#ifdef linear_interp
                ! Linear interpolation applies to float, double, complex
                do i = 1, nr - 1
                    if (var > rend(i) .and. var < rbeg(i + 1)) then
                        par = rpar(i) + (rpar(i + 1) - rpar(i))*(var - rend(i))/(rbeg(i + 1) - rend(i))
                        exit
                    end if
                end do
#endif
#ifdef nearest_interp
                ! Linear interpolation does not apply to int, logical, or string
                do i = 1, nr - 1
                    if (var > rend(i) .and. var < rend(i) + 0.5*(rbeg(i + 1) - rend(i))) then
                        par = rpar(i)
                        exit
                    else if (var >= rend(i) + 0.5*(rbeg(i + 1) - rend(i)) .and. var < rbeg(i + 1)) then
                        par = rpar(i + 1)
                        exit
                    end if
                end do
#endif
            end if

        end if

    else
        ! if not exists then set to default
        if (present(required)) then
            if (required) then
                ! if required then stop
                call warn(' Error: Parameter '//tidy(parname)//' is required but not set. ')
                stop
            else
                ! if not required then set to default
                par = default_value
            end if
        else
            ! if not required then set to default value
            par = default_value
        end if

    end if

end subroutine readxpar_

!
!> Read a parameter from command line arguments by name, e.g.,
!> parname = 10
!
subroutine getpar_(parname, par, default_value, required)

    character(len=*), intent(in) :: parname
    TT, intent(out) :: par
    TT, intent(in) :: default_value
    logical, intent(in), optional :: required

    integer :: i
    integer :: len_arg
    integer :: len_par
    logical :: getpar_success
    character(len=256) :: arg

    getpar_success = .false.
    len_par = len_trim(parname)

    do i = 1, command_argument_count()

        call get_command_argument(i, arg)
        len_arg = len_trim(arg)

        if (len_arg > len_par + 1) then
            ! For getpar, the parname must be followed by = without any space.
            if (to_lower(tidy(arg(1:len_par + 1))) == to_lower(tidy(parname))//'=') then
                par = extract_TT(arg(len_par + 2:len_arg))
                getpar_success = .true.
            end if
        end if

    end do

    if (present(required) .and. required .and. (.not. getpar_success)) then
        write (error_unit, *) ' Error: Parameter ', trim(parname), ' is required but not set. '
        stop
    end if

    if ((.not. present(required) .or. (present(required) .and. (.not. required))) &
            .and. (.not. getpar_success)) then
        par = default_value
    end if

end subroutine getpar_

!
!> Read a parameter from command line arguments by name, e.g.,
!> parname = 10, 20, 30, ...
!
subroutine getnpar_(parname, par, default_value, separator, required)

    character(len=*), intent(in) :: parname
    TT, allocatable, dimension(:), intent(out) :: par
    TT, dimension(:), intent(in) :: default_value
    character(len=*), intent(in), optional :: separator
    logical, intent(in), optional :: required

    integer :: i
    integer :: len_arg
    integer :: len_par
    logical :: getpar_success
    character(len=256) :: arg
    character(len=10) :: sepr

    getpar_success = .false.
    len_par = len_trim(parname)
    if (present(separator)) then
        sepr = separator
    else
        sepr = ','
    end if

    do i = 1, command_argument_count()

        call get_command_argument(i, arg)
        len_arg = len_trim(arg)

        if (len_arg > len_par + 1) then
            ! For getpar, the parname must be followed by = without any space.
            if (to_lower(tidy(arg(1:len_par + 1))) == to_lower(tidy(parname))//'=') then
                call extract_nTT(arg(len_par + 2:len_arg), trim(adjustl(sepr)), par)
                getpar_success = .true.
            end if
        end if

    end do

    if (present(required) .and. required .and. (.not. getpar_success)) then
        write (error_unit, *) ' Error: Parameter ', trim(parname), ' is required but not set. '
        stop
    end if

    if ((.not. present(required) .or. (present(required) .and. (.not. required))) &
            .and. (.not. getpar_success)) then
        i = size(default_value, 1)
        if (allocated(par)) then
            deallocate (par)
        end if
        allocate (par(1:i))
        par = default_value
    end if

end subroutine getnpar_

!
!> Read a parameter from command line arguments using the format
!> parname=range1:value1,range2:value2,...
!> e.g., par=1~10:0.0,11~20:1.0,21~30:2.0
!> where the variable ranges from 1~10, the parameter value is 10, ...
!
subroutine getxpar_(parname, par, default_value, var, required)

    character(len=*), intent(in) :: parname
    TT, intent(inout) :: par
    TT, intent(in) :: default_value
    real, intent(in) :: var
    logical, intent(in), optional :: required

    character(len=64), allocatable, dimension(:) :: tmpar
    real, allocatable, dimension(:) :: rbeg, rend
    TTT, allocatable, dimension(:) :: rpar
    integer :: i, l, nr
    integer :: len_arg
    integer :: len_par
    logical :: getpar_success
    character(len=256) :: arg

    ! this initialization is necessary
    ! otherwise, segmentation fault
    allocate (tmpar(1:nstring_max))

    getpar_success = .false.
    len_par = len_trim(parname)

    do i = 1, command_argument_count()

        call get_command_argument(i, arg)
        len_arg = len_trim(arg)

        if (len_arg > len_par + 1) then
            ! For getpar, the parname must be followed by = without any space.
            if (to_lower(tidy(arg(1:len_par + 1))) == to_lower(tidy(parname))//'=') then
                call extract_nstring(tidy(arg(len_par + 2:)), ',', tmpar)
                getpar_success = .true.
            end if
        end if

    end do

    if (getpar_success) then

        nr = size(tmpar)

        ! extracting range and values
        allocate (rbeg(1:nr))
        allocate (rend(1:nr))
        allocate (rpar(1:nr))

        do i = 1, nr
            call extract_xTT(tidy(tmpar(i)), rbeg(i), rend(i), rpar(i))
        end do

        ! find target value correponding to var
        if (var < rbeg(1)) then
            ! if var < min of given range
            par = rpar(1)

        else if (var > rend(nr)) then
            ! if var > max of given range
            par = rpar(nr)

        else
            ! if var between min and max of given range
            ! and in some given interval
            l = 0
            do i = 1, nr
                if (var >= rbeg(i) .and. var <= rend(i)) then
                    par = rpar(i)
                    l = l + 1
                end if
            end do

            ! if var is located somewhere between some non-defined interval [rend, rbeg_next]
            ! then linearly interpolate to obtain par value
            if (l == 0) then

#ifdef linear_interp
                ! Linear interpolation applies to float, double, complex
                do i = 1, nr - 1
                    if (var > rend(i) .and. var < rbeg(i + 1)) then
                        par = rpar(i) + (rpar(i + 1) - rpar(i))*(var - rend(i))/(rbeg(i + 1) - rend(i))
                        exit
                    end if
                end do
#endif
#ifdef nearest_interp
                ! Linear interpolation does not apply to int, logical, or string
                do i = 1, nr - 1
                    if (var > rend(i) .and. var < rend(i) + 0.5*(rbeg(i + 1) - rend(i))) then
                        par = rpar(i)
                        exit
                    else if (var >= rend(i) + 0.5*(rbeg(i + 1) - rend(i)) .and. var < rbeg(i + 1)) then
                        par = rpar(i + 1)
                        exit
                    end if
                end do
#endif
            end if

        end if

    else

        if (present(required) .and. required) then
            write (error_unit, *) ' Error: Parameter ', trim(parname), ' is required but not set. '
            stop
        end if

        if (.not. present(required) .or. (present(required) .and. (.not. required))) then
            par = default_value
        end if

    end if

end subroutine getxpar_

!
!> Read a parameter from a string by name, e.g.,
!> parname = 10
!
subroutine parsepar_(source, parname, par, default_value, required)

    character(len=*), dimension(:), intent(in) :: source
    character(len=*), intent(in) :: parname
    TT, intent(out) :: par
    TT, intent(in) :: default_value
    logical, intent(in), optional :: required

    integer :: i, j, eqindex
    logical :: parsepar_success
    character(len=256) :: arg

    parsepar_success = .false.

    do i = 1, size(source)

        arg = source(i)
        do j = 1, len(arg)
            if (arg(j:j) == achar(13) .or. arg(j:j) == achar(10)) then
                arg(j:j) = ' '
            end if
        end do
        eqindex = index(arg, '=')
        if (to_lower(tidy(arg(1:eqindex - 1))) == to_lower(tidy(parname))) then
            par = extract_TT(tidy(arg(eqindex + 1:)))
            parsepar_success = .true.
        end if

    end do

    if (present(required) .and. required .and. (.not. parsepar_success)) then
        write (error_unit, *) ' Error: Parameter ', trim(parname), ' is required but not set. '
        stop
    end if

    if ((.not. present(required) .or. (present(required) .and. (.not. required))) &
            .and. (.not. parsepar_success)) then
        par = default_value
    end if

end subroutine parsepar_

!
!> Read a parameter from a string by name, e.g.,
!> parname = 10, 20, 30, ...
!
subroutine parsenpar_(source, parname, par, default_value, separator, required)

    character(len=*), dimension(:), intent(in) :: source
    character(len=*), intent(in) :: parname
    TT, allocatable, dimension(:), intent(out) :: par
    TT, dimension(:), intent(in) :: default_value
    character(len=*), intent(in), optional :: separator
    logical, intent(in), optional :: required

    integer :: i, j, eqindex
    logical :: parsepar_success
    character(len=256) :: arg
    character(len=10) :: sepr

    parsepar_success = .false.
    if (present(separator)) then
        sepr = separator
    else
        sepr = ','
    end if

    do i = 1, size(source)

        arg = source(i)
        do j = 1, len(arg)
            if (arg(j:j) == achar(13) .or. arg(j:j) == achar(10)) then
                arg(j:j) = ' '
            end if
        end do
        eqindex = index(arg, '=')
        if (to_lower(tidy(arg(1:eqindex - 1))) == to_lower(tidy(parname))) then
            call extract_nTT(tidy(arg(eqindex + 1:)), trim(adjustl(sepr)), par)
            parsepar_success = .true.
        end if

    end do

    if (present(required) .and. required .and. (.not. parsepar_success)) then
        write (error_unit, *) ' Error: Parameter ', trim(parname), ' is required but not set. '
        stop
    end if

    if ((.not. present(required) .or. (present(required) .and. (.not. required))) &
            .and. (.not. parsepar_success)) then
        i = size(default_value, 1)
        if (allocated(par)) then
            deallocate (par)
        end if
        allocate (par(1:i))
        par = default_value
    end if

end subroutine parsenpar_

!
!> Read a parameter from a string using the format
!> parname=range1:value1,range2:value2,...
!> e.g., par=1~10:0.0,11~20:1.0,21~30:2.0
!> where the variable ranges from 1~10, the parameter value is 10, ...
!
subroutine parsexpar_(source, parname, par, default_value, var, required)

    character(len=*), dimension(:), intent(in) :: source
    character(len=*), intent(in) :: parname
    TT, allocatable, dimension(:), intent(out) :: par
    TT, dimension(:), intent(in) :: default_value
    real, intent(in) :: var
    logical, intent(in), optional :: required

    character(len=64), allocatable, dimension(:) :: tmpar
    real, allocatable, dimension(:) :: rbeg, rend
    TTT, allocatable, dimension(:) :: rpar
    integer :: i, j, l, nr, eqindex
    logical :: parsepar_success
    character(len=256) :: arg

    ! this initialization is necessary
    ! otherwise, segmentation fault
    allocate (tmpar(1:nstring_max))

    parsepar_success = .false.

    do i = 1, size(source)

        arg = source(i)
        do j = 1, len(arg)
            if (arg(j:j) == achar(13) .or. arg(j:j) == achar(10)) then
                arg(j:j) = ' '
            end if
        end do
        eqindex = index(arg, '=')
        if (to_lower(tidy(arg(1:eqindex - 1))) == to_lower(tidy(parname))) then
            call extract_nstring(tidy(arg(eqindex + 1:)), ',', tmpar)
            parsepar_success = .true.
        end if

    end do

    if (parsepar_success) then

        nr = size(tmpar)

        ! extracting range and values
        allocate (rbeg(1:nr))
        allocate (rend(1:nr))
        allocate (rpar(1:nr))

        do i = 1, nr
            call extract_xTT(tidy(tmpar(i)), rbeg(i), rend(i), rpar(i))
        end do

        ! find target value correponding to var
        if (var < rbeg(1)) then
            ! if var < min of given range
            par = rpar(1)

        else if (var > rend(nr)) then
            ! if var > max of given range
            par = rpar(nr)

        else
            ! if var between min and max of given range
            ! and in some given interval
            l = 0
            do i = 1, nr
                if (var >= rbeg(i) .and. var <= rend(i)) then
                    par = rpar(i)
                    l = l + 1
                end if
            end do

            ! if var is located somewhere between some non-defined interval [rend, rbeg_next]
            ! then linearly interpolate to obtain par value
            if (l == 0) then

#ifdef linear_interp
                ! Linear interpolation applies to float, double, complex
                do i = 1, nr - 1
                    if (var > rend(i) .and. var < rbeg(i + 1)) then
                        par = rpar(i) + (rpar(i + 1) - rpar(i))*(var - rend(i))/(rbeg(i + 1) - rend(i))
                        exit
                    end if
                end do
#endif
#ifdef nearest_interp
                ! Linear interpolation does not apply to int, logical, or string
                do i = 1, nr - 1
                    if (var > rend(i) .and. var < rend(i) + 0.5*(rbeg(i + 1) - rend(i))) then
                        par = rpar(i)
                        exit
                    else if (var >= rend(i) + 0.5*(rbeg(i + 1) - rend(i)) .and. var < rbeg(i + 1)) then
                        par = rpar(i + 1)
                        exit
                    end if
                end do
#endif
            end if

        end if

    else

        if (present(required) .and. required) then
            write (error_unit, *) ' Error: Parameter ', trim(parname), ' is required but not set. '
            stop
        end if

        if (.not. present(required) .or. (present(required) .and. (.not. required))) then
            par = default_value
        end if

    end if

end subroutine parsexpar_

#undef T
#undef TT
#undef TTT
#undef extract_TT
#undef extract_nTT
#undef extract_xTT
#undef linear_interp
#undef nearest_interp

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef readpar_
#undef readnpar_
#undef readxpar_
#undef getpar_
#undef getnpar_
#undef getxpar_
#undef parsepar_
#undef parsenpar_
#undef parsexpar_
