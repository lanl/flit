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

#define input_array_1d_      CONCAT(input_array_1d, T)
#define stdin_array_1d_      CONCAT(stdin_array_1d, T)

#define input_array_2d_      CONCAT(input_array_2d, T)
#define stdin_array_2d_      CONCAT(stdin_array_2d, T)

#define input_array_3d_      CONCAT(input_array_3d, T)
#define stdin_array_3d_      CONCAT(stdin_array_3d, T)

#define input_array_4d_      CONCAT(input_array_4d, T)
#define stdin_array_4d_      CONCAT(stdin_array_4d, T)

#define output_array_1d_      CONCAT(output_array_1d, T)
#define stdout_array_1d_      CONCAT(stdout_array_1d, T)

#define output_array_2d_      CONCAT(output_array_2d, T)
#define stdout_array_2d_      CONCAT(stdout_array_2d, T)

#define output_array_3d_      CONCAT(output_array_3d, T)
#define stdout_array_3d_      CONCAT(stdout_array_3d, T)

#define output_array_4d_      CONCAT(output_array_4d, T)
#define stdout_array_4d_      CONCAT(stdout_array_4d, T)

subroutine input_array_1d_(w, filename, pos, endian)

    TT, dimension(:), intent(inout) :: w
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: pos
    character(len=*), intent(in), optional :: endian

    integer :: funit, ppos
    character(len=24) :: endianness

    if (.not. present(pos)) then
        ppos = 1
    else
        ppos = pos
    end if

    if (present(endian)) then
        call assert(tidy(endian) == 'little' .or. tidy(endian) == 'big', &
            ' <intput_array_1d> Error: endian must be little or big. ')
        endianness = tidy(endian)//'_endian'
    else
        endianness = 'little_endian'
    end if

    call assert(file_exists(filename), &
        ' <input_array_1d> Error: '//tidy(filename)//' does not exist.')
    call assert(get_file_size(filename) >= size(w)*sizeof(w(1)), &
        ' <input_array_1d> Error: '//num2str(size(w)*sizeof(w(1))) &
        //' bytes expected while '//num2str(get_file_size(filename)) &
        //' bytes found. ')

    ! Read binary file
    open (newunit=funit, file=tidy(filename), &
        form='unformatted', action='read', access='stream', &
        status='old', convert=tidy(endianness))
    read (funit, pos=ppos) w
    close (funit)

end subroutine input_array_1d_

subroutine stdin_array_1d_(w, endian)

    TT, dimension(:), intent(inout) :: w
    character(len=*), intent(in), optional :: endian

    character(len=24) :: endianness

    if (present(endian)) then
        call assert(tidy(endian) == 'little' .or. tidy(endian) == 'big', &
            ' <stdin_array_1d> Error: endian must be little or big. ')
        endianness = tidy(endian)//'_endian'
    else
        endianness = 'little_endian'
    end if

    call assert(get_stdin_size() >= size(w)*sizeof(w(1)), &
        ' <stdin_array_1d> Error: '//num2str(size(w)*sizeof(w(1))) &
        //' bytes expected while '//num2str(get_stdin_size()) &
        //' bytes found. ')

    ! Read stream
    close (input_unit)
    open (unit=input_unit, form='unformatted', access='stream', convert=tidy(endianness))
    read (input_unit) w

end subroutine stdin_array_1d_

subroutine input_array_2d_(w, filename, transp, pos, endian)

    TT, dimension(:, :), intent(inout) :: w
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: transp
    integer, intent(in), optional :: pos
    character(len=*), intent(in), optional :: endian

    integer :: funit, ppos
    TT, allocatable, dimension(:, :) :: tw
    character(len=24) :: endianness

    if (.not. present(pos)) then
        ppos = 1
    else
        ppos = pos
    end if

    if (present(endian)) then
        call assert(tidy(endian) == 'little' .or. tidy(endian) == 'big', &
            ' <intput_array_2d> Error: endian must be little or big. ')
        endianness = tidy(endian)//'_endian'
    else
        endianness = 'little_endian'
    end if

    call assert(file_exists(filename), &
        ' <input_array_2d> Error: '//tidy(filename)//' does not exist.')
    call assert(get_file_size(filename) >= size(w)*sizeof(w(1, 1)), &
        ' <input_array_2d> Error: '//num2str(size(w)*sizeof(w(1, 1))) &
        //' bytes expected while '//num2str(get_file_size(filename)) &
        //' bytes found. ')

    ! if filename exists and file fize sufficient
    open (newunit=funit, file=tidy(filename), &
        form='unformatted', action='read', access='stream', &
        status='old', convert=tidy(endianness))
    if (present(transp)) then
        if (transp) then
            ! if data stored in different order
            allocate (tw(1:size(w, 2), 1:size(w, 1)))
            read (funit, pos=ppos) tw
            w = transpose(tw)
            deallocate (tw)
        end if
    else
        read (funit, pos=ppos) w
    end if
    close (funit)

end subroutine input_array_2d_

subroutine stdin_array_2d_(w, transp, endian)

    TT, dimension(:, :), intent(inout) :: w
    logical, intent(in), optional :: transp
    character(len=*), intent(in), optional :: endian

    TT, allocatable, dimension(:, :) :: tw
    character(len=24) :: endianness

    if (present(endian)) then
        call assert(tidy(endian) == 'little' .or. tidy(endian) == 'big', &
            ' <stdin_array_2d> Error: endian must be little or big. ')
        endianness = tidy(endian)//'_endian'
    else
        endianness = 'little_endian'
    end if

    call assert(get_stdin_size() >= size(w)*sizeof(w(1, 1)), &
        ' <stdin_array_2d> Error: '//num2str(size(w)*sizeof(w(1, 1))) &
        //' bytes expected while '//num2str(get_stdin_size()) &
        //' bytes found. ')

    close (input_unit)
    open (unit=input_unit, form='unformatted', access='stream', convert=tidy(endianness))
    if (present(transp)) then
        if (transp) then
            ! if data stored in different order
            allocate (tw(1:size(w, 2), 1:size(w, 1)))
            read (input_unit) tw
            w = transpose(tw)
            deallocate (tw)
        end if
    else
        read (input_unit) w
    end if

end subroutine stdin_array_2d_

subroutine input_array_3d_(w, filename, store, pos, endian)

    TT, dimension(:, :, :), intent(inout) :: w
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: store
    integer, intent(in), optional :: pos
    character(len=*), intent(in), optional :: endian

    integer :: funit, ppos
    TT, allocatable, dimension(:) :: tw
    integer, allocatable, dimension(:) :: sorder
    character(len=24) :: endianness

    if (.not. present(pos)) then
        ppos = 1
    else
        ppos = pos
    end if

    if (present(endian)) then
        call assert(tidy(endian) == 'little' .or. tidy(endian) == 'big', &
            ' <intput_array_3d> Error: endian must be little or big. ')
        endianness = tidy(endian)//'_endian'
    else
        endianness = 'little_endian'
    end if

    call assert(file_exists(filename), &
        ' <input_array_3d> Error: '//tidy(filename)//' does not exist.')
    call assert(get_file_size(filename) >= size(w)*sizeof(w(1, 1, 1)), &
        ' <input_array_3d> Error: '//num2str(size(w)*sizeof(w(1, 1, 1))) &
        //' bytes expected while '//num2str(get_file_size(filename)) &
        //' bytes found. ')

    ! if filename exists and file fize sufficient
    open (newunit=funit, file=tidy(filename), &
        form='unformatted', action='read', access='stream', &
        status='old', convert=tidy(endianness))
    if (present(store)) then
        ! if data stored in different order
        allocate (tw(1:size(w)))
        read (funit, pos=ppos) tw
        allocate (sorder(1:3))
        call get_integer_digits(store, 3, sorder)
        w = reshape(tw, shape=shape(w), order=sorder)
        deallocate (tw)
    else
        read (funit, pos=ppos) w
    end if
    close (funit)

end subroutine input_array_3d_

subroutine stdin_array_3d_(w, store, endian)

    TT, dimension(:, :, :), intent(inout) :: w
    integer, intent(in), optional :: store
    character(len=*), intent(in), optional :: endian

    TT, allocatable, dimension(:) :: tw
    integer, allocatable, dimension(:) :: sorder
    character(len=24) :: endianness

    if (present(endian)) then
        call assert(tidy(endian) == 'little' .or. tidy(endian) == 'big', &
            ' <stdin_array_3d> Error: endian must be little or big. ')
        endianness = tidy(endian)//'_endian'
    else
        endianness = 'little_endian'
    end if

    call assert(get_stdin_size() >= size(w)*sizeof(w(1, 1, 1)), &
        ' <stdin_array_3d> Error: '//num2str(size(w)*sizeof(w(1, 1, 1))) &
        //' bytes expected while '//num2str(get_stdin_size()) &
        //' bytes found. ')

    close (input_unit)
    open (unit=input_unit, form='unformatted', access='stream', convert=tidy(endianness))
    if (present(store)) then
        ! if data stored in different order
        allocate (tw(1:size(w)))
        read (input_unit, iostat=ioerr) tw
        allocate (sorder(1:3))
        call get_integer_digits(store, 3, sorder)
        w = reshape(tw, shape=shape(w), order=sorder)
        deallocate (tw)
    else
        read (input_unit) w
    end if

end subroutine stdin_array_3d_

subroutine input_array_4d_(w, filename, store, pos, endian)

    TT, dimension(:, :, :, :), intent(inout) :: w
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: store
    integer, intent(in), optional :: pos
    character(len=*), intent(in), optional :: endian

    integer :: funit, ppos
    TT, allocatable, dimension(:) :: tw
    integer, allocatable, dimension(:) :: sorder
    character(len=24) :: endianness

    if (.not. present(pos)) then
        ppos = 1
    else
        ppos = pos
    end if

    if (present(endian)) then
        call assert(tidy(endian) == 'little' .or. tidy(endian) == 'big', &
            ' <intput_array_4d> Error: endian must be little or big. ')
        endianness = tidy(endian)//'_endian'
    else
        endianness = 'little_endian'
    end if

    call assert(file_exists(filename), &
        ' <input_array_4d> Error: '//tidy(filename)//' does not exist.')
    call assert(get_file_size(filename) >= size(w)*sizeof(w(1, 1, 1, 1)), &
        ' <input_array_4d> Error: '//num2str(size(w)*sizeof(w(1, 1, 1, 1))) &
        //' bytes expected while '//num2str(get_file_size(filename)) &
        //' bytes found. ')

    ! if filename exists and file fize sufficient
    open (newunit=funit, file=tidy(filename), &
        form='unformatted', action='read', access='stream', &
        status='old', convert=tidy(endianness))
    if (present(store)) then
        ! if data stored in different order
        allocate (tw(1:size(w)))
        read (funit, pos=ppos) tw
        allocate (sorder(1:4))
        call get_integer_digits(store, 4, sorder)
        w = reshape(tw, shape=shape(w), order=sorder)
        deallocate (tw)
    else
        read (funit, pos=ppos) w
    end if
    close (funit)

end subroutine input_array_4d_

subroutine stdin_array_4d_(w, store, endian)

    TT, dimension(:, :, :, :), intent(inout) :: w
    integer, intent(in), optional :: store
    character(len=*), intent(in), optional :: endian

    TT, allocatable, dimension(:) :: tw
    integer, allocatable, dimension(:) :: sorder
    character(len=24) :: endianness

    if (present(endian)) then
        call assert(tidy(endian) == 'little' .or. tidy(endian) == 'big', &
            ' <stdin_array_4d> Error: endian must be little or big. ')
        endianness = tidy(endian)//'_endian'
    else
        endianness = 'little_endian'
    end if

    call assert(get_stdin_size() >= size(w)*sizeof(w(1, 1, 1, 1)), &
        ' <stdin_array_4d> Error: '//num2str(size(w)*sizeof(w(1, 1, 1, 1))) &
        //' bytes expected while '//num2str(get_stdin_size()) &
        //' bytes found. ')

    close (input_unit)
    open (unit=input_unit, form='unformatted', access='stream', convert=tidy(endianness))
    if (present(store)) then
        ! if data stored in different order
        allocate (tw(1:size(w)))
        read (input_unit, iostat=ioerr) tw
        allocate (sorder(1:4))
        call get_integer_digits(store, 4, sorder)
        w = reshape(tw, shape=shape(w), order=sorder)
        deallocate (tw)
    else
        read (input_unit) w
    end if

end subroutine stdin_array_4d_

subroutine output_array_1d_(w, filename, append, ascii, format)

    TT, dimension(:), intent(in) :: w
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: append, ascii
    character(len=*), optional :: format

    integer :: funit, i
    logical :: output_ascii
    character(len=12) :: output_format

    if (present(ascii)) then
        output_ascii = ascii
    else
        output_ascii = .false.
    end if

    if (present(format)) then
        output_format = tidy(format)
    else
        output_format = 'es'
    end if

    if (output_ascii) then
        ! Output ASCII

        if (present(append)) then
            if (append) then
                open (newunit=funit, file=tidy(filename), &
                    form='formatted', action='write', access='stream', &
                    status='old', position='append')
            end if
        else
            open (newunit=funit, file=tidy(filename), access='stream', &
                form='formatted', action='write', &
                status='replace')
        end if
        do i = 1, size(w)
            write (funit, '('//tidy(output_format)//')') w(i)
        end do

    else
        ! Output binary

        if (present(append)) then
            if (append) then
                open (newunit=funit, file=tidy(filename), &
                    form='unformatted', action='write', access='stream', &
                    status='old', position='append')
            end if
        else
            open (newunit=funit, file=tidy(filename), &
                form='unformatted', action='write', access='stream', &
                status='replace')
        end if
        write (funit) w

    end if

    close (funit)

end subroutine output_array_1d_

subroutine stdout_array_1d_(w, append, ascii, format)

    TT, dimension(:), intent(in) :: w
    logical, intent(in), optional :: append, ascii
    character(len=*), optional :: format

    logical :: output_ascii
    integer :: i
    character(len=12) :: output_format

    if (present(ascii)) then
        output_ascii = ascii
    else
        output_ascii = .false.
    end if

    if (present(format)) then
        output_format = tidy(format)
    else
        output_format = 'es'
    end if

    close (output_unit)

    if (output_ascii) then
        ! Output ASCII

        if (present(append)) then
            if (append) then
                open (unit=output_unit, &
                    form='formatted', action='write', access='stream', &
                    status='old', position='append')
            end if
        else
            open (unit=output_unit, &
                form='formatted', action='write', access='stream', &
                status='replace')
        end if
        do i = 1, size(w)
            write (output_unit, '('//tidy(output_format)//')') w(i)
        end do

    else

        if (present(append)) then
            if (append) then
                open (unit=output_unit, &
                    form='unformatted', action='write', access='stream', &
                    status='old', position='append')
            end if
        else
            open (unit=output_unit, &
                form='unformatted', action='write', access='stream', &
                status='replace')
        end if
        write (output_unit) w

    end if

end subroutine stdout_array_1d_

subroutine output_array_2d_(w, filename, transp, append, ascii, format)

    TT, dimension(:, :), intent(in) :: w
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: transp, append, ascii
    character(len=*), optional :: format

    integer :: funit, i
    logical :: output_ascii
    character(len=12) :: output_format

    if (present(ascii)) then
        output_ascii = ascii
    else
        output_ascii = .false.
    end if

    if (present(format)) then
        output_format = tidy(format)
    else
        output_format = 'es'
    end if

    if (output_ascii) then
        ! Output ASCII

        if (present(append)) then
            if (append) then
                open (newunit=funit, file=tidy(filename), &
                    form='formatted', action='write', access='stream', &
                    status='old', position='append')
            end if
        else
            open (newunit=funit, file=tidy(filename), access='stream', &
                form='formatted', action='write', &
                status='replace')
        end if

        if (present(transp)) then
            if (transp) then
                do i = 1, size(w, 2)
                    write (funit, '('//num2str(size(w, 1))//tidy(output_format)//')') w(:, i)
                end do
            end if
        else
            do i = 1, size(w, 1)
                write (funit, '('//num2str(size(w, 1))//tidy(output_format)//')') w(i, :)
            end do
        end if

    else
        ! Output binary

        if (present(append)) then
            if (append) then
                open (newunit=funit, file=tidy(filename), &
                    form='unformatted', action='write', access='stream', &
                    status='old', position='append')
            end if
        else
            open (newunit=funit, file=tidy(filename), &
                form='unformatted', action='write', access='stream', &
                status='replace')
        end if

        if (present(transp)) then
            if (transp) then
                write (funit) transpose(w)
            end if
        else
            write (funit) w
        end if

    end if

    close (funit)

end subroutine output_array_2d_

subroutine stdout_array_2d_(w, transp, append, ascii, format)

    TT, dimension(:, :), intent(in) :: w
    logical, intent(in), optional :: transp, append, ascii
    character(len=*), optional :: format

    logical :: output_ascii
    integer :: i
    character(len=12) :: output_format

    if (present(ascii)) then
        output_ascii = ascii
    else
        output_ascii = .false.
    end if

    if (present(format)) then
        output_format = tidy(format)
    else
        output_format = 'es'
    end if

    close (output_unit)
    if (output_ascii) then
        ! Output ASCII

        if (present(append)) then
            if (append) then
                open (unit=output_unit, &
                    form='formatted', action='write', access='stream', &
                    status='old', position='append')
            end if
        else
            open (unit=output_unit, &
                form='formatted', action='write', access='stream', &
                status='replace')
        end if
        if (present(transp)) then
            if (transp) then
                do i = 1, size(w, 2)
                    write (output_unit, '('//num2str(size(w, 1))//'es)') w(:, i)
                end do
            end if
        else
            do i = 1, size(w, 1)
                write (output_unit, '('//num2str(size(w, 1))//'es)') w(i, :)
            end do
        end if

    else
        ! Output binary

        if (present(append)) then
            if (append) then
                open (unit=output_unit, &
                    form='unformatted', action='write', access='stream', &
                    status='old', position='append')
            end if
        else
            open (unit=output_unit, &
                form='unformatted', action='write', access='stream', &
                status='replace')
        end if
        if (present(transp)) then
            if (transp) then
                write (output_unit) transpose(w)
            end if
        else
            write (output_unit) w
        end if

    end if

end subroutine stdout_array_2d_

subroutine output_array_3d_(w, filename, store, append)

    TT, dimension(:, :, :), intent(in) :: w
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: store
    logical, intent(in), optional :: append

    integer :: funit

    if (present(append)) then
        if (append) then
            open (newunit=funit, file=tidy(filename), &
                form='unformatted', action='write', access='stream', &
                status='old', position='append')
        end if
    else
        open (newunit=funit, file=tidy(filename), &
            form='unformatted', action='write', access='stream', &
            status='replace')
    end if
    if (present(store)) then
        write (funit) permute(w, store)
    else
        write (funit) w
    end if
    close (funit)

end subroutine output_array_3d_

subroutine stdout_array_3d_(w, store, append)

    TT, dimension(:, :, :), intent(in) :: w
    integer, intent(in), optional :: store
    logical, intent(in), optional :: append

    close (output_unit)
    if (present(append)) then
        if (append) then
            open (unit=output_unit, &
                form='unformatted', action='write', access='stream', &
                status='old', position='append')
        end if
    else
        open (unit=output_unit, &
            form='unformatted', action='write', access='stream', &
            status='replace')
    end if
    if (present(store)) then
        write (output_unit) permute(w, store)
    else
        write (output_unit) w
    end if

end subroutine stdout_array_3d_

subroutine output_array_4d_(w, filename, store, append)

    TT, dimension(:, :, :, :), intent(in) :: w
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: store
    logical, intent(in), optional :: append

    integer :: funit

    if (present(append)) then
        if (append) then
            open (newunit=funit, file=tidy(filename), &
                form='unformatted', action='write', access='stream', &
                status='old', position='append')
        end if
    else
        open (newunit=funit, file=tidy(filename), &
            form='unformatted', action='write', access='stream', &
            status='replace')
    end if
    if (present(store)) then
        write (funit) permute(w, store)
    else
        write (funit) w
    end if
    close (funit)

end subroutine output_array_4d_

subroutine stdout_array_4d_(w, store, append)

    TT, dimension(:, :, :, :), intent(in) :: w
    integer, intent(in), optional :: store
    logical, intent(in), optional :: append

    close (output_unit)
    if (present(append)) then
        if (append) then
            open (unit=output_unit, &
                form='unformatted', action='write', access='stream', &
                status='old', position='append')
        end if
    else
        open (unit=output_unit, &
            form='unformatted', action='write', access='stream', &
            status='replace')
    end if
    if (present(store)) then
        write (output_unit) permute(w, store)
    else
        write (output_unit) w
    end if

end subroutine stdout_array_4d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef input_array_1d_
#undef stdin_array_1d_

#undef input_array_2d_
#undef stdin_array_2d_

#undef input_array_3d_
#undef stdin_array_3d_

#undef input_array_4d_
#undef stdin_array_4d_

#undef output_array_1d_
#undef stdout_array_1d_

#undef output_array_2d_
#undef stdout_array_2d_

#undef output_array_3d_
#undef stdout_array_3d_

#undef output_array_4d_
#undef stdout_array_4d_
