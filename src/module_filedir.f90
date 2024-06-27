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


module libflit_filedir

    use libflit_string
    use libflit_random
    use libflit_error
    use iso_fortran_env
    use iso_c_binding
    use filesystem

contains

    !
    !> Get current directory
    !
    function get_current_directory() result(dir_name)

        character(len=:), allocatable :: dir_name

        dir_name = tidy(get_cwd())

    end function get_current_directory

    !
    !> Get the real/absolute path of a file
    !
    function get_real_path(filename) result(path)

        character(len=*), intent(in) :: filename
        character(len=:), allocatable :: path

        path = canonical(tidy(filename))

    end function get_real_path

    !
    !> Get the directory of a file
    !
    function get_file_directory(filename) result(path)

        character(len=*), intent(in) :: filename
        character(len=:), allocatable :: path

        character(len=:), allocatable :: fname

		fname = get_basename(filename)
        path = canonical(tidy(filename))
        path = path(1:len(path) - len(fname))

    end function get_file_directory

    !
    !> Make a directory
    !
    subroutine make_directory(dir)

        character(len=*), intent(in) :: dir

        if (.not. directory_exists(dir)) then
            call execute_command_line("mkdir -p "//tidy(dir))
        end if

    end subroutine make_directory

    !
    !> Delete a directory
    !
    subroutine delete_directory(dir)

        character(len=*), intent(in) :: dir

        call execute_command_line("rm -rf "//tidy(dir))

    end subroutine delete_directory

    !
    !> Rename a directory
    !
    subroutine move_directory(source_dir, target_dir)

        character(len=*), intent(in) :: source_dir, target_dir

        if (directory_exists(source_dir)) then
            call rename(tidy(source_dir), tidy(target_dir))
        end if

    end subroutine move_directory

    !
    !> Check if dir is a directory
    !
    function is_directory(dir) result(isadir)

        character(len=*), intent(in) :: dir
        logical :: isadir

        isadir = is_dir(get_real_path(tidy(dir)))

    end function is_directory

    !
    !> Check file existence
    !
    function file_exists(filename) result(fe)

        character(len=*), intent(in) :: filename
        logical :: fe

        inquire (file=tidy(filename), exist=fe)

    end function file_exists

    !
    !> Check file existence
    !
    function directory_exists(dirname) result(de)

        character(len=*), intent(in) :: dirname
        logical :: de

        inquire (directory=tidy(dirname), exist=de)

    end function directory_exists

    !
    !> Get file size
    !
    function get_file_size(filename) result(fs)

        character(len=*), intent(in) :: filename
        integer(8) :: fs

        if (file_exists(filename)) then
            inquire (file=tidy(filename), size=fs)
        else
            write (error_unit, *) ' <get_file_size> Warning: File '//tidy(filename)// &
                ' not found; Return size = 0. '
            fs = 0
        end if

    end function get_file_size

    !
    !> Get standard input size
    !
    function get_stdin_size() result(sse)

        integer(8) :: sse

        inquire (unit=input_unit, size=sse)

    end function get_stdin_size

    !
    !> Delete a file
    !
    subroutine delete_file(filename)

        character(len=*), intent(in) :: filename

        integer :: fileunit

        if (file_exists(filename)) then
            open (newunit=fileunit, file=tidy(filename))
            close (fileunit, status='delete')
        else
            write (error_unit, *) ' <delete_file> Warning: File '//tidy(filename)//' not found. '
        end if

    end subroutine delete_file

    !
    !> Copy a file
    !
    subroutine copy_file(source_file, target_file)

        character(len=*), intent(in) :: source_file, target_file

        if (file_exists(source_file)) then
            ! call execute_command_line("cp "//tidy(source_file)//' '//tidy(target_file))
            !call copyfile(tidy(source_file), tidy(target_file))
            call execute_command_line('rsync -azq '//tidy(source_file)//' '//tidy(target_file))
        else
            write (error_unit, *) ' <copy_file> Warning: Source file '//tidy(source_file)//' not found. '
        end if

    end subroutine copy_file

    !
    !> Move a file
    !
    subroutine move_file(source_file, target_file)

        character(len=*), intent(in) :: source_file, target_file

        if (file_exists(source_file)) then
            !call execute_command_line("mv "//tidy(source_file)//' '//tidy(target_file))
            call rename(tidy(source_file), tidy(target_file))
        else
            write (error_unit, *) ' <move_file> Warning: Source file '//tidy(source_file)//' not found. '
        end if

    end subroutine move_file

    !
    !> List all files within a directory
    !
    subroutine list_files(dir, filelist, pattern)

        character(len=*), intent(in) :: dir
        character(len=*), allocatable, dimension(:), intent(inout) :: filelist
        character(len=*), intent(in), optional :: pattern

        character(len=128) :: tmpfile
        integer :: fileunit
        integer :: i, nf

        ! first list files with system command
        tmpfile = 'tmplist_'//random_string()
        if (present(pattern)) then
            call execute_command_line("ls -1 -d "//tidy(dir)//'/'//tidy(pattern) &
                //' 2>/dev/null 1>'//tidy(tmpfile))
        else
            call execute_command_line("ls -1 -d "//tidy(dir)//'/*  2>/dev/null 1>' &
                //tidy(tmpfile))
        end if
        nf = count_lines(tmpfile)

        ! second read in files
        if (allocated(filelist)) then
            deallocate (filelist)
        end if
        allocate (filelist(1:nf))
        open (newunit=fileunit, file=tidy(tmpfile), status='old', action='read')
        do i = 1, nf
            read (fileunit, '(a)') filelist(i)
        end do
        close (fileunit, status='delete')

    end subroutine list_files

    !
    !> Count the number of all files or that matche a given pattern within a directory
    !
    function count_files(dir, pattern) result(nf)

        character(len=*), intent(in) :: dir
        character(len=*), intent(in), optional :: pattern

        integer :: nf

        character(len=128) :: tmpfile

        tmpfile = 'tmplist_'//random_string()

        if (present(pattern)) then
            call execute_command_line("ls -1 -d "//tidy(dir)//'/'//tidy(pattern) &
                //' 2>/dev/null 1>'//tidy(tmpfile))
        else
            call execute_command_line("ls -1 -d "//tidy(dir)//'/*  2>/dev/null 1>' &
                //tidy(tmpfile))
        end if

        nf = count_lines(tmpfile)

        call delete_file(tmpfile)

    end function count_files

    !
    !> Get extension from filename
    !
    function get_extension(filename) result(extension)

        character(len=*) :: filename
        character(len=:), allocatable :: extension

        integer :: ipos
        character(len=1024) :: tmpname

        ipos = index(filename, '.', back=.true.)
        tmpname = filename(ipos + 1:)

        allocate (character(len=len(tidy(tmpname))) :: extension)
        extension = tidy(tmpname)

    end function get_extension

    !
    !> Strip path from filename
    !
    function strip_extension(filename) result(extension)

        character(len=*) :: filename
        character(len=:), allocatable :: extension

        integer :: ipos
        character(len=1024) :: tmpname

        ipos = index(filename, '.', back=.true.)
        tmpname = filename(1:ipos - 1)

        allocate (character(len=len(tidy(tmpname))) :: extension)
        extension = tidy(tmpname)

    end function strip_extension

    !
    !> Strip path from filename
    !
    function get_basename(filename) result(basename)

        character(len=*) :: filename
        character(len=:), allocatable :: basename

        integer :: ipos
        character(len=1024) :: tmpname

        ipos = index(filename, '/', back=.true.)
        tmpname = filename(ipos + 1:)

        allocate (character(len=len(tidy(tmpname))) :: basename)
        basename = tidy(tmpname)

    end function get_basename

    !
    !> Count lines of an ascii file
    !
    function count_lines(infile) result(nl)

        character(len=*), intent(in) :: infile
        integer :: fileunit
        integer :: nl, ioerr

        open (newunit=fileunit, file=tidy(infile), status='old', action='read')
        nl = 0
        do
            read (fileunit, *, iostat=ioerr)
            if (ioerr /= 0) then
                exit
            end if
            nl = nl + 1
        end do
        close (fileunit)

    end function count_lines

    function count_nonempty_lines(infile) result(nl)

        character(len=*), intent(in) :: infile
        integer :: fileunit
        integer :: nl, ioerr
        character(len=1024) :: str

        open (newunit=fileunit, file=tidy(infile), status='old', action='read')
        nl = 0
        do
            read (fileunit, '(a)', iostat=ioerr) str
            if (ioerr /= 0) then
                exit
            end if
            if (str /= '') then
                nl = nl + 1
            end if
        end do
        close (fileunit)

    end function count_nonempty_lines

end module libflit_filedir
