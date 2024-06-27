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


!
! ................ Fortran reshape usage ................
!
!       reshape(w, shape=[], order=[])
!
! Note that the order indicates "where do the old dims go in new array"
! e.g., reshape the array w(35,30,40) to w(40,35,30):
!
!       reshape(w, shape=[40,35,30], order=[?])
!
! then order should be order=[2, 3, 1], because:
! 35 becomes dimension 2 in the new array
! 30 becomes dimension 3 in the new array
! 40 becomes dimension 1 in the new array
!
! The following operation is permitted:
!
!       w = reshape(w, shape=..., order=...)
!
! as long as w is allocatable array
!
module libflit_io

    use libflit_string
    use libflit_utility
    use libflit_date_time
    use libflit_filedir
    use libflit_array
    use libflit_error
    use iso_fortran_env
    use libflit_array_operation

    implicit none

    private

    interface input_array
        module procedure :: input_array_1d_int
        module procedure :: input_array_2d_int
        module procedure :: input_array_3d_int
        module procedure :: input_array_4d_int
        module procedure :: input_array_1d_float
        module procedure :: input_array_2d_float
        module procedure :: input_array_3d_float
        module procedure :: input_array_4d_float
        module procedure :: input_array_1d_double
        module procedure :: input_array_2d_double
        module procedure :: input_array_3d_double
        module procedure :: input_array_4d_double
        module procedure :: input_array_1d_complex
        module procedure :: input_array_2d_complex
        module procedure :: input_array_3d_complex
        module procedure :: input_array_4d_complex
        module procedure :: input_array_1d_dcomplex
        module procedure :: input_array_2d_dcomplex
        module procedure :: input_array_3d_dcomplex
        module procedure :: input_array_4d_dcomplex
    end interface input_array

    interface stdin_array
        module procedure :: stdin_array_1d_int
        module procedure :: stdin_array_2d_int
        module procedure :: stdin_array_3d_int
        module procedure :: stdin_array_4d_int
        module procedure :: stdin_array_1d_float
        module procedure :: stdin_array_2d_float
        module procedure :: stdin_array_3d_float
        module procedure :: stdin_array_4d_float
        module procedure :: stdin_array_1d_double
        module procedure :: stdin_array_2d_double
        module procedure :: stdin_array_3d_double
        module procedure :: stdin_array_4d_double
        module procedure :: stdin_array_1d_complex
        module procedure :: stdin_array_2d_complex
        module procedure :: stdin_array_3d_complex
        module procedure :: stdin_array_4d_complex
        module procedure :: stdin_array_1d_dcomplex
        module procedure :: stdin_array_2d_dcomplex
        module procedure :: stdin_array_3d_dcomplex
        module procedure :: stdin_array_4d_dcomplex
    end interface stdin_array

    interface output_array
        module procedure :: output_array_1d_int
        module procedure :: output_array_2d_int
        module procedure :: output_array_3d_int
        module procedure :: output_array_4d_int
        module procedure :: output_array_1d_float
        module procedure :: output_array_2d_float
        module procedure :: output_array_3d_float
        module procedure :: output_array_4d_float
        module procedure :: output_array_1d_double
        module procedure :: output_array_2d_double
        module procedure :: output_array_3d_double
        module procedure :: output_array_4d_double
        module procedure :: output_array_1d_complex
        module procedure :: output_array_2d_complex
        module procedure :: output_array_3d_complex
        module procedure :: output_array_4d_complex
        module procedure :: output_array_1d_dcomplex
        module procedure :: output_array_2d_dcomplex
        module procedure :: output_array_3d_dcomplex
        module procedure :: output_array_4d_dcomplex
    end interface output_array

    interface stdout_array
        module procedure :: stdout_array_1d_int
        module procedure :: stdout_array_2d_int
        module procedure :: stdout_array_3d_int
        module procedure :: stdout_array_4d_int
        module procedure :: stdout_array_1d_float
        module procedure :: stdout_array_2d_float
        module procedure :: stdout_array_3d_float
        module procedure :: stdout_array_4d_float
        module procedure :: stdout_array_1d_double
        module procedure :: stdout_array_2d_double
        module procedure :: stdout_array_3d_double
        module procedure :: stdout_array_4d_double
        module procedure :: stdout_array_1d_complex
        module procedure :: stdout_array_2d_complex
        module procedure :: stdout_array_3d_complex
        module procedure :: stdout_array_4d_complex
        module procedure :: stdout_array_1d_dcomplex
        module procedure :: stdout_array_2d_dcomplex
        module procedure :: stdout_array_3d_dcomplex
        module procedure :: stdout_array_4d_dcomplex
    end interface stdout_array

    integer :: ioerr
    logical :: file_exist
    integer :: file_size

    interface load
        module procedure :: load_array_1d_float
        module procedure :: load_array_2d_float
        module procedure :: load_array_3d_float
    end interface load

    public :: input_array
    public :: stdin_array
    public :: output_array
    public :: stdout_array
    public :: load

contains

#define T int
#define TT integer
#include "template_io.f90"

#define T float
#define TT real
#include "template_io.f90"

#define T double
#define TT double precision
#include "template_io.f90"

#define T complex
#define TT complex
#include "template_io.f90"

#define T dcomplex
#define TT double complex
#include "template_io.f90"

    function load_array_1d_float(filename, n1, ascii, pos, endian) result(w)

        integer, intent(in) :: n1
        character(len=*), intent(in) :: filename
        logical, intent(in), optional :: ascii
        integer, intent(in), optional :: pos
        character(len=*), intent(in), optional :: endian

        real, allocatable, dimension(:) :: w

        integer :: funit, fpos, i
        character(len=24) :: endianness
        logical :: file_ascii
        integer :: target_size

        if (present(ascii)) then
            file_ascii = ascii
        else
            file_ascii = .false.
        end if

        if (.not. present(pos)) then
            fpos = 1
        else
            fpos = pos
        end if

        if (.not. file_ascii) then
            if (.not. present(endian)) then
                endianness = 'little_endian'
            else
                endianness = tidy(endian)//'_endian'
                if (tidy(endianness) /= 'little_endian' &
                        .and. tidy(endianness) /= 'big_endian') then
                    endianness = 'little_endian'
                end if
            end if
        end if

        target_size = n1*4

        if (file_ascii) then
            ! ASCII file

            if (.not. file_exists(filename)) then
                ! if filename does not exist
                call warn(' <load_array_1d_float> Error: '//tidy(filename)//' does not exist.')
                stop
            else
                w = zeros(n1)
                open (newunit=funit, file=tidy(filename), action='read', status='old')
                do i = 1, n1
                    if (i >= fpos) then
                        read (funit, *) w(i)
                    end if
                end do
                close(funit)
            end if

        else
            ! If binary file

            if (.not. file_exists(filename)) then
                ! if filename does not exist
                call warn(' <load_array_1d_float> Error: '//tidy(filename)//' does not exist.')
                stop
            else
                if (get_file_size(filename) < target_size) then
                    ! if file size insufficient
                    call warn(' <load_array_1d_float> Error: '//num2str(get_file_size(filename)) &
                        //' bytes found while '//num2str(target_size) &
                        //' bytes required.')
                    stop
                else
                    ! if filename exists and file fize sufficient
                    open (newunit=funit, file=tidy(filename), &
                        form='unformatted', action='read', access='stream', &
                        status='old', convert=tidy(endianness))
                    w = zeros(n1)
                    read (funit, pos=fpos) w
                    close (funit)
                end if
            end if

        end if

    end function load_array_1d_float

    function load_array_2d_float(filename, n1, n2, transp, ascii, pos, endian) result(w)

        integer, intent(in) :: n1, n2
        character(len=*), intent(in) :: filename
        logical, intent(in), optional :: transp, ascii
        integer, intent(in), optional :: pos
        character(len=*), intent(in), optional :: endian

        real, allocatable, dimension(:, :) :: w

        integer :: funit, fpos, i
        character(len=24) :: endianness
        logical :: file_ascii
        integer :: target_size

        if (present(ascii)) then
            file_ascii = ascii
        else
            file_ascii = .false.
        end if

        if (.not. present(pos)) then
            fpos = 1
        else
            fpos = pos
        end if

        if (.not. file_ascii) then
            if (.not. present(endian)) then
                endianness = 'little_endian'
            else
                endianness = tidy(endian)//'_endian'
                if (tidy(endianness) /= 'little_endian' &
                        .and. tidy(endianness) /= 'big_endian') then
                    endianness = 'little_endian'
                end if
            end if
        end if

        target_size = n1*n2*4

        if (file_ascii) then
            ! ASCII file

            if (.not. file_exists(filename)) then
                ! if filename does not exist
                call warn(' <load_array_2d_float> Error: '//tidy(filename)//' does not exist.')
                stop
            else
                w = zeros(n1, n2)
                open (newunit=funit, file=tidy(filename), action='read', status='old')
                do i = 1, n1
                    if (i >= fpos) then
                        read (funit, *) w(i, :)
                    end if
                end do
                close(funit)
            end if

        else
            ! If binary file

            if (.not. file_exists(filename)) then
                ! if filename does not exist
                call warn(' <load_array_2d_float> Error: '//tidy(filename)//' does not exist.')
                stop
            else
                if (get_file_size(filename) < target_size) then
                    ! if file size insufficient
                    call warn(' <load_array_2d_float> Error: '//num2str(get_file_size(filename)) &
                        //' bytes found while '//num2str(target_size) &
                        //' bytes required.')
                    stop
                else
                    ! if filename exists and file fize sufficient
                    open (newunit=funit, file=tidy(filename), &
                        form='unformatted', action='read', access='stream', &
                        status='old', convert=tidy(endianness))
                    if (present(transp)) then
                        if (transp) then
                            ! if data stored in different order
                            w = zeros(n2, n1)
                            read (funit, pos=fpos) w
                            w = transpose(w)
                        end if
                    else
                        w = zeros(n1, n2)
                        read (funit, pos=fpos) w
                    end if
                    close (funit)
                end if
            end if

        end if

    end function load_array_2d_float

    function load_array_3d_float(filename, n1, n2, n3, store, ascii, pos, endian) result(w)

        integer, intent(in) :: n1, n2, n3
        character(len=*), intent(in) :: filename
        integer, intent(in), optional :: store
        logical, intent(in), optional :: ascii
        integer, intent(in), optional :: pos
        character(len=*), intent(in), optional :: endian

        real, allocatable, dimension(:, :, :) :: w
        real, allocatable, dimension(:) :: ww
        integer, allocatable, dimension(:) :: sorder

        integer :: funit, fpos, i, j, l
        character(len=24) :: endianness
        logical :: file_ascii
        integer :: target_size

        if (present(ascii)) then
            file_ascii = ascii
        else
            file_ascii = .false.
        end if

        if (.not. present(pos)) then
            fpos = 1
        else
            fpos = pos
        end if

        if (.not. file_ascii) then
            if (.not. present(endian)) then
                endianness = 'little_endian'
            else
                endianness = tidy(endian)//'_endian'
                if (tidy(endianness) /= 'little_endian' &
                        .and. tidy(endianness) /= 'big_endian') then
                    endianness = 'little_endian'
                end if
            end if
        end if

        target_size = n1*n2*n3*4

        if (file_ascii) then
            ! ASCII file

            if (.not. file_exists(filename)) then
                ! if filename does not exist
                call warn(' <load_array_3d_float> Error: '//tidy(filename)//' does not exist.')
                stop
            else
                w = zeros(n1, n2, n3)
                open (newunit=funit, file=tidy(filename), action='read', status='old')
                l = 1
                do i = 1, n1
                    do j = 1, n2
                        if (l >= fpos) then
                            read (funit, *) w(i, j, :)
                            l = l + 1
                        end if
                    end do
                end do
                close(funit)
            end if

        else
            ! If binary file

            if (.not. file_exists(filename)) then
                ! if filename does not exist
                call warn(' <load_array_3d_float> Error: '//tidy(filename)//' does not exist.')
                stop
            else
                if (get_file_size(filename) < target_size) then
                    ! if file size insufficient
                    call warn(' <load_array_3d_float> Error: '//num2str(get_file_size(filename)) &
                        //' bytes found while '//num2str(target_size) &
                        //' bytes required.')
                    stop
                else
                    ! if filename exists and file fize sufficient
                    open (newunit=funit, file=tidy(filename), &
                        form='unformatted', action='read', access='stream', &
                        status='old', convert=tidy(endianness))
                    if (present(store)) then
                        ! if data stored in different order
                        allocate (ww(1:n1*n2*n3))
                        read (funit, pos=fpos) ww
                        allocate (sorder(1:3))
                        call get_integer_digits(store, 3, sorder)
                        w = reshape(ww, shape=[n1, n2, n3], order=sorder)
                    else
                        w = zeros(n1, n2, n3)
                        read (funit, pos=fpos) w
                    end if
                    close (funit)
                end if
            end if

        end if

    end function load_array_3d_float

end module libflit_io
