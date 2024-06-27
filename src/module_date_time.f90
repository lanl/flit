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


module libflit_date_time

    implicit none

contains

    !
    !> Convert a time duration in seconds to
    !>        ? day ? hours ? minutes ? seconds
    !
    function time_duration(duration_time) result(duration)

        real, intent(in) :: duration_time

        integer :: nmin, nhrs, nday
        real :: nsec
        character(len=:), allocatable :: duration
        character(len=64) :: string_duration

        if (duration_time < 60.0d0) then
            ! shorter than a minute
            write (string_duration, '(f5.2,a)') duration_time, ' sec.'
        else if ((duration_time >= 60.0d0) .and. (duration_time < 60*60.0d0)) then
            ! shorter than a hour
            nmin = floor(duration_time/60.0d0)
            nsec = duration_time - nmin*60.0d0
            write (string_duration, '(i2,a,f5.2,a)') nmin, ' min ', nsec, ' sec.'
        else if ((duration_time >= 60*60.0d0) .and. (duration_time < 24*60*60.0d0)) then
            ! shorter than a day
            nhrs = floor(duration_time/(60*60.0d0))
            nmin = floor((duration_time - nhrs*60*60.0d0)/60.0d0)
            nsec = duration_time - nhrs*60*60.0d0 - nmin*60.0d0
            write (string_duration, '(i2,a,i2,a,f5.2,a)') nhrs, ' hrs ', nmin, ' min ', nsec, ' sec.'
        else
            ! longer than a day
            nday = floor(duration_time/(24*60*60.0d0))
            nhrs = floor((duration_time - nday*24*60*60.0d0)/(60*60.0d0))
            nmin = floor((duration_time - nday*24*60*60.0d0 - nhrs*60*60.0d0)/60.0d0)
            nsec = duration_time - nday*24*60*60.0d0 - nhrs*60*60.0d0 - nmin*60.0d0
            write (string_duration, '(i4,a,i2,a,i2,a,f5.2,a)') nday, ' day ', nhrs, ' hrs ', nmin, ' min ', nsec, ' sec.'
        end if

        ! print
        allocate (character(len=len_trim(string_duration)) :: duration)
        duration = trim(adjustl(string_duration))

    end function time_duration

    !
    !> Return a compact date-time string
    !
    function date_time_compact() result(w)

        integer :: date_time(8)
        character(len=12) :: real_clock(3)
        character(len=22) :: w

        call date_and_time(real_clock(1), real_clock(2), real_clock(3), date_time)
        write (w, '(a,i4.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a)') &
            ' [', date_time(1), '-', date_time(2), '-', date_time(3), &
            ' ', date_time(5), ':', date_time(6), ':', date_time(7), ']'

    end function date_time_compact

    !
    !> Return a date-time stamp
    !
    function date_time_stamp() result(w)

        integer :: date_time(8)
        character(len=12) :: real_clock(3)
        character(len=22) :: w

        call date_and_time(real_clock(1), real_clock(2), real_clock(3), date_time)
        write (w, '(i4.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') &
            date_time(1), '-', date_time(2), '-', date_time(3), &
            '_', date_time(5), ':', date_time(6), ':', date_time(7)

    end function date_time_stamp

    !
    !> Return a date-time string with only numbers and dashes
    !
    function date_time_string() result(w)

        integer :: date_time(8)
        character(len=12) :: real_clock(3)
        character(len=22) :: w

        call date_and_time(real_clock(1), real_clock(2), real_clock(3), date_time)
        write (w, '(i4.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') &
            date_time(1), '-', date_time(2), '-', date_time(3), &
            '-', date_time(5), '-', date_time(6), '-', date_time(7)

    end function date_time_string

    !
    !> Display date and time information
    !
    subroutine print_date_time()

        integer :: date_time(8)
        character(len=12) :: real_clock(3)

        call date_and_time(real_clock(1), real_clock(2), real_clock(3), date_time)

        write (*, '(a,i4.2,a,i2.2,a,i2.2)') ' @Date: ', date_time(1), '-', date_time(2), '-', date_time(3)
        write (*, '(a,a)') ' @Zone: ', real_clock(3)
        write (*, '(a,i2.2,a,i2.2,a,i2.2)') ' @Time: ', date_time(5), ':', date_time(6), ':', date_time(7)

    end subroutine print_date_time

end module libflit_date_time
