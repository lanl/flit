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


module libflit_error

    use iso_fortran_env

    implicit none

    private
    public :: assert
    public :: remind

contains

    !
    !> Ensure that a condition is satisfied,
    !>        otherwise, give error and stop
    !
    subroutine assert(condition, error_message)

        logical, intent(in) :: condition
        character(len=*), intent(in), optional :: error_message

        if (.not. condition) then
            if (present(error_message)) then
                write (error_unit, *) ' '//trim(adjustl(error_message))
            else
                write (error_unit, *) ' Error: Condition is violated. Exit.'
            end if
            stop
        end if

    end subroutine assert

    subroutine remind(condition, warn_message)

        logical, intent(in) :: condition
        character(len=*), intent(in), optional :: warn_message

        if (.not. condition) then
            if (present(warn_message)) then
                write (error_unit, *) ' '//trim(adjustl(warn_message))
            else
                write (error_unit, *) ' Warning: Condition is violated. Continue. '
            end if
        end if

    end subroutine remind

end module libflit_error
