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

#define any_in_1d_     CONCAT(any_in_1d, T)
#define all_in_1d_     CONCAT(all_in_1d, T)
#define remove_any_in_1d_     CONCAT(remove_any_in_1d, T)

!function intersect_1d_(a, b) result(c)

!    ! arguments
!    TT, dimension(:) :: a, b
!    TT, allocatable, dimension(:) :: c

!    ! local variables
!    integer :: i


!end subroutine intersect_1d_

function any_in_1d_(a, b) result(y)

    TT, dimension(:) :: a, b
    logical :: y

    integer :: i

    y = .false.
    do i = 1, size(a)
        y = y .or. any(a(i) == b)
    end do

end function any_in_1d_

function all_in_1d_(a, b) result(y)

    TT, dimension(:) :: a, b
    logical :: y

    integer :: i

    y = .true.
    do i = 1, size(a)
        y = y .and. any(a(i) == b)
        if (.not. y) then
            return
        end if
    end do

end function all_in_1d_

function remove_any_in_1d_(a, b) result(c)

    TT, dimension(:) :: a, b
    TTT, allocatable, dimension(:) :: c

    integer :: i

    c = a
    do i = 1, size(b)
        c = pack(c, c /= b(i))
    end do

end function remove_any_in_1d_

#undef T
#undef TT
#undef TTT
#undef DEFAULT_VALUE

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef any_in_1d_
#undef all_in_1d_
#undef remove_any_in_1d_

