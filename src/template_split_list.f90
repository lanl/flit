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

#define split_list_         CONCAT(split_list, T)

!
!> Split a 1D array into small chunks so that each
!> chunk has approximately similar number of elements
!> and get a specific chunck by index
!
function split_list_(list, nc, index) result(c)

    ! arguments
    TT, dimension(:), intent(in) :: list
    integer, intent(in) :: nc
    integer, intent(in) :: index
    TT, allocatable, dimension(:) :: c

    integer, allocatable, dimension(:, :) :: irange
    integer :: n

    irange = zeros(nc, 2)

    n = size(list)
    call assert(index >= 1 .and. index <= nc, ' <split_list> Error: Chunck index is out of range. ')

    call cut(1, n, nc, irange)
    c = list(irange(index, 1):irange(index, 2))

end function split_list_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef split_list_
