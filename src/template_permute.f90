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

#define permute_array_3d_     CONCAT(permute_array_3d, T)
#define permute_array_4d_     CONCAT(permute_array_4d, T)
#define permute_3d_     CONCAT(permute_3d, T)
#define permute_4d_     CONCAT(permute_4d, T)

subroutine permute_array_3d_(w, order)

    TT, allocatable, dimension(:, :, :), intent(inout) :: w
    integer, intent(in) :: order

    integer, allocatable, dimension(:) :: sorder
    integer :: n1, n2, n3

    allocate (sorder(1:3))
    sorder = breakint(order)

    n1 = size(w, sorder(1))
    n2 = size(w, sorder(2))
    n3 = size(w, sorder(3))

    w = reshape(w, shape=[n1, n2, n3], &
        order=[ &
        maxloc(sorder, mask=(sorder == 1)), &
        maxloc(sorder, mask=(sorder == 2)), &
        maxloc(sorder, mask=(sorder == 3))])

end subroutine permute_array_3d_

subroutine permute_array_4d_(w, order)

    TT, allocatable, dimension(:, :, :, :), intent(inout) :: w
    integer, intent(in) :: order

    integer, allocatable, dimension(:) :: sorder
    integer :: n1, n2, n3, n4

    allocate (sorder(1:4))
    sorder = breakint(order)

    n1 = size(w, sorder(1))
    n2 = size(w, sorder(2))
    n3 = size(w, sorder(3))
    n4 = size(w, sorder(4))

    w = reshape(w, shape=[n1, n2, n3, n4], &
        order=[ &
        maxloc(sorder, mask=(sorder == 1)), &
        maxloc(sorder, mask=(sorder == 2)), &
        maxloc(sorder, mask=(sorder == 3)), &
        maxloc(sorder, mask=(sorder == 4))])

end subroutine permute_array_4d_

function permute_3d_(w, order) result(pw)

    TT, dimension(:, :, :), intent(in) :: w
    integer, intent(in) :: order
    TT, allocatable, dimension(:, :, :) :: pw

    integer, allocatable, dimension(:) :: sorder
    integer :: n1, n2, n3

    allocate (sorder(1:3))
    sorder = breakint(order)

    n1 = size(w, sorder(1))
    n2 = size(w, sorder(2))
    n3 = size(w, sorder(3))

    allocate (pw(1:n1, 1:n2, 1:n3))

    pw = reshape(w, shape=shape(pw), &
        order=[ &
        maxloc(sorder, mask=(sorder == 1)), &
        maxloc(sorder, mask=(sorder == 2)), &
        maxloc(sorder, mask=(sorder == 3))])

end function permute_3d_

function permute_4d_(w, order) result(pw)

    TT, dimension(:, :, :, :), intent(in) :: w
    integer, intent(in) :: order
    TT, allocatable, dimension(:, :, :, :) :: pw

    integer, allocatable, dimension(:) :: sorder
    integer :: n1, n2, n3, n4

    allocate (sorder(1:4))
    sorder = breakint(order)

    n1 = size(w, sorder(1))
    n2 = size(w, sorder(2))
    n3 = size(w, sorder(3))
    n4 = size(w, sorder(4))

    allocate (pw(1:n1, 1:n2, 1:n3, 1:n4))

    pw = reshape(w, shape=shape(pw), &
        order=[ &
        maxloc(sorder, mask=(sorder == 1)), &
        maxloc(sorder, mask=(sorder == 2)), &
        maxloc(sorder, mask=(sorder == 3)), &
        maxloc(sorder, mask=(sorder == 4))])

    !            order=[ &
        !            findloc(sorder,1), &
        !            findloc(sorder,2), &
        !            findloc(sorder,3), &
        !            findloc(sorder,4)])

end function permute_4d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef permute_array_3d_
#undef permute_array_4d_
#undef permute_3d_
#undef permute_4d_
