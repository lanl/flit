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

#define const_1d_     CONCAT(const_1d, T)
#define const_2d_     CONCAT(const_2d, T)
#define const_3d_     CONCAT(const_3d, T)
#define const_4d_     CONCAT(const_4d, T)

#define const_like_1d_     CONCAT(const_like_1d, T)
#define const_like_2d_     CONCAT(const_like_2d, T)
#define const_like_3d_     CONCAT(const_like_3d, T)
#define const_like_4d_     CONCAT(const_like_4d, T)

! Constant array
function const_1d_(n1, val) result(w)

    TT, intent(in) :: val
    integer, intent(in) :: n1

    TTT, allocatable, dimension(:) :: w

    allocate (w(1:n1))
    w = val

end function const_1d_

function const_2d_(n1, n2, val) result(w)

    TT, intent(in) :: val
    integer, intent(in) :: n1, n2

    TTT, allocatable, dimension(:, :) :: w

    allocate (w(1:n1, 1:n2))
    w = val

end function const_2d_

function const_3d_(n1, n2, n3, val) result(w)

    TT, intent(in) :: val
    integer, intent(in) :: n1, n2, n3

    TTT, allocatable, dimension(:, :, :) :: w

    allocate (w(1:n1, 1:n2, 1:n3))
    w = val

end function const_3d_

function const_4d_(n1, n2, n3, n4, val) result(w)

    TT, intent(in) :: val
    integer, intent(in) :: n1, n2, n3, n4

    TTT, allocatable, dimension(:, :, :, :) :: w

    allocate (w(1:n1, 1:n2, 1:n3, 1:n4))
    w = val

end function const_4d_

! Constant-like array
function const_like_1d_(w, val) result(wr)

    TT, dimension(:), intent(in) :: w
    TT, intent(in) :: val
    TTT, allocatable, dimension(:) :: wr

    allocate (wr(1:size(w)))
    wr = val

end function const_like_1d_

function const_like_2d_(w, val) result(wr)

    TT, dimension(:, :), intent(in) :: w
    TT, intent(in) :: val
    TTT, allocatable, dimension(:, :) :: wr

    allocate (wr(1:size(w, 1), 1:size(w, 2)))
    wr = val

end function const_like_2d_

function const_like_3d_(w, val) result(wr)

    TT, dimension(:, :, :), intent(in) :: w
    TT, intent(in) :: val
    TTT, allocatable, dimension(:, :, :) :: wr

    allocate (wr(1:size(w, 1), 1:size(w, 2), 1:size(w, 3)))
    wr = val

end function const_like_3d_

function const_like_4d_(w, val) result(wr)

    TT, dimension(:, :, :, :), intent(in) :: w
    TT, intent(in) :: val
    TTT, allocatable, dimension(:, :, :, :) :: wr

    allocate (wr(1:size(w, 1), 1:size(w, 2), 1:size(w, 3), 1:size(w, 4)))
    wr = val

end function const_like_4d_

#undef T
#undef TT
#undef TTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef const_1d_
#undef const_2d_
#undef const_3d_
#undef const_4d_

#undef const_like_1d_
#undef const_like_2d_
#undef const_like_3d_
#undef const_like_4d_

