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

#define zeros_like_1d_     CONCAT(zeros_like_1d, T)
#define zeros_like_2d_     CONCAT(zeros_like_2d, T)
#define zeros_like_3d_     CONCAT(zeros_like_3d, T)
#define zeros_like_4d_     CONCAT(zeros_like_4d, T)

#define ones_like_1d_     CONCAT(ones_like_1d, T)
#define ones_like_2d_     CONCAT(ones_like_2d, T)
#define ones_like_3d_     CONCAT(ones_like_3d, T)
#define ones_like_4d_     CONCAT(ones_like_4d, T)

! Zeros-like array
function zeros_like_1d_(w) result(wr)

    TT, dimension(:), intent(in) :: w
    TT, allocatable, dimension(:) :: wr

    allocate (wr(1:size(w)))
    wr = 0

end function zeros_like_1d_

function zeros_like_2d_(w) result(wr)

    TT, dimension(:, :), intent(in) :: w
    TT, allocatable, dimension(:, :) :: wr

    allocate (wr(1:size(w, 1), 1:size(w, 2)))
    wr = 0

end function zeros_like_2d_

function zeros_like_3d_(w) result(wr)

    TT, dimension(:, :, :), intent(in) :: w
    TT, allocatable, dimension(:, :, :) :: wr

    allocate (wr(1:size(w, 1), 1:size(w, 2), 1:size(w, 3)))
    wr = 0

end function zeros_like_3d_

function zeros_like_4d_(w) result(wr)

    TT, dimension(:, :, :, :), intent(in) :: w
    TT, allocatable, dimension(:, :, :, :) :: wr

    allocate (wr(1:size(w, 1), 1:size(w, 2), 1:size(w, 3), 1:size(w, 4)))
    wr = 0

end function zeros_like_4d_

! Ones-like array
function ones_like_1d_(w) result(wr)

    TT, dimension(:), intent(in) :: w
    TT, allocatable, dimension(:) :: wr

    allocate (wr(1:size(w)))
    wr = 1

end function ones_like_1d_

function ones_like_2d_(w) result(wr)

    TT, dimension(:, :), intent(in) :: w
    TT, allocatable, dimension(:, :) :: wr

    allocate (wr(1:size(w, 1), 1:size(w, 2)))
    wr = 1

end function ones_like_2d_

function ones_like_3d_(w) result(wr)

    TT, dimension(:, :, :), intent(in) :: w
    TT, allocatable, dimension(:, :, :) :: wr

    allocate (wr(1:size(w, 1), 1:size(w, 2), 1:size(w, 3)))
    wr = 1

end function ones_like_3d_

function ones_like_4d_(w) result(wr)

    TT, dimension(:, :, :, :), intent(in) :: w
    TT, allocatable, dimension(:, :, :, :) :: wr

    allocate (wr(1:size(w, 1), 1:size(w, 2), 1:size(w, 3), 1:size(w, 4)))
    wr = 1

end function ones_like_4d_

#undef T
#undef TT
#undef TTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef zeros_like_1d_
#undef zeros_like_2d_
#undef zeros_like_3d_
#undef zeros_like_4d_

#undef ones_like_1d_
#undef ones_like_2d_
#undef ones_like_3d_
#undef ones_like_4d_

