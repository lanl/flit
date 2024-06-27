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

#define flatten_2d_      CONCAT(flatten_2d, T)
#define flatten_3d_      CONCAT(flatten_3d, T)
#define flatten_4d_      CONCAT(flatten_4d, T)

function flatten_2d_(w) result(wr)

    TT :: w(:, :)
    TT, allocatable :: wr(:)

    allocate (wr(1:size(w, 1)*size(w, 2)))
    wr = reshape(w, [size(wr)])

end function flatten_2d_

function flatten_3d_(w) result(wr)

    TT :: w(:, :, :)
    TT, allocatable :: wr(:)

    allocate (wr(1:size(w, 1)*size(w, 2)*size(w, 3)))
    wr = reshape(w, [size(wr)])

end function flatten_3d_

function flatten_4d_(w) result(wr)

    TT :: w(:, :, :, :)
    TT, allocatable :: wr(:)

    allocate (wr(1:size(w, 1)*size(w, 2)*size(w, 3)*size(w, 4)))
    wr = reshape(w, [size(wr)])

end function flatten_4d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef flatten_2d_
#undef flatten_3d_
#undef flatten_4d_
