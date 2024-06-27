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

#define binarize_1d_      CONCAT(binarize_1d, T)
#define binarize_2d_      CONCAT(binarize_2d, T)
#define binarize_3d_      CONCAT(binarize_3d, T)
#define binarize_4d_      CONCAT(binarize_4d, T)

function binarize_1d_(w, sep, value) result(wt)

    TT, dimension(:), intent(in) :: w
    TT, intent(in) :: sep
    TT, dimension(:), intent(in) :: value
    TT, allocatable, dimension(:) :: wt

    wt = zeros_like(w)
    where (w < sep)
        wt = value(1)
    end where
    where (w >= sep)
        wt = value(2)
    end where

end function binarize_1d_

function binarize_2d_(w, sep, value) result(wt)

    TT, dimension(:, :), intent(in) :: w
    TT, intent(in) :: sep
    TT, dimension(:), intent(in) :: value
    TT, allocatable, dimension(:, :) :: wt

    wt = zeros_like(w)
    where (w < sep)
        wt = value(1)
    end where
    where (w >= sep)
        wt = value(2)
    end where

end function binarize_2d_

function binarize_3d_(w, sep, value) result(wt)

    TT, dimension(:, :, :), intent(in) :: w
    TT, intent(in) :: sep
    TT, dimension(:), intent(in) :: value
    TT, allocatable, dimension(:, :, :) :: wt

    wt = zeros_like(w)
    where (w < sep)
        wt = value(1)
    end where
    where (w >= sep)
        wt = value(2)
    end where

end function binarize_3d_

function binarize_4d_(w, sep, value) result(wt)

    TT, dimension(:, :, :, :), intent(in) :: w
    TT, intent(in) :: sep
    TT, dimension(:), intent(in) :: value
    TT, allocatable, dimension(:, :, :, :) :: wt

    wt = zeros_like(w)
    where (w < sep)
        wt = value(1)
    end where
    where (w >= sep)
        wt = value(2)
    end where

end function binarize_4d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef binarize_1d_
#undef binarize_2d_
#undef binarize_3d_
#undef binarize_4d_
