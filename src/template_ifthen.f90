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

#define ifthen_         CONCAT(ifthen, T)
#define ifthen_1d_      CONCAT(ifthen_1d, T)
#define ifthen_2d_      CONCAT(ifthen_2d, T)
#define ifthen_3d_      CONCAT(ifthen_3d, T)
#define ifthen_4d_      CONCAT(ifthen_4d, T)

function ifthen_(condition, a, b) result(c)

    logical :: condition
    TT :: a, b, c

    if (condition) then
        c = a
    else
        c = b
    end if

end function ifthen_

function ifthen_1d_(condition, a, b) result(c)

    logical :: condition
    TT, dimension(:) :: a, b
    TT, allocatable, dimension(:) :: c

    if (condition) then
        c = a
    else
        c = b
    end if

end function ifthen_1d_

function ifthen_2d_(condition, a, b) result(c)

    logical :: condition
    TT, dimension(:, :) :: a, b
    TT, allocatable, dimension(:, :) :: c

    if (condition) then
        c = a
    else
        c = b
    end if

end function ifthen_2d_

function ifthen_3d_(condition, a, b) result(c)

    logical :: condition
    TT, dimension(:, :, :) :: a, b
    TT, allocatable, dimension(:, :, :) :: c

    if (condition) then
        c = a
    else
        c = b
    end if

end function ifthen_3d_

function ifthen_4d_(condition, a, b) result(c)

    logical :: condition
    TT, dimension(:, :, :, :) :: a, b
    TT, allocatable, dimension(:, :, :, :) :: c

    if (condition) then
        c = a
    else
        c = b
    end if

end function ifthen_4d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef ifthen_
#undef ifthen_1d_
#undef ifthen_2d_
#undef ifthen_3d_
#undef ifthen_4d_
