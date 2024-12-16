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

#define ifelse_         CONCAT(ifelse, T)
#define ifelse_1d_      CONCAT(ifelse_1d, T)
#define ifelse_2d_      CONCAT(ifelse_2d, T)
#define ifelse_3d_      CONCAT(ifelse_3d, T)
#define ifelse_4d_      CONCAT(ifelse_4d, T)

function ifelse_(condition, a, b) result(c)

    logical :: condition
    TT :: a, b, c

    if (condition) then
        c = a
    else
        c = b
    end if

end function ifelse_

function ifelse_1d_(condition, a, b) result(c)

    logical :: condition
    TT, dimension(:) :: a, b
    TT, allocatable, dimension(:) :: c

    if (condition) then
        c = a
    else
        c = b
    end if

end function ifelse_1d_

function ifelse_2d_(condition, a, b) result(c)

    logical :: condition
    TT, dimension(:, :) :: a, b
    TT, allocatable, dimension(:, :) :: c

    if (condition) then
        c = a
    else
        c = b
    end if

end function ifelse_2d_

function ifelse_3d_(condition, a, b) result(c)

    logical :: condition
    TT, dimension(:, :, :) :: a, b
    TT, allocatable, dimension(:, :, :) :: c

    if (condition) then
        c = a
    else
        c = b
    end if

end function ifelse_3d_

function ifelse_4d_(condition, a, b) result(c)

    logical :: condition
    TT, dimension(:, :, :, :) :: a, b
    TT, allocatable, dimension(:, :, :, :) :: c

    if (condition) then
        c = a
    else
        c = b
    end if

end function ifelse_4d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef ifelse_
#undef ifelse_1d_
#undef ifelse_2d_
#undef ifelse_3d_
#undef ifelse_4d_
