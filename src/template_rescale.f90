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

#define rescale_1d_      CONCAT(rescale_1d, T)
#define rescale_2d_      CONCAT(rescale_2d, T)
#define rescale_3d_      CONCAT(rescale_3d, T)
#define rescale_4d_      CONCAT(rescale_4d, T)

function rescale_1d_(w, range) result(ws)

    TT, dimension(:), intent(in) :: w
    TT, dimension(:), intent(in) :: range

    TT :: dr
    TT, allocatable, dimension(:) :: ws

    ws = w

    if (maxval(ws) /= minval(ws) .and. range(1) /= range(2)) then
        dr = range(2) - range(1)
        ws = ws - minval(ws)
        ws = ws/maxval(ws)
        ws = ws*dr + range(1)
    else
        ws = range(1)
    end if

end function rescale_1d_

function rescale_2d_(w, range) result(ws)

    TT, dimension(:, :), intent(in) :: w
    TT, dimension(:), intent(in) :: range

    TT :: dr
    TT, allocatable, dimension(:, :) :: ws

    ws = w

    if (maxval(ws) /= minval(ws) .and. range(1) /= range(2)) then
        dr = range(2) - range(1)
        ws = ws - minval(ws)
        ws = ws/maxval(ws)
        ws = ws*dr + range(1)
    else
        ws = range(1)
    end if

end function rescale_2d_

function rescale_3d_(w, range) result(ws)

    TT, dimension(:, :, :), intent(in) :: w
    TT, dimension(:), intent(in) :: range

    TT :: dr
    TT, allocatable, dimension(:, :, :) :: ws

    ws = w

    if (maxval(ws) /= minval(ws) .and. range(1) /= range(2)) then
        dr = range(2) - range(1)
        ws = ws - minval(ws)
        ws = ws/maxval(ws)
        ws = ws*dr + range(1)
    else
        ws = range(1)
    end if

end function rescale_3d_

function rescale_4d_(w, range) result(ws)

    TT, dimension(:, :, :, :), intent(in) :: w
    TT, dimension(:), intent(in) :: range

    TT :: dr
    TT, allocatable, dimension(:, :, :, :) :: ws

    ws = w

    if (maxval(ws) /= minval(ws) .and. range(1) /= range(2)) then
        dr = range(2) - range(1)
        ws = ws - minval(ws)
        ws = ws/maxval(ws)
        ws = ws*dr + range(1)
    else
        ws = range(1)
    end if

end function rescale_4d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef rescale_1d_
#undef rescale_2d_
#undef rescale_3d_
#undef rescale_4d_
