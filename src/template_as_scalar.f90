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

#define as_scalar_1d_     CONCAT(as_scalar_1d, T)
#define as_scalar_2d_     CONCAT(as_scalar_2d, T)
#define as_scalar_3d_     CONCAT(as_scalar_3d, T)
#define as_scalar_4d_     CONCAT(as_scalar_4d, T)

function as_scalar_1d_(w, i) result(s)

    TT, intent(in), dimension(:) :: w
    integer, intent(in), optional :: i
    TT :: s

    if (present(i)) then
        s = w(i)
    else
        s = w(1)
    end if

end function as_scalar_1d_

function as_scalar_2d_(w, i) result(s)

    TT, intent(in), dimension(:, :) :: w
    integer, dimension(:), intent(in), optional :: i
    TT :: s

    if (present(i)) then
        s = w(i(1), i(2))
    else
        s = w(1, 1)
    end if

end function as_scalar_2d_

function as_scalar_3d_(w, i) result(s)

    TT, intent(in), dimension(:, :, :) :: w
    integer, dimension(:), intent(in), optional :: i
    TT :: s

    if (present(i)) then
        s = w(i(1), i(2), i(3))
    else
        s = w(1, 1, 1)
    end if

end function as_scalar_3d_

function as_scalar_4d_(w, i) result(s)

    TT, intent(in), dimension(:, :, :, :) :: w
    integer, dimension(:), intent(in), optional :: i
    TT :: s

    if (present(i)) then
        s = w(i(1), i(2), i(3), i(4))
    else
        s = w(1, 1, 1, 1)
    end if

end function as_scalar_4d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef as_scalar_1d_
#undef as_scalar_2d_
#undef as_scalar_3d_
#undef as_scalar_4d_
