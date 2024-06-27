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

#define mean_1d_      CONCAT(mean_1d, T)
#define mean_2d_      CONCAT(mean_2d, T)
#define mean_3d_      CONCAT(mean_3d, T)
#define mean_4d_      CONCAT(mean_4d, T)

!
!> Compute the mean of a 1D array
!
function mean_1d_(w, power) result(m)

    TT, dimension(:), intent(in) :: w
    integer, intent(in), optional :: power
    TT :: m

    integer :: p

    if (present(power)) then
        p = power
        m = mTT((sum(nTT(w)**p, mask=(abs(w) /= 0))/size(w))**(1.0d0/p))
    else
        m = mTT(sum(nTT(w))/size(w))
    end if

end function mean_1d_

!
!> Compute the mean of a 2D array
!
function mean_2d_(w, power) result(m)

    TT, dimension(:, :), intent(in) :: w
    integer, intent(in), optional :: power
    TT :: m

    integer :: p

    if (present(power)) then
        p = power
        m = mTT((sum(nTT(w)**p, mask=(abs(w) /= 0))/size(w))**(1.0d0/p))
    else
        m = mTT(sum(nTT(w))/size(w))
    end if

end function mean_2d_

!
!> Compute the mean of a 3D array
!
function mean_3d_(w, power) result(m)

    TT, dimension(:, :, :), intent(in) :: w
    integer, intent(in), optional :: power
    TT :: m

    integer :: p

    if (present(power)) then
        p = power
        m = mTT((sum(nTT(w)**p, mask=(abs(w) /= 0))/size(w))**(1.0d0/p))
    else
        m = mTT(sum(nTT(w))/size(w))
    end if

end function mean_3d_

!
!> Compute the mean of a 4D array
!
function mean_4d_(w, power) result(m)

    TT, dimension(:, :, :, :), intent(in) :: w
    integer, intent(in), optional :: power
    TT :: m

    integer :: p

    if (present(power)) then
        p = power
        m = mTT((sum(nTT(w)**p, mask=(abs(w) /= 0))/size(w))**(1.0d0/p))
    else
        m = mTT(sum(nTT(w))/size(w))
    end if

end function mean_4d_

#undef T
#undef TT
#undef nTT
#undef mTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef mean_1d_
#undef mean_2d_
#undef mean_3d_
#undef mean_4d_
