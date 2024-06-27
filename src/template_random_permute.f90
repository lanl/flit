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

#define random_permute_     CONCAT(random_permute, T)

function random_permute_(w, seed) result(wr)

    TT, dimension(:) :: w
    integer, optional :: seed
    TT, allocatable, dimension(:) :: wr

    integer :: n, i, j

    n = size(w)
    wr = w

    if (n == 1) then
        return
    end if

    if (present(seed)) then
        do i = 1, n - 1
            j = i + irand(range=[0, n - i], seed=seed*i)
            call swap(wr(i), wr(j))
        end do
    else
        do i = 1, n - 1
            j = i + irand(range=[0, n - i])
            call swap(wr(i), wr(j))
        end do
    end if

end function random_permute_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef random_permute_
