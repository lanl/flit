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

#define cross_product_2d_      CONCAT(cross_product_2d, T)
#define cross_product_3d_      CONCAT(cross_product_3d, T)

function cross_product_2d_(a, b) result(c)

    TT, dimension(1:2), intent(in) :: a, b
    TT :: c

    c = a(1)*b(2) - a(2)*b(1)

end function cross_product_2d_

function cross_product_3d_(a, b) result(c)

    TT, dimension(1:3), intent(in) :: a, b
    TT, allocatable, dimension(:) :: c

    allocate(c(1:3))

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

end function cross_product_3d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef cross_product_2d_
#undef cross_product_3d_
