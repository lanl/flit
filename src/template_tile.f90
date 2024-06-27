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

#define tile_1d_      CONCAT(tile_1d, T)
#define tile_2d_      CONCAT(tile_2d, T)
#define tile_3d_      CONCAT(tile_3d, T)
#define tile_4d_      CONCAT(tile_4d, T)

!
!> Repeat a 1D araay
!
function tile_1d_(w, n) result(wr)

    TT, dimension(:) :: w
    integer :: n
    TT, allocatable, dimension(:) :: wr

    integer :: n1, i

    n1 = size(w)

    allocate(wr(1:n*n1))
    do i = 1, n
        wr((i - 1)*n1 + 1:i*n1) = w
    end do

end function tile_1d_

!
!> Repeat a 2D araay
!
function tile_2d_(w, n) result(wr)

    TT, dimension(:, :) :: w
    integer, dimension(1:2) :: n
    TT, allocatable, dimension(:, :) :: wr

    integer :: n1, n2, i, j

    n1 = size(w, 1)
    n2 = size(w, 2)

    allocate(wr(1:n(1)*n1, 1:n(2)*n2))
    !$omp parallel do private(i, j)
    do j = 1, n(2)
        do i = 1, n(1)
            wr((i - 1)*n1 + 1:i*n1, (j - 1)*n2 + 1:j*n2) = w
        end do
    end do
    !$omp end parallel do

end function tile_2d_

!
!> Repeat a 3D araay
!
function tile_3d_(w, n) result(wr)

    TT, dimension(:, :, :) :: w
    integer, dimension(1:3) :: n
    TT, allocatable, dimension(:, :, :) :: wr

    integer :: n1, n2, n3, i, j, k

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)

    allocate(wr(1:n(1)*n1, 1:n(2)*n2, 1:n(3)*n3))
    !$omp parallel do private(i, j, k)
    do k = 1, n(3)
        do j = 1, n(2)
            do i = 1, n(1)
                wr((i - 1)*n1 + 1:i*n1, (j - 1)*n2 + 1:j*n2, (k - 1)*n3 + 1:k*n3) = w
            end do
        end do
    end do
    !$omp end parallel do

end function tile_3d_

!
!> Repeat a 4D araay
!
function tile_4d_(w, n) result(wr)

    TT, dimension(:, :, :, :) :: w
    integer, dimension(1:4) :: n
    TT, allocatable, dimension(:, :, :, :) :: wr

    integer :: n1, n2, n3, n4, i, j, k, l

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)
    n4 = size(w, 4)

    allocate(wr(1:n(1)*n1, 1:n(2)*n2, 1:n(3)*n3, 1:n(4)*n4))
    !$omp parallel do private(i, j, k, l)
    do l = 1, n(4)
        do k = 1, n(3)
            do j = 1, n(2)
                do i = 1, n(1)
                    wr((i - 1)*n1 + 1:i*n1, (j - 1)*n2 + 1:j*n2, (k - 1)*n3 + 1:k*n3, (l - 1)*n4:l*n4) = w
                end do
            end do
        end do
    end do
    !$omp end parallel do

end function tile_4d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef tile_1d_
#undef tile_2d_
#undef tile_3d_
#undef tile_4d_
