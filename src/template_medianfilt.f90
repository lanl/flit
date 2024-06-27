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

#define median_filt_1d_      CONCAT(median_filt_1d, T)
#define median_filt_2d_      CONCAT(median_filt_2d, T)
#define median_filt_3d_      CONCAT(median_filt_3d, T)

function median_filt_1d_(w, nw) result(wr)

    TT, dimension(:), intent(in) :: w
    integer, intent(in), optional :: nw

    integer :: i, n1, nw1
    TT, allocatable, dimension(:) :: wr

    n1 = size(w)
    if (present(nw)) then
        nw1 = nw
    else
        nw1 = 1
    end if

    allocate(wr(1:n1))

    do i = 1, n1
        wr(i) = median(w(max(1, i - nw1):min(n1, i + nw1)))
    end do

end function median_filt_1d_

function median_filt_2d_(w, nw) result(wr)

    TT, dimension(:, :), intent(in) :: w
    integer, dimension(:), intent(in), optional :: nw

    integer :: i, j, n1, n2, nw1, nw2
    TT, allocatable, dimension(:, :) :: wr

    n1 = size(w, 1)
    n2 = size(w, 2)
    if (present(nw)) then
        nw1 = nw(1)
        nw2 = nw(2)
    else
        nw1 = 1
        nw2 = 1
    end if

    allocate(wr(1:n1, 1:n2))

    !$omp parallel do private(i, j)
    do j = 1, n2
        do i = 1, n1
            wr(i, j) = median(w(max(1, i - nw1):min(n1, i + nw1), &
                max(1, j - nw2):min(n2, j + nw2)))
        end do
    end do
    !$omp end parallel do

end function median_filt_2d_

function median_filt_3d_(w, nw) result(wr)

    TT, dimension(:, :, :), intent(in) :: w
    integer, dimension(:), intent(in), optional :: nw

    integer :: i, j, k, n1, n2, n3, nw1, nw2, nw3
    TT, allocatable, dimension(:, :, :) :: wr

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)
    if (present(nw)) then
        nw1 = nw(1)
        nw2 = nw(2)
        nw3 = nw(3)
    else
        nw1 = 1
        nw2 = 1
        nw3 = 1
    end if

    allocate (wr(1:n1, 1:n2, 1:n3))

    !$omp parallel do private(i, j, k)
    do k = 1, n3
        do j = 1, n2
            do i = 1, n1
                wr(i, j, k) = median(w( &
                    max(1, i - nw1):min(n1, i + nw1), &
                    max(1, j - nw2):min(n2, j + nw2), &
                    max(1, k - nw3):min(n3, k + nw3)))
            end do
        end do
    end do
    !$omp end parallel do

end function median_filt_3d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef median_filt_1d_
#undef median_filt_2d_
#undef median_filt_3d_
