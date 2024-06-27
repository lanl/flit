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

#define eye_2d_     CONCAT(eye_2d, T)
#define eye_2d_from_vector_     CONCAT(eye_2d_from_vector, T)
#define get_diagonal_2d_     CONCAT(get_diagonal_2d, T)

function eye_2d_(n, val) result(w)

    integer :: n
    TT :: val
    TT, allocatable, dimension(:, :) :: w

    integer :: i

    allocate (w(1:n, 1:n))
    w = DEFAULT_VALUE
    !$omp parallel do private(i)
    do i = 1, n
        w(i, i) = val
    end do
    !$omp end parallel do

end function eye_2d_

function eye_2d_from_vector_(v) result(w)

    TT, dimension(:) :: v
    TT, allocatable, dimension(:, :) :: w

    integer :: i, n

    n = size(v)

    allocate (w(1:n, 1:n))
    w = DEFAULT_VALUE
    !$omp parallel do private(i)
    do i = 1, n
        w(i, i) = v(i)
    end do
    !    forall(i = 1:n) w(i, i) = v(i)
    !$omp end parallel do

end function eye_2d_from_vector_

function get_diagonal_2d_(w) result(d)

    TT, dimension(:, :), intent(in) :: w
    TT, allocatable, dimension(:) :: d

    integer :: n, i

    n = min(size(w, 1), size(w, 2))
    allocate (d(1:n))
    !$omp parallel do private(i)
    do i = 1, n
        d(i) = w(i, i)
    end do
    !$omp end parallel do
    !    forall(i = 1:n) d(i) = w(i, i)

end function get_diagonal_2d_

#undef T
#undef TT
#undef DEFAULT_VALUE

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef eye_2d_
#undef eye_2d_from_vector_
#undef get_diagonal_2d_
