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

#define slice_1d_      CONCAT(slice_1d, T)
#define slice_2d_      CONCAT(slice_2d, T)
#define slice_3d_      CONCAT(slice_3d, T)

!
!> Take a slice from a 1D array
!
!> @param[in] w The input 1D array in
!> @param[in] index The index of slicing, must be 1 <= index <= dim_size
!> @return wr The scalar sliced from the 1D array
!
function slice_1d_(w, index) result(wr)

    TT, dimension(:) :: w
    integer :: index
    TT :: wr

    wr = w(index)

end function slice_1d_

!
!> Take a slice from a 2D array
!
!> @param[in] w The input 2D array in
!> @param[in] dim Which dimension to slice
!> @param[in] index The index of slicing, must be 1 <= index <= dim_size
!> @return wr The 1D array sliced from the 2D array
!
function slice_2d_(w, dim, index) result(wr)

    TT, dimension(:, :) :: w
    integer :: dim, index
    TT, allocatable, dimension(:) :: wr

    integer :: n1, n2

    n1 = size(w, 1)
    n2 = size(w, 2)

    select case (dim)
        case (1)
            allocate (wr(1:n2))
            wr(:) = w(index, :)
        case (2)
            allocate (wr(1:n1))
            wr(:) = w(:, index)
    end select

end function slice_2d_

!
!> Take a slice from a 3D array
!
!> @param[in] w The input 3D array in
!> @param[in] dim Which dimension to slice
!> @param[in] index The index of slicing, must be 1 <= index <= dim_size
!> @return wr The 2D array sliced from the 3D array
!
function slice_3d_(w, dim, index) result(wr)

    TT, dimension(:, :, :) :: w
    integer :: dim, index
    TT, allocatable, dimension(:, :) :: wr

    integer :: n1, n2, n3

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)

    select case (dim)
        case (1)
            allocate (wr(1:n2, 1:n3))
            wr(:, :) = w(index, :, :)
        case (2)
            allocate (wr(1:n1, 1:n3))
            wr(:, :) = w(:, index, :)
        case (3)
            allocate (wr(1:n1, 1:n2))
            wr(:, :) = w(:, :, index)
    end select

end function slice_3d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef slice_1d_
#undef slice_2d_
#undef slice_3d_
