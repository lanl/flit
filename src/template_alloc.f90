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

#define alloc_array_1d_     CONCAT(alloc_array_1d, T)
#define alloc_array_2d_     CONCAT(alloc_array_2d, T)
#define alloc_array_3d_     CONCAT(alloc_array_3d, T)
#define alloc_array_4d_     CONCAT(alloc_array_4d, T)
#define alloc_large_array_1d_     CONCAT(alloc_large_array_1d, T)

subroutine alloc_array_1d_(array, n, pad, source)

    TT, allocatable, dimension(:), intent(inout) :: array
    integer, dimension(1:2), intent(in) :: n
    integer, intent(in), optional :: pad
    TT, dimension(:), intent(in), optional :: source

    integer :: l
    TTT, allocatable, dimension(:) :: w

    if (present(source)) then
        allocate (w(size(source)), source=source)
    end if

    if (present(pad)) then
        l = pad
    else
        l = 0
    end if

    if (allocated(array)) then
        deallocate (array)
    end if

    allocate (array(n(1) - l:n(2) + l))

    if (present(source)) then
        array(:) = w(:)
    else
        array = DEFAULT_VALUE
    end if

end subroutine alloc_array_1d_

subroutine alloc_large_array_1d_(array, n, pad, source)

    TT, allocatable, dimension(:), intent(inout) :: array
    integer(kind=8), dimension(1:2), intent(in) :: n
    integer, intent(in), optional :: pad
    TT, dimension(:), intent(in), optional :: source

    integer :: l
    TTT, allocatable, dimension(:) :: w

    if (present(source)) then
        allocate (w(size(source)), source=source)
    end if

    if (present(pad)) then
        l = pad
    else
        l = 0
    end if

    if (allocated(array)) then
        deallocate (array)
    end if

    allocate (array(n(1) - l:n(2) + l))

    if (present(source)) then
        array(:) = w(:)
    else
        array = DEFAULT_VALUE
    end if

end subroutine alloc_large_array_1d_

subroutine alloc_array_2d_(array, n, pad, source)

    TT, allocatable, dimension(:, :), intent(inout) :: array
    integer, dimension(1:4), intent(in) :: n
    integer, intent(in), optional :: pad
    TT, dimension(:, :), intent(in), optional :: source

    integer :: l
    TTT, allocatable, dimension(:, :) :: w

    if (present(source)) then
        allocate (w(size(source, 1), size(source, 2)), source=source)
    end if

    if (present(pad)) then
        l = pad
    else
        l = 0
    end if

    if (allocated(array)) then
        deallocate (array)
    end if

    allocate (array(n(1) - l:n(2) + l, n(3) - l:n(4) + l))

    if (present(source)) then
        array(:, :) = w(:, :)
    else
        array = DEFAULT_VALUE
    end if

end subroutine alloc_array_2d_

subroutine alloc_array_3d_(array, n, pad, source)

    TT, allocatable, dimension(:, :, :), intent(inout) :: array
    integer, dimension(1:6), intent(in) :: n
    integer, intent(in), optional :: pad
    TT, dimension(:, :, :), intent(in), optional :: source

    integer :: l
    TTT, allocatable, dimension(:, :, :) :: w

    if (present(source)) then
        allocate (w(size(source, 1), size(source, 2), size(source, 3)), source=source)
    end if

    if (present(pad)) then
        l = pad
    else
        l = 0
    end if

    if (allocated(array)) then
        deallocate (array)
    end if

    allocate (array(n(1) - l:n(2) + l, n(3) - l:n(4) + l, n(5) - l:n(6) + l))

    if (present(source)) then
        array(:, :, :) = w(:, :, :)
    else
        array = DEFAULT_VALUE
    end if

end subroutine alloc_array_3d_


subroutine alloc_array_4d_(array, n, pad, source)

    TT, allocatable, dimension(:, :, :, :), intent(inout) :: array
    integer, dimension(1:8), intent(in) :: n
    integer, intent(in), optional :: pad
    TT, dimension(:, :, :, :), intent(in), optional :: source

    integer :: l
    TTT, allocatable, dimension(:, :, :, :) :: w

    if (present(source)) then
        allocate (w(size(source, 1), size(source, 2), size(source, 3), size(source, 4)), source=source)
    end if

    if (present(pad)) then
        l = pad
    else
        l = 0
    end if

    if (allocated(array)) then
        deallocate (array)
    end if

    allocate (array(n(1) - l:n(2) + l, n(3) - l:n(4) + l, n(5) - l:n(6) + l, n(7) - l:n(8) + l))

    if (present(source)) then
        array = w(:, :, :, :)
    else
        array = DEFAULT_VALUE
    end if

end subroutine alloc_array_4d_

#undef T
#undef TT
#undef TTT
#undef DEFAULT_VALUE

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef alloc_array_1d_
#undef alloc_array_2d_
#undef alloc_array_3d_
#undef alloc_array_4d_
#undef alloc_large_array_1d_
