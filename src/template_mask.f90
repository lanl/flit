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

#define mask_array_1d_      CONCAT(mask_array_1d, T)
#define mask_array_2d_      CONCAT(mask_array_2d, T)
#define mask_array_3d_      CONCAT(mask_array_3d, T)
#define mask_array_4d_      CONCAT(mask_array_4d, T)

#define mask_1d_      CONCAT(mask_1d, T)
#define mask_2d_      CONCAT(mask_2d, T)
#define mask_3d_      CONCAT(mask_3d, T)
#define mask_4d_      CONCAT(mask_4d, T)

subroutine mask_array_1d_(w, mask, zero_to_nan)

    TT, dimension(:), intent(inout) :: w
    TT, dimension(:), intent(in) :: mask
    logical, intent(in), optional :: zero_to_nan

    w = w*mask

    if (present(zero_to_nan)) then
        if (zero_to_nan) then
            where (mask == 0)
                w = nan()
            end where
        end if
    end if

end subroutine mask_array_1d_

subroutine mask_array_2d_(w, mask, zero_to_nan)

    TT, dimension(:, :), intent(inout) :: w
    TT, dimension(:, :), intent(in) :: mask
    logical, intent(in), optional :: zero_to_nan

    w = w*mask

    if (present(zero_to_nan)) then
        if (zero_to_nan) then
            where (mask == 0)
                w = nan()
            end where
        end if
    end if

end subroutine mask_array_2d_

subroutine mask_array_3d_(w, mask, zero_to_nan)

    TT, dimension(:, :, :), intent(inout) :: w
    TT, dimension(:, :, :), intent(in) :: mask
    logical, intent(in), optional :: zero_to_nan

    w = w*mask

    if (present(zero_to_nan)) then
        if (zero_to_nan) then
            where (mask == 0)
                w = nan()
            end where
        end if
    end if

end subroutine mask_array_3d_

subroutine mask_array_4d_(w, mask, zero_to_nan)

    TT, dimension(:, :, :, :), intent(inout) :: w
    TT, dimension(:, :, :, :), intent(in) :: mask
    logical, intent(in), optional :: zero_to_nan

    w = w*mask

    if (present(zero_to_nan)) then
        if (zero_to_nan) then
            where (mask == 0)
                w = nan()
            end where
        end if
    end if

end subroutine mask_array_4d_

function mask_1d_(w, mask, zero_to_nan) result(wm)

    TT, dimension(:), intent(in) :: w, mask
    logical, intent(in), optional :: zero_to_nan

    TT, allocatable, dimension(:) :: wm

    call alloc_array(wm, [1, size(w)])

    wm = w*mask

    if (present(zero_to_nan)) then
        if (zero_to_nan) then
            where (mask == 0)
                wm = nan()
            end where
        end if
    end if

end function mask_1d_

function mask_2d_(w, mask, zero_to_nan) result(wm)

    TT, dimension(:, :), intent(in) :: w, mask
    logical, intent(in), optional :: zero_to_nan

    TT, allocatable, dimension(:, :) :: wm

    call alloc_array(wm, [1, size(w), 1, size(w, 2)])

    wm = w*mask

    if (present(zero_to_nan)) then
        if (zero_to_nan) then
            where (mask == 0)
                wm = nan()
            end where
        end if
    end if

end function mask_2d_

function mask_3d_(w, mask, zero_to_nan) result(wm)

    TT, dimension(:, :, :), intent(in) :: w, mask
    logical, intent(in), optional :: zero_to_nan

    TT, allocatable, dimension(:, :, :) :: wm

    call alloc_array(wm, [1, size(w), 1, size(w, 2), 1, size(w, 3)])

    wm = w*mask

    if (present(zero_to_nan)) then
        if (zero_to_nan) then
            where (mask == 0)
                wm = nan()
            end where
        end if
    end if

end function mask_3d_

function mask_4d_(w, mask, zero_to_nan) result(wm)

    TT, dimension(:, :, :, :), intent(in) :: w, mask
    logical, intent(in), optional :: zero_to_nan

    TT, allocatable, dimension(:, :, :, :) :: wm

    call alloc_array(wm, [1, size(w), 1, size(w, 2), 1, size(w, 3), 1, size(w, 4)])

    wm = w*mask

    if (present(zero_to_nan)) then
        if (zero_to_nan) then
            where (mask == 0)
                wm = nan()
            end where
        end if
    end if

end function mask_4d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef mask_array_1d_
#undef mask_array_2d_
#undef mask_array_3d_
#undef mask_array_4d_

#undef mask_1d_
#undef mask_2d_
#undef mask_3d_
#undef mask_4d_
