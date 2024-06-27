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

#define taper_1d_      CONCAT(taper_1d, T)
#define taper_2d_      CONCAT(taper_2d, T)
#define taper_3d_      CONCAT(taper_3d, T)

!
!> Taper 1D array
!
function taper_1d_(w, len, method, alpha, protect) result(wt)

    TT, dimension(:), intent(in) :: w
    integer, dimension(:), intent(in), optional :: protect, len
    character(len=*), intent(in), dimension(:), optional :: method
    real, dimension(:), intent(in), optional :: alpha
    TT, allocatable, dimension(:) ::  wt

    real, allocatable, dimension(:) :: taper
    integer, dimension(1:2) :: taper_protect, taper_len
    integer :: n, np, i, taper_beg, taper_end
    character(len=24), dimension(1:2) :: taper_method
    real, dimension(1:2) :: taper_alpha

    ! Size of the array
    n = size(w)

    ! Taper parameters
    if (present(len)) then
        taper_len = len
    else
        taper_len = [max(5, nint(0.05*n)), max(5, nint(0.05*n))]
    end if

    if (present(protect)) then
        taper_protect = protect
    else
        taper_protect = [taper_len(1), n - taper_len(2) + 1]
    end if
    taper_protect = clip(taper_protect, 1, n)
    call assert(taper_protect(2) >= taper_protect(1), ' <taper_1d> Error: Protect range specification error.')

    if (present(method)) then
        taper_method = method
    else
        taper_method = ['hann', 'hann']
    end if

    if (present(alpha)) then
        taper_alpha = alpha
    else
        taper_alpha = [0.0, 0.0]
    end if

    ! Setup the taper
    np = taper_len(1) + taper_protect(2) - taper_protect(1) + 1 + taper_len(2) - count(taper_len /= 0)
    taper = taper_window(np, taper_len, taper_method, taper_alpha)

    ! Tapering the array by multiplication
    wt = w

    if (taper_len(1) == 0) then
        taper_beg = taper_protect(1)
    else
        taper_beg = taper_protect(1) - taper_len(1) + 1
    end if
    if (taper_len(2) == 0) then
        taper_end = taper_protect(2)
    else
        taper_end = taper_protect(2) + taper_len(2) - 1
    end if

    wt(1:taper_beg - 1) = 0.0
    do i = max(1, taper_beg), min(n, taper_end)
        wt(i) = wt(i)*taper(i - taper_beg + 1)
    end do
    wt(taper_end + 1:n) = 0.0

end function taper_1d_

!
!> Taper 2D array
!
function taper_2d_(w, len, method, alpha, protect) result(wt)

    TT, dimension(:, :) :: w
    integer, dimension(:), optional :: protect, len
    character(len=*), dimension(:), optional :: method
    real, dimension(:), optional :: alpha
    TT, allocatable, dimension(:, :) :: wt

    integer :: n1, n2, i, j
    integer, dimension(1:4) :: taper_protect, taper_len
    character(len=24), dimension(1:4) :: taper_method
    real, dimension(1:4) :: taper_alpha

    ! dimension of arrays
    n1 = size(w, 1)
    n2 = size(w, 2)

    ! Taper parameters
    if (present(len)) then
        taper_len = len
    else
        taper_len(1:2) = max(5, nint(0.05*n1))
        taper_len(3:4) = max(5, nint(0.05*n2))
    end if

    if (present(protect)) then
        taper_protect = protect
    else
        taper_protect(1) = taper_len(1)
        taper_protect(2) = n1 - taper_len(2) + 1
        taper_protect(3) = taper_len(3)
        taper_protect(4) = n2 - taper_len(4) + 1
    end if
    taper_protect(1:2) = clip(taper_protect(1:2), 1, n1)
    taper_protect(3:4) = clip(taper_protect(3:4), 1, n2)

    if (present(method)) then
        taper_method = method
    else
        taper_method = ['hann', 'hann', 'hann', 'hann']
    end if

    if (present(alpha)) then
        taper_alpha = real(alpha)
    else
        taper_alpha = [0.0, 0.0, 0.0, 0.0]
    end if

    ! Tapering the array by multiplication
    wt = w

    ! dim = 1 tapering
    !$omp parallel do private(j)
    do j = 1, n2
        wt(:, j) = taper_1d_(w(:, j), &
            taper_len(1:2), taper_method(1:2), taper_alpha(1:2), taper_protect(1:2))
    end do
    !$omp end parallel do

    ! dim = 2 tapering
    !$omp parallel do private(i)
    do i = 1, n1
        wt(i, :) = taper_1d_(wt(i, :), &
            taper_len(3:4), taper_method(3:4), taper_alpha(3:4), taper_protect(3:4))
    end do
    !$omp end parallel do

end function taper_2d_

!
!> Taper 3D array
!
function taper_3d_(w, len, method, alpha, protect) result(wt)

    TT, dimension(:, :, :) :: w
    integer, dimension(:), optional :: protect, len
    character(len=*), dimension(:), optional :: method
    real, dimension(:), optional :: alpha
    TT, allocatable, dimension(:, :, :) :: wt

    integer :: n1, n2, n3, i, j, k
    integer, dimension(1:6) :: taper_protect, taper_len
    character(len=24), dimension(1:6) :: taper_method
    real, dimension(1:6) :: taper_alpha

    ! dimension of arrays
    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)

    ! Taper parameters
    if (present(len)) then
        taper_len = len
    else
        taper_len(1:2) = max(5, nint(0.05*n1))
        taper_len(3:4) = max(5, nint(0.05*n2))
        taper_len(5:6) = max(5, nint(0.05*n3))
    end if

    if (present(protect)) then
        taper_protect = protect
    else
        taper_protect(1) = taper_len(1)
        taper_protect(2) = n1 - taper_len(2) + 1
        taper_protect(3) = taper_len(3)
        taper_protect(4) = n2 - taper_len(4) + 1
        taper_protect(5) = taper_len(5)
        taper_protect(6) = n3 - taper_len(6) + 1
    end if
    taper_protect(1:2) = clip(taper_protect(1:2), 1, n1)
    taper_protect(3:4) = clip(taper_protect(3:4), 1, n2)
    taper_protect(5:6) = clip(taper_protect(5:6), 1, n3)

    if (present(method)) then
        taper_method = method
    else
        taper_method = ['hann', 'hann', 'hann', 'hann', 'hann', 'hann']
    end if
    if (present(alpha)) then
        taper_alpha = real(alpha)
    else
        taper_alpha = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    end if

    ! Tapering the array by multiplication
    wt = w

    ! dim = 3 tapering
    !$omp parallel do private(i, j) collapse(2)
    do i = 1, n1
        do j = 1, n2
            wt(i, j, :) = taper_1d_(w(i, j, :), &
                taper_len(5:6), taper_method(5:6), taper_alpha(5:6), taper_protect(5:6))

        end do
    end do
    !$omp end parallel do

    ! dim = 2 tapering
    !$omp parallel do private(i, k) collapse(2)
    do i = 1, n1
        do k = 1, n3
            wt(i, :, k) = taper_1d_(wt(i, :, k), &
                taper_len(3:4), taper_method(3:4), taper_alpha(3:4), taper_protect(3:4))
        end do
    end do
    !$omp end parallel do

    ! dim = 1 tapering
    !$omp parallel do private(j, k) collapse(2)
    do j = 1, n2
        do k = 1, n3
            wt(:, j, k) = taper_1d_(wt(:, j, k), &
                taper_len(1:2), taper_method(1:2), taper_alpha(1:2), taper_protect(1:2))
        end do
    end do
    !$omp end parallel do

end function taper_3d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef taper_1d_
#undef taper_2d_
#undef taper_3d_
