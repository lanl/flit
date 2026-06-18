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

#define gst_local_dip_2d_      CONCAT(gst_local_dip_2d, T)
#define gst_local_dip_3d_      CONCAT(gst_local_dip_3d, T)
#define dip_filt_2d_      CONCAT(dip_filt_2d, T)
#define dip_filt_3d_      CONCAT(dip_filt_3d, T)

!
!> Estimate 2D local dip using general structural tensor
!
subroutine gst_local_dip_2d_(w, dip, k1, k2)

    TT, dimension(:, :), intent(in) :: w
    TT, allocatable, dimension(:, :), intent(inout) :: dip
    TT, allocatable, dimension(:, :), intent(inout), optional :: k1, k2

    TT, allocatable, dimension(:, :) :: s11, s12, s22
    type(TTT), allocatable, dimension(:, :) :: ev1, ev2
    TT, allocatable, dimension(:, :) :: ea1, ea2
    integer :: i, j
    integer :: n1, n2

    n1 = size(w, 1)
    n2 = size(w, 2)

    ! Compute structure tensor
    allocate (s11(1:n1, 1:n2))
    allocate (s12(1:n1, 1:n2))
    allocate (s22(1:n1, 1:n2))
    call compute_structure_tensor(w, s11, s12, s22)

    ! Smooth structure tensor with a sigma=1 Gaussian filter
    s11 = gauss_filt(s11, real([1.0, 1.0], fp))
    s12 = gauss_filt(s12, real([1.0, 1.0], fp))
    s22 = gauss_filt(s22, real([1.0, 1.0], fp))

    ! Compute eigenvectors
    allocate (ev1(1:n1, 1:n2))
    allocate (ev2(1:n1, 1:n2))
    allocate (ea1(1:n1, 1:n2))
    allocate (ea2(1:n1, 1:n2))
    call compute_structure_tensor_eigens(s11, s12, s22, ev1, ev2, ea1, ea2)

    call alloc_array(dip, [1, n1, 1, n2])
    if (present(k1)) then
        call alloc_array(k1, [1, n1, 1, n2])
    end if
    if (present(k2)) then
        call alloc_array(k2, [1, n1, 1, n2])
    end if

    ! Compute local dips
    !$omp parallel do private(i,j)
    do j = 1, size(w, 2)
        do i = 1, size(w, 1)
            dip(i, j) = -atan(ev1(i, j)%array(2)/(ev1(i, j)%array(1) + 1.0e-10))
            if (present(k1)) then
                k1(i, j) = ev1(i, j)%array(1)
            end if
            if (present(k2)) then
                k2(i, j) = ev1(i, j)%array(2)
            end if
        end do
    end do
    !$omp end parallel do

end subroutine

!
!> Estimate 3D local dip using general structural tensor
!
subroutine gst_local_dip_3d_(w, dip1, dip2, k1, k2, k3)

    TT, dimension(:, :, :), intent(in) :: w
    TT, allocatable, dimension(:, :, :), intent(inout) :: dip1, dip2
    TT, allocatable, dimension(:, :, :), intent(inout), optional :: k1, k2, k3

    TT, allocatable, dimension(:, :, :) :: s11, s12, s13, s22, s23, s33
    type(TTT), allocatable, dimension(:, :, :) :: ev1, ev2, ev3
    TT, allocatable, dimension(:, :, :) :: ea1, ea2, ea3
    integer :: i, j, k
    integer :: n1, n2, n3

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)

    ! Compute structure tensor
    allocate (s11(1:n1, 1:n2, 1:n3))
    allocate (s12(1:n1, 1:n2, 1:n3))
    allocate (s13(1:n1, 1:n2, 1:n3))
    allocate (s22(1:n1, 1:n2, 1:n3))
    allocate (s23(1:n1, 1:n2, 1:n3))
    allocate (s33(1:n1, 1:n2, 1:n3))
    call compute_structure_tensor(w, s11, s12, s13, s22, s23, s33)

    ! Smooth structure tensor with a sigma=1 Gaussian filter
    s11 = gauss_filt(s11, real([1.0, 1.0, 1.0], fp))
    s12 = gauss_filt(s12, real([1.0, 1.0, 1.0], fp))
    s13 = gauss_filt(s13, real([1.0, 1.0, 1.0], fp))
    s22 = gauss_filt(s22, real([1.0, 1.0, 1.0], fp))
    s23 = gauss_filt(s23, real([1.0, 1.0, 1.0], fp))
    s33 = gauss_filt(s33, real([1.0, 1.0, 1.0], fp))

    ! Compute eigenvectors
    allocate (ev1(1:n1, 1:n2, 1:n3))
    allocate (ev2(1:n1, 1:n2, 1:n3))
    allocate (ev3(1:n1, 1:n2, 1:n3))
    allocate (ea1(1:n1, 1:n2, 1:n3))
    allocate (ea2(1:n1, 1:n2, 1:n3))
    allocate (ea3(1:n1, 1:n2, 1:n3))
    call compute_structure_tensor_eigens(s11, s12, s13, s22, s23, s33, ev1, ev2, ev3, ea1, ea2, ea3)

    call alloc_array(dip1, [1, n1, 1, n2, 1, n3])
    call alloc_array(dip2, [1, n1, 1, n2, 1, n3])
    if (present(k1)) then
        call alloc_array(k1, [1, n1, 1, n2, 1, n3])
    end if
    if (present(k2)) then
        call alloc_array(k2, [1, n1, 1, n2, 1, n3])
    end if
    if (present(k3)) then
        call alloc_array(k3, [1, n1, 1, n2, 1, n3])
    end if

    ! Compute local dips
    !$omp parallel do private(i,j,k)
    do k = 1, size(w, 3)
        do j = 1, size(w, 2)
            do i = 1, size(w, 1)
                dip1(i, j, k) = -atan(ev1(i, j, k)%array(3)/(ev1(i, j, k)%array(1) + 1.0e-10))
                dip2(i, j, k) = -atan(ev1(i, j, k)%array(2)/(ev1(i, j, k)%array(1) + 1.0e-10))
                if (present(k1)) then
                    k1(i, j, k) = ev1(i, j, k)%array(1)
                end if
                if (present(k2)) then
                    k2(i, j, k) = ev1(i, j, k)%array(2)
                end if
                if (present(k3)) then
                    k3(i, j, k) = ev1(i, j, k)%array(3)
                end if
            end do
        end do
    end do
    !$omp end parallel do

end subroutine

!
!> 2D dip filtering based on estimated local dip angles
!
function dip_filt_2d_(w, d, slopes, amps) result(wr)

    TT, dimension(:), intent(in) :: slopes, amps
    TT, dimension(:), intent(in) :: d
    TT, dimension(:, :), intent(in) :: w
    TT, allocatable, dimension(:, :) :: wr

    integer :: n1, n2, i, j, nk1, nk2
    complex, allocatable, dimension(:, :) :: ww
    TT, allocatable, dimension(:) :: cf1, cf2
    TT :: slope(1), amp(1)
    integer :: nsl
    TT :: d1, d2

    n1 = size(w, 1)
    n2 = size(w, 2)
    d1 = d(1)
    d2 = d(2)

    nk1 = max(next_power_235(n1), 180)
    nk2 = max(next_power_235(n2), 180)

    ! Convert to f-k domain
    ww = zeros(nk1, nk2)
    ww(1:n1, 1:n2) = w
    ww = fft(ww)

    ! Slope filtering
    allocate (cf1(1:nk1))
    allocate (cf2(1:nk2))
    cf1 = fft_omega(nk1, d1)
    cf2 = fft_omega(nk2, d2)
    nsl = size(slopes)
    !$omp parallel do private(i, j, slope, amp)
    do j = 1, nk2
        do i = 1, nk1
            if (cf1(i) == 0) then
                slope(1) = cf2(j)/1.0e-10*const_pi
            else
                slope(1) = -cf2(j)/cf1(i)
            end if
            if (slope(1) > slopes(1) .and. slope(1) < slopes(nsl)) then
                amp = ginterp(slopes, amps, slope, method='linear')
            else if (slope(1) <= slopes(1)) then
                amp(1) = amps(1)
            else if (slope(1) >= slopes(nsl)) then
                amp(1) = amps(nsl)
            end if
            ww(i, j) = ww(i, j)*amp(1)
        end do
    end do
    !$omp end parallel do

    ! Convert to t-x domain
    ww = ifft(ww)
    wr = TF(ww(1:n1, 1:n2))

end function

!
!> 3D dip filtering based on estimated local dip angles
!
function dip_filt_3d_(w, d, slopes, amps, plane) result(wr)

    TT, dimension(:), intent(in) :: slopes, amps
    TT, dimension(:), intent(in) :: d
    TT, dimension(:, :, :), intent(in) :: w
    integer, intent(in) :: plane
    TT, allocatable, dimension(:, :, :) :: wr

    TT :: d1, d2, d3
    integer :: i

    d1 = d(1)
    d2 = d(2)
    d3 = d(3)
    wr = w

    select case (plane)
        case (12)
            do i = 1, size(w, 3)
                wr(:, :, i) = dip_filt_2d_(wr(:, :, i), [d1, d2], slopes, amps)
            end do
        case (13)
            do i = 1, size(w, 2)
                wr(:, i, :) = dip_filt_2d_(wr(:, i, :), [d1, d3], slopes, amps)
            end do
        case (23)
            do i = 1, size(w, 1)
                wr(i, :, :) = dip_filt_2d_(wr(i, :, :), [d2, d3], slopes, amps)
            end do
    end select

end function

#undef T
#undef TT
#undef TTT
#undef TF
#undef fp

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef gst_local_dip_2d_
#undef gst_local_dip_3d_
#undef dip_filt_2d_
#undef dip_filt_3d_
