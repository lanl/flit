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

module libflit_dip

    use libflit_array
    use libflit_array_extension
    use libflit_array_operation
    use libflit_gaussfilt
    use libflit_utility
    use libflit_andffilt
    use libflit_random
    use libflit_transform
    use libflit_interp

    implicit none

    interface gstdip
        module procedure :: gst_localdip2
        module procedure :: gst_localdip3
    end interface gstdip

    interface dip_filt
        module procedure :: dip_filt_2d_float
        module procedure :: dip_filt_3d_float
        module procedure :: dip_filt_2d_complex
    end interface dip_filt

    private
    public :: gstdip
    public :: dip_filt

contains

    !
    !> Estimate 2D local dip using general structural tensor
    !
    subroutine gst_localdip2(w, dp, k1, k2)

        ! arguments
        real, allocatable, dimension(:, :), intent(in) :: w
        real, allocatable, dimension(:, :), intent(inout) :: dp
        real, allocatable, dimension(:, :), intent(inout), optional :: k1, k2

        real, allocatable, dimension(:, :) :: s11, s12, s22
        type(meta_array1_real), allocatable, dimension(:, :) :: ev1, ev2
        real, allocatable, dimension(:, :) :: ea1, ea2
        integer :: i, j
        integer :: n1, n2

        n1 = size(w, 1)
        n2 = size(w, 2)

        ! compute structure tensor
        allocate (s11(1:n1, 1:n2))
        allocate (s12(1:n1, 1:n2))
        allocate (s22(1:n1, 1:n2))
        call compute_structure_tensor(w, s11, s12, s22)

        ! smooth structure tensor with a sigma=1 Gaussian filter
        s11 = gauss_filt(s11, [1.0, 1.0])
        s12 = gauss_filt(s12, [1.0, 1.0])
        s22 = gauss_filt(s22, [1.0, 1.0])

        ! compute eigenvectors
        allocate (ev1(1:n1, 1:n2))
        allocate (ev2(1:n1, 1:n2))
        allocate (ea1(1:n1, 1:n2))
        allocate (ea2(1:n1, 1:n2))
        call compute_structure_tensor_eigens(s11, s12, s22, ev1, ev2, ea1, ea2)

        call alloc_array(dp, [1, n1, 1, n2])
        if (present(k1)) then
            call alloc_array(k1, [1, n1, 1, n2])
        end if
        if (present(k2)) then
            call alloc_array(k2, [1, n1, 1, n2])
        end if

        ! compute local dips
        !$omp parallel do private(i,j)
        do j = 1, size(w, 2)
            do i = 1, size(w, 1)
                dp(i, j) = -atan(ev1(i, j)%array(2)/(ev1(i, j)%array(1) + 1.0e-10))
                if (present(k1)) then
                    k1(i, j) = ev1(i, j)%array(1)
                end if
                if (present(k2)) then
                    k2(i, j) = ev1(i, j)%array(2)
                end if
            end do
        end do
        !$omp end parallel do

    end subroutine gst_localdip2

    !
    !> Estimate 3D local dip using general structural tensor
    !
    !
    subroutine gst_localdip3(w, dp1, dp2, k1, k2, k3)

        ! arguments
        real, allocatable, dimension(:, :, :), intent(in) :: w
        real, allocatable, dimension(:, :, :), intent(inout) :: dp1, dp2
        real, allocatable, dimension(:, :, :), intent(inout), optional :: k1, k2, k3

        real, allocatable, dimension(:, :, :) :: s11, s12, s13, s22, s23, s33
        type(meta_array1_real), allocatable, dimension(:, :, :) :: ev1, ev2, ev3
        real, allocatable, dimension(:, :, :) :: ea1, ea2, ea3
        integer :: i, j, k
        integer :: n1, n2, n3

        n1 = size(w, 1)
        n2 = size(w, 2)
        n3 = size(w, 3)

        ! compute structure tensor
        call compute_structure_tensor(w, s11, s12, s13, s22, s23, s33)

        ! smooth structure tensor with a sigma=1 Gaussian filter
        s11 = gauss_filt(s11, [1.0, 1.0, 1.0])
        s12 = gauss_filt(s12, [1.0, 1.0, 1.0])
        s13 = gauss_filt(s13, [1.0, 1.0, 1.0])
        s22 = gauss_filt(s22, [1.0, 1.0, 1.0])
        s23 = gauss_filt(s23, [1.0, 1.0, 1.0])
        s33 = gauss_filt(s33, [1.0, 1.0, 1.0])

        ! compute eigenvectors
        call compute_structure_tensor_eigens(s11, s12, s13, s22, s23, s33, ev1, ev2, ev3, ea1, ea2, ea3)

        call alloc_array(dp1, [1, n1, 1, n2, 1, n3])
        call alloc_array(dp2, [1, n1, 1, n2, 1, n3])
        if (present(k1)) then
            call alloc_array(k1, [1, n1, 1, n2, 1, n3])
        end if
        if (present(k2)) then
            call alloc_array(k2, [1, n1, 1, n2, 1, n3])
        end if
        if (present(k3)) then
            call alloc_array(k3, [1, n1, 1, n2, 1, n3])
        end if

        ! compute local dips
        !$omp parallel do private(i,j,k)
        do k = 1, size(w, 3)
            do j = 1, size(w, 2)
                do i = 1, size(w, 1)
                    dp1(i, j, k) = -atan(ev1(i, j, k)%array(3)/(ev1(i, j, k)%array(1) + 1.0e-10))
                    dp2(i, j, k) = -atan(ev1(i, j, k)%array(2)/(ev1(i, j, k)%array(1) + 1.0e-10))
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

    end subroutine gst_localdip3

    function dip_filt_2d_float(w, d, slopes, amps) result(wr)

        real, dimension(:), intent(in) :: slopes, amps
        real, dimension(:), intent(in) :: d
        real, dimension(:, :), intent(in) :: w
        real, allocatable, dimension(:, :) :: wr

        integer :: n1, n2, i, j, nk1, nk2
        complex, allocatable, dimension(:, :) :: ww
        real, allocatable, dimension(:) :: cf1, cf2
        real :: slope(1), amp(1)
        integer :: nsl
        real :: d1, d2

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
                    !call interp_linear(nsl, slopes, amps, 1, slope, amp)
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
        wr = real(ww(1:n1, 1:n2))

    end function dip_filt_2d_float

    function dip_filt_3d_float(w, d, slopes, amps, plane) result(wr)

        real, dimension(:), intent(in) :: slopes, amps
        real, dimension(:), intent(in) :: d
        real, dimension(:, :, :), intent(in) :: w
        integer, intent(in) :: plane
        real, allocatable, dimension(:, :, :) :: wr

        real :: d1, d2, d3
        integer :: i

        d1 = d(1)
        d2 = d(2)
        d3 = d(3)
        wr = w

        select case (plane)
            case (12)
                do i = 1, size(w, 3)
                    wr(:, :, i) = dip_filt_2d_float(wr(:, :, i), [d1, d2], slopes, amps)
                end do
            case (13)
                do i = 1, size(w, 2)
                    wr(:, i, :) = dip_filt_2d_float(wr(:, i, :), [d1, d3], slopes, amps)
                end do
            case (23)
                do i = 1, size(w, 1)
                    wr(i, :, :) = dip_filt_2d_float(wr(i, :, :), [d2, d3], slopes, amps)
                end do
        end select

    end function dip_filt_3d_float

    function dip_filt_2d_complex(w, d, slopes, amps) result(ww)

        real, dimension(:), intent(in) :: slopes, amps
        real, dimension(:), intent(in) :: d
        complex, dimension(:, :), intent(in) :: w
        complex, allocatable, dimension(:, :) :: ww

        integer :: n1, n2, i, j, nk1, nk2
        real, allocatable, dimension(:) :: cf1, cf2
        real :: slope(1), amp(1)
        integer :: nsl
        real :: d1, d2

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
                    !call interp_linear(nsl, slopes, amps, 1, slope, amp)
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
        ww = ww(1:n1, 1:n2)

    end function dip_filt_2d_complex

end module libflit_dip
