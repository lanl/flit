
module iir_filter_mod

    use, intrinsic :: iso_fortran_env, only : real64, int32
    implicit none

    integer, parameter :: dp = real64

    type :: iir_filter
        integer :: n_sections = 0
        ! SOS coefficients, with a shape of (6, n_sections)
        ! Each SOS row is [b0, b1, b2, a0, a1, a2]; usually a0=1.
        real(dp), allocatable :: sos(:, :)
    contains
        procedure :: init              => iir_init
        procedure :: clear             => iir_clear
        procedure :: apply             => iir_apply
        procedure :: apply_zero_phase  => iir_apply_zero_phase
    end type iir_filter

    private
    public :: iir_filter

contains

    subroutine iir_init(this, sos_in)
        class(iir_filter), intent(inout) :: this
        real(dp), intent(in)             :: sos_in(:, :)
        integer :: nsec, j

        call this%clear()
        if (size(sos_in, 1) == 6) then
            nsec = size(sos_in, 2)
            allocate(this%sos(6, nsec))
            this%sos = sos_in
        else if (size(sos_in, 2) == 6) then
            nsec = size(sos_in, 1)
            allocate(this%sos(6, nsec))
            this%sos = transpose(sos_in)
        else
            error stop 'iir_init: SOS must have one dimension of length 6.'
        end if

        this%n_sections = nsec

        do j = 1, nsec
            if (abs(this%sos(4, j)) <= tiny(1.0_dp)) error stop 'iir_init: a0 is zero.'
            ! Normalize each section so a0=1. scipy/obspy SOS rows are normally already normalized.
            if (abs(this%sos(4, j) - 1.0_dp) > 10.0_dp * epsilon(1.0_dp)) then
                this%sos(1:3, j) = this%sos(1:3, j) / this%sos(4, j)
                this%sos(5:6, j) = this%sos(5:6, j) / this%sos(4, j)
                this%sos(4, j)   = 1.0_dp
            end if
        end do
    end subroutine iir_init

    subroutine iir_clear(this)
        class(iir_filter), intent(inout) :: this
        if (allocated(this%sos)) deallocate(this%sos)
        this%n_sections = 0
    end subroutine iir_clear

    subroutine iir_apply(this, x)
        class(iir_filter), intent(in) :: this
        real(dp), intent(inout)       :: x(:)
        call apply_sos(this%sos, x)
    end subroutine iir_apply

    subroutine iir_apply_zero_phase(this, x)
        class(iir_filter), intent(in) :: this
        real(dp), intent(inout)       :: x(:)
        call apply_sos_zero_phase(this%sos, x)
    end subroutine iir_apply_zero_phase

    subroutine apply_sos(sos, x)
        !! In-place SOS filtering with zero initial conditions.
        !! Direct Form II transposed, same difference equation as scipy.signal.sosfilt.
        real(dp), intent(in)    :: sos(:, :)
        real(dp), intent(inout) :: x(:)
        integer :: nsec, j, i, n
        real(dp) :: b0, b1, b2, a1, a2, z1, z2, y, xi

        if (.not. (size(sos, 1) == 6 .or. size(sos, 2) == 6)) &
            error stop 'apply_sos: SOS must have one dimension of length 6.'

        n = size(x)
        if (n == 0) return

        if (size(sos, 1) == 6) then
            nsec = size(sos, 2)
            do j = 1, nsec
                b0 = sos(1, j) / sos(4, j)
                b1 = sos(2, j) / sos(4, j)
                b2 = sos(3, j) / sos(4, j)
                a1 = sos(5, j) / sos(4, j)
                a2 = sos(6, j) / sos(4, j)
                z1 = 0.0_dp; z2 = 0.0_dp
                do i = 1, n
                    xi = x(i)
                    y  = b0 * xi + z1
                    z1 = b1 * xi - a1 * y + z2
                    z2 = b2 * xi - a2 * y
                    x(i) = y
                end do
            end do
        else
            nsec = size(sos, 1)
            do j = 1, nsec
                b0 = sos(j, 1) / sos(j, 4)
                b1 = sos(j, 2) / sos(j, 4)
                b2 = sos(j, 3) / sos(j, 4)
                a1 = sos(j, 5) / sos(j, 4)
                a2 = sos(j, 6) / sos(j, 4)
                z1 = 0.0_dp; z2 = 0.0_dp
                do i = 1, n
                    xi = x(i)
                    y  = b0 * xi + z1
                    z1 = b1 * xi - a1 * y + z2
                    z2 = b2 * xi - a2 * y
                    x(i) = y
                end do
            end do
        end if
    end subroutine apply_sos

    subroutine apply_sos_zero_phase(sos, x)
        !! ObsPy-style zero-phase IIR: sosfilt forward, reverse, sosfilt, reverse.
        real(dp), intent(in)    :: sos(:, :)
        real(dp), intent(inout) :: x(:)
        call apply_sos(sos, x)
        call reverse_inplace(x)
        call apply_sos(sos, x)
        call reverse_inplace(x)
    end subroutine apply_sos_zero_phase

    subroutine reverse_inplace(x)
        real(dp), intent(inout) :: x(:)
        integer :: i, j
        real(dp) :: tmp
        i = 1; j = size(x)
        do while (i < j)
            tmp = x(i); x(i) = x(j); x(j) = tmp
            i = i + 1; j = j - 1
        end do
    end subroutine reverse_inplace

end module iir_filter_mod
