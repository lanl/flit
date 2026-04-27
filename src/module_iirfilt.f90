!
! © 2024-2026. Triad National Security, LLC. All rights reserved.
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


module libflit_iirfilt

    use libflit_string
    use libflit_utility

    ! Butterworth IIR filtering
    use iir_design_butterworth_mod, only: design_butterworth_sos
    use iir_filter_mod, only: iir_filter

    interface iir_filt
        module procedure :: iir_filt_float
        module procedure :: iir_filt_double
    end interface iir_filt

    private
    public :: iir_filt

contains

    function iir_filt_float(x, dt, filter_type, flow, fhigh, order, zerophase) result(xx)

        real, dimension(:), intent(in) :: x
        real, intent(in) :: dt
        character(len=*), intent(in) :: filter_type
        real, intent(in) :: flow, fhigh
        integer, intent(in), optional :: order
        logical, intent(in), optional :: zerophase
        real, allocatable, dimension(:) :: xx

        type(iir_filter) :: filt
        double precision, allocatable, dimension(:, :) :: sos
        double precision, allocatable, dimension(:) :: data
        integer :: iir_order
        logical :: iir_zerophase

        if (present(order)) then
            iir_order = clip(order, 1, 10)
        else
            iir_order = 4
        end if

        if (present(zerophase)) then
            iir_zerophase = zerophase
        else
            iir_zerophase = .true.
        end if

        call design_butterworth_sos(tidy(filter_type), iir_order, dble(1.0/dt), dble(flow), dble(fhigh), sos)

        call filt%init(sos)

        data = x
        if (iir_zerophase) then
            call filt%apply_zero_phase(data)
        else
            call filt%apply(data)
        end if
        xx = data

    end function iir_filt_float


    function iir_filt_double(x, dt, filter_type, flow, fhigh, order, zerophase) result(xx)

        double precision, dimension(:), intent(in) :: x
        double precision, intent(in) :: dt
        character(len=*), intent(in) :: filter_type
        double precision, intent(in) :: flow, fhigh
        integer, intent(in), optional :: order
        logical, intent(in), optional :: zerophase
        double precision, allocatable, dimension(:) :: xx

        type(iir_filter) :: filt
        double precision, allocatable, dimension(:, :) :: sos
        double precision, allocatable, dimension(:) :: data
        integer :: iir_order
        logical :: iir_zerophase

        if (present(order)) then
            iir_order = clip(order, 1, 10)
        else
            iir_order = 4
        end if

        if (present(zerophase)) then
            iir_zerophase = zerophase
        else
            iir_zerophase = .true.
        end if

        call design_butterworth_sos(tidy(filter_type), iir_order, 1.0d0/dt, flow, fhigh, sos)

        call filt%init(sos)

        data = x
        if (iir_zerophase) then
            call filt%apply_zero_phase(data)
        else
            call filt%apply(data)
        end if
        xx = data

    end function iir_filt_double

end module libflit_iirfilt
