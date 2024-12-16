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


module libflit_interp

    use libflit_specialfunc
    use libflit_constants
    use libflit_array
    use libflit_array_operation
    use libflit_unique
    use libflit_linear_algebra
    use libflit_error
    use libflit_string
    use libflit_utility
    use libflit_calculus
    use libflit_taper

    implicit none

    ! External interfaces
    interface

        subroutine interp_cspline_1d_float(n, x, y, nn, xx, yy) bind(c, name='cspline1d')

            use iso_c_binding, only: c_int, c_float

            integer(kind=c_int), value :: n, nn

            real, dimension(*), intent(in) :: x, y, xx
            real, dimension(*), intent(out) :: yy

        end subroutine interp_cspline_1d_float

        subroutine interp_pchip_1d_float(n, x, y, nn, xx, yy) bind(c, name='pchip1d')

            use iso_c_binding, only: c_int, c_float

            integer(kind=c_int), value :: n, nn

            real, dimension(*), intent(in) :: x, y, xx
            real, dimension(*), intent(out) :: yy

        end subroutine interp_pchip_1d_float

        subroutine interp_quintic_1d_float(n, x, y, nn, xx, yy) bind(c, name='quintic1d')

            use iso_c_binding, only: c_int, c_float

            integer(kind=c_int), value :: n, nn

            real, dimension(*), intent(in) :: x, y, xx
            real, dimension(*), intent(out) :: yy

        end subroutine interp_quintic_1d_float

        subroutine interp_cspline_1d_double(n, x, y, nn, xx, yy) bind(c, name='cspline1d_double')

            use iso_c_binding, only: c_int, c_double

            integer(kind=c_int), value :: n, nn

            double precision, dimension(*), intent(in) :: x, y, xx
            double precision, dimension(*), intent(out) :: yy

        end subroutine interp_cspline_1d_double

        subroutine interp_pchip_1d_double(n, x, y, nn, xx, yy) bind(c, name='pchip1d_double')

            use iso_c_binding, only: c_int, c_double

            integer(kind=c_int), value :: n, nn

            double precision, dimension(*), intent(in) :: x, y, xx
            double precision, dimension(*), intent(out) :: yy

        end subroutine interp_pchip_1d_double

        subroutine interp_quintic_1d_double(n, x, y, nn, xx, yy) bind(c, name='quintic1d_double')

            use iso_c_binding, only: c_int, c_double

            integer(kind=c_int), value :: n, nn

            double precision, dimension(*), intent(in) :: x, y, xx
            double precision, dimension(*), intent(out) :: yy

        end subroutine interp_quintic_1d_double

        subroutine interp_mba_1d_float(n, x, y, nn, xx, yy) bind(c, name='mba1')

            use iso_c_binding, only: c_int, c_float

            integer(kind=c_int), value :: n, nn
            real(kind=c_float), dimension(*), intent(in) :: x, y, xx
            real(kind=c_float), dimension(*), intent(out) :: yy

        end subroutine interp_mba_1d_float

        subroutine interp_mba_2d_float(n, x, y, z, nn, xx, yy, zz) bind(c, name='mba2')

            use iso_c_binding, only: c_int, c_float

            integer(kind=c_int), value :: n, nn
            real(kind=c_float), dimension(*), intent(in) :: x, y, z, xx, yy
            real(kind=c_float), dimension(*), intent(out) :: zz

        end subroutine interp_mba_2d_float

        subroutine interp_mba_3d_float(n, x, y, z, v, nn, xx, yy, zz, vv) bind(c, name='mba3')

            use iso_c_binding, only: c_int, c_float

            integer(kind=c_int), value :: n, nn
            real(kind=c_float), dimension(*), intent(in) :: x, y, z, v, xx, yy, zz
            real(kind=c_float), dimension(*), intent(out) :: vv

        end subroutine interp_mba_3d_float

        subroutine interp_mba_1d_double(n, x, y, nn, xx, yy) bind(c, name='mba1_double')

            use iso_c_binding, only: c_int, c_double

            integer(kind=c_int), value :: n, nn
            real(kind=c_double), dimension(*), intent(in) :: x, y, xx
            real(kind=c_double), dimension(*), intent(out) :: yy

        end subroutine interp_mba_1d_double

        subroutine interp_mba_2d_double(n, x, y, z, nn, xx, yy, zz) bind(c, name='mba2_double')

            use iso_c_binding, only: c_int, c_double

            integer(kind=c_int), value :: n, nn
            real(kind=c_double), dimension(*), intent(in) :: x, y, z, xx, yy
            real(kind=c_double), dimension(*), intent(out) :: zz

        end subroutine interp_mba_2d_double

        subroutine interp_mba_3d_double(n, x, y, z, v, nn, xx, yy, zz, vv) bind(c, name='mba3_double')

            use iso_c_binding, only: c_int, c_double

            integer(kind=c_int), value :: n, nn
            real(kind=c_double), dimension(*), intent(in) :: x, y, z, v, xx, yy, zz
            real(kind=c_double), dimension(*), intent(out) :: vv

        end subroutine interp_mba_3d_double

    end interface

    ! Internal interfaces
    interface interp
        module procedure :: reg_to_reg_interp_1d_float
        module procedure :: reg_to_reg_interp_1d_double
        module procedure :: reg_to_reg_interp_2d_float
        module procedure :: reg_to_reg_interp_2d_double
        module procedure :: reg_to_reg_interp_3d_float
        module procedure :: reg_to_reg_interp_3d_double
    end interface interp

    interface resample
        module procedure :: resample_1d_float
        module procedure :: resample_1d_double
        module procedure :: resample_2d_float
        module procedure :: resample_2d_double
        module procedure :: resample_3d_float
        module procedure :: resample_3d_double
    end interface resample

    interface interp_to
        module procedure :: interp_to_1d_float
        module procedure :: interp_to_1d_double
        module procedure :: interp_to_2d_float
        module procedure :: interp_to_2d_double
        module procedure :: interp_to_3d_float
        module procedure :: interp_to_3d_double
    end interface interp_to

    interface interp_like
        module procedure :: interp_like_1d_float
        module procedure :: interp_like_1d_double
        module procedure :: interp_like_2d_float
        module procedure :: interp_like_2d_double
        module procedure :: interp_like_3d_float
        module procedure :: interp_like_3d_double
    end interface interp_like

    interface ginterp
        module procedure :: irreg_to_irreg_interp_1d_float
        module procedure :: irreg_to_irreg_interp_1d_double
        module procedure :: irreg_to_irreg_interp_2d_float
        module procedure :: irreg_to_irreg_interp_2d_double
        module procedure :: irreg_to_irreg_interp_3d_float
        module procedure :: irreg_to_irreg_interp_3d_double
        module procedure :: irreg_to_reg_interp_1d_float
        module procedure :: irreg_to_reg_interp_1d_double
        module procedure :: irreg_to_reg_interp_2d_float
        module procedure :: irreg_to_reg_interp_2d_double
        module procedure :: irreg_to_reg_interp_3d_float
        module procedure :: irreg_to_reg_interp_3d_double
    end interface ginterp

    interface inpaint
        module procedure :: inpaint_1d_float
        module procedure :: inpaint_1d_double
        module procedure :: inpaint_2d_float
        module procedure :: inpaint_2d_double
        module procedure :: inpaint_3d_float
        module procedure :: inpaint_3d_double
    end interface inpaint

    interface meshgrid
        module procedure :: meshgrid_float
        module procedure :: meshgrid_double
        module procedure :: meshgrid_1d_float
        module procedure :: meshgrid_1d_double
        module procedure :: meshgrid_2d_float
        module procedure :: meshgrid_2d_double
        module procedure :: meshgrid_3d_float
        module procedure :: meshgrid_3d_double
    end interface meshgrid

    interface point_interp_linear
        module procedure :: point_interp_linear_1d_float
        module procedure :: point_interp_linear_2d_float
        module procedure :: point_interp_linear_3d_float
        module procedure :: point_interp_linear_1d_double
        module procedure :: point_interp_linear_2d_double
        module procedure :: point_interp_linear_3d_double
    end interface point_interp_linear

    interface point_interp_barycentric
        module procedure :: point_interp_barycentric_2d_float
        module procedure :: point_interp_barycentric_3d_float
        module procedure :: point_interp_barycentric_2d_double
        module procedure :: point_interp_barycentric_3d_double
    end interface point_interp_barycentric

    private

    public :: interp
    public :: resample
    public :: interp_to
    public :: interp_like
    public :: ginterp
    public :: meshgrid
    public :: inpaint
    public :: point_interp_linear
    public :: point_interp_barycentric

contains

#define T float
#define TT real
#define TTT real
#define nTT real
#define nTTT real
#include "template_interp.f90"

#define T double
#define TT double precision
#define TTT double precision
#define nTT dble
#define nTTT dble
#include "template_interp.f90"

end module libflit_interp
