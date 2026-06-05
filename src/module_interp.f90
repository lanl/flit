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

    logical :: mba_verbose = .false.

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

    public :: mba_verbose

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
