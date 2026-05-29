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


module libflit_geometry

    use libflit_linear_algebra
    use libflit_constants
    use libflit_utility
    use libflit_array
    use libflit_calculus
    use libflit_error
    use libflit_fit
    use libflit_lowessfilt
    use libflit_statistics
    use libflit_interp
    use libflit_array_operation
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32

    implicit none

    interface rotation_matrix
        module procedure :: rotation_matrix_along_arbitrary_float
        module procedure :: rotation_matrix_along_axis_float
        module procedure :: rotation_matrix_2d_float
        module procedure :: rotation_matrix_3d_spherical_float
        module procedure :: rotation_matrix_3d_euler_float
        module procedure :: rotation_matrix_from_to_float
        module procedure :: rotation_matrix_along_arbitrary_double
        module procedure :: rotation_matrix_along_axis_double
        module procedure :: rotation_matrix_2d_double
        module procedure :: rotation_matrix_3d_spherical_double
        module procedure :: rotation_matrix_3d_euler_double
        module procedure :: rotation_matrix_from_to_double
    end interface

    interface rotate_point
        module procedure :: rotate_point_2d_float
        module procedure :: rotate_points_2d_float
        module procedure :: rotate_point_3d_float
        module procedure :: rotate_points_3d_float
        module procedure :: rotate_point_2d_double
        module procedure :: rotate_points_2d_double
        module procedure :: rotate_point_3d_double
        module procedure :: rotate_points_3d_double
    end interface

    interface convexhull
        module procedure :: convexhull_2d_float
        module procedure :: convexhull_3d_float
        module procedure :: convexhull_2d_double
        module procedure :: convexhull_3d_double
    end interface

    interface polygon_contains_point
        module procedure :: polygon_contains_point_float
        module procedure :: polygon_contains_point_double
    end interface

    interface distance_between_line_segments
        module procedure :: distance_between_line_segments_float
        module procedure :: distance_between_line_segments_double
    end interface

    interface point_distance_to_line_segement
        module procedure :: point_distance_to_line_segement_float
        module procedure :: point_distance_to_line_segement_double
    end interface

    interface point_distance_to_line
        module procedure :: point_distance_to_line_float
        module procedure :: point_distance_to_line_double
    end interface

    interface point_signed_distance_to_line
        module procedure :: point_signed_distance_to_line_float
        module procedure :: point_signed_distance_to_line_double
    end interface

    interface gaussian_curvature
        module procedure :: gaussian_curvature_float
        module procedure :: gaussian_curvature_double
    end interface

    interface spherical_to_cartesian
        module procedure :: spherical_to_cartesian_float
        module procedure :: spherical_to_cartesian_double
    end interface

    interface cartesian_to_spherical
        module procedure :: cartesian_to_spherical_float
        module procedure :: cartesian_to_spherical_double
    end interface

    interface fit_curve
        module procedure :: fit_curve_float
        module procedure :: fit_curve_double
    end interface

    interface fit_surface
        module procedure :: fit_surface_float
        module procedure :: fit_surface_double
    end interface

    private

    public :: convexhull
    public :: polygon_contains_point
    public :: rotate_point
    public :: rotation_matrix
    public :: distance_between_line_segments
    public :: point_distance_to_line_segement
    public :: point_distance_to_line
    public :: point_signed_distance_to_line
    public :: gaussian_curvature
    public :: spherical_to_cartesian
    public :: cartesian_to_spherical
    public :: fit_curve
    public :: fit_surface

contains

#define T float
#define TT real
#define fp sp
#include "template_rotate.f90"
#undef fp

#define T double
#define TT double precision
#define fp dp
#include "template_rotate.f90"
#undef fp

#define T float
#define TT real
#define TTT real
#define fp sp
#include "template_geometry.f90"
#undef fp

#define T double
#define TT double precision
#define TTT real
#define fp dp
#include "template_geometry.f90"
#undef fp

end module libflit_geometry
