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


module libflit_tvfilt

    use libflit_array
    use libflit_string
    use libflit_date_time
    use libflit_mpicomm
    use libflit_utility
    use libflit_constants
    use libflit_array_operation

    implicit none

    private

    !
    !> Isotropic TV filtering, implementing the algorithm in
    !>
    !>      Goldstein, Osher, 2009, The Split Bregman Method for L1-Regularized Problems,
    !>      SIAM Journal of Imaging Sciences, doi: 10.1137/080725891
    !
    interface tv_filt
        module procedure :: tv_iso_filt_1d_float
        module procedure :: tv_iso_filt_2d_float
        module procedure :: tv_iso_filt_3d_float
        module procedure :: tv_iso_filt_1d_double
        module procedure :: tv_iso_filt_2d_double
        module procedure :: tv_iso_filt_3d_double
    end interface tv_filt

    !
    !> Total generalized p-variation filtering, implementing the algorithm in
    !>
    !>      Knoll et al., 2011, Second order total generalized variation (TGV) for MRI,
    !>      DOI 10.1002/mrm.22595
    !>
    !>      Gao, Huang, 2019, Acoustic- and elastic-waveform inversion with
    !>      total generalized p-variation regularization
    !>      doi: 10.1093/gji/ggz203
    !
    interface tgpv_filt
        module procedure :: tgpv_filt_1d_float
        module procedure :: tgpv_filt_2d_float
        module procedure :: tgpv_filt_3d_float
        module procedure :: tgpv_filt_1d_double
        module procedure :: tgpv_filt_2d_double
        module procedure :: tgpv_filt_3d_double
    end interface tgpv_filt

    !
    !> MPI version of 2D and 3D TGpV filtering
    !
    interface tgpv_filt_mpi
        module procedure :: tgpv_filt_2d_mpi_float
        module procedure :: tgpv_filt_3d_mpi_float
        module procedure :: tgpv_filt_2d_mpi_double
        module procedure :: tgpv_filt_3d_mpi_double
    end interface tgpv_filt_mpi

    !
    !> Isotropic TV + sparsity filtering, implementing the algorithm in
    !>
    !>      Gao et al., 2022, SREMI: Super-resolution electromagnetic imaging with
    !>      single-channel ground-penetrating radar, Journal of Applied Geophysics,
    !>      doi: 10.1016/j.jappgeo.2022.104777
    !
    interface sparse_tv_filt
        module procedure :: sparse_tv_filt_1d_float
        module procedure :: sparse_tv_filt_2d_float
        module procedure :: sparse_tv_filt_3d_float
        module procedure :: sparse_tv_filt_1d_double
        module procedure :: sparse_tv_filt_2d_double
        module procedure :: sparse_tv_filt_3d_double
    end interface

    !
    !> Soft shrinkage, implementing Goldstein, Osher (2009)
    !
    interface soft_shrinkage
        module procedure :: soft_shrinkage_1d_float
        module procedure :: soft_shrinkage_2d_float
        module procedure :: soft_shrinkage_3d_float
        module procedure :: soft_shrinkage_1d_double
        module procedure :: soft_shrinkage_2d_double
        module procedure :: soft_shrinkage_3d_double
    end interface soft_shrinkage

    public :: tv_filt
    public :: sparse_tv_filt
    public :: tgpv_filt
    public :: tgpv_filt_mpi
    public :: soft_shrinkage

contains

#define T float
#define TT real
#include "template_tvfilt.f90"

#define T double
#define TT double precision
#include "template_tvfilt.f90"

end module libflit_tvfilt
