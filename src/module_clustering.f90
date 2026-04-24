!
! © 2024--2026. Triad National Security, LLC. All rights reserved.
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


module libflit_clustering

    use hdbscan_mod
    use libflit_error

    implicit none

    type hdbscan

        real, allocatable :: data(:, :)

        integer :: min_cluster_size = 5
        integer :: min_samples = 5
        real :: cluster_selection_epsilon = 0.0
        logical :: allow_single_cluster = .false.
        character(len=16) :: metric = 'euclidean'
        real :: minkowski_p = 2.0

        integer, allocatable :: labels(:)
        real, allocatable :: probabilities(:)
        real, allocatable :: outlier_scores(:)
        real, allocatable :: core_distances(:)
        integer :: n_clusters = 0
        integer :: n_noise = 0

        integer :: n_points = 0

    contains
        procedure :: cluster => hdbscan_clustering

    end type

    private
    public :: hdbscan

contains

    subroutine hdbscan_clustering(this)

        class(hdbscan), intent(inout) :: this

        type(hdbscan_params) :: params
        type(hdbscan_result) :: res

        ! Check data
        call assert(allocated(this%data), '<hdbscan_clustering> Error: data is not initialized. ')

        ! Copy parameter
        params%min_cluster_size = this%min_cluster_size
        params%min_samples = this%min_samples
        params%metric = this%metric
        params%cluster_selection_epsilon = this%cluster_selection_epsilon
        params%allow_single_cluster = this%allow_single_cluster

        ! Run HDBSCAN clustering
        call hdbscan_fit(transpose(dble(this%data)), params, res)

        ! Return results
        this%n_clusters = res%n_clusters
        this%n_noise = res%n_noise
        this%labels = res%labels
        this%probabilities = res%probabilities
        this%outlier_scores = res%outlier_scores

        this%n_points = size(this%data, 1)

        ! Free memory
        call hdbscan_free(res)

    end subroutine

end module
