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


module libflit_clustering

    use libflit_array
    use libflit_array_operation
    use libflit_utility
    use libflit_string
    use libflit_error

    implicit none

    interface

        !> Interface to C++ HDBSCAN to support the subroutines in this module,
        !> but could also be used as standalone subroutine for data of arbitrary dimension
        subroutine hdbscan_clustering(n, nd, x, min_sample, min_cluster_size, metric, p, ncluster, nnoisy, labels, probs) bind(c, name='hdbscan')

            use iso_c_binding, only: c_float, c_int

            !> Number of points
            integer(kind=c_int), value, intent(in) :: n
            !> Dimension of a point, nd >= 1
            integer(kind=c_int), value, intent(in) :: nd
            !> Input data organized as (x_1, y_1, ...), (x_2, y_2, ...), ..., (x_n, y_n, ...)
            real(kind=c_float), dimension(*), intent(in) :: x
            !> The number of samples in a neighborhood for a point to be considered as a core point, including the point itself.
            integer(kind=c_int), value, intent(in) :: min_sample
            !> The minimum number of samples in a group for that group to be considered a cluster;
            !> groupings smaller than this size will be left as noise.
            integer(kind=c_int), value, intent(in) :: min_cluster_size
            !> Distance metric, could be manhanttan (l1), euclidean (l2), and minkowski (l_p, with p >= 0)
            !> Minkowski distance is a generalization of Manhattan and Euclidean
            integer(kind=c_int), value, intent(in) :: metric
            !> If metric = minkowski, p is the power of the norm
            real(kind=c_float), value, intent(in) :: p
            !> Number of clusters
            integer(kind=c_int), intent(out) :: ncluster
            !> Number of noisy points (ungrouped)
            integer(kind=c_int), intent(out) :: nnoisy
            !> Labels of length n indicating which group a point belongs to; = 0 means a noisy point.
            integer(kind=c_int), dimension(*), intent(out) :: labels
            !> Probabilities of length n indicating the probability of a point belonging to the
            !> labeled group.
            real(kind=c_float), dimension(*), intent(out) :: probs

        end subroutine hdbscan_clustering

    end interface

    type hdbscan_param
        integer :: n
        integer :: nd
        real, allocatable, dimension(:, :) :: data
        integer :: min_sample = 5
        integer :: min_cluster_size = 5
        character(len=24) :: metric = 'euclidean'
        real :: p = 3.0
        integer :: ncluster
        integer :: nnoisy
        integer, allocatable, dimension(:) :: labels
        real, allocatable, dimension(:) :: probs
    end type hdbscan_param

    private
    public :: hdbscan_param
    public :: hdbscan
    public :: hdbscan_clustering

contains

    subroutine hdbscan(this)

        type(hdbscan_param), intent(inout) :: this

        call assert(this%n >= 2, ' <hdbscan> Error: Number of points must >= 2.')
        call assert(this%nd >= 1, ' <hdbscan> Error: Dimension of a point must >= 1.')
        call assert(allocated(this%data), ' <hdbscan> Error: Input data are not defined.')
        call assert(size(this%data, 1) == this%n .and. size(this%data, 2) == this%nd, &
            ' <hdbscan> Error: Dimensions of the data are incorrect.')

        this%labels = zeros(this%n)
        this%probs = zeros(this%n)

        select case (this%metric)
            case default
                call hdbscan_clustering(this%n, this%nd, flatten(transpose(this%data)), this%min_sample, this%min_cluster_size, &
                    1, this%p, this%ncluster, this%nnoisy, this%labels, this%probs)
            case ('euclidean')
                call hdbscan_clustering(this%n, this%nd, flatten(transpose(this%data)), this%min_sample, this%min_cluster_size, &
                    1, this%p, this%ncluster, this%nnoisy, this%labels, this%probs)
            case ('manhattan')
                call hdbscan_clustering(this%n, this%nd, flatten(transpose(this%data)), this%min_sample, this%min_cluster_size, &
                    2, this%p, this%ncluster, this%nnoisy, this%labels, this%probs)
            case ('minkowski')
                call hdbscan_clustering(this%n, this%nd, flatten(transpose(this%data)), this%min_sample, this%min_cluster_size, &
                    3, this%p, this%ncluster, this%nnoisy, this%labels, this%probs)
        end select

    end subroutine hdbscan

end module libflit_clustering
