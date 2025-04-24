!
! © 2025. Triad National Security, LLC. All rights reserved.
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


module libflit_mpicomm_group

    use mpi
    ! using mpi_f08 seems to cause mysterious errors in mpi_sendrecv
    use libflit_array
    use libflit_domain_decomposition
    use libflit_error
    use libflit_array_operation
    use libflit_mpicomm
    use iso_fortran_env

    implicit none

    interface bcast_array_group
        module procedure :: bcast_array_1d_group_string
        module procedure :: bcast_array_1d_group_int
        module procedure :: bcast_array_2d_group_int
        module procedure :: bcast_array_3d_group_int
        module procedure :: bcast_array_1d_group_double
        module procedure :: bcast_array_2d_group_double
        module procedure :: bcast_array_3d_group_double
        module procedure :: bcast_array_1d_group_float
        module procedure :: bcast_array_2d_group_float
        module procedure :: bcast_array_3d_group_float
        module procedure :: bcast_array_1d_group_complex
        module procedure :: bcast_array_2d_group_complex
        module procedure :: bcast_array_3d_group_complex
    end interface

    interface commute_array_group
        module procedure :: commute_array_1d_group_int
        module procedure :: commute_array_2d_group_int
        module procedure :: commute_array_3d_group_int
        module procedure :: commute_array_1d_group_double
        module procedure :: commute_array_2d_group_double
        module procedure :: commute_array_3d_group_double
        module procedure :: commute_array_1d_group_float
        module procedure :: commute_array_2d_group_float
        module procedure :: commute_array_3d_group_float
        module procedure :: commute_array_1d_group_complex
        module procedure :: commute_array_2d_group_complex
        module procedure :: commute_array_3d_group_complex
    end interface

    interface allreduce_group
        ! Scalar
        module procedure :: gather_distribute_group_int
        module procedure :: gather_distribute_group_double
        module procedure :: gather_distribute_group_float
        module procedure :: gather_distribute_group_complex
        module procedure :: gather_distribute_group_dcomplex
    end interface

    interface allreduce_array_group
        ! Array
        module procedure :: gather_distribute_array_1d_group_int
        module procedure :: gather_distribute_array_2d_group_int
        module procedure :: gather_distribute_array_3d_group_int
        module procedure :: gather_distribute_array_1d_group_double
        module procedure :: gather_distribute_array_2d_group_double
        module procedure :: gather_distribute_array_3d_group_double
        module procedure :: gather_distribute_array_1d_group_float
        module procedure :: gather_distribute_array_2d_group_float
        module procedure :: gather_distribute_array_3d_group_float
        module procedure :: gather_distribute_array_1d_group_complex
        module procedure :: gather_distribute_array_2d_group_complex
        module procedure :: gather_distribute_array_3d_group_complex
        module procedure :: gather_distribute_array_1d_group_dcomplex
        module procedure :: gather_distribute_array_2d_group_dcomplex
        module procedure :: gather_distribute_array_3d_group_dcomplex
        ! For array with size = int8
        module procedure :: gather_distribute_large_array_1d_group_int
        module procedure :: gather_distribute_large_array_1d_group_double
        module procedure :: gather_distribute_large_array_1d_group_float
        module procedure :: gather_distribute_large_array_1d_group_complex
        module procedure :: gather_distribute_large_array_1d_group_dcomplex
    end interface

    interface reduce_group
        module procedure :: gather_group_int
        module procedure :: gather_group_double
        module procedure :: gather_group_float
        module procedure :: gather_group_complex
        module procedure :: gather_group_dcomplex
    end interface

    interface reduce_array_group
        module procedure :: gather_array_1d_group_int
        module procedure :: gather_array_2d_group_int
        module procedure :: gather_array_3d_group_int
        module procedure :: gather_array_1d_group_double
        module procedure :: gather_array_2d_group_double
        module procedure :: gather_array_3d_group_double
        module procedure :: gather_array_1d_group_float
        module procedure :: gather_array_2d_group_float
        module procedure :: gather_array_3d_group_float
        module procedure :: gather_array_1d_group_complex
        module procedure :: gather_array_2d_group_complex
        module procedure :: gather_array_3d_group_complex
        module procedure :: gather_array_1d_group_dcomplex
        module procedure :: gather_array_2d_group_dcomplex
        module procedure :: gather_array_3d_group_dcomplex
    end interface

    interface domain_decomp_regular_group
        module procedure :: domain_decomp_regular_2d
        module procedure :: domain_decomp_regular_3d
    end interface

    !    interface mpi_group_min
    !        module procedure :: mpi_group_min_1d_int
    !        module procedure :: mpi_group_min_1d_float
    !        module procedure :: mpi_group_min_1d_double
    !    end interface mpi_group_min

    integer, public :: ngroup = 1
    integer, public :: groupid
    integer, public :: mpi_group_comm

    integer, public :: rank1_group, rank2_group, rank3_group
    integer, public :: block_x1left_group, block_x1right_group
    integer, public :: block_x2left_group, block_x2right_group
    integer, public :: block_x3left_group, block_x3right_group
    integer, public :: rankid_group
    integer, public :: blockid_group
    integer, public :: nrank_group
    type(mpi_status), public :: mpi_stats_group
    integer, public :: mpi_ierr_group

    public :: mpistart_group, mpiend_group, mpibarrier_group, mpistop_group
    public :: bcast_array_group
    public :: reduce_group
    public :: reduce_array_group
    public :: allreduce_group
    public :: allreduce_array_group
    public :: commute_array_group
    public :: domain_decomp_regular_group

contains

    subroutine mpistart_group

        integer, allocatable, dimension(:, :) :: r
        integer :: i

        call assert(ngroup >= 1 .and. ngroup <= nrank, ' <mpistart_group> Error: ngroup must >= 1 and <= nrank')

        call alloc_array(r, [0, ngroup - 1, 1, 2])
        call cut(0, nrank - 1, ngroup, r)

        do i = 0, ngroup - 1
            if (rankid >= r(i, 1) .and. rankid <= r(i, 2)) then
                groupid = i
                exit
            end if
        end do
        !        groupid = rankid/nrank_group
        call mpi_comm_split(mpi_comm_world, groupid, nrank, mpi_group_comm, mpi_ierr)

        call mpi_comm_size(mpi_group_comm, nrank_group, mpi_ierr_group)
        call mpi_comm_rank(mpi_group_comm, rankid_group, mpi_ierr_group)

        !        write(error_unit, '(a, i6, a, i6, a, i6)') ' @ group ', groupid, ' with group id ', rankid_group, '/', nrank_group

    end subroutine mpistart_group

    subroutine mpiend_group

        call mpi_barrier(mpi_group_comm, mpi_ierr_group)
        call mpi_comm_free(mpi_group_comm, mpi_ierr_group)

    end subroutine mpiend_group

    subroutine mpibarrier_group

        call mpi_barrier(mpi_group_comm, mpi_ierr_group)

    end subroutine mpibarrier_group

    subroutine mpistop_group

        call mpi_abort(mpi_group_comm, mpi_err_other, mpi_ierr_group)

    end subroutine mpistop_group

    !
    !> Decompose a matrix into MPI blocks
    !
    subroutine domain_decomp_regular_2d_group(n1, n2, n1beg, n1end, n2beg, n2end)

        integer, intent(in) :: n1, n2
        integer, intent(out) :: n1beg, n1end, n2beg, n2end

        integer, allocatable, dimension(:, :) :: blk1, blk2
        integer :: i, j

        call mpi_comm_rank(mpi_group_comm, rankid_group, mpi_ierr_group)

        call cut(1, n1, rank1_group, blk1)
        call cut(1, n2, rank2_group, blk2)

        if (rankid_group > rank1_group*rank2_group - 1) then
            n1beg = 1
            n1end = -1
            n2beg = 1
            n2end = -1
        else
            do j = 0, rank2_group - 1
                do i = 0, rank1_group - 1
                    if (rankid_group == j*rank1_group + i) then
                        n1beg = blk1(i + 1, 1)
                        n1end = blk1(i + 1, 2)
                        n2beg = blk2(j + 1, 1)
                        n2end = blk2(j + 1, 2)
                    end if
                end do
            end do

            ! Find adjacent blocks for each group
            do j = 0, rank2_group - 1
                do i = 0, rank1_group - 1

                    blockid = j*rank1_group + i

                    if (blockid == rankid_group) then

                        block_x1left_group = rankid_group - 1
                        block_x1right_group = rankid_group + 1
                        block_x2left_group = rankid_group - rank1_group
                        block_x2right_group = rankid_group + rank1_group

                        if (i == 0) then
                            block_x1left_group = mpi_proc_null
                        end if
                        if (i == rank1_group - 1) then
                            block_x1right_group = mpi_proc_null
                        end if
                        if (j == 0) then
                            block_x2left_group = mpi_proc_null
                        end if
                        if (j == rank2_group - 1) then
                            block_x2right_group = mpi_proc_null
                        end if

                    end if

                end do
            end do

        end if

    end subroutine domain_decomp_regular_2d_group

    !
    !> Decompose a volume into MPI blocks
    !
    subroutine domain_decomp_regular_3d_group(n1, n2, n3, n1beg, n1end, n2beg, n2end, n3beg, n3end)

        integer, intent(in) :: n1, n2, n3
        integer, intent(out) :: n1beg, n1end, n2beg, n2end, n3beg, n3end

        integer, allocatable, dimension(:, :) :: blk1, blk2, blk3
        integer :: i, j, k

        call mpi_comm_rank(mpi_group_comm, rankid_group, mpi_ierr_group)

        call cut(1, n1, rank1_group, blk1)
        call cut(1, n2, rank2_group, blk2)
        call cut(1, n3, rank3_group, blk3)

        if (rankid_group > rank1_group*rank2_group*rank3_group - 1) then
            n1beg = 1
            n1end = -1
            n2beg = 1
            n2end = -1
            n3beg = 1
            n3end = -1
        else
            do k = 0, rank3_group - 1
                do j = 0, rank2_group - 1
                    do i = 0, rank1_group - 1
                        if (rankid_group == k*rank1_group*rank2_group + j*rank1_group + i) then
                            n1beg = blk1(i + 1, 1)
                            n1end = blk1(i + 1, 2)
                            n2beg = blk2(j + 1, 1)
                            n2end = blk2(j + 1, 2)
                            n3beg = blk3(k + 1, 1)
                            n3end = blk3(k + 1, 2)
                        end if
                    end do
                end do
            end do

            ! Find adjacent blocks for each group
            do k = 0, rank3_group - 1
                do j = 0, rank2_group - 1
                    do i = 0, rank1_group - 1

                        blockid = k*rank2_group*rank1_group + j*rank1_group + i

                        if (blockid == rankid_group) then

                            block_x1left_group = rankid_group - 1
                            block_x1right_group = rankid_group + 1
                            block_x2left_group = rankid_group - rank1_group
                            block_x2right_group = rankid_group + rank1_group
                            block_x3left_group = rankid_group - rank1_group*rank2_group
                            block_x3right_group = rankid_group + rank1_group*rank2_group

                            if (i == 0) then
                                block_x1left_group = mpi_proc_null
                            end if
                            if (i == rank1_group - 1) then
                                block_x1right_group = mpi_proc_null
                            end if
                            if (j == 0) then
                                block_x2left_group = mpi_proc_null
                            end if
                            if (j == rank2_group - 1) then
                                block_x2right_group = mpi_proc_null
                            end if
                            if (k == 0) then
                                block_x3left_group = mpi_proc_null
                            end if
                            if (k == rank3_group - 1) then
                                block_x3right_group = mpi_proc_null
                            end if

                        end if

                    end do
                end do
            end do

        end if

    end subroutine domain_decomp_regular_3d_group

    ! MPI array operations
#define T int
#define TT integer
#define TTT mpi_integer
#include "template_mpicomm_group.f90"

#define T float
#define TT real
#define TTT mpi_real
#include "template_mpicomm_group.f90"

#define T double
#define TT double precision
#define TTT mpi_double_precision
#include "template_mpicomm_group.f90"

#define T complex
#define TT complex
#define TTT mpi_complex
#include "template_mpicomm_group.f90"

#define T dcomplex
#define TT double complex
#define TTT mpi_double_complex
#include "template_mpicomm_group.f90"

    ! For string
    subroutine bcast_array_1d_group_string(w, rankid)

        character(len=*), intent(inout) :: w
        integer, intent(in), optional :: rankid

        integer :: rid

        if (present(rankid)) then
            rid = rankid
        else
            rid = 0
        end if

        call mpi_barrier(mpi_group_comm, mpi_ierr_group)
        call mpi_bcast(w, len(w), mpi_character, rid, mpi_group_comm, mpi_ierr_group)
        call mpi_barrier(mpi_group_comm, mpi_ierr_group)

    end subroutine bcast_array_1d_group_string

end module libflit_mpicomm_group
