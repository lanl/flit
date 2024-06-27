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


module libflit_mpicomm

    use mpi
    ! using mpi_f08 seems to cause mysterious errors in mpi_sendrecv
    use libflit_array
    use libflit_domain_decomposition
    use libflit_error
    use libflit_array_operation

    implicit none

    interface bcast_array
        module procedure :: bcast_array_1d_global_string
        module procedure :: bcast_array_1d_global_int
        module procedure :: bcast_array_2d_global_int
        module procedure :: bcast_array_3d_global_int
        module procedure :: bcast_array_1d_global_double
        module procedure :: bcast_array_2d_global_double
        module procedure :: bcast_array_3d_global_double
        module procedure :: bcast_array_1d_global_float
        module procedure :: bcast_array_2d_global_float
        module procedure :: bcast_array_3d_global_float
        module procedure :: bcast_array_1d_global_complex
        module procedure :: bcast_array_2d_global_complex
        module procedure :: bcast_array_3d_global_complex
    end interface

    interface commute_array
        module procedure :: commute_array_1d_global_int
        module procedure :: commute_array_2d_global_int
        module procedure :: commute_array_3d_global_int
        module procedure :: commute_array_1d_global_double
        module procedure :: commute_array_2d_global_double
        module procedure :: commute_array_3d_global_double
        module procedure :: commute_array_1d_global_float
        module procedure :: commute_array_2d_global_float
        module procedure :: commute_array_3d_global_float
        module procedure :: commute_array_1d_global_complex
        module procedure :: commute_array_2d_global_complex
        module procedure :: commute_array_3d_global_complex
    end interface commute_array

    interface allreduce
        ! Scalar
        module procedure :: gather_distribute_global_int
        module procedure :: gather_distribute_global_double
        module procedure :: gather_distribute_global_float
        module procedure :: gather_distribute_global_complex
        module procedure :: gather_distribute_global_dcomplex
    end interface

    interface allreduce_array
        ! Array
        module procedure :: gather_distribute_array_1d_global_int
        module procedure :: gather_distribute_array_2d_global_int
        module procedure :: gather_distribute_array_3d_global_int
        module procedure :: gather_distribute_array_1d_global_double
        module procedure :: gather_distribute_array_2d_global_double
        module procedure :: gather_distribute_array_3d_global_double
        module procedure :: gather_distribute_array_1d_global_float
        module procedure :: gather_distribute_array_2d_global_float
        module procedure :: gather_distribute_array_3d_global_float
        module procedure :: gather_distribute_array_1d_global_complex
        module procedure :: gather_distribute_array_2d_global_complex
        module procedure :: gather_distribute_array_3d_global_complex
        module procedure :: gather_distribute_array_1d_global_dcomplex
        module procedure :: gather_distribute_array_2d_global_dcomplex
        module procedure :: gather_distribute_array_3d_global_dcomplex
        ! For array with size = int8
        module procedure :: gather_distribute_large_array_1d_global_int
        module procedure :: gather_distribute_large_array_1d_global_double
        module procedure :: gather_distribute_large_array_1d_global_float
        module procedure :: gather_distribute_large_array_1d_global_complex
        module procedure :: gather_distribute_large_array_1d_global_dcomplex
    end interface allreduce_array

    interface reduce
        module procedure :: gather_global_int
        module procedure :: gather_global_double
        module procedure :: gather_global_float
        module procedure :: gather_global_complex
        module procedure :: gather_global_dcomplex
    end interface reduce

    interface reduce_array
        module procedure :: gather_array_1d_global_int
        module procedure :: gather_array_2d_global_int
        module procedure :: gather_array_3d_global_int
        module procedure :: gather_array_1d_global_double
        module procedure :: gather_array_2d_global_double
        module procedure :: gather_array_3d_global_double
        module procedure :: gather_array_1d_global_float
        module procedure :: gather_array_2d_global_float
        module procedure :: gather_array_3d_global_float
        module procedure :: gather_array_1d_global_complex
        module procedure :: gather_array_2d_global_complex
        module procedure :: gather_array_3d_global_complex
        module procedure :: gather_array_1d_global_dcomplex
        module procedure :: gather_array_2d_global_dcomplex
        module procedure :: gather_array_3d_global_dcomplex
    end interface reduce_array

    interface domain_decomp_regular
        module procedure :: domain_decomp_regular_2d
        module procedure :: domain_decomp_regular_3d
    end interface domain_decomp_regular

    !    interface mpi_global_min
    !        module procedure :: mpi_global_min_1d_int
    !        module procedure :: mpi_global_min_1d_float
    !        module procedure :: mpi_global_min_1d_double
    !    end interface mpi_global_min

    integer, public :: rank1, rank2, rank3
    integer, public :: block_x1left, block_x1right
    integer, public :: block_x2left, block_x2right
    integer, public :: block_x3left, block_x3right
    integer, public :: groupid
    integer, public :: rankid
    integer, public :: blockid
    integer, public :: nrank
    integer, public :: ngroup
    type(mpi_status), public :: mpi_stats
    integer, public :: mpi_ierr

    public :: mpistart, mpiend, mpibarrier, mpistop
    public :: bcast_array
    public :: reduce
    public :: reduce_array
    public :: allreduce
    public :: allreduce_array
    public :: commute_array
    public :: domain_decomp_regular

contains

    subroutine mpistart

        call mpi_init(mpi_ierr)
        call mpi_comm_size(mpi_comm_world, nrank, mpi_ierr)
        call mpi_comm_rank(mpi_comm_world, rankid, mpi_ierr)

    end subroutine mpistart

    subroutine mpiend

        call mpi_barrier(mpi_comm_world, mpi_ierr)
        call mpi_finalize(mpi_ierr)

    end subroutine mpiend

    subroutine mpibarrier

        call mpi_barrier(mpi_comm_world, mpi_ierr)

    end subroutine mpibarrier

    subroutine mpistop

        call mpi_abort(mpi_comm_world, mpi_err_other, mpi_ierr)

    end subroutine mpistop

    !
    !> Decompose a matrix into MPI blocks
    !
    subroutine domain_decomp_regular_2d(n1, n2, n1beg, n1end, n2beg, n2end)

        integer, intent(in) :: n1, n2
        integer, intent(out) :: n1beg, n1end, n2beg, n2end

        integer, allocatable, dimension(:, :) :: blk1, blk2
        integer :: i, j

        call mpi_comm_rank(mpi_comm_world, rankid, mpi_ierr)

        call cut(1, n1, rank1, blk1)
        call cut(1, n2, rank2, blk2)

        if (rankid > rank1*rank2 - 1) then
            n1beg = 1
            n1end = -1
            n2beg = 1
            n2end = -1
        else
            do j = 0, rank2 - 1
                do i = 0, rank1 - 1
                    if (rankid == j*rank1 + i) then
                        n1beg = blk1(i + 1, 1)
                        n1end = blk1(i + 1, 2)
                        n2beg = blk2(j + 1, 1)
                        n2end = blk2(j + 1, 2)
                    end if
                end do
            end do

            ! Find adjacent blocks for each group
            do j = 0, rank2 - 1
                do i = 0, rank1 - 1

                    blockid = j*rank1 + i

                    if (blockid == rankid) then

                        block_x1left = rankid - 1
                        block_x1right = rankid + 1
                        block_x2left = rankid - rank1
                        block_x2right = rankid + rank1

                        if (i == 0) then
                            block_x1left = mpi_proc_null
                        end if
                        if (i == rank1 - 1) then
                            block_x1right = mpi_proc_null
                        end if
                        if (j == 0) then
                            block_x2left = mpi_proc_null
                        end if
                        if (j == rank2 - 1) then
                            block_x2right = mpi_proc_null
                        end if

                    end if

                end do
            end do

        end if

    end subroutine domain_decomp_regular_2d

    !
    !> Decompose a volume into MPI blocks
    !
    subroutine domain_decomp_regular_3d(n1, n2, n3, n1beg, n1end, n2beg, n2end, n3beg, n3end)

        integer, intent(in) :: n1, n2, n3
        integer, intent(out) :: n1beg, n1end, n2beg, n2end, n3beg, n3end

        integer, allocatable, dimension(:, :) :: blk1, blk2, blk3
        integer :: i, j, k

        call mpi_comm_rank(mpi_comm_world, rankid, mpi_ierr)

        call cut(1, n1, rank1, blk1)
        call cut(1, n2, rank2, blk2)
        call cut(1, n3, rank3, blk3)

        if (rankid > rank1*rank2*rank3 - 1) then
            n1beg = 1
            n1end = -1
            n2beg = 1
            n2end = -1
            n3beg = 1
            n3end = -1
        else
            do k = 0, rank3 - 1
                do j = 0, rank2 - 1
                    do i = 0, rank1 - 1
                        if (rankid == k*rank1*rank2 + j*rank1 + i) then
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
            do k = 0, rank3 - 1
                do j = 0, rank2 - 1
                    do i = 0, rank1 - 1

                        blockid = k*rank2*rank1 + j*rank1 + i

                        if (blockid == rankid) then

                            block_x1left = rankid - 1
                            block_x1right = rankid + 1
                            block_x2left = rankid - rank1
                            block_x2right = rankid + rank1
                            block_x3left = rankid - rank1*rank2
                            block_x3right = rankid + rank1*rank2

                            if (i == 0) then
                                block_x1left = mpi_proc_null
                            end if
                            if (i == rank1 - 1) then
                                block_x1right = mpi_proc_null
                            end if
                            if (j == 0) then
                                block_x2left = mpi_proc_null
                            end if
                            if (j == rank2 - 1) then
                                block_x2right = mpi_proc_null
                            end if
                            if (k == 0) then
                                block_x3left = mpi_proc_null
                            end if
                            if (k == rank3 - 1) then
                                block_x3right = mpi_proc_null
                            end if

                        end if

                    end do
                end do
            end do

        end if

    end subroutine domain_decomp_regular_3d

    ! MPI array operations
#define T int
#define TT integer
#define TTT mpi_integer
#include "template_mpicomm.f90"

#define T float
#define TT real
#define TTT mpi_real
#include "template_mpicomm.f90"

#define T double
#define TT double precision
#define TTT mpi_double_precision
#include "template_mpicomm.f90"

#define T complex
#define TT complex
#define TTT mpi_complex
#include "template_mpicomm.f90"

#define T dcomplex
#define TT double complex
#define TTT mpi_double_complex
#include "template_mpicomm.f90"

    ! For string
    subroutine bcast_array_1d_global_string(w, rankid)

        character(len=*), intent(inout) :: w
        integer, intent(in), optional :: rankid

        integer :: rid

        if (present(rankid)) then
            rid = rankid
        else
            rid = 0
        end if

        call mpi_barrier(mpi_comm_world, mpi_ierr)
        call mpi_bcast(w, len(w), mpi_character, rid, mpi_comm_world, mpi_ierr)
        call mpi_barrier(mpi_comm_world, mpi_ierr)

    end subroutine bcast_array_1d_global_string

end module libflit_mpicomm
