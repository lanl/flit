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


#define PASTE(X)            X
#define PASTE2(X)           PASTE(X)_
#define CONCATHELP(X, Y)    PASTE2(X)Y
#define CONCAT(X, Y)        CONCATHELP(X, Y)

#define bcast_array_1d_global_     CONCAT(bcast_array_1d_global, T)
#define bcast_array_2d_global_     CONCAT(bcast_array_2d_global, T)
#define bcast_array_3d_global_     CONCAT(bcast_array_3d_global, T)
#define gather_global_     CONCAT(gather_global, T)
#define gather_array_1d_global_     CONCAT(gather_array_1d_global, T)
#define gather_array_2d_global_     CONCAT(gather_array_2d_global, T)
#define gather_array_3d_global_     CONCAT(gather_array_3d_global, T)
#define commute_array_1d_global_     CONCAT(commute_array_1d_global, T)
#define commute_array_2d_global_     CONCAT(commute_array_2d_global, T)
#define commute_array_3d_global_     CONCAT(commute_array_3d_global, T)
#define gather_distribute_global_     CONCAT(gather_distribute_global, T)
#define gather_distribute_array_1d_global_     CONCAT(gather_distribute_array_1d_global, T)
#define gather_distribute_array_2d_global_     CONCAT(gather_distribute_array_2d_global, T)
#define gather_distribute_array_3d_global_     CONCAT(gather_distribute_array_3d_global, T)
#define gather_distribute_large_array_1d_global_     CONCAT(gather_distribute_large_array_1d_global, T)

!
! Broadcast
!
subroutine bcast_array_1d_global_(w, source)

    TT, dimension(:), intent(inout) :: w
    integer, intent(in), optional :: source

    integer :: rid

    if (present(source)) then
        rid = source
    else
        rid = 0
    end if

    call mpi_barrier(mpi_comm_world, mpi_ierr)
    call mpi_bcast(w, size(w), TTT, rid, mpi_comm_world, mpi_ierr)
    call mpi_barrier(mpi_comm_world, mpi_ierr)

end subroutine bcast_array_1d_global_

subroutine bcast_array_2d_global_(w, source)

    TT, dimension(:, :), intent(inout) :: w
    integer, intent(in), optional :: source

    integer :: rid

    if (present(source)) then
        rid = source
    else
        rid = 0
    end if

    call mpi_barrier(mpi_comm_world, mpi_ierr)
    call mpi_bcast(w, size(w), TTT, rid, mpi_comm_world, mpi_ierr)
    call mpi_barrier(mpi_comm_world, mpi_ierr)

end subroutine bcast_array_2d_global_

subroutine bcast_array_3d_global_(w, source)

    TT, dimension(:, :, :), intent(inout) :: w
    integer, intent(in), optional :: source

    integer :: rid

    if (present(source)) then
        rid = source
    else
        rid = 0
    end if

    call mpi_barrier(mpi_comm_world, mpi_ierr)
    call mpi_bcast(w, size(w), TTT, rid, mpi_comm_world, mpi_ierr)
    call mpi_barrier(mpi_comm_world, mpi_ierr)

end subroutine bcast_array_3d_global_

!
! Reduce
!
subroutine gather_global_(w, target)

    TT, intent(inout) :: w
    integer, intent(in), optional :: target

    TT :: wlocal
    integer :: rid

    if (present(target)) then
        rid = target
    else
        rid = 0
    end if

    wlocal = w

    call mpi_barrier(mpi_comm_world, mpi_ierr)
    call mpi_reduce(w, wlocal, 1, TTT, mpi_sum, rid, mpi_comm_world, mpi_ierr)
    call mpi_barrier(mpi_comm_world, mpi_ierr)

    w = wlocal
    call mpi_barrier(mpi_comm_world, mpi_ierr)

end subroutine gather_global_

subroutine gather_array_1d_global_(w, target)

    TT, dimension(:), intent(inout) :: w
    integer, intent(in), optional :: target

    TT, allocatable, dimension(:) :: wlocal
    integer :: rid

    if (present(target)) then
        rid = target
    else
        rid = 0
    end if

    wlocal = w

    call mpi_barrier(mpi_comm_world, mpi_ierr)
    call mpi_reduce(w, wlocal, size(w), TTT, mpi_sum, rid, mpi_comm_world, mpi_ierr)
    call mpi_barrier(mpi_comm_world, mpi_ierr)

    w = wlocal
    call mpi_barrier(mpi_comm_world, mpi_ierr)

end subroutine gather_array_1d_global_

subroutine gather_array_2d_global_(w, target)

    TT, dimension(:, :), intent(inout) :: w
    integer, intent(in), optional :: target

    TT, allocatable, dimension(:, :) :: wlocal
    integer :: rid

    if (present(target)) then
        rid = target
    else
        rid = 0
    end if

    wlocal = w

    call mpi_barrier(mpi_comm_world, mpi_ierr)
    call mpi_reduce(w, wlocal, size(w), TTT, mpi_sum, rid, mpi_comm_world, mpi_ierr)
    call mpi_barrier(mpi_comm_world, mpi_ierr)

    w = wlocal
    call mpi_barrier(mpi_comm_world, mpi_ierr)

end subroutine gather_array_2d_global_

subroutine gather_array_3d_global_(w, target)

    TT, dimension(:, :, :), intent(inout) :: w
    integer, intent(in), optional :: target

    TT, allocatable, dimension(:, :, :) :: wlocal
    integer :: rid

    if (present(target)) then
        rid = target
    else
        rid = 0
    end if

    wlocal = w

    call mpi_barrier(mpi_comm_world, mpi_ierr)
    call mpi_reduce(w, wlocal, size(w), TTT, mpi_sum, rid, mpi_comm_world, mpi_ierr)
    call mpi_barrier(mpi_comm_world, mpi_ierr)

    w = wlocal
    call mpi_barrier(mpi_comm_world, mpi_ierr)

end subroutine gather_array_3d_global_

!
! Send-recv
!
subroutine commute_array_1d_global_(w, nl)

    TT, allocatable, dimension(:), intent(inout) :: w
    integer, intent(in) :: nl

    integer :: blks1, n1l, n1u

    n1l = lbound(w, 1) + nl
    n1u = ubound(w, 1) - nl

    if (rank1 > 1) then

        blks1 = nl

        call mpi_sendrecv(w(n1l:n1l + (nl - 1)), blks1, TTT, block_x1left, 1, &
            w(n1u + 1:n1u + nl), blks1, TTT, block_x1right, 1, mpi_comm_world, mpi_stats, mpi_ierr)
        call mpi_sendrecv(w(n1u - (nl - 1):n1u), blks1, TTT, block_x1right, 2, &
            w(n1l - nl:n1l - 1), blks1, TTT, block_x1left, 2, mpi_comm_world, mpi_stats, mpi_ierr)

    end if

end subroutine commute_array_1d_global_

subroutine commute_array_2d_global_(w, nl, dim)

    TT, allocatable, dimension(:, :), intent(inout) :: w
    integer, intent(in) :: nl
    integer, intent(in), optional :: dim

    integer :: blks1, blks2, n1l, n1u, n2l, n2u
    integer :: axis

    n1l = lbound(w, 1) + nl
    n1u = ubound(w, 1) - nl
    n2l = lbound(w, 2) + nl
    n2u = ubound(w, 2) - nl

    if (present(dim)) then
        axis = dim
    else
        axis = 0
    end if

    if (rank1 > 1 .and. (axis == 0 .or. axis == 1)) then

        blks1 = (n2u - n2l + 1 + 2*nl)*nl

        call mpi_sendrecv(w(n1l:n1l + (nl - 1), n2l - nl:n2u + nl), blks1, TTT, block_x1left, 1, &
            w(n1u + 1:n1u + nl, n2l - nl:n2u + nl), blks1, TTT, block_x1right, 1, mpi_comm_world, mpi_stats, mpi_ierr)
        call mpi_sendrecv(w(n1u - (nl - 1):n1u, n2l - nl:n2u + nl), blks1, TTT, block_x1right, 2, &
            w(n1l - nl:n1l - 1, n2l - nl:n2u + nl), blks1, TTT, block_x1left, 2, mpi_comm_world, mpi_stats, mpi_ierr)

    end if

    if (rank2 > 1 .and. (axis == 0 .or. axis == 2)) then

        blks2 = (n1u - n1l + 1 + 2*nl)*nl

        call mpi_sendrecv(w(n1l - nl:n1u + nl, n2l:n2l + (nl - 1)), blks2, TTT, block_x2left, 3, &
            w(n1l - nl:n1u + nl, n2u + 1:n2u + nl), blks2, TTT, block_x2right, 3, mpi_comm_world, mpi_stats, mpi_ierr)
        call mpi_sendrecv(w(n1l - nl:n1u + nl, n2u - (nl - 1):n2u), blks2, TTT, block_x2right, 4, &
            w(n1l - nl:n1u + nl, n2l - nl:n2l - 1), blks2, TTT, block_x2left, 4, mpi_comm_world, mpi_stats, mpi_ierr)

    end if

end subroutine commute_array_2d_global_

subroutine commute_array_3d_global_(w, nl, dim)

    TT, allocatable, dimension(:, :, :), intent(inout) :: w
    integer, intent(in) :: nl
    integer, intent(in), optional :: dim

    integer :: blks1, blks2, blks3, n1l, n1u, n2l, n2u, n3l, n3u
    integer :: axis

    n1l = lbound(w, 1) + nl
    n1u = ubound(w, 1) - nl
    n2l = lbound(w, 2) + nl
    n2u = ubound(w, 2) - nl
    n3l = lbound(w, 3) + nl
    n3u = ubound(w, 3) - nl

    if (present(dim)) then
        axis = dim
    else
        axis = 0
    end if

    if (rank1 > 1 .and. (axis == 0 .or. axis == 1)) then

        blks1 = (n2u - n2l + 1 + 2*nl)*(n3u - n3l + 1 + 2*nl)*nl

        call mpi_sendrecv(w(n1l:n1l + (nl - 1), n2l - nl:n2u + nl, n3l - nl:n3u + nl), blks1, TTT, block_x1left, 1, &
            w(n1u + 1:n1u + nl, n2l - nl:n2u + nl, n3l - nl:n3u + nl), blks1, TTT, block_x1right, 1, mpi_comm_world, mpi_stats, mpi_ierr)
        call mpi_sendrecv(w(n1u - (nl - 1):n1u, n2l - nl:n2u + nl, n3l - nl:n3u + nl), blks1, TTT, block_x1right, 2, &
            w(n1l - nl:n1l - 1, n2l - nl:n2u + nl, n3l - nl:n3u + nl), blks1, TTT, block_x1left, 2, mpi_comm_world, mpi_stats, mpi_ierr)

    end if

    if (rank2 > 1 .and. (axis == 0 .or. axis == 2)) then

        blks2 = (n1u - n1l + 1 + 2*nl)*(n3u - n3l + 1 + 2*nl)*nl

        call mpi_sendrecv(w(n1l - nl:n1u + nl, n2l:n2l + (nl - 1), n3l - nl:n3u + nl), blks2, TTT, block_x2left, 3, &
            w(n1l - nl:n1u + nl, n2u + 1:n2u + nl, n3l - nl:n3u + nl), blks2, TTT, block_x2right, 3, mpi_comm_world, mpi_stats, mpi_ierr)
        call mpi_sendrecv(w(n1l - nl:n1u + nl, n2u - (nl - 1):n2u, n3l - nl:n3u + nl), blks2, TTT, block_x2right, 4, &
            w(n1l - nl:n1u + nl, n2l - nl:n2l - 1, n3l - nl:n3u + nl), blks2, TTT, block_x2left, 4, mpi_comm_world, mpi_stats, mpi_ierr)

    end if

    if (rank3 > 1 .and. (axis == 0 .or. axis == 3)) then

        blks3 = (n2u - n2l + 1 + 2*nl)*(n1u - n1l + 1 + 2*nl)*nl

        call mpi_sendrecv(w(n1l - nl:n1u + nl, n2l - nl:n2u + nl, n3l:n3l + (nl - 1)), blks3, TTT, block_x3left, 5, &
            w(n1l - nl:n1u + nl, n2l - nl:n2u + nl, n3u + 1:n3u + nl), blks3, TTT, block_x3right, 5, mpi_comm_world, mpi_stats, mpi_ierr)
        call mpi_sendrecv(w(n1l - nl:n1u + nl, n2l - nl:n2u + nl, n3u - (nl - 1):n3u), blks3, TTT, block_x3right, 6, &
            w(n1l - nl:n1u + nl, n2l - nl:n2u + nl, n3l - nl:n3l - 1), blks3, TTT, block_x3left, 6, mpi_comm_world, mpi_stats, mpi_ierr)

    end if

end subroutine commute_array_3d_global_

!
! Allreduce
!
subroutine gather_distribute_global_(w)

    TT, intent(inout) :: w

    call mpi_barrier(mpi_comm_world, mpi_ierr)
    call mpi_allreduce(mpi_in_place, w, 1, TTT, mpi_sum, mpi_comm_world, mpi_ierr)
    call mpi_barrier(mpi_comm_world, mpi_ierr)

end subroutine gather_distribute_global_

subroutine gather_distribute_array_1d_global_(w)

    TT, dimension(:), intent(inout) :: w

    call mpi_barrier(mpi_comm_world, mpi_ierr)
    call mpi_allreduce(mpi_in_place, w, size(w), TTT, mpi_sum, mpi_comm_world, mpi_ierr)
    call mpi_barrier(mpi_comm_world, mpi_ierr)

end subroutine gather_distribute_array_1d_global_

subroutine gather_distribute_large_array_1d_global_(w, nblock)

    TT, dimension(:), intent(inout) :: w
    integer, intent(in) :: nblock

    integer(kind=8) :: n1
    integer(kind=8), allocatable, dimension(:, :) :: blkrange
    integer :: i

    allocate (blkrange(0:nblock - 1, 1:2))
    n1 = size(w, 1)

    call cut(int(1, kind=8), n1, int(nblock, kind=8), blkrange)

    do i = 0, nblock - 1

        call mpi_barrier(mpi_comm_world, mpi_ierr)
        call mpi_allreduce(mpi_in_place, w(blkrange(i, 1):blkrange(i, 2)), int(blkrange(i, 2) - blkrange(i, 1) + 1, kind=4), &
            TTT, mpi_sum, mpi_comm_world, mpi_ierr)
        call mpi_barrier(mpi_comm_world, mpi_ierr)

    end do

end subroutine gather_distribute_large_array_1d_global_

subroutine gather_distribute_array_2d_global_(w)

    TT, dimension(:, :), intent(inout) :: w

    call mpi_barrier(mpi_comm_world, mpi_ierr)
    call mpi_allreduce(mpi_in_place, w, size(w), TTT, mpi_sum, mpi_comm_world, mpi_ierr)
    call mpi_barrier(mpi_comm_world, mpi_ierr)

end subroutine gather_distribute_array_2d_global_

subroutine gather_distribute_array_3d_global_(w)

    TT, dimension(:, :, :), intent(inout) :: w

    call mpi_barrier(mpi_comm_world, mpi_ierr)
    call mpi_allreduce(mpi_in_place, w, size(w), TTT, mpi_sum, mpi_comm_world, mpi_ierr)
    call mpi_barrier(mpi_comm_world, mpi_ierr)

end subroutine gather_distribute_array_3d_global_

#undef T
#undef TT
#undef TTT
#undef DEFAULT_VALUE

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef bcast_array_1d_global_
#undef bcast_array_2d_global_
#undef bcast_array_3d_global_
#undef gather_global_
#undef gather_array_1d_global_
#undef gather_array_2d_global_
#undef gather_array_3d_global_
#undef commute_array_1d_global_
#undef commute_array_2d_global_
#undef commute_array_3d_global_
#undef gather_distribute_global_
#undef gather_distribute_array_1d_global_
#undef gather_distribute_array_2d_global_
#undef gather_distribute_array_3d_global_
#undef gather_distribute_large_array_1d_global_
!#undef mpi_global_min_1d_
