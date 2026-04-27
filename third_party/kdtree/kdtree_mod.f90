!==============================================================================
! kdtree_mod.f90  —  Modern Fortran KD-tree package
!
! Provides fast spatial indexing for n-dimensional point clouds:
!   - Nearest neighbour          tree%query_nn(q, idx, dist)
!   - k nearest neighbours       tree%query_knn(q, k, idx, dists)
!   - Ball / range search        tree%query_ball(q, r, idx, n)
!
! Build:  O(n log n)   |   Query: O(log n) avg, O(n) worst case
!
! Accepted point-array types (resolved at compile time via generic build):
!   double precision, real, integer
!   All input is promoted to double precision internally.
!
!==============================================================================

!
!> K.G., 2026-04-25: All original integer(IP) were changed to integer
!
module kdtree_mod

    implicit none

    private

    !--------------------------------------------------------------------
    ! Public kind parameters — derived from Fortran native types.
    ! No iso_fortran_env dependency.
    !
    !   KD_DP  :  double precision kind   kind(1.0d0)
    !   KD_SP  :  single precision kind   kind(1.0e0)
    !   KD_IP  :  default integer kind    kind(0)
    !
    ! Import KD_DP / KD_IP for type-safe interop with query routines.
    !--------------------------------------------------------------------
    integer, parameter, public :: KD_DP = kind(1.0d0)
    integer, parameter, public :: KD_SP = kind(1.0e0)
    integer, parameter, public :: KD_IP = kind(0)

    ! Internal shorthand
    integer, parameter :: RP = KD_DP
    integer, parameter :: IP = KD_IP
    double precision, parameter :: KD_INF = huge(1.0d0)

    !--------------------------------------------------------------------
    ! Max-heap of capacity k.
    ! Stores the k *smallest* squared distances seen so far.
    ! Root always holds the current maximum -> O(1) rejection test.
    !--------------------------------------------------------------------
    type :: heap_t
        private
        double precision, allocatable :: d(:)
        integer, allocatable :: v(:)
        integer                   :: n = 0
    contains
        procedure :: h_init
        procedure :: h_push
        procedure :: h_top_dist
        procedure :: h_is_full
        procedure :: h_count
        procedure :: h_extract_sorted
    end type heap_t

    !--------------------------------------------------------------------
    ! KD-tree type
    !--------------------------------------------------------------------
    type, public :: kdtree
        !! KD-tree for spatial indexing of n-dimensional point clouds.
        !!
        !! build() accepts double precision, real, or integer point arrays.
        !! All arithmetic and returned distances are double precision.
        private
        integer :: ndim = 0
        integer :: npts = 0
        integer :: n_nodes = 0
        integer :: leaf_size = 10

        ! Flat node arrays, 1-indexed, DFS pre-order layout.
        ! nd_dim == 0  =>  leaf node
        ! nd_dim  > 0  =>  internal node, split along nd_dim
        integer, allocatable :: nd_dim(:)
        double precision, allocatable :: nd_val(:)
        integer, allocatable :: nd_left(:)
        integer, allocatable :: nd_right(:)
        integer, allocatable :: nd_lo(:)
        integer, allocatable :: nd_hi(:)

        double precision, allocatable :: pts(:, :)   ! (ndim, npts)
        integer, allocatable :: perm(:)

    contains
        ! Generic build — resolves to dp / real / integer overload at compile time
        generic, public :: build => build_dp, build_sp, build_int
        procedure, public :: query_nn => kd_query_nn
        procedure, public :: query_knn => kd_query_knn
        procedure, public :: query_ball => kd_query_ball
        procedure, public :: free => kd_free
        procedure, public :: built => kd_built
        procedure, public :: num_points => kd_num_points
        procedure, public :: num_dims => kd_num_dims

        procedure, private :: build_dp
        procedure, private :: build_sp
        procedure, private :: build_int
        procedure, private :: build_core
        procedure, private :: build_node
        procedure, private :: search_knn_node
        procedure, private :: search_ball_node
    end type kdtree

contains

    !==========================================================================
    ! Max-heap
    !==========================================================================

    subroutine h_init(h, k)
        class(heap_t), intent(inout) :: h
        integer, intent(in)    :: k
        if (allocated(h%d)) deallocate (h%d, h%v)
        allocate (h%d(k), h%v(k))
        h%n = 0
    end subroutine h_init

    logical function h_is_full(h)
        class(heap_t), intent(in) :: h
        h_is_full = (h%n == int(size(h%d), IP))
    end function h_is_full

    integer function h_count(h)
        class(heap_t), intent(in) :: h
        h_count = h%n
    end function h_count

    double precision function h_top_dist(h)
        class(heap_t), intent(in) :: h
        if (h%n > 0) then
            h_top_dist = h%d(1)
        else
            h_top_dist = KD_INF
        end if
    end function h_top_dist

    subroutine h_push(h, dist, idx)
        class(heap_t), intent(inout) :: h
        double precision, intent(in)    :: dist
        integer, intent(in)    :: idx
        integer      :: k, i, c
        double precision :: td
        integer      :: ti

        k = int(size(h%d), IP)
        if (h%n < k) then
            h%n = h%n + 1
            h%d(h%n) = dist
            h%v(h%n) = idx
            i = h%n
            do while (i > 1)
                if (h%d(i) > h%d(i/2)) then
                    td = h%d(i)
                    h%d(i) = h%d(i/2)
                    h%d(i/2) = td
                    ti = h%v(i)
                    h%v(i) = h%v(i/2)
                    h%v(i/2) = ti
                    i = i/2
                else
                    exit
                end if
            end do
        else if (dist < h%d(1)) then
            h%d(1) = dist
            h%v(1) = idx
            i = 1
            do
                c = 2*i
                if (c > h%n) exit
                if (c + 1 <= h%n .and. h%d(c + 1) > h%d(c)) c = c + 1
                if (h%d(c) > h%d(i)) then
                    td = h%d(i)
                    h%d(i) = h%d(c)
                    h%d(c) = td
                    ti = h%v(i)
                    h%v(i) = h%v(c)
                    h%v(c) = ti
                    i = c
                else
                    exit
                end if
            end do
        end if
    end subroutine h_push

    subroutine h_extract_sorted(h, dists, idxs)
        !! Copy heap into (dists, idxs), ascending by distance (insertion sort).
        class(heap_t), intent(in)  :: h
        double precision, intent(out) :: dists(:)
        integer, intent(out) :: idxs(:)
        integer      :: i, j
        double precision :: td
        integer      :: ti
        do i = 1, h%n
            dists(i) = h%d(i)
            idxs(i) = h%v(i)
        end do
        do i = 2, h%n
            td = dists(i)
            ti = idxs(i)
            j = i - 1
            do while (j >= 1 .and. dists(j) > td)
                dists(j + 1) = dists(j)
                idxs(j + 1) = idxs(j)
                j = j - 1
            end do
            dists(j + 1) = td
            idxs(j + 1) = ti
        end do
    end subroutine h_extract_sorted

    !==========================================================================
    ! Geometry utilities
    !==========================================================================

    integer function widest_dim(pts, perm, lo, hi, ndim)
        double precision, intent(in) :: pts(:, :)
        integer, intent(in) :: perm(:), lo, hi, ndim
        double precision :: mn, mx, rng, best
        integer      :: d, i
        widest_dim = 1
        best = -1.0d0
        do d = 1, ndim
            mn = pts(d, perm(lo))
            mx = mn
            do i = lo + 1, hi
                if (pts(d, perm(i)) < mn) mn = pts(d, perm(i))
                if (pts(d, perm(i)) > mx) mx = pts(d, perm(i))
            end do
            rng = mx - mn
            if (rng > best) then
                best = rng
                widest_dim = d
            end if
        end do
    end function widest_dim

    integer function med3(pts, perm, dim, a, b, c)
        double precision, intent(in) :: pts(:, :)
        integer, intent(in) :: perm(:), dim, a, b, c
        double precision :: va, vb, vc
        va = pts(dim, perm(a))
        vb = pts(dim, perm(b))
        vc = pts(dim, perm(c))
        if ((va <= vb .and. vb <= vc) .or. (vc <= vb .and. vb <= va)) then
            med3 = b
        else if ((vb <= va .and. va <= vc) .or. (vc <= va .and. va <= vb)) then
            med3 = a
        else
            med3 = c
        end if
    end function med3

    recursive subroutine quickselect(pts, perm, dim, lo, hi, k)
        double precision, intent(in)    :: pts(:, :)
        integer, intent(inout) :: perm(:)
        integer, intent(in)    :: dim, lo, hi, k
        integer      :: pi, i, j, tmp
        double precision :: pval
        if (lo >= hi) return
        pi = med3(pts, perm, dim, lo, (lo + hi)/2, hi)
        pval = pts(dim, perm(pi))
        tmp = perm(pi)
        perm(pi) = perm(hi)
        perm(hi) = tmp
        j = lo - 1
        do i = lo, hi - 1
            if (pts(dim, perm(i)) <= pval) then
                j = j + 1
                tmp = perm(i)
                perm(i) = perm(j)
                perm(j) = tmp
            end if
        end do
        j = j + 1
        tmp = perm(j)
        perm(j) = perm(hi)
        perm(hi) = tmp
        if (j > k) then
            call quickselect(pts, perm, dim, lo, j - 1, k)
        else if (j < k) then
            call quickselect(pts, perm, dim, j + 1, hi, k)
        end if
    end subroutine quickselect

    !==========================================================================
    ! Build overloads  (thin type-conversion wrappers -> build_core)
    !==========================================================================

    subroutine build_dp(self, points, leaf_size)
        !! Build from a double precision array, shape (ndim, npts).
        class(kdtree), intent(inout)        :: self
        double precision, intent(in)           :: points(:, :)
        integer, intent(in), optional :: leaf_size
        call self%build_core(points, leaf_size)
    end subroutine build_dp

    subroutine build_sp(self, points, leaf_size)
        !! Build from a real (single precision) array, shape (ndim, npts).
        !! Promoted to double precision; original array is unchanged.
        class(kdtree), intent(inout)        :: self
        real, intent(in)           :: points(:, :)
        integer, intent(in), optional :: leaf_size
        double precision, allocatable :: tmp(:, :)
        allocate (tmp(size(points, 1), size(points, 2)))
        tmp = dble(points)
        call self%build_core(tmp, leaf_size)
    end subroutine build_sp

    subroutine build_int(self, points, leaf_size)
        !! Build from an integer array, shape (ndim, npts).
        !! Converted to double precision; original array is unchanged.
        class(kdtree), intent(inout)        :: self
        integer, intent(in)           :: points(:, :)
        integer, intent(in), optional :: leaf_size
        double precision, allocatable :: tmp(:, :)
        allocate (tmp(size(points, 1), size(points, 2)))
        tmp = dble(points)
        call self%build_core(tmp, leaf_size)
    end subroutine build_int

    !==========================================================================
    ! Core build (double precision)
    !==========================================================================

    subroutine build_core(self, points, leaf_size)
        class(kdtree), intent(inout)        :: self
        double precision, intent(in)           :: points(:, :)
        integer, intent(in), optional :: leaf_size
        integer :: i, max_nodes

        call self%free()
        self%ndim = int(size(points, 1), IP)
        self%npts = int(size(points, 2), IP)
        if (self%ndim == 0 .or. self%npts == 0) return
        if (present(leaf_size)) self%leaf_size = max(1, leaf_size)

        allocate (self%pts(self%ndim, self%npts))
        self%pts = points

        allocate (self%perm(self%npts))
        do i = 1, self%npts
            self%perm(i) = i
        end do

        max_nodes = 4*self%npts/self%leaf_size + 16
        allocate (self%nd_dim(max_nodes))
        self%nd_dim = 0
        allocate (self%nd_val(max_nodes))
        self%nd_val = 0.0d0
        allocate (self%nd_left(max_nodes))
        self%nd_left = 0
        allocate (self%nd_right(max_nodes))
        self%nd_right = 0
        allocate (self%nd_lo(max_nodes))
        self%nd_lo = 0
        allocate (self%nd_hi(max_nodes))
        self%nd_hi = 0

        self%n_nodes = 0
        call self%build_node(1, self%npts)
    end subroutine build_core

    recursive subroutine build_node(self, lo, hi)
        class(kdtree), intent(inout) :: self
        integer, intent(in)    :: lo, hi
        integer :: node, mid, sdim, left_id, right_id

        self%n_nodes = self%n_nodes + 1
        node = self%n_nodes
        self%nd_lo(node) = lo
        self%nd_hi(node) = hi
        self%nd_left(node) = 0
        self%nd_right(node) = 0
        self%nd_dim(node) = 0

        if (hi - lo + 1 <= self%leaf_size) return

        sdim = widest_dim(self%pts, self%perm, lo, hi, self%ndim)
        mid = (lo + hi)/2
        call quickselect(self%pts, self%perm, sdim, lo, hi, mid)

        self%nd_dim(node) = sdim
        self%nd_val(node) = self%pts(sdim, self%perm(mid))

        left_id = self%n_nodes + 1
        call self%build_node(lo, mid)
        self%nd_left(node) = left_id

        right_id = self%n_nodes + 1
        call self%build_node(mid + 1, hi)
        self%nd_right(node) = right_id
    end subroutine build_node

    !==========================================================================
    ! Query routines
    !==========================================================================

    subroutine kd_query_nn(self, query, nn_idx, nn_dist)
        !! Single nearest neighbour.
        !! query: shape (ndim), double precision.
        !! nn_idx: 1-based index into original points array.
        !! nn_dist: Euclidean distance, double precision.
        class(kdtree), intent(in)  :: self
        double precision, intent(in)  :: query(:)
        integer, intent(out) :: nn_idx
        double precision, intent(out) :: nn_dist
        type(heap_t)     :: h
        double precision :: d(1)
        integer      :: v(1)
        if (.not. self%built()) then
            nn_idx = -1
            nn_dist = KD_INF
            return
        end if
        call h%h_init(1)
        call self%search_knn_node(1, query, h)
        call h%h_extract_sorted(d, v)
        nn_idx = int(v(1))
        nn_dist = sqrt(d(1))
    end subroutine kd_query_nn

    subroutine kd_query_knn(self, query, k, indices, distances)
        !! k nearest neighbours, sorted ascending by Euclidean distance.
        !! query: shape (ndim), double precision.
        !! If k > num_points(), only num_points() results are written.
        class(kdtree), intent(in)  :: self
        double precision, intent(in)  :: query(:)
        integer, intent(in)  :: k
        integer, intent(out) :: indices(:)
        double precision, intent(out) :: distances(:)
        type(heap_t) :: h
        integer  :: kk, i
        if (.not. self%built()) then
            indices = -1
            distances = KD_INF
            return
        end if
        kk = min(int(k, IP), self%npts)
        call h%h_init(kk)
        call self%search_knn_node(1, query, h)
        call h%h_extract_sorted(distances(1:kk), indices(1:kk))
        do i = 1, kk
            distances(i) = sqrt(distances(i))
        end do
    end subroutine kd_query_knn

    subroutine kd_query_ball(self, query, radius, indices, n_found)
        !! All points within Euclidean radius of query.
        !! query: shape (ndim), double precision.
        !! indices: allocatable output, 1-based. Caller deallocates.
        class(kdtree), intent(in)               :: self
        double precision, intent(in)               :: query(:)
        double precision, intent(in)               :: radius
        integer, intent(out), allocatable :: indices(:)
        integer, intent(out)              :: n_found
        integer, allocatable :: tmp(:)
        if (.not. self%built()) then
            n_found = 0
            allocate (indices(0))
            return
        end if
        allocate (tmp(self%npts))
        n_found = 0
        call self%search_ball_node(1, query, radius*radius, tmp, n_found)
        if (allocated(indices)) deallocate (indices)
        allocate (indices(n_found))
        if (n_found > 0) indices(1:n_found) = int(tmp(1:n_found))
        deallocate (tmp)
    end subroutine kd_query_ball

    pure logical function kd_built(self)
        class(kdtree), intent(in) :: self
        kd_built = (self%n_nodes > 0)
    end function kd_built

    pure integer function kd_num_points(self)
        class(kdtree), intent(in) :: self
        kd_num_points = int(self%npts)
    end function kd_num_points

    pure integer function kd_num_dims(self)
        class(kdtree), intent(in) :: self
        kd_num_dims = int(self%ndim)
    end function kd_num_dims

    subroutine kd_free(self)
        class(kdtree), intent(inout) :: self
        self%ndim = 0
        self%npts = 0
        self%n_nodes = 0
        self%leaf_size = 10
        if (allocated(self%pts)) deallocate (self%pts)
        if (allocated(self%perm)) deallocate (self%perm)
        if (allocated(self%nd_dim)) deallocate (self%nd_dim)
        if (allocated(self%nd_val)) deallocate (self%nd_val)
        if (allocated(self%nd_left)) deallocate (self%nd_left)
        if (allocated(self%nd_right)) deallocate (self%nd_right)
        if (allocated(self%nd_lo)) deallocate (self%nd_lo)
        if (allocated(self%nd_hi)) deallocate (self%nd_hi)
    end subroutine kd_free

    !==========================================================================
    ! Private search helpers
    !==========================================================================

    recursive subroutine search_knn_node(self, node, query, h)
        class(kdtree), intent(in)    :: self
        integer, intent(in)    :: node
        double precision, intent(in)    :: query(:)
        type(heap_t), intent(inout) :: h
        integer      :: i, pt, sdim, near_c, far_c
        double precision :: d2, diff

        if (self%nd_dim(node) == 0) then
            do i = self%nd_lo(node), self%nd_hi(node)
                pt = self%perm(i)
                d2 = sum((self%pts(:, pt) - query)**2)
                call h%h_push(d2, pt)
            end do
            return
        end if

        sdim = self%nd_dim(node)
        diff = query(sdim) - self%nd_val(node)
        if (diff <= 0.0d0) then
            near_c = self%nd_left(node)
            far_c = self%nd_right(node)
        else
            near_c = self%nd_right(node)
            far_c = self%nd_left(node)
        end if
        call self%search_knn_node(near_c, query, h)
        if (.not. h%h_is_full() .or. diff*diff < h%h_top_dist()) &
            call self%search_knn_node(far_c, query, h)
    end subroutine search_knn_node

    recursive subroutine search_ball_node(self, node, query, r2, buf, n)
        class(kdtree), intent(in)    :: self
        integer, intent(in)    :: node
        double precision, intent(in)    :: query(:)
        double precision, intent(in)    :: r2
        integer, intent(inout) :: buf(:)
        integer, intent(inout) :: n
        integer      :: i, pt
        double precision :: d2, diff

        if (self%nd_dim(node) == 0) then
            do i = self%nd_lo(node), self%nd_hi(node)
                pt = self%perm(i)
                d2 = sum((self%pts(:, pt) - query)**2)
                if (d2 <= r2) then
                    n = n + 1
                    buf(n) = pt
                end if
            end do
            return
        end if

        diff = query(self%nd_dim(node)) - self%nd_val(node)
        if (diff <= 0.0d0) then
            call self%search_ball_node(self%nd_left(node), query, r2, buf, n)
            if (diff*diff <= r2) &
                call self%search_ball_node(self%nd_right(node), query, r2, buf, n)
        else
            call self%search_ball_node(self%nd_right(node), query, r2, buf, n)
            if (diff*diff <= r2) &
                call self%search_ball_node(self%nd_left(node), query, r2, buf, n)
        end if
    end subroutine search_ball_node

end module kdtree_mod
