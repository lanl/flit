!>
!> Three-dimensional convex hull via the randomized incremental algorithm.
!>
!> Algorithm outline:
!>  1. Build an initial tetrahedron from four affinely independent points,
!>     chosen as extremes along each axis plus the most-distant points from
!>     the line, then the plane.
!>  2. Insert the remaining points in random order. For each new point p:
!>       a. Find all faces visible from p (signed distance > eps).
!>       b. If none, p is inside the current hull.
!>       c. Otherwise, locate horizon edges (edges between a visible face
!>          and a non-visible neighbor).
!>       d. Delete the visible faces and cap the hole with new triangles
!>          fanning from p to each horizon edge.
!>       e. Re-link face adjacency: new face `nbr(3)` points back at the
!>          retained face across the horizon edge; the other two point at
!>          the adjacent new faces around the fan.
!>  3. Compact the face list, removing logically deleted entries.
!>
!> Conventions:
!>  * Each face stores its three vertex indices CCW as seen from *outside*
!>    the hull.
!>  * nbr(k) is the neighbor across the edge *opposite* vertex v(k).
!>    That edge is (v(k+1), v(k+2)) with cyclic indexing.
!>  * Plane equation: n . x = off, with n the outward unit normal.
!>
!> Complexity: expected O(n log n) with randomization (O(n^2) worst case).
!>

module convex_hull_3d_mod

    use convex_hull_kinds, only: dp
    implicit none
    private
    public :: convex_hull_3d, face3d_t

    !> A triangular face on the 3-D hull.
    type :: face3d_t
        integer  :: v(3) = 0          !! vertex indices (CCW from outside)
        integer  :: nbr(3) = 0          !! neighbor face indices
        real(dp) :: n(3) = 0.0_dp     !! outward unit normal
        real(dp) :: off = 0.0_dp     !! plane offset: n . x = off
        logical  :: live = .true.     !! internal use; always .true. on output
    end type face3d_t

    !> Internal working state used while building the hull.
    type :: hull_state_t
        type(face3d_t), allocatable :: f(:)
        integer                     :: nf = 0   !! highest face index in use
        integer                     :: cap = 0   !! allocated capacity
    end type hull_state_t

    ! Status codes returned via the optional `stat` argument.
    integer, parameter, public :: HULL_OK = 0
    integer, parameter, public :: HULL_TOO_FEW_POINTS = -1
    integer, parameter, public :: HULL_ALL_COINCIDENT = 1
    integer, parameter, public :: HULL_ALL_COLLINEAR = 2
    integer, parameter, public :: HULL_ALL_COPLANAR = 3

contains

    !> Compute the 3-D convex hull of a point set.
    !>
    !> @param[in]  points     coordinates, shape (3, npts)
    !> @param[out] faces      triangular faces, CCW outward. Allocated here.
    !> @param[in]  tol        optional tolerance; defaults to a scale-aware
    !>                        epsilon (~1e-12 * max|coordinate|).
    !> @param[in]  randomize  optional; default .true. Randomizes the
    !>                        insertion order for expected O(n log n).
    !> @param[out] stat       optional status code (see HULL_* parameters).
    subroutine convex_hull_3d(points, faces, tol, randomize, stat)
        real(dp), intent(in)  :: points(:, :)
        type(face3d_t), allocatable, intent(out) :: faces(:)
        real(dp), optional, intent(in)  :: tol
        logical, optional, intent(in)  :: randomize
        integer, optional, intent(out) :: stat

        type(hull_state_t)   :: hs
        integer              :: npts, ip, i1, i2, i3, i4, degeneracy, i
        integer, allocatable :: order(:)
        logical              :: do_rand
        real(dp)             :: eps

        if (size(points, 1) /= 3) then
            error stop "convex_hull_3d: first extent of points must be 3"
        end if
        npts = size(points, 2)

        if (present(stat)) stat = HULL_OK

        if (npts < 4) then
            allocate (faces(0))
            if (present(stat)) stat = HULL_TOO_FEW_POINTS
            return
        end if

        do_rand = .true.
        if (present(randomize)) do_rand = randomize

        eps = default_eps(points)
        if (present(tol)) eps = tol

        call find_initial_simplex(points, eps, i1, i2, i3, i4, degeneracy)
        if (degeneracy /= HULL_OK) then
            allocate (faces(0))
            if (present(stat)) stat = degeneracy
            return
        end if

        call init_tetrahedron(hs, points, i1, i2, i3, i4)

        ! Insertion order: initial four vertices first (skipped during
        ! insertion), then the rest (optionally shuffled).
        allocate (order(npts))
        do i = 1, npts
            order(i) = i
        end do
        call swap_to_front(order, i1, 1)
        call swap_to_front(order, i2, 2)
        call swap_to_front(order, i3, 3)
        call swap_to_front(order, i4, 4)
        if (do_rand .and. npts > 5) call shuffle(order(5:))

        do ip = 5, npts
            call insert_point(hs, points, order(ip), eps)
        end do

        call compact_faces(hs, faces)
    end subroutine convex_hull_3d

    ! ==================================================================
    !  Initialization
    ! ==================================================================

    !> Scale-aware default tolerance for coplanarity/visibility tests.
    pure function default_eps(points) result(eps)
        real(dp), intent(in) :: points(:, :)
        real(dp)             :: eps
        eps = 1.0e-12_dp*max(maxval(abs(points)), 1.0_dp)
    end function default_eps

    !> Locate four affinely independent indices. Degeneracy is flagged when
    !> the input is coincident / collinear / coplanar.
    subroutine find_initial_simplex(points, eps, i1, i2, i3, i4, deg)
        real(dp), intent(in)  :: points(:, :)
        real(dp), intent(in)  :: eps
        integer, intent(out) :: i1, i2, i3, i4, deg

        integer  :: n, i, j, extremes(6)
        real(dp) :: dmax, d2, v1(3), v2(3), cr(3), nrm(3), nlen, off_, signed_d

        deg = HULL_OK
        n = size(points, 2)

        ! Extremes along the three coordinate axes: up to 6 candidates.
        extremes(1) = minloc(points(1, :), dim=1)
        extremes(2) = maxloc(points(1, :), dim=1)
        extremes(3) = minloc(points(2, :), dim=1)
        extremes(4) = maxloc(points(2, :), dim=1)
        extremes(5) = minloc(points(3, :), dim=1)
        extremes(6) = maxloc(points(3, :), dim=1)

        ! i1, i2 = farthest pair among the extremes.
        dmax = -1.0_dp
        i1 = extremes(1)
        i2 = extremes(2)
        do i = 1, 6
            do j = i + 1, 6
                d2 = sum((points(:, extremes(i)) - points(:, extremes(j)))**2)
                if (d2 > dmax) then
                    dmax = d2
                    i1 = extremes(i)
                    i2 = extremes(j)
                end if
            end do
        end do
        if (dmax <= eps*eps) then
            deg = HULL_ALL_COINCIDENT
            return
        end if

        ! i3 = point with maximum perpendicular distance to line (i1, i2).
        v1 = points(:, i2) - points(:, i1)
        dmax = -1.0_dp
        i3 = 0
        do i = 1, n
            if (i == i1 .or. i == i2) cycle
            v2 = points(:, i) - points(:, i1)
            cr = cross3d(v1, v2)
            d2 = sum(cr*cr)
            if (d2 > dmax) then
                dmax = d2
                i3 = i
            end if
        end do
        if (i3 == 0 .or. dmax <= eps*eps) then
            deg = HULL_ALL_COLLINEAR
            return
        end if

        ! i4 = point farthest from the plane of (i1, i2, i3).
        v1 = points(:, i2) - points(:, i1)
        v2 = points(:, i3) - points(:, i1)
        nrm = cross3d(v1, v2)
        nlen = norm2(nrm)
        if (nlen <= 0.0_dp) then
            deg = HULL_ALL_COLLINEAR
            return
        end if
        nrm = nrm/nlen
        off_ = dot_product(nrm, points(:, i1))

        dmax = -1.0_dp
        i4 = 0
        do i = 1, n
            if (i == i1 .or. i == i2 .or. i == i3) cycle
            signed_d = abs(dot_product(nrm, points(:, i)) - off_)
            if (signed_d > dmax) then
                dmax = signed_d
                i4 = i
            end if
        end do
        if (i4 == 0 .or. dmax <= eps) then
            deg = HULL_ALL_COPLANAR
            return
        end if
    end subroutine find_initial_simplex

    !> Build the initial four faces of a tetrahedron and wire up adjacency.
    subroutine init_tetrahedron(hs, points, i1, i2, i3, i4)
        type(hull_state_t), intent(inout) :: hs
        real(dp), intent(in)    :: points(:, :)
        integer, intent(in)    :: i1, i2, i3, i4

        integer :: inside_idx(4), k

        call ensure_capacity(hs, 16)
        hs%nf = 4

        ! Arbitrary initial winding; we flip below so all outward normals
        ! point away from the remaining tetrahedron vertex.
        hs%f(1)%v = [i1, i2, i3]; inside_idx(1) = i4
        hs%f(2)%v = [i1, i3, i4]; inside_idx(2) = i2
        hs%f(3)%v = [i1, i4, i2]; inside_idx(3) = i3
        hs%f(4)%v = [i2, i4, i3]; inside_idx(4) = i1

        do k = 1, 4
            hs%f(k)%live = .true.
            call compute_plane(hs%f(k), points)
            call orient_outward(hs%f(k), points, inside_idx(k))
        end do

        call build_tetra_adjacency(hs%f(1:4))
    end subroutine init_tetrahedron

    !> Flip vertex order (and plane) so that `inside_idx` lies on the negative
    !> side of the face.
    subroutine orient_outward(f, points, inside_idx)
        type(face3d_t), intent(inout) :: f
        real(dp), intent(in)    :: points(:, :)
        integer, intent(in)    :: inside_idx
        integer  :: tmp
        real(dp) :: d

        d = dot_product(f%n, points(:, inside_idx)) - f%off
        if (d > 0.0_dp) then
            tmp = f%v(2); f%v(2) = f%v(3); f%v(3) = tmp
            call compute_plane(f, points)
        end if
    end subroutine orient_outward

    !> Fill in adjacency pointers across the four faces of the initial
    !> tetrahedron by matching shared edges.
    subroutine build_tetra_adjacency(f)
        type(face3d_t), intent(inout) :: f(4)
        integer :: i, j, k, a, b

        do i = 1, 4
            do k = 1, 3
                a = f(i)%v(mod(k, 3) + 1)
                b = f(i)%v(mod(k + 1, 3) + 1)
                f(i)%nbr(k) = 0
                do j = 1, 4
                    if (j == i) cycle
                    if (any(f(j)%v == a) .and. any(f(j)%v == b)) then
                        f(i)%nbr(k) = j
                        exit
                    end if
                end do
            end do
        end do
    end subroutine build_tetra_adjacency

    ! ==================================================================
    !  Point insertion
    ! ==================================================================

    !> Insert a single new point into the current hull.
    subroutine insert_point(hs, points, ip, eps)
        type(hull_state_t), intent(inout) :: hs
        real(dp), intent(in)    :: points(:, :)
        integer, intent(in)    :: ip
        real(dp), intent(in)    :: eps

        integer, allocatable :: visible(:), new_faces(:)
        integer, allocatable :: horizon_face(:), horizon_a(:), horizon_b(:)
        integer  :: nvis, nhor
        integer  :: i, j, k, fi, fnew, nbr, old_count
        real(dp) :: p(3)

        p = points(:, ip)

        ! --- find visible faces -------------------------------------------
        old_count = hs%nf
        allocate (visible(old_count))
        nvis = 0
        do i = 1, old_count
            if (.not. hs%f(i)%live) cycle
            if (signed_dist(hs%f(i), p) > eps) then
                nvis = nvis + 1
                visible(nvis) = i
            end if
        end do
        if (nvis == 0) return   ! p is inside (or within tolerance of) the hull

        ! --- walk the horizon ---------------------------------------------
        allocate (horizon_face(3*nvis))
        allocate (horizon_a(3*nvis))
        allocate (horizon_b(3*nvis))
        nhor = 0
        do i = 1, nvis
            fi = visible(i)
            do k = 1, 3
                nbr = hs%f(fi)%nbr(k)
                if (nbr <= 0) cycle
                if (signed_dist(hs%f(nbr), p) > eps) cycle  ! nbr also visible
                nhor = nhor + 1
                horizon_face(nhor) = nbr
                horizon_a(nhor) = hs%f(fi)%v(mod(k, 3) + 1)
                horizon_b(nhor) = hs%f(fi)%v(mod(k + 1, 3) + 1)
            end do
        end do

        ! --- logically delete visible faces -------------------------------
        do i = 1, nvis
            hs%f(visible(i))%live = .false.
        end do

        ! --- create new fan faces -----------------------------------------
        allocate (new_faces(nhor))
        call ensure_capacity(hs, hs%nf + nhor)
        do i = 1, nhor
            hs%nf = hs%nf + 1
            fnew = hs%nf
            new_faces(i) = fnew

            hs%f(fnew)%v = [horizon_a(i), horizon_b(i), ip]
            hs%f(fnew)%live = .true.
            call compute_plane(hs%f(fnew), points)

            ! nbr(3) is opposite v(3) = p, i.e., edge (a, b): the retained
            ! neighbor across the horizon.
            hs%f(fnew)%nbr(3) = horizon_face(i)

            ! Flip the retained face's reciprocal pointer to the new face.
            call replace_neighbor(hs%f(horizon_face(i)), &
                horizon_a(i), horizon_b(i), fnew)
        end do

        ! --- link new faces to one another around the fan -----------------
        ! New face i has v = (a_i, b_i, p).
        !   nbr(1) opposite v(1)=a_i -> edge (b_i, p):
        !       shared with the new face j whose v(1) = b_i.
        !   nbr(2) opposite v(2)=b_i -> edge (p, a_i):
        !       shared with the new face j whose v(2) = a_i.
        do i = 1, nhor
            fi = new_faces(i)
            do j = 1, nhor
                if (j == i) cycle
                if (hs%f(new_faces(j))%v(1) == hs%f(fi)%v(2)) then
                    hs%f(fi)%nbr(1) = new_faces(j)
                    exit
                end if
            end do
            do j = 1, nhor
                if (j == i) cycle
                if (hs%f(new_faces(j))%v(2) == hs%f(fi)%v(1)) then
                    hs%f(fi)%nbr(2) = new_faces(j)
                    exit
                end if
            end do
        end do
    end subroutine insert_point

    !> In face f, find the slot k whose opposite edge is (a,b) or (b,a)
    !> and set f%nbr(k) = new_nbr.
    subroutine replace_neighbor(f, a, b, new_nbr)
        type(face3d_t), intent(inout) :: f
        integer, intent(in)    :: a, b, new_nbr
        integer :: k, v1, v2

        do k = 1, 3
            v1 = f%v(mod(k, 3) + 1)
            v2 = f%v(mod(k + 1, 3) + 1)
            if ((v1 == a .and. v2 == b) .or. (v1 == b .and. v2 == a)) then
                f%nbr(k) = new_nbr
                return
            end if
        end do
        error stop "replace_neighbor: shared edge not found (topology broken)"
    end subroutine replace_neighbor

    ! ==================================================================
    !  Geometry helpers
    ! ==================================================================

    pure function cross3d(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp)             :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross3d

    pure function signed_dist(f, p) result(d)
        type(face3d_t), intent(in) :: f
        real(dp), intent(in) :: p(3)
        real(dp)                   :: d
        d = dot_product(f%n, p) - f%off
    end function signed_dist

    !> Compute unit outward normal and offset for a face from its vertex
    !> indices. Assumes (v(1), v(2), v(3)) already wound CCW from outside;
    !> the caller may need to flip before calling.
    subroutine compute_plane(f, points)
        type(face3d_t), intent(inout) :: f
        real(dp), intent(in)    :: points(:, :)
        real(dp) :: a(3), b(3), n(3), nlen

        a = points(:, f%v(2)) - points(:, f%v(1))
        b = points(:, f%v(3)) - points(:, f%v(1))
        n = cross3d(a, b)
        nlen = norm2(n)
        if (nlen > 0.0_dp) then
            f%n = n/nlen
        else
            f%n = n   ! degenerate; leave caller to detect
        end if
        f%off = dot_product(f%n, points(:, f%v(1)))
    end subroutine compute_plane

    ! ==================================================================
    !  State / compaction
    ! ==================================================================

    !> Ensure `hs%f` has room for at least `min_cap` faces.
    subroutine ensure_capacity(hs, min_cap)
        type(hull_state_t), intent(inout) :: hs
        integer, intent(in)    :: min_cap
        type(face3d_t), allocatable       :: tmp(:)
        integer :: new_cap

        if (hs%cap >= min_cap) return

        if (hs%cap == 0) then
            new_cap = max(min_cap, 16)
        else
            new_cap = hs%cap
            do while (new_cap < min_cap)
                new_cap = new_cap*2
            end do
        end if

        if (.not. allocated(hs%f)) then
            allocate (hs%f(new_cap))
        else
            allocate (tmp(new_cap))
            tmp(1:hs%nf) = hs%f(1:hs%nf)
            call move_alloc(tmp, hs%f)
        end if
        hs%cap = new_cap
    end subroutine ensure_capacity

    !> Copy live faces into a fresh, tightly-sized array and remap the
    !> adjacency pointers.
    subroutine compact_faces(hs, faces)
        type(hull_state_t), intent(inout) :: hs
        type(face3d_t), allocatable, intent(out)   :: faces(:)
        integer, allocatable :: map_(:)
        integer :: i, j, nlive

        nlive = count(hs%f(1:hs%nf)%live)
        allocate (faces(nlive))
        allocate (map_(hs%nf))
        map_(:) = 0

        j = 0
        do i = 1, hs%nf
            if (hs%f(i)%live) then
                j = j + 1
                map_(i) = j
                faces(j) = hs%f(i)
            end if
        end do

        do i = 1, nlive
            do j = 1, 3
                if (faces(i)%nbr(j) > 0) then
                    faces(i)%nbr(j) = map_(faces(i)%nbr(j))
                end if
            end do
        end do
    end subroutine compact_faces

    ! ==================================================================
    !  Small utilities
    ! ==================================================================

    !> Move the entry equal to `val` to position `pos` by swapping.
    subroutine swap_to_front(order, val, pos)
        integer, intent(inout) :: order(:)
        integer, intent(in)    :: val, pos
        integer :: i, tmp

        do i = pos, size(order)
            if (order(i) == val) then
                tmp = order(pos); order(pos) = order(i); order(i) = tmp
                return
            end if
        end do
    end subroutine swap_to_front

    !> Fisher-Yates shuffle using the intrinsic random_number (Fortran's
    !> default RNG). Call random_seed() externally for reproducibility.
    subroutine shuffle(a)
        integer, intent(inout) :: a(:)
        integer  :: i, j, n, tmp
        real(dp) :: r

        n = size(a)
        do i = n, 2, -1
            call random_number(r)
            j = 1 + int(r*real(i, dp))
            if (j > i) j = i
            tmp = a(i); a(i) = a(j); a(j) = tmp
        end do
    end subroutine shuffle

end module convex_hull_3d_mod
