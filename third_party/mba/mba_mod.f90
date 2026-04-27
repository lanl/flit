!==============================================================================
! mba_mod.f90  --  Multilevel B-spline Approximation (Lee-Wolberg-Shin 1997)
!
! Scattered-data interpolation / approximation in 1-D, 2-D and 3-D using a
! hierarchy of uniform cubic-B-spline control lattices.
!
! Algorithm (per level h):
!   1. Lay out a uniform (m_h+3) [x (n_h+3)] [x (p_h+3)] control lattice
!      over the data bounding box, where m_h = m0 * 2^h.
!   2. For each data point at scaled coords (s,t[,u]) in [0,m_h] x ... :
!      - 4 (1-D), 4x4 (2-D) or 4x4x4 (3-D) cubic B-spline weights are computed.
!      - Accumulate per-control-point sums (BA averaging, Lee 1997 Eq. 4):
!            delta(c) += w(c)**3 * z_data / W
!            omega(c) += w(c)**2          where W = sum over the stencil of w**2
!   3. Compute phi(c) = delta(c) / omega(c)   (zero if omega(c) == 0).
!   4. Evaluate the level-h spline at each data point, subtract from
!      residual; the next level fits the residual ("MBA-I" refinement).
!
! Stops when residual RMS drops below tol_rms, or max_levels is reached.
!
! KD-tree integration:
!   For "extrapolation safety" the module can build an internal kdtree
!   over the data points (lazy, on first use).  The helper method
!   %nearest_data_distance(qpts, dists) returns each query's distance
!   to its nearest data point.  Users can mask unreliable regions
!   where the spline is extrapolating beyond the data's convex hull.
!
! Memory note (3-D):
!   Lattice at level h is (m0*2^h + 3)^3 doubles. With m0=4 and 8 bytes
!   per double that's level5=17 MB, level6=135 MB, level7=1 GB. The
!   build_3d routine warns if a level would exceed ~4 GB.
!
!==============================================================================

module mba_mod

    use kdtree_mod
    implicit none

    private
    public :: mba

    ! Internal precision shorthand
    integer, parameter :: DP = kind(1.0d0)
    integer, parameter :: IP = kind(0)

    !--------------------------------------------------------------------
    ! Per-level lattice storage (one of these is allocated per level)
    !--------------------------------------------------------------------
    type :: lat1d_t
        real(DP), allocatable :: phi(:)              ! phi(-1:m+1)
        integer :: m = 0
        real(DP)    :: dx = 0.0d0                   ! cell size at this level
    end type lat1d_t

    type :: lat2d_t
        real(DP), allocatable :: phi(:, :)            ! phi(-1:m+1, -1:n+1)
        integer :: m = 0, n = 0
        real(DP)    :: dx = 0.0d0, dy = 0.0d0        ! cell size at this level
    end type lat2d_t

    type :: lat3d_t
        real(DP), allocatable :: phi(:, :, :)          ! phi(-1:m+1, -1:n+1, -1:p+1)
        integer :: m = 0, n = 0, p = 0
        real(DP)    :: dx = 0.0d0, dy = 0.0d0, dz = 0.0d0
    end type lat3d_t

    !--------------------------------------------------------------------
    ! Public MBA type
    !--------------------------------------------------------------------
    type :: mba
        !! Multilevel B-spline approximation of scattered 1-D, 2-D or 3-D data.
        !!
        !! Build with %build(pts, vals [, m0, max_levels, tol_rms, verbose]).
        !! Evaluate at arbitrary query points with %eval(qpts, qvals).
        !!
        !! For extrapolation safety, %nearest_data_distance(qpts, dists)
        !! returns the distance from each query to the nearest data point,
        !! using a KD-tree built lazily on first call.
        private
        integer :: ndim = 0       ! 1, 2 or 3
        integer :: nlevels = 0       ! number of levels actually built
        integer :: m0 = 4       ! base lattice cells per axis
        real(DP)    :: bbox(2, 3) = 0.0d0  ! (min/max, dim) bounding box (with margin)

        type(lat1d_t), allocatable :: l1(:)
        type(lat2d_t), allocatable :: l2(:)
        type(lat3d_t), allocatable :: l3(:)

        ! Original data, kept for residual reporting and tree queries
        real(DP), allocatable :: data_pts(:, :)
        real(DP), allocatable :: data_val(:)

        ! Lazy KD-tree on data points (built on first nearest_data_distance call)
        type(kdtree) :: data_tree
        logical        :: tree_built = .false.

    contains
        procedure, public :: build => mba_build
        procedure, public :: eval => mba_eval
        procedure, public :: free => mba_free
        procedure, public :: built => mba_built
        procedure, public :: num_levels => mba_num_levels
        procedure, public :: num_dims => mba_num_dims
        procedure, public :: rms_residual => mba_rms_residual
        procedure, public :: max_residual => mba_max_residual
        procedure, public :: nearest_data_distance => mba_nearest_data_distance
    end type mba

    !--------------------------------------------------------------------
    ! Internal: cumulative residuals (filled during build, queried after)
    !--------------------------------------------------------------------
    ! We don't store these in mba to keep the type simple; they're
    ! computed on-demand from the stored data + spline.

contains

    !==========================================================================
    ! Cubic uniform B-spline basis functions
    !
    !   Given t in [0,1], the four overlapping basis values are:
    !     B0(t) = (1-t)^3 / 6
    !     B1(t) = (3 t^3 - 6 t^2 + 4) / 6
    !     B2(t) = (-3 t^3 + 3 t^2 + 3 t + 1) / 6
    !     B3(t) = t^3 / 6
    !==========================================================================
    pure subroutine bspline_w(t, w)
        real(DP), intent(in)  :: t
        real(DP), intent(out) :: w(0:3)
        real(DP) :: t2, t3, om
        om = 1.0d0 - t
        t2 = t*t
        t3 = t2*t
        w(0) = om*om*om/6.0d0
        w(1) = (3.0d0*t3 - 6.0d0*t2 + 4.0d0)/6.0d0
        w(2) = (-3.0d0*t3 + 3.0d0*t2 + 3.0d0*t + 1.0d0)/6.0d0
        w(3) = t3/6.0d0
    end subroutine bspline_w

    !==========================================================================
    ! Public: build
    !
    ! Dispatches on size(pts,1):  1 -> 1-D MBA, 2 -> 2-D MBA, 3 -> 3-D MBA.
    !==========================================================================
    subroutine mba_build(self, pts, vals, m0, max_levels, tol_rms, verbose)
        !! Build the multilevel B-spline.
        !!
        !! @param pts        Data point coordinates, shape (ndim, npts).
        !!                   ndim = 1, 2 or 3 (determined from size(pts,1)).
        !! @param vals       Scalar values at each data point, size npts.
        !! @param m0         (optional) Base lattice cells per axis.  Default 4.
        !!                   Must be >= 1.  Larger m0 captures higher-frequency
        !!                   structure at level 0 but costs more memory.
        !! @param max_levels (optional) Maximum number of refinement levels.
        !!                   Default 8 (1-D/2-D) or 5 (3-D).
        !! @param tol_rms    (optional) Stop refining when residual RMS at the
        !!                   data points falls below this.  Default 1.0d-8 *
        !!                   (max(vals) - min(vals)).
        !! @param verbose    (optional) If .true., print one line per level
        !!                   showing residual RMS / max.  Default .false.
        class(mba), intent(inout)        :: self
        real(DP), intent(in)           :: pts(:, :)
        real(DP), intent(in)           :: vals(:)
        integer, intent(in), optional :: m0
        integer, intent(in), optional :: max_levels
        real(DP), intent(in), optional :: tol_rms
        logical, intent(in), optional :: verbose

        integer :: nd, nlev_default
        real(DP)    :: tol
        logical     :: vb

        call self%free()

        nd = int(size(pts, 1), IP)
        if (nd /= 1 .and. nd /= 2 .and. nd /= 3) then
            write (*, '(a,i0)') '[mba] ERROR: build expects ndim = 1, 2 or 3, got ', nd
            return
        end if
        self%ndim = nd

        if (size(pts, 2) /= size(vals)) then
            write (*, '(a)') '[mba] ERROR: size(pts,2) /= size(vals)'
            return
        end if

        if (present(m0)) then
            self%m0 = max(1, m0)
        end if

        if (present(max_levels)) then
            nlev_default = max(1, max_levels)
        else
            nlev_default = merge(5, 8, nd == 3)
        end if

        if (present(tol_rms)) then
            tol = tol_rms
        else
            tol = 1.0d-8*(maxval(vals) - minval(vals))
            if (tol <= 0.0d0) then
                tol = 1.0d-12
            end if
        end if

        vb = .false.
        if (present(verbose)) then
            vb = verbose
        end if

        ! Stash data for later
        allocate (self%data_pts(nd, size(vals)))
        self%data_pts = pts(1:nd, :)
        allocate (self%data_val(size(vals)))
        self%data_val = vals

        ! Domain bounding box with a small margin so points never sit on edges
        call compute_bbox(pts, nd, self%bbox)

        if (nd == 1) then
            call build_1d_impl(self, nlev_default, tol, vb)
        else if (nd == 2) then
            call build_2d_impl(self, nlev_default, tol, vb)
        else
            call build_3d_impl(self, nlev_default, tol, vb)
        end if
    end subroutine mba_build

    !--------------------------------------------------------------------------
    pure subroutine compute_bbox(pts, nd, bbox)
        real(DP), intent(in)  :: pts(:, :)
        integer, intent(in)  :: nd
        real(DP), intent(out) :: bbox(2, 3)
        integer :: d
        real(DP)    :: lo, hi, margin
        bbox = 0.0d0
        do d = 1, nd
            lo = minval(pts(d, :))
            hi = maxval(pts(d, :))
            ! 1e-9 of range ensures floor() never returns m at the upper edge
            margin = max(1.0d-9*(hi - lo), 1.0d-15*max(abs(lo), abs(hi)))
            if (margin <= 0.0d0) then
                margin = 1.0d-15
            end if
            bbox(1, d) = lo - margin
            bbox(2, d) = hi + margin
        end do
    end subroutine compute_bbox

    !==========================================================================
    ! 1-D builder
    !==========================================================================
    subroutine build_1d_impl(self, max_lev, tol, verbose)
        type(mba), intent(inout) :: self
        integer, intent(in)    :: max_lev
        real(DP), intent(in)    :: tol
        logical, intent(in)    :: verbose

        integer :: h, np, m, lev
        real(DP)    :: rms, mx
        real(DP), allocatable :: resid(:), pred(:)

        allocate (self%l1(max_lev))
        np = size(self%data_val)
        allocate (resid(np), pred(np))
        resid = self%data_val

        if (verbose) then
            write (*, '(a)') &
                '[mba 1D] level   m              RMS resid       max resid'
        end if

        self%nlevels = 0
        do lev = 1, max_lev
            h = lev - 1
            m = self%m0*(2**h)
            call solve_level_1d(self, lev, m, resid)
            self%nlevels = lev

            ! Evaluate level lev at all data points, then fit the remaining residual
            call eval_level_1d_at_pts(self, lev, self%data_pts, pred)
            resid = resid - pred
            rms = sqrt(sum(resid**2)/real(np, DP))
            mx = maxval(abs(resid))

            if (verbose) then
                write (*, '(a,i3,2x,i6,11x,es14.6,3x,es14.6)') &
                    '         ', lev, m, rms, mx
            end if

            if (rms < tol) then
                exit
            end if
        end do

        deallocate (resid, pred)
    end subroutine build_1d_impl

    !-- Solve one 1-D level: BA averaging using current residuals as values
    subroutine solve_level_1d(self, lev, m, resid)
        type(mba), intent(inout) :: self
        integer, intent(in)    :: lev, m
        real(DP), intent(in)    :: resid(:)

        integer :: np, d, i, k
        real(DP)    :: dx, sx, s
        real(DP)    :: bx(0:3), w, W2, num
        real(DP), allocatable :: delta(:), omega(:)

        ! Allocate this level's lattice phi(-1:m+1)
        if (allocated(self%l1(lev)%phi)) then
            deallocate (self%l1(lev)%phi)
        end if
        allocate (self%l1(lev)%phi(-1:m + 1))
        self%l1(lev)%phi = 0.0d0
        self%l1(lev)%m = m

        dx = (self%bbox(2, 1) - self%bbox(1, 1))/real(m, DP)
        self%l1(lev)%dx = dx

        allocate (delta(-1:m + 1))
        delta = 0.0d0
        allocate (omega(-1:m + 1))
        omega = 0.0d0

        np = size(resid)
        do d = 1, np
            ! Scale to lattice coords [0, m]
            sx = (self%data_pts(1, d) - self%bbox(1, 1))/dx
            i = int(floor(sx), IP)
            s = sx - real(i, DP)
            call bspline_w(s, bx)

            ! Sum of squared weights over the 4-point stencil
            W2 = 0.0d0
            do k = 0, 3
                w = bx(k)
                W2 = W2 + w*w
            end do
            if (W2 <= 0.0d0) then
                cycle
            end if

            ! Accumulate BA contributions
            do k = 0, 3
                w = bx(k)
                num = w*w*w*resid(d)/W2
                delta(i - 1 + k) = delta(i - 1 + k) + num
                omega(i - 1 + k) = omega(i - 1 + k) + w*w
            end do
        end do

        ! phi = delta / omega   (zero where omega==0)
        where (omega > 0.0d0)
            self%l1(lev)%phi = delta/omega
        elsewhere
            self%l1(lev)%phi = 0.0d0
        end where

        deallocate (delta, omega)
    end subroutine solve_level_1d

    !-- Evaluate level lev's spline at given query points (1-D)
    subroutine eval_level_1d_at_pts(self, lev, qpts, qvals)
        type(mba), intent(in)    :: self
        integer, intent(in)    :: lev
        real(DP), intent(in)    :: qpts(:, :)
        real(DP), intent(out)   :: qvals(:)

        integer :: nq, d, i, k, m
        real(DP)    :: dx, sx, s
        real(DP)    :: bx(0:3), v

        m = self%l1(lev)%m
        dx = self%l1(lev)%dx
        nq = size(qpts, 2)

        do d = 1, nq
            sx = (qpts(1, d) - self%bbox(1, 1))/dx
            i = int(floor(sx), IP)
            if (i < 0) then
                i = 0
                s = 0.0d0
            else if (i >= m) then
                i = m - 1
                s = 1.0d0
            else
                s = sx - real(i, DP)
            end if
            call bspline_w(s, bx)
            v = 0.0d0
            do k = 0, 3
                v = v + bx(k)*self%l1(lev)%phi(i - 1 + k)
            end do
            qvals(d) = v
        end do
    end subroutine eval_level_1d_at_pts

    !==========================================================================
    ! 2-D builder
    !==========================================================================
    subroutine build_2d_impl(self, max_lev, tol, verbose)
        type(mba), intent(inout) :: self
        integer, intent(in)    :: max_lev
        real(DP), intent(in)    :: tol
        logical, intent(in)    :: verbose

        integer :: h, np, m, n, lev
        real(DP)    :: rms, mx
        real(DP), allocatable :: resid(:), pred(:)

        allocate (self%l2(max_lev))
        np = size(self%data_val)
        allocate (resid(np), pred(np))
        resid = self%data_val

        if (verbose) then
            write (*, '(a)') &
                '[mba 2D] level   m x n          RMS resid       max resid'
        end if

        self%nlevels = 0
        do lev = 1, max_lev
            h = lev - 1
            m = self%m0*(2**h)
            n = self%m0*(2**h)
            call solve_level_2d(self, lev, m, n, resid)
            self%nlevels = lev

            ! Evaluate level lev at all data points (cumulative add)
            call eval_level_2d_at_pts(self, lev, self%data_pts, pred)
            resid = resid - pred
            rms = sqrt(sum(resid**2)/real(np, DP))
            mx = maxval(abs(resid))

            if (verbose) then
                write (*, '(a,i3,2x,i6,a,i6,3x,es14.6,3x,es14.6)') &
                    '         ', lev, m, ' x ', n, rms, mx
            end if

            if (rms < tol) then
                exit
            end if
        end do

        deallocate (resid, pred)
    end subroutine build_2d_impl

    !-- Solve one 2-D level: BA averaging using current residuals as values
    subroutine solve_level_2d(self, lev, m, n, resid)
        type(mba), intent(inout) :: self
        integer, intent(in)    :: lev, m, n
        real(DP), intent(in)    :: resid(:)

        integer :: np, d, i, j, k, l
        real(DP)    :: dx, dy, sx, sy, s, t
        real(DP)    :: bx(0:3), by(0:3), w, W2, num
        real(DP), allocatable :: delta(:, :), omega(:, :)

        ! Allocate this level's lattice phi(-1:m+1, -1:n+1)
        if (allocated(self%l2(lev)%phi)) then
            deallocate (self%l2(lev)%phi)
        end if
        allocate (self%l2(lev)%phi(-1:m + 1, -1:n + 1))
        self%l2(lev)%phi = 0.0d0
        self%l2(lev)%m = m
        self%l2(lev)%n = n

        dx = (self%bbox(2, 1) - self%bbox(1, 1))/real(m, DP)
        dy = (self%bbox(2, 2) - self%bbox(1, 2))/real(n, DP)
        self%l2(lev)%dx = dx
        self%l2(lev)%dy = dy

        allocate (delta(-1:m + 1, -1:n + 1))
        delta = 0.0d0
        allocate (omega(-1:m + 1, -1:n + 1))
        omega = 0.0d0

        np = size(resid)
        do d = 1, np
            ! Scale to lattice coords [0, m] x [0, n]
            sx = (self%data_pts(1, d) - self%bbox(1, 1))/dx
            sy = (self%data_pts(2, d) - self%bbox(1, 2))/dy
            i = int(floor(sx), IP)
            s = sx - real(i, DP)
            j = int(floor(sy), IP)
            t = sy - real(j, DP)
            call bspline_w(s, bx)
            call bspline_w(t, by)

            ! Sum of squared weights over the 4x4 stencil
            W2 = 0.0d0
            do l = 0, 3
                do k = 0, 3
                    w = bx(k)*by(l)
                    W2 = W2 + w*w
                end do
            end do
            if (W2 <= 0.0d0) then
                cycle
            end if

            ! Accumulate BA contributions
            do l = 0, 3
                do k = 0, 3
                    w = bx(k)*by(l)
                    num = w*w*w*resid(d)/W2     ! w^3 * z / W
                    delta(i - 1 + k, j - 1 + l) = delta(i - 1 + k, j - 1 + l) + num
                    omega(i - 1 + k, j - 1 + l) = omega(i - 1 + k, j - 1 + l) + w*w
                end do
            end do
        end do

        ! phi = delta / omega   (zero where omega==0)
        where (omega > 0.0d0)
            self%l2(lev)%phi = delta/omega
        elsewhere
            self%l2(lev)%phi = 0.0d0
        end where

        deallocate (delta, omega)
    end subroutine solve_level_2d

    !-- Evaluate level lev's spline at given query points (2-D)
    subroutine eval_level_2d_at_pts(self, lev, qpts, qvals)
        type(mba), intent(in)    :: self
        integer, intent(in)    :: lev
        real(DP), intent(in)    :: qpts(:, :)
        real(DP), intent(out)   :: qvals(:)

        integer :: nq, d, i, j, k, l, m, n
        real(DP)    :: dx, dy, sx, sy, s, t
        real(DP)    :: bx(0:3), by(0:3), v

        m = self%l2(lev)%m
        n = self%l2(lev)%n
        dx = self%l2(lev)%dx
        dy = self%l2(lev)%dy
        nq = size(qpts, 2)

        !$omp parallel do private(d, sx, sy, i, j, l, k, s, t, bx, by, v)
        do d = 1, nq
            sx = (qpts(1, d) - self%bbox(1, 1))/dx
            sy = (qpts(2, d) - self%bbox(1, 2))/dy
            ! Robust integer-index clamping: if the query is on or beyond the
            ! upper boundary, fall back to the last valid cell.  Real-valued
            ! clamping (m - 1.0d-15) is fragile because eps near m=32 already
            ! exceeds 1e-15.
            i = int(floor(sx), IP)
            if (i < 0) then
                i = 0
                s = 0.0d0
            else if (i >= m) then
                i = m - 1
                s = 1.0d0
            else
                s = sx - real(i, DP)
            end if
            j = int(floor(sy), IP)
            if (j < 0) then
                j = 0
                t = 0.0d0
            else if (j >= n) then
                j = n - 1
                t = 1.0d0
            else
                t = sy - real(j, DP)
            end if
            call bspline_w(s, bx)
            call bspline_w(t, by)
            v = 0.0d0
            do l = 0, 3
                do k = 0, 3
                    v = v + bx(k)*by(l)*self%l2(lev)%phi(i - 1 + k, j - 1 + l)
                end do
            end do
            qvals(d) = v
        end do
        !$omp end parallel do
    end subroutine eval_level_2d_at_pts

    !==========================================================================
    ! 3-D builder  (mirrors 2-D, with one extra axis)
    !==========================================================================
    subroutine build_3d_impl(self, max_lev, tol, verbose)
        type(mba), intent(inout) :: self
        integer, intent(in)    :: max_lev
        real(DP), intent(in)    :: tol
        logical, intent(in)    :: verbose

        integer :: h, np, m, lev
        real(DP)    :: rms, mx, mem_gb
        real(DP), allocatable :: resid(:), pred(:)
        real(DP), parameter   :: GB = 1024.0d0**3

        allocate (self%l3(max_lev))
        np = size(self%data_val)
        allocate (resid(np), pred(np))
        resid = self%data_val

        if (verbose) then
            write (*, '(a)') &
                '[mba 3D] level     m            RMS resid       max resid'
        end if

        self%nlevels = 0
        do lev = 1, max_lev
            h = lev - 1
            m = self%m0*(2**h)

            ! Estimate this level's lattice memory (cubic)
            mem_gb = 8.0d0*real(m + 3, DP)**3/GB
            if (mem_gb > 4.0d0) then
                write (*, '(a,i0,a,f6.2,a)') &
                    '[mba 3D] WARNING: level ', lev, &
                    ' would need ~', mem_gb, ' GB of lattice memory; stopping.'
                exit
            end if

            call solve_level_3d(self, lev, m, m, m, resid)
            self%nlevels = lev

            call eval_level_3d_at_pts(self, lev, self%data_pts, pred)
            resid = resid - pred
            rms = sqrt(sum(resid**2)/real(np, DP))
            mx = maxval(abs(resid))

            if (verbose) then
                write (*, '(a,i3,4x,i6,3x,es14.6,3x,es14.6)') &
                    '         ', lev, m, rms, mx
            end if

            if (rms < tol) then
                exit
            end if
        end do

        deallocate (resid, pred)
    end subroutine build_3d_impl

    !-- Solve one 3-D level
    subroutine solve_level_3d(self, lev, m, n, p, resid)
        type(mba), intent(inout) :: self
        integer, intent(in)    :: lev, m, n, p
        real(DP), intent(in)    :: resid(:)

        integer :: np_data, d, i, j, kk, k, l, mm
        real(DP)    :: dx, dy, dz, sx, sy, sz, s, t, u
        real(DP)    :: bx(0:3), by(0:3), bz(0:3), w, W2, num
        real(DP), allocatable :: delta(:, :, :), omega(:, :, :)

        if (allocated(self%l3(lev)%phi)) then
            deallocate (self%l3(lev)%phi)
        end if
        allocate (self%l3(lev)%phi(-1:m + 1, -1:n + 1, -1:p + 1))
        self%l3(lev)%phi = 0.0d0
        self%l3(lev)%m = m
        self%l3(lev)%n = n
        self%l3(lev)%p = p

        dx = (self%bbox(2, 1) - self%bbox(1, 1))/real(m, DP)
        dy = (self%bbox(2, 2) - self%bbox(1, 2))/real(n, DP)
        dz = (self%bbox(2, 3) - self%bbox(1, 3))/real(p, DP)
        self%l3(lev)%dx = dx
        self%l3(lev)%dy = dy
        self%l3(lev)%dz = dz

        allocate (delta(-1:m + 1, -1:n + 1, -1:p + 1))
        delta = 0.0d0
        allocate (omega(-1:m + 1, -1:n + 1, -1:p + 1))
        omega = 0.0d0

        np_data = size(resid)
        do d = 1, np_data
            sx = (self%data_pts(1, d) - self%bbox(1, 1))/dx
            sy = (self%data_pts(2, d) - self%bbox(1, 2))/dy
            sz = (self%data_pts(3, d) - self%bbox(1, 3))/dz
            i = int(floor(sx), IP)
            s = sx - real(i, DP)
            j = int(floor(sy), IP)
            t = sy - real(j, DP)
            kk = int(floor(sz), IP)
            u = sz - real(kk, DP)
            call bspline_w(s, bx)
            call bspline_w(t, by)
            call bspline_w(u, bz)

            W2 = 0.0d0
            do mm = 0, 3
                do l = 0, 3
                    do k = 0, 3
                        w = bx(k)*by(l)*bz(mm)
                        W2 = W2 + w*w
                    end do
                end do
            end do
            if (W2 <= 0.0d0) then
                cycle
            end if

            do mm = 0, 3
                do l = 0, 3
                    do k = 0, 3
                        w = bx(k)*by(l)*bz(mm)
                        num = w*w*w*resid(d)/W2
                        delta(i - 1 + k, j - 1 + l, kk - 1 + mm) = delta(i - 1 + k, j - 1 + l, kk - 1 + mm) + num
                        omega(i - 1 + k, j - 1 + l, kk - 1 + mm) = omega(i - 1 + k, j - 1 + l, kk - 1 + mm) + w*w
                    end do
                end do
            end do
        end do

        where (omega > 0.0d0)
            self%l3(lev)%phi = delta/omega
        elsewhere
            self%l3(lev)%phi = 0.0d0
        end where

        deallocate (delta, omega)
    end subroutine solve_level_3d

    !-- Evaluate level lev's spline at given query points (3-D)
    subroutine eval_level_3d_at_pts(self, lev, qpts, qvals)
        type(mba), intent(in)  :: self
        integer, intent(in)  :: lev
        real(DP), intent(in)  :: qpts(:, :)
        real(DP), intent(out) :: qvals(:)

        integer :: nq, d, i, j, kk, k, l, mm, m, n, p
        real(DP)    :: dx, dy, dz, sx, sy, sz, s, t, u
        real(DP)    :: bx(0:3), by(0:3), bz(0:3), v

        m = self%l3(lev)%m
        n = self%l3(lev)%n
        p = self%l3(lev)%p
        dx = self%l3(lev)%dx
        dy = self%l3(lev)%dy
        dz = self%l3(lev)%dz
        nq = size(qpts, 2)

        !$omp parallel do private(d, sx, sy, sz, i, j, kk, l, k, mm, s, t, u, bx, by, bz, v)
        do d = 1, nq
            sx = (qpts(1, d) - self%bbox(1, 1))/dx
            sy = (qpts(2, d) - self%bbox(1, 2))/dy
            sz = (qpts(3, d) - self%bbox(1, 3))/dz
            i = int(floor(sx), IP)
            if (i < 0) then
                i = 0
                s = 0.0d0
            else if (i >= m) then
                i = m - 1
                s = 1.0d0
            else
                s = sx - real(i, DP)
            end if
            j = int(floor(sy), IP)
            if (j < 0) then
                j = 0
                t = 0.0d0
            else if (j >= n) then
                j = n - 1
                t = 1.0d0
            else
                t = sy - real(j, DP)
            end if
            kk = int(floor(sz), IP)
            if (kk < 0) then
                kk = 0
                u = 0.0d0
            else if (kk >= p) then
                kk = p - 1
                u = 1.0d0
            else
                u = sz - real(kk, DP)
            end if
            call bspline_w(s, bx)
            call bspline_w(t, by)
            call bspline_w(u, bz)
            v = 0.0d0
            do mm = 0, 3
                do l = 0, 3
                    do k = 0, 3
                        v = v + bx(k)*by(l)*bz(mm)*self%l3(lev)%phi(i - 1 + k, j - 1 + l, kk - 1 + mm)
                    end do
                end do
            end do
            qvals(d) = v
        end do
        !$omp end parallel do
    end subroutine eval_level_3d_at_pts

    !==========================================================================
    ! Public eval: cumulative sum across all built levels
    !==========================================================================
    subroutine mba_eval(self, qpts, qvals)
        !! Evaluate the full multilevel spline at arbitrary query points.
        !!
        !! @param qpts  Query points, shape (ndim, nq).  Must match build's ndim.
        !! @param qvals Output values, size nq.
        !!
        !! For points inside (or near) the data bounding box this returns the
        !! MBA approximation (essentially interpolating at data points).  For
        !! points well outside the data domain the spline extrapolates and
        !! results may be unreliable -- use %nearest_data_distance to mask them.
        class(mba), intent(inout) :: self
        real(DP), intent(in)    :: qpts(:, :)
        real(DP), intent(out)   :: qvals(:)
        integer :: lev, nq
        real(DP), allocatable :: tmp(:)

        if (.not. self%built()) then
            write (*, '(a)') '[mba] ERROR: eval called on unbuilt object'
            qvals = 0.0d0
            return
        end if
        if (int(size(qpts, 1), IP) /= self%ndim) then
            write (*, '(a,i0,a,i0)') '[mba] ERROR: query ndim=', size(qpts, 1), &
                ' but built ndim=', self%ndim
            qvals = 0.0d0
            return
        end if

        nq = int(size(qpts, 2), IP)
        allocate (tmp(nq))
        qvals = 0.0d0
        if (self%ndim == 1) then
            do lev = 1, self%nlevels
                call eval_level_1d_at_pts(self, lev, qpts, tmp)
                qvals = qvals + tmp
            end do
        else if (self%ndim == 2) then
            do lev = 1, self%nlevels
                call eval_level_2d_at_pts(self, lev, qpts, tmp)
                qvals = qvals + tmp
            end do
        else
            do lev = 1, self%nlevels
                call eval_level_3d_at_pts(self, lev, qpts, tmp)
                qvals = qvals + tmp
            end do
        end if
        deallocate (tmp)
    end subroutine mba_eval

    !==========================================================================
    ! KD-tree integration: nearest-data-point distance
    !
    ! Builds the tree lazily on first call, then reuses it.  Useful for
    ! masking query points that fall outside the data's effective domain
    ! where the spline extrapolates and results are unreliable.
    !==========================================================================
    subroutine mba_nearest_data_distance(self, qpts, dists)
        !! Distance from each query point to the nearest data point.
        !!
        !! @param qpts   Query points, shape (ndim, nq).
        !! @param dists  Output distances, size nq.
        class(mba), intent(inout) :: self
        real(DP), intent(in)    :: qpts(:, :)
        real(DP), intent(out)   :: dists(:)
        integer :: i, nq, idx

        if (.not. self%built()) then
            write (*, '(a)') '[mba] ERROR: nearest_data_distance on unbuilt object'
            dists = 0.0d0
            return
        end if

        if (.not. self%tree_built) then
            call self%data_tree%build(self%data_pts)
            self%tree_built = .true.
        end if

        nq = int(size(qpts, 2), IP)
        do i = 1, nq
            call self%data_tree%query_nn(qpts(:, i), idx, dists(i))
        end do
    end subroutine mba_nearest_data_distance

    !==========================================================================
    ! Status accessors
    !==========================================================================
    pure logical function mba_built(self)
        class(mba), intent(in) :: self
        mba_built = (self%nlevels > 0)
    end function mba_built

    pure integer function mba_num_levels(self)
        class(mba), intent(in) :: self
        mba_num_levels = self%nlevels
    end function mba_num_levels

    pure integer function mba_num_dims(self)
        class(mba), intent(in) :: self
        mba_num_dims = self%ndim
    end function mba_num_dims

    !-- Residual statistics at the original data points (final fit quality)
    function mba_rms_residual(self) result(rms)
        class(mba), intent(inout) :: self
        real(DP) :: rms
        real(DP), allocatable :: pred(:), r(:)
        integer :: np
        if (.not. self%built()) then
            rms = 0.0d0
            return
        end if
        np = size(self%data_val)
        allocate (pred(np), r(np))
        call self%eval(self%data_pts, pred)
        r = self%data_val - pred
        rms = sqrt(sum(r**2)/real(np, DP))
        deallocate (pred, r)
    end function mba_rms_residual

    function mba_max_residual(self) result(mx)
        class(mba), intent(inout) :: self
        real(DP) :: mx
        real(DP), allocatable :: pred(:)
        integer :: np
        if (.not. self%built()) then
            mx = 0.0d0
            return
        end if
        np = size(self%data_val)
        allocate (pred(np))
        call self%eval(self%data_pts, pred)
        mx = maxval(abs(self%data_val - pred))
        deallocate (pred)
    end function mba_max_residual

    !==========================================================================
    ! Free
    !==========================================================================
    subroutine mba_free(self)
        !! Release all memory and reset to uninitialised state.
        class(mba), intent(inout) :: self
        integer :: lev
        if (allocated(self%l1)) then
            do lev = 1, size(self%l1)
                if (allocated(self%l1(lev)%phi)) then
                    deallocate (self%l1(lev)%phi)
                end if
            end do
            deallocate (self%l1)
        end if
        if (allocated(self%l2)) then
            do lev = 1, size(self%l2)
                if (allocated(self%l2(lev)%phi)) then
                    deallocate (self%l2(lev)%phi)
                end if
            end do
            deallocate (self%l2)
        end if
        if (allocated(self%l3)) then
            do lev = 1, size(self%l3)
                if (allocated(self%l3(lev)%phi)) then
                    deallocate (self%l3(lev)%phi)
                end if
            end do
            deallocate (self%l3)
        end if
        if (allocated(self%data_pts)) then
            deallocate (self%data_pts)
        end if
        if (allocated(self%data_val)) then
            deallocate (self%data_val)
        end if
        if (self%tree_built) then
            call self%data_tree%free()
            self%tree_built = .false.
        end if
        self%ndim = 0
        self%nlevels = 0
        self%m0 = 4
        self%bbox = 0.0d0
    end subroutine mba_free

end module mba_mod
