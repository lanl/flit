!>
!> Two-dimensional convex hull via Andrew's monotone chain.
!>
!> Algorithm:
!>  1. Sort points lexicographically by (x, y).
!>  2. Scan left-to-right building the lower hull; pop from the stack as long
!>     as the last three points do not make a strict left turn.
!>  3. Scan right-to-left building the upper hull likewise.
!>  4. Concatenate (dropping the duplicated seam).
!>
!> Complexity: O(n log n) time (dominated by the sort), O(n) extra memory.
!>

module convex_hull_2d_mod

    use convex_hull_kinds, only: dp
    implicit none
    private
    public :: convex_hull_2d

contains

    !> Compute the 2-D convex hull of a point set.
    !>
    !> @param[in]  points  Point coordinates, shape (2, npts).
    !> @param[out] hull    Allocated on return; 1-based vertex indices in
    !>                     counter-clockwise order. Length is the number of
    !>                     hull vertices. The polygon is *not* closed
    !>                     (hull(1) is not repeated at the end).
    !> @param[in]  tol     Optional collinearity tolerance. A point is treated
    !>                     as "not strictly convex" when the signed-area test
    !>                     is <= tol. Default is 0 (strict convex hull: no
    !>                     interior collinear points are reported).
    subroutine convex_hull_2d(points, hull, tol)
        real(dp), intent(in)  :: points(:, :)
        integer, allocatable, intent(out) :: hull(:)
        real(dp), optional, intent(in)  :: tol

        integer              :: npts, i, k, lower_size
        integer, allocatable :: idx(:), stk(:)
        real(dp)             :: eps

        if (size(points, 1) /= 2) then
            error stop "convex_hull_2d: first extent of points must be 2"
        end if
        npts = size(points, 2)

        eps = 0.0_dp
        if (present(tol)) eps = tol

        ! Degenerate sizes: return the points in input order.
        if (npts <= 2) then
            allocate (hull(npts))
            do i = 1, npts
                hull(i) = i
            end do
            return
        end if

        allocate (idx(npts))
        do i = 1, npts
            idx(i) = i
        end do
        call sort_indices_lex(points, idx)

        allocate (stk(2*npts))
        k = 0

        ! Lower hull.
        do i = 1, npts
            do while (k >= 2)
                if (cross_z(points(:, stk(k - 1)), &
                    points(:, stk(k)), &
                    points(:, idx(i))) > eps) exit
                k = k - 1
            end do
            k = k + 1
            stk(k) = idx(i)
        end do
        lower_size = k + 1

        ! Upper hull (skip the rightmost point: already on the lower hull).
        do i = npts - 1, 1, -1
            do while (k >= lower_size)
                if (cross_z(points(:, stk(k - 1)), &
                    points(:, stk(k)), &
                    points(:, idx(i))) > eps) exit
                k = k - 1
            end do
            k = k + 1
            stk(k) = idx(i)
        end do

        ! The last vertex on the upper hull is the leftmost point, which also
        ! starts the lower hull -> drop the duplicated closing vertex.
        allocate (hull(k - 1))
        hull(:) = stk(1:k - 1)
    end subroutine convex_hull_2d

    ! ------------------------------------------------------------------
    !  Internal helpers
    ! ------------------------------------------------------------------

    !> z-component of (a-o) x (b-o). Positive = CCW turn.
    pure function cross_z(o, a, b) result(z)
        real(dp), intent(in) :: o(2), a(2), b(2)
        real(dp)             :: z
        z = (a(1) - o(1))*(b(2) - o(2)) &
            - (a(2) - o(2))*(b(1) - o(1))
    end function cross_z

    !> Lexicographic (x, y) comparison.
    pure function lex_less(p, q) result(r)
        real(dp), intent(in) :: p(2), q(2)
        logical              :: r
        if (p(1) /= q(1)) then
            r = p(1) < q(1)
        else
            r = p(2) < q(2)
        end if
    end function lex_less

    !> Sort indices so that points(:, idx(:)) is lexicographically ordered.
    !> Uses recursive merge sort: stable, guaranteed O(n log n).
    subroutine sort_indices_lex(points, idx)
        real(dp), intent(in)    :: points(:, :)
        integer, intent(inout) :: idx(:)
        integer, allocatable    :: work(:)

        allocate (work(size(idx)))
        call merge_sort(idx, work, 1, size(idx), points)
    end subroutine sort_indices_lex

    recursive subroutine merge_sort(a, w, lo, hi, pts)
        integer, intent(inout) :: a(:), w(:)
        integer, intent(in)    :: lo, hi
        real(dp), intent(in)    :: pts(:, :)
        integer                 :: mid

        if (hi <= lo) return
        mid = (lo + hi)/2
        call merge_sort(a, w, lo, mid, pts)
        call merge_sort(a, w, mid + 1, hi, pts)
        call merge_(a, w, lo, mid, hi, pts)
    end subroutine merge_sort

    subroutine merge_(a, w, lo, mid, hi, pts)
        integer, intent(inout) :: a(:), w(:)
        integer, intent(in)    :: lo, mid, hi
        real(dp), intent(in)    :: pts(:, :)
        integer                 :: i, j, k

        i = lo
        j = mid + 1
        k = lo
        do while (i <= mid .and. j <= hi)
            if (lex_less(pts(:, a(i)), pts(:, a(j)))) then
                w(k) = a(i); i = i + 1
            else
                w(k) = a(j); j = j + 1
            end if
            k = k + 1
        end do
        do while (i <= mid)
            w(k) = a(i); i = i + 1; k = k + 1
        end do
        do while (j <= hi)
            w(k) = a(j); j = j + 1; k = k + 1
        end do
        a(lo:hi) = w(lo:hi)
    end subroutine merge_

end module convex_hull_2d_mod
