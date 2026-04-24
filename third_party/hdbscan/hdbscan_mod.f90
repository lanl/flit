!==============================================================================
! MODULE: hdbscan_mod
!
! Hierarchical Density-Based Spatial Clustering of Applications with Noise
! (HDBSCAN) — Modern Fortran 2018 implementation
!
! Distance metrics supported (params%metric):
!   'euclidean'  — L2 distance (default)
!   'manhattan'  — L1 / taxicab distance
!   'chebyshev'  — L-inf distance
!   'minkowski'  — Lp distance (set params%minkowski_p; p=2 => euclidean)
!   'cosine'     — 1 - cos(angle) in [0,1]
!
! Reference:
!   Campello, R.J.G.B., Moulavi, D., Sander, J. (2013).
!   "Density-Based Clustering Based on Hierarchical Density Estimates."
!   PAKDD 2013, LNAI 7819, pp. 160-172.
!
! Standard: Fortran 2018
!==============================================================================

module hdbscan_mod
    use iso_fortran_env, only: rk => real64, ik => int32
    implicit none
    private

    public :: hdbscan_params, hdbscan_result
    public :: hdbscan_fit, hdbscan_free
    public :: HDBSCAN_NOISE

    integer(ik), parameter :: HDBSCAN_NOISE = -1_ik

    !----------------------------------------------------------------------------
    ! hdbscan_params
    !----------------------------------------------------------------------------
    type :: hdbscan_params
        integer(ik)       :: min_cluster_size          = 5_ik
        integer(ik)       :: min_samples               = 5_ik
        real(rk)          :: cluster_selection_epsilon = 0.0_rk
        logical           :: allow_single_cluster      = .false.
        character(len=16) :: metric                    = 'euclidean'
        real(rk)          :: minkowski_p               = 2.0_rk  ! used only for 'minkowski'
    end type hdbscan_params

    !----------------------------------------------------------------------------
    ! hdbscan_result
    !----------------------------------------------------------------------------
    type :: hdbscan_result
        integer(ik), allocatable :: labels(:)
        real(rk),    allocatable :: probabilities(:)
        real(rk),    allocatable :: outlier_scores(:)
        real(rk),    allocatable :: core_distances(:)
        integer(ik)              :: n_clusters = 0_ik
        integer(ik)              :: n_noise    = 0_ik
    end type hdbscan_result

    !-- Private derived types ---------------------------------------------------
    type :: mst_edge
        integer(ik) :: u = 0, v = 0
        real(rk)    :: weight = 0.0_rk
    end type mst_edge

    type :: dendro_node
        integer(ik) :: left_child  = 0
        integer(ik) :: right_child = 0
        real(rk)    :: lambda_val  = 0.0_rk
        integer(ik) :: sz          = 0
    end type dendro_node

    type :: ct_node
        integer(ik) :: parent    = 0
        integer(ik) :: child1    = 0
        integer(ik) :: child2    = 0
        real(rk)    :: lam_birth = 0.0_rk
        real(rk)    :: stability = 0.0_rk
        integer(ik) :: n_born    = 0
    end type ct_node

contains

    !=============================================================================
    ! Distance functions
    !=============================================================================

    !-- Dispatch: pairwise distance matrix for any supported metric
    subroutine compute_dist_matrix(pts, n, d, metric, mink_p, dmat)
        integer(ik),      intent(in)  :: n, d
        real(rk),         intent(in)  :: pts(d, n)
        character(len=*), intent(in)  :: metric
        real(rk),         intent(in)  :: mink_p
        real(rk),         intent(out) :: dmat(n, n)

        character(len=16) :: m
        integer(ik) :: i, j
        real(rk)    :: s

        m = trim(adjustl(metric))

        dmat = 0.0_rk

        select case (m)

            case ('euclidean')
                do j = 1, n
                    do i = j+1, n
                        s = sum((pts(1:d,i) - pts(1:d,j))**2)
                        dmat(i,j) = sqrt(s)
                        dmat(j,i) = dmat(i,j)
                    end do
                end do

            case ('manhattan')
                do j = 1, n
                    do i = j+1, n
                        s = sum(abs(pts(1:d,i) - pts(1:d,j)))
                        dmat(i,j) = s
                        dmat(j,i) = s
                    end do
                end do

            case ('chebyshev')
                do j = 1, n
                    do i = j+1, n
                        s = maxval(abs(pts(1:d,i) - pts(1:d,j)))
                        dmat(i,j) = s
                        dmat(j,i) = s
                    end do
                end do

            case ('minkowski')
                do j = 1, n
                    do i = j+1, n
                        s = sum(abs(pts(1:d,i) - pts(1:d,j))**mink_p)**(1.0_rk/mink_p)
                        dmat(i,j) = s
                        dmat(j,i) = s
                    end do
                end do

            case ('cosine')
                ! cosine distance = 1 - cos_similarity in [0,1]
                do j = 1, n
                    do i = j+1, n
                        block
                            real(rk) :: dot, na, nb, sim
                            dot = sum(pts(1:d,i) * pts(1:d,j))
                            na  = sqrt(sum(pts(1:d,i)**2))
                            nb  = sqrt(sum(pts(1:d,j)**2))
                            if (na > 0.0_rk .and. nb > 0.0_rk) then
                                sim = dot / (na * nb)
                                sim = min(1.0_rk, max(-1.0_rk, sim))
                                s = 1.0_rk - sim
                            else
                                s = 1.0_rk   ! undefined direction → max distance
                            end if
                            dmat(i,j) = s
                            dmat(j,i) = s
                        end block
                    end do
                end do

            case default
                ! Unknown metric: fall back to euclidean with a warning to stderr
                write(*, '(a,a,a)') 'hdbscan_mod WARNING: unknown metric "', &
                    trim(m), '"; falling back to euclidean'
                do j = 1, n
                    do i = j+1, n
                        s = sum((pts(1:d,i) - pts(1:d,j))**2)
                        dmat(i,j) = sqrt(s)
                        dmat(j,i) = dmat(i,j)
                    end do
                end do

        end select
    end subroutine compute_dist_matrix

    !=============================================================================
    ! PUBLIC: hdbscan_fit
    !   pts    — (n_dims, n_points)  column-major
    !   params — algorithm parameters (including metric)
    !   res    — allocated here; caller must call hdbscan_free
    !=============================================================================
    subroutine hdbscan_fit(pts, params, res)
        real(rk),             intent(in)  :: pts(:,:)
        type(hdbscan_params), intent(in)  :: params
        type(hdbscan_result), intent(out) :: res

        integer(ik) :: n_pts, n_dim, k, n_cl
        real(rk),          allocatable :: dist_mat(:,:)
        real(rk),          allocatable :: core_dist(:)
        real(rk),          allocatable :: mrd_mat(:,:)
        type(mst_edge),    allocatable :: mst_arr(:)
        type(dendro_node), allocatable :: dg(:)
        type(ct_node),     allocatable :: ct(:)
        integer(ik),       allocatable :: pt_cl(:)
        real(rk),          allocatable :: pt_lam(:)

        n_pts = size(pts, 2)
        n_dim = size(pts, 1)
        k     = max(1_ik, params%min_samples)

        call init_result(res, n_pts)
        if (n_pts < 2) then
            res%labels  = HDBSCAN_NOISE
            res%n_noise = n_pts
            return
        end if

        allocate(dist_mat(n_pts, n_pts))
        call compute_dist_matrix(pts, n_pts, n_dim, params%metric, &
            params%minkowski_p, dist_mat)

        allocate(core_dist(n_pts))
        call compute_core_dists(dist_mat, n_pts, k, core_dist)
        res%core_distances = core_dist

        allocate(mrd_mat(n_pts, n_pts))
        call compute_mrd(dist_mat, core_dist, n_pts, mrd_mat)
        deallocate(dist_mat)

        allocate(mst_arr(n_pts - 1))
        call prim_mst(mrd_mat, n_pts, mst_arr)
        deallocate(mrd_mat)

        call sort_edges(mst_arr, n_pts - 1)

        allocate(dg(n_pts - 1))
        call build_dendrogram(mst_arr, n_pts, dg)

        allocate(ct(2*n_pts))
        allocate(pt_cl(n_pts),  source=0_ik)
        allocate(pt_lam(n_pts), source=0.0_rk)
        call extract_condensed_tree(dg, n_pts, params%min_cluster_size, &
            n_cl, ct, pt_cl, pt_lam)

        block
            logical, allocatable :: is_sel(:)
            allocate(is_sel(n_cl), source=.false.)
            call eom_select(ct, n_cl, params, is_sel)
            call assign_labels(pt_cl, pt_lam, ct, is_sel, n_pts, n_cl, res)
        end block
        call compute_glosh(res, n_pts)

        res%n_clusters = maxval(res%labels, mask=(res%labels > 0))
        if (res%n_clusters < 0) then
            res%n_clusters = 0
        end if
        res%n_noise = count(res%labels == HDBSCAN_NOISE)

    end subroutine hdbscan_fit

    !=============================================================================
    subroutine hdbscan_free(res)
        type(hdbscan_result), intent(inout) :: res
        if (allocated(res%labels)) then
            deallocate(res%labels)
        end if
        if (allocated(res%probabilities)) then
            deallocate(res%probabilities)
        end if
        if (allocated(res%outlier_scores)) then
            deallocate(res%outlier_scores)
        end if
        if (allocated(res%core_distances)) then
            deallocate(res%core_distances)
        end if
        res%n_clusters = 0
        res%n_noise = 0
    end subroutine hdbscan_free

    !=============================================================================
    ! PRIVATE ROUTINES
    !=============================================================================

    subroutine init_result(res, n)
        type(hdbscan_result), intent(out) :: res
        integer(ik),          intent(in)  :: n
        allocate(res%labels(n),         source=HDBSCAN_NOISE)
        allocate(res%probabilities(n),  source=0.0_rk)
        allocate(res%outlier_scores(n), source=0.0_rk)
        allocate(res%core_distances(n), source=0.0_rk)
    end subroutine init_result

    subroutine compute_core_dists(dmat, n, k, cd)
        integer(ik), intent(in)  :: n, k
        real(rk),    intent(in)  :: dmat(n, n)
        real(rk),    intent(out) :: cd(n)
        integer(ik) :: i, j, ki, idx
        real(rk)    :: row(n), tmp

        do i = 1, n
            row(1:n) = dmat(i, 1:n)
            row(i)   = huge(1.0_rk)
            do ki = 1, k
                idx = ki
                do j = ki+1, n
                    if (row(j) < row(idx)) then
                        idx = j
                    end if
                end do
                tmp = row(ki)
                row(ki) = row(idx)
                row(idx) = tmp
            end do
            cd(i) = row(k)
        end do
    end subroutine compute_core_dists

    subroutine compute_mrd(dmat, cd, n, mrd_out)
        integer(ik), intent(in)  :: n
        real(rk),    intent(in)  :: dmat(n, n), cd(n)
        real(rk),    intent(out) :: mrd_out(n, n)
        integer(ik) :: i, j

        do j = 1, n
            do i = 1, n
                if (i == j) then
                    mrd_out(i,j) = 0.0_rk
                else
                    mrd_out(i,j) = max(cd(i), cd(j), dmat(i,j))
                end if
            end do
        end do
    end subroutine compute_mrd

    subroutine prim_mst(mrd_in, n, mst_out)
        integer(ik),    intent(in)  :: n
        real(rk),       intent(in)  :: mrd_in(n, n)
        type(mst_edge), intent(out) :: mst_out(n-1)

        real(rk),    allocatable :: key(:)
        integer(ik), allocatable :: par(:)
        logical,     allocatable :: intree(:)
        integer(ik) :: e, i, u, v
        real(rk)    :: mn

        allocate(key(n),    source=huge(1.0_rk))
        allocate(par(n),    source=0_ik)
        allocate(intree(n), source=.false.)

        intree(1) = .true.
        do v = 2, n
            key(v) = mrd_in(1, v)
            par(v) = 1
        end do

        do e = 1, n - 1
            mn = huge(1.0_rk)
            u = 0
            do i = 1, n
                if ((.not. intree(i)) .and. key(i) < mn) then
                    mn = key(i)
                    u = i
                end if
            end do
            if (u == 0) then
                exit
            end if

            intree(u) = .true.
            mst_out(e)%u      = par(u)
            mst_out(e)%v      = u
            mst_out(e)%weight = mrd_in(par(u), u)

            do v = 1, n
                if ((.not. intree(v)) .and. mrd_in(u, v) < key(v)) then
                    key(v) = mrd_in(u, v)
                    par(v) = u
                end if
            end do
        end do
    end subroutine prim_mst

    subroutine sort_edges(edges, m)
        integer(ik),    intent(in)    :: m
        type(mst_edge), intent(inout) :: edges(m)
        integer(ik)    :: i, j
        type(mst_edge) :: tmp

        do i = 2, m
            tmp = edges(i)
            j = i - 1
            do while (j >= 1)
                if (edges(j)%weight > tmp%weight) then
                    edges(j+1) = edges(j)
                    j = j - 1
                else
                    exit
                end if
            end do
            edges(j+1) = tmp
        end do
    end subroutine sort_edges

    pure function uf_find(par, x) result(root)
        integer(ik), intent(in) :: par(:), x
        integer(ik) :: root
        root = x
        do while (par(root) /= root)
            root = par(root)
        end do
    end function uf_find

    subroutine uf_union(par, rnk, x, y)
        integer(ik), intent(inout) :: par(:), rnk(:)
        integer(ik), intent(in)    :: x, y
        integer(ik) :: rx, ry
        rx = uf_find(par, x)
        ry = uf_find(par, y)
        if (rx == ry) then
            return
        end if
        if (rnk(rx) < rnk(ry)) then
            par(rx) = ry
        else if (rnk(rx) > rnk(ry)) then
            par(ry) = rx
        else
            par(ry) = rx
            rnk(rx) = rnk(rx) + 1
        end if
    end subroutine uf_union

    subroutine build_dendrogram(edges, n, dg)
        integer(ik),       intent(in)  :: n
        type(mst_edge),    intent(in)  :: edges(n-1)
        type(dendro_node), intent(out) :: dg(n-1)

        integer(ik), allocatable :: par(:), rnk(:), sub_sz(:)
        integer(ik) :: i, ru, rv, nid

        allocate(par(2*n),    source=0_ik)
        allocate(rnk(2*n),    source=0_ik)
        allocate(sub_sz(2*n), source=0_ik)
        do i = 1, 2*n
            par(i) = i
        end do
        sub_sz(1:n) = 1

        nid = n
        do i = 1, n-1
            ru = uf_find(par, edges(i)%u)
            rv = uf_find(par, edges(i)%v)
            nid = nid + 1

            dg(i)%left_child  = ru
            dg(i)%right_child = rv
            dg(i)%sz          = sub_sz(ru) + sub_sz(rv)

            if (edges(i)%weight > 0.0_rk) then
                dg(i)%lambda_val = 1.0_rk / edges(i)%weight
            else
                dg(i)%lambda_val = huge(1.0_rk)
            end if

            sub_sz(nid) = dg(i)%sz
            call uf_union(par, rnk, ru, rv)
            par(ru) = nid
            par(rv) = nid
        end do
    end subroutine build_dendrogram

    pure function node_size(nd, n, dg) result(sz)
        integer(ik),       intent(in) :: nd, n
        type(dendro_node), intent(in) :: dg(:)
        integer(ik) :: sz
        if (nd <= n) then
            sz = 1
        else
            sz = dg(nd - n)%sz
        end if
    end function node_size

    subroutine mark_subtree(nd_in, cl, lam_in, n, dg, pt_cl, pt_lam)
        integer(ik),       intent(in)    :: nd_in, cl, n
        real(rk),          intent(in)    :: lam_in
        type(dendro_node), intent(in)    :: dg(:)
        integer(ik),       intent(inout) :: pt_cl(:)
        real(rk),          intent(inout) :: pt_lam(:)

        integer(ik) :: stk(2*n), top, nd

        stk(1) = nd_in
        top = 1
        do while (top > 0)
            nd = stk(top)
            top = top - 1
            if (nd <= n) then
                pt_cl(nd)  = cl
                pt_lam(nd) = lam_in
            else
                top = top + 1
                stk(top) = dg(nd-n)%left_child
                top = top + 1
                stk(top) = dg(nd-n)%right_child
            end if
        end do
    end subroutine mark_subtree

    subroutine extract_condensed_tree(dg, n, min_cl_sz, n_cl, ct, pt_cl, pt_lam)
        integer(ik),       intent(in)    :: n, min_cl_sz
        type(dendro_node), intent(in)    :: dg(n-1)
        integer(ik),       intent(out)   :: n_cl
        type(ct_node),     intent(inout) :: ct(2*n)
        integer(ik),       intent(inout) :: pt_cl(n)
        real(rk),          intent(inout) :: pt_lam(n)

        integer(ik) :: stk_nd(4*n), stk_cl(4*n), stk_top
        real(rk)    :: stk_lb(4*n)
        integer(ik) :: nd, cl, i, lc, rc, sz_lc, sz_rc
        real(rk)    :: lam_b, lam_here

        ct    = ct_node()
        n_cl  = 1
        ct(1)%parent    = 0
        ct(1)%lam_birth = 0.0_rk
        ct(1)%n_born    = n

        stk_nd(1) = 2*n - 1
        stk_cl(1) = 1
        stk_lb(1) = 0.0_rk
        stk_top   = 1

        do while (stk_top > 0)
            nd    = stk_nd(stk_top)
            cl    = stk_cl(stk_top)
            lam_b = stk_lb(stk_top)
            stk_top = stk_top - 1

            if (nd <= n) then
                if (pt_cl(nd) == 0) then
                    pt_cl(nd)  = cl
                    pt_lam(nd) = lam_b
                end if
                cycle
            end if

            i        = nd - n
            lc       = dg(i)%left_child
            rc       = dg(i)%right_child
            sz_lc    = node_size(lc, n, dg)
            sz_rc    = node_size(rc, n, dg)
            lam_here = dg(i)%lambda_val

            if (sz_lc >= min_cl_sz .and. sz_rc >= min_cl_sz) then
                ct(cl)%stability = ct(cl)%stability + &
                    real(sz_lc + sz_rc, rk) * (lam_here - lam_b)

                n_cl = n_cl + 1
                ct(n_cl)%parent    = cl
                ct(n_cl)%lam_birth = lam_here
                ct(n_cl)%n_born    = sz_lc
                ct(cl)%child1      = n_cl
                if (lc > n) then
                    stk_top = stk_top + 1
                    stk_nd(stk_top) = lc
                    stk_cl(stk_top) = n_cl
                    stk_lb(stk_top) = lam_here
                else
                    pt_cl(lc) = n_cl
                    pt_lam(lc) = lam_here
                end if

                n_cl = n_cl + 1
                ct(n_cl)%parent    = cl
                ct(n_cl)%lam_birth = lam_here
                ct(n_cl)%n_born    = sz_rc
                ct(cl)%child2      = n_cl
                if (rc > n) then
                    stk_top = stk_top + 1
                    stk_nd(stk_top) = rc
                    stk_cl(stk_top) = n_cl
                    stk_lb(stk_top) = lam_here
                else
                    pt_cl(rc) = n_cl
                    pt_lam(rc) = lam_here
                end if

            else if (sz_lc >= min_cl_sz) then
                ct(cl)%stability = ct(cl)%stability + &
                    real(sz_rc, rk) * (lam_here - lam_b)
                call mark_subtree(rc, cl, lam_here, n, dg, pt_cl, pt_lam)
                if (lc > n) then
                    stk_top = stk_top + 1
                    stk_nd(stk_top) = lc
                    stk_cl(stk_top) = cl
                    stk_lb(stk_top) = lam_b
                else
                    pt_cl(lc) = cl
                    pt_lam(lc) = lam_here
                end if

            else if (sz_rc >= min_cl_sz) then
                ct(cl)%stability = ct(cl)%stability + &
                    real(sz_lc, rk) * (lam_here - lam_b)
                call mark_subtree(lc, cl, lam_here, n, dg, pt_cl, pt_lam)
                if (rc > n) then
                    stk_top = stk_top + 1
                    stk_nd(stk_top) = rc
                    stk_cl(stk_top) = cl
                    stk_lb(stk_top) = lam_b
                else
                    pt_cl(rc) = cl
                    pt_lam(rc) = lam_here
                end if

            else
                ct(cl)%stability = ct(cl)%stability + &
                    real(sz_lc + sz_rc, rk) * (lam_here - lam_b)
                call mark_subtree(lc, cl, lam_here, n, dg, pt_cl, pt_lam)
                call mark_subtree(rc, cl, lam_here, n, dg, pt_cl, pt_lam)
            end if

        end do
    end subroutine extract_condensed_tree

    !---------------------------------------------------------------------------
    ! Excess-of-Mass cluster selection — correct greedy bottom-up algorithm.
    !
    ! Key insight: we maintain a separate is_sel(:) array.
    ! When children beat parent: parent is deselected, children remain selected.
    ! When parent beats children: parent is selected, children deselected recursively.
    !---------------------------------------------------------------------------
    subroutine eom_select(ct, n_cl, params, is_sel)
        type(ct_node),        intent(in)  :: ct(:)
        integer(ik),          intent(in)  :: n_cl
        type(hdbscan_params), intent(in)  :: params
        logical,              intent(out) :: is_sel(n_cl)

        real(rk),    allocatable :: score(:)
        integer(ik) :: i, c1, c2
        real(rk)    :: child_sum

        allocate(score(n_cl), source=0.0_rk)

        ! Start: leaf clusters (no children) are always selected
        do i = 1, n_cl
            is_sel(i) = (ct(i)%child1 == 0 .and. ct(i)%child2 == 0)
            score(i)  = ct(i)%stability
        end do

        ! Bottom-up (children before parents; children have higher ids)
        do i = n_cl, 1, -1
            c1 = ct(i)%child1
            c2 = ct(i)%child2
            if (c1 <= 0 .and. c2 <= 0) then
                cycle   ! leaf: already handled above
            end if

            child_sum = 0.0_rk
            if (c1 > 0 .and. c1 <= n_cl) then
                child_sum = child_sum + score(c1)
            end if
            if (c2 > 0 .and. c2 <= n_cl) then
                child_sum = child_sum + score(c2)
            end if

            if (i == 1 .and. .not. params%allow_single_cluster) then
                is_sel(i) = .false.
                score(i)  = child_sum
            else if (ct(i)%stability >= child_sum) then
                ! Parent wins: select parent, deselect all descendants
                is_sel(i) = .true.
                score(i)  = ct(i)%stability
                call deselect_subtree(is_sel, ct, n_cl, c1)
                call deselect_subtree(is_sel, ct, n_cl, c2)
            else
                ! Children win: deselect parent, propagate children scores upward
                is_sel(i) = .false.
                score(i)  = child_sum
            end if
        end do

        ! If nothing is selected (can happen if all clusters have zero stability),
        ! fall back to the root — equivalent to allow_single_cluster=true
        if (.not. any(is_sel(1:n_cl))) then
            if (params%allow_single_cluster .or. n_cl == 1) then
                is_sel(1) = .true.
            end if
        end if

    end subroutine eom_select

    ! Recursively deselect a cluster and all its descendants
    recursive subroutine deselect_subtree(is_sel, ct, n_cl, c)
        logical,       intent(inout) :: is_sel(:)
        type(ct_node), intent(in)    :: ct(:)
        integer(ik),   intent(in)    :: n_cl, c
        if (c <= 0 .or. c > n_cl) then
            return
        end if
        is_sel(c) = .false.
        call deselect_subtree(is_sel, ct, n_cl, ct(c)%child1)
        call deselect_subtree(is_sel, ct, n_cl, ct(c)%child2)
    end subroutine deselect_subtree

    subroutine assign_labels(pt_cl, pt_lam, ct, is_sel, n_pts, n_cl, res)
        integer(ik),    intent(in)    :: n_pts, n_cl
        integer(ik),    intent(in)    :: pt_cl(n_pts)
        real(rk),       intent(in)    :: pt_lam(n_pts)
        type(ct_node),  intent(in)    :: ct(2*n_pts)
        logical,        intent(in)    :: is_sel(n_cl)
        type(hdbscan_result), intent(inout) :: res

        integer(ik) :: i, c, next_lbl
        integer(ik), allocatable :: cl_lbl(:)
        real(rk),    allocatable :: cl_max_lam(:)

        allocate(cl_lbl(n_cl),     source=0_ik)
        allocate(cl_max_lam(n_cl), source=0.0_rk)

        ! Assign output label ids to each selected condensed-tree cluster
        next_lbl = 1
        do i = 1, n_cl
            if (is_sel(i)) then
                cl_lbl(i) = next_lbl
                next_lbl = next_lbl + 1
            end if
        end do

        ! For each point, walk up the condensed tree to the nearest selected cluster
        do i = 1, n_pts
            c = pt_cl(i)
            if (c == 0) then
                cycle
            end if
            do while (c > 0 .and. c <= n_cl)
                if (is_sel(c)) then
                    res%labels(i) = cl_lbl(c)
                    exit
                end if
                c = ct(c)%parent
            end do
            ! Accumulate max lambda per selected cluster (for probability normalisation)
            if (res%labels(i) > 0) then
                c = pt_cl(i)
                do while (c > 0 .and. c <= n_cl)
                    if (is_sel(c)) then
                        exit
                    end if
                    c = ct(c)%parent
                end do
                if (c > 0 .and. c <= n_cl) then
                    cl_max_lam(c) = max(cl_max_lam(c), pt_lam(i))
                end if
            end if
        end do

        ! Membership probabilities = pt_lam / max_lam_in_cluster
        do i = 1, n_pts
            if (res%labels(i) <= 0) then
                cycle
            end if
            c = pt_cl(i)
            do while (c > 0 .and. c <= n_cl)
                if (is_sel(c)) then
                    exit
                end if
                c = ct(c)%parent
            end do
            if (c > 0 .and. c <= n_cl .and. cl_max_lam(c) > 0.0_rk) then
                res%probabilities(i) = min(1.0_rk, pt_lam(i) / cl_max_lam(c))
            else
                res%probabilities(i) = 1.0_rk
            end if
        end do

    end subroutine assign_labels

    subroutine compute_glosh(res, n)
        type(hdbscan_result), intent(inout) :: res
        integer(ik),          intent(in)    :: n

        integer(ik) :: i, lbl, mx_cl
        real(rk),    allocatable :: cl_max(:)

        mx_cl = maxval(res%labels, mask=(res%labels > 0))
        if (mx_cl <= 0) then
            res%outlier_scores = 1.0_rk
            return
        end if

        allocate(cl_max(mx_cl), source=0.0_rk)
        do i = 1, n
            lbl = res%labels(i)
            if (lbl < 1 .or. lbl > mx_cl) then
                cycle
            end if
            cl_max(lbl) = max(cl_max(lbl), res%probabilities(i))
        end do

        do i = 1, n
            lbl = res%labels(i)
            if (lbl < 1 .or. lbl > mx_cl) then
                res%outlier_scores(i) = 1.0_rk
            else if (cl_max(lbl) > 0.0_rk) then
                res%outlier_scores(i) = 1.0_rk - res%probabilities(i) / cl_max(lbl)
            else
                res%outlier_scores(i) = 1.0_rk
            end if
        end do

    end subroutine compute_glosh

end module hdbscan_mod
