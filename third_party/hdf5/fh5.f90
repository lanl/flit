!==============================================================================
! fh5 v2 — Minimal Fortran HDF5 interface
! Types : i8/i16/i32/i64, r32/r64, c32/c64 (complex stored as [2,...] real)
! Ranks : 0–6 scalars/arrays
! API   : fh5_open/close/finalize
!         fh5_write / fh5_read        (overloaded for all type×rank)
!         fh5_write_str / fh5_read_str
!         fh5_write_attr / fh5_read_attr / fh5_print_attr
!         fh5_exists / fh5_ndim / fh5_shape / fh5_dtype
!         fh5_print_tree
!==============================================================================
module fh5

    use hdf5
    use iso_fortran_env, only: int8, int16, int32, int64, real32, real64
    implicit none
    private

    ! Public kind aliases
    integer, parameter, public :: fh5_i8 = int8
    integer, parameter, public :: fh5_i16 = int16
    integer, parameter, public :: fh5_i32 = int32
    integer, parameter, public :: fh5_i64 = int64
    integer, parameter, public :: fh5_r32 = real32
    integer, parameter, public :: fh5_r64 = real64
    integer, parameter, public :: fh5_c32 = real32   ! COMPLEX(fh5_c32)
    integer, parameter, public :: fh5_c64 = real64   ! COMPLEX(fh5_c64)

    ! Private kind aliases (used inside module)
    integer, parameter :: i8 = int8
    integer, parameter :: i16 = int16
    integer, parameter :: i32 = int32
    integer, parameter :: i64 = int64
    integer, parameter :: r32 = real32
    integer, parameter :: r64 = real64

    logical, save :: lib_open = .false.

    integer(hid_t) :: fh5_fid

    ! Public API
    public :: fh5_open, fh5_close, fh5_finalize
    public :: fh5_write, fh5_read
    public :: fh5_write_str, fh5_read_str
    public :: fh5_write_attr, fh5_read_attr, fh5_print_attr
    public :: fh5_exists, fh5_ndim, fh5_shape, fh5_dtype
    public :: fh5_print_tree
    public :: fh5_write_chunk, fh5_delete, fh5_list, fh5_copy
    public :: fh5_fid

    interface fh5_write
        module procedure :: wr_i8_0, wr_i8_1, wr_i8_2, wr_i8_3, wr_i8_4, wr_i8_5, wr_i8_6
        module procedure :: wr_i16_0, wr_i16_1, wr_i16_2, wr_i16_3, wr_i16_4, wr_i16_5, wr_i16_6
        module procedure :: wr_i32_0, wr_i32_1, wr_i32_2, wr_i32_3, wr_i32_4, wr_i32_5, wr_i32_6
        module procedure :: wr_i64_0, wr_i64_1, wr_i64_2, wr_i64_3, wr_i64_4, wr_i64_5, wr_i64_6
        module procedure :: wr_r32_0, wr_r32_1, wr_r32_2, wr_r32_3, wr_r32_4, wr_r32_5, wr_r32_6
        module procedure :: wr_r64_0, wr_r64_1, wr_r64_2, wr_r64_3, wr_r64_4, wr_r64_5, wr_r64_6
        module procedure :: wr_c32_0, wr_c32_1, wr_c32_2, wr_c32_3, wr_c32_4, wr_c32_5, wr_c32_6
        module procedure :: wr_c64_0, wr_c64_1, wr_c64_2, wr_c64_3, wr_c64_4, wr_c64_5, wr_c64_6
    end interface fh5_write

    interface fh5_read
        module procedure :: rd_i8_0, rd_i8_1, rd_i8_2, rd_i8_3, rd_i8_4, rd_i8_5, rd_i8_6
        module procedure :: rd_i16_0, rd_i16_1, rd_i16_2, rd_i16_3, rd_i16_4, rd_i16_5, rd_i16_6
        module procedure :: rd_i32_0, rd_i32_1, rd_i32_2, rd_i32_3, rd_i32_4, rd_i32_5, rd_i32_6
        module procedure :: rd_i64_0, rd_i64_1, rd_i64_2, rd_i64_3, rd_i64_4, rd_i64_5, rd_i64_6
        module procedure :: rd_r32_0, rd_r32_1, rd_r32_2, rd_r32_3, rd_r32_4, rd_r32_5, rd_r32_6
        module procedure :: rd_r64_0, rd_r64_1, rd_r64_2, rd_r64_3, rd_r64_4, rd_r64_5, rd_r64_6
        module procedure :: rd_c32_0, rd_c32_1, rd_c32_2, rd_c32_3, rd_c32_4, rd_c32_5, rd_c32_6
        module procedure :: rd_c64_0, rd_c64_1, rd_c64_2, rd_c64_3, rd_c64_4, rd_c64_5, rd_c64_6
    end interface fh5_read

    ! Attribute interfaces
    interface fh5_write_attr
        module procedure :: wa_i32, wa_i64, wa_r32, wa_r64, wa_str
        module procedure :: wa_i32_1, wa_i64_1, wa_r32_1, wa_r64_1
    end interface
    interface fh5_read_attr
        module procedure :: ra_i32, ra_i64, ra_r32, ra_r64, ra_str
        module procedure :: ra_i32_1, ra_i64_1, ra_r32_1, ra_r64_1
    end interface

    interface fh5_write_chunk
        module procedure :: fh5_write_chunk_r64_1, fh5_write_chunk_r32_3
    end interface fh5_write_chunk

contains

    !================================================================
    ! FILE OPERATIONS
    !================================================================
    subroutine fh5_open(filename, fid, mode, error)
        character(*), intent(in) :: filename
        integer(hid_t), intent(out) :: fid
        character(*), intent(in), optional :: mode
        integer, intent(out), optional :: error
        character :: m
        integer :: e
        call lib_init(e)
        if (bail(e, "fh5_open", error)) return
        m = 'w'
        if (present(mode)) m = mode(1:1)
        select case (m)
            case ('r', 'R')
                call h5fopen_f(filename, H5F_ACC_RDONLY_F, fid, e)
            case ('a', 'A')
                call h5fopen_f(filename, H5F_ACC_RDWR_F, fid, e)
            case default
                call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, e)
        end select
        call chk(e, "fh5_open: "//trim(filename), error)
    end subroutine

    subroutine fh5_close(fid, error)
        integer(hid_t), intent(in) :: fid
        integer, intent(out), optional :: error
        integer :: e
        call h5fclose_f(fid, e)
        call chk(e, "fh5_close", error)
    end subroutine

    subroutine fh5_finalize(error)
        integer, intent(out), optional :: error
        integer :: e
        if (lib_open) then
            call h5close_f(e)
            lib_open = .false.
            call chk(e, "fh5_finalize", error)
        else
            if (present(error)) error = 0
        end if
    end subroutine

    !================================================================
    ! STRING DATASETS
    !================================================================
    subroutine fh5_write_str(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path, d
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did, tid
        integer(hsize_t) :: dims(1)
        integer(size_t) :: slen
        integer :: e
        call lib_init(e)
        if (bail(e, "fh5_write_str", error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, "fh5_write_str grp", error)) return
        call h5screate_f(H5S_SCALAR_F, sid, e)
        call h5tcopy_f(H5T_NATIVE_CHARACTER, tid, e)
        slen = max(1_size_t, int(len_trim(d), size_t))
        call h5tset_size_f(tid, slen, e)
        call h5dcreate_f(fid, trim(path), tid, sid, did, e)
        if (bail(e, "fh5_write_str dset", error)) then
            call h5tclose_f(tid, e)
            call h5sclose_f(sid, e)
            return
        end if
        dims = [1_hsize_t]
        call h5dwrite_f(did, tid, trim(d), dims, e)
        call h5dclose_f(did, e)
        call h5tclose_f(tid, e)
        call h5sclose_f(sid, e)
        call chk(e, "fh5_write_str", error)
    end subroutine

    subroutine fh5_read_str(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        character(*), intent(out) :: d
        integer, intent(out), optional :: error
        integer(hid_t) :: did, tid
        integer(hsize_t) :: dims(1)
        integer :: e
        d = ''
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, "fh5_read_str", error)) return
        call h5dget_type_f(did, tid, e)
        dims = [1_hsize_t]
        call h5dread_f(did, tid, d, dims, e)
        call h5tclose_f(tid, e)
        call h5dclose_f(did, e)
        call chk(e, "fh5_read_str", error)
    end subroutine

    !================================================================
    ! QUERY
    !================================================================
    subroutine fh5_exists(fid, path, exists, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        logical, intent(out) :: exists
        integer, intent(out), optional :: error
        integer :: e
        call h5lexists_f(fid, trim(path), exists, e)
        call chk(e, "fh5_exists", error)
    end subroutine

    subroutine fh5_ndim(fid, path, n, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer, intent(out) :: n
        integer, intent(out), optional :: error
        integer(hid_t) :: did, sid
        integer :: e
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, "fh5_ndim", error)) return
        call h5dget_space_f(did, sid, e)
        call h5sget_simple_extent_ndims_f(sid, n, e)
        call h5sclose_f(sid, e)
        call h5dclose_f(did, e)
        call chk(e, "fh5_ndim", error)
    end subroutine

    subroutine fh5_shape(fid, path, dims, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(hsize_t), allocatable, intent(out) :: dims(:)
        integer, intent(out), optional :: error
        integer(hid_t) :: did, sid
        integer(hsize_t), allocatable :: md(:)
        integer :: e, nd
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, "fh5_shape", error)) return
        call h5dget_space_f(did, sid, e)
        call h5sget_simple_extent_ndims_f(sid, nd, e)
        allocate (dims(nd), md(nd))
        call h5sget_simple_extent_dims_f(sid, dims, md, e)
        call h5sclose_f(sid, e)
        call h5dclose_f(did, e)
        call chk(e, "fh5_shape", error)
    end subroutine

    subroutine fh5_dtype(fid, path, label, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        character(*), intent(out) :: label
        integer, intent(out), optional :: error
        integer(hid_t) :: did, tid, sid
        integer(hsize_t) :: dims(16), md(16)
        integer :: e, nd, tc, i
        integer(size_t) :: tsz
        character(len=32) :: sh
        label = '?'
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, "fh5_dtype", error)) return
        call h5dget_type_f(did, tid, e)
        call h5tget_class_f(tid, tc, e)
        call h5tget_size_f(tid, tsz, e)
        call h5dget_space_f(did, sid, e)
        call h5sget_simple_extent_ndims_f(sid, nd, e)
        if (nd > 0) call h5sget_simple_extent_dims_f(sid, dims(1:nd), md(1:nd), e)
        call h5sclose_f(sid, e)
        call h5tclose_f(tid, e)
        call h5dclose_f(did, e)
        ! Use if/elseif because H5T_INTEGER_F etc. are variables, not parameters
        if (tc == H5T_INTEGER_F) then
            if (int(tsz) == 1) then
                label = "INT8"
            else if (int(tsz) == 2) then
                label = "INT16"
            else if (int(tsz) == 4) then
                label = "INT32"
            else if (int(tsz) == 8) then
                label = "INT64"
            else
                write (label, '(A,I0,A)') "INT[", tsz, "B]"
            end if
        else if (tc == H5T_FLOAT_F) then
            if (int(tsz) == 4) then
                label = "FLOAT32"
            else
                label = "FLOAT64"
            end if
        else if (tc == H5T_STRING_F) then
            write (label, '(A,I0,A)') "STRING[", tsz, "B]"
        else
            write (label, '(A,I0)') "TYPE_CLASS=", tc
        end if
        ! Detect complex: dim-0 == 2 and __complex__ attr present
        block
            logical :: cx
            cx = .false.
            if (nd > 0 .and. dims(1) == 2) then
                block
                    character(len=4) :: cattr
                    integer :: ae
                    call fh5_read_attr(fid, trim(path), '__complex__', cattr, error=ae)
                    if (ae == 0 .and. trim(cattr) == '1') cx = .true.
                end block
            end if
            if (cx) then
                if (trim(label) == 'FLOAT32') then
                    label = 'COMPLEX32'
                else if (trim(label) == 'FLOAT64') then
                    label = 'COMPLEX64'
                end if
                nd = nd - 1   ! hide the leading 2 from shape display
                dims(1:nd) = dims(2:nd + 1)
            end if
        end block
        if (nd == 0) then
            label = trim(label)//" (scalar)"
        else
            sh = " ["
            do i = 1, nd
                if (i > 1) sh = trim(sh)//"x"
                write (sh(len_trim(sh) + 1:), '(I0)') dims(i)
            end do
            label = trim(label)//trim(sh)//"]"
        end if
        if (present(error)) error = 0
    end subroutine

    !================================================================
    ! TREE PRINTER — uses h5gn_members_f + h5gget_obj_info_idx_f
    !================================================================
    subroutine fh5_print_tree(fid, path, iu)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in), optional :: path
        integer, intent(in), optional :: iu
        integer :: unit_, e, na, ia
        integer(hid_t) :: oid, aid
        character(len=512) :: root, aname
        unit_ = 6
        if (present(iu)) then
            unit_ = iu
        end if
        root = '/'
        if (present(path)) then
            root = trim(path)
        end if
        write (unit_, '(A)') trim(root)
        ! Print attributes attached to the root group itself
        call h5gopen_f(fid, trim(root), oid, e)
        if (e == 0) then
            call h5aget_num_attrs_f(oid, na, e)
            do ia = 0, na - 1
                aname = ''
                call h5aget_name_by_idx_f(fid, trim(root), H5_INDEX_NAME_F, &
                    H5_ITER_INC_F, int(ia, hsize_t), aname, e)
                if (e /= 0) then
                    cycle
                end if
                call h5aopen_f(oid, trim(aname), aid, e)
                if (e /= 0) then
                    cycle
                end if
                write (unit_, '(A,A)') "  @"//trim(aname)//" = ", trim(attr_val_str(aid))
                call h5aclose_f(aid, e)
            end do
            call h5gclose_f(oid, e)
        end if
        call tree_node(fid, trim(root), "", unit_)
    end subroutine fh5_print_tree

    recursive subroutine tree_node(fid, grp_path, prefix, iu)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: grp_path, prefix
        integer, intent(in) :: iu
        integer :: e, nmembers, i, obj_type, na, ia
        integer(hid_t) :: oid, aid
        character(len=512) :: cname, cpath, tlabel, dlabel, aname
        logical :: is_last

        ! Count members
        call h5gn_members_f(fid, trim(grp_path), nmembers, e)
        if (e /= 0 .or. nmembers <= 0) then
            return
        end if

        do i = 0, nmembers - 1
            is_last = (i == nmembers - 1)
            cname = ''
            ! Get name and type of member i
            call h5gget_obj_info_idx_f(fid, trim(grp_path), i, cname, obj_type, e)
            if (e /= 0) then
                cycle
            end if

            if (trim(grp_path) == '/') then
                cpath = '/'//trim(cname)
            else
                cpath = trim(grp_path)//'/'//trim(cname)
            end if

            if (is_last) then
                tlabel = prefix//"└── "//trim(cname)
            else
                tlabel = prefix//"├── "//trim(cname)
            end if

            if (obj_type == H5G_GROUP_F) then
                write (iu, '(A,A)') trim(tlabel), "/"
                ! Print attributes attached to this group
                call h5gopen_f(fid, trim(cpath), oid, e)
                if (e == 0) then
                    call h5aget_num_attrs_f(oid, na, e)
                    do ia = 0, na - 1
                        aname = ''
                        call h5aget_name_by_idx_f(fid, trim(cpath), H5_INDEX_NAME_F, &
                            H5_ITER_INC_F, int(ia, hsize_t), aname, e)
                        if (e /= 0) then
                            cycle
                        end if
                        call h5aopen_f(oid, trim(aname), aid, e)
                        if (e /= 0) then
                            cycle
                        end if
                        if (is_last) then
                            write (iu, '(A,A)') prefix//"       @"//trim(aname)//" = ", &
                                trim(attr_val_str(aid))
                        else
                            write (iu, '(A,A)') prefix//"│      @"//trim(aname)//" = ", &
                                trim(attr_val_str(aid))
                        end if
                        call h5aclose_f(aid, e)
                    end do
                    call h5gclose_f(oid, e)
                end if
                ! Recurse into group children
                if (is_last) then
                    call tree_node(fid, trim(cpath), prefix//"    ", iu)
                else
                    call tree_node(fid, trim(cpath), prefix//"│   ", iu)
                end if

            else if (obj_type == H5G_DATASET_F) then
                call fh5_dtype(fid, trim(cpath), dlabel)
                write (iu, '(A,A,A,A)') trim(tlabel), "  <", trim(dlabel), ">"
                ! Print attributes attached to this dataset
                call h5dopen_f(fid, trim(cpath), oid, e)
                if (e == 0) then
                    call h5aget_num_attrs_f(oid, na, e)
                    do ia = 0, na - 1
                        aname = ''
                        call h5aget_name_by_idx_f(fid, trim(cpath), H5_INDEX_NAME_F, &
                            H5_ITER_INC_F, int(ia, hsize_t), aname, e)
                        if (e /= 0) then
                            cycle
                        end if
                        call h5aopen_f(oid, trim(aname), aid, e)
                        if (e /= 0) then
                            cycle
                        end if
                        if (trim(aname) /= '__complex__') then
                            if (is_last) then
                                write (iu, '(A,A)') prefix//"       @"//trim(aname)//" = ", &
                                    trim(attr_val_str(aid))
                            else
                                write (iu, '(A,A)') prefix//"│      @"//trim(aname)//" = ", &
                                    trim(attr_val_str(aid))
                            end if
                        end if
                        call h5aclose_f(aid, e)
                    end do
                    call h5dclose_f(oid, e)
                end if

            else
                write (iu, '(A)') trim(tlabel)//" <other>"
            end if
        end do
    end subroutine tree_node

    !> Format an open attribute handle as a short string
    function attr_val_str(aid) result(s)
        integer(hid_t), intent(in) :: aid
        character(len=128) :: s
        integer(hid_t) :: tid, sid
        integer :: e, tc, nd, i
        integer(size_t) :: tsz
        integer(hsize_t) :: dims(8), md(8), ntot
        integer(i32), allocatable :: ai32(:)
        integer(i64), allocatable :: ai64(:)
        real(r32), allocatable :: ar32(:)
        real(r64), allocatable :: ar64(:)
        character(len=64) :: vs
        s = '?'
        call h5aget_type_f(aid, tid, e)
        call h5aget_space_f(aid, sid, e)
        call h5tget_class_f(tid, tc, e)
        call h5tget_size_f(tid, tsz, e)
        call h5sget_simple_extent_ndims_f(sid, nd, e)
        dims(1:max(nd, 1)) = 1_hsize_t
        if (nd > 0) then
            call h5sget_simple_extent_dims_f(sid, dims(1:nd), md(1:nd), e)
        end if
        call h5sclose_f(sid, e)
        ntot = product(dims(1:max(nd, 1)))
        if (tc == H5T_INTEGER_F) then
            if (tsz <= 4) then
                allocate (ai32(ntot))
                call h5aread_f(aid, h5kind_to_type(i32, H5_INTEGER_KIND), ai32, dims(1:max(nd, 1)), e)
                if (ntot == 1) then
                    write (s, '(I0)') ai32(1)
                else
                    s = "["
                    do i = 1, ntot
                        if (i > 1) then
                            s = trim(s)//","
                        end if
                        write (s(len_trim(s) + 1:), '(I0)') ai32(i)
                    end do
                    s = trim(s)//"]"
                end if
            else
                allocate (ai64(ntot))
                call h5aread_f(aid, h5kind_to_type(i64, H5_INTEGER_KIND), ai64, dims(1:max(nd, 1)), e)
                if (ntot == 1) then
                    write (s, '(I0)') ai64(1)
                else
                    s = "["
                    do i = 1, ntot
                        if (i > 1) then
                            s = trim(s)//","
                        end if
                        write (s(len_trim(s) + 1:), '(I0)') ai64(i)
                    end do
                    s = trim(s)//"]"
                end if
            end if
        else if (tc == H5T_FLOAT_F) then
            if (tsz == 4) then
                allocate (ar32(ntot))
                call h5aread_f(aid, h5kind_to_type(r32, H5_REAL_KIND), ar32, dims(1:max(nd, 1)), e)
                if (ntot == 1) then
                    write (s, '(G14.6)') ar32(1)
                else
                    s = "["
                    do i = 1, ntot
                        if (i > 1) then
                            s = trim(s)//","
                        end if
                        write (s(len_trim(s) + 1:), '(G12.5)') ar32(i)
                    end do
                    s = trim(s)//"]"
                end if
            else
                allocate (ar64(ntot))
                call h5aread_f(aid, h5kind_to_type(r64, H5_REAL_KIND), ar64, dims(1:max(nd, 1)), e)
                if (ntot == 1) then
                    write (s, '(G20.12)') ar64(1)
                else
                    s = "["
                    do i = 1, ntot
                        if (i > 1) then
                            s = trim(s)//","
                        end if
                        write (s(len_trim(s) + 1:), '(G14.6)') ar64(i)
                    end do
                    s = trim(s)//"]"
                end if
            end if
        else if (tc == H5T_STRING_F) then
            vs = ''
            dims(1) = 1_hsize_t
            call h5aread_f(aid, tid, vs, dims(1:1), e)
            s = '"'//trim(vs)//'"'
        else
            write (s, '(A,I0)') "type_class=", tc
        end if
        call h5tclose_f(tid, e)
    end function

    !> Print one attribute's value, with auto type detection (mirrors fh5_write_attr)
    subroutine fh5_print_attr(fid, obj_path, attr_name, iu)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj_path, attr_name
        integer, intent(in), optional :: iu
        integer(hid_t) :: oid, aid
        integer :: e, unit_
        unit_ = 6
        if (present(iu)) then
            unit_ = iu
        end if
        call h5oopen_f(fid, trim(obj_path), oid, e)
        if (e /= 0) then
            write (unit_, '(A)') "fh5_print_attr: cannot open "//trim(obj_path)
            return
        end if
        call h5aopen_f(oid, trim(attr_name), aid, e)
        if (e /= 0) then
            write (unit_, '(A)') "fh5_print_attr: '"//trim(attr_name)//"' not found in "//trim(obj_path)
            call h5oclose_f(oid, e)
            return
        end if
        write (unit_, '(A,A)') trim(obj_path)//"@"//trim(attr_name)//" = ", trim(attr_val_str(aid))
        call h5aclose_f(aid, e)
        call h5oclose_f(oid, e)
    end subroutine

    !================================================================
    ! INTERNALS
    !================================================================
    subroutine lib_init(e)
        integer, intent(out) :: e
        e = 0
        if (.not. lib_open) then
            call h5open_f(e)
            if (e == 0) lib_open = .true.
            ! Silence HDF5's own error stack so internal probes don't print noise
            call h5eset_auto_f(0, e)
            e = 0
        end if
    end subroutine

    subroutine chk(e, msg, error)
        integer, intent(in) :: e
        character(*), intent(in) :: msg
        integer, intent(out), optional :: error
        if (present(error)) then
            error = e
        else if (e /= 0) then
            write (*, '(A,I0)') "fh5 ERROR ["//trim(msg)//"] code=", e
            error stop "fh5: aborting."
        end if
    end subroutine

    logical function bail(e, msg, error)
        integer, intent(in) :: e
        character(*), intent(in) :: msg
        integer, intent(out), optional :: error
        bail = (e /= 0)
        if (bail) call chk(e, msg, error)
    end function

    subroutine ensure_groups(fid, dpath, e)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: dpath
        integer, intent(out) :: e
        integer(hid_t) :: gid
        logical :: exists
        integer :: i, last, n
        character(len=1024) :: pfx
        e = 0
        n = len_trim(dpath)
        last = 1
        do i = n, 1, -1
            if (dpath(i:i) == '/') then
                last = i
                exit
            end if
        end do
        if (last <= 1) return
        do i = 2, last
            if (dpath(i:i) == '/') then
                pfx = dpath(1:i - 1)
                call h5lexists_f(fid, trim(pfx), exists, e)
                if (e /= 0) return
                if (.not. exists) then
                    call h5gcreate_f(fid, trim(pfx), gid, e)
                    if (e /= 0) return
                    call h5gclose_f(gid, e)
                    if (e /= 0) return
                end if
            end if
        end do
        pfx = dpath(1:last - 1)
        call h5lexists_f(fid, trim(pfx), exists, e)
        if (e /= 0) return
        if (.not. exists) then
            call h5gcreate_f(fid, trim(pfx), gid, e)
            if (e /= 0) return
            call h5gclose_f(gid, e)
        end if
    end subroutine

    subroutine ensure_obj_exists(fid, path, e)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer, intent(out) :: e
        logical :: exists
        integer(hid_t) :: gid
        e = 0
        call h5lexists_f(fid, trim(path), exists, e)
        if (e /= 0) return
        if (.not. exists) then
            call h5gcreate_f(fid, trim(path), gid, e)
            if (e /= 0) return
            call h5gclose_f(gid, e)
        end if
    end subroutine

    subroutine get_dims(fid, path, nd, dims, e)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer, intent(in) :: nd
        integer(hsize_t), intent(out) :: dims(nd)
        integer, intent(out) :: e
        integer(hid_t) :: did, sid
        integer(hsize_t) :: md(nd)
        call h5dopen_f(fid, trim(path), did, e)
        if (e /= 0) return
        call h5dget_space_f(did, sid, e)
        call h5sget_simple_extent_dims_f(sid, dims, md, e)
        call h5sclose_f(sid, e)
        call h5dclose_f(did, e)
    end subroutine

    subroutine open_obj(fid, obj_path, oid, e)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj_path
        integer(hid_t), intent(out) :: oid
        integer, intent(out) :: e
        call h5oopen_f(fid, trim(obj_path), oid, e)
    end subroutine

    subroutine wr_i8_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i8), intent(in) :: d
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(1)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i8_0', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i8_0g', error)) return
        call h5screate_f(H5S_SCALAR_F, sid, e)
        if (bail(e, 'wr_i8_0s', error)) return
        dims = [1_hsize_t]
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i8, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i8_0d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i8, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i8_0', error)
    end subroutine wr_i8_0

    subroutine wr_i8_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i8), intent(in) :: d(:)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(1)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i8_1', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i8_1g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(1, dims, sid, e)
        if (bail(e, 'wr_i8_1s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i8, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i8_1d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i8, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i8_1', error)
    end subroutine wr_i8_1

    subroutine wr_i8_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i8), intent(in) :: d(:, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(2)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i8_2', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i8_2g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(2, dims, sid, e)
        if (bail(e, 'wr_i8_2s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i8, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i8_2d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i8, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i8_2', error)
    end subroutine wr_i8_2

    subroutine wr_i8_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i8), intent(in) :: d(:, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(3)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i8_3', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i8_3g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(3, dims, sid, e)
        if (bail(e, 'wr_i8_3s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i8, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i8_3d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i8, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i8_3', error)
    end subroutine wr_i8_3

    subroutine wr_i8_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i8), intent(in) :: d(:, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(4)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i8_4', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i8_4g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(4, dims, sid, e)
        if (bail(e, 'wr_i8_4s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i8, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i8_4d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i8, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i8_4', error)
    end subroutine wr_i8_4

    subroutine wr_i8_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i8), intent(in) :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(5)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i8_5', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i8_5g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(5, dims, sid, e)
        if (bail(e, 'wr_i8_5s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i8, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i8_5d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i8, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i8_5', error)
    end subroutine wr_i8_5

    subroutine wr_i8_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i8), intent(in) :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(6)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i8_6', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i8_6g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(6, dims, sid, e)
        if (bail(e, 'wr_i8_6s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i8, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i8_6d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i8, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i8_6', error)
    end subroutine wr_i8_6

    subroutine wr_i16_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i16), intent(in) :: d
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(1)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i16_0', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i16_0g', error)) return
        call h5screate_f(H5S_SCALAR_F, sid, e)
        if (bail(e, 'wr_i16_0s', error)) return
        dims = [1_hsize_t]
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i16, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i16_0d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i16, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i16_0', error)
    end subroutine wr_i16_0

    subroutine wr_i16_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i16), intent(in) :: d(:)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(1)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i16_1', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i16_1g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(1, dims, sid, e)
        if (bail(e, 'wr_i16_1s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i16, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i16_1d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i16, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i16_1', error)
    end subroutine wr_i16_1

    subroutine wr_i16_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i16), intent(in) :: d(:, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(2)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i16_2', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i16_2g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(2, dims, sid, e)
        if (bail(e, 'wr_i16_2s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i16, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i16_2d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i16, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i16_2', error)
    end subroutine wr_i16_2

    subroutine wr_i16_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i16), intent(in) :: d(:, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(3)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i16_3', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i16_3g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(3, dims, sid, e)
        if (bail(e, 'wr_i16_3s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i16, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i16_3d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i16, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i16_3', error)
    end subroutine wr_i16_3

    subroutine wr_i16_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i16), intent(in) :: d(:, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(4)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i16_4', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i16_4g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(4, dims, sid, e)
        if (bail(e, 'wr_i16_4s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i16, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i16_4d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i16, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i16_4', error)
    end subroutine wr_i16_4

    subroutine wr_i16_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i16), intent(in) :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(5)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i16_5', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i16_5g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(5, dims, sid, e)
        if (bail(e, 'wr_i16_5s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i16, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i16_5d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i16, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i16_5', error)
    end subroutine wr_i16_5

    subroutine wr_i16_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i16), intent(in) :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(6)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i16_6', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i16_6g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(6, dims, sid, e)
        if (bail(e, 'wr_i16_6s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i16, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i16_6d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i16, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i16_6', error)
    end subroutine wr_i16_6

    subroutine wr_i32_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i32), intent(in) :: d
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(1)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i32_0', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i32_0g', error)) return
        call h5screate_f(H5S_SCALAR_F, sid, e)
        if (bail(e, 'wr_i32_0s', error)) return
        dims = [1_hsize_t]
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i32, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i32_0d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i32, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i32_0', error)
    end subroutine wr_i32_0

    subroutine wr_i32_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i32), intent(in) :: d(:)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(1)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i32_1', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i32_1g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(1, dims, sid, e)
        if (bail(e, 'wr_i32_1s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i32, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i32_1d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i32, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i32_1', error)
    end subroutine wr_i32_1

    subroutine wr_i32_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i32), intent(in) :: d(:, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(2)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i32_2', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i32_2g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(2, dims, sid, e)
        if (bail(e, 'wr_i32_2s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i32, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i32_2d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i32, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i32_2', error)
    end subroutine wr_i32_2

    subroutine wr_i32_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i32), intent(in) :: d(:, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(3)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i32_3', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i32_3g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(3, dims, sid, e)
        if (bail(e, 'wr_i32_3s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i32, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i32_3d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i32, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i32_3', error)
    end subroutine wr_i32_3

    subroutine wr_i32_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i32), intent(in) :: d(:, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(4)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i32_4', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i32_4g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(4, dims, sid, e)
        if (bail(e, 'wr_i32_4s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i32, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i32_4d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i32, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i32_4', error)
    end subroutine wr_i32_4

    subroutine wr_i32_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i32), intent(in) :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(5)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i32_5', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i32_5g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(5, dims, sid, e)
        if (bail(e, 'wr_i32_5s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i32, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i32_5d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i32, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i32_5', error)
    end subroutine wr_i32_5

    subroutine wr_i32_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i32), intent(in) :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(6)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i32_6', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i32_6g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(6, dims, sid, e)
        if (bail(e, 'wr_i32_6s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i32, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i32_6d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i32, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i32_6', error)
    end subroutine wr_i32_6

    subroutine wr_i64_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i64), intent(in) :: d
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(1)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i64_0', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i64_0g', error)) return
        call h5screate_f(H5S_SCALAR_F, sid, e)
        if (bail(e, 'wr_i64_0s', error)) return
        dims = [1_hsize_t]
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i64, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i64_0d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i64, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i64_0', error)
    end subroutine wr_i64_0

    subroutine wr_i64_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i64), intent(in) :: d(:)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(1)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i64_1', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i64_1g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(1, dims, sid, e)
        if (bail(e, 'wr_i64_1s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i64, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i64_1d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i64, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i64_1', error)
    end subroutine wr_i64_1

    subroutine wr_i64_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i64), intent(in) :: d(:, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(2)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i64_2', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i64_2g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(2, dims, sid, e)
        if (bail(e, 'wr_i64_2s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i64, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i64_2d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i64, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i64_2', error)
    end subroutine wr_i64_2

    subroutine wr_i64_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i64), intent(in) :: d(:, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(3)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i64_3', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i64_3g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(3, dims, sid, e)
        if (bail(e, 'wr_i64_3s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i64, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i64_3d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i64, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i64_3', error)
    end subroutine wr_i64_3

    subroutine wr_i64_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i64), intent(in) :: d(:, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(4)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i64_4', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i64_4g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(4, dims, sid, e)
        if (bail(e, 'wr_i64_4s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i64, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i64_4d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i64, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i64_4', error)
    end subroutine wr_i64_4

    subroutine wr_i64_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i64), intent(in) :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(5)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i64_5', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i64_5g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(5, dims, sid, e)
        if (bail(e, 'wr_i64_5s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i64, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i64_5d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i64, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i64_5', error)
    end subroutine wr_i64_5

    subroutine wr_i64_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i64), intent(in) :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(6)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_i64_6', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_i64_6g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(6, dims, sid, e)
        if (bail(e, 'wr_i64_6s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(i64, H5_INTEGER_KIND), sid, did, e)
        if (bail(e, 'wr_i64_6d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(i64, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_i64_6', error)
    end subroutine wr_i64_6

    subroutine wr_r32_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(in) :: d
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(1)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_r32_0', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_r32_0g', error)) return
        call h5screate_f(H5S_SCALAR_F, sid, e)
        if (bail(e, 'wr_r32_0s', error)) return
        dims = [1_hsize_t]
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_r32_0d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_r32_0', error)
    end subroutine wr_r32_0

    subroutine wr_r32_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(in) :: d(:)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(1)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_r32_1', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_r32_1g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(1, dims, sid, e)
        if (bail(e, 'wr_r32_1s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_r32_1d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_r32_1', error)
    end subroutine wr_r32_1

    subroutine wr_r32_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(in) :: d(:, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(2)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_r32_2', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_r32_2g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(2, dims, sid, e)
        if (bail(e, 'wr_r32_2s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_r32_2d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_r32_2', error)
    end subroutine wr_r32_2

    subroutine wr_r32_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(in) :: d(:, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(3)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_r32_3', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_r32_3g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(3, dims, sid, e)
        if (bail(e, 'wr_r32_3s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_r32_3d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_r32_3', error)
    end subroutine wr_r32_3

    subroutine wr_r32_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(in) :: d(:, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(4)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_r32_4', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_r32_4g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(4, dims, sid, e)
        if (bail(e, 'wr_r32_4s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_r32_4d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_r32_4', error)
    end subroutine wr_r32_4

    subroutine wr_r32_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(in) :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(5)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_r32_5', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_r32_5g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(5, dims, sid, e)
        if (bail(e, 'wr_r32_5s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_r32_5d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_r32_5', error)
    end subroutine wr_r32_5

    subroutine wr_r32_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(in) :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(6)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_r32_6', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_r32_6g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(6, dims, sid, e)
        if (bail(e, 'wr_r32_6s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_r32_6d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_r32_6', error)
    end subroutine wr_r32_6

    subroutine wr_r64_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(in) :: d
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(1)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_r64_0', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_r64_0g', error)) return
        call h5screate_f(H5S_SCALAR_F, sid, e)
        if (bail(e, 'wr_r64_0s', error)) return
        dims = [1_hsize_t]
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_r64_0d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_r64_0', error)
    end subroutine wr_r64_0

    subroutine wr_r64_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(in) :: d(:)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(1)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_r64_1', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_r64_1g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(1, dims, sid, e)
        if (bail(e, 'wr_r64_1s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_r64_1d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_r64_1', error)
    end subroutine wr_r64_1

    subroutine wr_r64_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(in) :: d(:, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(2)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_r64_2', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_r64_2g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(2, dims, sid, e)
        if (bail(e, 'wr_r64_2s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_r64_2d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_r64_2', error)
    end subroutine wr_r64_2

    subroutine wr_r64_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(in) :: d(:, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(3)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_r64_3', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_r64_3g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(3, dims, sid, e)
        if (bail(e, 'wr_r64_3s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_r64_3d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_r64_3', error)
    end subroutine wr_r64_3

    subroutine wr_r64_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(in) :: d(:, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(4)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_r64_4', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_r64_4g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(4, dims, sid, e)
        if (bail(e, 'wr_r64_4s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_r64_4d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_r64_4', error)
    end subroutine wr_r64_4

    subroutine wr_r64_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(in) :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(5)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_r64_5', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_r64_5g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(5, dims, sid, e)
        if (bail(e, 'wr_r64_5s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_r64_5d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_r64_5', error)
    end subroutine wr_r64_5

    subroutine wr_r64_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(in) :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(6)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_r64_6', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_r64_6g', error)) return
        dims = shape(d, kind=hsize_t)
        call h5screate_simple_f(6, dims, sid, e)
        if (bail(e, 'wr_r64_6s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_r64_6d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        call chk(e, 'wr_r64_6', error)
    end subroutine wr_r64_6

    subroutine wr_c32_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r32), intent(in) :: d
        integer, intent(out), optional :: error
        real(r32) :: buf(2)
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(1)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_c32_0', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_c32_0g', error)) return
        buf(1) = real(d, r32)
        buf(2) = aimag(d)
        dims = [2_hsize_t]
        call h5screate_simple_f(1, dims, sid, e)
        if (bail(e, 'wr_c32_0s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_c32_0d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        if (e == 0) call fh5_write_attr(fid, trim(path), '__complex__', '1', error=e)
        call chk(e, 'wr_c32_0', error)
    end subroutine wr_c32_0

    subroutine wr_c32_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r32), intent(in) :: d(:)
        integer, intent(out), optional :: error
        real(r32) :: buf(2, size(d, 1))
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(2)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_c32_1', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_c32_1g', error)) return
        buf = reshape(transfer(d, buf), shape(buf))
        dims = [2_hsize_t, &
            int(size(d, 1), hsize_t)]
        call h5screate_simple_f(2, dims, sid, e)
        if (bail(e, 'wr_c32_1s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_c32_1d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        if (e == 0) call fh5_write_attr(fid, trim(path), '__complex__', '1', error=e)
        call chk(e, 'wr_c32_1', error)
    end subroutine wr_c32_1

    subroutine wr_c32_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r32), intent(in) :: d(:, :)
        integer, intent(out), optional :: error
        real(r32) :: buf(2, size(d, 1), size(d, 2))
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(3)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_c32_2', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_c32_2g', error)) return
        buf = reshape(transfer(d, buf), shape(buf))
        dims = [2_hsize_t, &
            int(size(d, 1), hsize_t), &
            int(size(d, 2), hsize_t)]
        call h5screate_simple_f(3, dims, sid, e)
        if (bail(e, 'wr_c32_2s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_c32_2d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        if (e == 0) call fh5_write_attr(fid, trim(path), '__complex__', '1', error=e)
        call chk(e, 'wr_c32_2', error)
    end subroutine wr_c32_2

    subroutine wr_c32_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r32), intent(in) :: d(:, :, :)
        integer, intent(out), optional :: error
        real(r32) :: buf(2, size(d, 1), size(d, 2), size(d, 3))
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(4)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_c32_3', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_c32_3g', error)) return
        buf = reshape(transfer(d, buf), shape(buf))
        dims = [2_hsize_t, &
            int(size(d, 1), hsize_t), &
            int(size(d, 2), hsize_t), &
            int(size(d, 3), hsize_t)]
        call h5screate_simple_f(4, dims, sid, e)
        if (bail(e, 'wr_c32_3s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_c32_3d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        if (e == 0) call fh5_write_attr(fid, trim(path), '__complex__', '1', error=e)
        call chk(e, 'wr_c32_3', error)
    end subroutine wr_c32_3

    subroutine wr_c32_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r32), intent(in) :: d(:, :, :, :)
        integer, intent(out), optional :: error
        real(r32) :: buf(2, size(d, 1), size(d, 2), size(d, 3), size(d, 4))
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(5)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_c32_4', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_c32_4g', error)) return
        buf = reshape(transfer(d, buf), shape(buf))
        dims = [2_hsize_t, &
            int(size(d, 1), hsize_t), &
            int(size(d, 2), hsize_t), &
            int(size(d, 3), hsize_t), &
            int(size(d, 4), hsize_t)]
        call h5screate_simple_f(5, dims, sid, e)
        if (bail(e, 'wr_c32_4s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_c32_4d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        if (e == 0) call fh5_write_attr(fid, trim(path), '__complex__', '1', error=e)
        call chk(e, 'wr_c32_4', error)
    end subroutine wr_c32_4

    subroutine wr_c32_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r32), intent(in) :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        real(r32) :: buf(2, size(d, 1), size(d, 2), size(d, 3), size(d, 4), size(d, 5))
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(6)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_c32_5', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_c32_5g', error)) return
        buf = reshape(transfer(d, buf), shape(buf))
        dims = [2_hsize_t, &
            int(size(d, 1), hsize_t), &
            int(size(d, 2), hsize_t), &
            int(size(d, 3), hsize_t), &
            int(size(d, 4), hsize_t), &
            int(size(d, 5), hsize_t)]
        call h5screate_simple_f(6, dims, sid, e)
        if (bail(e, 'wr_c32_5s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_c32_5d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        if (e == 0) call fh5_write_attr(fid, trim(path), '__complex__', '1', error=e)
        call chk(e, 'wr_c32_5', error)
    end subroutine wr_c32_5

    subroutine wr_c32_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r32), intent(in) :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        real(r32) :: buf(2, size(d, 1), size(d, 2), size(d, 3), size(d, 4), size(d, 5), size(d, 6))
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(7)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_c32_6', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_c32_6g', error)) return
        buf = reshape(transfer(d, buf), shape(buf))
        dims = [2_hsize_t, &
            int(size(d, 1), hsize_t), &
            int(size(d, 2), hsize_t), &
            int(size(d, 3), hsize_t), &
            int(size(d, 4), hsize_t), &
            int(size(d, 5), hsize_t), &
            int(size(d, 6), hsize_t)]
        call h5screate_simple_f(7, dims, sid, e)
        if (bail(e, 'wr_c32_6s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_c32_6d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        if (e == 0) call fh5_write_attr(fid, trim(path), '__complex__', '1', error=e)
        call chk(e, 'wr_c32_6', error)
    end subroutine wr_c32_6

    subroutine wr_c64_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r64), intent(in) :: d
        integer, intent(out), optional :: error
        real(r64) :: buf(2)
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(1)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_c64_0', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_c64_0g', error)) return
        buf(1) = real(d, r64)
        buf(2) = aimag(d)
        dims = [2_hsize_t]
        call h5screate_simple_f(1, dims, sid, e)
        if (bail(e, 'wr_c64_0s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_c64_0d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        if (e == 0) call fh5_write_attr(fid, trim(path), '__complex__', '1', error=e)
        call chk(e, 'wr_c64_0', error)
    end subroutine wr_c64_0

    subroutine wr_c64_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r64), intent(in) :: d(:)
        integer, intent(out), optional :: error
        real(r64) :: buf(2, size(d, 1))
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(2)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_c64_1', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_c64_1g', error)) return
        buf = reshape(transfer(d, buf), shape(buf))
        dims = [2_hsize_t, &
            int(size(d, 1), hsize_t)]
        call h5screate_simple_f(2, dims, sid, e)
        if (bail(e, 'wr_c64_1s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_c64_1d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        if (e == 0) call fh5_write_attr(fid, trim(path), '__complex__', '1', error=e)
        call chk(e, 'wr_c64_1', error)
    end subroutine wr_c64_1

    subroutine wr_c64_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r64), intent(in) :: d(:, :)
        integer, intent(out), optional :: error
        real(r64) :: buf(2, size(d, 1), size(d, 2))
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(3)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_c64_2', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_c64_2g', error)) return
        buf = reshape(transfer(d, buf), shape(buf))
        dims = [2_hsize_t, &
            int(size(d, 1), hsize_t), &
            int(size(d, 2), hsize_t)]
        call h5screate_simple_f(3, dims, sid, e)
        if (bail(e, 'wr_c64_2s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_c64_2d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        if (e == 0) call fh5_write_attr(fid, trim(path), '__complex__', '1', error=e)
        call chk(e, 'wr_c64_2', error)
    end subroutine wr_c64_2

    subroutine wr_c64_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r64), intent(in) :: d(:, :, :)
        integer, intent(out), optional :: error
        real(r64) :: buf(2, size(d, 1), size(d, 2), size(d, 3))
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(4)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_c64_3', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_c64_3g', error)) return
        buf = reshape(transfer(d, buf), shape(buf))
        dims = [2_hsize_t, &
            int(size(d, 1), hsize_t), &
            int(size(d, 2), hsize_t), &
            int(size(d, 3), hsize_t)]
        call h5screate_simple_f(4, dims, sid, e)
        if (bail(e, 'wr_c64_3s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_c64_3d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        if (e == 0) call fh5_write_attr(fid, trim(path), '__complex__', '1', error=e)
        call chk(e, 'wr_c64_3', error)
    end subroutine wr_c64_3

    subroutine wr_c64_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r64), intent(in) :: d(:, :, :, :)
        integer, intent(out), optional :: error
        real(r64) :: buf(2, size(d, 1), size(d, 2), size(d, 3), size(d, 4))
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(5)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_c64_4', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_c64_4g', error)) return
        buf = reshape(transfer(d, buf), shape(buf))
        dims = [2_hsize_t, &
            int(size(d, 1), hsize_t), &
            int(size(d, 2), hsize_t), &
            int(size(d, 3), hsize_t), &
            int(size(d, 4), hsize_t)]
        call h5screate_simple_f(5, dims, sid, e)
        if (bail(e, 'wr_c64_4s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_c64_4d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        if (e == 0) call fh5_write_attr(fid, trim(path), '__complex__', '1', error=e)
        call chk(e, 'wr_c64_4', error)
    end subroutine wr_c64_4

    subroutine wr_c64_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r64), intent(in) :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        real(r64) :: buf(2, size(d, 1), size(d, 2), size(d, 3), size(d, 4), size(d, 5))
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(6)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_c64_5', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_c64_5g', error)) return
        buf = reshape(transfer(d, buf), shape(buf))
        dims = [2_hsize_t, &
            int(size(d, 1), hsize_t), &
            int(size(d, 2), hsize_t), &
            int(size(d, 3), hsize_t), &
            int(size(d, 4), hsize_t), &
            int(size(d, 5), hsize_t)]
        call h5screate_simple_f(6, dims, sid, e)
        if (bail(e, 'wr_c64_5s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_c64_5d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        if (e == 0) call fh5_write_attr(fid, trim(path), '__complex__', '1', error=e)
        call chk(e, 'wr_c64_5', error)
    end subroutine wr_c64_5

    subroutine wr_c64_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r64), intent(in) :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        real(r64) :: buf(2, size(d, 1), size(d, 2), size(d, 3), size(d, 4), size(d, 5), size(d, 6))
        integer(hid_t) :: sid, did
        integer(hsize_t) :: dims(7)
        integer :: e
        call lib_init(e)
        if (bail(e, 'wr_c64_6', error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, 'wr_c64_6g', error)) return
        buf = reshape(transfer(d, buf), shape(buf))
        dims = [2_hsize_t, &
            int(size(d, 1), hsize_t), &
            int(size(d, 2), hsize_t), &
            int(size(d, 3), hsize_t), &
            int(size(d, 4), hsize_t), &
            int(size(d, 5), hsize_t), &
            int(size(d, 6), hsize_t)]
        call h5screate_simple_f(7, dims, sid, e)
        if (bail(e, 'wr_c64_6s', error)) return
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e)
        if (bail(e, 'wr_c64_6d', error)) then
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        call h5sclose_f(sid, e)
        if (e == 0) call fh5_write_attr(fid, trim(path), '__complex__', '1', error=e)
        call chk(e, 'wr_c64_6', error)
    end subroutine wr_c64_6

    subroutine rd_i8_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i8), intent(out) :: d
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(1)
        integer :: e
        dims = [1_hsize_t]
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i8_0o', error)) return
        call h5dread_f(did, h5kind_to_type(i8, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i8_0', error)
    end subroutine rd_i8_0

    subroutine rd_i8_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i8), intent(out), allocatable :: d(:)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(1)
        integer :: e
        call get_dims(fid, path, 1, dims, e)
        if (bail(e, 'rd_i8_1d', error)) return
        allocate (d(dims(1)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i8_1o', error)) return
        call h5dread_f(did, h5kind_to_type(i8, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i8_1', error)
    end subroutine rd_i8_1

    subroutine rd_i8_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i8), intent(out), allocatable :: d(:, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(2)
        integer :: e
        call get_dims(fid, path, 2, dims, e)
        if (bail(e, 'rd_i8_2d', error)) return
        allocate (d(dims(1), dims(2)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i8_2o', error)) return
        call h5dread_f(did, h5kind_to_type(i8, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i8_2', error)
    end subroutine rd_i8_2

    subroutine rd_i8_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i8), intent(out), allocatable :: d(:, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(3)
        integer :: e
        call get_dims(fid, path, 3, dims, e)
        if (bail(e, 'rd_i8_3d', error)) return
        allocate (d(dims(1), dims(2), dims(3)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i8_3o', error)) return
        call h5dread_f(did, h5kind_to_type(i8, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i8_3', error)
    end subroutine rd_i8_3

    subroutine rd_i8_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i8), intent(out), allocatable :: d(:, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(4)
        integer :: e
        call get_dims(fid, path, 4, dims, e)
        if (bail(e, 'rd_i8_4d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i8_4o', error)) return
        call h5dread_f(did, h5kind_to_type(i8, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i8_4', error)
    end subroutine rd_i8_4

    subroutine rd_i8_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i8), intent(out), allocatable :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(5)
        integer :: e
        call get_dims(fid, path, 5, dims, e)
        if (bail(e, 'rd_i8_5d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4), dims(5)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i8_5o', error)) return
        call h5dread_f(did, h5kind_to_type(i8, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i8_5', error)
    end subroutine rd_i8_5

    subroutine rd_i8_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i8), intent(out), allocatable :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(6)
        integer :: e
        call get_dims(fid, path, 6, dims, e)
        if (bail(e, 'rd_i8_6d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i8_6o', error)) return
        call h5dread_f(did, h5kind_to_type(i8, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i8_6', error)
    end subroutine rd_i8_6

    subroutine rd_i16_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i16), intent(out) :: d
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(1)
        integer :: e
        dims = [1_hsize_t]
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i16_0o', error)) return
        call h5dread_f(did, h5kind_to_type(i16, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i16_0', error)
    end subroutine rd_i16_0

    subroutine rd_i16_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i16), intent(out), allocatable :: d(:)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(1)
        integer :: e
        call get_dims(fid, path, 1, dims, e)
        if (bail(e, 'rd_i16_1d', error)) return
        allocate (d(dims(1)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i16_1o', error)) return
        call h5dread_f(did, h5kind_to_type(i16, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i16_1', error)
    end subroutine rd_i16_1

    subroutine rd_i16_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i16), intent(out), allocatable :: d(:, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(2)
        integer :: e
        call get_dims(fid, path, 2, dims, e)
        if (bail(e, 'rd_i16_2d', error)) return
        allocate (d(dims(1), dims(2)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i16_2o', error)) return
        call h5dread_f(did, h5kind_to_type(i16, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i16_2', error)
    end subroutine rd_i16_2

    subroutine rd_i16_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i16), intent(out), allocatable :: d(:, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(3)
        integer :: e
        call get_dims(fid, path, 3, dims, e)
        if (bail(e, 'rd_i16_3d', error)) return
        allocate (d(dims(1), dims(2), dims(3)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i16_3o', error)) return
        call h5dread_f(did, h5kind_to_type(i16, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i16_3', error)
    end subroutine rd_i16_3

    subroutine rd_i16_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i16), intent(out), allocatable :: d(:, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(4)
        integer :: e
        call get_dims(fid, path, 4, dims, e)
        if (bail(e, 'rd_i16_4d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i16_4o', error)) return
        call h5dread_f(did, h5kind_to_type(i16, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i16_4', error)
    end subroutine rd_i16_4

    subroutine rd_i16_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i16), intent(out), allocatable :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(5)
        integer :: e
        call get_dims(fid, path, 5, dims, e)
        if (bail(e, 'rd_i16_5d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4), dims(5)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i16_5o', error)) return
        call h5dread_f(did, h5kind_to_type(i16, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i16_5', error)
    end subroutine rd_i16_5

    subroutine rd_i16_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i16), intent(out), allocatable :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(6)
        integer :: e
        call get_dims(fid, path, 6, dims, e)
        if (bail(e, 'rd_i16_6d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i16_6o', error)) return
        call h5dread_f(did, h5kind_to_type(i16, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i16_6', error)
    end subroutine rd_i16_6

    subroutine rd_i32_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i32), intent(out) :: d
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(1)
        integer :: e
        dims = [1_hsize_t]
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i32_0o', error)) return
        call h5dread_f(did, h5kind_to_type(i32, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i32_0', error)
    end subroutine rd_i32_0

    subroutine rd_i32_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i32), intent(out), allocatable :: d(:)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(1)
        integer :: e
        call get_dims(fid, path, 1, dims, e)
        if (bail(e, 'rd_i32_1d', error)) return
        allocate (d(dims(1)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i32_1o', error)) return
        call h5dread_f(did, h5kind_to_type(i32, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i32_1', error)
    end subroutine rd_i32_1

    subroutine rd_i32_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i32), intent(out), allocatable :: d(:, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(2)
        integer :: e
        call get_dims(fid, path, 2, dims, e)
        if (bail(e, 'rd_i32_2d', error)) return
        allocate (d(dims(1), dims(2)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i32_2o', error)) return
        call h5dread_f(did, h5kind_to_type(i32, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i32_2', error)
    end subroutine rd_i32_2

    subroutine rd_i32_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i32), intent(out), allocatable :: d(:, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(3)
        integer :: e
        call get_dims(fid, path, 3, dims, e)
        if (bail(e, 'rd_i32_3d', error)) return
        allocate (d(dims(1), dims(2), dims(3)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i32_3o', error)) return
        call h5dread_f(did, h5kind_to_type(i32, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i32_3', error)
    end subroutine rd_i32_3

    subroutine rd_i32_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i32), intent(out), allocatable :: d(:, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(4)
        integer :: e
        call get_dims(fid, path, 4, dims, e)
        if (bail(e, 'rd_i32_4d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i32_4o', error)) return
        call h5dread_f(did, h5kind_to_type(i32, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i32_4', error)
    end subroutine rd_i32_4

    subroutine rd_i32_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i32), intent(out), allocatable :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(5)
        integer :: e
        call get_dims(fid, path, 5, dims, e)
        if (bail(e, 'rd_i32_5d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4), dims(5)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i32_5o', error)) return
        call h5dread_f(did, h5kind_to_type(i32, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i32_5', error)
    end subroutine rd_i32_5

    subroutine rd_i32_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i32), intent(out), allocatable :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(6)
        integer :: e
        call get_dims(fid, path, 6, dims, e)
        if (bail(e, 'rd_i32_6d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i32_6o', error)) return
        call h5dread_f(did, h5kind_to_type(i32, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i32_6', error)
    end subroutine rd_i32_6

    subroutine rd_i64_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i64), intent(out) :: d
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(1)
        integer :: e
        dims = [1_hsize_t]
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i64_0o', error)) return
        call h5dread_f(did, h5kind_to_type(i64, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i64_0', error)
    end subroutine rd_i64_0

    subroutine rd_i64_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i64), intent(out), allocatable :: d(:)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(1)
        integer :: e
        call get_dims(fid, path, 1, dims, e)
        if (bail(e, 'rd_i64_1d', error)) return
        allocate (d(dims(1)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i64_1o', error)) return
        call h5dread_f(did, h5kind_to_type(i64, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i64_1', error)
    end subroutine rd_i64_1

    subroutine rd_i64_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i64), intent(out), allocatable :: d(:, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(2)
        integer :: e
        call get_dims(fid, path, 2, dims, e)
        if (bail(e, 'rd_i64_2d', error)) return
        allocate (d(dims(1), dims(2)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i64_2o', error)) return
        call h5dread_f(did, h5kind_to_type(i64, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i64_2', error)
    end subroutine rd_i64_2

    subroutine rd_i64_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i64), intent(out), allocatable :: d(:, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(3)
        integer :: e
        call get_dims(fid, path, 3, dims, e)
        if (bail(e, 'rd_i64_3d', error)) return
        allocate (d(dims(1), dims(2), dims(3)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i64_3o', error)) return
        call h5dread_f(did, h5kind_to_type(i64, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i64_3', error)
    end subroutine rd_i64_3

    subroutine rd_i64_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i64), intent(out), allocatable :: d(:, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(4)
        integer :: e
        call get_dims(fid, path, 4, dims, e)
        if (bail(e, 'rd_i64_4d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i64_4o', error)) return
        call h5dread_f(did, h5kind_to_type(i64, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i64_4', error)
    end subroutine rd_i64_4

    subroutine rd_i64_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i64), intent(out), allocatable :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(5)
        integer :: e
        call get_dims(fid, path, 5, dims, e)
        if (bail(e, 'rd_i64_5d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4), dims(5)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i64_5o', error)) return
        call h5dread_f(did, h5kind_to_type(i64, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i64_5', error)
    end subroutine rd_i64_5

    subroutine rd_i64_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer(i64), intent(out), allocatable :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(6)
        integer :: e
        call get_dims(fid, path, 6, dims, e)
        if (bail(e, 'rd_i64_6d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_i64_6o', error)) return
        call h5dread_f(did, h5kind_to_type(i64, H5_INTEGER_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_i64_6', error)
    end subroutine rd_i64_6

    subroutine rd_r32_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(out) :: d
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(1)
        integer :: e
        dims = [1_hsize_t]
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_r32_0o', error)) return
        call h5dread_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_r32_0', error)
    end subroutine rd_r32_0

    subroutine rd_r32_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(out), allocatable :: d(:)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(1)
        integer :: e
        call get_dims(fid, path, 1, dims, e)
        if (bail(e, 'rd_r32_1d', error)) return
        allocate (d(dims(1)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_r32_1o', error)) return
        call h5dread_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_r32_1', error)
    end subroutine rd_r32_1

    subroutine rd_r32_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(out), allocatable :: d(:, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(2)
        integer :: e
        call get_dims(fid, path, 2, dims, e)
        if (bail(e, 'rd_r32_2d', error)) return
        allocate (d(dims(1), dims(2)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_r32_2o', error)) return
        call h5dread_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_r32_2', error)
    end subroutine rd_r32_2

    subroutine rd_r32_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(out), allocatable :: d(:, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(3)
        integer :: e
        call get_dims(fid, path, 3, dims, e)
        if (bail(e, 'rd_r32_3d', error)) return
        allocate (d(dims(1), dims(2), dims(3)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_r32_3o', error)) return
        call h5dread_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_r32_3', error)
    end subroutine rd_r32_3

    subroutine rd_r32_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(out), allocatable :: d(:, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(4)
        integer :: e
        call get_dims(fid, path, 4, dims, e)
        if (bail(e, 'rd_r32_4d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_r32_4o', error)) return
        call h5dread_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_r32_4', error)
    end subroutine rd_r32_4

    subroutine rd_r32_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(out), allocatable :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(5)
        integer :: e
        call get_dims(fid, path, 5, dims, e)
        if (bail(e, 'rd_r32_5d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4), dims(5)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_r32_5o', error)) return
        call h5dread_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_r32_5', error)
    end subroutine rd_r32_5

    subroutine rd_r32_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(out), allocatable :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(6)
        integer :: e
        call get_dims(fid, path, 6, dims, e)
        if (bail(e, 'rd_r32_6d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_r32_6o', error)) return
        call h5dread_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_r32_6', error)
    end subroutine rd_r32_6

    subroutine rd_r64_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(out) :: d
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(1)
        integer :: e
        dims = [1_hsize_t]
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_r64_0o', error)) return
        call h5dread_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_r64_0', error)
    end subroutine rd_r64_0

    subroutine rd_r64_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(out), allocatable :: d(:)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(1)
        integer :: e
        call get_dims(fid, path, 1, dims, e)
        if (bail(e, 'rd_r64_1d', error)) return
        allocate (d(dims(1)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_r64_1o', error)) return
        call h5dread_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_r64_1', error)
    end subroutine rd_r64_1

    subroutine rd_r64_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(out), allocatable :: d(:, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(2)
        integer :: e
        call get_dims(fid, path, 2, dims, e)
        if (bail(e, 'rd_r64_2d', error)) return
        allocate (d(dims(1), dims(2)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_r64_2o', error)) return
        call h5dread_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_r64_2', error)
    end subroutine rd_r64_2

    subroutine rd_r64_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(out), allocatable :: d(:, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(3)
        integer :: e
        call get_dims(fid, path, 3, dims, e)
        if (bail(e, 'rd_r64_3d', error)) return
        allocate (d(dims(1), dims(2), dims(3)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_r64_3o', error)) return
        call h5dread_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_r64_3', error)
    end subroutine rd_r64_3

    subroutine rd_r64_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(out), allocatable :: d(:, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(4)
        integer :: e
        call get_dims(fid, path, 4, dims, e)
        if (bail(e, 'rd_r64_4d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_r64_4o', error)) return
        call h5dread_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_r64_4', error)
    end subroutine rd_r64_4

    subroutine rd_r64_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(out), allocatable :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(5)
        integer :: e
        call get_dims(fid, path, 5, dims, e)
        if (bail(e, 'rd_r64_5d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4), dims(5)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_r64_5o', error)) return
        call h5dread_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_r64_5', error)
    end subroutine rd_r64_5

    subroutine rd_r64_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(out), allocatable :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        integer(hid_t) :: did
        integer(hsize_t) :: dims(6)
        integer :: e
        call get_dims(fid, path, 6, dims, e)
        if (bail(e, 'rd_r64_6d', error)) return
        allocate (d(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_r64_6o', error)) return
        call h5dread_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call chk(e, 'rd_r64_6', error)
    end subroutine rd_r64_6

    subroutine rd_c32_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r32), intent(out) :: d
        integer, intent(out), optional :: error
        real(r32), allocatable :: buf(:)
        integer(hid_t) :: did
        integer(hsize_t) :: dims(1)
        integer :: e
        call get_dims(fid, path, 1, dims, e)
        if (bail(e, 'rd_c32_0d', error)) return
        allocate (buf(dims(1)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_c32_0o', error)) return
        call h5dread_f(did, h5kind_to_type(r32, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        if (bail(e, 'rd_c32_0r', error)) return
        d = transfer(buf(1:2), d)
        call chk(0, 'rd_c32_0', error)
    end subroutine rd_c32_0

    subroutine rd_c32_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r32), intent(out), allocatable :: d(:)
        integer, intent(out), optional :: error
        real(r32), allocatable :: buf(:, :)
        integer(hid_t) :: did
        integer(hsize_t) :: dims(2)
        integer :: e
        call get_dims(fid, path, 2, dims, e)
        if (bail(e, 'rd_c32_1d', error)) return
        allocate (buf(dims(1), dims(2)))
        allocate (d(dims(2)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_c32_1o', error)) return
        call h5dread_f(did, h5kind_to_type(r32, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        if (bail(e, 'rd_c32_1r', error)) return
        d = reshape(transfer(buf, d), [dims(2)])
        call chk(0, 'rd_c32_1', error)
    end subroutine rd_c32_1

    subroutine rd_c32_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r32), intent(out), allocatable :: d(:, :)
        integer, intent(out), optional :: error
        real(r32), allocatable :: buf(:, :, :)
        integer(hid_t) :: did
        integer(hsize_t) :: dims(3)
        integer :: e
        call get_dims(fid, path, 3, dims, e)
        if (bail(e, 'rd_c32_2d', error)) return
        allocate (buf(dims(1), dims(2), dims(3)))
        allocate (d(dims(2), dims(3)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_c32_2o', error)) return
        call h5dread_f(did, h5kind_to_type(r32, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        if (bail(e, 'rd_c32_2r', error)) return
        d = reshape(transfer(buf, d), [dims(2), dims(3)])
        call chk(0, 'rd_c32_2', error)
    end subroutine rd_c32_2

    subroutine rd_c32_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r32), intent(out), allocatable :: d(:, :, :)
        integer, intent(out), optional :: error
        real(r32), allocatable :: buf(:, :, :, :)
        integer(hid_t) :: did
        integer(hsize_t) :: dims(4)
        integer :: e
        call get_dims(fid, path, 4, dims, e)
        if (bail(e, 'rd_c32_3d', error)) return
        allocate (buf(dims(1), dims(2), dims(3), dims(4)))
        allocate (d(dims(2), dims(3), dims(4)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_c32_3o', error)) return
        call h5dread_f(did, h5kind_to_type(r32, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        if (bail(e, 'rd_c32_3r', error)) return
        d = reshape(transfer(buf, d), [dims(2), dims(3), dims(4)])
        call chk(0, 'rd_c32_3', error)
    end subroutine rd_c32_3

    subroutine rd_c32_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r32), intent(out), allocatable :: d(:, :, :, :)
        integer, intent(out), optional :: error
        real(r32), allocatable :: buf(:, :, :, :, :)
        integer(hid_t) :: did
        integer(hsize_t) :: dims(5)
        integer :: e
        call get_dims(fid, path, 5, dims, e)
        if (bail(e, 'rd_c32_4d', error)) return
        allocate (buf(dims(1), dims(2), dims(3), dims(4), dims(5)))
        allocate (d(dims(2), dims(3), dims(4), dims(5)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_c32_4o', error)) return
        call h5dread_f(did, h5kind_to_type(r32, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        if (bail(e, 'rd_c32_4r', error)) return
        d = reshape(transfer(buf, d), [dims(2), dims(3), dims(4), dims(5)])
        call chk(0, 'rd_c32_4', error)
    end subroutine rd_c32_4

    subroutine rd_c32_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r32), intent(out), allocatable :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        real(r32), allocatable :: buf(:, :, :, :, :, :)
        integer(hid_t) :: did
        integer(hsize_t) :: dims(6)
        integer :: e
        call get_dims(fid, path, 6, dims, e)
        if (bail(e, 'rd_c32_5d', error)) return
        allocate (buf(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
        allocate (d(dims(2), dims(3), dims(4), dims(5), dims(6)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_c32_5o', error)) return
        call h5dread_f(did, h5kind_to_type(r32, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        if (bail(e, 'rd_c32_5r', error)) return
        d = reshape(transfer(buf, d), [dims(2), dims(3), dims(4), dims(5), dims(6)])
        call chk(0, 'rd_c32_5', error)
    end subroutine rd_c32_5

    subroutine rd_c32_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r32), intent(out), allocatable :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        real(r32), allocatable :: buf(:, :, :, :, :, :, :)
        integer(hid_t) :: did
        integer(hsize_t) :: dims(7)
        integer :: e
        call get_dims(fid, path, 7, dims, e)
        if (bail(e, 'rd_c32_6d', error)) return
        allocate (buf(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
        allocate (d(dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_c32_6o', error)) return
        call h5dread_f(did, h5kind_to_type(r32, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        if (bail(e, 'rd_c32_6r', error)) return
        d = reshape(transfer(buf, d), [dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)])
        call chk(0, 'rd_c32_6', error)
    end subroutine rd_c32_6

    subroutine rd_c64_0(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r64), intent(out) :: d
        integer, intent(out), optional :: error
        real(r64), allocatable :: buf(:)
        integer(hid_t) :: did
        integer(hsize_t) :: dims(1)
        integer :: e
        call get_dims(fid, path, 1, dims, e)
        if (bail(e, 'rd_c64_0d', error)) return
        allocate (buf(dims(1)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_c64_0o', error)) return
        call h5dread_f(did, h5kind_to_type(r64, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        if (bail(e, 'rd_c64_0r', error)) return
        d = transfer(buf(1:2), d)
        call chk(0, 'rd_c64_0', error)
    end subroutine rd_c64_0

    subroutine rd_c64_1(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r64), intent(out), allocatable :: d(:)
        integer, intent(out), optional :: error
        real(r64), allocatable :: buf(:, :)
        integer(hid_t) :: did
        integer(hsize_t) :: dims(2)
        integer :: e
        call get_dims(fid, path, 2, dims, e)
        if (bail(e, 'rd_c64_1d', error)) return
        allocate (buf(dims(1), dims(2)))
        allocate (d(dims(2)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_c64_1o', error)) return
        call h5dread_f(did, h5kind_to_type(r64, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        if (bail(e, 'rd_c64_1r', error)) return
        d = reshape(transfer(buf, d), [dims(2)])
        call chk(0, 'rd_c64_1', error)
    end subroutine rd_c64_1

    subroutine rd_c64_2(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r64), intent(out), allocatable :: d(:, :)
        integer, intent(out), optional :: error
        real(r64), allocatable :: buf(:, :, :)
        integer(hid_t) :: did
        integer(hsize_t) :: dims(3)
        integer :: e
        call get_dims(fid, path, 3, dims, e)
        if (bail(e, 'rd_c64_2d', error)) return
        allocate (buf(dims(1), dims(2), dims(3)))
        allocate (d(dims(2), dims(3)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_c64_2o', error)) return
        call h5dread_f(did, h5kind_to_type(r64, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        if (bail(e, 'rd_c64_2r', error)) return
        d = reshape(transfer(buf, d), [dims(2), dims(3)])
        call chk(0, 'rd_c64_2', error)
    end subroutine rd_c64_2

    subroutine rd_c64_3(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r64), intent(out), allocatable :: d(:, :, :)
        integer, intent(out), optional :: error
        real(r64), allocatable :: buf(:, :, :, :)
        integer(hid_t) :: did
        integer(hsize_t) :: dims(4)
        integer :: e
        call get_dims(fid, path, 4, dims, e)
        if (bail(e, 'rd_c64_3d', error)) return
        allocate (buf(dims(1), dims(2), dims(3), dims(4)))
        allocate (d(dims(2), dims(3), dims(4)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_c64_3o', error)) return
        call h5dread_f(did, h5kind_to_type(r64, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        if (bail(e, 'rd_c64_3r', error)) return
        d = reshape(transfer(buf, d), [dims(2), dims(3), dims(4)])
        call chk(0, 'rd_c64_3', error)
    end subroutine rd_c64_3

    subroutine rd_c64_4(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r64), intent(out), allocatable :: d(:, :, :, :)
        integer, intent(out), optional :: error
        real(r64), allocatable :: buf(:, :, :, :, :)
        integer(hid_t) :: did
        integer(hsize_t) :: dims(5)
        integer :: e
        call get_dims(fid, path, 5, dims, e)
        if (bail(e, 'rd_c64_4d', error)) return
        allocate (buf(dims(1), dims(2), dims(3), dims(4), dims(5)))
        allocate (d(dims(2), dims(3), dims(4), dims(5)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_c64_4o', error)) return
        call h5dread_f(did, h5kind_to_type(r64, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        if (bail(e, 'rd_c64_4r', error)) return
        d = reshape(transfer(buf, d), [dims(2), dims(3), dims(4), dims(5)])
        call chk(0, 'rd_c64_4', error)
    end subroutine rd_c64_4

    subroutine rd_c64_5(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r64), intent(out), allocatable :: d(:, :, :, :, :)
        integer, intent(out), optional :: error
        real(r64), allocatable :: buf(:, :, :, :, :, :)
        integer(hid_t) :: did
        integer(hsize_t) :: dims(6)
        integer :: e
        call get_dims(fid, path, 6, dims, e)
        if (bail(e, 'rd_c64_5d', error)) return
        allocate (buf(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
        allocate (d(dims(2), dims(3), dims(4), dims(5), dims(6)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_c64_5o', error)) return
        call h5dread_f(did, h5kind_to_type(r64, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        if (bail(e, 'rd_c64_5r', error)) return
        d = reshape(transfer(buf, d), [dims(2), dims(3), dims(4), dims(5), dims(6)])
        call chk(0, 'rd_c64_5', error)
    end subroutine rd_c64_5

    subroutine rd_c64_6(fid, path, d, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        complex(r64), intent(out), allocatable :: d(:, :, :, :, :, :)
        integer, intent(out), optional :: error
        real(r64), allocatable :: buf(:, :, :, :, :, :, :)
        integer(hid_t) :: did
        integer(hsize_t) :: dims(7)
        integer :: e
        call get_dims(fid, path, 7, dims, e)
        if (bail(e, 'rd_c64_6d', error)) return
        allocate (buf(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
        allocate (d(dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
        call h5dopen_f(fid, trim(path), did, e)
        if (bail(e, 'rd_c64_6o', error)) return
        call h5dread_f(did, h5kind_to_type(r64, H5_REAL_KIND), buf, dims, e)
        call h5dclose_f(did, e)
        if (bail(e, 'rd_c64_6r', error)) return
        d = reshape(transfer(buf, d), [dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)])
        call chk(0, 'rd_c64_6', error)
    end subroutine rd_c64_6

    subroutine wa_i32(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        integer(i32), intent(in) :: val
        integer, intent(out), optional::error
        integer(hid_t) :: oid, sid, aid
        integer(hsize_t) :: dims(1)
        integer::e
        call ensure_obj_exists(fid, obj, e)
        if (bail(e, "wa_i32_ensure", error)) return
        call open_obj(fid, obj, oid, e)
        if (bail(e, "wa_i32", error)) return
        call h5screate_f(H5S_SCALAR_F, sid, e)
        call h5acreate_f(oid, trim(name), h5kind_to_type(i32, H5_INTEGER_KIND), sid, aid, e)
        dims = [1_hsize_t]
        call h5awrite_f(aid, h5kind_to_type(i32, H5_INTEGER_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5sclose_f(sid, e)
        call h5oclose_f(oid, e)
        call chk(e, "wa_i32", error)
    end subroutine
    subroutine wa_i64(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        integer(i64), intent(in) :: val
        integer, intent(out), optional::error
        integer(hid_t) :: oid, sid, aid
        integer(hsize_t) :: dims(1)
        integer::e
        call ensure_obj_exists(fid, obj, e)
        if (bail(e, "wa_i64_ensure", error)) return
        call open_obj(fid, obj, oid, e)
        if (bail(e, "wa_i64", error)) return
        call h5screate_f(H5S_SCALAR_F, sid, e)
        call h5acreate_f(oid, trim(name), h5kind_to_type(i64, H5_INTEGER_KIND), sid, aid, e)
        dims = [1_hsize_t]
        call h5awrite_f(aid, h5kind_to_type(i64, H5_INTEGER_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5sclose_f(sid, e)
        call h5oclose_f(oid, e)
        call chk(e, "wa_i64", error)
    end subroutine
    subroutine wa_r32(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        real(r32), intent(in) :: val
        integer, intent(out), optional::error
        integer(hid_t) :: oid, sid, aid
        integer(hsize_t) :: dims(1)
        integer::e
        call ensure_obj_exists(fid, obj, e)
        if (bail(e, "wa_r32_ensure", error)) return
        call open_obj(fid, obj, oid, e)
        if (bail(e, "wa_r32", error)) return
        call h5screate_f(H5S_SCALAR_F, sid, e)
        call h5acreate_f(oid, trim(name), h5kind_to_type(r32, H5_REAL_KIND), sid, aid, e)
        dims = [1_hsize_t]
        call h5awrite_f(aid, h5kind_to_type(r32, H5_REAL_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5sclose_f(sid, e)
        call h5oclose_f(oid, e)
        call chk(e, "wa_r32", error)
    end subroutine
    subroutine wa_r64(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        real(r64), intent(in) :: val
        integer, intent(out), optional::error
        integer(hid_t) :: oid, sid, aid
        integer(hsize_t) :: dims(1)
        integer::e
        call ensure_obj_exists(fid, obj, e)
        if (bail(e, "wa_r64_ensure", error)) return
        call open_obj(fid, obj, oid, e)
        if (bail(e, "wa_r64", error)) return
        call h5screate_f(H5S_SCALAR_F, sid, e)
        call h5acreate_f(oid, trim(name), h5kind_to_type(r64, H5_REAL_KIND), sid, aid, e)
        dims = [1_hsize_t]
        call h5awrite_f(aid, h5kind_to_type(r64, H5_REAL_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5sclose_f(sid, e)
        call h5oclose_f(oid, e)
        call chk(e, "wa_r64", error)
    end subroutine
    subroutine wa_str(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name, val
        integer, intent(out), optional::error
        integer(hid_t) :: oid, sid, aid, tid
        integer(hsize_t) :: dims(1)
        integer(size_t) :: slen
        integer::e
        call ensure_obj_exists(fid, obj, e)
        if (bail(e, "wa_str_ensure", error)) return
        call open_obj(fid, obj, oid, e)
        if (bail(e, "wa_str", error)) return
        call h5screate_f(H5S_SCALAR_F, sid, e)
        call h5tcopy_f(H5T_NATIVE_CHARACTER, tid, e)
        slen = max(1_size_t, int(len_trim(val), size_t))
        call h5tset_size_f(tid, slen, e)
        call h5acreate_f(oid, trim(name), tid, sid, aid, e)
        dims = [1_hsize_t]
        call h5awrite_f(aid, tid, trim(val), dims, e)
        call h5aclose_f(aid, e)
        call h5tclose_f(tid, e)
        call h5sclose_f(sid, e)
        call h5oclose_f(oid, e)
        call chk(e, "wa_str", error)
    end subroutine
    subroutine wa_i32_1(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        integer(i32), intent(in) :: val(:)
        integer, intent(out), optional::error
        integer(hid_t) :: oid, sid, aid
        integer(hsize_t) :: dims(1)
        integer::e
        dims = shape(val, kind=hsize_t)
        call ensure_obj_exists(fid, obj, e)
        if (bail(e, "wa_i32_1_ensure", error)) return
        call open_obj(fid, obj, oid, e)
        if (bail(e, "wa_i32_1", error)) return
        call h5screate_simple_f(1, dims, sid, e)
        call h5acreate_f(oid, trim(name), h5kind_to_type(i32, H5_INTEGER_KIND), sid, aid, e)
        call h5awrite_f(aid, h5kind_to_type(i32, H5_INTEGER_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5sclose_f(sid, e)
        call h5oclose_f(oid, e)
        call chk(e, "wa_i32_1", error)
    end subroutine
    subroutine wa_i64_1(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        integer(i64), intent(in) :: val(:)
        integer, intent(out), optional::error
        integer(hid_t) :: oid, sid, aid
        integer(hsize_t) :: dims(1)
        integer::e
        dims = shape(val, kind=hsize_t)
        call ensure_obj_exists(fid, obj, e)
        if (bail(e, "wa_i64_1_ensure", error)) return
        call open_obj(fid, obj, oid, e)
        if (bail(e, "wa_i64_1", error)) return
        call h5screate_simple_f(1, dims, sid, e)
        call h5acreate_f(oid, trim(name), h5kind_to_type(i64, H5_INTEGER_KIND), sid, aid, e)
        call h5awrite_f(aid, h5kind_to_type(i64, H5_INTEGER_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5sclose_f(sid, e)
        call h5oclose_f(oid, e)
        call chk(e, "wa_i64_1", error)
    end subroutine
    subroutine wa_r32_1(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        real(r32), intent(in) :: val(:)
        integer, intent(out), optional::error
        integer(hid_t) :: oid, sid, aid
        integer(hsize_t) :: dims(1)
        integer::e
        dims = shape(val, kind=hsize_t)
        call ensure_obj_exists(fid, obj, e)
        if (bail(e, "wa_r32_1_ensure", error)) return
        call open_obj(fid, obj, oid, e)
        if (bail(e, "wa_r32_1", error)) return
        call h5screate_simple_f(1, dims, sid, e)
        call h5acreate_f(oid, trim(name), h5kind_to_type(r32, H5_REAL_KIND), sid, aid, e)
        call h5awrite_f(aid, h5kind_to_type(r32, H5_REAL_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5sclose_f(sid, e)
        call h5oclose_f(oid, e)
        call chk(e, "wa_r32_1", error)
    end subroutine
    subroutine wa_r64_1(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        real(r64), intent(in) :: val(:)
        integer, intent(out), optional::error
        integer(hid_t) :: oid, sid, aid
        integer(hsize_t) :: dims(1)
        integer::e
        dims = shape(val, kind=hsize_t)
        call ensure_obj_exists(fid, obj, e)
        if (bail(e, "wa_r64_1_ensure", error)) return
        call open_obj(fid, obj, oid, e)
        if (bail(e, "wa_r64_1", error)) return
        call h5screate_simple_f(1, dims, sid, e)
        call h5acreate_f(oid, trim(name), h5kind_to_type(r64, H5_REAL_KIND), sid, aid, e)
        call h5awrite_f(aid, h5kind_to_type(r64, H5_REAL_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5sclose_f(sid, e)
        call h5oclose_f(oid, e)
        call chk(e, "wa_r64_1", error)
    end subroutine

    !================================================================
    ! ATTRIBUTE READ
    !================================================================
    subroutine ra_i32(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        integer(i32), intent(out) :: val
        integer, intent(out), optional::error
        integer(hid_t) :: oid, aid
        integer(hsize_t) :: dims(1)
        integer::e
        dims = [1_hsize_t]
        call open_obj(fid, obj, oid, e)
        if (bail(e, "ra_i32", error)) return
        call h5aopen_f(oid, trim(name), aid, e)
        call h5aread_f(aid, h5kind_to_type(i32, H5_INTEGER_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5oclose_f(oid, e)
        call chk(e, "ra_i32", error)
    end subroutine
    subroutine ra_i64(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        integer(i64), intent(out) :: val
        integer, intent(out), optional::error
        integer(hid_t) :: oid, aid
        integer(hsize_t) :: dims(1)
        integer::e
        dims = [1_hsize_t]
        call open_obj(fid, obj, oid, e)
        if (bail(e, "ra_i64", error)) return
        call h5aopen_f(oid, trim(name), aid, e)
        call h5aread_f(aid, h5kind_to_type(i64, H5_INTEGER_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5oclose_f(oid, e)
        call chk(e, "ra_i64", error)
    end subroutine
    subroutine ra_r32(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        real(r32), intent(out) :: val
        integer, intent(out), optional::error
        integer(hid_t) :: oid, aid
        integer(hsize_t) :: dims(1)
        integer::e
        dims = [1_hsize_t]
        call open_obj(fid, obj, oid, e)
        if (bail(e, "ra_r32", error)) return
        call h5aopen_f(oid, trim(name), aid, e)
        call h5aread_f(aid, h5kind_to_type(r32, H5_REAL_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5oclose_f(oid, e)
        call chk(e, "ra_r32", error)
    end subroutine
    subroutine ra_r64(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        real(r64), intent(out) :: val
        integer, intent(out), optional::error
        integer(hid_t) :: oid, aid
        integer(hsize_t) :: dims(1)
        integer::e
        dims = [1_hsize_t]
        call open_obj(fid, obj, oid, e)
        if (bail(e, "ra_r64", error)) return
        call h5aopen_f(oid, trim(name), aid, e)
        call h5aread_f(aid, h5kind_to_type(r64, H5_REAL_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5oclose_f(oid, e)
        call chk(e, "ra_r64", error)
    end subroutine
    subroutine ra_str(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        character(*), intent(out) :: val
        integer, intent(out), optional::error
        integer(hid_t) :: oid, aid, tid
        integer(hsize_t) :: dims(1)
        integer::e
        dims = [1_hsize_t]
        val = ''
        call open_obj(fid, obj, oid, e)
        if (bail(e, "ra_str", error)) return
        call h5aopen_f(oid, trim(name), aid, e)
        call h5aget_type_f(aid, tid, e)
        call h5aread_f(aid, tid, val, dims, e)
        call h5tclose_f(tid, e)
        call h5aclose_f(aid, e)
        call h5oclose_f(oid, e)
        call chk(e, "ra_str", error)
    end subroutine
    subroutine ra_i32_1(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        integer(i32), intent(out), allocatable::val(:)
        integer, intent(out), optional::error
        integer(hid_t) :: oid, aid, sid
        integer(hsize_t) :: dims(1), md(1)
        integer::e, nd
        call open_obj(fid, obj, oid, e)
        if (bail(e, "ra_i32_1", error)) return
        call h5aopen_f(oid, trim(name), aid, e)
        call h5aget_space_f(aid, sid, e)
        call h5sget_simple_extent_ndims_f(sid, nd, e)
        call h5sget_simple_extent_dims_f(sid, dims, md, e)
        call h5sclose_f(sid, e)
        allocate (val(dims(1)))
        call h5aread_f(aid, h5kind_to_type(i32, H5_INTEGER_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5oclose_f(oid, e)
        call chk(e, "ra_i32_1", error)
    end subroutine
    subroutine ra_i64_1(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        integer(i64), intent(out), allocatable::val(:)
        integer, intent(out), optional::error
        integer(hid_t) :: oid, aid, sid
        integer(hsize_t) :: dims(1), md(1)
        integer::e, nd
        call open_obj(fid, obj, oid, e)
        if (bail(e, "ra_i64_1", error)) return
        call h5aopen_f(oid, trim(name), aid, e)
        call h5aget_space_f(aid, sid, e)
        call h5sget_simple_extent_ndims_f(sid, nd, e)
        call h5sget_simple_extent_dims_f(sid, dims, md, e)
        call h5sclose_f(sid, e)
        allocate (val(dims(1)))
        call h5aread_f(aid, h5kind_to_type(i64, H5_INTEGER_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5oclose_f(oid, e)
        call chk(e, "ra_i64_1", error)
    end subroutine
    subroutine ra_r32_1(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        real(r32), intent(out), allocatable::val(:)
        integer, intent(out), optional::error
        integer(hid_t) :: oid, aid, sid
        integer(hsize_t) :: dims(1), md(1)
        integer::e, nd
        call open_obj(fid, obj, oid, e)
        if (bail(e, "ra_r32_1", error)) return
        call h5aopen_f(oid, trim(name), aid, e)
        call h5aget_space_f(aid, sid, e)
        call h5sget_simple_extent_ndims_f(sid, nd, e)
        call h5sget_simple_extent_dims_f(sid, dims, md, e)
        call h5sclose_f(sid, e)
        allocate (val(dims(1)))
        call h5aread_f(aid, h5kind_to_type(r32, H5_REAL_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5oclose_f(oid, e)
        call chk(e, "ra_r32_1", error)
    end subroutine
    subroutine ra_r64_1(fid, obj, name, val, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: obj, name
        real(r64), intent(out), allocatable::val(:)
        integer, intent(out), optional::error
        integer(hid_t) :: oid, aid, sid
        integer(hsize_t) :: dims(1), md(1)
        integer::e, nd
        call open_obj(fid, obj, oid, e)
        if (bail(e, "ra_r64_1", error)) return
        call h5aopen_f(oid, trim(name), aid, e)
        call h5aget_space_f(aid, sid, e)
        call h5sget_simple_extent_ndims_f(sid, nd, e)
        call h5sget_simple_extent_dims_f(sid, dims, md, e)
        call h5sclose_f(sid, e)
        allocate (val(dims(1)))
        call h5aread_f(aid, h5kind_to_type(r64, H5_REAL_KIND), val, dims, e)
        call h5aclose_f(aid, e)
        call h5oclose_f(oid, e)
        call chk(e, "ra_r64_1", error)
    end subroutine

    !=================================================================
    ! CHUNKED / COMPRESSED WRITE
    ! fh5_write_chunk: writes any real/int array with optional gzip
    ! compression. chunk_dims defaults to the dataset shape if omitted.
    ! Useful for large geophysical arrays read slice-by-slice.
    !=================================================================

    !> Write a 1-D real64 array with chunking and optional gzip compression.
    !> chunk_dims: size of each chunk (defaults to full array if absent).
    !> level: gzip level 0-9 (0=no compression, default=6).
    subroutine fh5_write_chunk_r64_1(fid, path, d, chunk_dims, level, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r64), intent(in) :: d(:)
        integer(hsize_t), intent(in), optional :: chunk_dims(1)
        integer, intent(in), optional :: level
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did, pid
        integer(hsize_t) :: dims(1), cdims(1)
        integer :: e, lvl
        call lib_init(e)
        if (bail(e, "fh5_write_chunk_r64_1", error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, "fh5_write_chunk_r64_1 grp", error)) return
        dims = shape(d, kind=hsize_t)
        cdims = dims
        if (present(chunk_dims)) cdims = chunk_dims
        lvl = 6
        if (present(level)) lvl = level
        call h5screate_simple_f(1, dims, sid, e)
        if (bail(e, "fh5_write_chunk_r64_1 spc", error)) return
        call h5pcreate_f(H5P_DATASET_CREATE_F, pid, e)
        call h5pset_chunk_f(pid, 1, cdims, e)
        if (lvl > 0) call h5pset_deflate_f(pid, lvl, e)
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r64, H5_REAL_KIND), sid, did, e, &
            dcpl_id=pid)
        if (bail(e, "fh5_write_chunk_r64_1 dset", error)) then
            call h5pclose_f(pid, e)
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r64, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5pclose_f(pid, e)
        call h5sclose_f(sid, e)
        call chk(e, "fh5_write_chunk_r64_1", error)
    end subroutine fh5_write_chunk_r64_1

    !> Write a 3-D real32 array with chunking + optional gzip (common for seismic volumes).
    subroutine fh5_write_chunk_r32_3(fid, path, d, chunk_dims, level, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        real(r32), intent(in) :: d(:, :, :)
        integer(hsize_t), intent(in), optional :: chunk_dims(3)
        integer, intent(in), optional :: level
        integer, intent(out), optional :: error
        integer(hid_t) :: sid, did, pid
        integer(hsize_t) :: dims(3), cdims(3)
        integer :: e, lvl
        call lib_init(e)
        if (bail(e, "fh5_write_chunk_r32_3", error)) return
        call ensure_groups(fid, trim(path), e)
        if (bail(e, "fh5_write_chunk_r32_3 grp", error)) return
        dims = shape(d, kind=hsize_t)
        cdims = dims
        if (present(chunk_dims)) cdims = chunk_dims
        lvl = 6
        if (present(level)) lvl = level
        call h5screate_simple_f(3, dims, sid, e)
        if (bail(e, "fh5_write_chunk_r32_3 spc", error)) return
        call h5pcreate_f(H5P_DATASET_CREATE_F, pid, e)
        call h5pset_chunk_f(pid, 3, cdims, e)
        if (lvl > 0) call h5pset_deflate_f(pid, lvl, e)
        call h5dcreate_f(fid, trim(path), h5kind_to_type(r32, H5_REAL_KIND), sid, did, e, &
            dcpl_id=pid)
        if (bail(e, "fh5_write_chunk_r32_3 dset", error)) then
            call h5pclose_f(pid, e)
            call h5sclose_f(sid, e)
            return
        end if
        call h5dwrite_f(did, h5kind_to_type(r32, H5_REAL_KIND), d, dims, e)
        call h5dclose_f(did, e)
        call h5pclose_f(pid, e)
        call h5sclose_f(sid, e)
        call chk(e, "fh5_write_chunk_r32_3", error)
    end subroutine fh5_write_chunk_r32_3

    !=================================================================
    ! DELETE a dataset or group
    !=================================================================
    subroutine fh5_delete(fid, path, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: path
        integer, intent(out), optional :: error
        integer :: e
        call h5ldelete_f(fid, trim(path), e)
        call chk(e, "fh5_delete: "//trim(path), error)
    end subroutine fh5_delete

    !=================================================================
    ! LIST names of all direct children of a group
    !=================================================================
    subroutine fh5_list(fid, grp_path, names, error)
        integer(hid_t), intent(in) :: fid
        character(*), intent(in) :: grp_path
        character(len=256), allocatable, intent(out) :: names(:)
        integer, intent(out), optional :: error
        integer(hid_t) :: gid
        integer :: e, n, i
        character(len=256) :: nm
        integer :: obj_type
        call h5gopen_f(fid, trim(grp_path), gid, e)
        if (bail(e, "fh5_list open", error)) return
        call h5gn_members_f(fid, trim(grp_path), n, e)
        allocate (names(n))
        do i = 0, n - 1
            nm = ''
            call h5gget_obj_info_idx_f(fid, trim(grp_path), i, nm, obj_type, e)
            names(i + 1) = trim(nm)
        end do
        call h5gclose_f(gid, e)
        if (present(error)) error = 0
    end subroutine fh5_list

    !=================================================================
    ! COPY a dataset from one location to another (within same file
    ! or across files opened simultaneously)
    !=================================================================
    subroutine fh5_copy(src_fid, src_path, dst_fid, dst_path, error)
        integer(hid_t), intent(in) :: src_fid, dst_fid
        character(*), intent(in) :: src_path, dst_path
        integer, intent(out), optional :: error
        integer :: e
        ! Strip leading path to get just the object name for h5ocopy_f
        integer :: last_slash !, dst_parent_len
        character(len=512) :: dst_parent, dst_name
        last_slash = 0
        block
            integer :: k
            do k = len_trim(dst_path), 1, -1
                if (dst_path(k:k) == '/') then
                    last_slash = k
                    exit
                end if
            end do
        end block
        if (last_slash <= 1) then
            dst_parent = '/'
            dst_name = trim(dst_path)
        else
            dst_parent = dst_path(1:last_slash - 1)
            dst_name = dst_path(last_slash + 1:)
        end if
        call ensure_groups(dst_fid, trim(dst_path), e)
        if (bail(e, "fh5_copy ensure_groups", error)) return
        call h5ocopy_f(src_fid, trim(src_path), dst_fid, trim(dst_path), e)
        call chk(e, "fh5_copy: "//trim(src_path)//" -> "//trim(dst_path), error)
    end subroutine fh5_copy

end module fh5
