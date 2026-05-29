
program test

    use libflit

    implicit none

    integer :: i, j, k, l, m, n, cnt
    character(len=256), allocatable :: names(:)

    real, allocatable, dimension(:, : ,:) :: vp

    vp = random(64, 64, 32)

    call fh5_open("test.h5", fh5_fid, mode='w')

    call fh5_write_attr(fh5_fid, '/', 'version', '1.0.0')

    call fh5_write(fh5_fid, "/vp", vp)

    call fh5_write_attr(fh5_fid,"/vp","units","m/s")
    call fh5_write_attr(fh5_fid,"/vp", "d",  [10.0, 20.0, 30.0])
    call fh5_write_attr(fh5_fid,"/vp","origin", [0.0, 10.0, 20.0])


    call fh5_write(fh5_fid, "/vs", vp/sqrt(3.0))
    call fh5_write_attr(fh5_fid,"/vs","units","m/s")
    call fh5_write_attr(fh5_fid,"/vs", "d",  [10.0, 20.0, 30.0])
    call fh5_write_attr(fh5_fid,"/vs","origin", [0.0, 10.0, 20.0])


    call fh5_write(fh5_fid, "/complex", cmplx(vp/sqrt(3.0), vp/sqrt(5.0)))
    call fh5_write_attr(fh5_fid,"/complex","units","complex m/s")
    call fh5_write_attr(fh5_fid,"/complex", "d",  [10.0, 20.0, 30.0])
    call fh5_write_attr(fh5_fid,"/complex","origin", [0.0, 10.0, 20.0])

    call fh5_close(fh5_fid)

    call fh5_open("test.h5", fh5_fid, mode='r')
    call fh5_read(fh5_fid,"/vs", vp)
    print *, "vs(1,1,1)=", vp(1,1,1)
    print *, "vs(64,64,32)=", vp(64,64,32)
    call fh5_print_attr(fh5_fid,"/vs","units")
    call fh5_print_attr(fh5_fid,"/vs","d")
    call fh5_print_attr(fh5_fid,"/vs","origin")
    call fh5_print_tree(fh5_fid)
    call fh5_close(fh5_fid)

    !------------------------------------------------------------
    ! 1. Create file, write all 8 scalar types
    !------------------------------------------------------------
    call fh5_open("demo_v2.h5", fh5_fid, mode='w')
    print '(A)', "=== Writing scalars ==="
    call fh5_write(fh5_fid,"/scalars/i8",  int(127,   fh5_i8))
    call fh5_write(fh5_fid,"/scalars/i16", int(1024,  fh5_i16))
    call fh5_write(fh5_fid,"/scalars/i32", int(999999, fh5_i32))
    call fh5_write(fh5_fid,"/scalars/i64", 1234567890123_fh5_i64)
    call fh5_write(fh5_fid,"/scalars/r32", 3.14_fh5_r32)
    call fh5_write(fh5_fid,"/scalars/r64", 2.718281828459045_fh5_r64)
    call fh5_write(fh5_fid,"/scalars/c32", cmplx(1.0_fh5_c32,-2.0_fh5_c32,fh5_c32))
    call fh5_write(fh5_fid,"/scalars/c64", cmplx(3.0_fh5_c64, 4.0_fh5_c64,fh5_c64))
    print '(A)', "  all 8 types OK"

    !------------------------------------------------------------
    ! 2. String
    !------------------------------------------------------------
    call fh5_write_str(fh5_fid,"/meta/title","fh5 v2 demo")
    call fh5_write_str(fh5_fid,"/meta/version","2.0")

    !------------------------------------------------------------
    ! 3. Rank 1–6 arrays (mix of types)
    !------------------------------------------------------------
    print '(A)', "=== Writing arrays rank 1-6 ==="
    block
        integer(fh5_i8) :: v(8)
        v = int([10,20,30,40,50,60,70,80], fh5_i8)
        call fh5_write(fh5_fid,"/arrays/i8_1d",v)
    end block
    block
        real(fh5_r32) :: m2(4,5)
        do j=1,5
            do i=1,4
                m2(i,j)=real(i*j,fh5_r32)
            end do
        end do
        call fh5_write(fh5_fid,"/arrays/r32_2d",m2)
    end block

    block
        !real(fh5_r64) :: a3(4,4,4)

        double precision :: a3(4, 4, 4)
        do k=1,4
            do j=1,4
                do i=1,4
                    a3(i,j,k)=real(i+j+k,fh5_r64)*0.1d0
                end do
            end do
        end do
        call fh5_write(fh5_fid,"/arrays/r64_3d",a3)





    end block
    block
        complex(fh5_c32) :: a4(3,3,3,2)
        do l=1,2
            do k=1,3
                do j=1,3
                    do i=1,3
                        a4(i,j,k,l)=cmplx(real(i+j,fh5_r32),real(k+l,fh5_r32),fh5_c32)
                    end do
                end do
            end do
        end do
        call fh5_write(fh5_fid,"/arrays/c32_4d",a4)
    end block
    block
        integer(fh5_i16) :: a5(2,2,2,2,3)
        cnt=1
        do m=1,3
            do l=1,2
                do k=1,2
                    do j=1,2
                        do i=1,2
                            a5(i,j,k,l,m)=int(cnt,fh5_i16)
                            cnt=cnt+1
                        end do
                    end do
                end do
            end do
        end do
        call fh5_write(fh5_fid,"/arrays/i16_5d",a5)
    end block
    block
        complex(fh5_c64) :: a6(2,2,2,2,2,2)
        cnt=1
        do n=1,2
            do m=1,2
                do l=1,2
                    do k=1,2
                        do j=1,2
                            do i=1,2
                                a6(i,j,k,l,m,n)=cmplx(real(cnt,fh5_r64),real(-cnt,fh5_r64),fh5_c64)
                                cnt=cnt+1
                            end do
                        end do
                    end do
                end do
            end do
        end do
        call fh5_write(fh5_fid,"/arrays/c64_6d",a6)
    end block
    print '(A)', "  OK"

    !    !------------------------------------------------------------
    !    ! 4. Chunked + compressed write  (new feature)
    !    !------------------------------------------------------------
    !    print '(A)', "=== Chunked/compressed write ==="
    !    block
    !        real(fh5_r32)    :: vp(64,64,32)
    !        integer(hsize_t) :: chunks(3)
    !        do k=1,32
    !            do j=1,64
    !                do i=1,64
    !                    vp(i,j,k)=2000.0_fh5_r32+real(k,fh5_r32)*80.0_fh5_r32
    !                end do
    !            end do
    !        end do
    !        chunks = [64_hsize_t, 64_hsize_t, 1_hsize_t]   ! slab-per-depth
    !        call fh5_write_chunk(fh5_fid,"/model/Vp",vp,chunk_dims=chunks,level=6)
    !        call fh5_write_attr(fh5_fid,"/model/Vp","units","m/s")
    !        call fh5_write_attr(fh5_fid,"/model/Vp","dx",  50.0_fh5_r64)
    !        call fh5_write_attr(fh5_fid,"/model/Vp","origin",[0.0_fh5_r64,0.0_fh5_r64,0.0_fh5_r64])
    !        print '(A)', "  Vp [64x64x32] chunk+gzip OK"
    !    end block

    !------------------------------------------------------------
    ! 5. Complex wavefield with metadata
    !------------------------------------------------------------
    block
        complex(fh5_c64) :: snap(16,16,8)
        do k=1,8
            do j=1,16
                do i=1,16
                    snap(i,j,k)=cmplx(sin(real(i,fh5_r64)*0.2d0), &
                        cos(real(j,fh5_r64)*0.15d0),fh5_c64)
                end do
            end do
        end do
        call fh5_write(fh5_fid,"/wavefield/snap",snap)
        call fh5_write_attr(fh5_fid,"/wavefield/snap","time_step",100_fh5_i32)
        call fh5_write_attr(fh5_fid,"/wavefield/snap","dt",0.001_fh5_r64)
        call fh5_write_attr(fh5_fid,"/wavefield/snap","component","uz")
    end block

    call fh5_close(fh5_fid)

    !------------------------------------------------------------
    ! 6. Tree print (fixed indent for all nesting levels)
    !------------------------------------------------------------
    call fh5_open("demo_v2.h5",fh5_fid,mode='r')
    print '(/,A)', "=== fh5_print_tree ==="
    call fh5_print_tree(fh5_fid)

    !------------------------------------------------------------
    ! 7. fh5_print_attr  (auto type-detect, no need to know type)
    !------------------------------------------------------------
    print '(/,A)', "=== fh5_print_attr ==="
    call fh5_print_attr(fh5_fid,"/model/Vp","units")
    call fh5_print_attr(fh5_fid,"/model/Vp","dx")
    call fh5_print_attr(fh5_fid,"/model/Vp","origin")
    call fh5_print_attr(fh5_fid,"/wavefield/snap","time_step")
    call fh5_print_attr(fh5_fid,"/wavefield/snap","component")

    !------------------------------------------------------------
    ! 8. fh5_list
    !------------------------------------------------------------
    print '(/,A)', "=== fh5_list ==="
    call fh5_list(fh5_fid,"/scalars",names)
    print '(A,I0,A)', "  /scalars has ", size(names), " members:"
    do i=1,size(names)
        print '(A,A)', "    ", trim(names(i))
    end do

    call fh5_list(fh5_fid,"/arrays",names)
    print '(A,I0,A)', "  /arrays has ", size(names), " members:"
    do i=1,size(names)
        print '(A,A)', "    ", trim(names(i))
    end do
    call fh5_close(fh5_fid)

    !------------------------------------------------------------
    ! 9. fh5_delete + fh5_copy
    !------------------------------------------------------------
    print '(/,A)', "=== fh5_delete + fh5_copy ==="
    call fh5_open("demo_v2.h5",fh5_fid,mode='a')
    call fh5_delete(fh5_fid,"/arrays/r32_2d")
    print '(A)', "  deleted /arrays/r32_2d"
    call fh5_copy(fh5_fid,"/arrays/r64_3d",fh5_fid,"/arrays/r64_3d_backup")
    print '(A)', "  copied /arrays/r64_3d -> /arrays/r64_3d_backup"

    ! Verify copy
    block
        real(fh5_r64), allocatable :: orig(:,:,:), bak(:,:,:)
        call fh5_read(fh5_fid,"/arrays/r64_3d",orig)
        call fh5_read(fh5_fid,"/arrays/r64_3d_backup",bak)
        if (all(orig == bak)) then
            print '(A)', "  copy verified: original == backup"
        else
            print '(A)', "  ERROR: copy mismatch!"
        end if
    end block
    call fh5_close(fh5_fid)

    !------------------------------------------------------------
    ! 10. Read-back verification
    !------------------------------------------------------------
    print '(/,A)', "=== Read-back ==="
    call fh5_open("demo_v2.h5",fh5_fid,mode='r')
    block
        integer(fh5_i8)  :: a
        real(fh5_r64)    :: f
        complex(fh5_c32) :: g
        complex(fh5_c64) :: h
        call fh5_read(fh5_fid,"/scalars/i8",a)
        print '(A,I4)',    "  i8  =",a
        call fh5_read(fh5_fid,"/scalars/r64",f)
        print '(A,F20.15)',"  r64 =",f
        call fh5_read(fh5_fid,"/scalars/c32",g)
        print '(A,2F8.3)', "  c32 =",g
        call fh5_read(fh5_fid,"/scalars/c64",h)
        print '(A,2F8.3)', "  c64 =",h
    end block
    block
        complex(fh5_c32), allocatable :: a4(:,:,:,:)
        call fh5_read(fh5_fid,"/arrays/c32_4d",a4)
        print '(A,4I3)', "  c32_4d shape=",shape(a4)
        print '(A,2F7.3)',"  c32_4d(1,1,1,1)=",a4(1,1,1,1)
    end block
    block
        complex(fh5_c64), allocatable :: a6(:,:,:,:,:,:)
        call fh5_read(fh5_fid,"/arrays/c64_6d",a6)
        print '(A,6I3)', "  c64_6d shape=",shape(a6)
        print '(A,2F8.2)',"  c64_6d(2,2,2,2,2,2)=",a6(2,2,2,2,2,2)
    end block
    block
        real(fh5_r32), allocatable :: Vp(:,:,:)
        character(len=64) :: lbl
        call fh5_read(fh5_fid,"/model/Vp",Vp)
        call fh5_dtype(fh5_fid,"/model/Vp",lbl)
        print '(A,A)',    "  Vp dtype =",trim(lbl)
        print '(A,F8.1)', "  Vp(1,1,1)=",Vp(1,1,1)
    end block
    block
        character(len=128) :: title
        call fh5_read_str(fh5_fid,"/meta/title",title)
        print '(A,A)', "  title=",trim(title)
    end block
    call fh5_close(fh5_fid)
    call fh5_finalize()
    print '(/,A)', "=== ALL PASSED ==="

end program
