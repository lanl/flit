
program test

    use libflit

    character(len=512), allocatable, dimension(:) :: fs
    integer :: i

    print *, get_current_directory()
    print *, get_real_path('./')
    print *, is_file('./w2.bin')
    print *, get_file_directory('./w2.bin')
    print *, count_nonempty_lines('./test.sh')
    call make_directory('./test')
    print *, directory_exists('./test')
    call copy_file('./w2.bin', 'w2r.bin')
    print *, file_exists('./w2r.bin')
    call delete_file('./w2r.bin')
    print *, file_exists('./w2r.bin')
    print *, get_file_size('./w2.bin')
    call move_file('./w2.bin', 'w2m.bin')
    call move_file('./w2m.bin', 'w2.bin')
    call list_files('./', fs)
    do i = 1, size(fs)
        print *, tidy(fs(i))
    end do
    print *, get_extension('./w2.bin')
    print *, get_basename('w2.bin')

end program test
