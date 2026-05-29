
program test

    use libflit

    integer, allocatable, dimension(:) :: w
    
    w = irandom(23, range=[0, 100])
    print *, w
    print *, ''
    do i = 1, 5
		print *, i, split_list(w, 5, i)
	end do

end program test
