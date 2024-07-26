
program test

    use libflit

    real, allocatable, dimension(:) :: w
    integer, allocatable, dimension(:) :: id

    w = random(10)
    print *, w
    print *

    w = sort(w)
    print *, w
    print *

    w = random(10)
    print *, w
    print *

    call sort_index(w, id)
    print *, w
    print *, id
    print *

    call sort_index(w, id, -1)
    print *, w
    print *, id

end program test
