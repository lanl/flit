program test

    use libflit

    real, allocatable :: x(:)
    character(len=16) :: type
    real :: fs, flow, fhigh
    integer :: order

    x = load('./x.txt', count_nonempty_lines('./x.txt'), ascii=.true.)
    ! print *, x

    call get_command_argument(2, type)
    fs = extract_float(type)

    call get_command_argument(3, type)
    flow = extract_float(type)

    call get_command_argument(4, type)
    fhigh = extract_float(type)

    call get_command_argument(5, type)
    order = extract_int4(type)

    call get_command_argument(1, type)

    print *, type, fs, flow, fhigh, order

    x = iir_filt(x, 1.0/fs, type, flow, fhigh, order, zerophase=.true.)
    ! print *, x

    call output_array(x, 'y.txt', ascii=.true.)

end program test
