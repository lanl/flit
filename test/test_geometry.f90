
program test

    use libflit

    block

        real, allocatable, dimension(:, :) :: ps, p1, p2, p3, p4
        real, allocatable, dimension(:) :: pr1, pr2, pr3
        integer, allocatable, dimension(:) :: sindex
        integer, allocatable, dimension(:) :: unindex
        integer :: np, np1, np2, np3, np4, ic

        type(hdbscan_param) :: p

        np1 = 100
        p1 = zeros(np1, 3)
        p1(:, 1) = rescale(random(np1, dist='normal'), [1.0, 14.0])
        p1(:, 2) = rescale(random(np1, dist='normal'), [-1.0, 3.0])
        p1(:, 3) = rescale(random(np1, dist='normal'), [0.2, 4.0])
        p1 = rotate_point(p1, [0.0, real(30.0*const_deg2rad), real(50.0*const_deg2rad)], [0.0, 0.0, 0.0], 'zyx')

        np2 = 200
        p2 = zeros(np2, 3)
        p2(:, 1) = rescale(random(np2, dist='normal'), [1.0, 11.0]) + 4.0
        p2(:, 2) = rescale(random(np2, dist='normal'), [-1.0, 5.0]) - 2.0
        p2(:, 3) = rescale(random(np2, dist='normal'), [0.5, 6.0]) + 2.0
        p2 = rotate_point(p2, [real(20.0*const_deg2rad), real(-30.0*const_deg2rad), real(-10.0*const_deg2rad)], [0.0, 0.0, 0.0], 'zyx')


        np3 = 300
        p3 = zeros(np3, 3)
        p3(:, 1) = rescale(random(np3, dist='normal'), [2.0, 10.0])
        p3(:, 2) = rescale(random(np3, dist='normal'), [-5.0, -2.0]) - 2.0
        p3(:, 3) = rescale(random(np3, dist='normal'), [-3.0, -2.0]) + 2.0
        p3 = rotate_point(p3, [real(45.0*const_deg2rad), real(10.0*const_deg2rad), real(120.0*const_deg2rad)], [0.0, 0.0, 0.0], 'xyz')


        np4 = 400
        p4 = zeros(np4, 3)
        p4(:, 1) = rescale(random(np4, dist='normal'), [2.0, 10.0])
        p4(:, 2) = rescale(random(np4, dist='normal'), [-5.0, 5.0]) + 2.0
        p4(:, 3) = rescale(random(np4, dist='normal'), [-3.0, 1.0]) - 2.0
        p4 = rotate_point(p4, [real(45.0*const_deg2rad), real(90.0*const_deg2rad), real(20.0*const_deg2rad)], [0.0, 0.0, 0.0], 'yzx')

        !===========================================================================
        np = 200
        p%data = zeros(np1 + np2 + np3 + np4 + np, 3)
        p%data(1:np1, :) = p1
        p%data(np1 + 1:np1 + np2, :) = p2
        p%data(np1 + np2 + 1:np1 + np2 + np3, :) = p3
        p%data(np1 + np2 + np3 + 1:np1 + np2 + np3 + np4, :) = p4
        p%data(:np1 + np2 + np3 + np4, 1) = rescale(p%data(:np1 + np2 + np3 + np4, 1), [20.0, 280.0])
        p%data(:np1 + np2 + np3 + np4, 2) = rescale(p%data(:np1 + np2 + np3 + np4, 2), [50.0, 350.0])
        p%data(:np1 + np2 + np3 + np4, 3) = rescale(p%data(:np1 + np2 + np3 + np4, 3), [60.0, 440.0])
        p%data(np1 + np2 + np3 + np4 + 1:, 1) = random(np, range=[0.0, 300.0])
        p%data(np1 + np2 + np3 + np4 + 1:, 2) = random(np, range=[0.0, 400.0])
        p%data(np1 + np2 + np3 + np4 + 1:, 3) = random(np, range=[0.0, 500.0])

        open(3, file='./test_data_3d.txt')
        do i = 1, size(p%data, 1)
            write(3, *) p%data(i, :), rand()
        end do
        close(3)

        p%n = size(p%data, 1)
        p%nd = 3
        p%min_sample = 25
        p%min_cluster_size = 25

        call hdbscan(p)

        print *, p%ncluster, p%nnoisy

        open(3, file='./cluster3.txt')
        do i = 1, p%n
            write(3, *) p%data(i, :), p%labels(i)*1.0
        end do
        close(3)

        unindex = sort(unique(p%labels))

        open(3, file='./cluster3_multiple_fit.txt')
        do ic = 2, size(unindex)

            sindex = pack(regspace(1, 1, p%n), mask=(p%labels==unindex(ic)))
            np = size(sindex)
            print *, ic, np

            ps = zeros(np, 3)
            ps(:, 1) = p%data(sindex, 1)
            ps(:, 2) = p%data(sindex, 2)
            ps(:, 3) = p%data(sindex, 3)

            call fit_surface(ps(:, 1), ps(:, 2), ps(:, 3), pr1, pr2, pr3, method='polynomial', order=2)

            do i = 1, np
                write(3, *) pr1(i), pr2(i), pr3(i), ic*1.0
            end do

            print *, ic

        end do
        close(3)

    end block

    block

        real, allocatable, dimension(:, :) :: ps
        real, allocatable, dimension(:) :: pr1, pr2
        integer, allocatable, dimension(:) :: sindex
        integer, allocatable, dimension(:) :: unindex
        integer :: i, ic, np

        type(hdbscan_param) :: p

        p%data = load('./multiple_cluster_data.txt', 2309, 2, ascii=.true.)

        !===========================================================================
        p%n = size(p%data, 1)
        p%nd = 2
        p%min_sample = 15
        p%min_cluster_size = 15

        call hdbscan(p)

        open(3, file='./cluster_multiple.txt')
        do i = 1, p%n
            write(3, *) p%data(i, :), p%labels(i)*1.0, p%probs(i)
        end do
        close(3)

        unindex = sort(unique(p%labels))

        open(3, file='./cluster_multiple_fit.txt')
        do ic = 2, size(unindex)

            sindex = pack(regspace(1, 1, p%n), mask=(p%labels==unindex(ic)))
            np = size(sindex)

            ps = zeros(np, 2)
            ps(:, 1) = p%data(sindex, 1)
            ps(:, 2) = p%data(sindex, 2)

            call fit_curve(ps(:, 1), ps(:, 2), pr1, pr2, method='polynomial', order=3)

            do i = 1, np
                write(3, *) pr1(i), pr2(i), 100.0, 1.0
            end do

            print *, ic

        end do
        close(3)

        np = count(p%labels /= 0)

        ! This needs github.com/lanl/pymplot
        call  execute_command_line('x_showgraph -in=cluster_multiple.txt,cluster_multiple_fit.txt -ftype=ascii -n1='//num2str(p%n)//','//num2str(np)//' -ptype=3 -tick1beg=-0.6 -tick2beg=-0.6 -x1beg=-0.6 -x1end=0.6 -tick1d=0.3 -label1=X -x2beg=-0.6 -x2end=0.6 -tick2d=0.3 -mtick1=2 -mtick2=2 -select=1,2,3 -color=hsv -ctruncbeg=0.1 -cmin=0 -cmax='//num2str(p%ncluster + 10.0)//' -marker=o,v -markersizemin=5 -markersizemax=5 -markeredgecolor=gray,gray -out=cluster_multiple.png &')

    end block

    block

        real, allocatable, dimension(:) :: p1, p2
        integer :: n
        type(hdbscan_param) :: p
        character(len=1024) :: opts

        ! Raw data
        n = 400
        p%data = zeros(n, 2)
        p%data(:, 1) = random(n, range=[0.0, 1.0], dist='normal')*2 + 3
        p%data(:, 2) = random(n, range=[0.0, 1.0], dist='normal')*0.5 + 1

        open(3, file='./test_data.txt')
        do i = 1, size(p%data, 1)
            write(3, *) p%data(i, 1:2)
        end do
        close(3)

        call fit_curve(p%data(:, 1), p%data(:, 2), p1, p2, method='polynomial', order=2)
        open(3, file='./test_data_fit_p2.txt')
        do i = 1, size(p%data, 1)
            write(3, *) p1(i), p2(i)
        end do
        close(3)

        call fit_curve(p%data(:, 1), p%data(:, 2), p1, p2, method='lowess', smooth=0.25, order=2)
        open(3, file='./test_data_fit_l1.txt')
        do i = 1, size(p%data, 1)
            write(3, *) p1(i), p2(i)
        end do
        close(3)

        call fit_curve(p%data(:, 1), p%data(:, 2), p1, p2, method='lowess', smooth=0.5, order=2)
        open(3, file='./test_data_fit_l2.txt')
        do i = 1, size(p%data, 1)
            write(3, *) p1(i), p2(i)
        end do
        close(3)

        opts = '-x1beg=-5 -x1end=10 -x2beg=-5 -x2end=10 -ftype=ascii -n1=' &
            //num2str(n)//','//num2str(n)//' -ptype=2 -marker=o,* -markerfacecolor=b,r -markeredgecolor=none,none -markersize=5,5 ' &
            //' -label1=X -label2=Y '

        call  execute_command_line('x_showgraph -in=test_data.txt,test_data_fit_p2.txt '//tidy(opts)//' -out=raw_p2.png &')
        call  execute_command_line('x_showgraph -in=test_data.txt,test_data_fit_l1.txt '//tidy(opts)//' -out=raw_l1.png &')
        call  execute_command_line('x_showgraph -in=test_data.txt,test_data_fit_l2.txt '//tidy(opts)//' -out=raw_l2.png &')

        ! Rotate
        do i = 1, n
            p%data(i, :) = rotate_point(p%data(i, :), real(70.0*const_deg2rad), [0.0, 0.0])
        end do
        open(3, file='./test_data_rotation.txt')
        do i = 1, n
            write(3, *) p%data(i, 1:2)
        end do
        close(3)

        call fit_curve(p%data(:, 1), p%data(:, 2), p1, p2, method='polynomial', order=2)
        open(3, file='./test_data_rotation_fit_p2.txt')
        do i = 1, n
            write(3, *) p1(i), p2(i)
        end do
        close(3)

        call fit_curve(p%data(:, 1), p%data(:, 2), p1, p2, method='lowess', smooth=0.25, order=2)
        open(3, file='./test_data_rotation_fit_l1.txt')
        do i = 1, n
            write(3, *) p1(i), p2(i)
        end do
        close(3)

        call fit_curve(p%data(:, 1), p%data(:, 2), p1, p2, method='lowess', smooth=0.5, order=2)
        open(3, file='./test_data_rotation_fit_l2.txt')
        do i = 1, n
            write(3, *) p1(i), p2(i)
        end do
        close(3)

        call  execute_command_line('x_showgraph -in=test_data_rotation.txt '//tidy(opts)//' -n1='//num2str(n)//' -out=rotation.png &')
        call  execute_command_line('x_showgraph -in=test_data_rotation.txt,test_data_rotation_fit_p2.txt '//tidy(opts)//' -out=rotation_p2.png &')
        call  execute_command_line('x_showgraph -in=test_data_rotation.txt,test_data_rotation_fit_l1.txt '//tidy(opts)//' -out=rotation_l1.png &')
        call  execute_command_line('x_showgraph -in=test_data_rotation.txt,test_data_rotation_fit_l2.txt '//tidy(opts)//' -out=rotation_l2.png &')

    end block

    block

        real, allocatable, dimension(:, :) :: w
        real, allocatable, dimension(:, :) :: ps
        real, allocatable, dimension(:) :: pr1, pr2
        integer, allocatable, dimension(:) :: sindex
        integer, allocatable, dimension(:) :: unindex

        type(hdbscan_param) :: p


        p%data = zeros(1000, 2)

        p%data(1:400, 1) = random(400, range=[0.0, 1.0], dist='normal')*2 +0.7
        p%data(1:400, 2) = random(400, range=[0.0, 1.0], dist='normal')*0.5 + 3

        p%data(401:, 1) = random(600, range=[0.0, 1.0], dist='normal', sigma=2.0)/2 + 7.0
        p%data(401:, 2) = random(600, range=[0.0, 1.0], dist='normal', sigma=2.0)*2

        do i = 1, size(p%data, 1)
            p%data(i, :) = rotate_point(p%data(i, :), real(30.0*const_deg2rad), [0.0, 0.0])
        end do


        !===========================================================================
        p%n = size(p%data, 1)
        p%nd = 2

        open(3, file='./test_data.txt', status='replace')
        do i = 1, p%n
            write(3, *) p%data(i, :)
        end do
        close(3)

        p%min_sample = 35
        p%min_cluster_size = 35

        call hdbscan(p)

        open(3, file='./cluster.txt')
        do i = 1, p%n
            write(3, *) p%data(i, :), p%labels(i)*1.0
        end do
        close(3)

        unindex = sort(unique(p%labels))
        sindex = pack(regspace(1, 1, p%n), mask=(p%labels==unindex(2)))
        ps = zeros(size(sindex), 2)
        ps(:, 1) = p%data(sindex, 1)
        ps(:, 2) = p%data(sindex, 2)

        open(3, file='./select.txt')
        do i = 1, size(ps, 1)
            write(3, *) ps(i, 1:2)
        end do
        close(3)

        w = zeros(400, 600)

        !   w = fit_curve(ps(:, 1), ps(:, 2), shape(w), [0.025, 0.025], [-5.0, -5.0], 'lowess', 0.5, 2)
        !   w = fit_curve(ps(:, 1), ps(:, 2), shape(w), [0.025, 0.025], [-5.0, -5.0], 'polynomial', 0.75, 4)
        call fit_curve(ps(:, 1), ps(:, 2), pr1, pr2, w, shape(w), [0.025, 0.025], [-5.0, -5.0], 'lowess', 0.5, 2)
        call output_array(w, './b.bin')

        open(3, file='./select_r.txt')
        do i = 1, size(pr1)
            write(3, *) pr1(i), pr2(i)
        end do
        close(3)

        ! This needs github.com/lanl/pymplot
        call execute_command_line('x_showmatrix -in=b.bin -color=bwr -n1=400 -d1=0.025 -o1=-5 -d2=0.025 -o2=-5 -curve=select.txt,select_r.txt -curveselect=1,2 -curvesize=5,10 -curvefacecolor=w,yellow -curvestyle=scattero,scatter* -out=fit.png; vv fit.png ')

    end block

    stop

end program test
