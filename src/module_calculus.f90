!
! © 2024. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by
! Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by
! Triad National Security, LLC, and the U.S. Department of Energy/National
! Nuclear Security Administration. The Government is granted for itself and
! others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare. derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!


module libflit_calculus

    use libflit_array
    use libflit_error

    implicit none

    interface deriv
        module procedure :: differentiate_1d_float
        module procedure :: differentiate_2d_float
        module procedure :: differentiate_3d_float
        module procedure :: differentiate_1d_double
        module procedure :: differentiate_2d_double
        module procedure :: differentiate_3d_double
        module procedure :: differentiate_1d_complex
        module procedure :: differentiate_2d_complex
        module procedure :: differentiate_3d_complex
        module procedure :: differentiate_1d_dcomplex
        module procedure :: differentiate_2d_dcomplex
        module procedure :: differentiate_3d_dcomplex
    end interface deriv

    interface integ
        module procedure :: integrate_1d_float
        module procedure :: integrate_2d_float
        module procedure :: integrate_3d_float
        module procedure :: integrate_1d_double
        module procedure :: integrate_2d_double
        module procedure :: integrate_3d_double
        module procedure :: integrate_1d_complex
        module procedure :: integrate_2d_complex
        module procedure :: integrate_3d_complex
        module procedure :: integrate_1d_dcomplex
        module procedure :: integrate_2d_dcomplex
        module procedure :: integrate_3d_dcomplex
    end interface integ

    interface cumsum
        module procedure :: cumsum_1d_float
        module procedure :: cumsum_2d_float
        module procedure :: cumsum_3d_float
        module procedure :: cumsum_1d_double
        module procedure :: cumsum_2d_double
        module procedure :: cumsum_3d_double
        module procedure :: cumsum_1d_complex
        module procedure :: cumsum_2d_complex
        module procedure :: cumsum_3d_complex
        module procedure :: cumsum_1d_dcomplex
        module procedure :: cumsum_2d_dcomplex
        module procedure :: cumsum_3d_dcomplex
    end interface cumsum

    private
    public :: taylor_fd_coefs
    public :: deriv
    public :: integ
    public :: cumsum

    ! Center finite-difference coefficients
    double precision, dimension(1:3), parameter, public :: &
        fdcoef_1c10 = [ &
        -5.0000000000000000E-01, &
        0.0000000000000000E+00, &
        5.0000000000000000E-01]

    double precision, dimension(1:5), parameter, public :: &
        fdcoef_1c20 = [ &
        8.3333333333333329E-02, &
        -6.6666666666666663E-01, &
        0.0000000000000000E+00, &
        6.6666666666666663E-01, &
        -8.3333333333333329E-02]

    double precision, dimension(1:7), parameter, public :: &
        fdcoef_1c30 = [ &
        -1.6666666666666666E-02, &
        1.4999999999999999E-01, &
        -7.5000000000000000E-01, &
        -7.4014868308343765E-17, &
        7.5000000000000000E-01, &
        -1.5000000000000002E-01, &
        1.6666666666666666E-02]

    double precision, dimension(1:9), parameter, public :: &
        fdcoef_1c40 = [ &
        3.5714285714285718E-03, &
        -3.8095238095238092E-02, &
        1.9999999999999998E-01, &
        -8.0000000000000004E-01, &
        -3.0531133177191805E-16, &
        8.0000000000000016E-01, &
        -2.0000000000000001E-01, &
        3.8095238095238092E-02, &
        -3.5714285714285709E-03]

    double precision, dimension(1:11), parameter, public :: &
        fdcoef_1c50 = [ &
        -7.9365079365079365E-04, &
        9.9206349206349201E-03, &
        -5.9523809523809521E-02, &
        2.3809523809523808E-01, &
        -8.3333333333333337E-01, &
        -8.8817841970012528E-17, &
        8.3333333333333326E-01, &
        -2.3809523809523805E-01, &
        5.9523809523809521E-02, &
        -9.9206349206349201E-03, &
        7.9365079365079365E-04]

    double precision, dimension(1:13), parameter, public :: &
        fdcoef_1c60 = [ &
        1.8037518037518038E-04, &
        -2.5974025974025978E-03, &
        1.7857142857142856E-02, &
        -7.9365079365079361E-02, &
        2.6785714285714285E-01, &
        -8.5714285714285710E-01, &
        1.4802973661668753E-16, &
        8.5714285714285698E-01, &
        -2.6785714285714285E-01, &
        7.9365079365079361E-02, &
        -1.7857142857142860E-02, &
        2.5974025974025978E-03, &
        -1.8037518037518040E-04]

    double precision, dimension(1:15), parameter, public :: &
        fdcoef_1c70 = [ &
        -4.1625041625041618E-05, &
        6.7987567987567988E-04, &
        -5.3030303030303034E-03, &
        2.6515151515151516E-02, &
        -9.7222222222222238E-02, &
        2.9166666666666674E-01, &
        -8.7500000000000000E-01, &
        -3.6478756523397999E-16, &
        8.7500000000000000E-01, &
        -2.9166666666666663E-01, &
        9.7222222222222224E-02, &
        -2.6515151515151516E-02, &
        5.3030303030303034E-03, &
        -6.7987567987567988E-04, &
        4.1625041625041625E-05]

    double precision, dimension(1:17), parameter, public :: &
        fdcoef_1c80 = [ &
        9.7125097125097091E-06, &
        -1.7760017760017757E-04, &
        1.5540015540015540E-03, &
        -8.7024087024087024E-03, &
        3.5353535353535352E-02, &
        -1.1313131313131311E-01, &
        3.1111111111111112E-01, &
        -8.8888888888888884E-01, &
        -3.4694469519536142E-16, &
        8.8888888888888895E-01, &
        -3.1111111111111106E-01, &
        1.1313131313131312E-01, &
        -3.5353535353535352E-02, &
        8.7024087024087024E-03, &
        -1.5540015540015540E-03, &
        1.7760017760017760E-04, &
        -9.7125097125097125E-06]

    ! forward finite-difference coefficients
    double precision, dimension(1:2), parameter, public :: &
        fdcoef_1f10 = [ &
        -1.0000000000000000E+00, &
        1.0000000000000000E+00]

    double precision, dimension(1:3), parameter, public :: &
        fdcoef_1f20 = [ &
        -1.5000000000000000E+00, &
        2.0000000000000000E+00, &
        -5.0000000000000000E-01]

    double precision, dimension(1:4), parameter, public :: &
        fdcoef_1f30 = [ &
        -1.8333333333333333E+00, &
        3.0000000000000000E+00, &
        -1.5000000000000000E+00, &
        3.3333333333333331E-01]

    double precision, dimension(1:5), parameter, public :: &
        fdcoef_1f40 = [ &
        -2.0833333333333330E+00, &
        4.0000000000000000E+00, &
        -3.0000000000000000E+00, &
        1.3333333333333333E+00, &
        -2.5000000000000000E-01]

    double precision, dimension(1:6), parameter, public :: &
        fdcoef_1f50 = [ &
        -2.2833333333333328E+00, &
        5.0000000000000000E+00, &
        -5.0000000000000000E+00, &
        3.3333333333333330E+00, &
        -1.2500000000000000E+00, &
        2.0000000000000001E-01]

    double precision, dimension(1:7), parameter, public :: &
        fdcoef_1f60 = [ &
        -2.4499999999999993E+00, &
        6.0000000000000000E+00, &
        -7.5000000000000000E+00, &
        6.6666666666666670E+00, &
        -3.7500000000000000E+00, &
        1.2000000000000002E+00, &
        -1.6666666666666666E-01]

    double precision, dimension(1:8), parameter, public :: &
        fdcoef_1f70 = [ &
        -2.5928571428571421E+00, &
        7.0000000000000000E+00, &
        -1.0500000000000000E+01, &
        1.1666666666666668E+01, &
        -8.7500000000000000E+00, &
        4.2000000000000011E+00, &
        -1.1666666666666665E+00, &
        1.4285714285714285E-01]

    double precision, dimension(1:9), parameter, public :: &
        fdcoef_1f80 = [ &
        -2.7178571428571421E+00, &
        8.0000000000000000E+00, &
        -1.4000000000000000E+01, &
        1.8666666666666668E+01, &
        -1.7500000000000000E+01, &
        1.1200000000000003E+01, &
        -4.6666666666666661E+00, &
        1.1428571428571428E+00, &
        -1.2500000000000000E-01]

    ! Backward finite-difference coefficients
    double precision, dimension(1:2), parameter, public :: &
        fdcoef_1b10 = [ &
        -1.0000000000000000E+00, &
        1.0000000000000000E+00]

    double precision, dimension(1:3), parameter, public :: &
        fdcoef_1b20 = [ &
        5.0000000000000000E-01, &
        -2.0000000000000000E+00, &
        1.5000000000000000E+00]

    double precision, dimension(1:4), parameter, public :: &
        fdcoef_1b30 = [ &
        -3.3333333333333331E-01, &
        1.5000000000000000E+00, &
        -3.0000000000000000E+00, &
        1.8333333333333333E+00]

    double precision, dimension(1:5), parameter, public :: &
        fdcoef_1b40 = [ &
        2.5000000000000000E-01, &
        -1.3333333333333333E+00, &
        3.0000000000000000E+00, &
        -4.0000000000000000E+00, &
        2.0833333333333330E+00]

    double precision, dimension(1:6), parameter, public :: &
        fdcoef_1b50 = [ &
        -2.0000000000000001E-01, &
        1.2500000000000000E+00, &
        -3.3333333333333335E+00, &
        5.0000000000000000E+00, &
        -5.0000000000000000E+00, &
        2.2833333333333332E+00]

    double precision, dimension(1:7), parameter, public :: &
        fdcoef_1b60 = [ &
        1.6666666666666666E-01, &
        -1.2000000000000000E+00, &
        3.7500000000000000E+00, &
        -6.6666666666666670E+00, &
        7.5000000000000000E+00, &
        -6.0000000000000000E+00, &
        2.4500000000000002E+00]

    double precision, dimension(1:8), parameter, public :: &
        fdcoef_1b70 = [ &
        -1.4285714285714285E-01, &
        1.1666666666666667E+00, &
        -4.2000000000000002E+00, &
        8.7500000000000000E+00, &
        -1.1666666666666666E+01, &
        1.0500000000000000E+01, &
        -7.0000000000000000E+00, &
        2.5928571428571425E+00]

    double precision, dimension(1:9), parameter, public :: &
        fdcoef_1b80 = [ &
        1.2500000000000000E-01, &
        -1.1428571428571428E+00, &
        4.6666666666666670E+00, &
        -1.1199999999999999E+01, &
        1.7500000000000000E+01, &
        -1.8666666666666668E+01, &
        1.4000000000000000E+01, &
        -8.0000000000000000E+00, &
        2.7178571428571425E+00]

    ! Staggered-grid central finite-difference coefficients
    double precision, dimension(1:2), parameter, public :: &
        fdcoef_1s10 = [ &
        -1.0000000000000000E+00, &
        1.0000000000000000E+00]

    double precision, dimension(1:4), parameter, public :: &
        fdcoef_1s20 = [ &
        4.1666666666666664E-02, &
        -1.1250000000000000E+00, &
        1.1250000000000000E+00, &
        -4.1666666666666664E-02]

    double precision, dimension(1:6), parameter, public :: &
        fdcoef_1s30 = [ &
        -4.6874999999999998E-03, &
        6.5104166666666657E-02, &
        -1.1718750000000000E+00, &
        1.1718750000000000E+00, &
        -6.5104166666666644E-02, &
        4.6874999999999972E-03]

    double precision, dimension(1:8), parameter, public :: &
        fdcoef_1s40 = [ &
        6.9754464285714287E-04, &
        -9.5703124999999990E-03, &
        7.9752604166666671E-02, &
        -1.1962890625000000E+00, &
        1.1962890625000000E+00, &
        -7.9752604166666671E-02, &
        9.5703124999999972E-03, &
        -6.9754464285714233E-04]

    double precision, dimension(1:10), parameter, public :: &
        fdcoef_1s50 = [ &
        -1.1867947048611101E-04, &
        1.7656598772321430E-03, &
        -1.3842773437499999E-02, &
        8.9721679687500042E-02, &
        -1.2112426757812500E+00, &
        1.2112426757812500E+00, &
        -8.9721679687500000E-02, &
        1.3842773437500004E-02, &
        -1.7656598772321434E-03, &
        1.1867947048611115E-04]

    double precision, dimension(1:12), parameter, public :: &
        fdcoef_1s60 = [ &
        2.1847811612215930E-05, &
        -3.5900539822048704E-04, &
        2.9672895159040180E-03, &
        -1.7447662353515622E-02, &
        9.6931457519531278E-02, &
        -1.2213363647460938E+00, &
        1.2213363647460938E+00, &
        -9.6931457519531167E-02, &
        1.7447662353515608E-02, &
        -2.9672895159040119E-03, &
        3.5900539822048536E-04, &
        -2.1847811612215849E-05]

    double precision, dimension(1:14), parameter, public :: &
        fdcoef_1s70 = [ &
        -4.2365147517277731E-06, &
        7.6922503384676844E-05, &
        -6.8945354885525283E-04, &
        4.1789327348981577E-03, &
        -2.0476770401000974E-02, &
        1.0238385200500491E-01, &
        -1.2286062240600586E+00, &
        1.2286062240600586E+00, &
        -1.0238385200500486E-01, &
        2.0476770401000967E-02, &
        -4.1789327348981525E-03, &
        6.8945354885525077E-04, &
        -7.6922503384676709E-05, &
        4.2365147517277553E-06]

    double precision, dimension(1:16), parameter, public :: &
        fdcoef_1s80 = [ &
        8.5234642028808598E-07, &
        -1.7021711056049085E-05, &
        1.6641887751492568E-04, &
        -1.0772711700863316E-03, &
        5.3423855985913959E-03, &
        -2.3036366701126087E-02, &
        1.0664984583854678E-01, &
        -1.2340910732746124E+00, &
        1.2340910732746126E+00, &
        -1.0664984583854666E-01, &
        2.3036366701126097E-02, &
        -5.3423855985913968E-03, &
        1.0772711700863297E-03, &
        -1.6641887751492573E-04, &
        1.7021711056049027E-05, &
        -8.5234642028808460E-07]

    ! Center FD, 2nd order
    double precision, dimension(1:3), parameter, public :: &
        fdcoef_2c10 = [ &
        1.0000000000000000E+00, &
        -2.0000000000000000E+00, &
        1.0000000000000000E+00]

    double precision, dimension(1:5), parameter, public :: &
        fdcoef_2c20 = [ &
        -8.3333333333333329E-02, &
        1.3333333333333333E+00, &
        -2.5000000000000000E+00, &
        1.3333333333333335E+00, &
        -8.3333333333333343E-02]

    double precision, dimension(1:7), parameter, public :: &
        fdcoef_2c30 = [ &
        1.1111111111111108E-02, &
        -1.4999999999999999E-01, &
        1.5000000000000000E+00, &
        -2.7222222222222219E+00, &
        1.5000000000000000E+00, &
        -1.4999999999999999E-01, &
        1.1111111111111108E-02]

    double precision, dimension(1:9), parameter, public :: &
        fdcoef_2c40 = [ &
        -1.7857142857142863E-03, &
        2.5396825396825379E-02, &
        -1.9999999999999996E-01, &
        1.6000000000000001E+00, &
        -2.8472222222222214E+00, &
        1.5999999999999999E+00, &
        -1.9999999999999993E-01, &
        2.5396825396825373E-02, &
        -1.7857142857142833E-03]

    double precision, dimension(1:11), parameter, public :: &
        fdcoef_2c50 = [ &
        3.1746031746031789E-04, &
        -4.9603174603174618E-03, &
        3.9682539682539694E-02, &
        -2.3809523809523814E-01, &
        1.6666666666666670E+00, &
        -2.9272222222222219E+00, &
        1.6666666666666665E+00, &
        -2.3809523809523800E-01, &
        3.9682539682539666E-02, &
        -4.9603174603174566E-03, &
        3.1746031746031724E-04]

    double precision, dimension(1:13), parameter, public :: &
        fdcoef_2c60 = [ &
        -6.0125060125060141E-05, &
        1.0389610389610368E-03, &
        -8.9285714285714246E-03, &
        5.2910052910052983E-02, &
        -2.6785714285714285E-01, &
        1.7142857142857140E+00, &
        -2.9827777777777786E+00, &
        1.7142857142857142E+00, &
        -2.6785714285714296E-01, &
        5.2910052910052942E-02, &
        -8.9285714285714420E-03, &
        1.0389610389610407E-03, &
        -6.0125060125060222E-05]

    double precision, dimension(1:15), parameter, public :: &
        fdcoef_2c70 = [ &
        1.1892869035726177E-05, &
        -2.2662522662522704E-04, &
        2.1212121212121227E-03, &
        -1.3257575757575756E-02, &
        6.4814814814814797E-02, &
        -2.9166666666666669E-01, &
        1.7500000000000000E+00, &
        -3.0235941043083892E+00, &
        1.7499999999999998E+00, &
        -2.9166666666666663E-01, &
        6.4814814814814728E-02, &
        -1.3257575757575739E-02, &
        2.1212121212121175E-03, &
        -2.2662522662522606E-04, &
        1.1892869035726143E-05]

    double precision, dimension(1:17), parameter, public :: &
        fdcoef_2c80 = [ &
        -2.4281274281274256E-06, &
        5.0742907885765017E-05, &
        -5.1800051800051999E-04, &
        3.4809634809634810E-03, &
        -1.7676767676767680E-02, &
        7.5420875420875402E-02, &
        -3.1111111111111123E-01, &
        1.7777777777777781E+00, &
        -3.0548441043083896E+00, &
        1.7777777777777768E+00, &
        -3.1111111111111084E-01, &
        7.5420875420875347E-02, &
        -1.7676767676767652E-02, &
        3.4809634809634736E-03, &
        -5.1800051800051663E-04, &
        5.0742907885764868E-05, &
        -2.4281274281274197E-06]

    ! Forward FD, 2nd order
    double precision, dimension(1:3), parameter, public :: &
        fdcoef_2f10 = [ &
        1.0000000000000000E+00, &
        -2.0000000000000000E+00, &
        1.0000000000000000E+00]

    double precision, dimension(1:4), parameter, public :: &
        fdcoef_2f20 = [ &
        2.0000000000000000E+00, &
        -5.0000000000000000E+00, &
        4.0000000000000000E+00, &
        -1.0000000000000000E+00]

    double precision, dimension(1:5), parameter, public :: &
        fdcoef_2f30 = [ &
        2.9166666666666665E+00, &
        -8.6666666666666661E+00, &
        9.5000000000000000E+00, &
        -4.6666666666666670E+00, &
        9.1666666666666663E-01]

    double precision, dimension(1:6), parameter, public :: &
        fdcoef_2f40 = [ &
        3.7500000000000000E+00, &
        -1.2833333333333332E+01, &
        1.7833333333333332E+01, &
        -1.3000000000000002E+01, &
        5.0833333333333330E+00, &
        -8.3333333333333326E-01]

    double precision, dimension(1:7), parameter, public :: &
        fdcoef_2f50 = [ &
        4.5111111111111111E+00, &
        -1.7399999999999999E+01, &
        2.9250000000000000E+01, &
        -2.8222222222222229E+01, &
        1.6500000000000000E+01, &
        -5.4000000000000004E+00, &
        7.6111111111111107E-01]

    double precision, dimension(1:8), parameter, public :: &
        fdcoef_2f60 = [ &
        5.2111111111111104E+00, &
        -2.2299999999999997E+01, &
        4.3950000000000003E+01, &
        -5.2722222222222236E+01, &
        4.1000000000000000E+01, &
        -2.0100000000000001E+01, &
        5.6611111111111105E+00, &
        -6.9999999999999984E-01]

    double precision, dimension(1:9), parameter, public :: &
        fdcoef_2f70 = [ &
        5.8593253968253958E+00, &
        -2.7485714285714284E+01, &
        6.2100000000000001E+01, &
        -8.9022222222222240E+01, &
        8.6375000000000000E+01, &
        -5.6400000000000006E+01, &
        2.3811111111111110E+01, &
        -5.8857142857142843E+00, &
        6.4821428571428552E-01]

    double precision, dimension(1:10), parameter, public :: &
        fdcoef_2f80 = [ &
        6.4632936507936494E+00, &
        -3.2921428571428571E+01, &
        8.3842857142857142E+01, &
        -1.3975555555555559E+02, &
        1.6247499999999999E+02, &
        -1.3250000000000000E+02, &
        7.4544444444444437E+01, &
        -2.7628571428571423E+01, &
        6.0839285714285696E+00, &
        -6.0396825396825371E-01]

    ! Backward FD, 2nd order
    double precision, dimension(1:3), parameter, public :: &
        fdcoef_2b10 = [ &
        1.0000000000000000E+00, &
        -2.0000000000000000E+00, &
        1.0000000000000000E+00]

    double precision, dimension(1:4), parameter, public :: &
        fdcoef_2b20 = [ &
        -1.0000000000000000E+00, &
        4.0000000000000000E+00, &
        -5.0000000000000000E+00, &
        2.0000000000000000E+00]

    double precision, dimension(1:5), parameter, public :: &
        fdcoef_2b30 = [ &
        9.1666666666666663E-01, &
        -4.6666666666666670E+00, &
        9.5000000000000000E+00, &
        -8.6666666666666661E+00, &
        2.9166666666666665E+00]

    double precision, dimension(1:6), parameter, public :: &
        fdcoef_2b40 = [ &
        -8.3333333333333326E-01, &
        5.0833333333333330E+00, &
        -1.3000000000000000E+01, &
        1.7833333333333332E+01, &
        -1.2833333333333332E+01, &
        3.7500000000000000E+00]

    double precision, dimension(1:7), parameter, public :: &
        fdcoef_2b50 = [ &
        7.6111111111111107E-01, &
        -5.4000000000000004E+00, &
        1.6500000000000000E+01, &
        -2.8222222222222218E+01, &
        2.9250000000000000E+01, &
        -1.7400000000000002E+01, &
        4.5111111111111111E+00]

    double precision, dimension(1:8), parameter, public :: &
        fdcoef_2b60 = [ &
        -6.9999999999999996E-01, &
        5.6611111111111105E+00, &
        -2.0100000000000001E+01, &
        4.1000000000000000E+01, &
        -5.2722222222222221E+01, &
        4.3950000000000003E+01, &
        -2.2300000000000001E+01, &
        5.2111111111111104E+00]

    double precision, dimension(1:9), parameter, public :: &
        fdcoef_2b70 = [ &
        6.4821428571428563E-01, &
        -5.8857142857142852E+00, &
        2.3811111111111106E+01, &
        -5.6399999999999999E+01, &
        8.6375000000000000E+01, &
        -8.9022222222222226E+01, &
        6.2099999999999994E+01, &
        -2.7485714285714280E+01, &
        5.8593253968253958E+00]

    double precision, dimension(1:10), parameter, public :: &
        fdcoef_2b80 = [ &
        -6.0396825396825404E-01, &
        6.0839285714285722E+00, &
        -2.7628571428571430E+01, &
        7.4544444444444451E+01, &
        -1.3250000000000000E+02, &
        1.6247499999999999E+02, &
        -1.3975555555555556E+02, &
        8.3842857142857127E+01, &
        -3.2921428571428564E+01, &
        6.4632936507936494E+00]

    ! Staggered-grid center FD, 2nd order
    double precision, dimension(1:2), parameter, public :: &
        fdcoef_2s10 = [ &
        0.0000000000000000E+00, &
        0.0000000000000000E+00]

    double precision, dimension(1:4), parameter, public :: &
        fdcoef_2s20 = [ &
        5.0000000000000000E-01, &
        -5.0000000000000000E-01, &
        -5.0000000000000000E-01, &
        5.0000000000000000E-01]

    double precision, dimension(1:6), parameter, public :: &
        fdcoef_2s30 = [ &
        -1.0416666666666667E-01, &
        8.1249999999999989E-01, &
        -7.0833333333333337E-01, &
        -7.0833333333333315E-01, &
        8.1249999999999978E-01, &
        -1.0416666666666666E-01]

    double precision, dimension(1:8), parameter, public :: &
        fdcoef_2s40 = [ &
        2.2482638888888889E-02, &
        -2.1657986111111108E-01, &
        1.0148437500000003E+00, &
        -8.2074652777777724E-01, &
        -8.2074652777777768E-01, &
        1.0148437500000000E+00, &
        -2.1657986111111110E-01, &
        2.2482638888888885E-02]

    double precision, dimension(1:10), parameter, public :: &
        fdcoef_2s50 = [ &
        -5.0052703373015886E-03, &
        5.7519531250000006E-02, &
        -3.1668526785714285E-01, &
        1.1549913194444443E+00, &
        -8.9082031250000004E-01, &
        -8.9082031250000038E-01, &
        1.1549913194444446E+00, &
        -3.1668526785714285E-01, &
        5.7519531250000006E-02, &
        -5.0052703373015877E-03]

    double precision, dimension(1:12), parameter, public :: &
        fdcoef_2s60 = [ &
        1.1380537729414682E-03, &
        -1.5247754293774795E-02, &
        9.7351413302951384E-02, &
        -4.0203930082775285E-01, &
        1.2574161590091766E+00, &
        -9.3861857096354184E-01, &
        -9.3861857096354129E-01, &
        1.2574161590091761E+00, &
        -4.0203930082775280E-01, &
        9.7351413302951342E-02, &
        -1.5247754293774795E-02, &
        1.1380537729414678E-03]

    double precision, dimension(1:14), parameter, public :: &
        fdcoef_2s70 = [ &
        -2.6262464060010436E-04, &
        4.0269248195426172E-03, &
        -2.9429484886180437E-02, &
        1.3779560795536747E-01, &
        -4.7426107699278158E-01, &
        1.3354156772674073E+00, &
        -9.7328502352275470E-01, &
        -9.7328502352275603E-01, &
        1.3354156772674082E+00, &
        -4.7426107699278164E-01, &
        1.3779560795536747E-01, &
        -2.9429484886180447E-02, &
        4.0269248195426190E-03, &
        -2.6262464060010458E-04]

    double precision, dimension(1:16), parameter, public :: &
        fdcoef_2s80 = [ &
        6.1269042621576244E-05, &
        -1.0591221946805954E-03, &
        8.7446411014039868E-03, &
        -4.6155933521870758E-02, &
        1.7682398810531150E-01, &
        -5.3559138865697942E-01, &
        1.3967459889316052E+00, &
        -9.9956944280741133E-01, &
        -9.9956944280741311E-01, &
        1.3967459889316054E+00, &
        -5.3559138865697942E-01, &
        1.7682398810531155E-01, &
        -4.6155933521870751E-02, &
        8.7446411014039851E-03, &
        -1.0591221946805957E-03, &
        6.1269042621576230E-05]

contains

    !
    !> Compute weights for Taylor-expansion-based finite difference
    !
    function taylor_fd_coefs(m, x0, x) result(w)

        ! m is the order of the difference
        integer, intent(in) :: m
        ! x0 is the location of the difference to be determined
        real, intent(in) :: x0
        ! x is the locations of the difference coefficients to be determined
        real, dimension(:), intent(in) :: x
        ! w is the array of coefficients
        double precision, allocatable, dimension(:) :: w

        integer :: n, nu, mi, ni
        double precision, allocatable, dimension(:) :: a
        double precision, allocatable, dimension(:, :, :) :: delta
        double precision :: c1, c2, c3

        n = size(x) - 1
        call alloc_array(a, [0, n], source=dble(x(:)))
        call alloc_array(delta, [0, m, 0, n, 0, n])

        delta(0, 0, 0) = 1.0d0
        c1 = 1.0d0
        do ni = 1, n
            c2 = 1.0d0
            do nu = 0, ni - 1
                c3 = a(ni) - a(nu)
                c2 = c2*c3
                do mi = 0, min(ni, m)
                    delta(mi, ni, nu) = ((a(ni) - x0)*delta(mi, ni - 1, nu) &
                        - mi*delta(mi - 1, ni - 1, nu))/c3
                end do
            end do
            do mi = 0, min(ni, m)
                delta(mi, ni, ni) = c1/c2*(mi*delta(mi - 1, ni - 1, ni - 1) &
                    - (a(ni - 1) - x0)*delta(mi, ni - 1, ni - 1))
            end do
            c1 = c2
        end do

        call alloc_array(w, [1, n + 1], source=delta(m, n, :))

    end function taylor_fd_coefs

#define T float
#define TT real
#include "template_calculus.f90"

#define T double
#define TT double precision
#include "template_calculus.f90"

#define T complex
#define TT complex
#include "template_calculus.f90"

#define T dcomplex
#define TT double complex
#include "template_calculus.f90"

end module libflit_calculus
