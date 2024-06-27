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


module libflit_specialfunc

    use libflit_error
    use libflit_array

    implicit none

    interface sinc
        module procedure :: sinc_float
        module procedure :: sinc_double
        module procedure :: sinc_complex
        module procedure :: sinc_dcomplex
    end interface sinc

    interface dsinc
        module procedure :: dsinc_float
        module procedure :: dsinc_double
        module procedure :: dsinc_complex
        module procedure :: dsinc_dcomplex
    end interface dsinc

    interface dirac_delta
        module procedure :: dirac_delta_int
        module procedure :: dirac_delta_float
        module procedure :: dirac_delta_double
    end interface

    interface cgamma
        module procedure :: gamma_complex
        module procedure :: gamma_dcomplex
    end interface cgamma

    interface factorial
        module procedure :: factorial_int
        module procedure :: factorial_float
        module procedure :: factorial_double
    end interface factorial

    interface binomial
        module procedure :: binomial_int
        module procedure :: binomial_float
        module procedure :: binomial_double
    end interface binomial

    private

    ! The following special functions have been made standard in modern Fortran:
    !
    ! Gamma: gamma()
    ! Log(Gamma): log_gamma()
    ! Bessel function of the first kind, order 0: bessel_j0
    ! Bessel function of the first kind, order 1: bessel_j1
    ! Bessel function of the first kind, order n: bessel_jn(n, x)
    ! Bessel function of the second kind, order 0: bessel_y0
    ! Bessel function of the second kind, order 1: bessel_y1
    ! Bessel function of the second kind, order n: bessel_yn(n, x)
    !

    ! Sinc function = sin(x)/x
    public :: sinc

    ! Derivative of sinc function = cos(x)/x - sin(x)/x**2
    public :: dsinc

    ! Delta delta function
    public :: dirac_delta

    ! Modified Bessel function of the first kind of order 0
    public :: bessel_i0

    ! Modified Bessel function of the second kind of order 0
    public :: bessel_k0

    ! Factorial = n!
    public :: factorial

    ! Binomial coefficient C(n, k)
    public :: binomial

    ! Exponential integral: Ei = -int (-x, +inf, x /= 0) exp(-t)/t dt
    public :: ei

    ! Exponential integral: -ei(-x), where x > 0
    public :: expei

    ! Exponential integral: exp(-x)*ei(x), where x /= 0
    public :: eone

    ! Gamma function -- complex version
    public :: cgamma

contains

    !===================================================
    ! Binomial coefficient
    elemental function binomial_int(n, k) result(c)

        integer, intent(in) :: n, k
        real :: c

        c = factorial(n)/(factorial(k)*factorial(n - k))

    end function binomial_int

    elemental function binomial_float(n, k) result(c)

        real, intent(in) :: n, k
        real :: c

        c = gamma(n + 1.0d0)/(gamma(k + 1.0d0)*gamma(n - k + 1.0d0))

    end function binomial_float

    elemental function binomial_double(n, k) result(c)

        double precision, intent(in) :: n, k
        double precision :: c

        c = gamma(n + 1.0d0)/(gamma(k + 1.0d0)*gamma(n - k + 1.0d0))

    end function binomial_double

    !===================================================
    ! Factorial
    elemental function factorial_int(x) result(f)

        integer, intent(in) :: x
        double precision :: f

        f = gamma(x + 1.d0)

    end function factorial_int

    elemental function factorial_float(x) result(f)

        real, intent(in) :: x
        double precision :: f

        f = gamma(x + 1.d0)

    end function factorial_float

    elemental function factorial_double(x) result(f)

        double precision, intent(in) :: x
        double precision :: f

        f = gamma(x + 1.d0)

    end function factorial_double

    !=========================================================
    ! Dirac's delta function
    elemental function dirac_delta_int(x) result(d)

        integer, intent(in) :: x
        integer :: d

        if (x == 0) then
            d = 1
        else
            d = 0
        end if

    end function dirac_delta_int

    elemental function dirac_delta_float(x) result(d)

        real, intent(in) :: x
        real :: d

        if (x == 0) then
            d = 1.0
        else
            d = 0.0
        end if

    end function dirac_delta_float

    elemental function dirac_delta_double(x) result(d)

        double precision, intent(in) :: x
        double precision :: d

        if (x == 0) then
            d = 1.d0
        else
            d = 0.d0
        end if

    end function dirac_delta_double

    !=====================================================
    ! Sinc function
    elemental function sinc_float(x) result(xs)

        real, intent(in) :: x
        real :: xs

        if (x == 0) then
            xs = 1.0
        else
            xs = sin(x)/x
        end if

    end function sinc_float

    elemental function sinc_double(x) result(xs)

        double precision, intent(in) :: x
        double precision :: xs

        if (x == 0) then
            xs = 1.0
        else
            xs = sin(x)/x
        end if

    end function sinc_double

    elemental function sinc_complex(x) result(xs)

        complex, intent(in) :: x
        complex :: xs

        if (x == cmplx(0, 0)) then
            xs = 1.0
        else
            xs = sin(x)/x
        end if

    end function sinc_complex

    elemental function sinc_dcomplex(x) result(xs)

        double complex, intent(in) :: x
        double complex :: xs

        if (x == dcmplx(0, 0)) then
            xs = 1.0
        else
            xs = sin(x)/x
        end if

    end function sinc_dcomplex

    !=====================================================
    ! Derivative of since function
    elemental function dsinc_float(x) result(xs)

        real, intent(in) :: x
        real :: xs

        if (x == 0) then
            xs = 0.0
        else
            xs = cos(x)/x - sin(x)/x**2
        end if

    end function dsinc_float

    elemental function dsinc_double(x) result(xs)

        double precision, intent(in) :: x
        double precision :: xs

        if (x == 0) then
            xs = 0.0
        else
            xs = cos(x)/x - sin(x)/x**2
        end if

    end function dsinc_double

    elemental function dsinc_complex(x) result(xs)

        complex, intent(in) :: x
        complex :: xs

        if (x == cmplx(0, 0)) then
            xs = 0.0
        else
            xs = cos(x)/x - sin(x)/x**2
        end if

    end function dsinc_complex

    elemental function dsinc_dcomplex(x) result(xs)

        double complex, intent(in) :: x
        double complex :: xs

        if (x == dcmplx(0, 0)) then
            xs = 0.0
        else
            xs = cos(x)/x - sin(x)/x**2
        end if

    end function dsinc_dcomplex

    !=====================================================
    !> Compute the exponential integrals ei(x), e1(x), and  exp(-x)*ei(x) for real arguments x
    !
    !>  Originally written by Alan Miller
    !>  Modified by W. J. Cody and Laura Stoltz, Argonne National Laboratory
    !
    subroutine compute_exponetial_integral(arg, result, int)

        double precision, intent(in) :: arg
        double precision, intent(out) :: result
        integer, intent(in) :: int

        integer :: i
        double precision :: ei, frac, px(10), qx(10), sump, sumq, t, w, x, xmx0, y, ysq
        !  mathematical constants
        !   exp40 = exp(40)
        !   x0 = zero of ei
        !   x01/x11 + x02 = zero of ei to extra precision
        double precision, parameter :: zero = 0.0d0, p037 = 0.037d0, half = 0.5d0, &
            one = 1.0d0, two = 2.0d0, three = 3.0d0, &
            four = 4.0d0, six = 6.0d0, twelve = 12.0d0, &
            two4 = 24.0d0, fourty = 40.0d0, &
            exp40 = 2.3538526683701998541d17, &
            x01 = 381.5d0, x11 = 1024.0d0, &
            x02 = -5.1182968633365538008d-5, &
            x0 = 3.7250741078136663466d-1
        ! machine-dependent constants
        double precision, parameter :: xinf = 1.79d+308, xmax = 716.351d0, xbig = 701.84d0
        ! coefficients  for -1.0 <= x < 0.0
        double precision, parameter :: a(7) = &
            [1.1669552669734461083368d2, 2.1500672908092918123209d3, &
            1.5924175980637303639884d4, 8.9904972007457256553251d4, &
            1.5026059476436982420737d5, -1.4815102102575750838086d5, &
            5.0196785185439843791020d0]
        double precision, parameter :: b(6) = &
            [4.0205465640027706061433d1, 7.5043163907103936624165d2, &
            8.1258035174768735759855d3, 5.2440529172056355429883d4, &
            1.8434070063353677359298d5, 2.5666493484897117319268d5]
        ! coefficients for -4.0 <= x < -1.0
        double precision, parameter :: c(9) = &
            [3.828573121022477169108d-1, 1.107326627786831743809d+1, &
            7.246689782858597021199d+1, 1.700632978311516129328d+2, &
            1.698106763764238382705d+2, 7.633628843705946890896d+1, &
            1.487967702840464066613d+1, 9.999989642347613068437d-1, &
            1.737331760720576030932d-8]
        double precision, parameter :: d(9) = &
            [8.258160008564488034698d-2, 4.344836335509282083360d+0, &
            4.662179610356861756812d+1, 1.775728186717289799677d+2, &
            2.953136335677908517423d+2, 2.342573504717625153053d+2, &
            9.021658450529372642314d+1, 1.587964570758947927903d+1, &
            1.000000000000000000000d0]
        ! coefficients for x < -4.0
        double precision, parameter :: e(10) = &
            [1.3276881505637444622987d+2, 3.5846198743996904308695d+4, &
            1.7283375773777593926828d+5, 2.6181454937205639647381d+5, &
            1.7503273087497081314708d+5, 5.9346841538837119172356d+4, &
            1.0816852399095915622498d+4, 1.0611777263550331766871d03, &
            5.2199632588522572481039d+1, 9.9999999999999999087819d-1]
        double precision, parameter :: f(10) = &
            [3.9147856245556345627078d+4, 2.5989762083608489777411d+5, &
            5.5903756210022864003380d+5, 5.4616842050691155735758d+5, &
            2.7858134710520842139357d+5, 7.9231787945279043698718d+4, &
            1.2842808586627297365998d+4, 1.1635769915320848035459d+3, &
            5.4199632588522559414924d+1, 1.0d0]
        !  coefficients for rational approximation to ln(x/a), |1-x/a| < .1
        double precision, parameter :: plg(4) = &
            [-2.4562334077563243311d+01, 2.3642701335621505212d+02, &
            -5.4989956895857911039d+02, 3.5687548468071500413d+02]
        double precision, parameter :: qlg(4) = &
            [-3.5553900764052419184d+01, 1.9400230218539473193d+02, &
            -3.3442903192607538956d+02, 1.7843774234035750207d+02]
        ! coefficients for  0.0 < x < 6.0,
        !  ratio of chebyshev polynomials
        double precision, parameter :: p(10) = &
            [-1.2963702602474830028590d01, -1.2831220659262000678155d03, &
            -1.4287072500197005777376d04, -1.4299841572091610380064d06, &
            -3.1398660864247265862050d05, -3.5377809694431133484800d08, &
            3.1984354235237738511048d08, -2.5301823984599019348858d10, &
            1.2177698136199594677580d10, -2.0829040666802497120940d11]
        double precision, parameter :: q(10) = &
            [7.6886718750000000000000d01, -5.5648470543369082846819d03, &
            1.9418469440759880361415d05, -4.2648434812177161405483d06, &
            6.4698830956576428587653d07, -7.0108568774215954065376d08, &
            5.4229617984472955011862d09, -2.8986272696554495342658d10, &
            9.8900934262481749439886d10, -8.9673749185755048616855d10]
        ! j-fraction coefficients for 6.0 <= x < 12.0
        double precision, parameter :: r(10) = &
            [-2.645677793077147237806d0, -2.378372882815725244124d0, &
            -2.421106956980653511550d01, 1.052976392459015155422d01, &
            1.945603779539281810439d01, -3.015761863840593359165d01, &
            1.120011024227297451523d01, -3.988850730390541057912d0, &
            9.565134591978630774217d0, 9.981193787537396413219d-1]
        double precision, parameter :: s(9) = &
            [1.598517957704779356479d-4, 4.644185932583286942650d0, &
            3.697412299772985940785d02, -8.791401054875438925029d0, &
            7.608194509086645763123d02, 2.852397548119248700147d01, &
            4.731097187816050252967d02, -2.369210235636181001661d02, &
            1.249884822712447891440d0]
        ! j-fraction coefficients for 12.0 <= x < 24.0
        double precision, parameter :: p1(10) = &
            [-1.647721172463463140042d0, -1.860092121726437582253d01, &
            -1.000641913989284829961d01, -2.105740799548040450394d01, &
            -9.134835699998742552432d-1, -3.323612579343962284333d01, &
            2.495487730402059440626d01, 2.652575818452799819855d01, &
            -1.845086232391278674524d0, 9.999933106160568739091d-1]
        double precision, parameter :: q1(9) = &
            [9.792403599217290296840d01, 6.403800405352415551324d01, &
            5.994932325667407355255d01, 2.538819315630708031713d02, &
            4.429413178337928401161d01, 1.192832423968601006985d03, &
            1.991004470817742470726d02, -1.093556195391091143924d01, &
            1.001533852045342697818d0]
        ! j-fraction coefficients for  x .ge. 24.0
        double precision, parameter :: p2(10) = &
            [1.75338801265465972390d02, -2.23127670777632409550d02, &
            -1.81949664929868906455d01, -2.79798528624305389340d01, &
            -7.63147701620253630855d0, -1.52856623636929636839d01, &
            -7.06810977895029358836d0, -5.00006640413131002475d0, &
            -3.00000000320981265753d0, 1.00000000000000485503d0]
        double precision, parameter :: q2(9) = &
            [3.97845977167414720840d04, 3.97277109100414518365d0, &
            1.37790390235747998793d02, 1.17179220502086455287d02, &
            7.04831847180424675988d01, -1.20187763547154743238d01, &
            -7.99243595776339741065d0, -2.99999894040324959612d0, &
            1.99999999999048104167d0]

        x = arg
        if (x == zero) then
            ei = -xinf
            if (int == 2) then
                ei = -ei
            end if
        else if (x < zero .or. int == 2) then
            ! calculate ei for negative argument or for e1.
            y = abs(x)
            if (y <= one) then
                sump = a(7)*y + a(1)
                sumq = y + b(1)
                do i = 2, 6
                    sump = sump*y + a(i)
                    sumq = sumq*y + b(i)
                end do
                ei = log(y) - sump/sumq
                if (int == 3) then
                    ei = ei*exp(y)
                end if
            else if (y <= four) then
                w = one/y
                sump = c(1)
                sumq = d(1)
                do i = 2, 9
                    sump = sump*w + c(i)
                    sumq = sumq*w + d(i)
                end do
                ei = -sump/sumq
                if (int /= 3) then
                    ei = ei*exp(-y)
                end if
            else
                if (y > xbig .and. int < 3) then
                    ei = zero
                else
                    w = one/y
                    sump = e(1)
                    sumq = f(1)
                    do i = 2, 10
                        sump = sump*w + e(i)
                        sumq = sumq*w + f(i)
                    end do
                    ei = -w*(one - w*sump/sumq)
                    if (int /= 3) then
                        ei = ei*exp(-y)
                    end if
                end if
            end if
            if (int == 2) then
                ei = -ei
            end if
        else if (x < six) then
            !  To improve conditioning, rational approximations are expressed
            !    in terms of chebyshev polynomials for 0 <= x < 6, and in
            !    continued fraction form for larger x.
            t = x + x
            t = t/three - two
            px(1) = zero
            qx(1) = zero
            px(2) = p(1)
            qx(2) = q(1)
            do i = 2, 9
                px(i + 1) = t*px(i) - px(i - 1) + p(i)
                qx(i + 1) = t*qx(i) - qx(i - 1) + q(i)
            end do
            sump = half*t*px(10) - px(9) + p(10)
            sumq = half*t*qx(10) - qx(9) + q(10)
            frac = sump/sumq
            xmx0 = (x - x01/x11) - x02
            if (abs(xmx0) >= p037) then
                ei = log(x/x0) + xmx0*frac
                if (int == 3) then
                    ei = exp(-x)*ei
                end if
            else
                ! Special approximation to  ln(x/x0)  for x close to x0
                y = xmx0/(x + x0)
                ysq = y*y
                sump = plg(1)
                sumq = ysq + qlg(1)
                do i = 2, 4
                    sump = sump*ysq + plg(i)
                    sumq = sumq*ysq + qlg(i)
                end do
                ei = (sump/(sumq*(x + x0)) + frac)*xmx0
                if (int == 3) then
                    ei = exp(-x)*ei
                end if
            end if
        else if (x < twelve) then
            frac = zero
            do i = 1, 9
                frac = s(i)/(r(i) + x + frac)
            end do
            ei = (r(10) + frac)/x
            if (int /= 3) then
                ei = ei*exp(x)
            end if
        else if (x <= two4) then
            frac = zero
            do i = 1, 9
                frac = q1(i)/(p1(i) + x + frac)
            end do
            ei = (p1(10) + frac)/x
            if (int /= 3) then
                ei = ei*exp(x)
            end if
        else
            if (x >= xmax .and. int < 3) then
                ei = xinf
            else
                y = one/x
                frac = zero
                do i = 1, 9
                    frac = q2(i)/(p2(i) + x + frac)
                end do
                frac = p2(10) + frac
                ei = y + y*y*frac
                if (int /= 3) then
                    if (x <= xmax - two4) then
                        ei = ei*exp(x)
                    else
                        ! Calculation reformulated to avoid premature overflow
                        ei = (ei*exp(x - fourty))*exp40
                    end if
                end if
            end if
        end if

        result = ei

    end subroutine compute_exponetial_integral

    function ei(x) result(fn_val)

        double precision, intent(in) :: x
        double precision :: fn_val

        call compute_exponetial_integral(x, fn_val, 1)

    end function ei

    function expei(x) result(fn_val)

        double precision, intent(in) :: x
        double precision :: fn_val

        call compute_exponetial_integral(x, fn_val, 3)

    end function expei

    function eone(x) result(fn_val)

        double precision, intent(in) :: x
        double precision :: fn_val

        call compute_exponetial_integral(x, fn_val, 2)

    end function eone

    !=====================================================
    !
    !> Compute modified Bessel functions of the first kind and order zero Bessel_in(0, x)
    !
    !>  Originally written by Alan Miller
    !>  Modified by W. J. Cody and Laura Stoltz, Argonne National Laboratory
    !
    pure subroutine compute_bessel_i0(arg, result, jint)

        double precision, intent(in) :: arg
        double precision, intent(out) :: result
        integer, intent(in) :: jint

        integer :: i
        double precision :: a, b, sump, sumq, x, xx
        !  Mathematical constants
        double precision, parameter :: one = 1.0d0, one5 = 15.0d0, forty = 40.0d0, &
            exp40 = 2.353852668370199854d17, two25 = 225.0d0, &
            rec15 = 6.6666666666666666666d-2
        !  Machine-dependent constants
        double precision, parameter :: xsmall = 5.55d-17, xinf = 1.79d308, xmax = 713.986d0
        !  Coefficients for xsmall .le. abs(arg) .lt. 15.0
        double precision, parameter :: p(15) = &
            [-5.2487866627945699800d-18, -1.5982226675653184646d-14, &
            -2.6843448573468483278d-11, -3.0517226450451067446d-08, &
            -2.5172644670688975051d-05, -1.5453977791786851041d-02, &
            -7.0935347449210549190d+00, -2.4125195876041896775d+03, &
            -5.9545626019847898221d+05, -1.0313066708737980747d+08, &
            -1.1912746104985237192d+10, -8.4925101247114157499d+11, &
            -3.2940087627407749166d+13, -5.5050369673018427753d+14, &
            -2.2335582639474375249d+15]
        double precision, parameter :: q(5) = &
            [-3.7277560179962773046d+03, 6.5158506418655165707d+06, &
            -6.5626560740833869295d+09, 3.7604188704092954661d+12, &
            -9.7087946179594019126d+14]
        !  Coefficients for 15.0 .le. abs(arg)
        double precision, parameter :: pp(8) = &
            [-3.9843750000000000000d-01, 2.9205384596336793945d+00, &
            -2.4708469169133954315d+00, 4.7914889422856814203d-01, &
            -3.7384991926068969150d-03, -2.6801520353328635310d-03, &
            9.9168777670983678974d-05, -2.1877128189032726730d-06]
        double precision, parameter :: qq(7) = &
            [-3.1446690275135491500d+01, 8.5539563258012929600d+01, &
            -6.0228002066743340583d+01, 1.3982595353892851542d+01, &
            -1.1151759188741312645d+00, 3.2547697594819615062d-02, &
            -5.5194330231005480228d-04]

        x = abs(arg)
        if (x < xsmall) then
            result = one
        else if (x < one5) then
            !  xsmall .le.  abs(arg)  .lt. 15.0
            xx = x*x
            sump = p(1)
            do i = 2, 15
                sump = sump*xx + p(i)
            end do
            xx = xx - two25
            sumq = ((((xx + q(1))*xx + q(2))*xx + q(3))*xx + q(4))*xx + q(5)
            result = sump/sumq
            if (jint == 2) result = result*exp(-x)
        else if (x >= one5) then
            if ((jint == 1) .and. (x > xmax)) then
                result = xinf
            else
                !  15.0  .le.  abs(arg)
                xx = one/x - rec15
                sump = ((((((pp(1)*xx + pp(2))*xx + pp(3))*xx + pp(4))*xx + pp(5))*xx + &
                    pp(6))*xx + pp(7))*xx + pp(8)
                sumq = ((((((xx + qq(1))*xx + qq(2))*xx + qq(3))*xx + qq(4))*xx + &
                    qq(5))*xx + qq(6))*xx + qq(7)
                result = sump/sumq
                if (jint == 2) then
                    result = (result - pp(1))/sqrt(x)
                else
                    !  calculation reformulated to avoid premature overflow
                    if (x <= (xmax - one5)) then
                        a = exp(x)
                        b = one
                    else
                        a = exp(x - forty)
                        b = exp40
                    end if
                    result = ((result*a - pp(1)*a)/sqrt(x))*b
                end if
            end if
        end if

        ! Return for abs(arg) .lt. xsmall

    end subroutine compute_bessel_i0

    elemental function bessel_i0(x) result(fn_val)

        double precision, intent(in) :: x
        double precision :: fn_val

        call compute_bessel_i0(x, fn_val, 1)

    end function bessel_i0

    !=====================================================
    !
    !> Compute modified Bessel functions of the second kind and order zero Bessel_kn(0, x)
    !
    !>  Originally written by Alan Miller
    !>  Modified by W. J. Cody and Laura Stoltz, Argonne National Laboratory
    !
    pure subroutine compute_bessel_k0(arg, result, jint)

        double precision, intent(in) :: arg
        double precision, intent(out) :: result
        integer, intent(in) :: jint

        integer :: i
        double precision :: sumf, sumg, sump, sumq, temp, x, xx
        ! Constants
        double precision, parameter :: one = 1.0d0, zero = 0.0d0
        ! Machine-dependent constants
        double precision, parameter :: xsmall = 1.11d-16, xinf = 1.79d+308, xmax = 705.342d0
        ! Coefficients for xsmall .le.  arg  .le. 1.0
        double precision, parameter :: p(6) = &
            [5.8599221412826100000d-04, 1.3166052564989571850d-01, &
            1.1999463724910714109d+01, 4.6850901201934832188d+02, &
            5.9169059852270512312d+03, 2.4708152720399552679d+03]
        double precision, parameter :: q(2) = &
            [-2.4994418972832303646d+02, 2.1312714303849120380d+04]
        double precision, parameter :: f(4) = &
            [-1.6414452837299064100d+00, -2.9601657892958843866d+02, &
            -1.7733784684952985886d+04, -4.0320340761145482298d+05]
        double precision, parameter :: g(3) = &
            [-2.5064972445877992730d+02, 2.9865713163054025489d+04, &
            -1.6128136304458193998d+06]
        ! Coefficients for  1.0 .lt. arg
        double precision, parameter :: pp(10) = &
            [1.1394980557384778174d+02, 3.6832589957340267940d+03, &
            3.1075408980684392399d+04, 1.0577068948034021957d+05, &
            1.7398867902565686251d+05, 1.5097646353289914539d+05, &
            7.1557062783764037541d+04, 1.8321525870183537725d+04, &
            2.3444738764199315021d+03, 1.1600249425076035558d+02]
        double precision, parameter :: qq(10) = &
            [2.0013443064949242491d+02, 4.4329628889746408858d+03, &
            3.1474655750295278825d+04, 9.7418829762268075784d+04, &
            1.5144644673520157801d+05, 1.2689839587977598727d+05, &
            5.8824616785857027752d+04, 1.4847228371802360957d+04, &
            1.8821890840982713696d+03, 9.2556599177304839811d+01]

        x = arg
        if (x > zero) then
            if (x <= one) then
                ! 0.0 .lt.  arg  .le. 1.0
                temp = log(x)
                if (x < xsmall) then
                    ! return for small arg
                    result = p(6)/q(2) - temp
                else
                    xx = x*x
                    sump = ((((p(1)*xx + p(2))*xx + p(3))*xx + p(4))*xx + p(5))*xx + p(6)
                    sumq = (xx + q(1))*xx + q(2)
                    sumf = ((f(1)*xx + f(2))*xx + f(3))*xx + f(4)
                    sumg = ((xx + g(1))*xx + g(2))*xx + g(3)
                    result = sump/sumq - xx*sumf*temp/sumg - temp
                    if (jint == 2) result = result*exp(x)
                end if
            else if ((jint == 1) .and. (x > xmax)) then
                ! error return for arg .gt. xmax
                result = zero
            else
                ! 1.0 .lt. arg
                xx = one/x
                sump = pp(1)
                do i = 2, 10
                    sump = sump*xx + pp(i)
                end do
                sumq = xx
                do i = 1, 9
                    sumq = (sumq + qq(i))*xx
                end do
                sumq = sumq + qq(10)
                result = sump/sumq/sqrt(x)
                if (jint == 1) result = result*exp(-x)
            end if
        else
            ! error return for arg .le. 0.0
            result = xinf
        end if

    end subroutine compute_bessel_k0

    elemental function bessel_k0(x) result(fn_val)

        double precision, intent(in) :: x
        double precision :: fn_val

        call compute_bessel_k0(x, fn_val, 1)

    end function bessel_k0

    !=====================================================
    !> Compute the complex gamma and loggamma functions, accurate to within 5 in the 14th significant
    !
    !>  Originally written by Alfred H. Morris, Jr., Naval Surface Warfare Center
    !>  Modified by Alan Miller
    !
    pure subroutine cgamma_(mo, z, w)

        integer, intent(in) :: mo
        double complex, intent(in) :: z
        double complex, intent(out) :: w

        double complex :: eta, eta2, sum
        double precision, parameter :: c0(12) = [0.833333333333333d-01, &
            -0.277777777777778d-02, 0.793650793650794d-03, &
            -0.595238095238095d-03, 0.841750841750842d-03, &
            -0.191752691752692d-02, 0.641025641025641d-02, &
            -0.295506535947712d-01, 0.179644372368831d0, &
            -1.39243221690590d0, 13.4028640441684d0, &
            -156.848284626002d0], pi = 3.14159265358979d0, &
            pi2 = 6.28318530717959d0, alpi = 1.14472988584940d0, &
            hl2p = 0.918938533204673d0, half = 0.5d0
        double precision :: a, a1, a2, c, cn, cut, d, eps, et, e2t, h1, h2, s, sn, &
            s1, s2, t, t1, t2, u, u1, u2, v1, v2, w1, w2, x, y, y2
        integer :: j, k, l, m, max, n, nm1

        !---------------------------
        ! alpi = log(pi)
        ! hl2p = 0.5 * log(2*pi)
        !---------------------------

        ! ****** max and eps are machine dependent constants.
        ! max is the largest positive integer that may
        ! be used, and eps is the smallest real number
        ! such that 1.0 + eps > 1.0.

        ! max = ipmpar(3)
        max = huge(3)
        eps = epsilon(1.0)

        !---------------------------
        x = z
        y = aimag(z)
        if (x < 0.0) then
            ! case when the real part of z is negative
            y = abs(y)
            t = -pi*y
            et = exp(t)
            e2t = et*et
            ! set  a1 = (1 + e2t)/2  and  a2 = (1 - e2t)/2
            a1 = half*(1.0 + e2t)
            t2 = t + t
            if (t2 >= -0.15) then
                a2 = -half*rexp(t2)
            else
                a2 = half*(half + (half - e2t))
            end if
            ! compute sin(pi*x) and cos(pi*x)
            if (abs(x) >= min(real(max), 1.0/eps)) goto 70
            k = abs(x)
            u = x + k
            k = mod(k, 2)
            if (u <= -half) then
                u = half + (half + u)
                k = k + 1
            end if
            u = pi*u
            sn = sin(u)
            cn = cos(u)
            if (k == 1) then
                sn = -sn
                cn = -cn
            end if
            ! set  h1 + h2*i  to  pi/sin(pi*z)  or  log(pi/sin(pi*z))
            a1 = sn*a1
            a2 = cn*a2
            a = a1*a1 + a2*a2
            if (a == 0.0) then
                goto 70
            end if
            if (mo == 0) then

                h1 = a1/a
                h2 = -a2/a
                c = pi*et
                h1 = c*h1
                h2 = c*h2
            else

                h1 = (alpi + t) - half*log(a)
                h2 = -atan2(a2, a1)
            end if
            if (aimag(z) >= 0.0) then
                x = 1.0 - x
                y = -y
            else
                h2 = -h2
                x = 1.0 - x
            end if
        end if
        !  case when the real part of z is nonnegative
        w1 = 0.0
        w2 = 0.0
        n = 0
        t = x
        y2 = y*y
        a = t*t + y2
        cut = 36.0
        if (eps > 1.d-8) then
            cut = 16.0
        end if
        if (a < cut) then
            if (a == 0.0) then
                goto 70
            end if
            10 continue
            n = n + 1
            t = t + 1.0
            a = t*t + y2
            if (a < cut) goto 10
            ! let s1 + s2*i be the product of the terms (z+j)/(z+n)
            u1 = (x*t + y2)/a
            u2 = y/a
            s1 = u1
            s2 = n*u2
            if (n >= 2) then
                u = t/a
                nm1 = n - 1
                do j = 1, nm1
                    v1 = u1 + j*u
                    v2 = (n - j)*u2
                    c = s1*v1 - s2*v2
                    d = s1*v2 + s2*v1
                    s1 = c
                    s2 = d
                end do
            end if
            ! set  w1 + w2*i = log(s1 + s2*i)  when mo is nonzero
            s = s1*s1 + s2*s2
            if (mo /= 0) then
                w1 = half*log(s)
                w2 = atan2(s2, s1)
            end if
        end if
        ! set  v1 + v2*i = (z - 0.5) * log(z + n) - z
        t1 = half*log(a) - 1.0
        t2 = atan2(y, t)
        u = x - half
        v1 = (u*t1 - half) - y*t2
        v2 = u*t2 + y*t1
        ! let a1 + a2*i be the asymptotic sum
        eta = dcmplx(t/a, -y/a)
        eta2 = eta*eta
        m = 12
        if (a >= 289.0) then
            m = 6
        end if
        if (eps > 1.d-8) then
            m = m/2
        end if
        sum = dcmplx(c0(m), 0.0)
        l = m
        do j = 2, m
            l = l - 1
            sum = dcmplx(c0(l), 0.0) + sum*eta2
        end do
        sum = sum*eta
        a1 = real(sum)
        a2 = aimag(sum)
        ! gathering together the results
        w1 = (((a1 + hl2p) - w1) + v1) - n
        w2 = (a2 - w2) + v2
        if (real(z) < 0.0) then
            goto 50
        end if
        if (mo == 0) then
            ! case when the real part of z is nonnegative and mo = 0
            a = exp(w1)
            w1 = a*cos(w2)
            w2 = a*sin(w2)
            if (n == 0) then
                goto 60
            end if
            c = (s1*w1 + s2*w2)/s
            d = (s1*w2 - s2*w1)/s
            w1 = c
            w2 = d
            goto 60
        end if
        ! case when the real part of z is nonnegative and mo is nonzero.
        ! the angle w2 is reduced to the interval -pi < w2 <= pi.
        40 continue
        if (w2 <= pi) then
            k = half - w2/pi2
            w2 = w2 + pi2*k
            goto 60
        end if
        k = w2/pi2 - half
        w2 = w2 - pi2*real(k + 1)
        if (w2 <= -pi) then
            w2 = pi
        end if
        goto 60
        ! case when the real part of z is negative and mo is nonzero
        50 continue
        if (mo /= 0) then
            w1 = h1 - w1
            w2 = h2 - w2
            goto 40
        end if
        ! case when the real part of z is negative and mo = 0
        a = exp(-w1)
        t1 = a*cos(-w2)
        t2 = a*sin(-w2)
        w1 = h1*t1 - h2*t2
        w2 = h1*t2 + h2*t1
        if (n /= 0) then
            c = w1*s1 - w2*s2
            d = w1*s2 + w2*s1
            w1 = c
            w2 = d
        end if
        ! termination
        60 continue
        w = dcmplx(w1, w2)
        return
        !  the requested value cannot be computed
        70 continue
        w = (0.0, 0.0)
        return

    contains

        pure function rexp(x) result(fn_val)

            ! evaluation of the function exp(x) - 1
            double precision, intent(in) :: x
            double precision :: fn_val

            double precision, parameter :: p1 = 0.914041914819518d-09, &
                p2 = 0.238082361044469d-01, q1 = -0.499999999085958d0, &
                q2 = 0.107141568980644d0, q3 = -0.119041179760821d-01, &
                q4 = 0.595130811860248d-03
            double precision :: e

            if (abs(x) <= 0.15) then
                fn_val = x*(((p2*x + p1)*x + 1.0)/ &
                    ((((q4*x + q3)*x + q2)*x + q1)*x + 1.0))
                return
            end if

            if (x >= 0.0) then
                e = exp(x)
                fn_val = e*(half + (half - 1.0/e))
                return
            end if

            if (x >= -37.0) then
                fn_val = (exp(x) - half) - half
                return
            end if

            fn_val = -1.0

        end function rexp

    end subroutine cgamma_

    elemental function gamma_complex(x) result(y)

        complex, intent(in) :: x
        complex :: y

        double complex :: yy

        call cgamma_(0, dcmplx(x%re, x%im), yy)
        y = cmplx(yy%re, yy%im)

    end function gamma_complex

    elemental function gamma_dcomplex(x) result(y)

        double complex, intent(in) :: x
        double complex :: y

        call cgamma_(0, x, y)

    end function gamma_dcomplex

end module libflit_specialfunc
