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


!
!> Constants for modeling, inversion and imaging
!
module libflit_constants

    use iso_fortran_env
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

    implicit none

    ! Constant unit imaginary number
    complex, parameter :: const_i = cmplx(0.0, 1.0)

    ! Constant complex zero
    complex, parameter :: const_complex_zero = cmplx(0.0, 0.0)
    complex, parameter :: const_dcomplex_zero = dcmplx(0.0d0, 0.0d0)

    ! Constant Pi
    double precision, parameter :: const_pi = 3.141592653589793238462643
    double precision, parameter :: const_pi_half = 1.570796326794896619231321
    double precision, parameter :: const_pi_inv = 0.318309886183790671537767
    double precision, parameter :: const_pi_sqr = 9.869604401089358618834491
    double precision, parameter :: const_pi_sqrt = 1.772453850905516027298167
    double precision, parameter :: const_pi_ln = 1.144729885849400174143427
    double precision, parameter :: const_pi_log10 = 0.497149872694133854351268

    ! Euler's number e
    double precision, parameter :: const_ee = 2.71828182845904523560287
    double precision, parameter :: const_ee_inv = 0.367879441171442321595523
    double precision, parameter :: const_ee_sqr = 7.389056098930650227230427
    double precision, parameter :: const_ee_log10 = 0.434294481903251827651129

    ! sqrt(2)
    double precision, parameter :: const_sqrt2 = 1.4142135623730951
    double precision, parameter :: const_sqrt2_inv = 0.7071067811865475

    ! sqrt(3)
    double precision, parameter :: const_sqrt3 = 1.7320508075688772
    double precision, parameter :: const_sqrt3_inv = 0.5773502691896258

    ! sqrt(5)
    double precision, parameter :: const_sqrt5 = 2.23606797749979
    double precision, parameter :: const_sqrt5_inv = 0.4472135954999579

    ! degree radian
    double precision, parameter :: const_deg2rad = 0.017453292519943295
    double precision, parameter :: const_rad2deg = 57.29577951308232

    ! Volumes
    double precision, parameter :: const_oz2cbm = 0.0000295703125
    double precision, parameter :: const_cbm2oz = 33817.70145310435931307794

    ! Length
    real, parameter :: const_ft2m = 0.3048
    real, parameter :: const_m2ft = 3.28084

    ! Area
    double precision, parameter :: const_sqft2sqm = 0.09290304
    double precision, parameter :: const_sqm2sqft = 10.76391041670972230833

    ! golden ratio
    double precision, parameter :: const_golden = 1.618033988749894

    ! light in vacuum
    real, parameter :: const_lightc = 299752458.0

    ! Planck constant
    double precision, parameter :: const_planckh = 6.62607004081d-34

    ! Vacuum permitivity
    double precision, parameter :: const_vacuum_permittivity = 8.8541878128d-12

    ! Vacuum permeability
    double precision, parameter :: const_vacuum_permeability = 1.25663706212d-6

    ! Electron charge
    double precision, parameter :: const_elementary_charge = 1.602176634d-19

    ! Fine structure constant
    double precision, parameter :: const_fine_structure = 0.00729735256278713575
    double precision, parameter :: const_fine_structure_reciprocial = 137.035999206

    ! Scales
    double precision, parameter :: const_yotta = 1.0d24
    double precision, parameter :: const_zetta = 1.0d21
    double precision, parameter :: const_exa = 1.0d18
    double precision, parameter :: const_peta = 1.0d15
    double precision, parameter :: const_tera = 1.0d12
    double precision, parameter :: const_giga = 1.0d9
    double precision, parameter :: const_mega = 1.0d6
    double precision, parameter :: const_kilo = 1.0d3
    double precision, parameter :: const_milli = 1.0d-3
    double precision, parameter :: const_micro = 1.0d-6
    double precision, parameter :: const_nano = 1.0d-9
    double precision, parameter :: const_pico = 1.0d-12
    double precision, parameter :: const_femto = 1.0d-15
    double precision, parameter :: const_atto = 1.0d-18
    double precision, parameter :: const_zepto = 1.0d-21
    double precision, parameter :: const_yocto = 1.0d-24

    ! tiny, small and huge value
    real, parameter :: float_tiny = tiny(0.0)
    real, parameter :: float_small = 1.0e-6
    real, parameter :: float_large = 1.0e+6
    real, parameter :: float_huge = huge(0.0)

    double precision, parameter :: double_tiny = tiny(0.0d0)
    double precision, parameter :: double_small = 1.0d-6
    double precision, parameter :: double_large = 1.0d+6
    double precision, parameter :: double_huge = huge(0.0d0)

    integer(4), parameter :: int4_tiny = 0
    integer(4), parameter :: int4_huge = huge(int(0, kind=4))

    integer(2), parameter :: int2_tiny = 0
    integer(2), parameter :: int2_huge = huge(int(0, kind=2))

    integer(1), parameter :: int1_tiny = 0
    integer(1), parameter :: int1_huge = huge(int(0, kind=1))

    !    ! The following definitions don't work... so have to define a nan function which is nan = sqrt(-1.0)
    !    integer, parameter :: int_nan = huge(0) + huge(0)
    !    real, parameter :: float_nan = transfer(-4194304_int32, 1._real32)
    !    double precision, parameter :: double_nan = transfer(-2251799813685248_int64, 1._real64)
    !    complex, parameter :: complex_nan = cmplx(float_nan, float_nan)
    !    double complex, parameter :: dcomplex_nan = dcmplx(double_nan, double_nan)

    character(len=1), parameter :: char_space = ' '
    character(len=1), parameter :: char_newline = achar(10)
    character(len=1), parameter :: char_vrectangle = achar(219)
    character(len=1), parameter :: char_hrectangle = achar(220)
    character(len=1), parameter :: char_square = achar(254)
    character(len=1), parameter :: char_degree = achar(248)
    character(len=1), parameter :: char_inf = achar(236)
    character(len=1), parameter :: char_delta = achar(235)
    character(len=1), parameter :: char_approx = achar(247)
    character(len=1), parameter :: char_ge = achar(242)
    character(len=1), parameter :: char_le = achar(243)
    character(len=1), parameter :: char_pm = achar(241)
    character(len=1), parameter :: char_define = achar(240)
    character(len=1), parameter :: char_pi = achar(227)
    character(len=1), parameter :: char_sum = achar(228)

    ! For reference, in C++, Greek letters can be defined through unicodes.
    ! But at the same time, if one want to print out a Greek letter,
    ! it is always more convenient to directly copy the letters directly
    ! to the code and let it print...
    !!    Letter   Description  Escape-Sequence
    !    -------------------------------------
    !    A        Alpha        \u0391
    !    B        Beta         \u0392
    !    Γ        Gamma        \u0393
    !    Δ        Delta        \u0394
    !    Ε        Epsilon      \u0395
    !    Ζ        Zeta         \u0396
    !    Η        Eta          \u0397
    !    Θ        Theta        \u0398
    !    Ι        Iota         \u0399
    !    Κ        Kappa        \u039A
    !    Λ        Lambda       \u039B
    !    Μ        Mu           \u039C
    !    Ν        Nu           \u039D
    !    Ξ        Xi           \u039E
    !    Ο        Omicron      \u039F
    !    Π        Pi           \u03A0
    !    Ρ        Rho          \u03A1
    !    Σ        Sigma        \u03A3
    !    Τ        Tau          \u03A4
    !    Υ        Upsilon      \u03A5
    !    Φ        Phi          \u03A6
    !    Χ        Chi          \u03A7
    !    Ψ        Psi          \u03A8
    !    Ω        Omega        \u03A9
    !    α        Alpha        \u03B1
    !    β        Beta         \u03B2
    !    γ        Gamma        \u03B3
    !    δ        Delta        \u03B4
    !    ε        Epsilon      \u03B5
    !    ζ        Zeta         \u03B6
    !    η        Eta          \u03B7
    !    θ        Theta        \u03B8
    !    ι        Iota         \u03B9
    !    κ        Kappa        \u03BA
    !    λ        Lambda       \u03BB
    !    μ        Mu           \u03BC
    !    ν        Nu           \u03BD
    !    ξ        Xi           \u03BE
    !    ο        Omicron      \u03BF
    !    π        Pi           \u03C0
    !    ρ        Rho          \u03C1
    !    σ        Sigma        \u03C3
    !    τ        Tau          \u03C4
    !    υ        Upsilon      \u03C5
    !    φ        Phi          \u03C6
    !    χ        Chi          \u03C7
    !    ψ        Psi          \u03C8
    !    ω        Omega        \u03C9

contains

    function nan() result(x)

        real :: x

        x = sqrt(-1.0)

    end function nan

end module libflit_constants
