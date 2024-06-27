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


module libflit_random

    use libflit_array
    use libflit_utility
    use libflit_error
    use libflit_array_operation
    use libflit_sort

    implicit none

    ! C++ random number generators using the PCG random engine
    interface

        ! random number following uniform distribution
        subroutine uniform_rand_int(w, nw, lb, ub) bind(c, name='uniform_rand_int')
            use iso_c_binding, only: c_int
            integer(kind=c_int), dimension(*) :: w
            integer(kind=c_int), value :: nw
            integer(kind=c_int), value :: lb, ub
        end subroutine uniform_rand_int

        ! random number following uniform distribution, given seed
        subroutine uniform_rand_seed_int(w, nw, lb, ub, seed) bind(c, name='uniform_rand_seed_int')
            use iso_c_binding, only: c_int
            integer(kind=c_int), dimension(*) :: w
            integer(kind=c_int), value :: nw
            integer(kind=c_int), value :: lb, ub
            integer(kind=c_int), value :: seed
        end subroutine uniform_rand_seed_int

        ! random number following normal distribution
        subroutine normal_rand_int(w, nw, mu, sigma) bind(c, name='normal_rand_int')
            use iso_c_binding, only: c_float, c_int
            integer(kind=c_int), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: mu, sigma
        end subroutine normal_rand_int

        ! random number following normal distribution, given seed
        subroutine normal_rand_seed_int(w, nw, mu, sigma, seed) bind(c, name='normal_rand_seed_int')
            use iso_c_binding, only: c_float, c_int
            integer(kind=c_int), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: mu, sigma
            integer(kind=c_int), value :: seed
        end subroutine normal_rand_seed_int

        ! random number following uniform distribution
        subroutine uniform_rand(w, nw, lb, ub) bind(c, name='uniform_rand')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: lb, ub
        end subroutine uniform_rand

        ! random number following uniform distribution, given seed
        subroutine uniform_rand_seed(w, nw, lb, ub, seed) bind(c, name='uniform_rand_seed')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: lb, ub
            integer(kind=c_int), value :: seed
        end subroutine uniform_rand_seed

        ! random number following normal distribution
        subroutine normal_rand(w, nw, mu, sigma) bind(c, name='normal_rand')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: mu, sigma
        end subroutine normal_rand

        ! random number following normal distribution, given seed
        subroutine normal_rand_seed(w, nw, mu, sigma, seed) bind(c, name='normal_rand_seed')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: mu, sigma
            integer(kind=c_int), value :: seed
        end subroutine normal_rand_seed

        ! random number following Cauchy distribution
        subroutine cauchy_rand(w, nw, a, b) bind(c, name='cauchy_rand')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: a, b
        end subroutine cauchy_rand

        ! random number following Cauchy distribution, given seed
        subroutine cauchy_rand_seed(w, nw, a, b, seed) bind(c, name='cauchy_rand_seed')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: a, b
            integer(kind=c_int), value :: seed
        end subroutine cauchy_rand_seed

        ! random number following Poisson distribution
        subroutine poisson_rand(w, nw, mean) bind(c, name='poisson_rand')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: mean
        end subroutine poisson_rand

        ! random number following Poisson distribution, given seed
        subroutine poisson_rand_seed(w, nw, mean, seed) bind(c, name='poisson_rand_seed')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: mean
            integer(kind=c_int), value :: seed
        end subroutine poisson_rand_seed

        ! random number following exponential distribution
        subroutine exponential_rand(w, nw, lambda) bind(c, name='exponential_rand')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: lambda
        end subroutine exponential_rand

        ! random number following exponential distribution, given seed
        subroutine exponential_rand_seed(w, nw, lambda, seed) bind(c, name='exponential_rand_seed')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: lambda
            integer(kind=c_int), value :: seed
        end subroutine exponential_rand_seed

        ! random number following uniform distribution
        subroutine uniform_rand_double(w, nw, lb, ub) bind(c, name='uniform_rand_double')
            use iso_c_binding, only: c_double, c_int
            real(kind=c_double), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_double), value :: lb, ub
        end subroutine uniform_rand_double

        ! random number following uniform distribution, given seed
        subroutine uniform_rand_seed_double(w, nw, lb, ub, seed) bind(c, name='uniform_rand_seed_double')
            use iso_c_binding, only: c_double, c_int
            real(kind=c_double), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_double), value :: lb, ub
            integer(kind=c_int), value :: seed
        end subroutine uniform_rand_seed_double

        ! random number following normal distribution
        subroutine normal_rand_double(w, nw, mu, sigma) bind(c, name='normal_rand_double')
            use iso_c_binding, only: c_double, c_int
            real(kind=c_double), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_double), value :: mu, sigma
        end subroutine normal_rand_double

        ! random number following normal distribution, given seed
        subroutine normal_rand_seed_double(w, nw, mu, sigma, seed) bind(c, name='normal_rand_seed_double')
            use iso_c_binding, only: c_double, c_int
            real(kind=c_double), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_double), value :: mu, sigma
            integer(kind=c_int), value :: seed
        end subroutine normal_rand_seed_double

        ! random number following Cauchy distribution
        subroutine cauchy_rand_double(w, nw, a, b) bind(c, name='cauchy_rand_double')
            use iso_c_binding, only: c_double, c_int
            real(kind=c_double), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_double), value :: a, b
        end subroutine cauchy_rand_double

        ! random number following Cauchy distribution, given seed
        subroutine cauchy_rand_seed_double(w, nw, a, b, seed) bind(c, name='cauchy_rand_seed_double')
            use iso_c_binding, only: c_double, c_int
            real(kind=c_double), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_double), value :: a, b
            integer(kind=c_int), value :: seed
        end subroutine cauchy_rand_seed_double

        ! random number following Poisson distribution
        subroutine poisson_rand_double(w, nw, mean) bind(c, name='poisson_rand_double')
            use iso_c_binding, only: c_double, c_int
            real(kind=c_double), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_double), value :: mean
        end subroutine poisson_rand_double

        ! random number following Poisson distribution, given seed
        subroutine poisson_rand_seed_double(w, nw, mean, seed) bind(c, name='poisson_rand_seed_double')
            use iso_c_binding, only: c_double, c_int
            real(kind=c_double), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_double), value :: mean
            integer(kind=c_int), value :: seed
        end subroutine poisson_rand_seed_double

        ! random number following exponential distribution
        subroutine exponential_rand_double(w, nw, lambda) bind(c, name='exponential_rand_double')
            use iso_c_binding, only: c_double, c_int
            real(kind=c_double), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_double), value :: lambda
        end subroutine exponential_rand_double

        ! random number following exponential distribution, given seed
        subroutine exponential_rand_seed_double(w, nw, lambda, seed) bind(c, name='exponential_rand_seed_double')
            use iso_c_binding, only: c_double, c_int
            real(kind=c_double), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_double), value :: lambda
            integer(kind=c_int), value :: seed
        end subroutine exponential_rand_seed_double

    end interface

    interface irandom
        module procedure :: irand
        module procedure :: rand_array_1d_int
        module procedure :: rand_array_2d_int
        module procedure :: rand_array_3d_int
    end interface irandom

    interface random
        module procedure :: rand
        module procedure :: rand_array_1d_float
        module procedure :: rand_array_2d_float
        module procedure :: rand_array_3d_float
    end interface random

    interface drandom
        module procedure :: drand
        module procedure :: rand_array_1d_double
        module procedure :: rand_array_2d_double
        module procedure :: rand_array_3d_double
    end interface drandom

    interface randomize
        module procedure :: randomize_1d_float
        module procedure :: randomize_2d_float
        module procedure :: randomize_3d_float
    end interface randomize

    interface random_permute
        module procedure :: random_permute_int
        module procedure :: random_permute_float
        module procedure :: random_permute_double
        module procedure :: random_permute_complex
        module procedure :: random_permute_dcomplex
        module procedure :: random_permute_logical
    end interface random_permute

    interface random_mask
        module procedure :: random_mask_1d_float
        module procedure :: random_mask_2d_float
        module procedure :: random_mask_3d_float
    end interface random_mask

    private
    public :: irand, rand, drand
    public :: irandom, random, drandom
    public :: randomize
    public :: random_order
    public :: random_string
    public :: random_permute
    public :: random_mask

contains

    !
    !> Return a random integer
    !
    function irand(dist, seed, range, mu, sigma, a, b, mean, lambda) result(r)

        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        integer, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        real, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        real, optional :: a, b
        ! Possion distribution -- mean
        real, optional :: mean
        ! Exponential distribution -- lambda
        real, optional :: lambda
        ! Output
        integer :: r

        integer, dimension(1:1) :: random_value
        integer, dimension(1:2) :: value_range
        real :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        real :: exponential_lambda
        character(len=32) :: distribution

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0, 9]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    call uniform_rand_seed_int(random_value, 1, value_range(1), value_range(2), seed)
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed_int(random_value, 1, normal_mu, normal_sigma, seed)
                    !                case ('cauchy', 'Cauchy')
                    !                    call cauchy_rand_seed(random_value, 1, cauchy_a, cauchy_b, seed)
                    !                case ('poisson', 'Poisson')
                    !                    call poisson_rand_seed(random_value, 1, poisson_mean, seed)
                    !                case ('exp', 'exponential')
                    !                    call exponential_rand_seed(random_value, 1, exponential_lambda, seed)
            end select
        else
            select case (distribution)
                case ('uniform')
                    call uniform_rand_int(random_value, 1, value_range(1), value_range(2))
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_int(random_value, 1, normal_mu, normal_sigma)
                    !                case ('cauchy', 'Cauchy')
                    !                    call cauchy_rand(random_value, 1, cauchy_a, cauchy_b)
                    !                case ('poisson', 'Poisson')
                    !                    call poisson_rand(random_value, 1, poisson_mean)
                    !                case ('exp', 'exponential')
                    !                    call exponential_rand(random_value, 1, exponential_lambda)
            end select
        end if

        r = random_value(1)

    end function irand

    !
    !> Return a 1D array of random integers
    !
    function rand_array_1d_int(n1, dist, seed, range, mu, sigma, a, b, mean, lambda) result(r)

        integer :: n1
        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        integer, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        real, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        real, optional :: a, b
        ! Possion distribution -- mean
        real, optional :: mean
        ! Exponential distribution -- lambda
        real, optional :: lambda
        ! Minimum distance
        !        integer, optional :: spacing
        ! Output
        integer, allocatable, dimension(:) :: r

        integer, allocatable, dimension(:) :: random_value
        integer, dimension(1:2) :: value_range
        real :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        real :: exponential_lambda
        character(len=32) :: distribution
        integer :: n
        !        integer :: random_spacing
        !        integer :: empty_space
        !        integer, allocatable, dimension(:) :: index

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0, 9]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        !                if (present(spacing)) then
        !            random_spacing = spacing
        !        else
        !            random_spacing = 0
        !        end if

        n = n1
        !        if (present(spacing)) then
        !            empty_space = value_range(2) - value_range(1) - (n - 1)*random_spacing
        !            call assert(empty_space >= 0, ' <rand_array_1d> Error: The expected minimum spacing is too big. ')
        !        end if

        allocate (random_value(1:n))

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    call uniform_rand_seed_int(random_value, n, value_range(1), value_range(2), seed)
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed_int(random_value, n, normal_mu, normal_sigma, seed)
                    !                case ('cauchy', 'Cauchy')
                    !                    call cauchy_rand_seed(random_value, n, cauchy_a, cauchy_b, seed)
                    !                case ('poisson', 'Poisson')
                    !                    call poisson_rand_seed(random_value, n, poisson_mean, seed)
                    !                case ('exp', 'exponential')
                    !                    call exponential_rand_seed(random_value, n, exponential_lambda, seed)
            end select
        else
            select case (distribution)
                case ('uniform')
                    call uniform_rand_int(random_value, n, value_range(1), value_range(2))
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_int(random_value, n, normal_mu, normal_sigma)
                    !                case ('cauchy', 'Cauchy')
                    !                    call cauchy_rand(random_value, n, cauchy_a, cauchy_b)
                    !                case ('poisson', 'Poisson')
                    !                    call poisson_rand(random_value, n, poisson_mean)
                    !                case ('exp', 'exponential')
                    !                    call exponential_rand(random_value, n, exponential_lambda)
            end select
        end if

        r = random_value

    end function rand_array_1d_int

    !
    !> Return a 2D array of random floats
    !
    function rand_array_2d_int(n1, n2, dist, seed, range, mu, sigma, a, b, mean, lambda) result(r)

        integer :: n1, n2
        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        real, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        real, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        real, optional :: a, b
        ! Possion distribution -- mean
        real, optional :: mean
        ! Exponential distribution -- lambda
        real, optional :: lambda
        ! Output
        integer, allocatable, dimension(:, :) :: r

        integer, allocatable, dimension(:) :: random_value
        integer, dimension(1:2) :: value_range
        real :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        real :: exponential_lambda
        character(len=32) :: distribution
        integer :: n
        !        integer :: random_spacing
        !        integer :: empty_space
        !        integer, allocatable, dimension(:) :: index

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0.0, 1.0]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        n = n1*n2
        allocate (random_value(1:n))

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    call uniform_rand_seed_int(random_value, n, value_range(1), value_range(2), seed)
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed_int(random_value, n, normal_mu, normal_sigma, seed)
                    !                case ('cauchy', 'Cauchy')
                    !                    call cauchy_rand_seed(random_value, n, cauchy_a, cauchy_b, seed)
                    !                case ('poisson', 'Poisson')
                    !                    call poisson_rand_seed(random_value, n, poisson_mean, seed)
                    !                case ('exp', 'exponential')
                    !                    call exponential_rand_seed(random_value, n, exponential_lambda, seed)
            end select
        else
            select case (distribution)
                case ('uniform')
                    call uniform_rand_int(random_value, n, value_range(1), value_range(2))
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_int(random_value, n, normal_mu, normal_sigma)
                    !                case ('cauchy', 'Cauchy')
                    !                    call cauchy_rand(random_value, n, cauchy_a, cauchy_b)
                    !                case ('poisson', 'Poisson')
                    !                    call poisson_rand(random_value, n, poisson_mean)
                    !                case ('exp', 'exponential')
                    !                    call exponential_rand(random_value, n, exponential_lambda)
            end select
        end if

        r = reshape(random_value, [n1, n2])

    end function rand_array_2d_int

    !
    !> Return a 3D array of random floats
    !
    function rand_array_3d_int(n1, n2, n3, dist, seed, range, mu, sigma, a, b, mean, lambda) result(r)

        integer :: n1, n2, n3
        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        real, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        real, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        real, optional :: a, b
        ! Possion distribution -- mean
        real, optional :: mean
        ! Exponential distribution -- lambda
        real, optional :: lambda
        ! Output
        integer, allocatable, dimension(:, :, :) :: r

        integer, allocatable, dimension(:) :: random_value
        integer, dimension(1:2) :: value_range
        real :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        real :: exponential_lambda
        character(len=32) :: distribution
        integer :: n
        !        integer :: random_spacing
        !        integer :: empty_space
        !        integer, allocatable, dimension(:) :: index

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0.0, 1.0]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        n = n1*n2*n3
        allocate (random_value(1:n))

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    call uniform_rand_seed_int(random_value, n, value_range(1), value_range(2), seed)
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed_int(random_value, n, normal_mu, normal_sigma, seed)
                    !                case ('cauchy', 'Cauchy')
                    !                    call cauchy_rand_seed(random_value, n, cauchy_a, cauchy_b, seed)
                    !                case ('poisson', 'Poisson')
                    !                    call poisson_rand_seed(random_value, n, poisson_mean, seed)
                    !                case ('exp', 'exponential')
                    !                    call exponential_rand_seed(random_value, n, exponential_lambda, seed)
            end select
        else
            select case (distribution)
                case ('uniform')
                    call uniform_rand_int(random_value, n, value_range(1), value_range(2))
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_int(random_value, n, normal_mu, normal_sigma)
                    !                case ('cauchy', 'Cauchy')
                    !                    call cauchy_rand(random_value, n, cauchy_a, cauchy_b)
                    !                case ('poisson', 'Poisson')
                    !                    call poisson_rand(random_value, n, poisson_mean)
                    !                case ('exp', 'exponential')
                    !                    call exponential_rand(random_value, n, exponential_lambda)
            end select
        end if

        !        allocate (r(1:n1, 1:n2, 1:n3), source=reshape(random_value, [n1, n2, n3]))
        r = reshape(random_value, [n1, n2, n3])

    end function rand_array_3d_int


    !
    !> Return a random float
    !
    function rand(dist, seed, range, mu, sigma, a, b, mean, lambda) result(r)

        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        real, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        real, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        real, optional :: a, b
        ! Possion distribution -- mean
        real, optional :: mean
        ! Exponential distribution -- lambda
        real, optional :: lambda
        ! Output
        real :: r

        real, dimension(1:1) :: random_value
        real, dimension(1:2) :: value_range
        real :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        real :: exponential_lambda
        character(len=32) :: distribution

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0.0, 1.0]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    call uniform_rand_seed(random_value, 1, value_range(1), value_range(2), seed)
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed(random_value, 1, normal_mu, normal_sigma, seed)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_seed(random_value, 1, cauchy_a, cauchy_b, seed)
                case ('poisson', 'Poisson')
                    call poisson_rand_seed(random_value, 1, poisson_mean, seed)
                case ('exp', 'exponential')
                    call exponential_rand_seed(random_value, 1, exponential_lambda, seed)
            end select
        else
            ! When seed = -1, the random generator uses null seed
            select case (distribution)
                case ('uniform')
                    call uniform_rand(random_value, 1, value_range(1), value_range(2))
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand(random_value, 1, normal_mu, normal_sigma)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand(random_value, 1, cauchy_a, cauchy_b)
                case ('poisson', 'Poisson')
                    call poisson_rand(random_value, 1, poisson_mean)
                case ('exp', 'exponential')
                    call exponential_rand(random_value, 1, exponential_lambda)
            end select
        end if

        r = random_value(1)

    end function rand

    !
    !> Return a 1D array of random floats
    !
    function rand_array_1d_float(n1, dist, seed, range, mu, sigma, a, b, mean, lambda, spacing) result(r)

        integer :: n1
        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        real, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        real, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        real, optional :: a, b
        ! Possion distribution -- mean
        real, optional :: mean
        ! Exponential distribution -- lambda
        real, optional :: lambda
        ! Minimum distance
        real, optional :: spacing
        ! Output
        real, allocatable, dimension(:) :: r

        real, allocatable, dimension(:) :: random_value
        real, dimension(1:2) :: value_range
        real :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        real :: exponential_lambda
        character(len=32) :: distribution
        integer :: n
        real :: random_spacing
        real :: empty_space
        integer, allocatable, dimension(:) :: index

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0.0, 1.0]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        if (present(spacing)) then
            random_spacing = spacing
        else
            random_spacing = 0.0
        end if

        n = n1
        if (present(spacing)) then
            empty_space = value_range(2) - value_range(1) - (n - 1.0)*random_spacing
            call assert(empty_space >= 0.0, ' <rand_array_1d> Error: The expected minimum spacing is too big. ')
        end if

        allocate (random_value(1:n))

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    if (present(spacing)) then
                        call uniform_rand_seed(random_value, n, 0.0, 1.0, seed)
                        call sort_index(random_value, index)
                        random_value = empty_space*random_value + value_range(1) + random_spacing*linspace(0.0, n - 1.0, n)
                        random_value = random_value(index)
                    else
                        call uniform_rand_seed(random_value, n, value_range(1), value_range(2), seed)
                    end if
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed(random_value, n, normal_mu, normal_sigma, seed)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_seed(random_value, n, cauchy_a, cauchy_b, seed)
                case ('poisson', 'Poisson')
                    call poisson_rand_seed(random_value, n, poisson_mean, seed)
                case ('exp', 'exponential')
                    call exponential_rand_seed(random_value, n, exponential_lambda, seed)
            end select
        else
            select case (distribution)
                case ('uniform')
                    if (present(spacing)) then
                        call uniform_rand(random_value, n, 0.0, 1.0)
                        call sort_index(random_value, index)
                        random_value = empty_space*random_value + value_range(1) + random_spacing*linspace(0.0, n - 1.0, n)
                        random_value = random_value(index)
                    else
                        call uniform_rand(random_value, n, value_range(1), value_range(2))
                    end if
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand(random_value, n, normal_mu, normal_sigma)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand(random_value, n, cauchy_a, cauchy_b)
                case ('poisson', 'Poisson')
                    call poisson_rand(random_value, n, poisson_mean)
                case ('exp', 'exponential')
                    call exponential_rand(random_value, n, exponential_lambda)
            end select
        end if

        r = random_value

    end function rand_array_1d_float

    !
    !> Return a 2D array of random floats
    !
    function rand_array_2d_float(n1, n2, dist, seed, range, mu, sigma, a, b, mean, lambda, spacing) result(r)

        integer :: n1, n2
        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        real, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        real, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        real, optional :: a, b
        ! Possion distribution -- mean
        real, optional :: mean
        ! Exponential distribution -- lambda
        real, optional :: lambda
        ! Minimum distance
        real, optional :: spacing
        ! Output
        real, allocatable, dimension(:, :) :: r

        real, allocatable, dimension(:) :: random_value
        real, dimension(1:2) :: value_range
        real :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        real :: exponential_lambda
        character(len=32) :: distribution
        integer :: n
        real :: random_spacing
        real :: empty_space
        integer, allocatable, dimension(:) :: index

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0.0, 1.0]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        if (present(spacing)) then
            random_spacing = spacing
        else
            random_spacing = 0.0
        end if

        n = n1*n2
        if (present(spacing)) then
            empty_space = value_range(2) - value_range(1) - (n - 1.0)*random_spacing
            call assert(empty_space >= 0.0, ' <rand_array_1d> Error: The expected minimum spacing is too big. ')
        end if

        allocate (random_value(1:n))

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    if (present(spacing)) then
                        call uniform_rand_seed(random_value, n, 0.0, 1.0, seed)
                        call sort_index(random_value, index)
                        random_value = empty_space*random_value + value_range(1) + random_spacing*linspace(0.0, n - 1.0, n)
                        random_value = random_value(index)
                    else
                        call uniform_rand_seed(random_value, n, value_range(1), value_range(2), seed)
                    end if
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed(random_value, n, normal_mu, normal_sigma, seed)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_seed(random_value, n, cauchy_a, cauchy_b, seed)
                case ('poisson', 'Poisson')
                    call poisson_rand_seed(random_value, n, poisson_mean, seed)
                case ('exp', 'exponential')
                    call exponential_rand_seed(random_value, n, exponential_lambda, seed)
            end select
        else
            select case (distribution)
                case ('uniform')
                    if (present(spacing)) then
                        call uniform_rand(random_value, n, 0.0, 1.0)
                        call sort_index(random_value, index)
                        random_value = empty_space*random_value + value_range(1) + random_spacing*linspace(0.0, n - 1.0, n)
                        random_value = random_value(index)
                    else
                        call uniform_rand(random_value, n, value_range(1), value_range(2))
                    end if
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand(random_value, n, normal_mu, normal_sigma)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand(random_value, n, cauchy_a, cauchy_b)
                case ('poisson', 'Poisson')
                    call poisson_rand(random_value, n, poisson_mean)
                case ('exp', 'exponential')
                    call exponential_rand(random_value, n, exponential_lambda)
            end select
        end if

        r = reshape(random_value, [n1, n2])

    end function rand_array_2d_float

    !
    !> Return a 3D array of random floats
    !
    function rand_array_3d_float(n1, n2, n3, dist, seed, range, mu, sigma, a, b, mean, lambda, spacing) result(r)

        integer :: n1, n2, n3
        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        real, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        real, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        real, optional :: a, b
        ! Possion distribution -- mean
        real, optional :: mean
        ! Exponential distribution -- lambda
        real, optional :: lambda
        ! Output
        real, allocatable, dimension(:, :, :) :: r
        ! Minimum distance
        real, optional :: spacing

        real, allocatable, dimension(:) :: random_value
        real, dimension(1:2) :: value_range
        real :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        real :: exponential_lambda
        character(len=32) :: distribution
        integer :: n
        real :: random_spacing
        real :: empty_space
        integer, allocatable, dimension(:) :: index

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0.0, 1.0]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        if (present(spacing)) then
            random_spacing = spacing
        else
            random_spacing = 0.0
        end if

        n = n1*n2*n3
        if (present(spacing)) then
            empty_space = value_range(2) - value_range(1) - (n - 1.0)*random_spacing
            call assert(empty_space >= 0.0, ' <rand_array_1d> Error: The expected minimum spacing is too big. ')
        end if

        allocate (random_value(1:n))

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    if (present(spacing)) then
                        call uniform_rand_seed(random_value, n, 0.0, 1.0, seed)
                        call sort_index(random_value, index)
                        random_value = empty_space*random_value + value_range(1) + random_spacing*linspace(0.0, n - 1.0, n)
                        random_value = random_value(index)
                    else
                        call uniform_rand_seed(random_value, n, value_range(1), value_range(2), seed)
                    end if
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed(random_value, n, normal_mu, normal_sigma, seed)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_seed(random_value, n, cauchy_a, cauchy_b, seed)
                case ('poisson', 'Poisson')
                    call poisson_rand_seed(random_value, n, poisson_mean, seed)
                case ('exp', 'exponential')
                    call exponential_rand_seed(random_value, n, exponential_lambda, seed)
            end select
        else
            select case (distribution)
                case ('uniform')
                    if (present(spacing)) then
                        call uniform_rand(random_value, n, 0.0, 1.0)
                        call sort_index(random_value, index)
                        random_value = empty_space*random_value + value_range(1) + random_spacing*linspace(0.0, n - 1.0, n)
                        random_value = random_value(index)
                    else
                        call uniform_rand(random_value, n, value_range(1), value_range(2))
                    end if
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand(random_value, n, normal_mu, normal_sigma)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand(random_value, n, cauchy_a, cauchy_b)
                case ('poisson', 'Poisson')
                    call poisson_rand(random_value, n, poisson_mean)
                case ('exp', 'exponential')
                    call exponential_rand(random_value, n, exponential_lambda)
            end select
        end if

        r = reshape(random_value, [n1, n2, n3])

    end function rand_array_3d_float

    !
    !> Return a random integer
    !
    function drand(dist, seed, range, mu, sigma, a, b, mean, lambda) result(r)

        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        integer, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        double precision, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        double precision, optional :: a, b
        ! Possion distribution -- mean
        double precision, optional :: mean
        ! Exponential distribution -- lambda
        double precision, optional :: lambda
        ! Output
        double precision :: r

        double precision, dimension(1:1) :: random_value
        double precision, dimension(1:2) :: value_range
        double precision :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        double precision :: exponential_lambda
        character(len=32) :: distribution

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0, 9]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    call uniform_rand_seed_double(random_value, 1, value_range(1), value_range(2), seed)
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed_double(random_value, 1, normal_mu, normal_sigma, seed)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_seed_double(random_value, 1, cauchy_a, cauchy_b, seed)
                case ('poisson', 'Poisson')
                    call poisson_rand_seed_double(random_value, 1, poisson_mean, seed)
                case ('exp', 'exponential')
                    call exponential_rand_seed_double(random_value, 1, exponential_lambda, seed)
            end select
        else
            select case (distribution)
                case ('uniform')
                    call uniform_rand_double(random_value, 1, value_range(1), value_range(2))
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_double(random_value, 1, normal_mu, normal_sigma)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_double(random_value, 1, cauchy_a, cauchy_b)
                case ('poisson', 'Poisson')
                    call poisson_rand_double(random_value, 1, poisson_mean)
                case ('exp', 'exponential')
                    call exponential_rand_double(random_value, 1, exponential_lambda)
            end select
        end if

        r = random_value(1)

    end function drand

    !
    !> Return a 1D array of random integers
    !
    function rand_array_1d_double(n1, dist, seed, range, mu, sigma, a, b, mean, lambda, spacing) result(r)

        integer :: n1
        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        integer, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        double precision, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        double precision, optional :: a, b
        ! Possion distribution -- mean
        double precision, optional :: mean
        ! Exponential distribution -- lambda
        double precision, optional :: lambda
        ! Minimum distance
        double precision, optional :: spacing
        ! Output
        double precision, allocatable, dimension(:) :: r

        double precision, allocatable, dimension(:) :: random_value
        double precision, dimension(1:2) :: value_range
        double precision :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        double precision :: exponential_lambda
        character(len=32) :: distribution
        integer :: n
        real :: random_spacing
        real :: empty_space
        integer, allocatable, dimension(:) :: index

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0, 9]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        if (present(spacing)) then
            random_spacing = spacing
        else
            random_spacing = 0.0
        end if

        n = n1
        if (present(spacing)) then
            empty_space = value_range(2) - value_range(1) - (n - 1.0)*random_spacing
            call assert(empty_space >= 0.0, ' <rand_array_1d> Error: The expected minimum spacing is too big. ')
        end if

        allocate (random_value(1:n))

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    if (present(spacing)) then
                        call uniform_rand_seed_double(random_value, n, 0.0d0, 1.0d0, seed)
                        call sort_index(random_value, index)
                        random_value = empty_space*random_value + value_range(1) + random_spacing*linspace(0.0, n - 1.0, n)
                        random_value = random_value(index)
                    else
                        call uniform_rand_seed_double(random_value, n, value_range(1), value_range(2), seed)
                    end if
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed_double(random_value, n, normal_mu, normal_sigma, seed)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_seed_double(random_value, n, cauchy_a, cauchy_b, seed)
                case ('poisson', 'Poisson')
                    call poisson_rand_seed_double(random_value, n, poisson_mean, seed)
                case ('exp', 'exponential')
                    call exponential_rand_seed_double(random_value, n, exponential_lambda, seed)
            end select
        else
            select case (distribution)
                case ('uniform')
                    if (present(spacing)) then
                        call uniform_rand_double(random_value, n, 0.0d0, 1.0d0)
                        call sort_index(random_value, index)
                        random_value = empty_space*random_value + value_range(1) + random_spacing*linspace(0.0, n - 1.0, n)
                        random_value = random_value(index)
                    else
                        call uniform_rand_double(random_value, n, value_range(1), value_range(2))
                    end if
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_double(random_value, n, normal_mu, normal_sigma)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_double(random_value, n, cauchy_a, cauchy_b)
                case ('poisson', 'Poisson')
                    call poisson_rand_double(random_value, n, poisson_mean)
                case ('exp', 'exponential')
                    call exponential_rand_double(random_value, n, exponential_lambda)
            end select
        end if

        r = random_value

    end function rand_array_1d_double

    !
    !> Return a 2D array of random floats
    !
    function rand_array_2d_double(n1, n2, dist, seed, range, mu, sigma, a, b, mean, lambda, spacing) result(r)

        integer :: n1, n2
        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        double precision, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        double precision, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        double precision, optional :: a, b
        ! Possion distribution -- mean
        double precision, optional :: mean
        ! Exponential distribution -- lambda
        double precision, optional :: lambda
        ! Minimum distance
        double precision, optional :: spacing
        ! Output
        double precision, allocatable, dimension(:, :) :: r

        double precision, allocatable, dimension(:) :: random_value
        double precision, dimension(1:2) :: value_range
        double precision :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        double precision :: exponential_lambda
        character(len=32) :: distribution
        integer :: n
        real :: random_spacing
        real :: empty_space
        integer, allocatable, dimension(:) :: index

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0.0, 1.0]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        if (present(spacing)) then
            random_spacing = spacing
        else
            random_spacing = 0.0
        end if

        n = n1*n2
        if (present(spacing)) then
            empty_space = value_range(2) - value_range(1) - (n - 1.0)*random_spacing
            call assert(empty_space >= 0.0, ' <rand_array_1d> Error: The expected minimum spacing is too big. ')
        end if

        allocate (random_value(1:n))

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    if (present(spacing)) then
                        call uniform_rand_seed_double(random_value, n, 0.0d0, 1.0d0, seed)
                        call sort_index(random_value, index)
                        random_value = empty_space*random_value + value_range(1) + random_spacing*linspace(0.0, n - 1.0, n)
                        random_value = random_value(index)
                    else
                        call uniform_rand_seed_double(random_value, n, value_range(1), value_range(2), seed)
                    end if
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed_double(random_value, n, normal_mu, normal_sigma, seed)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_seed_double(random_value, n, cauchy_a, cauchy_b, seed)
                case ('poisson', 'Poisson')
                    call poisson_rand_seed_double(random_value, n, poisson_mean, seed)
                case ('exp', 'exponential')
                    call exponential_rand_seed_double(random_value, n, exponential_lambda, seed)
            end select
        else
            select case (distribution)
                case ('uniform')
                    if (present(spacing)) then
                        call uniform_rand_double(random_value, n, 0.0d0, 1.0d0)
                        call sort_index(random_value, index)
                        random_value = empty_space*random_value + value_range(1) + random_spacing*linspace(0.0, n - 1.0, n)
                        random_value = random_value(index)
                    else
                        call uniform_rand_double(random_value, n, value_range(1), value_range(2))
                    end if
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_double(random_value, n, normal_mu, normal_sigma)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_double(random_value, n, cauchy_a, cauchy_b)
                case ('poisson', 'Poisson')
                    call poisson_rand_double(random_value, n, poisson_mean)
                case ('exp', 'exponential')
                    call exponential_rand_double(random_value, n, exponential_lambda)
            end select
        end if

        r = reshape(random_value, [n1, n2])

    end function rand_array_2d_double

    !
    !> Return a 3D array of random floats
    !
    function rand_array_3d_double(n1, n2, n3, dist, seed, range, mu, sigma, a, b, mean, lambda, spacing) result(r)

        integer :: n1, n2, n3
        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        double precision, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        double precision, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        double precision, optional :: a, b
        ! Possion distribution -- mean
        double precision, optional :: mean
        ! Exponential distribution -- lambda
        double precision, optional :: lambda
        ! Minimum distance
        double precision, optional :: spacing
        ! Output
        double precision, allocatable, dimension(:, :, :) :: r

        double precision, allocatable, dimension(:) :: random_value
        double precision, dimension(1:2) :: value_range
        double precision :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        double precision :: exponential_lambda
        character(len=32) :: distribution
        integer :: n
        real :: random_spacing
        real :: empty_space
        integer, allocatable, dimension(:) :: index

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0.0, 1.0]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        if (present(spacing)) then
            random_spacing = spacing
        else
            random_spacing = 0.0
        end if

        n = n1*n2*n3
        if (present(spacing)) then
            empty_space = value_range(2) - value_range(1) - (n - 1.0)*random_spacing
            call assert(empty_space >= 0.0, ' <rand_array_1d> Error: The expected minimum spacing is too big. ')
        end if

        allocate (random_value(1:n))

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    if (present(spacing)) then
                        call uniform_rand_seed_double(random_value, n, 0.0d0, 1.0d0, seed)
                        call sort_index(random_value, index)
                        random_value = empty_space*random_value + value_range(1) + random_spacing*linspace(0.0, n - 1.0, n)
                        random_value = random_value(index)
                    else
                        call uniform_rand_seed_double(random_value, n, value_range(1), value_range(2), seed)
                    end if
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed_double(random_value, n, normal_mu, normal_sigma, seed)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_seed_double(random_value, n, cauchy_a, cauchy_b, seed)
                case ('poisson', 'Poisson')
                    call poisson_rand_seed_double(random_value, n, poisson_mean, seed)
                case ('exp', 'exponential')
                    call exponential_rand_seed_double(random_value, n, exponential_lambda, seed)
            end select
        else
            select case (distribution)
                case ('uniform')
                    if (present(spacing)) then
                        call uniform_rand_double(random_value, n, 0.0d0, 1.0d0)
                        call sort_index(random_value, index)
                        random_value = empty_space*random_value + value_range(1) + random_spacing*linspace(0.0, n - 1.0, n)
                        random_value = random_value(index)
                    else
                        call uniform_rand_double(random_value, n, value_range(1), value_range(2))
                    end if
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_double(random_value, n, normal_mu, normal_sigma)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_double(random_value, n, cauchy_a, cauchy_b)
                case ('poisson', 'Poisson')
                    call poisson_rand_double(random_value, n, poisson_mean)
                case ('exp', 'exponential')
                    call exponential_rand_double(random_value, n, exponential_lambda)
            end select
        end if

        r = reshape(random_value, [n1, n2, n3])

    end function rand_array_3d_double


    !
    !> Randomize a 1D float array
    !
    subroutine randomize_1d_float(r, dist, seed, range, mu, sigma, a, b, mean, lambda)

        real, dimension(:), intent(inout) :: r

        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        real, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        real, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        real, optional :: a, b
        ! Possion distribution -- mean
        real, optional :: mean
        ! Exponential distribution -- lambda
        real, optional :: lambda
        ! Output
        real, allocatable, dimension(:) :: random_value

        real, dimension(1:2) :: value_range
        real :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        real :: exponential_lambda
        character(len=32) :: distribution
        integer :: n, n1

        n1 = size(r)

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0.0, 1.0]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        n = n1
        allocate (random_value(1:n))

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    call uniform_rand_seed(random_value, n, value_range(1), value_range(2), seed)
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed(random_value, n, normal_mu, normal_sigma, seed)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_seed(random_value, n, cauchy_a, cauchy_b, seed)
                case ('poisson', 'Poisson')
                    call poisson_rand_seed(random_value, n, poisson_mean, seed)
                case ('exp', 'exponential')
                    call exponential_rand_seed(random_value, n, exponential_lambda, seed)
            end select
        else
            select case (distribution)
                case ('uniform')
                    call uniform_rand(random_value, n, value_range(1), value_range(2))
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand(random_value, n, normal_mu, normal_sigma)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand(random_value, n, cauchy_a, cauchy_b)
                case ('poisson', 'Poisson')
                    call poisson_rand(random_value, n, poisson_mean)
                case ('exp', 'exponential')
                    call exponential_rand(random_value, n, exponential_lambda)
            end select
        end if

        r = random_value

    end subroutine randomize_1d_float

    !
    !> Randomize a 2D float array
    !
    subroutine randomize_2d_float(r, dist, seed, range, mu, sigma, a, b, mean, lambda)

        real, dimension(:, :), intent(inout) :: r

        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        real, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        real, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        real, optional :: a, b
        ! Possion distribution -- mean
        real, optional :: mean
        ! Exponential distribution -- lambda
        real, optional :: lambda
        ! Output
        real, allocatable, dimension(:) :: random_value

        real, dimension(1:2) :: value_range
        real :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        real :: exponential_lambda
        character(len=32) :: distribution
        integer :: n, n1, n2

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0.0, 1.0]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        n1 = size(r, 1)
        n2 = size(r, 2)
        n = n1*n2
        allocate (random_value(1:n))

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    call uniform_rand_seed(random_value, n, value_range(1), value_range(2), seed)
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed(random_value, n, normal_mu, normal_sigma, seed)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_seed(random_value, n, cauchy_a, cauchy_b, seed)
                case ('poisson', 'Poisson')
                    call poisson_rand_seed(random_value, n, poisson_mean, seed)
                case ('exp', 'exponential')
                    call exponential_rand_seed(random_value, n, exponential_lambda, seed)
            end select
        else
            select case (distribution)
                case ('uniform')
                    call uniform_rand(random_value, n, value_range(1), value_range(2))
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand(random_value, n, normal_mu, normal_sigma)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand(random_value, n, cauchy_a, cauchy_b)
                case ('poisson', 'Poisson')
                    call poisson_rand(random_value, n, poisson_mean)
                case ('exp', 'exponential')
                    call exponential_rand(random_value, n, exponential_lambda)
            end select
        end if

        r = reshape(random_value, [n1, n2])

    end subroutine randomize_2d_float

    !
    !> Randomize a 3D float array
    !
    subroutine randomize_3d_float(r, dist, seed, range, mu, sigma, a, b, mean, lambda)

        real, dimension(:, :, :), intent(inout) :: r

        ! Given seed or not
        integer, optional :: seed
        ! Distribution
        character(len=*), optional :: dist
        ! Uniform distribution -- value range
        real, dimension(:), optional :: range
        ! Gaussian distribution -- mean and standard deviation
        real, optional :: mu, sigma
        ! Cauchy distribution -- parameters a and b
        real, optional :: a, b
        ! Possion distribution -- mean
        real, optional :: mean
        ! Exponential distribution -- lambda
        real, optional :: lambda
        ! Output
        real, allocatable, dimension(:) :: random_value

        real, dimension(1:2) :: value_range
        real :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
        real :: exponential_lambda
        character(len=32) :: distribution
        integer :: n, n1, n2, n3

        if (present(dist)) then
            distribution = trim(adjustl(dist))
        else
            distribution = 'uniform'
        end if

        if (present(range)) then
            value_range = range
        else
            value_range = [0.0, 1.0]
        end if

        if (present(mu)) then
            normal_mu = mu
        else
            normal_mu = 0.0
        end if

        if (present(sigma)) then
            normal_sigma = sigma
        else
            normal_sigma = 1.0
        end if

        if (present(a)) then
            cauchy_a = a
        else
            cauchy_a = 0.0
        end if

        if (present(b)) then
            cauchy_b = b
        else
            cauchy_b = 1.0
        end if

        if (present(mean)) then
            poisson_mean = mean
        else
            poisson_mean = 0.0
        end if

        if (present(lambda)) then
            exponential_lambda = lambda
        else
            exponential_lambda = 1.0
        end if

        n1 = size(r, 1)
        n2 = size(r, 2)
        n3 = size(r, 3)
        n = n1*n2*n3
        allocate (random_value(1:n))

        if (present(seed) .and. seed > 0) then
            select case (distribution)
                case ('uniform')
                    call uniform_rand_seed(random_value, n, value_range(1), value_range(2), seed)
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand_seed(random_value, n, normal_mu, normal_sigma, seed)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand_seed(random_value, n, cauchy_a, cauchy_b, seed)
                case ('poisson', 'Poisson')
                    call poisson_rand_seed(random_value, n, poisson_mean, seed)
                case ('exp', 'exponential')
                    call exponential_rand_seed(random_value, n, exponential_lambda, seed)
            end select
        else
            select case (distribution)
                case ('uniform')
                    call uniform_rand(random_value, n, value_range(1), value_range(2))
                case ('normal', 'gaussian', 'Gaussian')
                    call normal_rand(random_value, n, normal_mu, normal_sigma)
                case ('cauchy', 'Cauchy')
                    call cauchy_rand(random_value, n, cauchy_a, cauchy_b)
                case ('poisson', 'Poisson')
                    call poisson_rand(random_value, n, poisson_mean)
                case ('exp', 'exponential')
                    call exponential_rand(random_value, n, exponential_lambda)
            end select
        end if

        r = reshape(random_value, [n1, n2, n3])

    end subroutine randomize_3d_float

    !
    !> Generate a random ordering of the integers 1 ... n.
    !
    function random_order(n, seed) result(order)

        integer :: n
        integer, optional :: seed
        integer, allocatable, dimension(:) :: order

        integer :: i, j, k
        real, allocatable, dimension(:) :: wkr
        real :: wk

        order = regspace(1, 1, n)
        if (present(seed)) then
            wkr = random(n, seed=seed)
        else
            wkr = random(n)
        end if

        ! Starting at the end, swap the current last indicator with one
        ! randomly chosen from those preceeding it.
        if (present(seed)) then
            do i = n, 2, -1
                wk = wkr(i) !rand(range=[0.0, 1.0], seed=seed*i)
                j = 1 + i*wk
                if (j < i) then
                    k = order(i)
                    order(i) = order(j)
                    order(j) = k
                end if
            end do
        else
            do i = n, 2, -1
                wk = wkr(i) !rand(range=[0.0, 1.0])
                j = 1 + i*wk
                if (j < i) then
                    k = order(i)
                    order(i) = order(j)
                    order(j) = k
                end if
            end do
        end if

    end function random_order

    !
    !> Generate random string of A-Z, a-z, 0-9
    !
    function random_string(nc) result(str)

        integer, optional :: nc

        integer :: nchar, i
        character(len=1), dimension(1:72) :: pool
        character(len=:), allocatable :: str

        if (present(nc)) then
            nchar = max(1, nc)
        else
            nchar = 10
        end if

        pool = [ &
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', &
            'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', &
            'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', &
            'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', &
            '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', &
            '0', '9', '8', '7', '6', '5', '4', '3', '2', '1' &
            ]

        allocate (character(len=nchar)::str)
        do i = 1, nchar
            str(i:i) = pool(nint(rand(range=[1.0, 72.0])))
        end do

    end function random_string


    function random_mask_1d_float(m, n, seed) result(mm)

        real, dimension(:), intent(in) :: m
        integer, intent(in) :: n
        integer, intent(in), optional :: seed
        real, allocatable, dimension(:) :: mm

        integer, allocatable, dimension(:) :: i
        integer :: l

        mm = m

        call assert(n <= size(m), ' <random_mask_1d_float> Error: n must <= size(m)')

        if (present(seed)) then
            i = random_order(size(m), seed=seed)
        else
            i = random_order(size(m))
        end if

        do l = 1, n
            mm(i(l)) = 0.0
        end do

    end function random_mask_1d_float

    function random_mask_2d_float(m, n, seed) result(mm)

        real, dimension(:, :), intent(in) :: m
        integer, intent(in) :: n
        integer, intent(in), optional :: seed
        real, allocatable, dimension(:, :) :: mm

        integer, allocatable, dimension(:) :: ij
        integer, allocatable, dimension(:, :) :: index
        integer :: i, j, l

        mm = m

        call assert(n <= size(m), ' <random_mask_2d_float> Error: n must <= size(m)')

        if (present(seed)) then
            ij = random_order(size(m), seed=seed)
        else
            ij = random_order(size(m))
        end if

        index = zeros(size(m), 3)
        l = 1
        do i = 1, size(m, 1)
            do j = 1, size(m, 2)
                index(l, :) = [l, i, j]
                l = l + 1
            end do
        end do

        do l = 1, n
            mm(index(ij(l), 2), index(ij(l), 3)) = 0.0
        end do

    end function random_mask_2d_float

    function random_mask_3d_float(m, n, seed) result(mm)

        real, dimension(:, :, :), intent(in) :: m
        integer, intent(in) :: n
        integer, intent(in), optional :: seed
        real, allocatable, dimension(:, :, :) :: mm

        integer, allocatable, dimension(:) :: ijk
        integer, allocatable, dimension(:, :) :: index
        integer :: i, j, k, l

        mm = m

        call assert(n <= size(m), ' <random_mask_3d_float> Error: n must <= size(m)')

        if (present(seed)) then
            ijk = random_order(size(m), seed=seed)
        else
            ijk = random_order(size(m))
        end if

        index = zeros(size(m), 4)
        l = 1
        do i = 1, size(m, 1)
            do j = 1, size(m, 2)
                do k = 1, size(m, 3)
                    index(l, :) = [l, i, j, k]
                    l = l + 1
                end do
            end do
        end do

        do l = 1, n
            mm(index(ijk(l), 2), index(ijk(l), 3), index(ijk(l), 4)) = 0.0
        end do

    end function random_mask_3d_float

    !
    !> Random pertutation of a 1D array of integers
    !
#define T int
#define TT integer
#include "template_random_permute.f90"

#define T float
#define TT real
#include "template_random_permute.f90"

#define T double
#define TT double precision
#include "template_random_permute.f90"

#define T complex
#define TT complex
#include "template_random_permute.f90"

#define T dcomplex
#define TT double complex
#include "template_random_permute.f90"

#define T logical
#define TT logical
#include "template_random_permute.f90"

end module libflit_random
