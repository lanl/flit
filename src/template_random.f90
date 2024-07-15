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

#define PASTE(X)            X
#define PASTE2(X)           PASTE(X)_
#define CONCATHELP(X, Y)    PASTE2(X)Y
#define CONCAT(X, Y)        CONCATHELP(X, Y)

#define uniform_rand_seed_      CONCAT(uniform_rand_seed, T)
#define normal_rand_seed_      CONCAT(normal_rand_seed, T)
#define cauchy_rand_seed_      CONCAT(cauchy_rand_seed, T)
#define poisson_rand_seed_      CONCAT(poisson_rand_seed, T)
#define exponential_rand_seed_      CONCAT(exponential_rand_seed, T)
#define uniform_rand_      CONCAT(uniform_rand, T)
#define normal_rand_      CONCAT(normal_rand, T)
#define cauchy_rand_      CONCAT(cauchy_rand, T)
#define poisson_rand_      CONCAT(poisson_rand, T)
#define exponential_rand_      CONCAT(exponential_rand, T)

#define rand_      CONCAT(rand, T)
#define rand_array_1d_      CONCAT(rand_array_1d, T)
#define rand_array_2d_      CONCAT(rand_array_2d, T)
#define rand_array_3d_      CONCAT(rand_array_3d, T)

!
!> Return a random number
!
function rand_(dist, seed, range, mu, sigma, a, b, mean, lambda) result(r)

    ! Given seed or not
    integer, optional :: seed
    ! Distribution
    character(len=*), optional :: dist
    ! Uniform distribution -- value range
    TT, dimension(:), optional :: range
    ! Gaussian distribution -- mean and standard deviation
    TTT, optional :: mu, sigma
    ! Cauchy distribution -- parameters a and b
    TTT, optional :: a, b
    ! Possion distribution -- mean
    TTT, optional :: mean
    ! Exponential distribution -- lambda
    TTT, optional :: lambda
    ! Output
    TT :: r

    TTT, dimension(1:1) :: random_value
    TTT, dimension(1:2) :: value_range
    TTT :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
    TTT :: exponential_lambda
    character(len=32) :: distribution
    integer :: rs

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

    if (present(seed)) then
        rs = seed
    else
        rs = -1
    end if

    random_value = rand_array_1d_(1, distribution, rs, nTT(value_range), normal_mu, normal_sigma, &
        cauchy_a, cauchy_b, poisson_mean, exponential_lambda)

    r = nTT(random_value(1))

end function rand_

!
!> Return a 1D array of random numbers
!
function rand_array_1d_(n1, dist, seed, range, mu, sigma, a, b, mean, lambda, spacing) result(r)

    integer :: n1
    ! Given seed or not
    integer, optional :: seed
    ! Distribution
    character(len=*), optional :: dist
    ! Uniform distribution -- value range
    TT, dimension(:), optional :: range
    ! Gaussian distribution -- mean and standard deviation
    TTT, optional :: mu, sigma
    ! Cauchy distribution -- parameters a and b
    TTT, optional :: a, b
    ! Possion distribution -- mean
    TTT, optional :: mean
    ! Exponential distribution -- lambda
    TTT, optional :: lambda
    ! Minimum distance
    TT, optional :: spacing
    ! Output
    TT, allocatable, dimension(:) :: r

    TTT, allocatable, dimension(:) :: random_value
    TTT, dimension(1:2) :: value_range
    TTT :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
    TTT :: exponential_lambda
    character(len=32) :: distribution
    integer :: n
    TTT :: random_spacing
    TTT :: empty_space
    integer, allocatable, dimension(:) :: index
    integer :: rs

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

    if (present(seed)) then
        rs = seed
    else
        rs = -1
    end if

    n = n1
    if (random_spacing > 0) then
        empty_space = value_range(2) - value_range(1) - (n - 1.0)*random_spacing
        call assert(empty_space >= 0.0, ' <rand_array_1d> Error: The expected minimum spacing is too big. ')
    end if

    allocate (random_value(1:n))

    if (rs >= 0) then
        select case (distribution)
            case ('uniform')
                if (random_spacing > 0) then
                    call uniform_rand_seed_(random_value, n, 0.0, 1.0, rs)
                    call sort_index(random_value, index)
                    random_value = empty_space*random_value + value_range(1) + random_spacing*linspace(0.0, n - 1.0, n)
                    random_value = random_value(index)
                else
                    call uniform_rand_seed_(random_value, n, value_range(1), value_range(2), rs)
                end if
            case ('normal', 'gaussian', 'Gaussian')
                call normal_rand_seed_(random_value, n, normal_mu, normal_sigma, rs)
            case ('cauchy', 'Cauchy')
                call cauchy_rand_seed_(random_value, n, cauchy_a, cauchy_b, rs)
            case ('poisson', 'Poisson')
                call poisson_rand_seed_(random_value, n, poisson_mean, rs)
            case ('exp', 'exponential')
                call exponential_rand_seed_(random_value, n, exponential_lambda, rs)
        end select
    else
        select case (distribution)
            case ('uniform')
                if (random_spacing > 0) then
                    call uniform_rand_(random_value, n, 0.0, 1.0)
                    call sort_index(random_value, index)
                    random_value = empty_space*random_value + value_range(1) + random_spacing*linspace(0.0, n - 1.0, n)
                    random_value = random_value(index)
                else
                    call uniform_rand_(random_value, n, value_range(1), value_range(2))
                end if
            case ('normal', 'gaussian', 'Gaussian')
                call normal_rand_(random_value, n, normal_mu, normal_sigma)
            case ('cauchy', 'Cauchy')
                call cauchy_rand_(random_value, n, cauchy_a, cauchy_b)
            case ('poisson', 'Poisson')
                call poisson_rand_(random_value, n, poisson_mean)
            case ('exp', 'exponential')
                call exponential_rand_(random_value, n, exponential_lambda)
        end select
    end if

    r = nTT(random_value)

end function rand_array_1d_

!
!> Return a 2D array of random numbers
!
function rand_array_2d_(n1, n2, dist, seed, range, mu, sigma, a, b, mean, lambda, spacing) result(r)

    integer :: n1, n2
    ! Given seed or not
    integer, optional :: seed
    ! Distribution
    character(len=*), optional :: dist
    ! Uniform distribution -- value range
    TT, dimension(:), optional :: range
    ! Gaussian distribution -- mean and standard deviation
    TTT, optional :: mu, sigma
    ! Cauchy distribution -- parameters a and b
    TTT, optional :: a, b
    ! Possion distribution -- mean
    TTT, optional :: mean
    ! Exponential distribution -- lambda
    TTT, optional :: lambda
    ! Minimum distance
    TT, optional :: spacing
    ! Output
    TT, allocatable, dimension(:, :) :: r

    TTT, allocatable, dimension(:) :: random_value
    TTT, dimension(1:2) :: value_range
    TTT :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
    TTT :: exponential_lambda
    character(len=32) :: distribution
    TTT :: random_spacing
    integer :: rs

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

    if (present(seed)) then
        rs = seed
    else
        rs = -1
    end if

    random_value = rand_array_1d_(n1*n2, distribution, rs, nTT(value_range), normal_mu, normal_sigma, &
        cauchy_a, cauchy_b, poisson_mean, exponential_lambda, nTT(random_spacing))

    r = nTT(reshape(random_value, [n1, n2]))

end function rand_array_2d_

!
!> Return a 3D array of random numbers
!
function rand_array_3d_(n1, n2, n3, dist, seed, range, mu, sigma, a, b, mean, lambda, spacing) result(r)

    integer :: n1, n2, n3
    ! Given seed or not
    integer, optional :: seed
    ! Distribution
    character(len=*), optional :: dist
    ! Uniform distribution -- value range
    TT, dimension(:), optional :: range
    ! Gaussian distribution -- mean and standard deviation
    TTT, optional :: mu, sigma
    ! Cauchy distribution -- parameters a and b
    TTT, optional :: a, b
    ! Possion distribution -- mean
    TTT, optional :: mean
    ! Exponential distribution -- lambda
    TTT, optional :: lambda
    ! Minimum distance
    TT, optional :: spacing
    ! Output
    TT, allocatable, dimension(:, :, :) :: r

    TTT, allocatable, dimension(:) :: random_value
    TTT, dimension(1:2) :: value_range
    TTT :: normal_mu, normal_sigma, cauchy_a, cauchy_b, poisson_mean
    TTT :: exponential_lambda
    character(len=32) :: distribution
    TTT :: random_spacing
    integer :: rs

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

    if (present(seed)) then
        rs = seed
    else
        rs = -1
    end if

    random_value = rand_array_1d_(n1*n2*n3, distribution, rs, nTT(value_range), normal_mu, normal_sigma, &
        cauchy_a, cauchy_b, poisson_mean, exponential_lambda, nTT(random_spacing))

    r = nTT(reshape(random_value, [n1, n2, n3]))

end function rand_array_3d_

#undef uniform_rand_seed_
#undef normal_rand_seed_
#undef cauchy_rand_seed_
#undef poisson_rand_seed_
#undef exponential_rand_seed_
#undef uniform_rand_
#undef normal_rand_
#undef cauchy_rand_
#undef poisson_rand_
#undef exponential_rand_

#undef rand_
#undef rand_array_1d_
#undef rand_array_2d_
#undef rand_array_3d_

#undef T
#undef TT
#undef TTT
#undef nTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef median_filt_1d_
#undef median_filt_2d_
#undef median_filt_3d_
