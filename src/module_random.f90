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
        subroutine uniform_rand_int(w, nw, lb, ub) bind(c, name='uniform_rand')
            use iso_c_binding, only: c_int, c_int
            real(kind=c_int), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_int), value :: lb, ub
        end subroutine uniform_rand_int

        ! random number following uniform distribution, given seed
        subroutine uniform_rand_seed_int(w, nw, lb, ub, seed) bind(c, name='uniform_rand_seed')
            use iso_c_binding, only: c_int, c_int
            real(kind=c_int), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_int), value :: lb, ub
            integer(kind=c_int), value :: seed
        end subroutine uniform_rand_seed_int

        ! random number following normal distribution
        subroutine normal_rand_int(w, nw, mu, sigma) bind(c, name='normal_rand')
            use iso_c_binding, only: c_int, c_int
            real(kind=c_int), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_int), value :: mu, sigma
        end subroutine normal_rand_int

        ! random number following normal distribution, given seed
        subroutine normal_rand_seed_int(w, nw, mu, sigma, seed) bind(c, name='normal_rand_seed')
            use iso_c_binding, only: c_int, c_int
            real(kind=c_int), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_int), value :: mu, sigma
            integer(kind=c_int), value :: seed
        end subroutine normal_rand_seed_int

        ! random number following Cauchy distribution
        subroutine cauchy_rand_int(w, nw, a, b) bind(c, name='cauchy_rand')
            use iso_c_binding, only: c_int, c_int
            real(kind=c_int), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_int), value :: a, b
        end subroutine cauchy_rand_int

        ! random number following Cauchy distribution, given seed
        subroutine cauchy_rand_seed_int(w, nw, a, b, seed) bind(c, name='cauchy_rand_seed')
            use iso_c_binding, only: c_int, c_int
            real(kind=c_int), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_int), value :: a, b
            integer(kind=c_int), value :: seed
        end subroutine cauchy_rand_seed_int

        ! random number following Poisson distribution
        subroutine poisson_rand_int(w, nw, mean) bind(c, name='poisson_rand')
            use iso_c_binding, only: c_int, c_int
            real(kind=c_int), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_int), value :: mean
        end subroutine poisson_rand_int

        ! random number following Poisson distribution, given seed
        subroutine poisson_rand_seed_int(w, nw, mean, seed) bind(c, name='poisson_rand_seed')
            use iso_c_binding, only: c_int, c_int
            real(kind=c_int), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_int), value :: mean
            integer(kind=c_int), value :: seed
        end subroutine poisson_rand_seed_int

        ! random number following exponential distribution
        subroutine exponential_rand_int(w, nw, lambda) bind(c, name='exponential_rand')
            use iso_c_binding, only: c_int, c_int
            real(kind=c_int), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_int), value :: lambda
        end subroutine exponential_rand_int

        ! random number following exponential distribution, given seed
        subroutine exponential_rand_seed_int(w, nw, lambda, seed) bind(c, name='exponential_rand_seed')
            use iso_c_binding, only: c_int, c_int
            real(kind=c_int), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_int), value :: lambda
            integer(kind=c_int), value :: seed
        end subroutine exponential_rand_seed_int

        ! random number following uniform distribution
        subroutine uniform_rand_float(w, nw, lb, ub) bind(c, name='uniform_rand')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: lb, ub
        end subroutine uniform_rand_float

        ! random number following uniform distribution, given seed
        subroutine uniform_rand_seed_float(w, nw, lb, ub, seed) bind(c, name='uniform_rand_seed')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: lb, ub
            integer(kind=c_int), value :: seed
        end subroutine uniform_rand_seed_float

        ! random number following normal distribution
        subroutine normal_rand_float(w, nw, mu, sigma) bind(c, name='normal_rand')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: mu, sigma
        end subroutine normal_rand_float

        ! random number following normal distribution, given seed
        subroutine normal_rand_seed_float(w, nw, mu, sigma, seed) bind(c, name='normal_rand_seed')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: mu, sigma
            integer(kind=c_int), value :: seed
        end subroutine normal_rand_seed_float

        ! random number following Cauchy distribution
        subroutine cauchy_rand_float(w, nw, a, b) bind(c, name='cauchy_rand')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: a, b
        end subroutine cauchy_rand_float

        ! random number following Cauchy distribution, given seed
        subroutine cauchy_rand_seed_float(w, nw, a, b, seed) bind(c, name='cauchy_rand_seed')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: a, b
            integer(kind=c_int), value :: seed
        end subroutine cauchy_rand_seed_float

        ! random number following Poisson distribution
        subroutine poisson_rand_float(w, nw, mean) bind(c, name='poisson_rand')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: mean
        end subroutine poisson_rand_float

        ! random number following Poisson distribution, given seed
        subroutine poisson_rand_seed_float(w, nw, mean, seed) bind(c, name='poisson_rand_seed')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: mean
            integer(kind=c_int), value :: seed
        end subroutine poisson_rand_seed_float

        ! random number following exponential distribution
        subroutine exponential_rand_float(w, nw, lambda) bind(c, name='exponential_rand')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: lambda
        end subroutine exponential_rand_float

        ! random number following exponential distribution, given seed
        subroutine exponential_rand_seed_float(w, nw, lambda, seed) bind(c, name='exponential_rand_seed')
            use iso_c_binding, only: c_float, c_int
            real(kind=c_float), dimension(*) :: w
            integer(kind=c_int), value :: nw
            real(kind=c_float), value :: lambda
            integer(kind=c_int), value :: seed
        end subroutine exponential_rand_seed_float

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

    interface irand
        module procedure :: rand_int
    end interface irand

    interface rand
        module procedure :: rand_float
    end interface

    interface drand
        module procedure :: rand_double
    end interface

    interface irandom
        module procedure :: rand_int
        module procedure :: rand_array_1d_int
        module procedure :: rand_array_2d_int
        module procedure :: rand_array_3d_int
    end interface irandom

    interface random
        module procedure :: rand_float
        module procedure :: rand_array_1d_float
        module procedure :: rand_array_2d_float
        module procedure :: rand_array_3d_float
    end interface random

    interface drandom
        module procedure :: rand_double
        module procedure :: rand_array_1d_double
        module procedure :: rand_array_2d_double
        module procedure :: rand_array_3d_double
    end interface drandom

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
    public :: random_order
    public :: random_string
    public :: random_permute
    public :: random_mask

contains

#define T int
#define TT integer
#define TTT real
#define nTT nint
#include "template_random.f90"

#define T float
#define TT real
#define TTT real
#define nTT
#include "template_random.f90"

#define T double
#define TT double precision
#define TTT double precision
#define nTT
#include "template_random.f90"

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
    function random_string(nc, seed) result(str)

        integer, optional :: nc, seed

        integer :: nchar, i
        character(len=1), dimension(1:72) :: pool
        character(len=:), allocatable :: str
        integer, allocatable, dimension(:) :: ri

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
        if (present(seed)) then
            ri = irandom(nchar, range=[1, 72], seed=seed)
        else
            ri = irandom(nchar, range=[1, 72])
        end if

        do i = 1, nchar
            str(i:i) = pool(ri(i))
        end do

    end function random_string

    !
    !> Generate random mask
    !
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
    !> Random permute an array
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
