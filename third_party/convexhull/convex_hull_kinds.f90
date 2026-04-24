!>
!> Numerical kind parameters for the convex_hull package.
!>

module convex_hull_kinds

    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private
    public :: dp
    integer, parameter :: dp = real64

end module convex_hull_kinds
