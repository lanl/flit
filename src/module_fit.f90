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


module libflit_fit

    use libflit_linear_algebra
    use libflit_string
    use libflit_array
    use iso_fortran_env
    use libflit_array_operation

    implicit none

    private

    type :: polyfit1
        integer :: order = 1
        real :: l1err, l2err, linferr
        real, allocatable, dimension(:) :: coef, misfit
    contains
        procedure :: build => polyfit1_build
        procedure :: eval => polyfit1_evaluate
    end type

    type :: polyfit2
        integer :: order = 1
        real :: l1err, l2err, linferr
        real, allocatable, dimension(:) :: coef, misfit
    contains
        procedure :: build => polyfit2_build
        procedure :: eval => polyfit2_evaluate
    end type

    public :: polyfit1
    public :: polyfit2

contains

    subroutine polyfit1_build(this, x, v, verbose)

        class(polyfit1), intent(inout) :: this
        real, dimension(:), intent(in) :: x, v
        logical, intent(in), optional :: verbose

        real, allocatable, dimension(:, :) :: a
        integer :: n, i, j
        logical :: fit_verbose

        if (present(verbose)) then
            fit_verbose = verbose
        else
            fit_verbose = .false.
        end if

        n = size(x)
        if (size(v) /= n) then
            write (error_unit, *) ' Error: size(x) /= size(v); Exit'
        end if

        ! Form the Vandermonde matrix
        allocate (a(1:n, 1:this%order + 1))
        do i = 1, n
            do j = 1, this%order + 1
                a(i, j) = x(i)**(j - 1)
            end do
        end do

        ! Solve the least-squares problem
        call alloc_array(this%coef, [1, this%order + 1])
        this%coef = lsqsolve(a, v)

        ! Misfit vectors
        call alloc_array(this%misfit, [1, n])
        this%misfit = v - matx(a, this%coef)

        ! Errors
        this%l1err = sum(abs(this%misfit))
        this%l2err = sqrt(sum(this%misfit**2))
        this%linferr = maxval(abs(this%misfit))

        if (fit_verbose) then
            write (error_unit, *)
            write (error_unit, '(a)') ' ------------------- Polynomial Fitting ------------------- '
            write (error_unit, '(a)') '   Form: '
            write (error_unit, '(a)') '      y = sum_{i=0}^{'//num2str(this%order)//'} c_i * x^i'
            write (error_unit, *)
            write (error_unit, '(a)') '   Coefficients: '
            do i = 1, this%order + 1
                write (error_unit, '(a)') '      c_'//num2str(i - 1)//' = '//num2str(this%coef(i), '(es)')
            end do
            write (error_unit, *)
            write (error_unit, '(a)') '   Errors: '
            write (error_unit, '(a)') '      err_l1 = '//num2str(this%l1err, '(es)')
            write (error_unit, '(a)') '      err_l2 = '//num2str(this%l2err, '(es)')
            write (error_unit, '(a)') '      err_linf = '//num2str(this%linferr, '(es)')
            write (error_unit, *)
            write (error_unit, '(a14, x, a14, x, a14, x, a14)') 'X', 'V', 'Fit', 'Misfit'
            write (error_unit, '(a14, x, a14, x, a14, x, a14)') '--------------', '--------------', &
                '--------------', '--------------'
            do i = 1, n
                write (error_unit, '(es14.7, x, es14.7, x, es14.7, x, es14.7)') &
                    x(i), v(i), v(i) - this%misfit(i), this%misfit(i)
            end do
            write (error_unit, '(a14, x, a14, x, a14, x, a14)') '--------------', '--------------', &
                '--------------', '--------------'
            write (error_unit, *)
        end if

        deallocate(a)

    end subroutine polyfit1_build

    function polyfit1_evaluate(this, x) result(y)

        class(polyfit1), intent(in) :: this
        real, intent(in) :: x

        real :: y
        integer :: i

        y = 0.0
        do i = 1, this%order + 1
            y = y + this%coef(i)*x**(i - 1)
        end do

    end function polyfit1_evaluate

    subroutine polyfit2_build(this, x, y, v, verbose)

        class(polyfit2), intent(inout) :: this
        real, dimension(:), intent(in) :: x, y, v
        logical, intent(in), optional :: verbose

        real, allocatable, dimension(:, :) :: a
        integer :: n, i, j, k, l
        logical :: fit_verbose

        if (present(verbose)) then
            fit_verbose = verbose
        else
            fit_verbose = .false.
        end if

        n = size(x)
        if (size(y) /= n) then
            write (error_unit, *) ' Error: size(x) /= size(y); Exit'
        end if
        if (size(v) /= n) then
            write (error_unit, *) ' Error: size(x) /= size(v); Exit'
        end if

        ! Form the Vandermonde matrix
        allocate (a(1:n, 1:(this%order + 1)**2))
        l = 1
        do j = 1, this%order + 1
            do k = 1, this%order + 1
                if (j + k <= this%order + 2) then
                    do i = 1, n
                        a(i, l) = x(i)**(j - 1)*y(i)**(k - 1)
                    end do
                    l = l + 1
                end if
            end do
        end do
        l = l - 1
        call alloc_array(a, [1, n, 1, l], source=a(:, 1:l))

        ! Solve the least-squares problem
        call alloc_array(this%coef, [1, l])
        this%coef = lsqsolve(a, v)

        ! Misfit vectors
        call alloc_array(this%misfit, [1, n])
        this%misfit = v - matx(a, this%coef)

        ! Errors
        this%l1err = sum(abs(this%misfit))
        this%l2err = sqrt(sum(this%misfit**2))
        this%linferr = maxval(abs(this%misfit))

        if (fit_verbose) then
            write (error_unit, *)
            write (error_unit, '(a)') ' ------------------- Polynomial Fitting ------------------- '
            write (error_unit, '(a)') '   Form: '
            write (error_unit, '(a)') '      y = sum_{i=0}^{'//num2str(this%order)//'} sum_{j=0}^{' &
                //num2str(this%order)//'} c_ij * x^i * y^j (i + j <= '//num2str(this%order)//')'
            write (error_unit, *)
            write (error_unit, '(a)') '   Coefficients: '
            l = 1
            do i = 1, this%order + 1
                do j = 1, this%order + 1
                    if (i + j <= this%order + 2) then
                        write (error_unit, '(a)') '      c_'//num2str(i - 1)//'' &
                            //num2str(j - 1)//' = '//num2str(this%coef(l), '(es)')
                        l = l + 1
                    end if
                end do
            end do
            write (error_unit, *)
            write (error_unit, '(a)') '   Errors: '
            write (error_unit, '(a)') '      err_l1 = '//num2str(this%l1err, '(es)')
            write (error_unit, '(a)') '      err_l2 = '//num2str(this%l2err, '(es)')
            write (error_unit, '(a)') '      err_linf = '//num2str(this%linferr, '(es)')
            write (error_unit, *)
            write (error_unit, '(a14, x, a14, x, a14, x, a14, x, a14)') 'X', 'Y', 'V', 'Fit', 'Misfit'
            write (error_unit, '(a14, x, a14, x, a14, x, a14, x, a14)') '--------------', '--------------', &
                '--------------', '--------------', '--------------'
            do i = 1, n
                write (error_unit, '(es14.7, x, es14.7, x, es14.7, x, es14.7, x, es14.7)') &
                    x(i), y(i), v(i), v(i) - this%misfit(i), this%misfit(i)
            end do
            write (error_unit, '(a14, x, a14, x, a14, x, a14, x, a14)') '--------------', '--------------', &
                '--------------', '--------------', '--------------'
            write (error_unit, *)
        end if

        deallocate(a)

    end subroutine polyfit2_build

    function polyfit2_evaluate(this, x, y) result(z)

        class(polyfit2), intent(in) :: this
        real, intent(in) :: x, y

        real :: z
        integer :: i, j, l

        z = 0.0
        l = 1
        do i = 1, this%order + 1
            do j = 1, this%order + 1
                if (i + j <= this%order + 2) then
                    z = z + this%coef(l)*x**(i - 1)*y**(j - 1)
                    l = l + 1
                end if
            end do
        end do

    end function polyfit2_evaluate

end module libflit_fit
