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

#define libflit_meta_array_     CONCAT(libflit_meta_array, T)

#define meta_array1_     CONCAT(meta_array1, T)
#define meta_array2_     CONCAT(meta_array2, T)
#define meta_array3_     CONCAT(meta_array3, T)
#define meta_array4_     CONCAT(meta_array4, T)

#define list_meta_array1_     CONCAT(list_meta_array1, T)
#define list_meta_array2_     CONCAT(list_meta_array2, T)
#define list_meta_array3_     CONCAT(list_meta_array3, T)
#define list_meta_array4_     CONCAT(list_meta_array4, T)

#define get_meta_array1_     CONCAT(get_meta_array1, T)
#define get_meta_array2_     CONCAT(get_meta_array2, T)
#define get_meta_array3_     CONCAT(get_meta_array3, T)
#define get_meta_array4_     CONCAT(get_meta_array4, T)

#define get_meta_array_core1_     CONCAT(get_meta_array_core1, T)
#define get_meta_array_core2_     CONCAT(get_meta_array_core2, T)
#define get_meta_array_core3_     CONCAT(get_meta_array_core3, T)
#define get_meta_array_core4_     CONCAT(get_meta_array_core4, T)

#define set_meta_array_core1_     CONCAT(set_meta_array_core1, T)
#define set_meta_array_core2_     CONCAT(set_meta_array_core2, T)
#define set_meta_array_core3_     CONCAT(set_meta_array_core3, T)
#define set_meta_array_core4_     CONCAT(set_meta_array_core4, T)

module libflit_meta_array_

    use libflit_array
    use libflit_error
    use libflit_string

    implicit none

    type meta_array1_

        character(len=24) :: name
        integer :: id
        integer :: n1
        TT, allocatable, dimension(:) :: array

    contains
        generic, public :: assignment(=) => copy_meta_array1_
        procedure, private :: copy_meta_array1_
        procedure, public :: init => initialize_meta_array1_

    end type meta_array1_

    type meta_array2_

        character(len=24) :: name
        integer :: id
        integer :: n1, n2
        TT, allocatable, dimension(:, :) :: array

    contains
        generic, public :: assignment(=) => copy_meta_array2_
        procedure, private :: copy_meta_array2_
        procedure, public :: init => initialize_meta_array2_

    end type meta_array2_

    type meta_array3_

        character(len=24) :: name
        integer :: id
        integer :: n1, n2, n3
        TT, allocatable, dimension(:, :, :) :: array

    contains
        generic, public :: assignment(=) => copy_meta_array3_
        procedure, private :: copy_meta_array3_
        procedure, public :: init => initialize_meta_array3_

    end type meta_array3_

    type meta_array4_

        character(len=24) :: name
        integer :: id
        integer :: n1, n2, n3, n4
        TT, allocatable, dimension(:, :, :, :) :: array

    contains
        generic, public :: assignment(=) => copy_meta_array4_
        procedure, private :: copy_meta_array4_
        procedure, public :: init => initialize_meta_array4_

    end type meta_array4_

    type list_meta_array1_

        integer :: n
        type(meta_array1_), allocatable, dimension(:) :: meta_array

    contains
        procedure, public :: get => get_from_list_meta_array1_
        procedure, public :: init => initialize_list_meta_array1_
        procedure, public :: size => size_of_list_meta_array1_
        procedure, public :: has => probe_existence_list_meta_array1_
        procedure, public :: get_core => get_core_from_list_meta_array1_

    end type list_meta_array1_

    type list_meta_array2_

        integer :: n
        type(meta_array2_), allocatable, dimension(:) :: meta_array

    contains
        procedure, public :: get => get_from_list_meta_array2_
        procedure, public :: init => initialize_list_meta_array2_
        procedure, public :: size => size_of_list_meta_array2_
        procedure, public :: has => probe_existence_list_meta_array2_
        procedure, public :: get_core => get_core_from_list_meta_array2_

    end type list_meta_array2_

    type list_meta_array3_

        integer :: n
        type(meta_array3_), allocatable, dimension(:) :: meta_array

    contains
        procedure, public :: get => get_from_list_meta_array3_
        procedure, public :: init => initialize_list_meta_array3_
        procedure, public :: size => size_of_list_meta_array3_
        procedure, public :: has => probe_existence_list_meta_array3_
        procedure, public :: get_core => get_core_from_list_meta_array3_

    end type list_meta_array3_

    type list_meta_array4_

        integer :: n
        type(meta_array4_), allocatable, dimension(:) :: meta_array

    contains
        procedure, public :: get => get_from_list_meta_array4_
        procedure, public :: init => initialize_list_meta_array4_
        procedure, public :: size => size_of_list_meta_array4_
        procedure, public :: has => probe_existence_list_meta_array4_
        procedure, public :: get_core => get_core_from_list_meta_array4_

    end type list_meta_array4_

    private
    public :: meta_array1_
    public :: meta_array2_
    public :: meta_array3_
    public :: meta_array4_
    public :: list_meta_array1_
    public :: list_meta_array2_
    public :: list_meta_array3_
    public :: list_meta_array4_
    public :: get_meta_array1_
    public :: get_meta_array2_
    public :: get_meta_array3_
    public :: get_meta_array4_
    public :: get_meta_array_core1_
    public :: get_meta_array_core2_
    public :: get_meta_array_core3_
    public :: get_meta_array_core4_
    public :: set_meta_array_core1_
    public :: set_meta_array_core2_
    public :: set_meta_array_core3_
    public :: set_meta_array_core4_

contains

    subroutine set_meta_array_core1_(w, name, x)

        type(meta_array1_), dimension(:), intent(inout) :: w
        character(len=*), intent(in) :: name
        TT, dimension(:), intent(in) :: x

        integer :: i

        do i = 1, size(w)
            if (tidy(w(i)%name) == tidy(name)) then
                w(i)%array = x
                exit
            end if
        end do

    end subroutine set_meta_array_core1_

    subroutine set_meta_array_core2_(w, name, x)

        type(meta_array2_), dimension(:), intent(inout) :: w
        character(len=*), intent(in) :: name
        TT, dimension(:, :), intent(in) :: x

        integer :: i

        do i = 1, size(w)
            if (tidy(w(i)%name) == tidy(name)) then
                w(i)%array = x
                exit
            end if
        end do

    end subroutine set_meta_array_core2_

    subroutine set_meta_array_core3_(w, name, x)

        type(meta_array3_), dimension(:), intent(inout) :: w
        character(len=*), intent(in) :: name
        TT, dimension(:, :, :), intent(in) :: x

        integer :: i

        do i = 1, size(w)
            if (tidy(w(i)%name) == tidy(name)) then
                w(i)%array = x
                exit
            end if
        end do

    end subroutine set_meta_array_core3_

    subroutine set_meta_array_core4_(w, name, x)

        type(meta_array4_), dimension(:), intent(inout) :: w
        character(len=*), intent(in) :: name
        TT, dimension(:, :, :, :), intent(in) :: x

        integer :: i

        do i = 1, size(w)
            if (tidy(w(i)%name) == tidy(name)) then
                w(i)%array = x
                exit
            end if
        end do

    end subroutine set_meta_array_core4_

    function get_meta_array_core1_(w, name) result(x)

        type(meta_array1_), dimension(:), intent(in) :: w
        character(len=*), intent(in) :: name
        TT, allocatable, dimension(:) :: x

        integer :: i

        do i = 1, size(w)
            if (tidy(w(i)%name) == tidy(name)) then
                x = w(i)%array
                exit
            end if
        end do

    end function get_meta_array_core1_

    function get_meta_array_core2_(w, name) result(x)

        type(meta_array2_), dimension(:), intent(in) :: w
        character(len=*), intent(in) :: name
        TT, allocatable, dimension(:, :) :: x

        integer :: i

        do i = 1, size(w)
            if (tidy(w(i)%name) == tidy(name)) then
                x = w(i)%array
                exit
            end if
        end do

    end function get_meta_array_core2_

    function get_meta_array_core3_(w, name) result(x)

        type(meta_array3_), dimension(:), intent(in) :: w
        character(len=*), intent(in) :: name
        TT, allocatable, dimension(:, :, :) :: x

        integer :: i

        do i = 1, size(w)
            if (tidy(w(i)%name) == tidy(name)) then
                x = w(i)%array
                exit
            end if
        end do

    end function get_meta_array_core3_

    function get_meta_array_core4_(w, name) result(x)

        type(meta_array4_), dimension(:), intent(in) :: w
        character(len=*), intent(in) :: name
        TT, allocatable, dimension(:, :, :, :) :: x

        integer :: i

        do i = 1, size(w)
            if (tidy(w(i)%name) == tidy(name)) then
                x = w(i)%array
                exit
            end if
        end do

    end function get_meta_array_core4_

    function get_meta_array1_(w, name) result(x)

        type(meta_array1_), dimension(:), intent(in) :: w
        character(len=*), intent(in) :: name
        type(meta_array1_) :: x

        integer :: i

        do i = 1, size(w)
            if (tidy(w(i)%name) == tidy(name)) then
                x = w(i)
                exit
            end if
        end do

    end function get_meta_array1_

    function get_meta_array2_(w, name) result(x)

        type(meta_array2_), dimension(:), intent(in) :: w
        character(len=*), intent(in) :: name
        type(meta_array2_) :: x

        integer :: i

        do i = 1, size(w)
            if (tidy(w(i)%name) == tidy(name)) then
                x = w(i)
                exit
            end if
        end do

    end function get_meta_array2_

    function get_meta_array3_(w, name) result(x)

        type(meta_array3_), dimension(:), intent(in) :: w
        character(len=*), intent(in) :: name
        type(meta_array3_) :: x

        integer :: i

        do i = 1, size(w)
            if (tidy(w(i)%name) == tidy(name)) then
                x = w(i)
                exit
            end if
        end do

    end function get_meta_array3_

    function get_meta_array4_(w, name) result(x)

        type(meta_array4_), dimension(:), intent(in) :: w
        character(len=*), intent(in) :: name
        type(meta_array4_) :: x

        integer :: i

        do i = 1, size(w)
            if (tidy(w(i)%name) == tidy(name)) then
                x = w(i)
                exit
            end if
        end do

    end function get_meta_array4_

    !
    !> Initialize meta array
    !
    subroutine initialize_meta_array1_(this, n1, name, id, value)

        class(meta_array1_), intent(inout) :: this
        integer, intent(in) :: n1
        character(len=*), intent(in), optional :: name
        integer, intent(in), optional :: id
        TT, intent(in), optional :: value

        if (allocated(this%array)) then
            deallocate(this%array)
        end if

        this%n1 = n1
        allocate(this%array(1:this%n1))

        if (present(name)) then
            this%name = tidy(name)
        else
            this%name = 'meta array'
        end if

        if (present(id)) then
            this%id = id
        else
            this%id = 0
        end if

        if (present(value)) then
            this%array = value
        else
            this%array = 0
        end if

    end subroutine initialize_meta_array1_

    subroutine initialize_meta_array2_(this, n1, n2, name, id, value)

        class(meta_array2_), intent(inout) :: this
        integer, intent(in) :: n1, n2
        character(len=*), intent(in), optional :: name
        integer, intent(in), optional :: id
        TT, intent(in), optional :: value

        if (allocated(this%array)) then
            deallocate(this%array)
        end if

        this%n1 = n1
        this%n2 = n2
        allocate(this%array(1:this%n1, 1:this%n2))

        if (present(name)) then
            this%name = tidy(name)
        else
            this%name = 'meta array'
        end if

        if (present(id)) then
            this%id = id
        else
            this%id = 0
        end if

        if (present(value)) then
            this%array = value
        else
            this%array = 0
        end if

    end subroutine initialize_meta_array2_

    subroutine initialize_meta_array3_(this, n1, n2, n3, name, id, value)

        class(meta_array3_), intent(inout) :: this
        integer, intent(in) :: n1, n2, n3
        character(len=*), intent(in), optional :: name
        integer, intent(in), optional :: id
        TT, intent(in), optional :: value

        if (allocated(this%array)) then
            deallocate(this%array)
        end if

        this%n1 = n1
        this%n2 = n2
        this%n3 = n3
        allocate(this%array(1:this%n1, 1:this%n2, 1:this%n3))

        if (present(name)) then
            this%name = tidy(name)
        else
            this%name = 'meta array'
        end if

        if (present(id)) then
            this%id = id
        else
            this%id = 0
        end if

        if (present(value)) then
            this%array = value
        else
            this%array = 0
        end if

    end subroutine initialize_meta_array3_

    subroutine initialize_meta_array4_(this, n1, n2, n3, n4, name, id, value)

        class(meta_array4_), intent(inout) :: this
        integer, intent(in) :: n1, n2, n3, n4
        character(len=*), intent(in), optional :: name
        integer, intent(in), optional :: id
        TT, intent(in), optional :: value

        if (allocated(this%array)) then
            deallocate(this%array)
        end if

        this%n1 = n1
        this%n2 = n2
        this%n3 = n3
        this%n4 = n4
        allocate(this%array(1:this%n1, 1:this%n2, 1:this%n3, 1:this%n4))

        if (present(name)) then
            this%name = tidy(name)
        else
            this%name = 'meta array'
        end if

        if (present(id)) then
            this%id = id
        else
            this%id = 0
        end if

        if (present(value)) then
            this%array = value
        else
            this%array = 0
        end if

    end subroutine initialize_meta_array4_

    !
    !> Copy (=) meta array
    !
    subroutine copy_meta_array1_(this, from)

        class(meta_array1_), intent(inout) :: this
        type(meta_array1_), intent(in) :: from

        this%name = from%name
        this%id = from%id
        this%n1 = from%n1

        if (allocated(this%array)) then
            deallocate(this%array)
        end if
        allocate(this%array(lbound(from%array, 1):ubound(from%array, 1)))
        this%array = from%array

    end subroutine copy_meta_array1_

    subroutine copy_meta_array2_(this, from)

        class(meta_array2_), intent(inout) :: this
        type(meta_array2_), intent(in) :: from

        this%name = from%name
        this%id = from%id
        this%n1 = from%n1
        this%n2 = from%n2

        if (allocated(this%array)) then
            deallocate(this%array)
        end if
        allocate(this%array(lbound(from%array, 1):ubound(from%array, 1), &
            lbound(from%array, 2):ubound(from%array, 2)))
        this%array = from%array

    end subroutine copy_meta_array2_

    subroutine copy_meta_array3_(this, from)

        class(meta_array3_), intent(inout) :: this
        type(meta_array3_), intent(in) :: from

        this%name = from%name
        this%id = from%id
        this%n1 = from%n1
        this%n2 = from%n2
        this%n3 = from%n3

        if (allocated(this%array)) then
            deallocate(this%array)
        end if
        allocate(this%array(lbound(from%array, 1):ubound(from%array, 1), &
            lbound(from%array, 2):ubound(from%array, 2), &
            lbound(from%array, 3):ubound(from%array, 3)))
        this%array = from%array

    end subroutine copy_meta_array3_

    subroutine copy_meta_array4_(this, from)

        class(meta_array4_), intent(inout) :: this
        type(meta_array4_), intent(in) :: from

        this%name = from%name
        this%id = from%id
        this%n1 = from%n1
        this%n2 = from%n2
        this%n3 = from%n3
        this%n4 = from%n4

        if (allocated(this%array)) then
            deallocate(this%array)
        end if
        allocate(this%array(lbound(from%array, 1):ubound(from%array, 1), &
            lbound(from%array, 2):ubound(from%array, 2), &
            lbound(from%array, 3):ubound(from%array, 3), &
            lbound(from%array, 4):ubound(from%array, 4)))
        this%array = from%array

    end subroutine copy_meta_array4_

    !
    !> Initialize a list of meta array structure
    !
    subroutine initialize_list_meta_array1_(this, n, source)

        class(list_meta_array1_), intent(inout) :: this
        integer, intent(in) :: n
        type(meta_array1_), intent(in), optional :: source

        integer :: i

        this%n = n

        if (allocated(this%meta_array)) then
            deallocate(this%meta_array)
        end if

        allocate(this%meta_array(1:this%n))

        if (present(source)) then
            do i = 1, this%n
                this%meta_array(i) = source
            end do
        end if

    end subroutine initialize_list_meta_array1_

    subroutine initialize_list_meta_array2_(this, n, source)

        class(list_meta_array2_), intent(inout) :: this
        integer, intent(in) :: n
        type(meta_array2_), intent(in), optional :: source

        integer :: i

        this%n = n

        if (allocated(this%meta_array)) then
            deallocate(this%meta_array)
        end if

        allocate(this%meta_array(1:this%n))

        if (present(source)) then
            do i = 1, this%n
                this%meta_array(i) = source
            end do
        end if

    end subroutine initialize_list_meta_array2_

    subroutine initialize_list_meta_array3_(this, n, source)

        class(list_meta_array3_), intent(inout) :: this
        integer, intent(in) :: n
        type(meta_array3_), intent(in), optional :: source

        integer :: i

        this%n = n

        if (allocated(this%meta_array)) then
            deallocate(this%meta_array)
        end if

        allocate(this%meta_array(1:this%n))

        if (present(source)) then
            do i = 1, this%n
                this%meta_array(i) = source
            end do
        end if

    end subroutine initialize_list_meta_array3_

    subroutine initialize_list_meta_array4_(this, n, source)

        class(list_meta_array4_), intent(inout) :: this
        integer, intent(in) :: n
        type(meta_array4_), intent(in), optional :: source

        integer :: i

        this%n = n

        if (allocated(this%meta_array)) then
            deallocate(this%meta_array)
        end if

        allocate(this%meta_array(1:this%n))

        if (present(source)) then
            do i = 1, this%n
                this%meta_array(i) = source
            end do
        end if

    end subroutine initialize_list_meta_array4_

    !
    !> Get a meta array from a list of meta array structure
    !
    function get_from_list_meta_array1_(this, name) result(m)

        class(list_meta_array1_), intent(in) :: this
        character(len=*), intent(in) :: name
        type(meta_array1_) :: m

        integer :: i

        do i = 1, this%n
            if (tidy(this%meta_array(i)%name) == tidy(name)) then
                m = this%meta_array(i)
                return
            end if
        end do

    end function get_from_list_meta_array1_

    function get_from_list_meta_array2_(this, name) result(m)

        class(list_meta_array2_), intent(in) :: this
        character(len=*), intent(in) :: name
        type(meta_array2_) :: m

        integer :: i

        do i = 1, this%n
            if (tidy(this%meta_array(i)%name) == tidy(name)) then
                m = this%meta_array(i)
                return
            end if
        end do

    end function get_from_list_meta_array2_

    function get_from_list_meta_array3_(this, name) result(m)

        class(list_meta_array3_), intent(in) :: this
        character(len=*), intent(in) :: name
        type(meta_array3_) :: m

        integer :: i

        do i = 1, this%n
            if (tidy(this%meta_array(i)%name) == tidy(name)) then
                m = this%meta_array(i)
                return
            end if
        end do

    end function get_from_list_meta_array3_

    function get_from_list_meta_array4_(this, name) result(m)

        class(list_meta_array4_), intent(in) :: this
        character(len=*), intent(in) :: name
        type(meta_array4_) :: m

        integer :: i

        do i = 1, this%n
            if (tidy(this%meta_array(i)%name) == tidy(name)) then
                m = this%meta_array(i)
                return
            end if
        end do

    end function get_from_list_meta_array4_

    !
    !> Probe the shape of a meta array
    !
    function size_of_list_meta_array1_(this, name, dim) result(n)

        class(list_meta_array1_), intent(in) :: this
        character(len=*), intent(in) :: name
        integer, intent(in) :: dim
        integer :: n

        integer :: i

        do i = 1, this%n
            if (this%meta_array(i)%name == name) then
                n = size(this%meta_array(i)%array, dim=dim)
                return
            end if
        end do

    end function size_of_list_meta_array1_

    function size_of_list_meta_array2_(this, name, dim) result(n)

        class(list_meta_array2_), intent(in) :: this
        character(len=*), intent(in) :: name
        integer, intent(in) :: dim
        integer :: n

        integer :: i

        do i = 1, this%n
            if (this%meta_array(i)%name == name) then
                n = size(this%meta_array(i)%array, dim=dim)
                return
            end if
        end do

    end function size_of_list_meta_array2_

    function size_of_list_meta_array3_(this, name, dim) result(n)

        class(list_meta_array3_), intent(in) :: this
        character(len=*), intent(in) :: name
        integer, intent(in) :: dim
        integer :: n

        integer :: i

        do i = 1, this%n
            if (this%meta_array(i)%name == name) then
                n = size(this%meta_array(i)%array, dim=dim)
                return
            end if
        end do

    end function size_of_list_meta_array3_

    function size_of_list_meta_array4_(this, name, dim) result(n)

        class(list_meta_array4_), intent(in) :: this
        character(len=*), intent(in) :: name
        integer, intent(in) :: dim
        integer :: n

        integer :: i

        do i = 1, this%n
            if (this%meta_array(i)%name == name) then
                n = size(this%meta_array(i)%array, dim=dim)
                return
            end if
        end do

    end function size_of_list_meta_array4_

    !
    !> Probe if has a meta array with a specific name
    !
    function probe_existence_list_meta_array1_(this, name) result(n)

        class(list_meta_array1_), intent(in) :: this
        character(len=*), intent(in) :: name
        logical :: n

        if (allocated(this%meta_array)) then
            n = any(this%meta_array(:)%name == name)
        else
            n = .false.
        end if

    end function probe_existence_list_meta_array1_

    function probe_existence_list_meta_array2_(this, name) result(n)

        class(list_meta_array2_), intent(in) :: this
        character(len=*), intent(in) :: name
        logical :: n

        if (allocated(this%meta_array)) then
            n = any(this%meta_array(:)%name == name)
        else
            n = .false.
        end if

    end function probe_existence_list_meta_array2_

    function probe_existence_list_meta_array3_(this, name) result(n)

        class(list_meta_array3_), intent(in) :: this
        character(len=*), intent(in) :: name
        logical :: n

        if (allocated(this%meta_array)) then
            n = any(this%meta_array(:)%name == name)
        else
            n = .false.
        end if

    end function probe_existence_list_meta_array3_

    function probe_existence_list_meta_array4_(this, name) result(n)

        class(list_meta_array4_), intent(in) :: this
        character(len=*), intent(in) :: name
        logical :: n

        if (allocated(this%meta_array)) then
            n = any(this%meta_array(:)%name == name)
        else
            n = .false.
        end if

    end function probe_existence_list_meta_array4_

    !
    !> Get a meta array from a list of meta array structure
    !
    function get_core_from_list_meta_array1_(this, name) result(m)

        class(list_meta_array1_), intent(in) :: this
        character(len=*), intent(in) :: name
        TT, allocatable, dimension(:) :: m

        integer :: i

        do i = 1, this%n
            if (tidy(this%meta_array(i)%name) == tidy(name)) then
                allocate(m(lbound(this%meta_array(i)%array, 1):ubound(this%meta_array(i)%array, 1)))
                m = this%meta_array(i)%array
                return
            end if
        end do

    end function get_core_from_list_meta_array1_

    function get_core_from_list_meta_array2_(this, name) result(m)

        class(list_meta_array2_), intent(in) :: this
        character(len=*), intent(in) :: name
        TT, allocatable, dimension(:, :) :: m

        integer :: i

        do i = 1, this%n
            if (tidy(this%meta_array(i)%name) == tidy(name)) then
                allocate(m(lbound(this%meta_array(i)%array, 1):ubound(this%meta_array(i)%array, 1), &
                    lbound(this%meta_array(i)%array, 2):ubound(this%meta_array(i)%array, 2)))
                m = this%meta_array(i)%array
                return
            end if
        end do

    end function get_core_from_list_meta_array2_

    function get_core_from_list_meta_array3_(this, name) result(m)

        class(list_meta_array3_), intent(in) :: this
        character(len=*), intent(in) :: name
        TT, allocatable, dimension(:, :, :) :: m

        integer :: i

        do i = 1, this%n
            if (tidy(this%meta_array(i)%name) == tidy(name)) then
                allocate(m(lbound(this%meta_array(i)%array, 1):ubound(this%meta_array(i)%array, 1), &
                    lbound(this%meta_array(i)%array, 2):ubound(this%meta_array(i)%array, 2), &
                    lbound(this%meta_array(i)%array, 3):ubound(this%meta_array(i)%array, 3)))
                m = this%meta_array(i)%array
                return
            end if
        end do

    end function get_core_from_list_meta_array3_

    function get_core_from_list_meta_array4_(this, name) result(m)

        class(list_meta_array4_), intent(in) :: this
        character(len=*), intent(in) :: name
        TT, allocatable, dimension(:, :, :, :) :: m

        integer :: i

        do i = 1, this%n
            if (tidy(this%meta_array(i)%name) == tidy(name)) then
                allocate(m(lbound(this%meta_array(i)%array, 1):ubound(this%meta_array(i)%array, 1), &
                    lbound(this%meta_array(i)%array, 2):ubound(this%meta_array(i)%array, 2), &
                    lbound(this%meta_array(i)%array, 3):ubound(this%meta_array(i)%array, 3), &
                    lbound(this%meta_array(i)%array, 4):ubound(this%meta_array(i)%array, 4)))
                m = this%meta_array(i)%array
                return
            end if
        end do

    end function get_core_from_list_meta_array4_

end module libflit_meta_array_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef libflit_meta_array_

#undef meta_array1_
#undef meta_array2_
#undef meta_array3_
#undef meta_array4_

#undef list_meta_array1_
#undef list_meta_array2_
#undef list_meta_array3_
#undef list_meta_array4_

#undef get_meta_array1_
#undef get_meta_array2_
#undef get_meta_array3_
#undef get_meta_array4_

#undef get_meta_array_core1_
#undef get_meta_array_core2_
#undef get_meta_array_core3_
#undef get_meta_array_core4_

#undef set_meta_array_core1_
#undef set_meta_array_core2_
#undef set_meta_array_core3_
#undef set_meta_array_core4_
