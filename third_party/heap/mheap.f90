!
! Copyright (c) 2014, Daniel Pena
! All rights reserved.

! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:

! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.

! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!
! Note by Kai Gao:
!   This is a reformatted and modified version of the original code
!   to accommodate libflit's template programming.
!   It must be inserted into libflit compilation flow to work.
!
module mheap_

    implicit none

    abstract interface
        function heapfunc_(node1, node2) result(res)
            TT, intent(in) :: node1(:)
            TT, intent(in) :: node2(:)
            logical :: res
        end function heapfunc_
    end interface

    type :: heap_
        integer :: nmax ! max size
        integer :: n ! current heap size
        integer :: m ! current tree size
        integer :: nlen ! node size in real  units
        TT, allocatable :: data(:, :) ! node data
        integer, allocatable :: indx(:) ! nodes index
        procedure(heapfunc_), nopass, pointer :: fun ! heap function to find root node
    contains
        procedure :: init => heap_init_
        procedure :: insert => heap_insert_
        procedure :: peek => heap_peek_
        procedure :: pop => heap_pop_
        procedure :: reheap => heap_reheap_
        procedure :: size => heap_size_
        procedure :: free => heap_free_
    end type heap_

    private
    public :: heap_

contains

    ! returns the heap current size
    integer function heap_size_(heap)

        class(heap_) :: heap

        heap_size_ = heap%n

    end function heap_size_

    ! initializes the heap
    ! nmax  -  max size of the heap
    ! nlen  -  size of each node
    ! hpfun -  the heap function (provides comparison between two nodes' data)
    subroutine heap_init_(heap, nmax, nlen, hpfun)

        class(heap_) :: heap
        integer, intent(in) :: nmax, nlen
        procedure(heapfunc_) :: hpfun
        integer :: i

        heap%nmax = nmax
        heap%n = 0
        heap%m = 0
        heap%nlen = nlen
        heap%fun => hpfun

        allocate (heap%indx(nmax))
        allocate (heap%data(nlen, nmax))

        do i = 1, nmax
            heap%indx(i) = i
        end do

    end subroutine heap_init_

    ! releases all the allocated memory and resets the heap
    subroutine heap_free_(heap)

        class(heap_) :: heap

        deallocate (heap%indx)
        deallocate (heap%data)
        heap%n = 0
        heap%m = 0
        heap%nmax = 0
        heap%fun => null()

    end subroutine heap_free_

    ! insert a node into a heap. the resulting tree is re-heaped.
    !  input
    !        heap - the heap
    !        node - a real  array, nlen long, which
    !               contains the node's information to be inserted.
    subroutine heap_insert_(heap, node)

        class(heap_) :: heap
        TT, intent(in) :: node(heap%nlen)

        integer :: k1, k2, il, ir

        if (heap%n == heap%nmax) then
            return
        end if

        ! add one element and copy node data to new element
        heap%n = heap%n + 1
        heap%m = heap%m + 1
        heap%data(:, heap%indx(heap%n)) = node(:)

        ! re-index the heap from the bottom up
        k2 = heap%n
        do while (k2 /= 1)
            k1 = k2/2
            ir = heap%indx(k2)
            il = heap%indx(k1)
            if (heap%fun(heap%data(:, il), heap%data(:, ir))) then
                return
            end if
            call swapint(heap%indx(k2), heap%indx(k1))
            k2 = k2/2
        end do

    end subroutine heap_insert_

    ! retrieve the root element off the heap. the resulting tree is re-heaped.
    ! no data is deleted, thus the original
    !   input
    !        heap - the heap
    !   output
    !        node - the deleted node
    subroutine heap_pop_(heap, node)

        class(heap_) :: heap
        TT, optional :: node(heap%nlen)

        if (heap%n == 0) then
            return
        end if

        if (present(node)) then
            node(:) = heap%data(:, heap%indx(1))
        end if

        call swapint(heap%indx(1), heap%indx(heap%n))

        heap%n = heap%n - 1

        call heap_grow_(heap, 1)

    end subroutine heap_pop_

    ! access the k-th node of the heap
    subroutine heap_peek_(heap, k, node)

        class(heap_) :: heap
        integer, intent(in) :: k
        TT, intent(out) :: node(heap%nlen)

        if (k < 1 .or. k > heap%n .or. heap%n > heap%nmax) then
            return
        end if

        node(:) = heap%data(:, heap%indx(k))

    end subroutine heap_peek_

    ! forms a heap out of a tree. used privately by heap_reheap.
    ! the root node of the tree is stored in the location indx(ktemp).
    ! the first child node is in location indx(2*ktemp)...
    ! the next child node is in location indx(2*ktemp+1).
    ! this subroutines assumes each branch of the tree is itself a heap.
    subroutine heap_grow_(heap, ktemp)

        integer :: i, k, il, ir
        type(heap_) :: heap
        integer :: ktemp

        if (heap%n > heap%nmax) return

        k = ktemp
        do while (2*k <= heap%n)

            i = 2*k

            ! if there is more than one child node, find which is the smallest.
            if (2*k /= heap%n) then
                il = heap%indx(2*k + 1)
                ir = heap%indx(2*k)
                if (heap%fun(heap%data(:, il), heap%data(:, ir))) then
                    i = i + 1
                end if
            end if

            ! if a child is larger than its parent, interchange them... this destroys
            ! the heap property, so the remaining elements must be re-heaped.
            il = heap%indx(k)
            ir = heap%indx(i)
            if (heap%fun(heap%data(:, il), heap%data(:, ir))) then
                return
            end if

            call swapint(heap%indx(i), heap%indx(k))

            k = i

        end do

    end subroutine heap_grow_

    !
    ! builds the heap from the element data using the provided heap function.
    ! at exit, the root node satisfies the heap condition:
    !   hpfun( root_node, node ) = .true. for any other node
    !
    subroutine heap_reheap_(heap, hpfun)

        class(heap_) :: heap
        procedure(heapfunc_), optional :: hpfun
        integer :: k

        heap%n = heap%m
        if (present(hpfun)) then
            heap%fun => hpfun
        end if

        if (heap%nmax < heap%n) then
            return
        end if

        do k = heap%n/2, 1, -1
            call heap_grow_(heap, k)
        end do

    end subroutine heap_reheap_

    subroutine swapint(i, k)

        integer :: i, k, t

        t = i
        i = k
        k = t

    end subroutine swapint

end module mheap_

