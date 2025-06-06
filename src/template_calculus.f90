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

#define differentiate_1d_      CONCAT(differentiate_1d, T)
#define differentiate_2d_      CONCAT(differentiate_2d, T)
#define differentiate_3d_      CONCAT(differentiate_3d, T)
#define integrate_1d_      CONCAT(integrate_1d, T)
#define integrate_2d_      CONCAT(integrate_2d, T)
#define integrate_3d_      CONCAT(integrate_3d, T)
#define cumsum_1d_      CONCAT(cumsum_1d, T)
#define cumsum_2d_      CONCAT(cumsum_2d, T)
#define cumsum_3d_      CONCAT(cumsum_3d, T)

!
!> Differentiate 1D array using finite difference
!
function differentiate_1d_(w, order, method, accuracy) result(wt)

    TT, dimension(:), intent(in) :: w
    integer, intent(in), optional :: order, accuracy
    character(len=*), intent(in), optional :: method
    TT, allocatable, dimension(:) :: wt

    integer :: i, n
    character(len=12) :: deriv_method
    integer :: deriv_order, deriv_accuracy

    if (present(order)) then
        deriv_order = order
    else
        deriv_order = 1
    end if

    if (present(accuracy)) then
        deriv_accuracy = accuracy
    else
        deriv_accuracy = 1
    end if

    if (present(method)) then
        deriv_method = trim(adjustl((method)))
    else
        deriv_method = 'center'
    end if

    n = size(w)
    allocate (wt(1:n))

    select case (deriv_order)

        case (1)

            select case (deriv_method)

                case ('forward')
                    select case (deriv_accuracy)
                        case (1)
                            do i = 1, n - 1
                                wt(i) = sum(fdcoef_1f10*w(i:i + 1))
                            end do
                            wt(n) = w(n) - w(n - 1)
                        case (2)
                            do i = 1, n - 2
                                wt(i) = sum(fdcoef_1f20*w(i:i + 2))
                            end do
                            wt(n - 1) = sum(fdcoef_1f10*w(n - 1:n))
                            wt(n) = w(n) - w(n - 1)
                        case (3)
                            do i = 1, n - 3
                                wt(i) = sum(fdcoef_1f30*w(i:i + 3))
                            end do
                            wt(n - 2) = sum(fdcoef_1f20*w(n - 2:n))
                            wt(n - 1) = sum(fdcoef_1f10*w(n - 1:n))
                            wt(n) = w(n) - w(n - 1)
                        case (4)
                            do i = 1, n - 4
                                wt(i) = sum(fdcoef_1f40*w(i:i + 4))
                            end do
                            wt(n - 3) = sum(fdcoef_1f30*w(n - 3:n))
                            wt(n - 2) = sum(fdcoef_1f20*w(n - 2:n))
                            wt(n - 1) = sum(fdcoef_1f10*w(n - 1:n))
                            wt(n) = w(n) - w(n - 1)
                        case (5)
                            do i = 1, n - 5
                                wt(i) = sum(fdcoef_1f50*w(i:i + 5))
                            end do
                            wt(n - 4) = sum(fdcoef_1f40*w(n - 4:n))
                            wt(n - 3) = sum(fdcoef_1f30*w(n - 3:n))
                            wt(n - 2) = sum(fdcoef_1f20*w(n - 2:n))
                            wt(n - 1) = sum(fdcoef_1f10*w(n - 1:n))
                            wt(n) = w(n) - w(n - 1)
                        case (6)
                            do i = 1, n - 6
                                wt(i) = sum(fdcoef_1f60*w(i:i + 6))
                            end do
                            wt(n - 5) = sum(fdcoef_1f50*w(n - 5:n))
                            wt(n - 4) = sum(fdcoef_1f40*w(n - 4:n))
                            wt(n - 3) = sum(fdcoef_1f30*w(n - 3:n))
                            wt(n - 2) = sum(fdcoef_1f20*w(n - 2:n))
                            wt(n - 1) = sum(fdcoef_1f10*w(n - 1:n))
                            wt(n) = w(n) - w(n - 1)
                        case (7)
                            do i = 1, n - 7
                                wt(i) = sum(fdcoef_1f70*w(i:i + 7))
                            end do
                            wt(n - 6) = sum(fdcoef_1f60*w(n - 6:n))
                            wt(n - 5) = sum(fdcoef_1f50*w(n - 5:n))
                            wt(n - 4) = sum(fdcoef_1f40*w(n - 4:n))
                            wt(n - 3) = sum(fdcoef_1f30*w(n - 3:n))
                            wt(n - 2) = sum(fdcoef_1f20*w(n - 2:n))
                            wt(n - 1) = sum(fdcoef_1f10*w(n - 1:n))
                            wt(n) = w(n) - w(n - 1)
                        case (8)
                            do i = 1, n - 8
                                wt(i) = sum(fdcoef_1f80*w(i:i + 8))
                            end do
                            wt(n - 7) = sum(fdcoef_1f70*w(n - 7:n))
                            wt(n - 6) = sum(fdcoef_1f60*w(n - 6:n))
                            wt(n - 5) = sum(fdcoef_1f50*w(n - 5:n))
                            wt(n - 4) = sum(fdcoef_1f40*w(n - 4:n))
                            wt(n - 3) = sum(fdcoef_1f30*w(n - 3:n))
                            wt(n - 2) = sum(fdcoef_1f20*w(n - 2:n))
                            wt(n - 1) = sum(fdcoef_1f10*w(n - 1:n))
                            wt(n) = w(n) - w(n - 1)
                    end select

                case ('center')
                    select case (deriv_accuracy)
                        case (1)
                            do i = 2, n - 1
                                wt(i) = sum(fdcoef_1c10(1:3)*w(i - 1:i + 1))
                            end do
                            wt(1) = w(2) - w(1)
                            wt(n) = w(n) - w(n - 1)
                        case (2)
                            do i = 3, n - 2
                                wt(i) = sum(fdcoef_1c20(1:5)*w(i - 2:i + 2))
                            end do
                            wt(2) = sum(fdcoef_1c10(1:3)*w(1:3))
                            wt(n - 1) = sum(-fdcoef_1c10(3:1:-1)*w(n - 2:n))
                            wt(1) = w(2) - w(1)
                            wt(n) = w(n) - w(n - 1)
                        case (3)
                            do i = 4, n - 3
                                wt(i) = sum(fdcoef_1c30(1:7)*w(i - 3:i + 3))
                            end do
                            wt(3) = sum(fdcoef_1c20(1:5)*w(1:5))
                            wt(n - 2) = sum(-fdcoef_1c20(5:1:-1)*w(n - 4:n))
                            wt(2) = sum(fdcoef_1c10(1:3)*w(1:3))
                            wt(n - 1) = sum(-fdcoef_1c10(3:1:-1)*w(n - 2:n))
                            wt(1) = w(2) - w(1)
                            wt(n) = w(n) - w(n - 1)
                        case (4)
                            do i = 5, n - 4
                                wt(i) = sum(fdcoef_1c40(1:9)*w(i - 4:i + 4))
                            end do
                            wt(4) = sum(fdcoef_1c30(1:7)*w(1:7))
                            wt(n - 3) = sum(-fdcoef_1c30(7:1:-1)*w(n - 6:n))
                            wt(3) = sum(fdcoef_1c20(1:5)*w(1:5))
                            wt(n - 2) = sum(-fdcoef_1c20(5:1:-1)*w(n - 4:n))
                            wt(2) = sum(fdcoef_1c10(1:3)*w(1:3))
                            wt(n - 1) = sum(-fdcoef_1c10(3:1:-1)*w(n - 2:n))
                            wt(1) = w(2) - w(1)
                            wt(n) = w(n) - w(n - 1)
                        case (5)
                            do i = 6, n - 5
                                wt(i) = sum(fdcoef_1c50(1:11)*w(i - 5:i + 5))
                            end do
                            wt(5) = sum(fdcoef_1c40(1:9)*w(1:9))
                            wt(n - 4) = sum(-fdcoef_1c40(9:1:-1)*w(n - 8:n))
                            wt(4) = sum(fdcoef_1c30(1:7)*w(1:7))
                            wt(n - 3) = sum(-fdcoef_1c30(7:1:-1)*w(n - 6:n))
                            wt(3) = sum(fdcoef_1c20(1:5)*w(1:5))
                            wt(n - 2) = sum(-fdcoef_1c20(5:1:-1)*w(n - 4:n))
                            wt(2) = sum(fdcoef_1c10(1:3)*w(1:3))
                            wt(n - 1) = sum(-fdcoef_1c10(3:1:-1)*w(n - 2:n))
                            wt(1) = w(2) - w(1)
                            wt(n) = w(n) - w(n - 1)
                        case (6)
                            do i = 7, n - 6
                                wt(i) = sum(fdcoef_1c60(1:13)*w(i - 6:i + 6))
                            end do
                            wt(6) = sum(fdcoef_1c50(1:11)*w(1:11))
                            wt(n - 5) = sum(-fdcoef_1c50(11:1:-1)*w(n - 10:n))
                            wt(5) = sum(fdcoef_1c40(1:9)*w(1:9))
                            wt(n - 4) = sum(-fdcoef_1c40(9:1:-1)*w(n - 8:n))
                            wt(4) = sum(fdcoef_1c30(1:7)*w(1:7))
                            wt(n - 3) = sum(-fdcoef_1c30(7:1:-1)*w(n - 6:n))
                            wt(3) = sum(fdcoef_1c20(1:5)*w(1:5))
                            wt(n - 2) = sum(-fdcoef_1c20(5:1:-1)*w(n - 4:n))
                            wt(2) = sum(fdcoef_1c10(1:3)*w(1:3))
                            wt(n - 1) = sum(-fdcoef_1c10(3:1:-1)*w(n - 2:n))
                            wt(1) = w(2) - w(1)
                            wt(n) = w(n) - w(n - 1)
                        case (7)
                            do i = 8, n - 7
                                wt(i) = sum(fdcoef_1c70(1:15)*w(i - 7:i + 7))
                            end do
                            wt(7) = sum(fdcoef_1c60(1:13)*w(1:13))
                            wt(n - 6) = sum(-fdcoef_1c60(13:1:-1)*w(n - 12:n))
                            wt(6) = sum(fdcoef_1c50(1:11)*w(1:11))
                            wt(n - 5) = sum(-fdcoef_1c50(11:1:-1)*w(n - 10:n))
                            wt(5) = sum(fdcoef_1c40(1:9)*w(1:9))
                            wt(n - 4) = sum(-fdcoef_1c40(9:1:-1)*w(n - 8:n))
                            wt(4) = sum(fdcoef_1c30(1:7)*w(1:7))
                            wt(n - 3) = sum(-fdcoef_1c30(7:1:-1)*w(n - 6:n))
                            wt(3) = sum(fdcoef_1c20(1:5)*w(1:5))
                            wt(n - 2) = sum(-fdcoef_1c20(5:1:-1)*w(n - 4:n))
                            wt(2) = sum(fdcoef_1c10(1:3)*w(1:3))
                            wt(n - 1) = sum(-fdcoef_1c10(3:1:-1)*w(n - 2:n))
                            wt(1) = w(2) - w(1)
                            wt(n) = w(n) - w(n - 1)
                        case (8)
                            do i = 9, n - 8
                                wt(i) = sum(fdcoef_1c80(1:17)*w(i - 8:i + 8))
                            end do
                            wt(8) = sum(fdcoef_1c70(1:15)*w(1:15))
                            wt(n - 7) = sum(-fdcoef_1c70(15:1:-1)*w(n - 14:n))
                            wt(7) = sum(fdcoef_1c60(1:13)*w(1:13))
                            wt(n - 6) = sum(-fdcoef_1c60(13:1:-1)*w(n - 12:n))
                            wt(6) = sum(fdcoef_1c50(1:11)*w(1:11))
                            wt(n - 5) = sum(-fdcoef_1c50(11:1:-1)*w(n - 10:n))
                            wt(5) = sum(fdcoef_1c40(1:9)*w(1:9))
                            wt(n - 4) = sum(-fdcoef_1c40(9:1:-1)*w(n - 8:n))
                            wt(4) = sum(fdcoef_1c30(1:7)*w(1:7))
                            wt(n - 3) = sum(-fdcoef_1c30(7:1:-1)*w(n - 6:n))
                            wt(3) = sum(fdcoef_1c20(1:5)*w(1:5))
                            wt(n - 2) = sum(-fdcoef_1c20(5:1:-1)*w(n - 4:n))
                            wt(2) = sum(fdcoef_1c10(1:3)*w(1:3))
                            wt(n - 1) = sum(-fdcoef_1c10(3:1:-1)*w(n - 2:n))
                            wt(1) = w(2) - w(1)
                            wt(n) = w(n) - w(n - 1)
                    end select

                case ('backward')
                    select case (deriv_accuracy)
                        case (1)
                            do i = 2, n
                                wt(i) = sum(fdcoef_1b10*w(i - 1:i))
                            end do
                            wt(1) = w(2) - w(1)
                        case (2)
                            do i = 3, n
                                wt(i) = sum(fdcoef_1b20*w(i - 2:i))
                            end do
                            wt(2) = sum(fdcoef_1b10*w(1:2))
                            wt(1) = w(2) - w(1)
                        case (3)
                            do i = 4, n
                                wt(i) = sum(fdcoef_1b30*w(i - 3:i))
                            end do
                            wt(3) = sum(fdcoef_1b20*w(1:3))
                            wt(2) = sum(fdcoef_1b10*w(1:2))
                            wt(1) = w(2) - w(1)
                        case (4)
                            do i = 5, n
                                wt(i) = sum(fdcoef_1b40*w(i - 4:i))
                            end do
                            wt(4) = sum(fdcoef_1b30*w(1:4))
                            wt(3) = sum(fdcoef_1b20*w(1:3))
                            wt(2) = sum(fdcoef_1b10*w(1:2))
                            wt(1) = w(2) - w(1)
                        case (5)
                            do i = 6, n
                                wt(i) = sum(fdcoef_1b50*w(i - 5:i))
                            end do
                            wt(5) = sum(fdcoef_1b40*w(1:5))
                            wt(4) = sum(fdcoef_1b30*w(1:4))
                            wt(3) = sum(fdcoef_1b20*w(1:3))
                            wt(2) = sum(fdcoef_1b10*w(1:2))
                            wt(1) = w(2) - w(1)
                        case (6)
                            do i = 7, n
                                wt(i) = sum(fdcoef_1b60*w(i - 6:i))
                            end do
                            wt(6) = sum(fdcoef_1b50*w(1:6))
                            wt(5) = sum(fdcoef_1b40*w(1:5))
                            wt(4) = sum(fdcoef_1b30*w(1:4))
                            wt(3) = sum(fdcoef_1b20*w(1:3))
                            wt(2) = sum(fdcoef_1b10*w(1:2))
                            wt(1) = w(2) - w(1)
                        case (7)
                            do i = 8, n
                                wt(i) = sum(fdcoef_1b70*w(i - 7:i))
                            end do
                            wt(7) = sum(fdcoef_1b60*w(1:7))
                            wt(6) = sum(fdcoef_1b50*w(1:6))
                            wt(5) = sum(fdcoef_1b40*w(1:5))
                            wt(4) = sum(fdcoef_1b30*w(1:4))
                            wt(3) = sum(fdcoef_1b20*w(1:3))
                            wt(2) = sum(fdcoef_1b10*w(1:2))
                            wt(1) = w(2) - w(1)
                        case (8)
                            do i = 9, n
                                wt(i) = sum(fdcoef_1b80*w(i - 8:i))
                            end do
                            wt(8) = sum(fdcoef_1b70*w(1:8))
                            wt(7) = sum(fdcoef_1b60*w(1:7))
                            wt(6) = sum(fdcoef_1b50*w(1:6))
                            wt(5) = sum(fdcoef_1b40*w(1:5))
                            wt(4) = sum(fdcoef_1b30*w(1:4))
                            wt(3) = sum(fdcoef_1b20*w(1:3))
                            wt(2) = sum(fdcoef_1b10*w(1:2))
                            wt(1) = w(2) - w(1)
                    end select

                case ('staggered-center')
                    select case (deriv_accuracy)
                        case (1)
                            do i = 1, n - 1
                                wt(i) = sum(fdcoef_1s10*w(i:i + 1))
                            end do
                        case (2)
                            do i = 2, n - 2
                                wt(i) = sum(fdcoef_1s20*w(i - 1:i + 2))
                            end do
                            wt(1) = sum(fdcoef_1s10*w(1:2))
                            wt(n - 1) = sum(fdcoef_1s10*w(n - 1:n))
                        case (3)
                            do i = 3, n - 3
                                wt(i) = sum(fdcoef_1s30*w(i - 2:i + 3))
                            end do
                            wt(1) = sum(fdcoef_1s10*w(1:2))
                            wt(n - 1) = sum(fdcoef_1s10*w(n - 1:n))
                            wt(2) = sum(fdcoef_1s20*w(1:4))
                            wt(n - 2) = sum(fdcoef_1s20*w(n - 3:n))
                        case (4)
                            do i = 4, n - 4
                                wt(i) = sum(fdcoef_1s40*w(i - 3:i + 4))
                            end do
                            wt(1) = sum(fdcoef_1s10*w(1:2))
                            wt(n - 1) = sum(fdcoef_1s10*w(n - 1:n))
                            wt(2) = sum(fdcoef_1s20*w(1:4))
                            wt(n - 2) = sum(fdcoef_1s20*w(n - 3:n))
                            wt(3) = sum(fdcoef_1s30*w(1:6))
                            wt(n - 3) = sum(fdcoef_1s30*w(n - 5:n))

                    end select

            end select

        case (2)

            select case (deriv_method)

                case ('forward')
                    do i = 1, n - 2
                        wt(i) = w(i + 2) - 2*w(i + 1) + w(i)
                    end do
                    wt(n - 1) = wt(n - 2) + (wt(n - 2) - wt(n - 3))
                    wt(n) = wt(n - 1) + (wt(n - 1) - wt(n - 2))

                case ('center')
                    select case (deriv_accuracy)
                        case (1)
                            do i = 2, n - 1
                                wt(i) = sum(fdcoef_2c10*w(i - 1:i + 1))
                            end do
                            wt(1) = w(1) - 2*w(2) + w(3)
                            wt(n) = w(n) - 2*w(n - 1) + w(n - 2)
                        case (2)
                            do i = 3, n - 2
                                wt(i) = sum(fdcoef_2c20*w(i - 2:i + 2))
                            end do
                            wt(2) = sum(fdcoef_2c10*w(1:3))
                            wt(n - 1) = sum(fdcoef_2c10*w(n - 2:n))
                            wt(1) = w(1) - 2*w(2) + w(3)
                            wt(n) = w(n) - 2*w(n - 1) + w(n - 2)
                        case (3)
                            do i = 4, n - 3
                                wt(i) = sum(fdcoef_2c30*w(i - 3:i + 3))
                            end do
                            wt(3) = sum(fdcoef_2c20*w(1:5))
                            wt(n - 2) = sum(fdcoef_2c20*w(n - 4:n))
                            wt(2) = sum(fdcoef_2c10*w(1:3))
                            wt(n - 1) = sum(fdcoef_2c10*w(n - 2:n))
                            wt(1) = w(1) - 2*w(2) + w(3)
                            wt(n) = w(n) - 2*w(n - 1) + w(n - 2)
                        case (4)
                            do i = 5, n - 4
                                wt(i) = sum(fdcoef_2c40*w(i - 4:i + 4))
                            end do
                            wt(4) = sum(fdcoef_2c30*w(1:7))
                            wt(n - 3) = sum(fdcoef_2c30*w(n - 6:n))
                            wt(3) = sum(fdcoef_2c20*w(1:5))
                            wt(n - 2) = sum(fdcoef_2c20*w(n - 4:n))
                            wt(2) = sum(fdcoef_2c10*w(1:3))
                            wt(n - 1) = sum(fdcoef_2c10*w(n - 2:n))
                            wt(1) = w(1) - 2*w(2) + w(3)
                            wt(n) = w(n) - 2*w(n - 1) + w(n - 2)
                        case (5)
                            do i = 6, n - 5
                                wt(i) = sum(fdcoef_2c50*w(i - 5:i + 5))
                            end do
                            wt(5) = sum(fdcoef_2c40*w(1:9))
                            wt(n - 4) = sum(fdcoef_2c40*w(n - 8:n))
                            wt(4) = sum(fdcoef_2c30*w(1:7))
                            wt(n - 3) = sum(fdcoef_2c30*w(n - 6:n))
                            wt(3) = sum(fdcoef_2c20*w(1:5))
                            wt(n - 2) = sum(fdcoef_2c20*w(n - 4:n))
                            wt(2) = sum(fdcoef_2c10*w(1:3))
                            wt(n - 1) = sum(fdcoef_2c10*w(n - 2:n))
                            wt(1) = w(1) - 2*w(2) + w(3)
                            wt(n) = w(n) - 2*w(n - 1) + w(n - 2)
                        case (6)
                            do i = 7, n - 6
                                wt(i) = sum(fdcoef_2c60*w(i - 6:i + 6))
                            end do
                            wt(6) = sum(fdcoef_2c50*w(1:11))
                            wt(n - 5) = sum(fdcoef_2c50*w(n - 10:n))
                            wt(5) = sum(fdcoef_2c40*w(1:9))
                            wt(n - 4) = sum(fdcoef_2c40*w(n - 8:n))
                            wt(4) = sum(fdcoef_2c30*w(1:7))
                            wt(n - 3) = sum(fdcoef_2c30*w(n - 6:n))
                            wt(3) = sum(fdcoef_2c20*w(1:5))
                            wt(n - 2) = sum(fdcoef_2c20*w(n - 4:n))
                            wt(2) = sum(fdcoef_2c10*w(1:3))
                            wt(n - 1) = sum(fdcoef_2c10*w(n - 2:n))
                            wt(1) = w(1) - 2*w(2) + w(3)
                            wt(n) = w(n) - 2*w(n - 1) + w(n - 2)
                        case (7)
                            do i = 8, n - 7
                                wt(i) = sum(fdcoef_2c70*w(i - 7:i + 7))
                            end do
                            wt(7) = sum(fdcoef_2c60*w(1:13))
                            wt(n - 6) = sum(fdcoef_2c60*w(n - 12:n))
                            wt(6) = sum(fdcoef_2c50*w(1:11))
                            wt(n - 5) = sum(fdcoef_2c50*w(n - 10:n))
                            wt(5) = sum(fdcoef_2c40*w(1:9))
                            wt(n - 4) = sum(fdcoef_2c40*w(n - 8:n))
                            wt(4) = sum(fdcoef_2c30*w(1:7))
                            wt(n - 3) = sum(fdcoef_2c30*w(n - 6:n))
                            wt(3) = sum(fdcoef_2c20*w(1:5))
                            wt(n - 2) = sum(fdcoef_2c20*w(n - 4:n))
                            wt(2) = sum(fdcoef_2c10*w(1:3))
                            wt(n - 1) = sum(fdcoef_2c10*w(n - 2:n))
                            wt(1) = w(1) - 2*w(2) + w(3)
                            wt(n) = w(n) - 2*w(n - 1) + w(n - 2)
                        case (8)
                            do i = 9, n - 8
                                wt(i) = sum(fdcoef_2c80*w(i - 8:i + 8))
                            end do
                            wt(8) = sum(fdcoef_2c70*w(1:15))
                            wt(n - 7) = sum(fdcoef_2c70*w(n - 14:n))
                            wt(7) = sum(fdcoef_2c60*w(1:13))
                            wt(n - 6) = sum(fdcoef_2c60*w(n - 12:n))
                            wt(6) = sum(fdcoef_2c50*w(1:11))
                            wt(n - 5) = sum(fdcoef_2c50*w(n - 10:n))
                            wt(5) = sum(fdcoef_2c40*w(1:9))
                            wt(n - 4) = sum(fdcoef_2c40*w(n - 8:n))
                            wt(4) = sum(fdcoef_2c30*w(1:7))
                            wt(n - 3) = sum(fdcoef_2c30*w(n - 6:n))
                            wt(3) = sum(fdcoef_2c20*w(1:5))
                            wt(n - 2) = sum(fdcoef_2c20*w(n - 4:n))
                            wt(2) = sum(fdcoef_2c10*w(1:3))
                            wt(n - 1) = sum(fdcoef_2c10*w(n - 2:n))
                            wt(1) = w(1) - 2*w(2) + w(3)
                            wt(n) = w(n) - 2*w(n - 1) + w(n - 2)
                    end select

                case ('backward')
                    do i = 3, n
                        wt(i) = w(i - 2) - 2*w(i - 1) + w(i)
                    end do
                    wt(1) = wt(2) - (wt(3) - wt(2))
                    wt(2) = wt(3) - (wt(4) - wt(3))

                case ('staggered-center')
                    select case (deriv_accuracy)
                        case (1)
                            do i = 1, n
                            end do
                    end select

            end select

    end select

end function differentiate_1d_

!
!> Differentiate 2D array using finite difference
!
function differentiate_2d_(w, dim, order, method, accuracy) result(wt)

    TT, dimension(:, :), intent(in) :: w
    integer, intent(in), optional :: dim, order, accuracy
    character(len=*), intent(in), optional :: method
    TT, allocatable, dimension(:, :) :: wt

    integer :: along, n1, n2, i, j
    character(len=12) :: deriv_method
    integer :: deriv_order, deriv_accuracy

    if (present(dim)) then
        along = dim
    else
        along = 1
    end if

    if (present(order)) then
        deriv_order = order
    else
        deriv_order = 1
    end if

    if (present(accuracy)) then
        deriv_accuracy = accuracy
    else
        deriv_accuracy = 1
    end if

    if (present(method)) then
        deriv_method = trim(adjustl((method)))
    else
        deriv_method = 'center'
    end if

    n1 = size(w, 1)
    n2 = size(w, 2)
    wt = w

    select case (along)
        case (1)
            !$omp parallel do private(j)
            do j = 1, n2
                wt(:, j) = differentiate_1d_(wt(:, j), deriv_order, deriv_method, deriv_accuracy)
            end do
            !$omp end parallel do
        case (2)
            !$omp parallel do private(i)
            do i = 1, n1
                wt(i, :) = differentiate_1d_(wt(i, :), deriv_order, deriv_method, deriv_accuracy)
            end do
            !$omp end parallel do
    end select

end function differentiate_2d_

!
!> Differentiate 3D array using finite difference
!
function differentiate_3d_(w, dim, order, method, accuracy) result(wt)

    TT, dimension(:, :, :), intent(in) :: w
    integer, intent(in), optional :: dim, order, accuracy
    character(len=*), intent(in), optional :: method
    TT, allocatable, dimension(:, :, :) :: wt

    integer :: along, n1, n2, n3, i, j, k
    character(len=12) :: deriv_method
    integer :: deriv_order, deriv_accuracy

    if (present(dim)) then
        along = dim
    else
        along = 1
    end if

    if (present(order)) then
        deriv_order = order
    else
        deriv_order = 1
    end if

    if (present(accuracy)) then
        deriv_accuracy = accuracy
    else
        deriv_accuracy = 1
    end if

    if (present(method)) then
        deriv_method = trim(adjustl((method)))
    else
        deriv_method = 'center'
    end if

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)
    wt = w

    select case (along)
        case (1)
            !$omp parallel do private(j, k)
            do k = 1, n3
                do j = 1, n2
                    wt(:, j, k) = differentiate_1d_(wt(:, j, k), deriv_order, deriv_method, deriv_accuracy)
                end do
            end do
            !$omp end parallel do
        case (2)
            !$omp parallel do private(i, k)
            do k = 1, n3
                do i = 1, n1
                    wt(i, :, k) = differentiate_1d_(wt(i, :, k), deriv_order, deriv_method, deriv_accuracy)
                end do
            end do
            !$omp end parallel do
        case (3)
            !$omp parallel do private(i, j)
            do j = 1, n2
                do i = 1, n1
                    wt(i, j, :) = differentiate_1d_(wt(i, j, :), deriv_order, deriv_method, deriv_accuracy)
                end do
            end do
            !$omp end parallel do
    end select

end function differentiate_3d_

!
!> Integrate 1D array using cumulative trapezoidal rule
!
function integrate_1d_(w) result(wt)

    TT, dimension(:) :: w
    TT, allocatable, dimension(:) :: wt

    integer :: i, n

    n = size(w)
    allocate (wt(1:n))

    wt(1) = 0
    do i = 2, n
        wt(i) = wt(i - 1) + 0.5d0*(w(i - 1) + w(i))
    end do

end function integrate_1d_

!
!> Integrate 2D array using cumulative trapezoidal rule
!
function integrate_2d_(w, dim) result(wt)

    TT, dimension(:, :) :: w
    integer, optional :: dim
    TT, allocatable, dimension(:, :) :: wt

    integer :: along, n1, n2, i, j

    if (present(dim)) then
        along = dim
    else
        along = 1
    end if

    n1 = size(w, 1)
    n2 = size(w, 2)
    allocate (wt(1:n1, 1:n2))

    select case (along)

        case (1)
            !$omp parallel do private(j)
            do j = 1, n2
                wt(:, j) = integrate_1d_(w(:, j))
            end do
            !$omp end parallel do

        case (2)
            !$omp parallel do private(i)
            do i = 1, n1
                wt(i, :) = integrate_1d_(w(i, :))
            end do
            !$omp end parallel do

    end select

end function integrate_2d_

!
!> Integrate 3D array using cumulative trapezoidal rule
!
function integrate_3d_(w, dim) result(wt)

    TT, dimension(:, :, :) :: w
    integer, optional :: dim
    TT, allocatable, dimension(:, :, :) :: wt

    integer :: along, n1, n2, n3, i, j, k

    if (present(dim)) then
        along = dim
    else
        along = 1
    end if

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)
    allocate (wt(1:n1, 1:n2, 1:n3))

    select case (along)

        case (1)
            !$omp parallel do private(j, k)
            do k = 1, n3
                do j = 1, n2
                    wt(:, j, k) = integrate_1d_(w(:, j, k))
                end do
            end do
            !$omp end parallel do

        case (2)
            !$omp parallel do private(i, k)
            do k = 1, n3
                do i = 1, n1
                    wt(i, :, k) = integrate_1d_(w(i, :, k))
                end do
            end do
            !$omp end parallel do

        case (3)
            !$omp parallel do private(i, j)
            do j = 1, n2
                do i = 1, n1
                    wt(i, j, :) = integrate_1d_(w(i, j, :))
                end do
            end do
            !$omp end parallel do

    end select

end function integrate_3d_

!
!> Integrate 1D array using cumulative summation
!
function cumsum_1d_(w) result(wt)

    TT, dimension(:) :: w
    TT, allocatable, dimension(:) :: wt

    integer :: i, n

    n = size(w)
    allocate (wt(1:n))

    wt(1) = w(1)
    do i = 2, n
        wt(i) = wt(i - 1) + w(i)
    end do

end function cumsum_1d_

!
!> Integrate 2D array using cumulative summation
!
function cumsum_2d_(w, dim) result(wt)

    TT, dimension(:, :) :: w
    integer, optional :: dim
    TT, allocatable, dimension(:, :) :: wt

    integer :: along, n1, n2, i, j

    if (present(dim)) then
        along = dim
    else
        along = 1
    end if

    n1 = size(w, 1)
    n2 = size(w, 2)
    allocate (wt(1:n1, 1:n2))

    select case (along)
        case (1)
            wt(1, :) = w(1, :)
            !$omp parallel do private(i)
            do i = 2, n1
                wt(i, :) = wt(i - 1, :) + w(i, :)
            end do
            !$omp end parallel do
        case (2)
            wt(:, 1) = w(:, 1)
            !$omp parallel do private(j)
            do j = 2, n2
                wt(:, j) = wt(:, j - 1) + w(:, j)
            end do
            !$omp end parallel do
    end select

end function cumsum_2d_

!
!> Integrate 3D array using cumulative summation
!
function cumsum_3d_(w, dim) result(wt)

    TT, dimension(:, :, :) :: w
    integer, optional :: dim
    TT, allocatable, dimension(:, :, :) :: wt

    integer :: along, n1, n2, n3, i, j, k

    if (present(dim)) then
        along = dim
    else
        along = 1
    end if

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)
    allocate (wt(1:n1, 1:n2, 1:n3))

    select case (along)
        case (1)
            wt(1, :, :) = w(1, :, :)
            !$omp parallel do private(i)
            do i = 2, n1
                wt(i, :, :) = wt(i - 1, :, :) + w(i, :, :)
            end do
            !$omp end parallel do
        case (2)
            wt(:, 1, :) = w(:, 1, :)
            !$omp parallel do private(j)
            do j = 2, n2
                wt(:, j, :) = wt(:, j - 1, :) + w(:, j, :)
            end do
            !$omp end parallel do
        case (3)
            wt(:, :, 1) = w(:, :, 1)
            !$omp parallel do private(k)
            do k = 2, n3
                wt(:, :, k) = wt(:, :, k - 1) + w(:, :, k)
            end do
            !$omp end parallel do
    end select

end function cumsum_3d_

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef differentiate_1d_
#undef differentiate_2d_
#undef differentiate_3d_
#undef integrate_1d_
#undef integrate_2d_
#undef integrate_3d_
#undef cumsum_1d_
#undef cumsum_2d_
#undef cumsum_3d_
