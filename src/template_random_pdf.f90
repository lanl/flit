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

#define random_pdf_1d_      CONCAT(random_pdf_1d, T)
#define random_pdf_2d_      CONCAT(random_pdf_2d, T)
#define random_pdf_3d_      CONCAT(random_pdf_3d, T)

function random_pdf_1d_(n, pdf, seed, range) result(y)

    integer, intent(in) :: n
    TT, dimension(:), intent(in) :: pdf
    integer, intent(in), optional :: seed
    TT, dimension(:), intent(in), optional :: range
    TT, allocatable, dimension(:) :: y

    TT, allocatable, dimension(:) :: r, x, xx, yy, cdf
    integer :: np

    if (present(range)) then
        r = range
    else
        r = [0.0, 1.0]
    end if

    np = size(pdf)

    x = linspace(r(1), r(2), np)
    xx = linspace(r(1), r(2), max(np, 1000))
    yy = ginterp(x, pdf, xx, 'linear')

    cdf = integ(yy)
    cdf = cdf/maxval(cdf)

    if (present(seed)) then
        y = _random_(n, seed=seed)
    else
        y = _random_(n)
    end if

    y = ginterp(cdf, xx, y, 'linear')

end function

function random_pdf_2d_(n1, n2, pdf, seed, range) result(y)

    integer, intent(in) :: n1, n2
    TT, dimension(:), intent(in) :: pdf
    integer, intent(in), optional :: seed
    TT, dimension(:), intent(in), optional :: range
    TT, allocatable, dimension(:, :) :: y

    TT, allocatable, dimension(:) :: r, x, xx, yy, cdf
    integer :: np

    if (present(range)) then
        r = range
    else
        r = [0.0, 1.0]
    end if

    np = size(pdf)

    x = linspace(r(1), r(2), np)
    xx = linspace(r(1), r(2), max(np, 1000))
    yy = ginterp(x, pdf, xx, 'linear')

    cdf = integ(yy)
    cdf = cdf/maxval(cdf)

    if (present(seed)) then
        y = _random_(n1, n2, seed=seed)
    else
        y = _random_(n1, n2)
    end if

    y = reshape(ginterp(cdf, xx, flatten(y), 'linear'), shape(y))

end function

function random_pdf_3d_(n1, n2, n3, pdf, seed, range) result(y)

    integer, intent(in) :: n1, n2, n3
    TT, dimension(:), intent(in) :: pdf
    integer, intent(in), optional :: seed
    TT, dimension(:), intent(in), optional :: range
    TT, allocatable, dimension(:, :, :) :: y

    TT, allocatable, dimension(:) :: r, x, xx, yy, cdf
    integer :: np

    if (present(range)) then
        r = range
    else
        r = [0.0, 1.0]
    end if

    np = size(pdf)

    x = linspace(r(1), r(2), np)
    xx = linspace(r(1), r(2), max(np, 1000))
    yy = ginterp(x, pdf, xx, 'linear')

    cdf = integ(yy)
    cdf = cdf/maxval(cdf)

    if (present(seed)) then
        y = _random_(n1, n2, n3, seed=seed)
    else
        y = _random_(n1, n2, n3)
    end if

    y = reshape(ginterp(cdf, xx, flatten(y), 'linear'), shape(y))

end function

#undef T
#undef TT
#undef _random_

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef random_pdf_1d_
#undef random_pdf_2d_
#undef random_pdf_3d_
