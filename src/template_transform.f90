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

#define convfast_1d_     CONCAT(convfast_1d, T)
#define convsum_1d_     CONCAT(convsum_1d, T)
#define discrete_conv_1d_     CONCAT(discrete_conv_1d, T)
#define discrete_conv_2d_     CONCAT(discrete_conv_2d, T)
#define discrete_conv_3d_     CONCAT(discrete_conv_3d, T)
#define conv_1d_     CONCAT(conv_1d, T)
#define conv_2d_     CONCAT(conv_2d, T)
#define conv_3d_     CONCAT(conv_3d, T)
#define fftshift_1d_     CONCAT(fftshift_1d, T)
#define fftshift_2d_     CONCAT(fftshift_2d, T)
#define fftshift_3d_     CONCAT(fftshift_3d, T)
#define deconv_1d_     CONCAT(deconv_1d, T)
#define deconv_2d_     CONCAT(deconv_2d, T)
#define deconv_3d_     CONCAT(deconv_3d, T)

!
!> Fast 1D discrete convolution
!
subroutine convfast_1d_(lx, kx, x, ly, ky, y, lz, kz, z)

    integer, intent(in) :: lx, ly, lz
    integer, intent(in) :: kx, ky, kz
    TT, dimension(0:), intent(in) :: x, y
    TT, dimension(0:), intent(inout) :: z

    ! Bounds for index i.
    integer :: imin
    integer :: imax

    ! Variables that we expect to reside in registers.
    integer :: i, ilo, ihi, j, jlo, jhi, iz
    TT :: sa, sb, xa, xb, ya, yb

    imin = kz - kx - ky
    imax = imin + lz - 1

    ! Off left: imin <= i <= -1
    ilo = imin
    ihi = min(-1, imax)
    !    for (i=ilo,iz=i-imin i<=ihi ++i,++iz)
    do i = ilo, ihi
        iz = i - imin
        z(iz) = 0.0
    end do

    ! Rolling on: 0 <= i <= lx-2 and 0 <= j <= i
    ilo = max(0, imin)
    ihi = min(lx - 2, imax)
    jlo = 0
    jhi = ilo
    !    for (i=ilo,iz=i-imin i<ihi i+=2,iz+=2,jhi+=2) {
    do i = ilo, ihi - 1, 2

        iz = i - imin

        sa = 0.0
        sb = 0.0
        yb = y(i - jlo + 1)
        !      for (j=jlo j<jhi j+=2) {
        do j = jlo, jhi - 1, 2
            xa = x(j)
            sb = sb + xa*yb
            ya = y(i - j)
            sa = sa + xa*ya
            xb = x(j + 1)
            sb = sb + xb*ya
            yb = y(i - j - 1)
            sa = sa + xb*yb
        end do
        xa = x(j)
        sb = sb + xa*yb
        if (j == jhi) then
            ya = y(i - j)
            sa = sa + xa*ya
            xb = x(j + 1)
            sb = sb + xb*ya
        end if
        z(iz) = sa
        z(iz + 1) = sb

        jhi = jhi + 2

    end do

    if (i == ihi) then
        iz = i - imin
        jlo = 0
        jhi = i
        sa = 0.0
        !      for (j=jlo j<=jhi ++j)
        do j = jlo, jhi
            sa = sa + x(j)*y(i - j)
        end do
        z(iz) = sa
    end if

    ! Middle: lx-1 <= i <= ly-1 and 0 <= j <= lx-1
    ilo = max(lx - 1, imin)
    ihi = min(ly - 1, imax)
    jlo = 0
    jhi = lx - 1
    !    for (i=ilo,iz=i-imin i<ihi i+=2,iz+=2) {
    do i = ilo, ihi - 1, 2

        iz = i - imin

        sa = 0.0
        sb = 0.0
        yb = y(i - jlo + 1)
        !      for (j=jlo j<jhi j+=2) {
        do j = jlo, jhi - 1, 2
            xa = x(j)
            sb = sb + xa*yb
            ya = y(i - j)
            sa = sa + xa*ya
            xb = x(j + 1)
            sb = sb + xb*ya
            yb = y(i - j - 1)
            sa = sa + xb*yb
        end do
        if (j == jhi) then
            xa = x(j)
            sb = sb + xa*yb
            ya = y(i - j)
            sa = sa + xa*ya
        end if
        z(iz) = sa
        z(iz + 1) = sb

    end do

    if (i == ihi) then
        iz = i - imin
        sa = 0.0
        !      for (j=jlo j<=jhi ++j)
        do j = jlo, jhi
            sa = sa + x(j)*y(i - j)
        end do
        z(iz) = sa
    end if

    ! Rolling off: ly <= i <= lx+ly-2 and i-ly+1 <= j <= lx-1
    ilo = max(ly, imin)
    ihi = min(lx + ly - 2, imax)
    jlo = ihi - ly + 1
    jhi = lx - 1
    !    for (i=ihi,iz=i-imin i>ilo i-=2,iz-=2,jlo-=2) {
    do i = ihi, ilo + 1, -2

        iz = i - imin

        sa = 0.0
        sb = 0.0
        yb = y(i - jhi - 1)
        !      for (j=jhi j>jlo j-=2) {
        do j = jhi, jlo + 1, -2
            xa = x(j)
            sb = sb + xa*yb
            ya = y(i - j)
            sa = sa + xa*ya
            xb = x(j - 1)
            sb = sb + xb*ya
            yb = y(i - j + 1)
            sa = sa + xb*yb
        end do
        xa = x(j)
        sb = sb + xa*yb
        if (j == jlo) then
            ya = y(i - j)
            sa = sa + xa*ya
            xb = x(j - 1)
            sb = sb + xb*ya
        end if
        z(iz) = sa
        z(iz - 1) = sb

        jlo = jlo - 2

    end do

    if (i == ilo) then
        iz = i - imin
        jlo = i - ly + 1
        jhi = lx - 1
        sa = 0.0
        !      for (j=jhi j>=jlo --j)
        do j = jhi, jlo, -1
            sa = sa + x(j)*y(i - j)
        end do
        z(iz) = sa
    end if

    ! Off right: lx+ly-1 <= i <= imax
    ilo = max(lx + ly - 1, imin)
    ihi = imax
    !    for (i=ilo,iz=i-imin i<=ihi ++i,++iz)
    do i = ilo, ihi
        iz = i - imin
        z(iz) = 0.0
    end do

end subroutine convfast_1d_

!
!> Fast 1D discrete convolution
!
subroutine convsum_1d_(lx, kx, x, ly, ky, y, lz, kz, z)

    integer, intent(in) :: lx, ly, lz
    integer, intent(in) :: kx, ky, kz
    TT, dimension(0:), intent(in) :: x, y
    TT, dimension(0:), intent(inout) :: z

    ! Bounds for index i.
    integer :: imin
    integer :: imax

    ! Variables that we expect to reside in registers.
    integer :: i, ilo, ihi, j, jlo, jhi, iz
    TT :: sa, sb, xa, xb, ya, yb

    imin = kz - kx - ky
    imax = imin + lz - 1

    ! Rolling on: 0 <= i <= lx-2 and 0 <= j <= i
    ilo = max(0, imin)
    ihi = min(lx - 2, imax)
    jlo = 0
    jhi = ilo
    !    for (i=ilo,iz=i-imin i<ihi i+=2,iz+=2,jhi+=2) {
    do i = ilo, ihi - 1, 2

        iz = i - imin

        sa = z(iz)
        sb = z(iz + 1)
        yb = y(i - jlo + 1)
        !      for (j=jlo j<jhi j+=2) {
        do j = jlo, jhi - 1, 2
            xa = x(j)
            sb = sb + xa*yb
            ya = y(i - j)
            sa = sa + xa*ya
            xb = x(j + 1)
            sb = sb + xb*ya
            yb = y(i - j - 1)
            sa = sa + xb*yb
        end do
        xa = x(j)
        sb = sb + xa*yb
        if (j == jhi) then
            ya = y(i - j)
            sa = sa + xa*ya
            xb = x(j + 1)
            sb = sb + xb*ya
        end if
        z(iz) = sa
        z(iz + 1) = sb

        jhi = jhi + 2

    end do

    if (i == ihi) then
        iz = i - imin
        jlo = 0
        jhi = i
        sa = z(iz)
        !      for (j=jlo j<=jhi ++j)
        do j = jlo, jhi
            sa = sa + x(j)*y(i - j)
        end do
        z(iz) = sa
    end if

    ! Middle: lx-1 <= i <= ly-1 and 0 <= j <= lx-1
    ilo = max(lx - 1, imin)
    ihi = min(ly - 1, imax)
    jlo = 0
    jhi = lx - 1
    !    for (i=ilo,iz=i-imin i<ihi i+=2,iz+=2) {
    do i = ilo, ihi - 1, 2

        iz = i - imin

        sa = z(iz)
        sb = z(iz + 1)
        yb = y(i - jlo + 1)
        !      for (j=jlo j<jhi j+=2) {
        do j = jlo, jhi - 1, 2
            xa = x(j)
            sb = sb + xa*yb
            ya = y(i - j)
            sa = sa + xa*ya
            xb = x(j + 1)
            sb = sb + xb*ya
            yb = y(i - j - 1)
            sa = sa + xb*yb
        end do
        if (j == jhi) then
            xa = x(j)
            sb = sb + xa*yb
            ya = y(i - j)
            sa = sa + xa*ya
        end if
        z(iz) = sa
        z(iz + 1) = sb

    end do

    if (i == ihi) then
        iz = i - imin
        sa = z(iz)
        !      for (j=jlo j<=jhi ++j)
        do j = jlo, jhi
            sa = sa + x(j)*y(i - j)
        end do
        z(iz) = sa
    end if

    ! Rolling off: ly <= i <= lx+ly-2 and i-ly+1 <= j <= lx-1
    ilo = max(ly, imin)
    ihi = min(lx + ly - 2, imax)
    jlo = ihi - ly + 1
    jhi = lx - 1
    !    for (i=ihi,iz=i-imin i>ilo i-=2,iz-=2,jlo-=2) {
    do i = ihi, ilo + 1, -2

        iz = i - imin

        sa = z(iz)
        sb = z(iz - 1)
        yb = y(i - jhi - 1)
        !      for (j=jhi j>jlo j-=2) {
        do j = jhi, jlo + 1, -2
            xa = x(j)
            sb = sb + xa*yb
            ya = y(i - j)
            sa = sa + xa*ya
            xb = x(j - 1)
            sb = sb + xb*ya
            yb = y(i - j + 1)
            sa = sa + xb*yb
        end do
        xa = x(j)
        sb = sb + xa*yb
        if (j == jlo) then
            ya = y(i - j)
            sa = sa + xa*ya
            xb = x(j - 1)
            sb = sb + xb*ya
        end if
        z(iz) = sa
        z(iz - 1) = sb

        jlo = jlo - 2

    end do

    if (i == ilo) then
        iz = i - imin
        jlo = i - ly + 1
        jhi = lx - 1
        sa = z(iz)
        !      for (j=jhi j>=jlo --j)
        do j = jhi, jlo, -1
            sa = sa + x(j)*y(i - j)
        end do
        z(iz) = sa
    end if

end subroutine convsum_1d_

!
!> Discrete 1D convolution translated from convFast in minesjtk
!
subroutine discrete_conv_1d_(x, y, z, kx, ky, kz)

    TT, dimension(:), intent(in) :: x, y
    TT, dimension(:), intent(inout) :: z
    integer, intent(in) :: kx, ky, kz

    integer :: lx, ly, lz

    lx = size(x, 1)
    ly = size(y, 1)
    lz = size(z, 1)

    z = 0.0

    if (lx > ly) then
        call convfast_1d_( &
            ly, ky - 1, y, &
            lx, kx - 1, x, &
            lz, kz - 1, z)
    else
        call convfast_1d_( &
            lx, kx - 1, x, &
            ly, ky - 1, y, &
            lz, kz - 1, z)
    end if

end subroutine discrete_conv_1d_

!
!> Discrete 2D convolution translated from convFast in minesjtk
!
subroutine discrete_conv_2d_(x, y, z, &
        kx1, kx2, &
        ky1, ky2, &
        kz1, kz2)

    TT, dimension(:, :), intent(in) :: x, y
    TT, dimension(:, :), intent(inout) :: z
    integer, intent(in) :: kx1, ky1, kz1
    integer, intent(in) :: kx2, ky2, kz2

    integer :: lx1, ly1, lz1
    integer :: lx2, ly2, lz2
    integer :: ilo2, ihi2, jlo2, jhi2
    integer :: i2, j2

    lx1 = size(x, 1)
    lx2 = size(x, 2)
    ly1 = size(y, 1)
    ly2 = size(y, 2)
    lz1 = size(z, 1)
    lz2 = size(z, 2)

    z = 0.0

    ilo2 = (kz2 - 1) - (kx2 - 1) - (ky2 - 1)
    ihi2 = ilo2 + lz2 - 1

    if (lx1 > ly1) then
        !            !$omp parallel do private(i2,j2,jlo2,jhi2)
        do i2 = ilo2, ihi2
            jlo2 = max(0, i2 - ly2 + 1)
            jhi2 = min(lx2 - 1, i2)
            do j2 = jlo2, jhi2
                call convsum_1d_( &
                    ly1, ky1 - 1, y(:, i2 - j2 + 1), &
                    lx1, kx1 - 1, x(:, j2 + 1), &
                    lz1, kz1 - 1, z(:, i2 - ilo2 + 1))
            end do
        end do
        !            !$omp end parallel do
    else
        !            !$omp parallel do private(i2,j2,jlo2,jhi2)
        do i2 = ilo2, ihi2
            jlo2 = max(0, i2 - ly2 + 1)
            jhi2 = min(lx2 - 1, i2)
            do j2 = jlo2, jhi2
                call convsum_1d_( &
                    lx1, kx1 - 1, x(:, j2 + 1), &
                    ly1, ky1 - 1, y(:, i2 - j2 + 1), &
                    lz1, kz1 - 1, z(:, i2 - ilo2 + 1))
            end do
        end do
        !            !$omp end parallel do
    end if

end subroutine discrete_conv_2d_

!
!> Discrete 3D convolution
!
subroutine discrete_conv_3d_(x, y, z, &
        kx1, kx2, kx3, &
        ky1, ky2, ky3, &
        kz1, kz2, kz3)

    integer, intent(in) :: kx1, ky1, kz1
    integer, intent(in) :: kx2, ky2, kz2
    integer, intent(in) :: kx3, ky3, kz3
    TT, dimension(:, :, :), intent(in) :: x, y
    TT, dimension(:, :, :), intent(inout) :: z

    integer :: lx1, ly1, lz1
    integer :: lx2, ly2, lz2
    integer :: lx3, ly3, lz3
    integer :: ilo2, ihi2, ilo3, ihi3
    integer :: jlo2, jhi2, jlo3, jhi3
    integer :: i2, j2, i3, j3

    lx1 = size(x, 1)
    lx2 = size(x, 2)
    lx3 = size(x, 3)
    ly1 = size(y, 1)
    ly2 = size(y, 2)
    ly3 = size(y, 3)
    lz1 = size(z, 1)
    lz2 = size(z, 2)
    lz3 = size(z, 3)

    z = 0.0

    ilo2 = (kz2 - 1) - (kx2 - 1) - (ky2 - 1)
    ilo3 = (kz3 - 1) - (kx3 - 1) - (ky3 - 1)
    ihi2 = ilo2 + lz2 - 1
    ihi3 = ilo3 + lz3 - 1

    if (lx1 > ly1) then
        !            !$omp parallel do private(i2,i3,j2,j3,jlo2,jhi2,jlo3,jhi3)
        do i3 = ilo3, ihi3
            jlo3 = max(0, i3 - ly3 + 1)
            jhi3 = min(lx3 - 1, i3)
            do j3 = jlo3, jhi3
                do i2 = ilo2, ihi2
                    jlo2 = max(0, i2 - ly2 + 1)
                    jhi2 = min(lx2 - 1, i2)
                    do j2 = jlo2, jhi2
                        call convsum_1d_( &
                            ly1, ky1 - 1, y(:, i2 - j2 + 1, i3 - j3 + 1), &
                            lx1, kx1 - 1, x(:, j2 + 1, j3 + 1), &
                            lz1, kz1 - 1, z(:, i2 - ilo2 + 1, i3 - ilo3 + 1))
                    end do
                end do
            end do
        end do
        !            !$omp end parallel do
    else
        !            !$omp parallel do private(i2,i3,j2,j3,jlo2,jhi2,jlo3,jhi3)
        do i3 = ilo3, ihi3
            jlo3 = max(0, i3 - ly3 + 1)
            jhi3 = min(lx3 - 1, i3)
            do j3 = jlo3, jhi3
                do i2 = ilo2, ihi2
                    jlo2 = max(0, i2 - ly2 + 1)
                    jhi2 = min(lx2 - 1, i2)
                    do j2 = jlo2, jhi2
                        call convsum_1d_( &
                            lx1, kx1 - 1, x(:, j2 + 1, j3 + 1), &
                            ly1, ky1 - 1, y(:, i2 - j2 + 1, i3 - j3 + 1), &
                            lz1, kz1 - 1, z(:, i2 - ilo2 + 1, i3 - ilo3 + 1))
                    end do
                end do
            end do
        end do
        !            !$omp end parallel do
    end if

end subroutine discrete_conv_3d_

function conv_1d_(x, y, method) result(z)

    TT, dimension(:), intent(in) :: x, y
    character(len=*), intent(in), optional :: method

    TT, allocatable, dimension(:) :: z

    integer :: nx, ny, n
    character(len=12) :: conv_method

    if (present(method)) then
        conv_method = method
    else
        conv_method = 'full'
    end if

    nx = size(x)
    ny = size(y)
    n = next_power_235(nx + ny - 1)

    allocate (z(1:n))
    z = TTT(ifft(fft(pad(x, [0, n - nx], ['const', 'const']))*fft(pad(y, [0, n - ny], ['const', 'const']))))
    select case (conv_method)
        case ('full')
            call alloc_array(z, [1, nx + ny - 1], source=z(1:nx + ny - 1))
        case ('same')
            n = nint((ny + 1)/2.0)
            call alloc_array(z, [1, nx], source=z(n:n + nx - 1))
    end select

end function conv_1d_

function conv_2d_(x, y, method) result(z)

    TT, dimension(:, :), intent(in) :: x, y
    character(len=*), intent(in), optional :: method

    TT, allocatable, dimension(:, :) :: z

    integer :: nx1, nx2, ny1, ny2, n1, n2
    character(len=12) :: conv_method

    if (present(method)) then
        conv_method = method
    else
        conv_method = 'full'
    end if

    nx1 = size(x, 1)
    nx2 = size(x, 2)
    ny1 = size(y, 1)
    ny2 = size(y, 2)
    n1 = next_power_235(nx1 + ny1 - 1)
    n2 = next_power_235(nx2 + ny2 - 1)

    allocate (z(1:n1, 1:n2))
    z = TTT(ifft( &
        fft(pad(x, [0, n1 - nx1, 0, n2 - nx2], ['const', 'const', 'const', 'const']))* &
        fft(pad(y, [0, n1 - ny1, 0, n2 - ny2], ['const', 'const', 'const', 'const']))))
    select case (conv_method)
        case ('full')
            call alloc_array(z, [1, nx1 + ny1 - 1, 1, nx2 + ny2 - 1], &
                source=z(1:nx1 + ny1 - 1, 1:nx2 + ny2 - 1))
        case ('same')
            n1 = nint((ny1 + 1)/2.0)
            n2 = nint((ny2 + 1)/2.0)
            call alloc_array(z, [1, nx1, 1, nx2], &
                source=z(n1:n1 + nx1 - 1, n2:n2 + nx2 - 1))
    end select

end function conv_2d_

function conv_3d_(x, y, method) result(z)

    TT, dimension(:, :, :), intent(in) :: x, y
    character(len=*), intent(in), optional :: method

    TT, allocatable, dimension(:, :, :) :: z

    integer :: nx1, nx2, nx3, ny1, ny2, ny3, n1, n2, n3
    character(len=12) :: conv_method

    if (present(method)) then
        conv_method = method
    else
        conv_method = 'full'
    end if

    nx1 = size(x, 1)
    nx2 = size(x, 2)
    nx3 = size(x, 3)
    ny1 = size(y, 1)
    ny2 = size(y, 2)
    ny3 = size(y, 3)
    n1 = next_power_235(nx1 + ny1 - 1)
    n2 = next_power_235(nx2 + ny2 - 1)
    n3 = next_power_235(nx3 + ny3 - 1)

    allocate (z(1:n1, 1:n2, 1:n3))
    z = TTT(ifft( &
        fft(pad(x, [0, n1 - nx1, 0, n2 - nx2, 0, n3 - nx3], &
        ['const', 'const', 'const', 'const', 'const', 'const']))* &
        fft(pad(y, [0, n1 - ny1, 0, n2 - ny2, 0, n3 - ny3], &
        ['const', 'const', 'const', 'const', 'const', 'const']))))
    select case (conv_method)
        case ('full')
            call alloc_array(z, &
                [1, nx1 + ny1 - 1, 1, nx2 + ny2 - 1, 1, nx3 + ny3 - 1], &
                source=z(1:nx1 + ny1 - 1, 1:nx2 + ny2 - 1, 1:nx3 + ny3 - 1))
        case ('same')
            n1 = nint((ny1 + 1)/2.0)
            n2 = nint((ny2 + 1)/2.0)
            n3 = nint((ny3 + 1)/2.0)
            call alloc_array(z, [1, nx1, 1, nx2, 1, nx3], &
                source=z(n1:n1 + nx1 - 1, n2:n2 + nx2 - 1, n3:n3 + nx3 - 1))
    end select

end function conv_3d_

!
!> Shift 1D FFT result so that low-high-high-low becomes high-low-low-high,
!>        or vice versa
!
function fftshift_1d_(w) result(wr)

    TT, dimension(:) :: w

    integer :: n1
    TT, allocatable, dimension(:) :: wr

    n1 = size(w)

    wr = w
    wr = cshift(wr, nint(n1/2.0))

end function fftshift_1d_

!
!> Shift 2D FFT result so that low-high-high-low becomes high-low-low-high,
!>        or vice versa
!
function fftshift_2d_(w) result(wr)

    TT, dimension(:, :) :: w

    integer :: n1, n2
    TT, allocatable, dimension(:, :) :: wr

    n1 = size(w, 1)
    n2 = size(w, 2)

    wr = w
    wr = cshift(wr, nint(n1/2.0), dim=1)
    wr = cshift(wr, nint(n2/2.0), dim=2)

end function fftshift_2d_

!
!> Shift 2D FFT result so that low-high-high-low becomes high-low-low-high,
!>        or vice versa
!
function fftshift_3d_(w) result(wr)

    TT, dimension(:, :, :) :: w

    integer :: n1, n2, n3
    TT, allocatable, dimension(:, :, :) :: wr

    n1 = size(w, 1)
    n2 = size(w, 2)
    n3 = size(w, 3)

    wr = w
    wr = cshift(wr, nint(n1/2.0), dim=1)
    wr = cshift(wr, nint(n2/2.0), dim=2)
    wr = cshift(wr, nint(n3/2.0), dim=3)

end function fftshift_3d_

!
!> Simple deconvolution (a.k.a. inverse filtering) using FFT
!
function deconv_1d_(x, y, maxlag, eps) result(z)

    TT, dimension(:), intent(in) :: x, y
    integer, intent(in), optional :: maxlag
    TT, intent(in), optional :: eps
    TT, allocatable, dimension(:) :: z

    integer :: nx, ny, n, nlag
    TTTT, allocatable, dimension(:) :: a, b

    ! Dimensions
    nx = size(x)
    ny = size(y)
    if (present(maxlag)) then
        nlag = min(maxlag, max(nx, ny))
    else
        nlag = 0
    end if

    ! Pad to next power 235
    n = next_power_235(max(nx, ny) + nlag)

    ! Prepare arrays
    allocate (a(1:n))
    allocate (b(1:n))
    b = fft(pad(y, [0, n - ny], ['const', 'const']))
    a = fft(pad(x, [0, n - nx], ['const', 'const']))*conjg(b)
    b = b*conjg(b)
    if (present(eps)) then
        b = b + eps*maxval(abs(b))
    else
        ! b = b + 1.0e-8*sum(b)/n
        b = b + 1.0e-3*maxval(abs(b))
    end if

    ! Do the deconvolution
    ! Traditionally, deconvolution is u/d where u is the ground truth and
    ! d is the filter. u/d is to recover the groud truth from blurred measurements.
    ! Translated into the variables used here,
    ! a/b is a like a stabilized deconvolution.
    ! If a maximum time lag is given, then the following first option
    ! is more like a normalized cross-correlation, because the results
    ! embraces both positive and negative time lags.
    ! The second option is more of the traditional deconvolution, whose length
    ! is from 1 to nt (the same with the origional signal)
    allocate (z(1:n))
    if (present(maxlag)) then
        z = fftshift(TTT(ifft(a/b)))
        call alloc_array(z, [-nlag, nlag], &
            source=z(nint((n + 1)/2.0) - nlag:nint((n + 1)/2.0) + nlag))
    else
        z = TTT(ifft(a/b))
        call alloc_array(z, [1, nx], source=z(1:nx))
    end if

end function deconv_1d_

!
!> Simple deconvolution (a.k.a. inverse filtering) using FFT
!
function deconv_2d_(x, y, maxlag, eps) result(z)

    TT, dimension(:, :), intent(in) :: x, y
    integer, dimension(:), intent(in), optional :: maxlag
    TT, intent(in), optional :: eps
    TT, allocatable, dimension(:, :) :: z

    integer :: n1x, n1y, n2x, n2y, n1, n2, nlag1, nlag2
    TTTT, allocatable, dimension(:, :) :: a, b

    ! Dimensions
    n1x = size(x, 1)
    n2x = size(x, 2)
    n1y = size(y, 1)
    n2y = size(y, 2)
    if (present(maxlag)) then
        nlag1 = min(maxlag(1), max(n1x, n1y))
        nlag2 = min(maxlag(2), max(n2x, n2y))
    else
        nlag1 = 0
        nlag2 = 0
    end if

    ! Pad to next power 235
    n1 = next_power_235(max(n1x, n1y) + nlag1)
    n2 = next_power_235(max(n2x, n2y) + nlag2)

    ! Prepare arrays
    allocate (a(1:n1, 1:n2))
    allocate (b(1:n1, 1:n2))
    b = fft(pad(y, [0, n1 - n1y, 0, n2 - n2y], ['const', 'const', 'const', 'const']))
    a = fft(pad(x, [0, n1 - n1x, 0, n2 - n2x], ['const', 'const', 'const', 'const']))*conjg(b)
    b = b*conjg(b)
    if (present(eps)) then
        b = b + eps*maxval(abs(b))
    else
        !b = b + 1.0e-8*sum(b)/(n1*n2)
        b = b + 1.0e-3*maxval(abs(b))
    end if

    ! Do the deconvolution
    ! Traditionally, deconvolution is u/d where u is the ground truth and
    ! d is the filter. u/d is to recover the groud truth from blurred measurements.
    ! Translated into the variables used here,
    ! a/b is a like a stabilized deconvolution.
    ! If a maximum time lag is given, then the following first option
    ! is more like a normalized cross-correlation, because the results
    ! embraces both positive and negative time lags.
    ! The second option is more of the traditional deconvolution, whose length
    ! is from 1 to nt (the same with the origional signal)
    allocate (z(1:n1, 1:n2))
    if (present(maxlag)) then
        z = fftshift(TTT(ifft(a/b)))
        call alloc_array(z, [-nlag1, nlag1, -nlag2, nlag2], &
            source=z(nint((n1 + 1)/2.0) - nlag1:nint((n1 + 1)/2.0) + nlag1, &
            nint((n2 + 1)/2.0) - nlag2:nint((n2 + 1)/2.0) + nlag2))
    else
        z = TTT(ifft(a/b))
        call alloc_array(z, [1, n1x, 1, n2x], source=z(1:n1x, 1:n2x))
    end if

end function deconv_2d_

!
!> Simple deconvolution (a.k.a. inverse filtering) using FFT
!
function deconv_3d_(x, y, maxlag, eps) result(z)

    TT, dimension(:, :, :), intent(in) :: x, y
    integer, dimension(:), intent(in), optional :: maxlag
    TT, intent(in), optional :: eps
    TT, allocatable, dimension(:, :, :) :: z

    integer :: n1x, n1y, n2x, n2y, n3x, n3y, n1, n2, n3, nlag1, nlag2, nlag3
    TTTT, allocatable, dimension(:, :, :) :: a, b

    ! Dimensions
    n1x = size(x, 1)
    n2x = size(x, 2)
    n3x = size(x, 3)
    n1y = size(y, 1)
    n2y = size(y, 2)
    n3y = size(y, 3)
    if (present(maxlag)) then
        nlag1 = min(maxlag(1), max(n1x, n1y))
        nlag2 = min(maxlag(2), max(n2x, n2y))
        nlag3 = min(maxlag(3), max(n3x, n3y))
    else
        nlag1 = 0
        nlag2 = 0
        nlag3 = 0
    end if

    ! Pad to next power 235
    n1 = next_power_235(max(n1x, n1y) + nlag1)
    n2 = next_power_235(max(n2x, n2y) + nlag2)
    n3 = next_power_235(max(n3x, n3y) + nlag3)

    ! Prepare arrays
    allocate (a(1:n1, 1:n2, 1:n3))
    allocate (b(1:n1, 1:n2, 1:n3))
    b = fft(pad(y, [0, n1 - n1y, 0, n2 - n2y, 0, n3 - n3y], ['const', 'const', 'const', 'const', 'const', 'const']))
    a = fft(pad(x, [0, n1 - n1x, 0, n2 - n2x, 0, n3 - n3x], ['const', 'const', 'const', 'const', 'const', 'const']))*conjg(b)
    b = b*conjg(b)
    if (present(eps)) then
        b = b + eps*maxval(abs(b))
    else
        !b = b + 1.0e-8*sum(b)/(n1*n2)
        b = b + 1.0e-3*maxval(abs(b))
    end if

    ! Do the deconvolution
    ! Traditionally, deconvolution is u/d where u is the ground truth and
    ! d is the filter. u/d is to recover the groud truth from blurred measurements.
    ! Translated into the variables used here,
    ! a/b is a like a stabilized deconvolution.
    ! If a maximum time lag is given, then the following first option
    ! is more like a normalized cross-correlation, because the results
    ! embraces both positive and negative time lags.
    ! The second option is more of the traditional deconvolution, whose length
    ! is from 1 to nt (the same with the origional signal)
    allocate (z(1:n1, 1:n2, 1:n3))
    if (present(maxlag)) then
        z = fftshift(TTT(ifft(a/b)))
        call alloc_array(z, [-nlag1, nlag1, -nlag2, nlag2, -nlag3, nlag3], &
            source=z(nint((n1 + 1)/2.0) - nlag1:nint((n1 + 1)/2.0) + nlag1, &
            nint((n2 + 1)/2.0) - nlag2:nint((n2 + 1)/2.0) + nlag2, &
            nint((n3 + 1)/2.0) - nlag3:nint((n3 + 1)/2.0) + nlag3))
    else
        z = TTT(ifft(a/b))
        call alloc_array(z, [1, n1x, 1, n2x, 1, n3x], source=z(1:n1x, 1:n2x, 1:n3x))
    end if

end function deconv_3d_

#undef T
#undef TT
#undef TTT
#undef TTTT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef convfast_1d_
#undef convsum_1d_
#undef discrete_conv_1d_
#undef discrete_conv_2d_
#undef discrete_conv_3d_
#undef conv_1d_
#undef conv_2d_
#undef conv_3d_
#undef fftshift_1d_
#undef fftshift_2d_
#undef fftshift_3d_
#undef deconv_1d_
#undef deconv_2d_
#undef deconv_3d_

