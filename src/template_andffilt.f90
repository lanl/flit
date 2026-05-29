!
! © 2024-2026. Triad National Security, LLC. All rights reserved.
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

#define compute_structure_tensor_2d_     CONCAT(compute_structure_tensor_2d, T)
#define compute_structure_tensor_eigens_2d_     CONCAT(compute_structure_tensor_eigens_2d, T)
#define compute_diffusion_tensor_2d_     CONCAT(compute_diffusion_tensor_2d, T)
#define compute_diffusion_tensor_2d_     CONCAT(compute_diffusion_tensor_2d, T)
#define andf_fed_2d_     CONCAT(andf_fed_2d, T)
#define andf_fed_2d_mpi_     CONCAT(andf_fed_2d_mpi, T)

#define compute_structure_tensor_3d_     CONCAT(compute_structure_tensor_3d, T)
#define compute_structure_tensor_eigens_3d_     CONCAT(compute_structure_tensor_eigens_3d, T)
#define compute_diffusion_tensor_3d_     CONCAT(compute_diffusion_tensor_3d, T)
#define compute_diffusion_tensor_3d_     CONCAT(compute_diffusion_tensor_3d, T)
#define andf_fed_3d_     CONCAT(andf_fed_3d, T)
#define andf_fed_3d_mpi_     CONCAT(andf_fed_3d_mpi, T)

#define isprime_     CONCAT(isprime, T)
#define find_fed_steps_     CONCAT(find_fed_steps, T)
#define fine_fed_steps_by_process_time_     CONCAT(fine_fed_steps_by_process_time, T)

!
!> Compute general strcture tensor for 2D
!
subroutine compute_structure_tensor_2d_(u, s11, s12, s22)

    TT, dimension(:, :), intent(in) :: u
    TT, dimension(:, :), intent(inout) :: s11, s12, s22

    integer :: i, j
    TT :: d1u, d2u
    integer :: n1, n2

    n1 = size(u, 1)
    n2 = size(u, 2)

    !$omp parallel do private(i, j, d1u, d2u)
    do j = 1, n2
        do i = 1, n1

            if (i == 1) then
                d1u = u(i + 1, j) - u(i, j)
            else if (i == n1) then
                d1u = u(i, j) - u(i - 1, j)
            else
                d1u = (u(i + 1, j) - u(i - 1, j))*0.5
            end if

            if (j == 1) then
                d2u = u(i, j + 1) - u(i, j)
            else if (j == n2) then
                d2u = u(i, j) - u(i, j - 1)
            else
                d2u = (u(i, j + 1) - u(i, j - 1))*0.5
            end if

            ! Structure tensor components
            s11(i, j) = d1u**2
            s12(i, j) = d1u*d2u
            s22(i, j) = d2u**2

        end do
    end do
    !$omp end parallel do

end subroutine

!
!> Compute eigenvectors and eigenvalues from structure tensor for 2D
!
subroutine compute_structure_tensor_eigens_2d_(s11, s12, s22, ev1, ev2, ea1, ea2)

    TT, dimension(:, :), intent(in) :: s11, s12, s22
    type(TTT), dimension(:, :), intent(inout) :: ev1, ev2
    TT, dimension(:, :), intent(inout) :: ea1, ea2

    integer :: i, j
    TT, allocatable, dimension(:, :) :: strc, eigvec
    TT, allocatable, dimension(:) :: eigval
    integer :: n1, n2

    n1 = size(s11, 1)
    n2 = size(s11, 2)

    strc = zeros(2, 2)
    eigvec = zeros(2, 2)
    eigval = zeros(2)

    !$omp parallel do private(i, j, strc, eigvec, eigval)
    do j = 1, n2
        do i = 1, n1

            ! Structure tensor
            strc(1, 1) = s11(i, j)
            strc(1, 2) = s12(i, j)
            strc(2, 1) = s12(i, j)
            strc(2, 2) = s22(i, j)

            ! Solve for eigenvectors and eigenvalues
            call eigen_symm2x2(strc, eigval, eigvec)

            ! Assign to arrays
            ev1(i, j)%array = eigvec(:, 1)
            ev2(i, j)%array = eigvec(:, 2)

            ea1(i, j) = eigval(1)
            ea2(i, j) = eigval(2)

        end do
    end do
    !$omp end parallel do

end subroutine

!
!> Compute diffusion tensor for 2D
!
subroutine compute_diffusion_tensor_2d_(param, ev1, ev2, ea1, ea2, df11, df12, df22)

    type(andf_param), intent(in) :: param
    type(TTT), dimension(:, :), intent(in) :: ev1, ev2
    TT, dimension(:, :), intent(in) :: ea1, ea2
    TT, dimension(:, :), intent(inout) :: df11, df12, df22

    TT, dimension(1:2, 1:2) :: df
    TT, dimension(1:2, 1) :: ev
    integer :: i, j
    TT :: lambda1, lambda2
    integer :: n1, n2

    n1 = size(ev1, 1)
    n2 = size(ev2, 2)

    ! Anisotropic diffusion tensor
    !$omp parallel do private(i, j, df, ev, lambda1, lambda2)
    do j = 1, n2
        do i = 1, n1

            lambda1 = param%lambda1
            if (param%lambda2 >= 0) then
                lambda2 = param%lambda2
            else
                if (ea1(i, j) == ea2(i, j)) then
                    lambda2 = param%lambda1
                else
                    lambda2 = param%lambda1 + (1.0_fp - param%lambda1) &
                        *exp(-param%coefc/(ea1(i, j) - ea2(i, j))**(2.0_fp*param%powerm))
                end if
            end if

            ev(:, 1) = ev1(i, j)%array
            df = lambda1*matmul(ev, transpose(ev))
            ev(:, 1) = ev2(i, j)%array
            df = df + lambda2*matmul(ev, transpose(ev))

            df11(i, j) = df(1, 1)
            df12(i, j) = df(1, 2)
            df22(i, j) = df(2, 2)

        end do
    end do
    !$omp end parallel do

end subroutine

!
!> Comptue general structure tensor for 3D
!
subroutine compute_structure_tensor_3d_(u, s11, s12, s13, s22, s23, s33)

    TT, dimension(:, :, :), intent(in) :: u
    TT, dimension(:, :, :), intent(inout) :: s11, s12, s13, s22, s23, s33

    integer :: i, j, k
    TT :: d1u, d2u, d3u
    integer :: n1, n2, n3

    n1 = size(u, 1)
    n2 = size(u, 2)
    n3 = size(u, 3)

    !$omp parallel do private(i, j, k, d1u, d2u, d3u)
    do k = 1, n3
        do j = 1, n2
            do i = 1, n1

                if (i == 1) then
                    d1u = u(i + 1, j, k) - u(i, j, k)
                else if (i == n1) then
                    d1u = u(i, j, k) - u(i - 1, j, k)
                else
                    d1u = (u(i + 1, j, k) - u(i - 1, j, k))*0.5
                end if

                if (j == 1) then
                    d2u = u(i, j + 1, k) - u(i, j, k)
                else if (j == n2) then
                    d2u = u(i, j, k) - u(i, j - 1, k)
                else
                    d2u = (u(i, j + 1, k) - u(i, j - 1, k))*0.5
                end if

                if (k == 1) then
                    d3u = u(i, j, k + 1) - u(i, j, k)
                else if (k == n3) then
                    d3u = u(i, j, k) - u(i, j, k - 1)
                else
                    d3u = (u(i, j, k + 1) - u(i, j, k - 1))*0.5
                end if

                ! Structure tensor
                s11(i, j, k) = d1u**2
                s12(i, j, k) = d1u*d2u
                s13(i, j, k) = d1u*d3u
                s22(i, j, k) = d2u**2
                s23(i, j, k) = d2u*d3u
                s33(i, j, k) = d3u**2

            end do
        end do
    end do
    !$omp end parallel do

end subroutine

!
!> Compute eigenvectors and eigenvalues from structure tensor for 3D
!
subroutine compute_structure_tensor_eigens_3d_(s11, s12, s13, s22, s23, s33, ev1, ev2, ev3, ea1, ea2, ea3)

    TT, dimension(:, :, :), intent(in) :: s11, s12, s13, s22, s23, s33
    type(TTT), dimension(:, :, :), intent(inout) :: ev1, ev2, ev3
    TT, dimension(:, :, :), intent(inout) :: ea1, ea2, ea3

    integer :: i, j, k
    TT, allocatable, dimension(:, :) :: strc, eigvec
    TT, allocatable, dimension(:) :: eigval
    integer :: n1, n2, n3

    n1 = size(s11, 1)
    n2 = size(s11, 2)
    n3 = size(s11, 3)

    strc = zeros(3, 3)
    eigvec = zeros(3, 3)
    eigval = zeros(3)

    !$omp parallel do private(i, j, k, strc, eigvec, eigval)
    do k = 1, n3
        do j = 1, n2
            do i = 1, n1

                ! Structure tensor
                strc(1, 1) = s11(i, j, k)
                strc(1, 2) = s12(i, j, k)
                strc(1, 3) = s13(i, j, k)
                strc(2, 1) = s12(i, j, k)
                strc(2, 2) = s22(i, j, k)
                strc(2, 3) = s23(i, j, k)
                strc(3, 1) = s13(i, j, k)
                strc(3, 2) = s23(i, j, k)
                strc(3, 3) = s33(i, j, k)

                ! Solve for eigenvectors and eigenvalues
                call eigen_symm3x3(strc, eigval, eigvec)

                ! Assign to arrays
                ev1(i, j, k)%array = eigvec(:, 1)
                ev2(i, j, k)%array = eigvec(:, 2)
                ev3(i, j, k)%array = eigvec(:, 3)

                ea1(i, j, k) = eigval(1)
                ea2(i, j, k) = eigval(2)
                ea3(i, j, k) = eigval(3)

            end do
        end do
    end do
    !$omp end parallel do

end subroutine

!
!> Compute diffusion tensor for 3D
!
subroutine compute_diffusion_tensor_3d_(param, ev1, ev2, ev3, ea1, ea2, ea3, df11, df12, df13, df22, df23, df33)

    type(andf_param), intent(in) :: param
    type(TTT), dimension(:, :, :), intent(in) :: ev1, ev2, ev3
    TT, dimension(:, :, :), intent(in) :: ea1, ea2, ea3
    TT, dimension(:, :, :), intent(inout) :: df11, df12, df13, df22, df23, df33

    TT, dimension(1:3, 1:3) :: df
    TT, dimension(1:3, 1) :: ev
    integer :: i, j, k
    TT :: lambda1, lambda2, lambda3
    integer :: n1, n2, n3

    n1 = size(ev1, 1)
    n2 = size(ev1, 2)
    n3 = size(ev1, 3)

    ! Anisotropic diffusion tensor
    !$omp parallel do private(i, j, k, df, ev, lambda1, lambda2, lambda3)
    do k = 1, n3
        do j = 1, n2
            do i = 1, n1

                ! lambda 1
                lambda1 = param%lambda1
                ! lambda 2
                if (param%lambda2 >= 0) then
                    lambda2 = param%lambda2
                else
                    if (ea2(i, j, k) == ea3(i, j, k)) then
                        lambda2 = param%lambda1
                    else
                        lambda2 = param%lambda1 + (1.0_fp - param%lambda1) &
                            *exp(-param%coefc/(ea2(i, j, k) - ea3(i, j, k))**2)
                    end if
                end if
                ! lambda 3
                if (param%lambda3 >= 0) then
                    lambda3 = param%lambda3
                else
                    if (ea1(i, j, k) == ea3(i, j, k)) then
                        lambda3 = param%lambda1
                    else
                        lambda3 = param%lambda1 + (1.0_fp - param%lambda1) &
                            *exp(-param%coefc/(ea1(i, j, k) - ea3(i, j, k))**2)
                    end if
                end if

                ev(:, 1) = ev1(i, j, k)%array
                df = lambda1*matmul(ev, transpose(ev))
                ev(:, 1) = ev2(i, j, k)%array
                df = df + lambda2*matmul(ev, transpose(ev))
                ev(:, 1) = ev3(i, j, k)%array
                df = df + lambda3*matmul(ev, transpose(ev))

                df11(i, j, k) = df(1, 1)
                df12(i, j, k) = df(1, 2)
                df13(i, j, k) = df(1, 3)
                df22(i, j, k) = df(2, 2)
                df23(i, j, k) = df(2, 3)
                df33(i, j, k) = df(3, 3)

            end do
        end do
    end do
    !$omp end parallel do

end subroutine

function isprime_(w) result(tof)

    integer :: w
    logical :: tof

    integer :: s, i

    s = int(sqrt(w + 1.0))

    tof = .true.
    do i = 2, s
        if (mod(w, i) == 0) then
            tof = .false.
            exit
        end if
    end do

end function

subroutine find_fed_steps_(n, scale, tau_max, tau)

    integer, intent(in) :: n
    TT, intent(in) :: scale
    TT, intent(in) :: tau_max
    TT, dimension(:), intent(inout) :: tau

    integer :: k, kappa, prime, l, i
    TT :: c, d

    TT, allocatable, dimension(:) :: tauh

    ! Compute time saver
    c = 1.0_fp/(4.0_fp*n + 2.0_fp)
    d = scale*tau_max/2.0_fp

    allocate (tauh(1:n))

    ! Setup original ordered tau vector
    do k = 1, n
        tauh(k) = d/(cos(const_pi*(2.0_fp*(k - 1) + 1.0_fp)*c))**2
    end do

    ! Permute list of time steps
    if (n > fed_maxkappa + 1) then
        kappa = int(n/4.0_fp)
    else
        kappa = fed_kappalookup(n + 1)
    end if

    ! Get modulus for permutation
    prime = n + 1
    do while (.not. isprime_(prime))
        prime = prime + 1
    end do

    ! Perform permutation
    k = 0
    l = 0
    do while (l < n)

        i = mod((k + 1)*kappa, prime) - 1

        do while (i >= n)
            k = k + 1
            i = mod((k + 1)*kappa, prime) - 1
        end do

        tau(l + 1) = tauh(i + 1)

        k = k + 1
        l = l + 1

    end do

end subroutine

subroutine fine_fed_steps_by_process_time_(t, m, tau_max, tau, n)

    TT, intent(in) :: t
    integer, intent(in) :: m
    TT, intent(in) :: tau_max
    TT, allocatable, dimension(:), intent(inout) :: tau
    integer, intent(inout) :: n

    ! Compute necessary number of time steps
    n = int(ceiling(sqrt(3.0_fp*(t/m)/tau_max + 0.25_fp) - 0.5_fp - 1.0e-8) + 0.5_fp)

    ! Allocate memory
    if (allocated(tau)) then
        deallocate (tau)
    end if
    allocate (tau(1:n))

    ! Find FED steps
    call find_fed_steps_(n, 3.0_fp*(t/m)/(tau_max*(n*(n + 1))), tau_max, tau)

end subroutine

!
!> 2D nonlinear anisotropic diffusion filtering by fast explicit diffusion (FED)
!
function andf_fed_2d_(v, param, aux, df11, df12, df22, acoh) result(u)

    TT, dimension(:, :) :: v
    type(andf_param) :: param
    TT, dimension(:, :), optional :: aux, df11, df12, df22, acoh
    TT, allocatable, dimension(:, :) :: u

    integer :: n1, n2
    TT, allocatable, dimension(:, :) :: d1u, d2u
    TT, allocatable, dimension(:, :) :: s11, s12, s22
    TT, allocatable, dimension(:, :) :: coh
    TT, allocatable, dimension(:, :) :: d11, d12, d22, ea1, ea2
    type(TTT), allocatable, dimension(:, :) :: ev1, ev2
    integer :: nstep
    integer :: t, istep, i, j
    TT :: tmp1, tmp2
    TT, allocatable, dimension(:) :: tau
    TT :: scalar

    u = v

    n1 = size(u, 1)
    n2 = size(u, 2)
    scalar = maxval(abs(u))
    if (scalar == 0) then
        return
    else
        u = u/scalar
    end if

    allocate (d1u(1:n1, 1:n2))
    allocate (d2u(1:n1, 1:n2))
    allocate (d11(1:n1, 1:n2))
    allocate (d12(1:n1, 1:n2))
    allocate (d22(1:n1, 1:n2))
    allocate (s11(1:n1, 1:n2))
    allocate (s12(1:n1, 1:n2))
    allocate (s22(1:n1, 1:n2))
    allocate (ev1(1:n1, 1:n2))
    allocate (ev2(1:n1, 1:n2))
    allocate (ea1(1:n1, 1:n2))
    allocate (ea2(1:n1, 1:n2))
    allocate (coh(1:n1, 1:n2))

    ! Compute FED steps
    call fine_fed_steps_by_process_time_(param%sigma**2/2.0_fp, param%niter, 0.5_fp, tau, nstep)

    do t = 1, param%niter

        d1u = 0.0
        d2u = 0.0

        if (present(df11) .and. present(df12) .and. present(df22)) then
            ! If given diffusion tensor, then use
            d11 = df11
            d12 = df12
            d22 = df22

        else
            ! Otherwise compute diffusion tensor

            ! Compute structure tensor
            if (present(aux)) then
                if (maxval(abs(aux)) /= 0) then
                    call compute_structure_tensor_2d_(aux, s11, s12, s22)
                else
                    call compute_structure_tensor_2d_(u, s11, s12, s22)
                end if
            else
                call compute_structure_tensor_2d_(u, s11, s12, s22)
            end if

            ! Smooth structure tensor
            s11 = gauss_filt(s11, real([param%smooth1, param%smooth2], fp))
            s12 = gauss_filt(s12, real([param%smooth1, param%smooth2], fp))
            s22 = gauss_filt(s22, real([param%smooth1, param%smooth2], fp))

            ! Compute eigenvectors tensor
            call compute_structure_tensor_eigens_2d_(s11, s12, s22, ev1, ev2, ea1, ea2)

            ! Compute diffusion tensor
            call compute_diffusion_tensor_2d_(param, ev1, ev2, ea1, ea2, d11, d12, d22)

        end if

        if (param%powerm > 0) then

            if (present(acoh)) then
                if (maxval(abs(acoh)) /= 0) then
                    coh = acoh
                else
                    coh = ea1
                    where (coh == 0)
                        coh = float_huge
                    end where
                    coh = (ea1 - ea2)/coh
                end if
            else
                coh = ea1
                where (coh == 0)
                    coh = float_huge
                end where
                coh = (ea1 - ea2)/coh
            end if
            coh = clip(coh, 0.0_fp, 1.0_fp)

            !$omp parallel do private(i, j)
            do j = 1, n2
                do i = 1, n1
                    d11(i, j) = d11(i, j)*coh(i, j)**param%powerm
                    d12(i, j) = d12(i, j)*coh(i, j)**param%powerm
                    d22(i, j) = d22(i, j)*coh(i, j)**param%powerm
                end do
            end do
            !$omp end parallel do

        end if

        ! FED-based nonlinear anisotropic diffusion cycles
        do istep = 1, nstep

            ! Reflective boundary condition
            u(1, :) = u(2, :)
            u(n1, :) = u(n1 - 1, :)
            u(:, 1) = u(:, 2)
            u(:, n2) = u(:, n2 - 1)

            !$omp parallel do private(i, j, tmp1, tmp2)
            do j = 1, n2
                do i = 1, n1

                    if (i == 1) then
                        tmp1 = u(i + 1, j) - u(i, j)
                    else if (i == n1) then
                        tmp1 = u(i, j) - u(i - 1, j)
                    else
                        tmp1 = (u(i + 1, j) - u(i - 1, j))*0.5
                    end if

                    if (j == 1) then
                        tmp2 = u(i, j + 1) - u(i, j)
                    else if (j == n2) then
                        tmp2 = u(i, j) - u(i, j - 1)
                    else
                        tmp2 = (u(i, j + 1) - u(i, j - 1))*0.5
                    end if

                    d1u(i, j) = d11(i, j)*tmp1 + d12(i, j)*tmp2
                    d2u(i, j) = d12(i, j)*tmp1 + d22(i, j)*tmp2

                end do
            end do
            !$omp end parallel do

            !$omp parallel do private(i, j, tmp1, tmp2)
            do j = 2, n2 - 1
                do i = 2, n1 - 1

                    tmp1 = (d1u(i + 1, j) - d1u(i - 1, j))*0.5
                    tmp2 = (d2u(i, j + 1) - d2u(i, j - 1))*0.5

                    u(i, j) = u(i, j) + tau(istep)*(tmp1 + tmp2)

                end do
            end do
            !$omp end parallel do

        end do

        call warn(date_time_compact()//' Fast explicit anisotropic diffusion iteration '//num2str(t))

    end do

    u = u*scalar

end function

!
!> 3D nonlinear anisotropic diffusion filtering by fast explicit diffusion (FED)
!
function andf_fed_3d_(v, param, aux, df11, df12, df13, df22, df23, df33, acoh) result(u)

    TT, dimension(:, :, :) :: v
    type(andf_param) :: param
    TT, dimension(:, :, :), optional :: aux, df11, df12, df13, df22, df23, df33, acoh
    TT, allocatable, dimension(:, :, :) :: u

    integer :: n1, n2, n3
    TT, allocatable, dimension(:, :, :) :: d1u, d2u, d3u
    TT, allocatable, dimension(:, :, :) :: s11, s12, s13, s22, s23, s33
    TT, allocatable, dimension(:, :, :) :: coh
    TT, allocatable, dimension(:, :, :) :: d11, d12, d13, d22, d23, d33, ea1, ea2, ea3
    type(TTT), allocatable, dimension(:, :, :) :: ev1, ev2, ev3
    integer :: nstep
    integer :: t, istep, i, j, k
    TT :: tmp1, tmp2, tmp3
    TT, allocatable, dimension(:) :: tau
    TT :: scalar

    u = v

    n1 = size(u, 1)
    n2 = size(u, 2)
    n3 = size(u, 3)
    scalar = maxval(abs(u))
    if (scalar == 0) then
        return
    else
        u = u/scalar
    end if

    allocate (d1u(1:n1, 1:n2, 1:n3))
    allocate (d2u(1:n1, 1:n2, 1:n3))
    allocate (d3u(1:n1, 1:n2, 1:n3))
    allocate (d11(1:n1, 1:n2, 1:n3))
    allocate (d12(1:n1, 1:n2, 1:n3))
    allocate (d13(1:n1, 1:n2, 1:n3))
    allocate (d22(1:n1, 1:n2, 1:n3))
    allocate (d23(1:n1, 1:n2, 1:n3))
    allocate (d33(1:n1, 1:n2, 1:n3))
    allocate (s11(1:n1, 1:n2, 1:n3))
    allocate (s12(1:n1, 1:n2, 1:n3))
    allocate (s13(1:n1, 1:n2, 1:n3))
    allocate (s22(1:n1, 1:n2, 1:n3))
    allocate (s23(1:n1, 1:n2, 1:n3))
    allocate (s33(1:n1, 1:n2, 1:n3))
    allocate (ev1(1:n1, 1:n2, 1:n3))
    allocate (ev2(1:n1, 1:n2, 1:n3))
    allocate (ev3(1:n1, 1:n2, 1:n3))
    allocate (ea1(1:n1, 1:n2, 1:n3))
    allocate (ea2(1:n1, 1:n2, 1:n3))
    allocate (ea3(1:n1, 1:n2, 1:n3))
    allocate (coh(1:n1, 1:n2, 1:n3))

    ! Compute FED step sizes
    call fine_fed_steps_by_process_time_(param%sigma**2/2.0_fp, param%niter, 0.5_fp, tau, nstep)

    do t = 1, param%niter

        d1u = 0.0
        d2u = 0.0
        d3u = 0.0

        if (present(df11) .and. present(df12) .and. present(df13) &
                .and. present(df22) .and. present(df23) .and. present(df33)) then
            ! If given diffusion tensor, then use
            d11 = df11
            d12 = df12
            d13 = df13
            d22 = df22
            d23 = df23
            d33 = df33

        else

            ! Compute structure tensor
            if (present(aux)) then
                if (maxval(abs(aux)) /= 0) then
                    call compute_structure_tensor_3d_(aux, s11, s12, s13, s22, s23, s33)
                else
                    call compute_structure_tensor_3d_(u, s11, s12, s13, s22, s23, s33)
                end if
            else
                call compute_structure_tensor_3d_(u, s11, s12, s13, s22, s23, s33)
            end if

            ! Smooth structure tensor
            s11 = gauss_filt(s11, real([param%smooth1, param%smooth2, param%smooth3], fp))
            s12 = gauss_filt(s12, real([param%smooth1, param%smooth2, param%smooth3], fp))
            s13 = gauss_filt(s13, real([param%smooth1, param%smooth2, param%smooth3], fp))
            s22 = gauss_filt(s22, real([param%smooth1, param%smooth2, param%smooth3], fp))
            s23 = gauss_filt(s23, real([param%smooth1, param%smooth2, param%smooth3], fp))
            s33 = gauss_filt(s33, real([param%smooth1, param%smooth2, param%smooth3], fp))

            ! Compute eigenvectors tensor
            call compute_structure_tensor_eigens_3d_(s11, s12, s13, s22, s23, s33, ev1, ev2, ev3, ea1, ea2, ea3)

            ! Compute diffusion tensor
            call compute_diffusion_tensor_3d_(param, ev1, ev2, ev3, ea1, ea2, ea3, d11, d12, d13, d22, d23, d33)

            if (param%powerm > 0) then

                if (present(acoh)) then
                    if (maxval(abs(acoh)) /= 0) then
                        coh = acoh
                    else
                        coh = (ea1 + ea2)*(ea2 + ea3)
                        where (coh == 0)
                            coh = float_huge
                        end where
                        coh = 1.0_fp - 2.0_fp*ea2*(ea2 - ea3)/coh
                    end if
                else
                    coh = (ea1 + ea2)*(ea2 + ea3)
                    where (coh == 0)
                        coh = float_huge
                    end where
                    coh = 1.0_fp - 2.0_fp*ea2*(ea2 - ea3)/coh
                end if
                coh = clip(coh, 0.0_fp, 1.0_fp)

                !$omp parallel do private(i, j, k)
                do k = 1, n3
                    do j = 1, n2
                        do i = 1, n1
                            d11(i, j, k) = d11(i, j, k)*coh(i, j, k)**param%powerm
                            d12(i, j, k) = d12(i, j, k)*coh(i, j, k)**param%powerm
                            d13(i, j, k) = d13(i, j, k)*coh(i, j, k)**param%powerm
                            d22(i, j, k) = d22(i, j, k)*coh(i, j, k)**param%powerm
                            d23(i, j, k) = d23(i, j, k)*coh(i, j, k)**param%powerm
                            d33(i, j, k) = d33(i, j, k)*coh(i, j, k)**param%powerm
                        end do
                    end do
                end do
                !$omp end parallel do

            end if

        end if

        ! FED-based nonlinear anisotropic diffusion cycles
        do istep = 1, nstep

            ! Reflective boundary condition
            u(1, :, :) = u(2, :, :)
            u(n1, :, :) = u(n1 - 1, :, :)
            u(:, 1, :) = u(:, 2, :)
            u(:, n2, :) = u(:, n2 - 1, :)
            u(:, :, 1) = u(:, :, 2)
            u(:, :, n3) = u(:, :, n3 - 1)

            !$omp parallel do private(i, j, k, tmp1, tmp2, tmp3)
            do k = 1, n3
                do j = 1, n2
                    do i = 1, n1

                        if (i == 1) then
                            tmp1 = u(i + 1, j, k) - u(i, j, k)
                        else if (i == n1) then
                            tmp1 = u(i, j, k) - u(i - 1, j, k)
                        else
                            tmp1 = (u(i + 1, j, k) - u(i - 1, j, k))*0.5
                        end if

                        if (j == 1) then
                            tmp2 = u(i, j + 1, k) - u(i, j, k)
                        else if (j == n2) then
                            tmp2 = u(i, j, k) - u(i, j - 1, k)
                        else
                            tmp2 = (u(i, j + 1, k) - u(i, j - 1, k))*0.5
                        end if

                        if (k == 1) then
                            tmp3 = u(i, j, k + 1) - u(i, j, k)
                        else if (k == n3) then
                            tmp3 = u(i, j, k) - u(i, j, k - 1)
                        else
                            tmp3 = (u(i, j, k + 1) - u(i, j, k - 1))*0.5
                        end if

                        d1u(i, j, k) = d11(i, j, k)*tmp1 + d12(i, j, k)*tmp2 + d13(i, j, k)*tmp3
                        d2u(i, j, k) = d12(i, j, k)*tmp1 + d22(i, j, k)*tmp2 + d23(i, j, k)*tmp3
                        d3u(i, j, k) = d13(i, j, k)*tmp1 + d23(i, j, k)*tmp2 + d33(i, j, k)*tmp3

                    end do
                end do
            end do
            !$omp end parallel do

            !$omp parallel do private(i, j, k, tmp1, tmp2, tmp3)
            do k = 2, n3 - 1
                do j = 2, n2 - 1
                    do i = 2, n1 - 1

                        tmp1 = (d1u(i + 1, j, k) - d1u(i - 1, j, k))*0.5
                        tmp2 = (d2u(i, j + 1, k) - d2u(i, j - 1, k))*0.5
                        tmp3 = (d3u(i, j, k + 1) - d3u(i, j, k - 1))*0.5

                        u(i, j, k) = u(i, j, k) + tau(istep)*(tmp1 + tmp2 + tmp3)

                    end do
                end do
            end do
            !$omp end parallel do

        end do

        call warn(date_time_compact()//' Fast explicit anisotropic diffusion iteration '//num2str(t))

    end do

    u = u*scalar

end function

!
!> 2D nonlinear anisotropic diffusion filtering by fast explicit diffusion (FED)
!
function andf_fed_2d_mpi_(v, param, aux, df11, df12, df22, acoh) result(u)

    TT, dimension(:, :) :: v
    type(andf_param) :: param
    TT, dimension(:, :), optional :: aux, df11, df12, df22, acoh
    TT, allocatable, dimension(:, :) :: u

    integer :: n1, n2
    TT, allocatable, dimension(:, :) :: d1u, d2u
    TT, allocatable, dimension(:, :) :: s11, s12, s22
    TT, allocatable, dimension(:, :) :: coh
    TT, allocatable, dimension(:, :) :: d11, d12, d22, ea1, ea2
    type(TTT), allocatable, dimension(:, :) :: ev1, ev2
    integer :: nstep
    integer :: t, istep, i, j
    TT :: tmp1, tmp2
    TT, allocatable, dimension(:) :: tau
    integer :: n1beg, n1end, n2beg, n2end
    TT, allocatable, dimension(:, :) :: uu
    TT, allocatable, dimension(:, :) :: s11a, s12a, s22a
    TT, allocatable, dimension(:, :) :: auxu
    integer :: m1beg, m1end, m2beg, m2end
    TT :: scalar

    u = v

    n1 = size(u, 1)
    n2 = size(u, 2)
    scalar = maxval(abs(u))
    if (scalar == 0) then
        return
    else
        u = u/scalar
    end if

    ! divide domain
    call domain_decomp_regular(n1, n2, n1beg, n1end, n2beg, n2end)

    ! allocate memory
    call alloc_array(d1u, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(d2u, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(uu, [n1beg, n1end, n2beg, n2end], pad=1)
    allocate (d11(n1beg:n1end, n2beg:n2end))
    allocate (d12(n1beg:n1end, n2beg:n2end))
    allocate (d22(n1beg:n1end, n2beg:n2end))
    call alloc_array(s11, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(s12, [n1beg, n1end, n2beg, n2end], pad=1)
    call alloc_array(s22, [n1beg, n1end, n2beg, n2end], pad=1)
    allocate (ev1(n1beg:n1end, n2beg:n2end))
    allocate (ev2(n1beg:n1end, n2beg:n2end))
    allocate (ea1(n1beg:n1end, n2beg:n2end))
    allocate (ea2(n1beg:n1end, n2beg:n2end))
    allocate (coh(n1beg:n1end, n2beg:n2end))
    allocate (s11a(1:n1, 1:n2))
    allocate (s12a(1:n1, 1:n2))
    allocate (s22a(1:n1, 1:n2))
    call alloc_array(auxu, [n1beg, n1end, n2beg, n2end], pad=1)

    ! copy original
    m1beg = max(1, n1beg - 1)
    m1end = min(n1, n1end + 1)
    m2beg = max(1, n2beg - 1)
    m2end = min(n2, n2end + 1)

    uu(m1beg:m1end, m2beg:m2end) = u(m1beg:m1end, m2beg:m2end)
    u = 0.0

    ! Compute FED step sizes
    call fine_fed_steps_by_process_time_(param%sigma**2/2.0_fp, param%niter, 0.5_fp, tau, nstep)

    do t = 1, param%niter

        d1u = 0.0
        d2u = 0.0
        s11a = 0.0
        s12a = 0.0
        s22a = 0.0

        if (present(df11) .and. present(df12) .and. present(df22)) then
            ! If given diffusion tensor, then use
            d11 = df11(n1beg:n1end, n2beg:n2end)
            d12 = df12(n1beg:n1end, n2beg:n2end)
            d22 = df22(n1beg:n1end, n2beg:n2end)

        else

            ! Compute structure tensor
            if (present(aux)) then
                if (maxval(abs(aux)) /= 0) then
                    auxu(m1beg:m1end, m2beg:m2end) = aux(m1beg:m1end, m2beg:m2end)
                    call compute_structure_tensor_2d_( &
                        auxu(m1beg:m1end, m2beg:m2end), &
                        s11(m1beg:m1end, m2beg:m2end), &
                        s12(m1beg:m1end, m2beg:m2end), &
                        s22(m1beg:m1end, m2beg:m2end))
                else
                    call compute_structure_tensor_2d_( &
                        uu(m1beg:m1end, m2beg:m2end), &
                        s11(m1beg:m1end, m2beg:m2end), &
                        s12(m1beg:m1end, m2beg:m2end), &
                        s22(m1beg:m1end, m2beg:m2end))
                end if
            else
                call compute_structure_tensor_2d_( &
                    uu(m1beg:m1end, m2beg:m2end), &
                    s11(m1beg:m1end, m2beg:m2end), &
                    s12(m1beg:m1end, m2beg:m2end), &
                    s22(m1beg:m1end, m2beg:m2end))
            end if

            ! Smooth structure tensor
            ! The structure tensor is heterogeneous in space, and may be very heterogeneous.
            ! Therefore, to ensure sufficient accuracy,
            ! the smoothing is applied to the whole array, instead of blocks.
            s11a(n1beg:n1end, n2beg:n2end) = s11(n1beg:n1end, n2beg:n2end)
            s12a(n1beg:n1end, n2beg:n2end) = s12(n1beg:n1end, n2beg:n2end)
            s22a(n1beg:n1end, n2beg:n2end) = s22(n1beg:n1end, n2beg:n2end)

            call allreduce_array(s11a)
            call allreduce_array(s12a)
            call allreduce_array(s22a)

            s11a = gauss_filt(s11a, real([param%smooth1, param%smooth2], fp))
            s12a = gauss_filt(s12a, real([param%smooth1, param%smooth2], fp))
            s22a = gauss_filt(s22a, real([param%smooth1, param%smooth2], fp))

            s11(n1beg:n1end, n2beg:n2end) = s11a(n1beg:n1end, n2beg:n2end)
            s12(n1beg:n1end, n2beg:n2end) = s12a(n1beg:n1end, n2beg:n2end)
            s22(n1beg:n1end, n2beg:n2end) = s22a(n1beg:n1end, n2beg:n2end)

            ! Compute eigenvectors tensor
            call compute_structure_tensor_eigens_2d_( &
                s11(n1beg:n1end, n2beg:n2end), &
                s12(n1beg:n1end, n2beg:n2end), &
                s22(n1beg:n1end, n2beg:n2end), &
                ev1, ev2, ea1, ea2)

            ! Compute diffusion tensor
            call compute_diffusion_tensor_2d_(param, ev1, ev2, ea1, ea2, d11, d12, d22)

            if (param%powerm > 0) then

                if (present(acoh)) then
                    if (maxval(abs(acoh)) /= 0) then
                        coh = acoh(n1beg:n1end, n2beg:n2end)
                    else
                        coh = ea1
                        where (coh == 0)
                            coh = float_huge
                        end where
                        coh = (ea1 - ea2)/coh
                    end if
                else
                    coh = ea1
                    where (coh == 0)
                        coh = float_huge
                    end where
                    coh = (ea1 - ea2)/coh
                end if
                coh = clip(coh, 0.0_fp, 1.0_fp)

                !$omp parallel do private(i, j)
                do j = n2beg, n2end
                    do i = n1beg, n1end
                        d11(i, j) = d11(i, j)*coh(i, j)**param%powerm
                        d12(i, j) = d12(i, j)*coh(i, j)**param%powerm
                        d22(i, j) = d22(i, j)*coh(i, j)**param%powerm
                    end do
                end do
                !$omp end parallel do

            end if

        end if

        ! FED-based nonlinear anisotropic diffusion cycles
        do istep = 1, nstep

            ! Reflective boundary condition
            if (n1beg == 1) then
                uu(1, :) = uu(2, :)
            end if
            if (n1end == n1) then
                uu(n1, :) = uu(n1 - 1, :)
            end if
            if (n2beg == 1) then
                uu(:, 1) = uu(:, 2)
            end if
            if (n2end == n2) then
                uu(:, n2) = uu(:, n2 - 1)
            end if

            !$omp parallel do private(i, j, tmp1, tmp2)
            do j = n2beg, n2end
                do i = n1beg, n1end

                    if (i == 1) then
                        tmp1 = uu(i + 1, j) - uu(i, j)
                    else if (i == n1) then
                        tmp1 = uu(i, j) - uu(i - 1, j)
                    else
                        tmp1 = (uu(i + 1, j) - uu(i - 1, j))*0.5
                    end if

                    if (j == 1) then
                        tmp2 = uu(i, j + 1) - uu(i, j)
                    else if (j == n2) then
                        tmp2 = uu(i, j) - uu(i, j - 1)
                    else
                        tmp2 = (uu(i, j + 1) - uu(i, j - 1))*0.5
                    end if

                    d1u(i, j) = d11(i, j)*tmp1 + d12(i, j)*tmp2
                    d2u(i, j) = d12(i, j)*tmp1 + d22(i, j)*tmp2

                end do
            end do
            !$omp end parallel do

            if (n1beg <= n1end .and. n2beg <= n2end) then
                call commute_array(d1u, 1)
                call commute_array(d2u, 1)
            end if

            !$omp parallel do private(i, j, tmp1, tmp2)
            do j = max(2, n2beg), min(n2 - 1, n2end)
                do i = max(2, n1beg), min(n1 - 1, n1end)

                    tmp1 = (d1u(i + 1, j) - d1u(i - 1, j))*0.5
                    tmp2 = (d2u(i, j + 1) - d2u(i, j - 1))*0.5

                    uu(i, j) = uu(i, j) + tau(istep)*(tmp1 + tmp2)

                end do
            end do
            !$omp end parallel do

            if (n1beg <= n1end .and. n2beg <= n2end) then
                call commute_array(uu, 1)
            end if

            call mpi_barrier(mpi_comm_world, mpi_ierr)

        end do

        if (rankid == 0) then
            call warn(date_time_compact()//' Fast explicit anisotropic diffusion iteration '//num2str(t))
        end if

    end do

    ! copy back
    u(n1beg:n1end, n2beg:n2end) = uu(n1beg:n1end, n2beg:n2end)*scalar

    ! sum array
    call allreduce_array(u)

end function

!
!> 3D nonlinear anisotropic diffusion filtering by fast explicit diffusion (FED)
!
function andf_fed_3d_mpi_(v, param, aux, df11, df12, df13, df22, df23, df33, acoh) result(u)

    TT, dimension(:, :, :) :: v
    type(andf_param) :: param
    TT, dimension(:, :, :), optional :: aux, df11, df12, df13, df22, df23, df33, acoh
    TT, allocatable, dimension(:, :, :) :: u

    integer :: n1, n2, n3
    TT, allocatable, dimension(:, :, :) :: d1u, d2u, d3u
    TT, allocatable, dimension(:, :, :) :: s11, s12, s13, s22, s23, s33
    TT, allocatable, dimension(:, :, :) :: coh
    TT, allocatable, dimension(:, :, :) :: d11, d12, d13, d22, d23, d33, ea1, ea2, ea3
    type(TTT), allocatable, dimension(:, :, :) :: ev1, ev2, ev3
    integer :: nstep
    integer :: t, istep, i, j, k
    TT :: tmp1, tmp2, tmp3
    TT, allocatable, dimension(:) :: tau
    integer :: n1beg, n1end, n2beg, n2end, n3beg, n3end
    TT, allocatable, dimension(:, :, :) :: uu
    TT, allocatable, dimension(:, :, :) :: s11a, s12a, s13a, s22a, s23a, s33a
    TT, allocatable, dimension(:, :, :) :: auxu
    integer :: m1beg, m1end, m2beg, m2end, m3beg, m3end
    TT :: scalar

    u = v

    n1 = size(u, 1)
    n2 = size(u, 2)
    n3 = size(u, 3)
    scalar = maxval(abs(u))
    if (scalar == 0) then
        return
    else
        u = u/scalar
    end if

    ! divide domain
    call domain_decomp_regular(n1, n2, n3, n1beg, n1end, n2beg, n2end, n3beg, n3end)

    ! allocate memory
    call alloc_array(d1u, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(d2u, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(d3u, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(uu, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    allocate (d11(n1beg:n1end, n2beg:n2end, n3beg:n3end))
    allocate (d12(n1beg:n1end, n2beg:n2end, n3beg:n3end))
    allocate (d13(n1beg:n1end, n2beg:n2end, n3beg:n3end))
    allocate (d22(n1beg:n1end, n2beg:n2end, n3beg:n3end))
    allocate (d23(n1beg:n1end, n2beg:n2end, n3beg:n3end))
    allocate (d33(n1beg:n1end, n2beg:n2end, n3beg:n3end))
    call alloc_array(s11, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s12, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s13, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s22, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s23, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    call alloc_array(s33, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)
    allocate (ev1(n1beg:n1end, n2beg:n2end, n3beg:n3end))
    allocate (ev2(n1beg:n1end, n2beg:n2end, n3beg:n3end))
    allocate (ev3(n1beg:n1end, n2beg:n2end, n3beg:n3end))
    allocate (ea1(n1beg:n1end, n2beg:n2end, n3beg:n3end))
    allocate (ea2(n1beg:n1end, n2beg:n2end, n3beg:n3end))
    allocate (ea3(n1beg:n1end, n2beg:n2end, n3beg:n3end))
    allocate (coh(n1beg:n1end, n2beg:n2end, n3beg:n3end))
    allocate (s11a(1:n1, 1:n2, 1:n3))
    allocate (s12a(1:n1, 1:n2, 1:n3))
    allocate (s13a(1:n1, 1:n2, 1:n3))
    allocate (s22a(1:n1, 1:n2, 1:n3))
    allocate (s23a(1:n1, 1:n2, 1:n3))
    allocate (s33a(1:n1, 1:n2, 1:n3))
    call alloc_array(auxu, [n1beg, n1end, n2beg, n2end, n3beg, n3end], pad=1)

    ! copy original
    m1beg = max(1, n1beg - 1)
    m1end = min(n1, n1end + 1)
    m2beg = max(1, n2beg - 1)
    m2end = min(n2, n2end + 1)
    m3beg = max(1, n3beg - 1)
    m3end = min(n3, n3end + 1)

    uu(m1beg:m1end, m2beg:m2end, m3beg:m3end) = u(m1beg:m1end, m2beg:m2end, m3beg:m3end)
    u = 0.0

    ! Compute FED step sizes
    call fine_fed_steps_by_process_time_(param%sigma**2/2.0_fp, param%niter, 0.5_fp, tau, nstep)

    do t = 1, param%niter

        d1u = 0.0
        d2u = 0.0
        d3u = 0.0
        s11a = 0.0
        s12a = 0.0
        s13a = 0.0
        s22a = 0.0
        s23a = 0.0
        s33a = 0.0

        if (present(df11) .and. present(df12) .and. present(df13) &
                .and. present(df22) .and. present(df23) .and. present(df33)) then
            ! If given diffusion tensor, then use
            d11 = df11(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            d12 = df12(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            d13 = df13(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            d22 = df22(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            d23 = df23(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            d33 = df33(n1beg:n1end, n2beg:n2end, n3beg:n3end)

        else

            ! Compute structure tensor
            if (present(aux)) then
                if (maxval(abs(aux)) /= 0) then
                    auxu(m1beg:m1end, m2beg:m2end, m3beg:m3end) = aux(m1beg:m1end, m2beg:m2end, m3beg:m3end)
                    call compute_structure_tensor_3d_( &
                        auxu(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                        s11(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                        s12(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                        s13(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                        s22(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                        s23(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                        s33(m1beg:m1end, m2beg:m2end, m3beg:m3end))
                else
                    call compute_structure_tensor_3d_( &
                        uu(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                        s11(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                        s12(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                        s13(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                        s22(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                        s23(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                        s33(m1beg:m1end, m2beg:m2end, m3beg:m3end))
                end if
            else
                call compute_structure_tensor_3d_( &
                    uu(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                    s11(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                    s12(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                    s13(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                    s22(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                    s23(m1beg:m1end, m2beg:m2end, m3beg:m3end), &
                    s33(m1beg:m1end, m2beg:m2end, m3beg:m3end))
            end if

            ! Smooth structure tensor
            ! The structure tensor is heterogeneous in space, and may be very heterogeneous.
            ! Therefore, to ensure sufficient accuracy,
            ! the smoothing is applied to the whole array, instead of blocks.
            s11a(n1beg:n1end, n2beg:n2end, n3beg:n3end) = s11(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            s12a(n1beg:n1end, n2beg:n2end, n3beg:n3end) = s12(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            s13a(n1beg:n1end, n2beg:n2end, n3beg:n3end) = s13(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            s22a(n1beg:n1end, n2beg:n2end, n3beg:n3end) = s22(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            s23a(n1beg:n1end, n2beg:n2end, n3beg:n3end) = s23(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            s33a(n1beg:n1end, n2beg:n2end, n3beg:n3end) = s33(n1beg:n1end, n2beg:n2end, n3beg:n3end)

            call allreduce_array(s11a)
            call allreduce_array(s12a)
            call allreduce_array(s13a)
            call allreduce_array(s22a)
            call allreduce_array(s23a)
            call allreduce_array(s33a)

            s11a = gauss_filt(s11a, real([param%smooth1, param%smooth2, param%smooth3], fp))
            s12a = gauss_filt(s12a, real([param%smooth1, param%smooth2, param%smooth3], fp))
            s13a = gauss_filt(s13a, real([param%smooth1, param%smooth2, param%smooth3], fp))
            s22a = gauss_filt(s22a, real([param%smooth1, param%smooth2, param%smooth3], fp))
            s23a = gauss_filt(s23a, real([param%smooth1, param%smooth2, param%smooth3], fp))
            s33a = gauss_filt(s33a, real([param%smooth1, param%smooth2, param%smooth3], fp))

            s11(n1beg:n1end, n2beg:n2end, n3beg:n3end) = s11a(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            s12(n1beg:n1end, n2beg:n2end, n3beg:n3end) = s12a(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            s13(n1beg:n1end, n2beg:n2end, n3beg:n3end) = s13a(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            s22(n1beg:n1end, n2beg:n2end, n3beg:n3end) = s22a(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            s23(n1beg:n1end, n2beg:n2end, n3beg:n3end) = s23a(n1beg:n1end, n2beg:n2end, n3beg:n3end)
            s33(n1beg:n1end, n2beg:n2end, n3beg:n3end) = s33a(n1beg:n1end, n2beg:n2end, n3beg:n3end)

            ! Compute eigenvectors tensor
            call compute_structure_tensor_eigens_3d_( &
                s11(n1beg:n1end, n2beg:n2end, n3beg:n3end), &
                s12(n1beg:n1end, n2beg:n2end, n3beg:n3end), &
                s13(n1beg:n1end, n2beg:n2end, n3beg:n3end), &
                s22(n1beg:n1end, n2beg:n2end, n3beg:n3end), &
                s23(n1beg:n1end, n2beg:n2end, n3beg:n3end), &
                s33(n1beg:n1end, n2beg:n2end, n3beg:n3end), &
                ev1, ev2, ev3, ea1, ea2, ea3)

            ! Compute diffusion tensor
            call compute_diffusion_tensor_3d_(param, ev1, ev2, ev3, ea1, ea2, ea3, d11, d12, d13, d22, d23, d33)

            if (param%powerm > 0) then

                if (present(acoh)) then
                    if (maxval(abs(acoh)) /= 0) then
                        coh = acoh(n1beg:n1end, n2beg:n2end, n3beg:n3end)
                    else
                        coh = (ea1 + ea2)*(ea2 + ea3)
                        where (coh == 0)
                            coh = float_huge
                        end where
                        coh = 1.0_fp - 2.0_fp*ea2*(ea2 - ea3)/coh
                    end if
                else
                    coh = (ea1 + ea2)*(ea2 + ea3)
                    where (coh == 0)
                        coh = float_huge
                    end where
                    coh = 1.0_fp - 2.0_fp*ea2*(ea2 - ea3)/coh
                end if
                coh = clip(coh, 0.0_fp, 1.0_fp)

                !$omp parallel do private(i, j, k)
                do k = n3beg, n3end
                    do j = n2beg, n2end
                        do i = n1beg, n1end
                            d11(i, j, k) = d11(i, j, k)*coh(i, j, k)**param%powerm
                            d12(i, j, k) = d12(i, j, k)*coh(i, j, k)**param%powerm
                            d13(i, j, k) = d13(i, j, k)*coh(i, j, k)**param%powerm
                            d22(i, j, k) = d22(i, j, k)*coh(i, j, k)**param%powerm
                            d23(i, j, k) = d23(i, j, k)*coh(i, j, k)**param%powerm
                            d33(i, j, k) = d33(i, j, k)*coh(i, j, k)**param%powerm
                        end do
                    end do
                end do
                !$omp end parallel do

            end if

        end if

        ! FED-based nonlinear anisotropic diffusion cycles
        do istep = 1, nstep

            ! Reflective boundary condition
            if (n1beg == 1) then
                uu(1, :, :) = uu(2, :, :)
            end if
            if (n1end == n1) then
                uu(n1, :, :) = uu(n1 - 1, :, :)
            end if
            if (n2beg == 1) then
                uu(:, 1, :) = uu(:, 2, :)
            end if
            if (n2end == n2) then
                uu(:, n2, :) = uu(:, n2 - 1, :)
            end if
            if (n3beg == 1) then
                uu(:, :, 1) = uu(:, :, 2)
            end if
            if (n3end == n3) then
                uu(:, :, n3) = uu(:, :, n3 - 1)
            end if

            !$omp parallel do private(i, j, k, tmp1, tmp2, tmp3)
            do k = n3beg, n3end
                do j = n2beg, n2end
                    do i = n1beg, n1end

                        if (i == 1) then
                            tmp1 = uu(i + 1, j, k) - uu(i, j, k)
                        else if (i == n1) then
                            tmp1 = uu(i, j, k) - uu(i - 1, j, k)
                        else
                            tmp1 = (uu(i + 1, j, k) - uu(i - 1, j, k))*0.5
                        end if

                        if (j == 1) then
                            tmp2 = uu(i, j + 1, k) - uu(i, j, k)
                        else if (j == n2) then
                            tmp2 = uu(i, j, k) - uu(i, j - 1, k)
                        else
                            tmp2 = (uu(i, j + 1, k) - uu(i, j - 1, k))*0.5
                        end if

                        if (k == 1) then
                            tmp3 = uu(i, j, k + 1) - uu(i, j, k)
                        else if (k == n3) then
                            tmp3 = uu(i, j, k) - uu(i, j, k - 1)
                        else
                            tmp3 = (uu(i, j, k + 1) - uu(i, j, k - 1))*0.5
                        end if

                        d1u(i, j, k) = d11(i, j, k)*tmp1 + d12(i, j, k)*tmp2 + d13(i, j, k)*tmp3
                        d2u(i, j, k) = d12(i, j, k)*tmp1 + d22(i, j, k)*tmp2 + d23(i, j, k)*tmp3
                        d3u(i, j, k) = d13(i, j, k)*tmp1 + d23(i, j, k)*tmp2 + d33(i, j, k)*tmp3

                    end do
                end do
            end do
            !$omp end parallel do

            if (n1beg <= n1end .and. n2beg <= n2end .and. n3beg <= n3end) then
                call commute_array(d1u, 1)
                call commute_array(d2u, 1)
                call commute_array(d3u, 1)
            end if

            !$omp parallel do private(i, j, k, tmp1, tmp2, tmp3)
            do k = max(2, n3beg), min(n3 - 1, n3end)
                do j = max(2, n2beg), min(n2 - 1, n2end)
                    do i = max(2, n1beg), min(n1 - 1, n1end)

                        tmp1 = (d1u(i + 1, j, k) - d1u(i - 1, j, k))*0.5
                        tmp2 = (d2u(i, j + 1, k) - d2u(i, j - 1, k))*0.5
                        tmp3 = (d3u(i, j, k + 1) - d3u(i, j, k - 1))*0.5

                        uu(i, j, k) = uu(i, j, k) + tau(istep)*(tmp1 + tmp2 + tmp3)

                    end do
                end do
            end do
            !$omp end parallel do

            if (n1beg <= n1end .and. n2beg <= n2end .and. n3beg <= n3end) then
                call commute_array(uu, 1)
            end if

            call mpi_barrier(mpi_comm_world, mpi_ierr)

        end do

        if (rankid == 0) then
            call warn(date_time_compact()//' Fast explicit anisotropic diffusion iteration '//num2str(t))
        end if

    end do

    ! copy back
    u(n1beg:n1end, n2beg:n2end, n3beg:n3end) = uu(n1beg:n1end, n2beg:n2end, n3beg:n3end)*scalar

    ! sum array
    call allreduce_array(u)

end function

#undef T
#undef TT
#undef TTT
#undef fp

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef compute_structure_tensor_2d_
#undef compute_structure_tensor_eigens_2d_
#undef compute_diffusion_tensor_2d_
#undef compute_diffusion_tensor_2d_
#undef andf_fed_2d_
#undef andf_fed_2d_mpi_

#undef compute_structure_tensor_3d_
#undef compute_structure_tensor_eigens_3d_
#undef compute_diffusion_tensor_3d_
#undef compute_diffusion_tensor_3d_
#undef andf_fed_3d_
#undef andf_fed_3d_mpi_

#undef isprime_
#undef find_fed_steps_
#undef fine_fed_steps_by_process_time_

