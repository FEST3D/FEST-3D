module muscl
    !-----------------------------------------------------------------
    ! MUSCL (Monotone Upwing Schemes for Scalar Conservation Laws is
    ! a scheme which replaces the piecewise constant approximation by
    ! reconstructing the states at the left and right side of each face.
    ! This is a one parameter upwind scheme which results in at most 3rd
    ! order accuracy.
    !
    ! The MUSCL scheme alone creates non-physical oscillations near 
    ! discontinuities like shocks. Hence, MUSCL is combined with
    ! some TVD (Total Variation Diminishing) to reduce such oscillations.
    ! TVD schemes also ensure that no new extrema of the state variables
    ! is created at the faces.
    !-----------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg

    use global_vars, only : imx
    use global_vars, only : jmx
    use global_vars, only : kmx
    use global_vars, only : grid_x
    use global_vars, only : grid_y
    use global_vars, only : grid_z

    use global_vars, only : qp
    use global_vars, only : n_var
    use global_vars, only : pressure
    use global_vars, only : pressure_inf
    use global_vars, only : ilimiter_switch
    use global_vars, only : PB_switch

    implicit none
    private

    ! Private variables
    real, dimension(:, :, :, :), allocatable, target :: x_qp_left, &
        x_qp_right, y_qp_left, y_qp_right, z_qp_left, z_qp_right
    real :: phi, kappa
    real, dimension(:, :, :, :), pointer :: f_qp_left, f_qp_right
    real, dimension(:, :, :), allocatable :: pdif
    integer :: switch_L=1
    integer :: switch_P=1
    !TODO: Convert to system of flags to write all 3 directions in a single subroutine

!   character(len=30) :: TVD_scheme

    ! Public members
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_muscl_states
    public :: x_qp_left, x_qp_right
    public :: y_qp_left, y_qp_right
    public :: z_qp_left, z_qp_right

 !  TVD_scheme = trim('koren')

    contains
        
        subroutine setup_scheme()

        implicit none

        call dmsg(1, 'muscl', 'setup_muscl')

        phi = 1.0

        kappa = 1./3.

        call alloc(x_qp_left, 0, imx+1, 1, jmx-1, 1, kmx-1, 1, n_var, &
            errmsg='Error: Unable to allocate memory for ' // &
                'x_qp_left.')
        call alloc(x_qp_right, 0, imx+1, 1, jmx-1, 1, kmx-1, 1, n_var, &
            errmsg='Error: Unable to allocate memory for ' // &
                'x_qp_right.')

        call alloc(y_qp_left, 1, imx-1, 0, jmx+1, 1, kmx-1, 1, n_var, &
            errmsg='Error: Unable to allocate memory for ' // &
                'y_qp_left.')
        call alloc(y_qp_right, 1, imx-1, 0, jmx+1, 1, kmx-1, 1, n_var, &
            errmsg='Error: Unable to allocate memory for ' // &
                'y_qp_right.')

        call alloc(z_qp_left, 1, imx-1, 1, jmx-1, 0, kmx+1, 1, n_var, &
            errmsg='Error: Unable to allocate memory for ' // &
                'z_qp_left.')
        call alloc(z_qp_right, 1, imx-1, 1, jmx-1, 0, kmx+1, 1, n_var, &
            errmsg='Error: Unable to allocate memory for ' // &
                'z_qp_right.')

        call alloc(pdif, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for' // &
                    'pdif')

        end subroutine setup_scheme


        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'muscl', 'destroy_muscl')

            call dealloc(x_qp_left)
            call dealloc(x_qp_right)
            call dealloc(y_qp_left)
            call dealloc(y_qp_right)
            call dealloc(z_qp_left)
            call dealloc(z_qp_right)
            call dealloc(pdif)

        end subroutine destroy_scheme


!       function min_mod(r)
!           
!           implicit none

!           real, intent(in) :: r
!           real :: min_mod
!           min_mod = min(1., (3 - kappa) * r / (1 - kappa))
!       
!       end function min_mod

!       function koren(r)

!           implicit none
!           
!           real, intent(in) :: r
!           real :: koren
!           koren = max(0., min(2*r, (2 + r)/3., 2.))

!       end function koren(r)


        subroutine compute_xi_face_states()

            implicit none

            integer :: i, j, k, l
            real :: psi1, psi2, fd, bd, r

            phi = 1.0
            kappa = 1./3.
            switch_L=ilimiter_switch

            do l = 1, n_var
              if(l==6  .or. l==7 )then
                switch_L=1
              end if
             do k = 1, kmx - 1
              do j = 1, jmx - 1
               do i = 0, imx 
                ! Cell based
                ! All faces interior only (even at boundaries)
                ! Hence (i=1, left) and (i=imx, right) will be dealt separately
                ! Koren limiter for now
                ! From paper: delta: forward difference 'fd'
                !             nabla: backward difference 'bd'
                fd = qp(i+1, j, k, l) - qp(i, j, k, l)
                bd = qp(i, j, k, l) - qp(i-1, j, k, l)
                r = fd / bd
!                psi1 = min(1., (3 - kappa) * r / (1 - kappa)) !minmod
                psi1 = max(0., min(2*r, (2 + r)/3., 2.))  !koren limiter
!                psi1 = max(0., min(2*r,1.), min(r,2.))    ! superbee
                r = bd / fd
!                psi2 = min(1., (3 - kappa) * r / (1 - kappa))
                psi2 = max(0., min(2*r, (2 + r)/3., 2.))
!                psi2 = max(0., min(2*r,1.), min(r,2.))
                psi1 = (1 - (1 - psi1)*switch_L )
                psi2 = (1 - (1 - psi2)*switch_L )

                x_qp_left(i+1, j, k, l) = qp(i, j, k, l) + 0.25*phi* &
                    (((1-kappa) * psi1 * bd) + ((1+kappa) * psi2 * fd))
                x_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
                    (((1+kappa) * psi1 * bd) + ((1-kappa) * psi2 * fd))
                !check for turbulent variables
               end do
              end do
             end do
            end do

        !   do k = 1, kmx - 1
        !    do j = 1, jmx - 1
        !       ! Exterior boundaries
        !       x_qp_left(1, j, k, :) = 0.5 * (qp(0, j, k, :) + qp(1, j, k, :))
        !       x_qp_right(imx, j, k, :) = 0.5 * (qp(imx-1, j, k, :) + qp(imx, j, k, :))
        !    end do
        !   end do

        end subroutine compute_xi_face_states


        subroutine compute_eta_face_states()

            implicit none

            integer :: i, j, k, l
            real :: psi1, psi2, fd, bd, r
            switch_L=ilimiter_switch

            do l = 1, n_var
              if(l==6  .or. l==7 )then
                switch_L=1
              end if
             do k = 1, kmx - 1
              do j = 0, jmx 
               do i = 1, imx - 1
                ! Cell based
                ! All faces interior only (even at boundaries)
                ! Hence (j=1, left) and (j=jmx, right) will be dealt separately
                ! Koren limiter for now
                ! From paper: delta: forward difference 'fd'
                !             nabla: backward difference 'bd'
                fd = qp(i, j+1, k, l) - qp(i, j, k, l)
                bd = qp(i, j, k, l) - qp(i, j-1, k, l)
                r = fd / bd
                psi1 = max(0., min(2*r, (2 + r)/3., 2.))
!                psi1 = max(0., min(2*r,1.), min(r,2.))
!                psi1 = min(1., (3 - kappa) * r / (1 - kappa))
                !psi1 = (1 - (1 - psi1)*ilimiter_switch )
                r = bd / fd
                psi2 = max(0., min(2*r, (2 + r)/3., 2.))
!                psi2 = max(0., min(2*r,1.), min(r,2.))
!                psi2 = min(1., (3 - kappa) * r / (1 - kappa))
                !psi2 = (1 - (1 - psi2)*ilimiter_switch )
                !check for turbulent variables
                psi1 = (1 - (1 - psi1)*switch_L )
                psi2 = (1 - (1 - psi2)*switch_L )

                y_qp_left(i, j+1, k, l) = qp(i, j, k, l) + 0.25*phi* &
                    (((1-kappa) * psi1 * bd) + ((1+kappa) * psi2 * fd))
                y_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
                    (((1+kappa) * psi1 * bd) + ((1-kappa) * psi2 * fd))
               end do
              end do
             end do
            end do

!           do k = 1, kmx - 1
!            do i = 1, imx - 1
!               ! Exterior boundaries
!               y_qp_left(i, 1, k, :) = 0.5 * (qp(i, 0, k, :) + qp(i, 1, k, :))
!               y_qp_right(i, jmx, k, :) = 0.5 * (qp(i, jmx-1, k, :) + qp(i, jmx, k, :))
!            end do
!           end do

        end subroutine compute_eta_face_states


        subroutine compute_zeta_face_states()

            implicit none

            real :: psi1, psi2, fd, bd, r
            integer :: i, j, k, l
            switch_L=ilimiter_switch
            !TODO: Figure out why in muscl, compute_zeta_face_states(), the order of looping matters

            do k = 0, kmx
             do l = 1, n_var
              if(l==6  .or. l==7 )then
                switch_L=1
              end if
              do j = 1, jmx - 1
               do i = 1, imx - 1
                ! Cell based
                ! All faces interior only (even at boundaries)
                ! Hence (k=1, left) and (k=kmx, right) will be dealt separately
                ! Koren limiter for now
                ! From paper: delta: forward difference 'fd'
                !             nabla: backward difference 'bd'
                fd = qp(i, j, k+1, l) - qp(i, j, k, l)
                bd = qp(i, j, k, l) - qp(i, j, k-1, l)
                r = fd / bd
                psi1 = max(0., min(2*r, (2 + r)/3., 2.))
!                psi1 = max(0., min(2*r,1.), min(r,2.))
!                psi1 = min(1., (3 - kappa) * r / (1 - kappa))
                !psi1 = (1 - (1 - psi1)*ilimiter_switch )
                r = bd / fd
                psi2 = max(0., min(2*r, (2 + r)/3., 2.))
!                psi2 = max(0., min(2*r,1.), min(r,2.))
!                psi2 = min(1., (3 - kappa) * r / (1 - kappa))
                !psi2 = (1 - (1 - psi2)*ilimiter_switch )
                !check for turbulent variables
                psi1 = (1 - (1 - psi1)*switch_L )
                psi2 = (1 - (1 - psi2)*switch_L )

                z_qp_left(i, j, k+1, l) = qp(i, j, k, l) + 0.25*phi* &
                    (((1-kappa) * psi1 * bd) + ((1+kappa) * psi2 * fd))
                z_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
                    (((1+kappa) * psi1 * bd) + ((1-kappa) * psi2 * fd))
               end do
              end do
             end do
            end do

!           do j = 1, jmx - 1
!            do i = 1, imx - 1
!               ! Exterior boundaries
!               z_qp_left(i, j, 1, :) = 0.5 * (qp(i, k, 0, :) + qp(i, j, 1, :))
!               z_qp_right(i, j, kmx, :) = 0.5 * (qp(i, j, kmx-1, :) + qp(i, j, kmx, :))
!            end do
!           end do

        end subroutine compute_zeta_face_states

        subroutine pressure_based_switching(f_dir)

            implicit none
            ! Character can be x or y or z
            character, intent(in) :: f_dir
            integer :: i, j, k, i_end, j_end, k_end
            integer :: i_f, j_f, k_f  ! Flags to determine face direction
            real :: pd2

            call dmsg(1, 'muscl', 'pressure_based_switching')

            select case (f_dir)
                case ('x')
                    f_qp_left => x_qp_left
                    f_qp_right => x_qp_right
                    i_f = 1
                    j_f = 0
                    k_f = 0
                    i_end = imx
                    j_end = jmx - 1
                    k_end = kmx - 1
                case ('y')
                    f_qp_left => y_qp_left
                    f_qp_right => y_qp_right
                    i_f = 0
                    j_f = 1
                    k_f = 0
                    i_end = imx - 1
                    j_end = jmx 
                    k_end = kmx - 1
                case ('z')
                    f_qp_left => z_qp_left
                    f_qp_right => z_qp_right
                    i_f = 0
                    j_f = 0
                    k_f = 1
                    i_end = imx - 1
                    j_end = jmx - 1 
                    k_end = kmx
                case default
                    call dmsg(5, 'ppm', 'pressure_based_switching', &
                            'Direction not recognised')
                    stop
            end select

            ! i_end and j_end denote number of faces
            ! Total number of cells including ghost_cells is
            ! (i_end+1) * j_end for xi faces and i_end*(j_end+1) for
            ! eta faces. 

            ! Loop over cells (physical)
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                pd2 = abs(pressure(i + i_f*1, j + j_f*1, k + k_f*1) - &
                          pressure(i - i_f*1, j - j_f*1, k - k_f*1))
                pdif(i, j, k) = 1 - (pd2/(pd2 + pressure_inf))
              end do
             end do
            end do

            ! Update at ghost cells
            pdif((1-i_f):(1-i_f)*(imx-1), (1-j_f):(1-j_f)*(jmx-1), &
                 (1-k_f):(1-k_f)*(kmx-1)) = &
                pdif(1:imx-1 - i_f*(imx-2), 1:jmx-1 - j_f*(jmx-2), &
                     1:kmx-1 - k_f*(kmx-2))
                    
            pdif((1-i_f*(-imx+1)):(i_f)+(imx-1), (1-j_f*(-jmx+1)):(j_f)+(jmx-1), &
                 (1-k_f*(-kmx+1)):(k_f)+(kmx-1)) = &
                pdif(1+i_f*(imx-2):imx-1 , 1+j_f*(jmx-2):jmx-1, &
                     1+k_f*(kmx-2):kmx-1)
             

            ! Loop over faces
            do k = 1, kmx - (1 - k_f)            
             do j = 1, jmx - (1 - j_f)
              do i = 1, imx - (1 - i_f)
                f_qp_left(i, j, k, :) = qp(i - i_f*1, j - j_f*1, k - k_f*1, :) + (&
                    pdif(i - i_f*1, j - j_f*1, k - k_f*1) * ( &
                    f_qp_left(i, j, k, :) - qp(i - i_f*1, j - j_f*1, k - k_f*1, :)))

                f_qp_right(i, j, k, :) = qp(i, j, k, :) - (&
                    pdif(i, j, k) * ( &
                    qp(i, j, k, :) - f_qp_right(i, j, k, :)))
              end do
             end do
            end do

        end subroutine pressure_based_switching
        
        subroutine compute_muscl_states()
            !---------------------------------------------------------
            ! Implement MUSCL scheme to get left and right states at
            ! each face. The computation is done through all cells
            ! and first level ghost cells
            !---------------------------------------------------------
            
            call compute_xi_face_states()
            if(PB_switch==1)then
              call pressure_based_switching('x')
            end if
            call compute_eta_face_states()
            if(PB_switch==1)then
              call pressure_based_switching('y')
            end if
            call compute_zeta_face_states()
            if(PB_switch==1)then
              call pressure_based_switching('z')
            end if

        end subroutine compute_muscl_states

end module muscl
