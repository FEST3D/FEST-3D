    !< Higher order face state reconstruction method: MUSCL. 
module muscl
    !< MUSCL (Monotone Upwing Schemes for Scalar Conservation Laws is
    !< a scheme which replaces the piecewise constant approximation by
    !< reconstructing the states at the left and right side of each face.
    !< This is a one parameter upwind scheme which results in at most 3rd
    !< order accuracy.
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
    use global_vars, only : jlimiter_switch
    use global_vars, only : klimiter_switch
    use global_vars, only : itlimiter_switch
    use global_vars, only : jtlimiter_switch
    use global_vars, only : ktlimiter_switch
    use global_vars, only : iPB_switch
    use global_vars, only : jPB_switch
    use global_vars, only : kPB_switch

    implicit none
    private

    ! Private variables
    real, dimension(:, :, :, :), allocatable, target :: x_qp_left
      !< Store primitive state at the I-face left side
    real, dimension(:, :, :, :), allocatable, target :: x_qp_right
      !< Store primitive state at the I-face right side
    real, dimension(:, :, :, :), allocatable, target :: y_qp_left
      !< Store primitive state at the J-face left side
    real, dimension(:, :, :, :), allocatable, target :: y_qp_right
      !< Store primitive state at the J-face right side
    real, dimension(:, :, :, :), allocatable, target :: z_qp_left
      !< Store primitive state at the K-face left side
    real, dimension(:, :, :, :), allocatable, target :: z_qp_right
      !< Store primitive state at the K-face right side
    real :: phi, kappa
    
    real, dimension(:, :, :, :), pointer :: f_qp_left
    !< Generalized pointer for any I-J-K direction> f_qp_left can 
    !< either point to x_qp_left, y_qp_left or z_qp_left
    real, dimension(:, :, :, :), pointer :: f_qp_right
    !< Generalized pointer for any I-J-K direction> f_qp_right can 
    !< either point to x_qp_right, y_qp_right or z_qp_right
    real, dimension(:, :, :), allocatable :: pdif
    !< Used for pressure based witch
    integer :: switch_L=1
    !< Limiter switch 

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
          !< Allocate memoery to all array which store state
          !< the face.

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
          !< Deallocate all the array used 

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


        subroutine pressure_based_switching(f_dir)
          !< Pressure based switching. 
          !< User x,y, or z for I,J,or K face respectively
          !----------------------------------------------

            implicit none
            character, intent(in) :: f_dir
            !< Character can be x or y or z
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


        subroutine compute_face_state(f_dir, lam_switch, turb_switch)
          !< Subroutine to calculate state at the face, generalized for
          !< all direction : I,J, and K.
            implicit none
            ! Character can be x or y or z
            character, intent(in) :: f_dir
            !< Input direction x,y,or, z for which subroutine is called
            integer, intent(in) :: lam_switch
            !< Limiter switch for laminar variables
            integer, intent(in) :: turb_switch
            !< Limiter switch for turbulent variables
            integer :: i, j, k, l
            !< integer used for DO loop
            integer :: ii, jj, kk
            !< Variable for ALFA family limiter
            real :: alpha
            !< Flags to determine face direction
            real :: psi1, psi2
            !< limiters
            real :: fd
            !< forward difference
            real :: bd
            !< backward difference
            real :: r
            !< ratio of differences
            real, dimension(:, :, :, :), pointer :: f_qp_left
            !< Generalized pointer for any I-J-K direction> f_qp_left can 
            !< either point to x_qp_left, y_qp_left or z_qp_left
            real, dimension(:, :, :, :), pointer :: f_qp_right
            !< Generalized pointer for any I-J-K direction> f_qp_right can 
            !< either point to x_qp_right, y_qp_right or z_qp_right

            call dmsg(1, 'muscl', 'compute_face_state')

            select case (f_dir)
                case ('x')
                    f_qp_left(0:imx+1, 1:jmx-1, 1:kmx-1, 1:n_var) => x_qp_left
                    f_qp_right(0:imx+1, 1:jmx-1, 1:kmx-1, 1:n_var) => x_qp_right
                    ii = 1
                    jj = 0
                    kk = 0
                case ('y')
                    f_qp_left(1:imx-1, 0:jmx+1, 1:kmx-1, 1:n_var) => y_qp_left
                    f_qp_right(1:imx-1, 0:jmx+1, 1:kmx-1, 1:n_var) => y_qp_right
                    ii = 0
                    jj = 1
                    kk = 0
                case ('z')
                    f_qp_left(1:imx-1, 1:jmx-1, 0:kmx+1, 1:n_var) => z_qp_left
                    f_qp_right(1:imx-1, 1:jmx-1, 0:kmx+1, 1:n_var) => z_qp_right
                    ii = 0
                    jj = 0
                    kk = 1
                case default
                    call dmsg(5, 'muscl', 'compute_face_state', &
                            'Direction not recognised')
                    stop
            end select

            alpha = 2./3. !Koren limiter 
            phi = 1.0
            kappa = 1./3.
            switch_L=lam_switch

            do l = 1, n_var
              if(l>=6)then
                switch_L=turb_switch
              end if
             do k = 1-kk, kmx - 1 + kk
              do j = 1-jj, jmx - 1 + jj
               do i = 1-ii, imx - 1 + ii
                ! Cell based
                ! Koren limiter for now
                ! From paper: delta: forward difference 'fd'
                !             nabla: backward difference 'bd'
                fd = qp(i+ii, j+jj, k+kk, l) - qp(i, j, k, l)
                bd = qp(i, j, k, l) - qp(i-ii, j-jj, k-kk, l)

                r = fd / max(bd,1e-10)
                psi1 = max(0., min(2*r, alpha*(r-1.0) + 1.0, 2.))  !alpha limiter
!                psi1 = max(0., min(2*r,1.), min(r,2.))    ! superbee
!                psi1 = ((r*r) + r)/((r*r) + 1.0)          ! Van-Albda 
!                psi1 = (abs(r) + r)/(abs(r) + 1.0)          ! Van-Leer

                r = bd / max(fd, 1e-10)
                psi2 = max(0., min(2*r, alpha*(r-1.0) + 1.0, 2.))
!                psi2 = max(0., min(2*r,1.), min(r,2.))
!                psi2 = ((r*r) + r)/((r*r) + 1.0)          ! Van-Albda 
!                psi2 = (abs(r) + r)/(abs(r) + 1.0)          ! Van-Leer

                psi1 = (1 - (1 - psi1)*switch_L )
                psi2 = (1 - (1 - psi2)*switch_L )

                f_qp_left(i+ii, j+jj, k+kk, l) = qp(i, j, k, l) + 0.25*phi* &
                    (((1.-kappa) * psi1 * bd) + ((1.+kappa) * psi2 * fd))
                f_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
                    (((1.+kappa) * psi1 * bd) + ((1.-kappa) * psi2 * fd))
               end do
              end do
             end do
            end do


        end subroutine compute_face_state
        
        
        subroutine compute_muscl_states()
            !< Implement MUSCL scheme to get left and right states at
            !< each face. The computation is done through all cells
            !< and first level ghost cells
            !---------------------------------------------------------
            
            !call compute_xi_face_states()
            call compute_face_state('x', ilimiter_switch, itlimiter_switch)
            if(iPB_switch==1)then
              call pressure_based_switching('x')
            end if
            !call compute_eta_face_states()
            call compute_face_state('y', jlimiter_switch, jtlimiter_switch)
            if(jPB_switch==1)then
              call pressure_based_switching('y')
            end if
            !call compute_zeta_face_states()
            call compute_face_state('z', klimiter_switch, ktlimiter_switch)
            if(kPB_switch==1)then
              call pressure_based_switching('z')
            end if

        end subroutine compute_muscl_states

end module muscl
