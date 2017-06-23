module weno
    !-----------------------------------------------------------------
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
    use global_vars, only : process_id

    implicit none
    private

    ! Private variables
    real, dimension(:, :, :, :), allocatable, target :: x_qp_left, &
        x_qp_right, y_qp_left, y_qp_right, z_qp_left, z_qp_right
    real :: phi, kappa
    real, dimension(:, :, :, :), pointer :: f_qp_left, f_qp_right
    real, dimension(:, :, :), allocatable :: pdif
    !TODO: Convert to system of flags to write all 3 directions in a single subroutine

!   character(len=30) :: TVD_scheme

    ! Public members
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_weno_states
    public :: x_qp_left, x_qp_right
    public :: y_qp_left, y_qp_right
    public :: z_qp_left, z_qp_right

 !  TVD_scheme = trim('koren')

    contains
        
        subroutine setup_scheme()

        implicit none

        call dmsg(1, 'weno', 'setup_weno')

        phi = 1.0

        kappa = -1.0

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

            call dmsg(1, 'weno', 'destroy_weno')

            call dealloc(x_qp_left)
            call dealloc(x_qp_right)
            call dealloc(y_qp_left)
            call dealloc(y_qp_right)
            call dealloc(z_qp_left)
            call dealloc(z_qp_right)
            call dealloc(pdif)

        end subroutine destroy_scheme


        subroutine compute_face_states(dir)
            implicit none

            character(len=*), intent(in) :: dir
            integer :: i, j, k, l
            integer :: i_f=0, j_f=0, k_f=0
            real, dimension(3) :: P ! polynomial approximation
            real, dimension(3) :: B ! smoothness factor
            real, dimension(3) :: w ! wieght
            real, dimension(3) :: g ! linear wieght
            real, dimension(-2:2) :: u !state_variable
            real               :: eps=1e-6

            g(1) = 1./10.
            g(2) = 6./10.
            g(3) = 3./10.

            select case (dir)
              case('x')
              i_f=1
              j_f=0
              k_f=0
              f_qp_left  => x_qp_left
              f_qp_right => x_qp_right
              case('y')
              i_f=0
              j_f=1
              k_f=0
              f_qp_left  => y_qp_left
              f_qp_right => y_qp_right
              case('z')
              i_f=0
              j_f=0
              k_f=1
              f_qp_left  => z_qp_left
              f_qp_right => z_qp_right
            end select

            do l = 1, n_var
             do k = 1-k_f, kmx-1+k_f
              do j = 1-j_f, jmx-1+j_f
               do i = 1-i_f, imx-1+i_f
                 U(-2) = qp(i-2*i_f,j-2*j_f,k-2*k_f,l)  !u_{i-2}
                 U(-1) = qp(i-1*i_f,j-1*j_f,k-1*k_f,l)  !u_{i-1}
                 u( 0) = qp(i      ,j      ,k      ,l)  !u_{i}
                 U( 1) = qp(i+1*i_f,j+1*j_f,k+1*k_f,l)  !u_{i+1}
                 U( 2) = qp(i+2*i_f,j+2*j_f,k+2*k_f,l)  !u_{i+2}

                 P(1) = ( 2*U(-2) -  7*U(-1) + 11*U(0))/6
                 P(2) = (-1*U(-1) +  5*U( 0) +  2*U(1))/6
                 P(3) = ( 2*U( 0) +  5*U( 1) -  1*U(2))/6

                 B(1) =(13/12)*(U(-2)-2*U(-1)+U(0))**2 + (1/4)*(  U(-2)-4*U(-1)+3*U(0))**2
                 B(2) =(13/12)*(U(-1)-2*U( 0)+U(1))**2 + (1/4)*(  U(-1)-          U(1))**2
                 B(3) =(13/12)*(U( 0)-2*U( 1)+U(2))**2 + (1/4)*(3*U( 0)-4*U( 1)+  U(2))**2

                 w(:) = g(:)/(eps + B(:))**2 
                 f_qp_left(i+i_f,j+j_f,k+k_f,l)  = SUM(w*P)/SUM(w)

                 P(1) = ( 2*U(2) -  7*U( 1) + 11*U( 0))/6
                 P(2) = (-1*U(1) +  5*U( 0) +  2*U(-1))/6
                 P(3) = ( 2*U(0) +  5*U(-1) -  1*U(-2))/6

                 B(1) = (13/12)*(U( 2)-2*U( 1)+U( 0))**2 + (1/4)*(  U(2)-4*U( 1)+3*U( 0))**2
                 B(2) = (13/12)*(U( 1)-2*U( 0)+U(-1))**2 + (1/4)*(  U(1)-          U(-1))**2
                 B(3) = (13/12)*(U( 0)-2*U(-1)+U(-2))**2 + (1/4)*(3*U(0)-4*U(-1)+  U(-2))**2

                 w(:) = g(:)/(eps + B(:))**2 
                 f_qp_right(i,j,k,l) = SUM(w*P)/SUM(w)
               end do
              end do
             end do
            end do

        end subroutine compute_face_states


        subroutine compute_weno_states()

            call compute_face_states('x')
            call compute_face_states('y')
            call compute_face_states('z')

        end subroutine compute_weno_states

end module weno
