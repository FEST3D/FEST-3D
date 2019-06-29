    !< Higher face state reconstuction method: WENO
module weno
    !<
    !< Reference: 1 Shu, C.-W., “High-order Finite Difference and Finite Volume 
    !< WENO Schemes and Discontinuous Galerkin Methods for CFD,” 
    !< International Journal of Computational Fluid Dynamics, vol. 17, 2003, pp. 107–118.
    !< Reference: 2 Huang, W. F., Ren, Y. X., and Jiang, X., 
    !<“A simple algorithm to improve the performance of the WENO scheme on non-uniform grids,” 
    !<Acta Mechanica Sinica/Lixue Xuebao, 2017, pp. 1–11.
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
    use global_vars, only : process_id

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
    
    real, dimension(:, :, :, :), pointer :: f_qp_left
    !< Generalized pointer for any I-J-K direction> f_qp_left can 
    !< either point to x_qp_left, y_qp_left or z_qp_left
    real, dimension(:, :, :, :), pointer :: f_qp_right
    !< Generalized pointer for any I-J-K direction> f_qp_right can 
    !< either point to x_qp_right, y_qp_right or z_qp_right

    ! Public members
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_weno_states
    public :: x_qp_left, x_qp_right
    public :: y_qp_left, y_qp_right
    public :: z_qp_left, z_qp_right

    contains
        
        subroutine setup_scheme()
          !< Allocate memoery to all array which store state
          !< the face.

        implicit none

        call dmsg(1, 'weno', 'setup_weno')

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


        end subroutine setup_scheme


        subroutine destroy_scheme()
          !< Deallocate all the array used 

            implicit none

            call dmsg(1, 'weno', 'destroy_weno')

            call dealloc(x_qp_left)
            call dealloc(x_qp_right)
            call dealloc(y_qp_left)
            call dealloc(y_qp_right)
            call dealloc(z_qp_left)
            call dealloc(z_qp_right)

        end subroutine destroy_scheme


        subroutine compute_face_states(dir)
          !< Subroutine to calculate state at the face, generalized for
          !< all direction : I,J, and K.
            implicit none

            character(len=*), intent(in) :: dir
            integer :: i, j, k, l
            integer :: i_f=0, j_f=0, k_f=0
            real, dimension(3) :: P !< polynomial approximation
            real, dimension(3) :: B !< smoothness factor
            real, dimension(3) :: w !< wieght
            real, dimension(3) :: g !< linear wieght
            real, dimension(-2:2) :: u !< state_variable
            real               :: eps=1e-6

            g(1) = 1.0/10.0
            g(2) = 6.0/10.0
            g(3) = 3.0/10.0

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

                 P(1) = ( 2.0*U(-2) -  7.0*U(-1) + 11.0*U(0))/6.0
                 P(2) = (-1.0*U(-1) +  5.0*U( 0) +  2.0*U(1))/6.0
                 P(3) = ( 2.0*U( 0) +  5.0*U( 1) -  1.0*U(2))/6.0

                 B(1) = (13.0/12.0)*(U(-2)-2.0*U(-1)+U(0))**2 + (1.0/4.0)*(    U(-2)-4.0*U(-1)+3.0*U(0))**2
                 B(2) = (13.0/12.0)*(U(-1)-2.0*U( 0)+U(1))**2 + (1.0/4.0)*(    U(-1)-              U(1))**2
                 B(3) = (13.0/12.0)*(U( 0)-2.0*U( 1)+U(2))**2 + (1.0/4.0)*(3.0*U( 0)-4.0*U( 1)+    U(2))**2

                 w(:) = g(:)/(eps + B(:))**2 
                 f_qp_left(i+i_f,j+j_f,k+k_f,l)  = SUM(w*P)/SUM(w)

                 P(1) = ( 2.0*U(2) -  7.0*U( 1) + 11.0*U( 0))/6.0 
                 P(2) = (-1.0*U(1) +  5.0*U( 0) +  2.0*U(-1))/6.0
                 P(3) = ( 2.0*U(0) +  5.0*U(-1) -  1.0*U(-2))/6.0

                 !B(1) = (13.0/12.0)*(U( 2)-2.0*U( 1)+U( 0))**2 + (1.0/4.0)*(    U(2)-4.0*U( 1)+3.0*U( 0))**2
                 !B(2) = (13.0/12.0)*(U( 1)-2.0*U( 0)+U(-1))**2 + (1.0/4.0)*(    U(1)-              U(-1))**2
                 !B(3) = (13.0/12.0)*(U( 0)-2.0*U(-1)+U(-2))**2 + (1.0/4.0)*(3.0*U(0)-4.0*U(-1)+    U(-2))**2

                 w(1) = g(1)/(eps + B(3))**2 
                 w(2) = g(2)/(eps + B(2))**2 
                 w(3) = g(3)/(eps + B(1))**2 
                 f_qp_right(i,j,k,l) = SUM(w*P)/SUM(w)
               end do
              end do
             end do
            end do

        end subroutine compute_face_states


        subroutine compute_weno_states()
          !< Call Weno scheme for all the three direction I,J, and K

            call compute_face_states('x')
            call compute_face_states('y')
            call compute_face_states('z')

        end subroutine compute_weno_states

end module weno
