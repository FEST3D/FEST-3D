    !< Higher face state reconstuction method: WENO for non uniform grids
module weno_NM
    !<
    !<Reference: 2 Huang, W. F., Ren, Y. X., and Jiang, X., 
    !<“A simple algorithm to improve the performance of the WENO scheme on non-uniform grids,” 
    !<Acta Mechanica Sinica/Lixue Xuebao, 2017, pp. 1–11.
    !-----------------------------------------------------------------

#include "../../debug.h"

    use vartypes
    use utils, only: alloc
!    use global_vars, only : volume
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

    integer :: imx, jmx, kmx, n_var
    ! Public members
    public :: setup_scheme
!    public :: destroy_scheme
    public :: compute_weno_NM_states
    public :: x_qp_left, x_qp_right
    public :: y_qp_left, y_qp_right
    public :: z_qp_left, z_qp_right

 !  TVD_scheme = trim('koren')

    contains
        
        subroutine setup_scheme(control, dims)
          !< Allocate memoery to all array which store state
          !< the face.

        implicit none
        type(controltype), intent(in) :: control
        type(extent), intent(in) :: dims

        DebugCall('setup_weno_NM')

        imx = dims%imx
        jmx = dims%jmx
        kmx = dims%kmx

        n_var = control%n_var

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


!        subroutine destroy_scheme()
!          !< Deallocate all the array used 
!
!            implicit none
!
!            DebugCall('destroy_weno_NM')
!
!            call dealloc(x_qp_left)
!            call dealloc(x_qp_right)
!            call dealloc(y_qp_left)
!            call dealloc(y_qp_right)
!            call dealloc(z_qp_left)
!            call dealloc(z_qp_right)
!
!        end subroutine destroy_scheme
!

        subroutine compute_face_states(qp, cells, dims, dir)
          !< Subroutine to calculate state at the face, generalized for
          !< all direction : I,J, and K.
            implicit none

            type(extent), intent(in) :: dims
            real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in):: qp
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
            character(len=*), intent(in) :: dir
            integer :: i, j, k, l
            integer :: i_f=0, j_f=0, k_f=0
            real, dimension(3) :: P !< polynomial approximation
            real, dimension(3) :: B !< smoothness factor
            real, dimension(3) :: w !< wieght
            real, dimension(3) :: g !< linear wieght
            real, dimension(-2:2) :: u !<state_variable
            real               :: eps=1e-6

            real, dimension(-2:2) :: vol
            real               :: U11
            real               :: U00
            real               :: U21
            real               :: U10
            real               :: U01
            real               :: U12
            real               :: alpha12
            real               :: alpha01
            real               :: alpha10
            real               :: alpha21

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

            do l = 1, dims%n_var
             do k = 1-k_f, dims%kmx-1+k_f
              do j = 1-j_f, dims%jmx-1+j_f
               do i = 1-i_f, dims%imx-1+i_f
                 U(-2) = qp(i-2*i_f,j-2*j_f,k-2*k_f,l)  !u_{i-2}
                 U(-1) = qp(i-1*i_f,j-1*j_f,k-1*k_f,l)  !u_{i-1}
                 u( 0) = qp(i      ,j      ,k      ,l)  !u_{i}
                 U( 1) = qp(i+1*i_f,j+1*j_f,k+1*k_f,l)  !u_{i+1}
                 U( 2) = qp(i+2*i_f,j+2*j_f,k+2*k_f,l)  !u_{i+2}

                 Vol(-2) = cells(i-2*i_f,j-2*j_f,k-2*k_f)%volume  !volume_{i-2}
                 Vol(-1) = cells(i-1*i_f,j-1*j_f,k-1*k_f)%volume  !volume_{i-1}
                 Vol( 0) = cells(i      ,j      ,k      )%volume  !volume_{i}
                 Vol( 1) = cells(i+1*i_f,j+1*j_f,k+1*k_f)%volume  !volume_{i+1}
                 Vol( 2) = cells(i+2*i_f,j+2*j_f,k+2*k_f)%volume  !volume_{i+2}

                 alpha12 = Vol( 2)/(Vol( 1) + Vol( 2))
                 alpha01 = Vol( 1)/(Vol( 0) + Vol( 1))
                 alpha10 = Vol( 0)/(Vol(-1) + Vol( 0))
                 alpha21 = vol(-1)/(Vol(-2) + Vol(-1))

                 U01     = (1.0-alpha01)*U(0) + alpha01*U(1)
                 U12     = (1.0-alpha12)*U(1) + alpha12*U(2)
                 U10     = (1.0-alpha10)*U(-1) + alpha10*U(0)
                 U21     = (1.0-alpha21)*U(-2) + alpha21*U(-1)
                 U00     = U(-1) + (1.0-alpha21)*(U(-1) - U(-2))
                 U11     = U( 1) + alpha12*(U(1) - U(2))

                 P(1) = ( 6.0*U(0) -  1.0*U10  - 2.0*U00)/3.0
                 P(2) = (-1.0*U10  +  2.0*U(0) + 2.0*U01)/3.0
                 P(3) = ( 2.0*U01  +  2.0*U(1) - 1.0*U12)/3.0

                 B(1) = (13.0/12.0)*(2*U10-2.0*U00       )**2 + (1.0/4.0)*(4*U(0)-2.0*U10  -2.0*U00)**2
                 B(2) = (13.0/12.0)*(2*U10-4.0*U(0)+2*U01)**2 + (1.0/4.0)*(-2*u10          +2.0*U01)**2
                 B(3) = (13.0/12.0)*(2*U01-4.0*U(1)+2*U12)**2 + (1.0/4.0)*(-6*U01+8.0*U( 1)-2.0*U12)**2

                 w(:) = g(:)/(eps + B(:))**2 
                 f_qp_left(i+i_f,j+j_f,k+k_f,l)  = SUM(w*P)/SUM(w)

                 P(1) = ( 6.0*U(0) -  1.0*U01  - 2.0*U11)/3.0
                 P(2) = (-1.0*U01  +  2.0*U(0) + 2.0*U10)/3.0
                 P(3) = ( 2.0*U10  +  2.0*U(-1) - 1.0*U21)/3.0

                 B(1) = (13.0/12.0)*(2*U01-2.0*U11        )**2 + (1.0/4.0)*(4*U(0)-2.0*U01  -2.0*U11)**2
                 B(2) = (13.0/12.0)*(2*U01-4.0*U( 0)+2*U10)**2 + (1.0/4.0)*(-2*u01          +2.0*U10)**2
                 B(3) = (13.0/12.0)*(2*U10-4.0*U(-1)+2*U21)**2 + (1.0/4.0)*(-6*U10+8.0*U(-1)-2.0*U21)**2

                 w(:) = g(:)/(eps + B(:))**2 
                 f_qp_right(i,j,k,l) = SUM(w*P)/SUM(w)
               end do
              end do
             end do
            end do

        end subroutine compute_face_states


        subroutine compute_weno_NM_states(qp, cells, dims)
          !< Call Weno scheme for all the three direction I,J, and K

            implicit none
            type(extent), intent(in) :: dims
            real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in):: qp
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            call compute_face_states(qp, cells, dims, 'x')
            call compute_face_states(qp, cells, dims, 'y')
            call compute_face_states(qp, cells, dims, 'z')

        end subroutine compute_weno_NM_states

end module weno_NM
