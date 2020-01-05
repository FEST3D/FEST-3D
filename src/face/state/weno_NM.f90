    !< Higher face state reconstuction method: WENO for non uniform grids
module weno_NM
    !<
    !<Reference: 2 Huang, W. F., Ren, Y. X., and Jiang, X., 
    !<“A simple algorithm to improve the performance of the WENO scheme on non-uniform grids,” 
    !<Acta Mechanica Sinica/Lixue Xuebao, 2017, pp. 1–11.
    !-----------------------------------------------------------------

#include "../../debug.h"

    use vartypes
    implicit none
    private

    ! Public members
    public :: compute_weno_NM_states

    contains
        

        subroutine compute_face_states(qp, f_qp_left, f_qp_right, flags, cells, dims)
          !< Subroutine to calculate state at the face, generalized for
          !< all direction : I,J, and K.
            implicit none

            type(extent), intent(in) :: dims
            !< Extent of the domain:imx,jmx,kmx
            integer, dimension(3), intent(in) :: flags
            !< flags for direction switch
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in):: qp
            !< Store primitive variable at cell center
            real(wp), dimension(1-flags(1):dims%imx-1+2*flags(1), 1-flags(2):dims%jmx-1+2*flags(2), 1-flags(3):dims%kmx-1+2*flags(3), 1:dims%n_var), intent(inout) :: f_qp_left, f_qp_right
            !< primitive state at faces
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
            integer :: i, j, k, l
            integer :: i_f=0, j_f=0, k_f=0
            real(wp), dimension(3) :: P !< polynomial approximation
            real(wp), dimension(3) :: B !< smoothness factor
            real(wp), dimension(3) :: w !< wieght
            real(wp), dimension(3) :: g !< linear wieght
            real(wp), dimension(-2:2) :: u !<state_variable
            real(wp)               :: eps=1e-6

            real(wp), dimension(-2:2) :: vol
            real(wp)               :: U11
            real(wp)               :: U00
            real(wp)               :: U21
            real(wp)               :: U10
            real(wp)               :: U01
            real(wp)               :: U12
            real(wp)               :: alpha12
            real(wp)               :: alpha01
            real(wp)               :: alpha10
            real(wp)               :: alpha21

            g(1) = 1.0/10.0
            g(2) = 6.0/10.0
            g(3) = 3.0/10.0


            i_f = flags(1)
            j_f = flags(2)
            k_f = flags(3)

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


        subroutine compute_weno_NM_states(qp, x_qp_l, x_qp_r, y_qp_l, y_qp_r, z_qp_l, z_qp_r, cells, dims)
          !< Call Weno scheme for all the three direction I,J, and K

            implicit none
            type(extent), intent(in) :: dims
            !< Extent of the domain:imx,jmx,kmx
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in):: qp
            !< Store primitive variable at cell center
            real(wp), dimension(0:dims%imx+1,1:dims%jmx-1,1:dims%kmx-1,1:dims%n_var), intent(inout) :: x_qp_l, x_qp_r
            !< Store primitive state at the I-face 
            real(wp), dimension(1:dims%imx-1,0:dims%jmx+1,1:dims%kmx-1,1:dims%n_var), intent(inout) :: y_qp_l, y_qp_r
            !< Store primitive state at the J-face 
            real(wp), dimension(1:dims%imx-1,1:dims%jmx-1,0:dims%kmx+1,1:dims%n_var), intent(inout) :: z_qp_l, z_qp_r
            !< Store primitive state at the K-face 
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            integer, dimension(3) :: flags
            !< flags for different directions
            flags=(/1,0,0/)
            call compute_face_states(qp, x_qp_l, x_qp_r, flags, cells, dims)
            flags=(/0,1,0/)
            call compute_face_states(qp, y_qp_l, y_qp_r, flags, cells, dims)
            flags=(/0,0,1/)
            call compute_face_states(qp, z_qp_l, z_qp_r, flags, cells, dims)

        end subroutine compute_weno_NM_states

end module weno_NM
