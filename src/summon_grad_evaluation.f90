 !< Calculate gradients of any primitive variables and temperature
module summon_grad_evaluation
 !< Calculate gradients of any primitive variables and temperature
  !----------------------------------------------------------
  !170608  -Jatinder Pal Singh Sandhu
  ! Aim : call is made to all the required gradients
  !                         based on input conditions
  !----------------------------------------------------------
#include "error.inc"
  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx
  use global_vars, only : x_speed
  use global_vars, only : y_speed
  use global_vars, only : z_speed
  use global_vars, only : tk
  use global_vars, only : tw
  use global_vars, only : tv
  use global_vars, only : te
  use global_vars, only : tkl
  use global_vars, only : tgm
  use global_vars, only : turbulence
  use global_vars, only : transition
  use global_vars, only : gradu_x
  use global_vars, only : gradu_y
  use global_vars, only : gradu_z
  use global_vars, only : gradv_x
  use global_vars, only : gradv_y
  use global_vars, only : gradv_z
  use global_vars, only : gradw_x
  use global_vars, only : gradw_y
  use global_vars, only : gradw_z
  use global_vars, only : gradT_x
  use global_vars, only : gradT_y
  use global_vars, only : gradT_z
  use global_vars, only : gradtk_x
  use global_vars, only : gradtk_y
  use global_vars, only : gradtk_z
  use global_vars, only : gradtw_x
  use global_vars, only : gradtw_y
  use global_vars, only : gradtw_z
  use global_vars, only : gradtv_x
  use global_vars, only : gradtv_y
  use global_vars, only : gradtv_z
  use global_vars, only : gradte_x
  use global_vars, only : gradte_y
  use global_vars, only : gradte_z
  use global_vars, only : gradtkl_x
  use global_vars, only : gradtkl_y
  use global_vars, only : gradtkl_z
  use global_vars, only : gradtgm_x
  use global_vars, only : gradtgm_y
  use global_vars, only : gradtgm_z
  use global_vars, only : gradqp_z
  use global_vars, only : process_id
  use global_vars, only : xn
  use global_vars, only : yn
  use global_vars, only : zn
  use global_vars, only : xnx
  use global_vars, only : xny
  use global_vars, only : xnz
  use global_vars, only : ynx
  use global_vars, only : yny
  use global_vars, only : ynz
  use global_vars, only : znx
  use global_vars, only : zny
  use global_vars, only : znz
  use global_vars, only : xA
  use global_vars, only : yA
  use global_vars, only : zA
  use global_vars, only : volume
  use global_vars, only : density
  use global_vars, only : pressure
  use global_vars, only : R_gas
  use global_vars, only : gm
  use utils      , only : alloc
  use utils      , only : dealloc
  use utils      , only : dmsg
  use utils      , only : turbulence_read_error
  use string
  use ghost_gradients, only : apply_gradient_bc 

  implicit none
  private

  real, dimension(6)               :: T
  !< Temperaure array for six neighbours
  real                             :: cell_T
  !< Temperature at cell center
  integer :: i,j,k
  !< integer for DO loop
  public :: evaluate_all_gradients

  contains

    subroutine evaluate_all_gradients()
      !< Call to all the required gradients and 
      !< apply boundary condition for ghost cell
      !< gradients

      implicit none

      call dmsg(1, 'summon_grad_evaluation', 'evaluate_all_gradients')

      call compute_gradient_G(gradu_x, x_speed, 'x')
      call compute_gradient_G(gradv_x, y_speed, 'x')
      call compute_gradient_G(gradw_x, z_speed, 'x')
      call compute_gradient_T(gradT_x         , 'x')
      call compute_gradient_G(gradu_y, x_speed, 'y')
      call compute_gradient_G(gradv_y, y_speed, 'y')
      call compute_gradient_G(gradw_y, z_speed, 'y')
      call compute_gradient_T(gradT_y         , 'y')
      if(kmx>2) then
      call compute_gradient_G(gradu_z, x_speed, 'z')
      call compute_gradient_G(gradv_z, y_speed, 'z')
      call compute_gradient_G(gradw_z, z_speed, 'z')
      call compute_gradient_T(gradT_z         , 'z')
      else
       gradqp_z=0.0
      end if

!      include "compute_gradu_x.inc"
!      include "compute_gradu_y.inc"
!      include "compute_gradv_x.inc"
!      include "compute_gradv_y.inc"
!      include "compute_gradw_x.inc"
!      include "compute_gradw_y.inc"
!      include "compute_gradT_x.inc"
!      include "compute_gradT_y.inc"
!      if(kmx>2) then
!      include "compute_gradu_z.inc"
!      include "compute_gradv_z.inc"
!      include "compute_gradw_z.inc"
!      include "compute_gradT_z.inc"
!      else
!      gradqp_z=0.0
!      end if
      select case (trim(turbulence))

        case ('none')
          !do nothing
          continue

        case ('sa', 'saBC')
          call compute_gradient_G(gradtv_x, tv, 'x')
          call compute_gradient_G(gradtv_y, tv, 'y')
          if(kmx>2)then
          call compute_gradient_G(gradtv_z, tv, 'z')
          end if

        case ('sst', 'sst2003')
          call compute_gradient_G(gradtk_x, tk, 'x')
          call compute_gradient_G(gradtw_x, tw, 'x')
          call compute_gradient_G(gradtk_y, tk, 'y')
          call compute_gradient_G(gradtw_y, tw, 'y')
          if(kmx>2)then
          call compute_gradient_G(gradtk_z, tk, 'z')
          call compute_gradient_G(gradtw_z, tw, 'z')
          end if

        case ('kkl')
          call compute_gradient_G(gradtk_x , tk , 'x')
          call compute_gradient_G(gradtkl_x, tkl, 'x')
          call compute_gradient_G(gradtk_y , tk , 'y')
          call compute_gradient_G(gradtkl_y, tkl, 'y')
          if(kmx>2)then
          call compute_gradient_G(gradtk_z , tk , 'z')
          call compute_gradient_G(gradtkl_z, tkl, 'z')
          end if

        case DEFAULT
          !call turbulence_read_error()
          Fatal_error

      end select


      select case(trim(transition))
        case('lctm2015')
          call compute_gradient_G(gradtgm_x, tgm, 'x')
          call compute_gradient_G(gradtgm_y, tgm, 'y')
          if(kmx>2)then
          call compute_gradient_G(gradtgm_z, tgm, 'z')
          end if

        case('bc', 'none')
          !do nothing
          continue

        case DEFAULT
          Fatal_error

      end Select


      !gradqp_z=0.0
      !applying boundary condition to gradients
      call apply_gradient_bc()

    end subroutine evaluate_all_gradients


    subroutine compute_gradient_G(grad, var, dir)
      !<  Compute gradient of any input scalar
      implicit none
      real, dimension( 0:imx  , 0:jmx  , 0:kmx  ), intent(out) :: grad
      !< Output variable storing the graident of var
      real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2), intent(in) :: var
      !< Input variable of which graident is required
      character(len=*)                           , intent(in) :: dir
      !< Direction with respect to which gradients are calculated
      
      real, dimension(:,:,:), pointer  :: nx
      real, dimension(:,:,:), pointer  :: ny
      real, dimension(:,:,:), pointer  :: nz

      integer :: i
      integer :: j
      integer :: k

      select case(dir)
        case('x')
          nx(-2:imx+3,-2:jmx+2,-2:kmx+2) => xn(:,:,:,1)
          ny(-2:imx+2,-2:jmx+3,-2:kmx+2) => yn(:,:,:,1)
          nz(-2:imx+2,-2:jmx+2,-2:kmx+3) => zn(:,:,:,1)
        case('y')
          nx(-2:imx+3,-2:jmx+2,-2:kmx+2) => xn(:,:,:,2)
          ny(-2:imx+2,-2:jmx+3,-2:kmx+2) => yn(:,:,:,2)
          nz(-2:imx+2,-2:jmx+2,-2:kmx+3) => zn(:,:,:,2)
        case('z')
          nx(-2:imx+3,-2:jmx+2,-2:kmx+2) => xn(:,:,:,3)
          ny(-2:imx+2,-2:jmx+3,-2:kmx+2) => yn(:,:,:,3)
          nz(-2:imx+2,-2:jmx+2,-2:kmx+3) => zn(:,:,:,3)
        case DEFAULT
          nx(-2:imx+3,-2:jmx+2,-2:kmx+2) => xn(:,:,:,1)
          ny(-2:imx+2,-2:jmx+3,-2:kmx+2) => yn(:,:,:,1)
          nz(-2:imx+2,-2:jmx+2,-2:kmx+3) => zn(:,:,:,1)
          print*, "ERROR: gradient direction error"
      end select
      grad = 0.0

      do k=0,kmx
        do j=0,jmx
          do i=0,imx
            grad(i,j,k) =(-(var(i-1,j  ,k  )+var(i,j,k))*nx(i,j,k)*xA(i,j,k) &
                          -(var(i  ,j-1,k  )+var(i,j,k))*ny(i,j,k)*yA(i,j,k) &
                          -(var(i  ,j  ,k-1)+var(i,j,k))*nz(i,j,k)*zA(i,j,k) &
                          +(var(i+1,j  ,k  )+var(i,j,k))*nx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                          +(var(i  ,j+1,k  )+var(i,j,k))*ny(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                          +(var(i  ,j  ,k+1)+var(i,j,k))*nz(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                         )/(2*volume(i,j,k))
          end do
        end do
      end do
      if(any(isnan(grad)))then
        Fatal_error
      end if

    end subroutine compute_gradient_G

    subroutine compute_gradient_T(grad, dir)
      !< Calculate gradient of temperature

      implicit none
      real, dimension( 0:imx  , 0:jmx  , 0:kmx  ), intent(out) :: grad
      !< Output gradient of termperature
      character(len=*)                           , intent(in) :: dir
      !< Direction with respect to which gradients are calculated
      
      real, dimension(6)               :: T
      real                             :: cell_T
      real, dimension(:,:,:), pointer  :: nx
      real, dimension(:,:,:), pointer  :: ny
      real, dimension(:,:,:), pointer  :: nz

      integer :: i
      integer :: j
      integer :: k

      select case(dir)
        case('x')
          nx(-2:imx+3,-2:jmx+2,-2:kmx+2) => xn(:,:,:,1)
          ny(-2:imx+2,-2:jmx+3,-2:kmx+2) => yn(:,:,:,1)
          nz(-2:imx+2,-2:jmx+2,-2:kmx+3) => zn(:,:,:,1)
        case('y')
          nx(-2:imx+3,-2:jmx+2,-2:kmx+2) => xn(:,:,:,2)
          ny(-2:imx+2,-2:jmx+3,-2:kmx+2) => yn(:,:,:,2)
          nz(-2:imx+2,-2:jmx+2,-2:kmx+3) => zn(:,:,:,2)
        case('z')
          nx(-2:imx+3,-2:jmx+2,-2:kmx+2) => xn(:,:,:,3)
          ny(-2:imx+2,-2:jmx+3,-2:kmx+2) => yn(:,:,:,3)
          nz(-2:imx+2,-2:jmx+2,-2:kmx+3) => zn(:,:,:,3)
        case DEFAULT
          nx(-2:imx+3,-2:jmx+2,-2:kmx+2) => xn(:,:,:,1)
          ny(-2:imx+2,-2:jmx+3,-2:kmx+2) => yn(:,:,:,1)
          nz(-2:imx+2,-2:jmx+2,-2:kmx+3) => zn(:,:,:,1)
          print*, "ERROR: gradient direction error"
      end select
      grad = 0.0

      do k=0,kmx
        do j=0,jmx
          do i=0,imx

            cell_T = (pressure(i,j,k)/density(i,j,k))/R_gas

            T(1)   = (pressure(i-1,j,k)/density(i-1,j,k))/R_gas + cell_T
            T(2)   = (pressure(i,j-1,k)/density(i,j-1,k))/R_gas + cell_T
            T(3)   = (pressure(i,j,k-1)/density(i,j,k-1))/R_gas + cell_T
            T(4)   = (pressure(i+1,j,k)/density(i+1,j,k))/R_gas + cell_T
            T(5)   = (pressure(i,j+1,k)/density(i,j+1,k))/R_gas + cell_T
            T(6)   = (pressure(i,j,k+1)/density(i,j,k+1))/R_gas + cell_T

            grad(i,j,k) =(-T(1)*nx(i,j,k)*xA(i,j,k) &
                          -T(2)*ny(i,j,k)*yA(i,j,k) &
                          -T(3)*nz(i,j,k)*zA(i,j,k) &
                          +T(4)*nx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                          +T(5)*ny(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                          +T(6)*nz(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                         )/(2*volume(i,j,k))
          end do
        end do
      end do
      if(any(isnan(grad)))then
        Fatal_error
      end if

    end subroutine compute_gradient_T


end module summon_grad_evaluation
