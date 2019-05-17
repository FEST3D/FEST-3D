  !< Calculate Cell-center and normal through them for transition model
module CC
  !< In order to calculate pressure gradient in the transition model, two
  !< quantities are required, the distance of the cell-center from the wall 
  !< the normal made the distance vector (from wall to cell-center).
  !< This module calucate both with gradient of V.n also.

#include "debug.h"
#include "error.h"

  use global_vars, only: imx
  use global_vars, only: jmx
  use global_vars, only: kmx
  use global_vars, only: process_id
  use global_vars, only: xn
  use global_vars, only: yn
  use global_vars, only: zn
  use global_vars, only: xA
  use global_vars, only: yA
  use global_vars, only: zA
  use global_vars, only: volume
  use global_vars, only: CCnormalX
  use global_vars, only: CCnormalY
  use global_vars, only: CCnormalZ
  use global_vars, only: CCVn
  use global_vars, only: DCCVnX
  use global_vars, only: DCCVnY
  use global_vars, only: DCCVnZ
  use global_vars, only: dist
  use global_vars, only: x_speed
  use global_vars, only: y_speed
  use global_vars, only: z_speed
  use global_vars, only: transition
  use global_vars, only: turbulence

  use utils, only : alloc
  use utils, only : dealloc
  implicit none
  private
  public :: find_DCCVn
  public :: setupCC
  public :: destroyCC
  
  contains

    subroutine setupCC()
      !< Allocate memory for the cell center variable only in case of transition model
      implicit none

      DebugCall("Setup CC")

      if((transition=='lctm2015') .and. turbulence/='none')then
        call alloc(CCnormalX, -2, imx+2, -2, jmx+2, -2, kmx+2, AErrMsg("CCnormalX"))
        call alloc(CCnormalY, -2, imx+2, -2, jmx+2, -2, kmx+2, AErrMsg("CCnormalY"))
        call alloc(CCnormalZ, -2, imx+2, -2, jmx+2, -2, kmx+2, AErrMsg("CCnormalZ"))
        call alloc(CCVn,   -2, imx+2, -2, jmx+2, -2, kmx+2, AErrMsg("CCVn"))
        call alloc(DCCVnX, -2, imx+2, -2, jmx+2, -2, kmx+2, AErrMsg("DCCVnZ"))
        call alloc(DCCVnY, -2, imx+2, -2, jmx+2, -2, kmx+2, AErrMsg("DCCVnY"))
        call alloc(DCCVnZ, -2, imx+2, -2, jmx+2, -2, kmx+2, AErrMsg("DCCVnZ"))
        call find_CCnormal()
      end if

    end subroutine setupCC


    subroutine destroyCC()
      !< Deallocate memory from the cell-center variables.
      implicit none

      DebugCall("Destroy CC")

      call dealloc(CCnormalX)
      call dealloc(CCnormalY)
      call dealloc(CCnormalZ)
      call dealloc(CCVn)
      call dealloc(DCCVnX)
      call dealloc(DCCVnY)
      call dealloc(DCCVnZ)

    end subroutine destroyCC


    subroutine find_CCnormal()
      !< Find the cell-center unit normal
      implicit none
      call compute_gradient(CCnormalX, dist, 'x')
      call compute_gradient(CCnormalY, dist, 'y')
      call compute_gradient(CCnormalZ, dist, 'z')
      !using already allocated memeory for storing magnitude
      CCVn = sqrt(CCnormalX**2 + CCnormalY**2 + CCnormalZ**2)
      !CCVn hold the magnitude of CCnormal temporaraly and can be 
      !overwritten after next three lines of code.
      CCnormalX = CCnormalX/(CCVn + 1e-12)
      CCnormalY = CCnormalY/(CCVn + 1e-12)
      CCnormalZ = CCnormalZ/(CCVn + 1e-12)
    end subroutine find_CCnormal


    subroutine find_CCVn()
      !< Taking a dot product between Cell-center velocity and unit normal
      implicit none
      CCVn = CCnormalX*x_speed + CCnormalY*y_speed + CCnormalZ*z_speed
    end subroutine find_CCVn


    subroutine find_DCCVn()
      !< Find gradient of the dot product between cell velocity and unit normal
      implicit none
      call find_CCVn()
      call compute_gradient(DCCVnX, dist, 'x')
      call compute_gradient(DCCVnY, dist, 'y')
      call compute_gradient(DCCVnZ, dist, 'z')
    end subroutine find_DCCVn


    subroutine compute_gradient(grad, var, dir)
      !< Generalized subroutine to calculate gradients
      implicit none
      real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2), intent(out) :: grad
      real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2), intent(in) :: var
      character(len=*)                           , intent(in) :: dir
      
      real, dimension(:,:,:), pointer  :: nx
      real, dimension(:,:,:), pointer  :: ny
      real, dimension(:,:,:), pointer  :: nz

      integer :: i
      integer :: j
      integer :: k

      ! initialize
      grad = 0.0

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

    end subroutine compute_gradient

end module CC
