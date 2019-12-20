  !< Calculate Cell-center and normal through them for transition model
module CC
  !< In order to calculate pressure gradient in the transition model, two
  !< quantities are required: the distance of the cell-center from the wall 
  !< andn the normal made by the distance vector field (from wall to cell-center).
  !< This module calucate both with gradient of V.n also.

#include "debug.h"
#include "error.h"

   use vartypes
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

  use utils, only : alloc
  implicit none
  private
  public :: find_DCCVn
  public :: setupCC
!  public :: destroyCC
  
  contains

    subroutine setupCC(scheme, dims)
      !< Allocate memory for the cell center variable only in case of transition model
      implicit none
      type(schemetype), intent(in) :: scheme
      type(extent), intent(in) :: dims

      DebugCall("Setup CC")

      if((scheme%transition=='lctm2015') .and. scheme%turbulence/='none')then
        call alloc(CCnormalX, -2, dims%imx+2, -2, dims%jmx+2, -2, dims%kmx+2, AErrMsg("CCnormalX"))
        call alloc(CCnormalY, -2, dims%imx+2, -2, dims%jmx+2, -2, dims%kmx+2, AErrMsg("CCnormalY"))
        call alloc(CCnormalZ, -2, dims%imx+2, -2, dims%jmx+2, -2, dims%kmx+2, AErrMsg("CCnormalZ"))
        call alloc(CCVn,   -2, dims%imx+2, -2, dims%jmx+2, -2, dims%kmx+2, AErrMsg("CCVn"))
        call alloc(DCCVnX, -2, dims%imx+2, -2, dims%jmx+2, -2, dims%kmx+2, AErrMsg("DCCVnZ"))
        call alloc(DCCVnY, -2, dims%imx+2, -2, dims%jmx+2, -2, dims%kmx+2, AErrMsg("DCCVnY"))
        call alloc(DCCVnZ, -2, dims%imx+2, -2, dims%jmx+2, -2, dims%kmx+2, AErrMsg("DCCVnZ"))
        call find_CCnormal(dims)
      end if

    end subroutine setupCC


!    subroutine destroyCC()
!      !< Deallocate memory from the cell-center variables
!      implicit none
!
!      DebugCall("Destroy CC")
!
!      call dealloc(CCnormalX)
!      call dealloc(CCnormalY)
!      call dealloc(CCnormalZ)
!      call dealloc(CCVn)
!      call dealloc(DCCVnX)
!      call dealloc(DCCVnY)
!      call dealloc(DCCVnZ)
!
!    end subroutine destroyCC


    subroutine find_CCnormal(dims)
      !< Find the cell-center unit normal
      implicit none
      type(extent), intent(in) :: dims
      call compute_gradient(CCnormalX, dist, 'x', dims)
      call compute_gradient(CCnormalY, dist, 'y', dims)
      call compute_gradient(CCnormalZ, dist, 'z', dims)
      !using already allocated memeory for storing magnitude
      CCVn = sqrt(CCnormalX**2 + CCnormalY**2 + CCnormalZ**2)
      !CCVn hold the magnitude of CCnormal temporaraly and can be 
      !overwritten after next three lines of code.
      CCnormalX = CCnormalX/(CCVn + 1e-12)
      CCnormalY = CCnormalY/(CCVn + 1e-12)
      CCnormalZ = CCnormalZ/(CCVn + 1e-12)
    end subroutine find_CCnormal


    subroutine find_CCVn(qp, dims)
      !< Taking a dot product between Cell-center velocity and unit normal
      implicit none
      type(extent), intent(in) :: dims
      real, dimension(-2:dims%imx,-2:dims%jmx,-2:dims%kmx,-2:dims%n_var), intent(in) :: qp
      CCVn = CCnormalX*qp(:,:,:,2) + CCnormalY*qp(:,:,:,3) + CCnormalZ*qp(:,:,:,4) ! (nx,ny,nz).(u,v,w)
    end subroutine find_CCVn


    subroutine find_DCCVn(qp, dims)
      !< Find gradient of the dot product between cell velocity and unit normal
      implicit none
      type(extent), intent(in) :: dims
      real, dimension(-2:dims%imx,-2:dims%jmx,-2:dims%kmx,-2:dims%n_var), intent(in) :: qp
      call find_CCVn(qp, dims)
      call compute_gradient(DCCVnX, dist, 'x', dims)
      call compute_gradient(DCCVnY, dist, 'y', dims)
      call compute_gradient(DCCVnZ, dist, 'z', dims)
    end subroutine find_DCCVn


    subroutine compute_gradient(grad, var, dir, dims)
      !< Generalized subroutine to calculate gradients
      implicit none
      type(extent), intent(in) :: dims
      real, dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: grad
      real, dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: var
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
          nx(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2) => xn(:,:,:,1)
          ny(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2) => yn(:,:,:,1)
          nz(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3) => zn(:,:,:,1)
        case('y')
          nx(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2) => xn(:,:,:,2)
          ny(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2) => yn(:,:,:,2)
          nz(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3) => zn(:,:,:,2)
        case('z')
          nx(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2) => xn(:,:,:,3)
          ny(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2) => yn(:,:,:,3)
          nz(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3) => zn(:,:,:,3)
        case DEFAULT
          nx(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2) => xn(:,:,:,1)
          ny(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2) => yn(:,:,:,1)
          nz(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3) => zn(:,:,:,1)
          print*, "ERROR: gradient direction error"
      end select
      grad = 0.0

      do k=0,dims%kmx
        do j=0,dims%jmx
          do i=0,dims%imx
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
