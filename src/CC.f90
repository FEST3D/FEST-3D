  !< Calculate Cell-center and normal through them for transition model
module CC
  !< In order to calculate pressure gradient in the transition model, two
  !< quantities are required: the distance of the cell-center from the wall 
  !< andn the normal made by the distance vector field (from wall to cell-center).
  !< This module calucate both with gradient of V.n also.

#include "debug.h"
#include "error.h"

   use vartypes
  use wall_dist, only: dist

  use utils, only : alloc
  implicit none
  private
  real(wp), dimension(:, :, :), allocatable             :: CCnormalX 
   !< Cell-Center normal nx with respect to wall; used for transition model (pressure gradient calcualtion)
  real(wp), dimension(:, :, :), allocatable             :: CCnormalY
   !< Cell-Center normal ny with respect to wall; used for transiton model (pressure gradient calculation)
  real(wp), dimension(:, :, :), allocatable             :: CCnormalZ
   !< Cell-Center normal nz with respect to wall; used for transiton model (pressure gradient calculation)
  real(wp), dimension(:, :, :), allocatable             :: CCVn 
  !< Store value at Cell-Center of dot product between velocity vector and cell-center normal. {vec(Velocity).normal}
  real(wp), dimension(:, :, :), allocatable             :: DCCVnX
  !< Store Derivative of Cell-Center CCVn with respect to x
  real(wp), dimension(:, :, :), allocatable             :: DCCVnY
  !< Store Derivative of Cell-Center CCVn with respect to y
  real(wp), dimension(:, :, :), allocatable             :: DCCVnZ
  !< Store Derivative of Cell-Center CCVn with respect to z
  public :: find_DCCVn
  public :: setupCC
  public :: CCnormalX
  public :: CCnormalY
  public :: CCnormalZ

  public :: DCCVnX
  public :: DCCVnY
  public :: DCCVnZ
  
  contains

    subroutine setupCC(scheme, cells, Ifaces, Jfaces, Kfaces, dims)
      !< Allocate memory for the cell center variable only in case of transition model
      implicit none
      type(schemetype), intent(in) :: scheme
      type(extent), intent(in) :: dims
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal

      DebugCall("Setup CC")

      if((scheme%transition=='lctm2015') .and. scheme%turbulence/='none')then
        call alloc(CCnormalX, -2, dims%imx+2, -2, dims%jmx+2, -2, dims%kmx+2, AErrMsg("CCnormalX"))
        call alloc(CCnormalY, -2, dims%imx+2, -2, dims%jmx+2, -2, dims%kmx+2, AErrMsg("CCnormalY"))
        call alloc(CCnormalZ, -2, dims%imx+2, -2, dims%jmx+2, -2, dims%kmx+2, AErrMsg("CCnormalZ"))
        call alloc(CCVn,   -2, dims%imx+2, -2, dims%jmx+2, -2, dims%kmx+2, AErrMsg("CCVn"))
        call alloc(DCCVnX, -2, dims%imx+2, -2, dims%jmx+2, -2, dims%kmx+2, AErrMsg("DCCVnZ"))
        call alloc(DCCVnY, -2, dims%imx+2, -2, dims%jmx+2, -2, dims%kmx+2, AErrMsg("DCCVnY"))
        call alloc(DCCVnZ, -2, dims%imx+2, -2, dims%jmx+2, -2, dims%kmx+2, AErrMsg("DCCVnZ"))
        call find_CCnormal(cells, Ifaces, Jfaces, Kfaces, dims)
      end if

    end subroutine setupCC


    subroutine find_CCnormal(cells, Ifaces, Jfaces, Kfaces,dims)
      !< Find the cell-center unit normal
      implicit none
      type(extent), intent(in) :: dims
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      call compute_gradient(CCnormalX, dist, cells, Ifaces, Jfaces, Kfaces, 'x', dims)
      call compute_gradient(CCnormalY, dist, cells, Ifaces, Jfaces, Kfaces, 'y', dims)
      call compute_gradient(CCnormalZ, dist, cells, Ifaces, Jfaces, Kfaces, 'z', dims)
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
      real(wp), dimension(-2:dims%imx,-2:dims%jmx,-2:dims%kmx,-2:dims%n_var), intent(in) :: qp
      CCVn = CCnormalX*qp(:,:,:,2) + CCnormalY*qp(:,:,:,3) + CCnormalZ*qp(:,:,:,4) ! (nx,ny,nz).(u,v,w)
    end subroutine find_CCVn


    subroutine find_DCCVn(qp,cells, Ifaces, Jfaces, Kfaces,dims)
      !< Find gradient of the dot product between cell velocity and unit normal
      implicit none
      type(extent), intent(in) :: dims
      real(wp), dimension(-2:dims%imx,-2:dims%jmx,-2:dims%kmx,-2:dims%n_var), intent(in) :: qp
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      call find_CCVn(qp, dims)
      call compute_gradient(DCCVnX, dist, cells, Ifaces, Jfaces, Kfaces, 'x', dims)
      call compute_gradient(DCCVnY, dist, cells, Ifaces, Jfaces, Kfaces, 'y', dims)
      call compute_gradient(DCCVnZ, dist, cells, Ifaces, Jfaces, Kfaces, 'z', dims)
    end subroutine find_DCCVn


    subroutine compute_gradient(grad, var, cells, Ifaces, Jfaces, Kfaces, dir, dims)
      !< Generalized subroutine to calculate gradients
      implicit none
      type(extent), intent(in) :: dims
      real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: grad
      real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: var
      character(len=*)                                          , intent(in) :: dir
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      
      integer :: i
      integer :: j
      integer :: k

      ! initialize
      grad = 0.0

      select case(dir)
        case('x')
          do k=0,dims%kmx
            do j=0,dims%jmx
              do i=0,dims%imx
                grad(i,j,k) =(-(var(i-1,j  ,k  )+var(i,j,k))*Ifaces(i,j,k)%nx*Jfaces(i,j,k)%A &
                              -(var(i  ,j-1,k  )+var(i,j,k))*Ifaces(i,j,k)%ny*Jfaces(i,j,k)%A &
                              -(var(i  ,j  ,k-1)+var(i,j,k))*Ifaces(i,j,k)%nz*Jfaces(i,j,k)%A &
                              +(var(i+1,j  ,k  )+var(i,j,k))*Ifaces(i+1,j  ,k  )%nx*Ifaces(i+1,j  ,k  )%A &
                              +(var(i  ,j+1,k  )+var(i,j,k))*Ifaces(i  ,j+1,k  )%ny*Ifaces(i  ,j+1,k  )%A &
                              +(var(i  ,j  ,k+1)+var(i,j,k))*Ifaces(i  ,j  ,k+1)%nz*Ifaces(i  ,j  ,k+1)%A &
                             )/(2*cells(i,j,k)%volume)
              end do
            end do
          end do
        case('y')
          do k=0,dims%kmx
            do j=0,dims%jmx
              do i=0,dims%imx
                grad(i,j,k) =(-(var(i-1,j  ,k  )+var(i,j,k))*Jfaces(i,j,k)%nx*Jfaces(i,j,k)%A &
                              -(var(i  ,j-1,k  )+var(i,j,k))*Jfaces(i,j,k)%ny*Jfaces(i,j,k)%A &
                              -(var(i  ,j  ,k-1)+var(i,j,k))*Jfaces(i,j,k)%nz*Jfaces(i,j,k)%A &
                              +(var(i+1,j  ,k  )+var(i,j,k))*Jfaces(i+1,j  ,k  )%nx*Jfaces(i+1,j  ,k  )%A &
                              +(var(i  ,j+1,k  )+var(i,j,k))*Jfaces(i  ,j+1,k  )%ny*Jfaces(i  ,j+1,k  )%A &
                              +(var(i  ,j  ,k+1)+var(i,j,k))*Jfaces(i  ,j  ,k+1)%nz*Jfaces(i  ,j  ,k+1)%A &
                             )/(2*cells(i,j,k)%volume)
              end do
            end do
          end do
        case('z')
          do k=0,dims%kmx
            do j=0,dims%jmx
              do i=0,dims%imx
                grad(i,j,k) =(-(var(i-1,j  ,k  )+var(i,j,k))*Kfaces(i,j,k)%nx*Kfaces(i,j,k)%A &
                              -(var(i  ,j-1,k  )+var(i,j,k))*Kfaces(i,j,k)%ny*Kfaces(i,j,k)%A &
                              -(var(i  ,j  ,k-1)+var(i,j,k))*Kfaces(i,j,k)%nz*Kfaces(i,j,k)%A &
                              +(var(i+1,j  ,k  )+var(i,j,k))*Kfaces(i+1,j  ,k  )%nx*Kfaces(i+1,j  ,k  )%A &
                              +(var(i  ,j+1,k  )+var(i,j,k))*Kfaces(i  ,j+1,k  )%ny*Kfaces(i  ,j+1,k  )%A &
                              +(var(i  ,j  ,k+1)+var(i,j,k))*Kfaces(i  ,j  ,k+1)%nz*Kfaces(i  ,j  ,k+1)%A &
                             )/(2*cells(i,j,k)%volume)
              end do
            end do
          end do
        case DEFAULT
          print*, "ERROR: gradient direction error"
          Fatal_error
      end select
      if(any(isnan(grad)))then
        Fatal_error
      end if

    end subroutine compute_gradient

end module CC
