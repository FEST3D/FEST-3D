module residual_smoothing

  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx
  use global_vars, only : residue
  use global_vars, only : eps_x
  use global_vars, only : eps_y
  use global_vars, only : eps_z
  use global_vars, only : process_id
  use utils

#include "error.inc"
  implicit none
  private

  ! LU[Rb]=Residue
  ! L[Rt] = Residue
  ! U[Rb] = Rt
  real, dimension(:,:), allocatable :: Rb
  real, dimension(:,:), allocatable :: Rt
  real, dimension(:), allocatable :: Ax
  real, dimension(:), allocatable :: Bx
  real, dimension(:), allocatable :: Cx
  real, dimension(:), allocatable :: Ay
  real, dimension(:), allocatable :: By
  real, dimension(:), allocatable :: Cy
  real, dimension(:), allocatable :: Az
  real, dimension(:), allocatable :: Bz
  real, dimension(:), allocatable :: Cz
  real  :: t1
  public :: get_smoothen_residue
  public :: setup_implicit_residual_smoothing
  public :: destroy_implicit_residual_smoothing

  contains

    subroutine setup_implicit_residual_smoothing()
      implicit none
      integer :: i,j,k
      eps_x=0.4
      eps_y=0.4
      eps_z=0.4

      !I-direction
      call alloc(Ax,2,imx-1)
      call alloc(Bx,1,imx-1)
      call alloc(Cx,1,imx-2)
      Cx=-eps_x
      Bx(1) = 1+2*eps_x
      do i = 2,imx-1
        Ax(i)=-eps_x/Bx(i-1)
        Bx(i)=1+2*eps_x - Ax(i)*Cx(i-1)
      end do
      if(any(Bx==0.))then
        Fatal_error
      end if


      !J-direction
      call alloc(Ay,2,jmx-1)
      call alloc(By,1,jmx-1)
      call alloc(Cy,1,jmx-2)
      Cy=-eps_y
      By(1) = 1+2*eps_y
      do j = 2,jmx-1
        Ay(j)=-eps_y/By(j-1)
        By(j)=1+2*eps_y - Ay(j)*Cy(j-1)
      end do
      if(any(By==0.))then
        Fatal_error
      end if


      if(kmx>2)then
      !I-direction
      call alloc(Az,2,kmx-1)
      call alloc(Bz,1,kmx-1)
      call alloc(Cz,1,kmx-2)
      Cz=-eps_z
      Bz(1) = 1+2*eps_z
      do k = 2,kmx-1
        Az(k)=-eps_z/Bz(k-1)
        Bz(k)=1+2*eps_z - Az(k)*Cz(k-1)
      end do
      end if
    end subroutine setup_implicit_residual_smoothing



    subroutine destroy_implicit_residual_smoothing()
      implicit none

      !I-direction
      call dealloc(Ax)
      call dealloc(Bx)
      call dealloc(Cx)

      !J-direction
      call dealloc(Ay)
      call dealloc(By)
      call dealloc(Cy)

      if(kmx>2)then
      !I-direction
      call dealloc(Az)
      call dealloc(Bz)
      call dealloc(Cz)
      end if
    end subroutine destroy_implicit_residual_smoothing



    subroutine get_smoothen_residue
      implicit none
      integer :: i,j,k

      ! I-implicit smoothing
      ! forward sweep
      !residue(1,:,:) = residue(1,:,:)
      do i=2,imx-1
        residue(i,:,:,:)=residue(i,:,:,:)-Ax(i)*residue(i-1,:,:,:)
        if(isnan(residue(i,1,1,1))) then
          Fatal_error
        end if
      end do

      ! backward subsitution
      residue(imx-1,:,:,:)=residue(imx-1,:,:,:)/Bx(imx-1)
      do i=imx-2,1,-1
        residue(i,:,:,:)=(residue(i,:,:,:)-Cx(i)*residue(i+1,:,:,:))/Bx(i)
        if(isnan(residue(i,1,1,1))) then
          Fatal_error
        end if
      end do

      ! J-implicit smoothing
      ! forward sweep
      !residue(1,:,:) = residue(1,:,:)
      do j=2,jmx-1
        residue(:,j,:,:)=residue(:,j,:,:)-Ay(j)*residue(:,j-1,:,:)
        if(isnan(residue(1,j,1,1))) then
          Fatal_error
        end if
      end do

      ! backward subsitution
      residue(:,jmx-1,:,:)=residue(:,jmx-1,:,:)/By(jmx-1)
      do j=jmx-2,1,-1
        residue(:,j,:,:)=(residue(:,j,:,:)-Cy(j)*residue(:,j+1,:,:))/By(j)
        if(isnan(residue(1,j,1,1))) then
          print*, "BY", By
          print*, "CY", Cy
          print*, "Bx", Bx
          print*, "Cx", Cx
          !print*, residue(:,j,:,:)
          Fatal_error
        end if
      end do




      ! Jacobi iteration from fluent
      ! 3D
!      do i=1,2
!      residue(2:imx-2,2:jmx-2,2:kmx-2,:) = (residue(2:imx-2,2:jmx-2,2:kmx-2,:)&
!        +0.5*(residue(1:imx-3,2:jmx-2,2:kmx-2,:)&
!             +residue(3:imx-1,2:jmx-2,2:kmx-2,:)&
!             +residue(2:imx-2,1:jmx-3,2:kmx-2,:)&
!             +residue(2:imx-2,3:jmx-1,2:kmx-2,:)&
!             +residue(2:imx-2,2:jmx-2,1:kmx-3,:)&
!             +residue(2:imx-2,2:jmx-2,3:kmx-1,:)&
!             ))/4.
!      end do
      !2D
!      do i=1,2
!      residue(2:imx-2,2:jmx-2,2:kmx-2,:) = (residue(2:imx-2,2:jmx-2,2:kmx-2,:)&
!        +0.5*(residue(1:imx-3,2:jmx-2,:,:)&
!             +residue(3:imx-1,2:jmx-2,:,:)&
!             +residue(2:imx-2,1:jmx-3,:,:)&
!             +residue(2:imx-2,3:jmx-1,:,:)&
!             ))/4.
!      end do
    end subroutine get_smoothen_residue


end module residual_smoothing
