module dissipation
  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx
  use global_vars, only : pressure
  use global_vars, only : density
  use global_vars, only : x_speed
  use global_vars, only : y_speed
  use global_vars, only : z_speed
  use global_vars, only : qp
  use global_vars, only : Diss
  use global_vars, only : Volume
  use global_vars, only : n_var
  use global_vars, only : process_id
  use utils
  
#include "error.inc"
  implicit none
  private

  real, dimension(:,:,:), allocatable :: epsi2
  real, dimension(:,:,:), allocatable :: epsj2
  real, dimension(:,:,:), allocatable :: epsk2
  real, parameter :: k2=0.25
  real, parameter :: k4 =1./256.

  public :: setup_dissipation
  public :: destroy_dissipation
  public :: get_dissipation

  contains

    subroutine setup_dissipation()
      implicit none
      call alloc(epsi2,1,imx,1,jmx-1,1,kmx-1)
      call alloc(epsj2,1,imx-1,1,jmx,1,kmx-1)
      call alloc(epsk2,1,imx-1,1,jmx-1,1,kmx)
      call alloc(Diss,1,imx-1,1,jmx-1,1,kmx-1,1,n_var)
    end subroutine setup_dissipation

    subroutine destroy_dissipation()
      implicit none
      call dealloc(epsi2)
      call dealloc(epsj2)
      call dealloc(epsk2)
      call dealloc(Diss)
    end subroutine destroy_dissipation

    subroutine get_dissipation()
      implicit none
      integer :: i,j,k,l
      real :: mu1, mu2
      real :: ei2p, ej2p, ek2p
      real :: ei2m, ej2m, ek2m
      real :: ei4p, ej4p, ek4p
      real :: ei4m, ej4m, ek4m

      ! loop through I faces
      do k=1,kmx-1
        do j=1,jmx-1
          do i=1,imx
            mu1 = abs(pressure(i+1,j,k)-2*pressure(i,j,k)+pressure(i-1,j,k))
            mu1 = mu1/(abs(pressure(i+1,j,k))+abs(2*pressure(i,j,k))+abs(pressure(i-1,j,k)))

            mu2 = abs(pressure(i,j,k)-2*pressure(i-1,j,k)+pressure(i-2,j,k))
            mu2 = mu2/(abs(pressure(i,j,k))+abs(2*pressure(i-1,j,k))+abs(pressure(i-2,j,k)))

            epsi2(i,j,k) = k2*max(mu1,mu2)
          end do
        end do
      end do

      ! loop through J faces
      do k=1,kmx-1
        do j=1,jmx
          do i=1,imx-1
            mu1 = abs(pressure(i,j+1,k)-2*pressure(i,j,k)+pressure(i,j-1,k))
            mu1 = mu1/(abs(pressure(i,j+1,k))+abs(2*pressure(i,j,k))+abs(pressure(i,j-1,k)))

            mu2 = abs(pressure(i,j,k)-2*pressure(i,j-1,k)+pressure(i,j-2,k))
            mu2 = mu2/(abs(pressure(i,j,k))+abs(2*pressure(i,j-1,k))+abs(pressure(i,j-2,k)))

            epsj2(i,j,k) = k2*max(mu1,mu2)
          end do
        end do
      end do

      ! loop through K faces
      do k=1,kmx
        do j=1,jmx-1
          do i=1,imx-1
            mu1 = abs(pressure(i,j,k+1)-2*pressure(i,j,k)+pressure(i,j,k-1))
            mu1 = mu1/(abs(pressure(i,j,k+1))+abs(2*pressure(i,j,k))+abs(pressure(i,j,k-1)))

            mu2 = abs(pressure(i,j,k)-2*pressure(i,j,k-1)+pressure(i,j,k-2))
            mu2 = mu2/(abs(pressure(i,j,k))+abs(2*pressure(i,j,k-1))+abs(pressure(i,j,k-2)))

            epsk2(i,j,k) = k2*max(mu1,mu2)
          end do
        end do
      end do

      ! loop through all the cell to add dissipation to residue array
      do l=1,n_var
        do k=1,kmx-1
          do j=1,jmx-1
            do i=1,imx-1
              ei2p=epsi2(i+1,j,k)
              ei2m=epsi2(i,j,k)
              ej2p=epsj2(i,j+1,k)
              ej2m=epsj2(i,j,k)
              ek2p=epsk2(i,j,k+1)
              ek2m=epsk2(i,j,k)
              ei4p=max(0.0,(k4-ei2p))
              ei4m=max(0.0,(k4-ei2m))
              ej4p=max(0.0,(k4-ej2p))
              ej4m=max(0.0,(k4-ej2m))
              ek4p=max(0.0,(k4-ek2p))
              ek4m=max(0.0,(k4-ek2m))
              Diss(i,j,k,l)=ei2p*(qp(i+1,j,k,l)-qp(i,j,k,l))&
                           +ej2p*(qp(i,j+1,k,l)-qp(i,j,k,l))&
                           +ek2p*(qp(i,j,k+1,l)-qp(i,j,k,l))&
                           -ei2m*(qp(i,j,k,l)-qp(i-1,j,k,l))&
                           -ej2m*(qp(i,j,k,l)-qp(i,j-1,k,l))&
                           -ek2m*(qp(i,j,k,l)-qp(i,j,k-1,l))&
                           -ei4p*(qp(i+2,j,k,l)-3*qp(i+1,j,k,l)+3*qp(i,j,k,l)-qp(i-1,j,k,l))&
                           -ej4p*(qp(i,j+2,k,l)-3*qp(i,j+1,k,l)+3*qp(i,j,k,l)-qp(i,j-1,k,l))&
                           -ek4p*(qp(i,j,k+2,l)-3*qp(i,j,k+1,l)+3*qp(i,j,k,l)-qp(i,j,k-1,l))&
                           +ei4m*(qp(i+1,j,k,l)-3*qp(i,j,k,l)+3*qp(i-1,j,k,l)-qp(i-2,j,k,l))&
                           +ej4m*(qp(i,j+1,k,l)-3*qp(i,j,k,l)+3*qp(i,j-1,k,l)-qp(i,j-2,k,l))&
                           +ek4m*(qp(i,j,k+1,l)-3*qp(i,j,k,l)+3*qp(i,j,k-1,l)-qp(i,j,k-2,l))
              Diss(i,j,k,l)=Diss(i,j,k,l)*Volume(i,j,k)
              if(isnan(Diss(i,j,k,l)))then
                Fatal_error
              end if
            end do
          end do
        end do
      end do
    end subroutine get_dissipation

end module dissipation
