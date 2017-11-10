module summon_grad_evaluation
  !----------------------------------------------------------
  !170608  -Jatinder Pal Singh Sandhu
  ! Aim : call is made to all the required gradients
  !                         based on input conditions
  !----------------------------------------------------------
#include "error.inc"
  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx
  use global_vars, only : n_var
  use global_vars, only : n_grad
  use global_vars, only : qp
  use global_vars, only : x_speed
  use global_vars, only : y_speed
  use global_vars, only : z_speed
  use global_vars, only : tk
  use global_vars, only : tw
  use global_vars, only : tv
  use global_vars, only : te
  use global_vars, only : tkl
  use global_vars, only : turbulence
  use global_vars, only : gradqp_x
  use global_vars, only : gradqp_y
  use global_vars, only : gradqp_z
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

  public :: evaluate_all_gradients

  contains

    subroutine evaluate_all_gradients()

      implicit none
      integer            :: i,j,k,l,f
      real, dimension(6) :: area=1.0
      real, dimension(6) :: qface=0.0
      real, dimension(6) :: x_factor=0.0
      real, dimension(6) :: y_factor=0.0
      real, dimension(6) :: z_factor=0.0
      real               :: prefix=1.0
      real               :: qcurrent
      real               :: GradqpSumX
      real               :: GradqpSumY
      real               :: GradqpSumZ

      call dmsg(1, 'summon_grad_evaluation', 'evaluate_all_gradients')

      do k=0,kmx
        do j=0,jmx
          do i=0,imx
            !prefix (2*volume of cell)
            prefix = 1.0/(2*volume(i,j,k))
            ! area of all the faces of cell
            area(1) = xA(i  ,j  ,k  )
            area(4) = xA(i+1,j  ,k  )
            area(2) = yA(i  ,j  ,k  )
            area(5) = yA(i  ,j+1,k  )
            area(3) = zA(i  ,j  ,k  )
            area(6) = zA(i  ,j  ,k+1)

            !x direction loop costant
            x_factor(1)=-xn(i  ,j  ,k  ,1)*area(1)
            x_factor(4)=+xn(i+1,j  ,k  ,1)*area(4)
            x_factor(2)=-yn(i  ,j  ,k  ,1)*area(2)
            x_factor(5)=+yn(i  ,j+1,k  ,1)*area(5)
            x_factor(3)=-zn(i  ,j  ,k  ,1)*area(3)
            x_factor(6)=+zn(i  ,j  ,k+1,1)*area(6)

            !y direction loop costant
            y_factor(1)=-xn(i  ,j  ,k  ,2)*area(1)
            y_factor(4)=+xn(i+1,j  ,k  ,2)*area(4)
            y_factor(2)=-yn(i  ,j  ,k  ,2)*area(2)
            y_factor(5)=+yn(i  ,j+1,k  ,2)*area(5)
            y_factor(3)=-zn(i  ,j  ,k  ,2)*area(3)
            y_factor(6)=+zn(i  ,j  ,k+1,2)*area(6)

            !z direction loop costant
            z_factor(1)=-xn(i  ,j  ,k  ,3)*area(1)
            z_factor(4)=+xn(i+1,j  ,k  ,3)*area(4)
            z_factor(2)=-yn(i  ,j  ,k  ,3)*area(2)
            z_factor(5)=+yn(i  ,j+1,k  ,3)*area(5)
            z_factor(3)=-zn(i  ,j  ,k  ,3)*area(3)
            z_factor(6)=+zn(i  ,j  ,k+1,3)*area(6)

            do l=2,4 ! all component of velocity
              qcurrent = qp(i,j,k,l)
              qface(1)=qp(i-1,j  ,k  ,l)+qcurrent
              qface(4)=qp(i+1,j  ,k  ,l)+qcurrent
              qface(2)=qp(i  ,j-1,k  ,l)+qcurrent
              qface(5)=qp(i  ,j+1,k  ,l)+qcurrent
              qface(3)=qp(i  ,j  ,k-1,l)+qcurrent
              qface(6)=qp(i  ,j  ,k+1,l)+qcurrent
              GradqpSumX=0.0
              GradqpSumY=0.0
              GradqpSumZ=0.0
              do f=1,6
                GradqpSumX=GradqpSumX+(qface(f)*x_factor(f))
                GradqpSumY=GradqpSumY+(qface(f)*y_factor(f))
                GradqpSumZ=GradqpSumZ+(qface(f)*z_factor(f))
              end do
              gradqp_x(i,j,k,l-1)=GradqpSumX
              gradqp_y(i,j,k,l-1)=GradqpSumY
              gradqp_z(i,j,k,l-1)=GradqpSumZ
            end do

            do l=5,5 ! all component of temperature
              qcurrent = (qp(i,j,k,l)/qp(i,j,k,1))
              qface(1)=((qp(i-1,j  ,k  ,l)/qp(i-1,j  ,k  ,1))+qcurrent)/R_gas
              qface(4)=((qp(i+1,j  ,k  ,l)/qp(i+1,j  ,k  ,1))+qcurrent)/R_gas
              qface(2)=((qp(i  ,j-1,k  ,l)/qp(i  ,j-1,k  ,1))+qcurrent)/R_gas
              qface(5)=((qp(i  ,j+1,k  ,l)/qp(i  ,j+1,k  ,1))+qcurrent)/R_gas
              qface(3)=((qp(i  ,j  ,k-1,l)/qp(i  ,j  ,k-1,1))+qcurrent)/R_gas
              qface(6)=((qp(i  ,j  ,k+1,l)/qp(i  ,j  ,k+1,1))+qcurrent)/R_gas
              GradqpSumX=0.0
              GradqpSumY=0.0
              GradqpSumZ=0.0
              do f=1,6
                GradqpSumX=GradqpSumX+(qface(f)*x_factor(f))
                GradqpSumY=GradqpSumY+(qface(f)*y_factor(f))
                GradqpSumZ=GradqpSumZ+(qface(f)*z_factor(f))
              end do
              gradqp_x(i,j,k,l-1)=GradqpSumX
              gradqp_y(i,j,k,l-1)=GradqpSumY
              gradqp_z(i,j,k,l-1)=GradqpSumZ
            end do

            do l=6,n_var ! all turbulent variables
              qcurrent = qp(i,j,k,l)
              qface(1)=qp(i-1,j  ,k  ,l)+qcurrent
              qface(4)=qp(i+1,j  ,k  ,l)+qcurrent
              qface(2)=qp(i  ,j-1,k  ,l)+qcurrent
              qface(5)=qp(i  ,j+1,k  ,l)+qcurrent
              qface(3)=qp(i  ,j  ,k-1,l)+qcurrent
              qface(6)=qp(i  ,j  ,k+1,l)+qcurrent
              GradqpSumX=0.0
              GradqpSumY=0.0
              GradqpSumZ=0.0
              do f=1,6
                GradqpSumX=GradqpSumX+(qface(f)*x_factor(f))
                GradqpSumY=GradqpSumY+(qface(f)*y_factor(f))
                GradqpSumZ=GradqpSumZ+(qface(f)*z_factor(f))
              end do
              gradqp_x(i,j,k,l-1)=GradqpSumX
              gradqp_y(i,j,k,l-1)=GradqpSumY
              gradqp_z(i,j,k,l-1)=GradqpSumZ
            end do
          end do
        end do
      end do

      call apply_gradient_bc()

    end subroutine evaluate_all_gradients


end module summon_grad_evaluation
