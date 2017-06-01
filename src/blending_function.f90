module blending_function

  use global_sst , only : sst_F1
  use global_sst , only : sigma_w2
  use global_sst , only : bstar

  use global_vars , only : mu
  use global_vars , only : tk
  use global_vars , only : tw
  use global_vars , only : density
  use global_vars , only : dist
  
  ! gradients

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

  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx
  use global_vars, only : turbulence

  use utils

  implicit none
  private

  public :: setup_sst_F1
  public :: destroy_sst_F1
  public :: calculate_sst_F1

  contains

    subroutine setup_sst_F1()
      implicit none
      call alloc(sst_F1, -2,imx+2, -2,jmx+2, -2,kmx+2)
      sst_F1=0.
    end subroutine setup_sst_F1

    subroutine destroy_sst_F1()
      implicit none
      call dealloc(sst_F1)
    end subroutine destroy_sst_F1

    subroutine calculate
      implicit none
      integer :: i,j,k
      real :: arg1
      real :: CD
      real :: var1
      real :: var2
      real :: right
      real :: left

      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1

            CD = max(2*density(i,j,k)*sigma_w2*(                             & 
                                                gradtk_x(i,j,k)*gradtw_x(i,j,k)&
                                              + gradtk_y(i,j,k)*gradtw_y(i,j,k)&
                                              + gradtk_z(i,j,k)*gradtw_z(i,j,k)&
                                               )/tw(i,j,k),                  &
                     1e-20)

            var1 = sqrt(tk(i,j,k))/(bstar*tw(i,j,k)*dist(i,j,k))
            var2 = 500*(mu(i,j,k)/density(i,j,k))/((dist(i,j,k)**2)*tw(i,j,k))
            right = 4*(density(i,j,k)*sigma_w2*tk(i,j,k))/(CD*(dist(i,j,k)**2))
            left = max(var1, var2)
            arg1 = min(left, right)
            sst_F1(i,j,k) = tanh(arg1**4)

          end do
        end do
      end do


    end subroutine calculate

    subroutine calculate_sst_F1
      implicit none

      select case (turbulence)

        case ('none')
          !do nothing
          continue

        case ('sst')
          call calculate()

        case Default
          call turbulence_read_error()

      end select
  

    end subroutine calculate_sst_F1

end module blending_function

