module vortex_transport_state
  use global_vars, only: density
  use global_vars, only: x_speed
  use global_vars, only: y_speed
  use global_vars, only: z_speed
  use global_vars, only: pressure
  use global_vars, only : density_inf
  use global_vars, only : x_speed_inf
  use global_vars, only : y_speed_inf
  use global_vars, only : z_speed_inf
  use global_vars, only : pressure_inf
  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx
  use global_vars, only : R_gas
  use global_vars, only : gm
  use geometry   , only : CellCenter

  implicit none
  private

  public :: set_initial_vortex_state

  contains

    subroutine set_initial_vortex_state()
      implicit none
      real :: radius
      real :: Cp
      real :: Temp
      real, parameter  :: Radius_vortex=0.005
      real, parameter  :: Strength=1./5.
      real, parameter  :: xc=0.05
      real, parameter  :: yc=0.05
      real, parameter  :: Temp_inf = 300.
      real, parameter  :: Mach = 0.5

      integer :: i,j,k

      pressure_inf = 1.e+5
      density_inf  = pressure_inf/(R_gas*Temp_inf)
      x_speed_inf = Mach*sqrt(gm*R_gas*Temp_inf)
      Cp = R_gas*gm/(gm-1.)

      do k = -2,kmx+2
        do j = -2,jmx+2
          do i = -2,imx+2
            radius = sqrt((CellCenter(i,j,k,1)-xc)**2 + (CellCenter(i,j,k,2)-yc)**2)/Radius_vortex
            x_speed(i,j,k) = X_speed_inf*(1.-(Strength*(CellCenter(i,j,k,2)-yc)*exp(-0.5*radius**2)/Radius_vortex))
            y_speed(i,j,k) = X_speed_inf*Strength*(CellCenter(i,j,k,1)-xc)*exp(-0.5*radius**2)/Radius_vortex
            Temp = Temp_inf - 0.5*((X_speed_inf**2) * (Strength**2) * exp(-radius**2)/Cp)
            density(i,j,k) = density_inf*((Temp/Temp_inf)**(1./(gm-1.)))
            pressure(i,j,k) = density(i,j,k)*R_gas*Temp
          end do
        end do
      end do
      z_speed = 0.0


    end subroutine set_initial_vortex_state
end module vortex_transport_state
