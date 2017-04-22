module molecular_viscosity

  !------------------------------------------------------------------
  ! author    - Jatinder Pal Singh Sandhu
  ! objective - setup  (memory allocation and intialization),
  !             destroy (memory deallocation), and
  !             calcualte (using sutherlad law)
  ! molecular viscosity 
  !------------------------------------------------------------------

  use global_vars  , only : mu
  use global_vars , only : mu_ref
  use global_vars , only : Sutherland_temp
  use global_vars , only : T_ref
  use global_vars , only : R_gas
  use global_vars , only : mu_variation

  use global_vars , only : pressure
  use global_vars , only : density

  use global_vars , only : imx
  use global_vars , only : jmx
  use global_vars , only : kmx

  use utils       , only : dmsg
  use utils       , only :   alloc
  use utils       , only : dealloc

  implicit none
  private

  public :: setup_molecular_viscosity
  public :: destroy_molecular_viscosity
  public :: calculate_molecular_viscosity

  contains

    subroutine setup_molecular_viscosity()
      implicit none

      call alloc(mu, -2, imx+2, -2, jmx+2, -2, kmx+2)
      mu = mu_ref !intialize

    end subroutine setup_molecular_viscosity

    subroutine destroy_molecular_viscosity()
      implicit none

      call dealloc(mu)

    end subroutine destroy_molecular_viscosity

    subroutine calculate_molecular_viscosity()

      select case (trim(mu_variation))
        case ('sutherland_law')
          call apply_sutherland_law()

        case ('constant')
          !do nothing
          !mu will be equal to mu_ref

        case DEFAULT
          print*,"mu_variation not recognized:"
          print*, "   found '",trim(mu_variation),"'"
          print*, "accepted values: 1) sutherland_law"
          print*, "                 2) constant"
      end select

    end subroutine calculate_molecular_viscosity

    subroutine apply_sutherland_law()
      implicit none
      integer :: i,j,k
      real :: T,ST

      ST = Sutherland_temp

      do k = 1,kmx
        do j = 1,jmx
          do i = 1,imx
            T = pressure(i,j,k)/(density(i,j,k)*R_gas)
            mu(i,j,k) = mu_ref * ((T/T_ref)**(1.5)) *((T_ref + ST)/(T + ST))
          end do
        end do
      end do

    end subroutine apply_sutherland_law

end module molecular_viscosity
