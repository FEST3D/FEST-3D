module turbulent_source


  use global_vars , only : turbulence

  use sst_grad, only : setup_sst_grad
  use sst_grad, only : destroy_sst_grad
  use sst_grad, only : calculate_sst_grad
  use sst_source , only : add_sst_source
  use utils

  implicit none
  private

  public :: setup_turbulent_grad
  public :: destroy_turbulent_grad
  public :: compute_turbulent_grad
  
  contains

    subroutine setup_turbulent_grad()
      implicit none
      select case (turbulence)
        
        case ('none')
          !do nothing
          continue

        case ('sst')
          call setup_sst_grad()

        case DEFAULT
          call turbulence_read_error()

      end select

    end subroutine setup_turbulent_grad

    subroutine destroy_turbulent_grad()
      implicit none
      select case (turbulence)
        
        case ('none')
          !do nothing
          continue

        case ('sst')
          call destroy_sst_grad()

        case DEFAULT
          call turbulence_read_error()

      end select

    end subroutine destroy_turbulent_grad

    subroutine compute_turbulent_grad()
      implicit none
      select case (turbulence)
        
        case ('none')
          !do nothing
          continue

        case ('sst')
          call calculate_sst_grad()
          call add_sst_source()

        case DEFAULT
          call turbulence_read_error()

      end select

    end subroutine compute_turbulent_grad

end module turbulent_source
