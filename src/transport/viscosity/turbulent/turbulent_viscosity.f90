module turbulent_viscosity

  !------------------------------------------------------------------
  ! author    - Jatinder Pal Singh Sandhu
  ! objective - setup   (memory allocation, link pointer)
  !             destroy (deallocate memory, nullify pointer)
  !             calcualte (call only)
  ! turbulent viscosity.
  !------------------------------------------------------------------
#include "../../../error.inc"
  use global_vars    , only : imx
  use global_vars    , only : jmx
  use global_vars    , only : kmx
  use global_vars    , only : mu_t
  use global_vars    , only : sst_mu
  use global_vars    , only : kkl_mu
  use global_vars    , only : turbulence
  use sst_viscosity  , only : calculate_sst_mu
  use kkl_viscosity  , only : calculate_kkl_mu

  use utils

  implicit none

  private

  public :: setup_turbulent_viscosity
  public :: destroy_turbulent_viscosity
  public :: calculate_turbulent_viscosity

  contains

    subroutine setup_turbulent_viscosity()
      implicit none

      call alloc(mu_t, -2,imx+2, -2,jmx+2, -2,kmx+2)

      select case (trim(turbulence))

        case ('none')
          !do nothing
          continue

        case ('sst')
          sst_mu(-2:imx+2,-2:jmx+2,-2:kmx+2) => mu_t(:,:,:)

        case ('kkl')
          kkl_mu(-2:imx+2,-2:jmx+2,-2:kmx+2) => mu_t(:,:,:)

        case DEFAULT 
          !call turbulence_read_error()
          Error

      end select

    end subroutine setup_turbulent_viscosity

  
    subroutine destroy_turbulent_viscosity()
      implicit none

      select case (trim(turbulence))

        case ('none')
          !do nothing
          continue

        case ('sst')
          nullify(sst_mu)

        case ('kkl')
          nullify(kkl_mu)

        case DEFAULT 
          !call turbulence_read_error()
          Error

      end select
      call dealloc(mu_t)

    end subroutine destroy_turbulent_viscosity


    subroutine calculate_turbulent_viscosity()
      implicit none

      select case (trim(turbulence))

        case ('none')
          !do nothing
          continue

        case ('sst')
          call calculate_sst_mu()

        case ('kkl')
          call calculate_kkl_mu()

        case DEFAULT 
          !call turbulence_read_error()
          Error

      end select

    end subroutine calculate_turbulent_viscosity


end module turbulent_viscosity

