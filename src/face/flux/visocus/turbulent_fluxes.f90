module turbulent_fluxes
  !-----------------------------------------------------------------
  !170609  - jatinder Pal Singh Sandhu
  ! AIM: to add source term residue to already calculate residuals
  !      if requires
  !-----------------------------------------------------------------
#include "../../../error.inc"
  use global_vars, only : turbulence
  use utils      , only : dmsg
  use utils      , only : turbulence_read_error
  use sst_turbulent_flux , only : compute_sst_fluxes
  use kkl_turbulent_flux , only : compute_kkl_fluxes

  implicit none
  private
  public :: compute_turbulent_fluxes

  contains

    
    subroutine compute_turbulent_fluxes(F,G,H)

      implicit none
      real, dimension(:, :, :, :), pointer :: F, G, H

      call dmsg(1, 'turbulent_fluexes', 'compute_turbulent_fluxes')

      select case (trim(turbulence))

        case ('none')
          !do nothing
          continue

        case ('sst')
          call compute_sst_fluxes(F,G,H)

        case ('kkl')
          call compute_kkl_fluxes(F,G,H)

        case DEFAULT
          !call turbulence_read_error()
          Error

      end select

    end subroutine compute_turbulent_fluxes


end module turbulent_fluxes


