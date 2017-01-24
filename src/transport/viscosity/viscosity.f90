module viscosity
  !----------------------------------------------------
  ! author   - Jatinder Pal Singh Sandhu
  ! obective - setup     (call only),
  !            destroy   (call only), and 
  !            calculate (call only)
  ! molecular and turbulence viscosity based on switch
  !-----------------------------------------------------

  use global_vars          , only : mu_ref
  use global_vars          , only : turbulence
  use molecular_viscosity  , only : calculate_molecular_viscosity
  use turbulent_viscosity  , only : calculate_turbulent_viscosity
  use molecular_viscosity  , only : setup_molecular_viscosity
  use turbulent_viscosity  , only : setup_turbulent_viscosity
  use molecular_viscosity  , only : destroy_molecular_viscosity
  use turbulent_viscosity  , only : destroy_turbulent_viscosity

  implicit none
  private

  public :: setup_viscosity
  public :: destroy_viscosity
  public :: calculate_viscosity

  contains

    subroutine calculate_viscosity()

      if (mu_ref/=0.) then
        call calculate_molecular_viscosity()
      end if

      if (turbulence/='none') then
        call calculate_turbulent_viscosity()
      end if

    end subroutine calculate_viscosity

    subroutine setup_viscosity()

      if (mu_ref/=0.) then
        call setup_molecular_viscosity()
      end if

      if (turbulence/='none') then
        call setup_turbulent_viscosity()
      end if

    end subroutine setup_viscosity

    subroutine destroy_viscosity()

      if (mu_ref/=0.) then
        call destroy_molecular_viscosity()
      end if

      if (turbulence/='none') then
        call destroy_turbulent_viscosity()
      end if

    end subroutine destroy_viscosity

end module viscosity
