module transport

  !------------------------------------------------------------------
  ! author    - Jatinder Pal Singh Sandhu
  ! objective - setup     (call only)
  !             destroy   (call only)
  !             calculate (call only)
  ! all transport properties that are required by solver 
  !------------------------------------------------------------------

  use viscosity, only : setup_viscosity
  use viscosity, only : destroy_viscosity
  use viscosity, only : calculate_viscosity

  implicit none
  private

  public :: setup_transport
  public :: destroy_transport
  public :: calculate_transport

  contains

    subroutine setup_transport()

      call setup_viscosity

    end subroutine setup_transport

    subroutine destroy_transport()

      call destroy_viscosity()

    end subroutine destroy_transport

    subroutine calculate_transport()

      call calculate_viscosity()

    end subroutine calculate_transport

end module transport
