module source
  !-----------------------------------------------------------------
  !170609  - jatinder Pal Singh Sandhu
  ! AIM: to add source term residue to already calculate residuals
  !      if requires
  !-----------------------------------------------------------------
#include "../error.inc"
  use global_vars, only : turbulence
  use global_vars, only : process_id
  use utils      , only : dmsg
  use utils      , only : turbulence_read_error
  use sst_source , only : add_sst_source
  use kkl_source , only : add_kkl_source

  implicit none
  private
  public :: add_source_term_residue

  contains

    
    subroutine add_source_term_residue()

      implicit none

      call dmsg(1, 'source', 'add_source_term_residue')

      select case (turbulence)

        case ('none')
          !do nothing
          continue

        case ('sst')
          call add_sst_source()

        case ('kkl')
          call add_kkl_source()

        case DEFAULT
          !call turbulence_read_error()
          Fatal_error

      end select

    end subroutine add_source_term_residue


end module source


