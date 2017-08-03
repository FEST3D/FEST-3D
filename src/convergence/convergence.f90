module convergence
  !-------------------------------------
  ! 170803  -Jatinder Pal Singh Sandhu
  !           find if solution converged
  !-------------------------------------
  use global_vars, only: Res_abs
  use global_vars, only: Res_rel
  use global_vars, only: tolerance
  use global_vars, only: tolerance_type
#include "../error.inc"
  implicit none
  private
  public :: converged

  contains
    function converged() result(c)
        !-----------------------------------------------------------
        ! Check if the solution seems to have converged
        !
        ! The solution is said to have converged if the change in 
        ! the residue norm is "negligible".
        !-----------------------------------------------------------

        implicit none
        logical :: c
        real    :: check=10.

        select case(trim(tolerance_type))
        include "convergence_select.inc"
        end select

        if (check < tolerance) then
          c = .TRUE.
        else
          c = .FALSE.
        end if

    end function converged

end module convergence

