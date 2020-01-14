  !< Check for solution's convergence
module convergence
  !< Check for solution's convergence
  use vartypes
  use resnorm    , only: Res_abs
  use resnorm    , only: Res_rel
#include "error.h"
  implicit none
  private
  public :: converged

  contains
    function converged(control) result(c)
        !< Check if the solution seems to have converged
        !< The solution is said to have converged if the change in 
        !< the residue norm is "negligible".
        !-----------------------------------------------------------

        implicit none
        type(controltype), intent(in) :: control
        !< control parameters
        logical :: c
        !< convergence result:True or false
        real(wp)    :: check=10.

        select case(trim(control%tolerance_type))
          case('Mass_abs')
            check =  Res_abs(0)

          case('Resnorm_abs')
            check =  sqrt(sum(Res_abs(1:)**2))

          case('Viscous_abs')
            check =  sqrt(sum(Res_abs(1:5)**2))

          case('Turbulent_abs')
            check =  sqrt(sum(Res_abs(6:)**2))

          case('Continuity_abs')
            check =  Res_abs(1)

          case('X-mom_abs')
            check =  Res_abs(2)

          case('Z-mom_abs')
            check =  Res_abs(3)

          case('Y-mom_abs')
            check =  Res_abs(4)

          case('Energy_abs')
            check =  Res_abs(5)

          case('Mass_rel')
            check =  Res_rel(0)

          case('Resnorm_rel')
            check =  sqrt(sum(Res_rel(1:)**2))

          case('Viscous_rel')
            check =  sqrt(sum(Res_rel(1:5)**2))

          case('Turbulent_rel')
            check =  sqrt(sum(Res_rel(6:)**2))

          case('Continuity_rel')
            check =  Res_rel(1)

          case('X-mom_rel')
            check =  Res_rel(2)

          case('Z-mom_rel')
            check =  Res_rel(3)

          case('Y-mom_rel')
            check =  Res_rel(4)

          case('Energy_rel')
            check =  Res_rel(5)

          case('TKE_abs')
            check =  Res_abs(6)

          case('tv_abs')
            check =  Res_abs(6)

          case('Dissipation_abs')
            check =  Res_abs(7)

          case('Omega_abs')
            check =  Res_abs(7)

          case('Kl_abs')
            check =  Res_abs(7)

          case('TKE_rel')
            check =  Res_rel(6)

          case('tv_rel')
            check =  Res_rel(6)

          case('Dissipation_rel')
            check =  Res_rel(7)

          case('Omega_rel')
            check =  Res_rel(7)

          case('Kl_rel')
            check =  Res_rel(7)

          case DEFAULT
            ! making absolute resnorm default
            check =  sqrt(sum(Res_abs(1:)**2))
            Issue_warning
        end select

        if (check < control%tolerance .and. control%current_iter>10) then
          c = .TRUE.
        else
          c = .FALSE.
        end if

    end function converged

end module convergence

