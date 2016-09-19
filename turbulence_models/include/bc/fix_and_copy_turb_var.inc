!module boundary_conditions

!!! SELECT COMMENTED OUT DUE TO ERROR ENCOUNTER

!      select case (turbulence)
!
!          case ("none")
!              !include nothing
!              continue
!
!          case ("sst")
!              include "turbulence_models/sst/bc/fix_and_copy_turb_var.inc"
!
!          case DEFAULT
!              call dmsg(5, "turbulence_models/include/bc", "fix_and_copy_turb_var", &
!                        "ERROR: Turbulence model not recognised")
!              STOP
!
!      end select
              include "turbulence_models/sst/bc/fix_and_copy_turb_var.inc"
