
!module boundary_conditions

      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/bc/set_wall_bc_at_faces.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/bc", "set_wall_bc_at_faces", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
