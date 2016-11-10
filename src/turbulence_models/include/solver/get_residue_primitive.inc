!include in solver.f90
            
      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/solver/get_residue_primitive.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/solver", "get_residue_primitive", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
