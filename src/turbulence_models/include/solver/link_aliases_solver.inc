!include in solver.f90/link_aliases_solver

      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/solver/link_aliases_solver.inc"

          case ("kkl")
              include "turbulence_models/kkl/solver/link_aliases_solver.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/state", "link_aliases_solver", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
