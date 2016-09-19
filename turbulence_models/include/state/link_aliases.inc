!should be included in state state, in subroutine with same name.

      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/state/link_aliases.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/state", "link_aliases", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
