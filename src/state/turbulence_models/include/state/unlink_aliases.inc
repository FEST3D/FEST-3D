!should be included in state state, in subroutine with same name.

      select case (trim(turbulence))

          case ("none")
              !include nothing
              continue
          
          case ("sst")
              include "turbulence_models/sst/state/unlink_aliases.inc"

          case ("kkl")
              include "turbulence_models/kkl/state/unlink_aliases.inc"

          case ("sa")
              include "turbulence_models/sa/state/unlink_aliases.inc"

          case ("kw")
              include "turbulence_models/kw/state/unlink_aliases.inc"

          case ("ke")
              include "turbulence_models/ke/state/unlink_aliases.inc"

          case ("des-sst")
              include "turbulence_models/des-sst/state/unlink_aliases.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/state", "unlink_aliases", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
