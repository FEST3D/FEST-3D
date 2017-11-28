!should be included in state state, in subroutine with same name.

      select case (trim(turbulence))

          case ("none")
              !include nothing
              continue

          case ("sst", "bsl")
              include "turbulence_models/sst/state/link_aliases.inc"

          case ("kkl")
              include "turbulence_models/kkl/state/link_aliases.inc"

          case ("sa")
              include "turbulence_models/sa/state/link_aliases.inc"

          case ("kw")
              include "turbulence_models/kw/state/link_aliases.inc"

          case ("ke")
              include "turbulence_models/ke/state/link_aliases.inc"

          case ("des-sst")
              include "turbulence_models/des-sst/state/link_aliases.inc"

          case ("les")
              include "turbulence_models/les/state/link_aliases.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/state", "link_aliases", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
