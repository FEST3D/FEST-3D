!should be included in state state, in subroutine with same name.

      select case (trim(turbulence))
          
          case ("none")
              !include nothing
              continue

          case ("sst", "bsl")
              include "turbulence_models/sst/state/init_infinity_values.inc"

          case ("kkl")
              include "turbulence_models/kkl/state/init_infinity_values.inc"

          case ("sa")
              include "turbulence_models/sa/state/init_infinity_values.inc"

          case ("saBC")
              include "turbulence_models/saBC/state/init_infinity_values.inc"

          case ("kw")
              include "turbulence_models/kw/state/init_infinity_values.inc"

          case ("ke")
              include "turbulence_models/ke/state/init_infinity_values.inc"

          case ("des-sst")
              include "turbulence_models/des-sst/state/init_infinity_values.inc"

          case ("les")
              include "turbulence_models/les/state/init_infinity_values.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/state", "init_infinity_values", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
