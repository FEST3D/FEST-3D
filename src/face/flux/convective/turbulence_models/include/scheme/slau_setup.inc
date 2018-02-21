!to  include in schemee.f90/setup_scheme/slau


      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/scheme/slau_setup.inc"

          case ("kkl")
              include "turbulence_models/kkl/scheme/slau_setup.inc"

          case ("sa")
              include "turbulence_models/sa/scheme/slau_setup.inc"

          case ("ke")
              include "turbulence_models/ke/scheme/slau_setup.inc"

          case ("kw")
              include "turbulence_models/kw/scheme/slau_setup.inc"

          case ("des-sst")
              include "turbulence_models/des-sst/scheme/slau_setup.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/scheme", "slau_setup", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
