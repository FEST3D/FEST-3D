!to  include in schemee.f90/setup_scheme/vanleer


      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/scheme/van_leer_setup.inc"

          case ("kkl")
              include "turbulence_models/kkl/scheme/van_leer_setup.inc"

          case ("sa")
              include "turbulence_models/sa/scheme/van_leer_setup.inc"

          case ("ke")
              include "turbulence_models/ke/scheme/van_leer_setup.inc"

          case ("kw")
              include "turbulence_models/kw/scheme/van_leer_setup.inc"

          case ("des-sst")
              include "turbulence_models/des-sst/scheme/van_leer_setup.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/scheme", "van_leer_setup", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
