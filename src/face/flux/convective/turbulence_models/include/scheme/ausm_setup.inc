!to  include in schemee.f90/setup_scheme/ausm


      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/scheme/ausm_setup.inc"

          case ("kkl")
              include "turbulence_models/kkl/scheme/ausm_setup.inc"

          case ("sa", "saBC")
              include "turbulence_models/sa/scheme/ausm_setup.inc"

          case ("ke")
              include "turbulence_models/ke/scheme/ausm_setup.inc"

          case ("kw")
              include "turbulence_models/kw/scheme/ausm_setup.inc"

          case ("des-sst")
              include "turbulence_models/des-sst/scheme/ausm_setup.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/scheme", "ausm_setup", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
