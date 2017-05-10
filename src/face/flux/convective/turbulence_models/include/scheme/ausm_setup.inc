!to  include in schemee.f90/setup_scheme/ausm


      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/scheme/ausm_setup.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/scheme", "ausm_setup", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
