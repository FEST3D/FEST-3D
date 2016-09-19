!to  include in schemee.f90/setup_scheme/ldfss0   


      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/scheme/ldfss0_setup.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/scheme", "ldfss0_setup", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
