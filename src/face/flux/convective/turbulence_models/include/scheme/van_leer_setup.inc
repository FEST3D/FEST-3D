!to  include in schemee.f90/setup_scheme/vanleer


      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/scheme/van_leer_setup.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/scheme", "van_leer_setup", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
