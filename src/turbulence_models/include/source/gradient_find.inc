!should be inlcuded in source/gradient_find/ 


      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/source/gradient_find.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/source", "gradient_find", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
