!should be inlcuded in source/setup_source/ also initialize here


      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/source/setup_source.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/source", "setup_source", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
