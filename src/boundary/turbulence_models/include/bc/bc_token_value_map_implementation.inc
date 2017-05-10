!should be included in inc file : "bc_token_value_map_implementation.inc"

      print*, "turbulence_token"
      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              print*, "inside case"
              include "turbulence_models/sst/bc/bc_token_value_map_implementation.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/bc", "bc_token_value_map_implementation.inc", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
