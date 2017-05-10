! should be include in 'include' file : bclist_definition in main file
!  use slect
    
      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/bc/bclist_definition.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/bc", "bclist_definition", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
