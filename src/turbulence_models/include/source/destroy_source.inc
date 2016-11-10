! to be include in source/destroy_source/


        call dmsg(1, 'turbulence_models/sst/source', 'deallocate_memory.inc')

      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/source/destroy_source.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/source", "destroy_source", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
