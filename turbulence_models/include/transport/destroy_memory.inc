! include in viscous module (transport) in setup_viscous

      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/transport/destroy_memory.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/transport", "destroy_memory", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
