!should be included in ausm file in compte flux subroutine
            

      select case (turbulence)
          
          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/ausm/compute_flux_select.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/ausm", "compute_flux_select", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
