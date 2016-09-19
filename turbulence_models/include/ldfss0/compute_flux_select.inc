!should be included in ldfss0 file in compte flux subroutine
            

      select case (turbulence)
          
          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/ldfss0/compute_flux_select.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/ldfss0", "compute_flux_select", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
