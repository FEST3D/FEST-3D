!should be included in ldfss0 file in compte flux subroutine
            

      select case (turbulence)
          
          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/ldfss0/compute_flux_select.inc"

          case ("kkl")
              include "turbulence_models/kkl/ldfss0/compute_flux_select.inc"

          case ("sa", "saBC")
              include "turbulence_models/sa/ldfss0/compute_flux_select.inc"

          case ("kw")
              include "turbulence_models/kw/ldfss0/compute_flux_select.inc"

          case ("ke")
              include "turbulence_models/ke/ldfss0/compute_flux_select.inc"

          case ("des-sst")
              include "turbulence_models/des-sst/ldfss0/compute_flux_select.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/ldfss0", "compute_flux_select", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
