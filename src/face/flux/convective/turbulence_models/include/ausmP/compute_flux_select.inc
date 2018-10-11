!should be include in ausmP.f90 in compute_flux after F_minus ahs been computed


      select case (trim(turbulence))
          
          case ("none")
              !include nothing
              continue

          case ("sst", "sst2003")
              include "turbulence_models/sst/ausmP/compute_flux_select.inc"

          case ("kkl")

              include "turbulence_models/kkl/ausmP/compute_flux_select.inc"

          case ("sa", "saBC")
              include "turbulence_models/sa/ausmP/compute_flux_select.inc"

          case ("kw")
              include "turbulence_models/kw/ausmP/compute_flux_select.inc"

          case ("ke")
              include "turbulence_models/ke/ausmP/compute_flux_select.inc"

          case ("des-sst")
              include "turbulence_models/des-sst/ausmP/compute_flux_select.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/ausmP", "compute_flux_select", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
