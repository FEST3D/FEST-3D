!should be include in ausm.f90 in compute_flux after F_minus ahs been computed


      select case (turbulence)
          
          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/ausm/Fcompute_flux.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/ausm", "Fcompute_flux", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
