!should be include in ausm.f90 in compute_flux after F_minus ahs been computed


      select case (turbulence)
          
          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/ausm/Fcompute_flux.inc"

          case ("kkl")
              include "turbulence_models/kkl/ausm/Fcompute_flux.inc"

          case ("sa", "saBC")
              include "turbulence_models/sa/ausm/Fcompute_flux.inc"

          case ("kw")
              include "turbulence_models/kw/ausm/Fcompute_flux.inc"

          case ("ke")
              include "turbulence_models/ke/ausm/Fcompute_flux.inc"

          case ("des-sst")
              include "turbulence_models/des-sst/ausm/Fcompute_flux.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/ausm", "Fcompute_flux", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
