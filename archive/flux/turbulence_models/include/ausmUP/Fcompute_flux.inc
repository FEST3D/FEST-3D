!should be include in ausmUP.f90 in compute_flux after F_minus ahs been computed


      select case (turbulence)
          
          case ("none")
              !include nothing
              continue

          case ("sst", "sst2003")
              include "turbulence_models/sst/ausmUP/Fcompute_flux.inc"

          case ("kkl")
              include "turbulence_models/kkl/ausmUP/Fcompute_flux.inc"

          case ("sa", "saBC")
              include "turbulence_models/sa/ausmUP/Fcompute_flux.inc"

          case ("kw")
              include "turbulence_models/kw/ausmUP/Fcompute_flux.inc"

          case ("ke")
              include "turbulence_models/ke/ausmUP/Fcompute_flux.inc"

          case ("des-sst")
              include "turbulence_models/des-sst/ausmUP/Fcompute_flux.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/ausmUP", "Fcompute_flux", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
