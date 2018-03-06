! inlcude in scheme.f90/deallocate-memory
! nullifying the pointers 


      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/scheme/deallocate_memory.inc"

          case ("kkl")
              include "turbulence_models/kkl/scheme/deallocate_memory.inc"

          case ("sa", "saBC")
              include "turbulence_models/sa/scheme/deallocate_memory.inc"

          case ("ke")
              include "turbulence_models/ke/scheme/deallocate_memory.inc"

          case ("kw")
              include "turbulence_models/kw/scheme/deallocate_memory.inc"

          case ("des-sst")
              include "turbulence_models/des-sst/scheme/deallocate_memory.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/scheme", "deallocate_memory_setup", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
            nullify(TKE_residue)
            nullify(omega_residue)
