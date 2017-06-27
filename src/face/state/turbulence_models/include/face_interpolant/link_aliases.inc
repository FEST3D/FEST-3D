!should be included in face_interpolant.f90 file
        

!should be included in state state, in subroutine with same name.

      select case (turbulence)
 
          case ("none")
              !include nothing
              continue
         
          case ("sst")
              include "turbulence_models/sst/face_interpolant/link_aliases.inc"

          case ("kkl")
              include "turbulence_models/kkl/face_interpolant/link_aliases.inc"

          case ("sa")
              include "turbulence_models/sa/face_interpolant/link_aliases.inc"

          case ("ke")
              include "turbulence_models/ke/face_interpolant/link_aliases.inc"

          case ("kw")
              include "turbulence_models/kw/face_interpolant/link_aliases.inc"

          case ("des-sst")
              include "turbulence_models/des-sst/face_interpolant/link_aliases.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/face_interpolant", "link_aliases", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
