! include in paralle module, send_recv_j subroutine

      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/parallel/recv_jmax.inc"

          case ("kkl")
              include "turbulence_models/kkl/parallel/recv_jmax.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/parallel", "recv_jmax", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
