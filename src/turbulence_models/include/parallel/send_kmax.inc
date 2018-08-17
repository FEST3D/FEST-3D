! include in paralle module, send_recv_i subroutine

      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst", "sst2003")
              include "turbulence_models/sst/parallel/send_kmax.inc"

          case ("kkl")
              include "turbulence_models/kkl/parallel/send_kmax.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/parallel", "send_kmax", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select
