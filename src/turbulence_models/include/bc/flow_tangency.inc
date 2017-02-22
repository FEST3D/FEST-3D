!should be included in subroutine - flow_tangency 
!    in module boundary conditons

  if(turbulence == "none")then 
    !do nothing
    continue
  else if(turbulence == "sst")then
    include "turbulence_models/sst/bc/flow_tangency.inc"
  else 
    print*, "Turbulence nodel not recongnised, tbc, flow_tangency, sst"
  end if
