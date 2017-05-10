!should be included in subroutine - no_slip 
!    in module boundary conditons

  if(turbulence == "none")then 
    !do nothing
    continue
  else if(turbulence == "sst")then
    include "turbulence_models/sst/bc/no_slip.inc"
  else 
    print*, "Turbulence nodel not recongnised, tbc, no_slip, sst"
  end if
