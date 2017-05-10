!module boundary_conditions
!TODO Index 6 and 7 are wrong , check logrith definintion in main BC file

      if(turbulence == "none")then
        !do nothing
        continue
      else if (turbulence == "sst")then
        include "turbulence_models/sst/bc/apply_boundary_condition_jmax.inc"
      else
          print*, "***** Turbulence models not recognised ***"
      end if
