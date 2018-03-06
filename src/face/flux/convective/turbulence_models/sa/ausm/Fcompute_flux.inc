!should be include in ausm.f90 in compute_flux after F_minus ahs been computed

            
  F_plus(6) = (F_plus(1) * f_tv_left(i, j, k))
  
  ! Construct other fluxes in terms of the F mass flux
  F_minus(6) = (F_minus(1) * f_tv_right(i, j, k))

