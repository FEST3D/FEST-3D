!should be include in ausm.f90 in compute_flux after F_minus ahs been computed

            
!            call dmsg(1, 'turbulence_models/sst/ausm', 'compute_flux.inc')
            
            select case (f_dir)
                case ('x')
                    f_tk_left => x_tk_left
                    f_tk_right => x_tk_right
                    f_tw_left => x_tw_left
                    f_tw_right => x_tw_right
                case ('y')
                    f_tk_left => y_tk_left
                    f_tk_right => y_tk_right
                    f_tw_left => y_tw_left
                    f_tw_right => y_tw_right
                case ('z')
                    f_tk_left => z_tk_left
                    f_tk_right => z_tk_right
                    f_tw_left => z_tw_left
                    f_tw_right => z_tw_right
                case default
                    call dmsg(5, 'turbulence_models/sst/ausm', 'compute_flux.inc', &
                            'Direction not recognised')
                    stop
            end select
                ! Construct other fluxes in terms of the F mass flux
                F_plus(5) =  F_plus(5) + (F_plus(1) * f_tk_left(i, j, k))
                F_plus(6) = (F_plus(1) * f_tk_left(i, j, k)) 
                F_plus(7) = (F_plus(1) * f_tw_left(i, j, k))

                
                ! Construct other fluxes in terms of the F mass flux
                F_minus(5) =  F_minus(5) + (F_minus(1) * f_tk_right(i, j, k))
                F_minus(6) = (F_minus(1) * f_tk_right(i, j, k))
                F_minus(7) = (F_minus(1) * f_tw_right(i, j, k)) 

