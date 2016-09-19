!should be included in state file, in subroutine name same as this file name 

            !-----------------------------------------------------------
            ! Initialize the state based on the infinity values.
            !-----------------------------------------------------------
            
            
            call dmsg(1, 'turbulence/sst/t_state', 'init_state_with_infinity_values')
            
            do i = 1,sst_n_var
                tqp(:, :, :, i) = tqp_inf(i)
            end do  
