!should be included in state file, in subroutine name same as this file name

            !-----------------------------------------------------------
            ! Allocate memory for the tstate variables
            !
            ! This assumes that imx and jmx and Kmx (the grid size) has been set
            ! within the state module.
            !-----------------------------------------------------------


            call dmsg(1, 'turbulence_models/sst/state', 'allocate_memory.inc')


            call alloc(tqp, -2, imx+2, -2, jmx+2, -2, kmx+2, 1, sst_n_var, &
                    errmsg='Error: Unable to allocate memory for tstate ' // &
                        'variable tqp.')
            call alloc(tqp_inf, 1, sst_n_var, &
                    errmsg='Error: Unable to allocate memory for tstate ' // &
                        'variable tqp_inf.')
