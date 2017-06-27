!should be included in ausm file in compte flux subroutine
            
            select case (f_dir)
                case ('x')
                    f_tv_left  => x_tv_left
                    f_tv_right => x_tv_right
                case ('y')
                    f_tv_left  => y_tv_left
                    f_tv_right => y_tv_right
                case ('z')
                    f_tv_left  => z_tv_left
                    f_tv_right => z_tv_right
                case default
                    call dmsg(5, 'turbulence_models/sa/ldfss0', 'compute_flux.inc', &
                            'Direction not recognised')
                    stop
            end select

