!should be included in ausm file in compte flux subroutine
            select case (f_dir)
                case ('x')
                    f_tk_left  => x_tk_left
                    f_tk_right => x_tk_right
                    f_tkl_left  => x_tkl_left
                    f_tkl_right => x_tkl_right
                case ('y')
                    f_tk_left  => y_tk_left
                    f_tk_right => y_tk_right
                    f_tkl_left  => y_tw_left
                    f_tkl_right => y_tw_right
                case ('z')
                    f_tk_left  => z_tk_left
                    f_tk_right => z_tk_right
                    f_tkl_left  => z_tkl_left
                    f_tkl_right => z_tkl_right
                case default
                    call dmsg(5, 'turbulence_models/kkl/ldfss0', 'compute_flux.inc', &
                            'Direction not recognised')
                    stop
            end select

