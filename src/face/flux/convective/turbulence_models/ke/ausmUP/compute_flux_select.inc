!should be included in ausm file in compte flux subroutine
            select case (f_dir)
                case ('x')
                    f_tk_left => x_tk_left
                    f_tk_right => x_tk_right
                    f_te_left => x_te_left
                    f_te_right => x_te_right
                case ('y')
                    f_tk_left => y_tk_left
                    f_tk_right => y_tk_right
                    f_te_left => y_te_left
                    f_te_right => y_te_right
                case ('z')
                    f_tk_left => z_tk_left
                    f_tk_right => z_tk_right
                    f_te_left => z_te_left
                    f_te_right => z_te_right
                case default
                    call dmsg(5, 'turbulence_models/ke/ausm', 'compute_flux.inc', &
                            'Direction not recognised')
                    stop
            end select
