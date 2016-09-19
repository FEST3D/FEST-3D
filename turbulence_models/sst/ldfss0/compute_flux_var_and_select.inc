!should be included in ausm file in compte flux subroutine

            real, dimension(:, :, :), pointer :: f_tk_left, f_tk_right, &
                                                 f_tw_left, f_tw_right

            call dmsg(1, 'turbulence_models/sst/ldfss0', 'compute_flux.inc')
            
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
                    call dmsg(5, 'turbulence_models/sst/ldfss0', 'compute_flux.inc', &
                            'Direction not recognised')
                    stop
            end select
