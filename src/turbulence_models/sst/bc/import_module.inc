! include in boundary condition.f90

!    use state, only:  tk, tw, tk_inf, tw_inf, Sutherland_temp, T_ref, mu_ref,&
!                      R_gas

    use face_interpolant, only: x_tk_left, x_tk_right,&
                y_tk_left, y_tk_right, z_tk_left, z_tk_right, &
                x_tw_left, x_tw_right, y_tw_left, y_tw_right, &
                z_tw_left, z_tw_right, x_density_right, x_density_left,&
                y_density_left, y_density_right, z_density_left, z_density_right, &
                 x_pressure_right, x_pressure_left,&
                y_pressure_left, y_pressure_right, z_pressure_left, z_pressure_right
!    use wall_dist, only: dist
    use source, only: beta1
