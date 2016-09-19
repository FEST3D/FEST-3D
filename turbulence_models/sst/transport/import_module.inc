!should be included in ausm before variable deceleration

    use state, only:  tk, tw
    use face_interpolant, only: x_tk_left, x_tw_left, x_tk_right, x_tw_right, &
                                y_tk_left, y_tw_left, y_tk_right, y_tw_right, &
                                z_tk_left, z_tw_left, z_tk_right, z_tw_right
    use source, only: gradtk_x, gradtk_y, gradtk_z, gradtw_x, gradtw_y, gradtw_z,&
                      F1, sst_mu
