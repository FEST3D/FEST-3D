!include in source module

  use state, only: tk, tw
 ! use viscous, only:gradu_x, gradu_y, gradu_z, gradv_x, gradv_y, gradv_z,&
 !                   gradw_x, gradw_y, gradw_z, gradT_x, gradT_y, gradT_z,&
 !                   gradtk_x, gradtk_y, gradtk_z, gradtw_x, gradtw_y, gradtw_z
  use scheme, only: TKE_residue, omega_residue

  use face_interpolant, only: x_tk_left, x_tw_left, x_tk_right, x_tw_right, &
                                y_tk_left, y_tw_left, y_tk_right, y_tw_right, &
                                z_tk_left, z_tw_left, z_tk_right, z_tw_right

  use wall_dist, only: dist
