!include in viscous (transport) module

        real :: tk_face, tw_face

        call dmsg(1, 'turbuence_models/sst/source', 'compute_cell_center_grad')
        if (turbulence == "sst")then
          F1 = 0.0
          sst_mu =0.0
          gradtk_x = 0.0
          gradtk_y = 0.0
          gradtk_z = 0.0
          gradtw_x = 0.0
          gradtw_y = 0.0
          gradtw_z = 0.0
        end if

