! include in viscous module (transport) in setup_viscous


        call dmsg(1, 'turbulence_models/sst/transport', 'setup_viscous')

        call alloc(gradtk_x, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradtk_x - sst/transport')
        call alloc(gradtk_y, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradtk_y - sst/transport')
        call alloc(gradtk_z, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradtk_z - sst/transport')

        call alloc(gradtw_x, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradtw_x - sst/transport')
        call alloc(gradtw_y, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradtw_y - sst/transport')
        call alloc(gradtw_z, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradtw_z - sst/transport')
