! include in viscous (transport module) 


        call dmsg(1, 'turbulence_models/sst/transport', 'deallocate_memory.inc')

        call dealloc(gradtk_x)
        call dealloc(gradtk_y)
        call dealloc(gradtk_z)

        call dealloc(gradtw_x)
        call dealloc(gradtw_y)
        call dealloc(gradtw_z)
