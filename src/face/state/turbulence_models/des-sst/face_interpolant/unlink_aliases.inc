!should be included in face_nterpolant.f90 file
            
            ! Unlink xi faces left pointers
            nullify(x_tk_left)
            nullify(x_tw_left)

            ! Unlink xi faces right pointers
            nullify(x_tk_right)
            nullify(x_tw_right)

            ! Unlink eta faces left pointers
            nullify(y_tk_left)
            nullify(y_tw_left)

            ! Unlink eta faces right pointers
            nullify(y_tk_right)
            nullify(y_tw_right)
            
            ! Unlink tau faces left pointers
            nullify(z_tk_left)
            nullify(z_tw_left)

            ! Unlink tau faces right pointers
            nullify(z_tk_right)
            nullify(z_tw_right)


