!include in viscous (transport) module 

! turbulent kinectic energy gradient
            ! x_face(i, j, k)
            tk_face = (x_tk_left(i,j,k)+ x_tk_right(i,j,k))
            gradtk_x(i, j, k) = gradtk_x(i, j, k) - &
                               (tk_face * xnx(i, j, k) * xA(i, j, k))
            gradtk_y(i, j, k) = gradtk_y(i, j, k) - &
                               (tk_face * xny(i, j, k) * xA(i, j, k))
            gradtk_z(i, j, k) = gradtk_z(i, j, k) - &
                               (tk_face * xnz(i, j, k) * xA(i, j, k))

            ! x_face(i+1, j, k)
            tk_face = (x_tk_left(i+1,j,k)+ x_tk_right(i+1,j,k))
            gradtk_x(i, j, k) = gradtk_x(i, j, k) + &
                               (tk_face * xnx(i+1, j, k) * xA(i+1, j, k))
            gradtk_y(i, j, k) = gradtk_y(i, j, k) + &
                               (tk_face * xny(i+1, j, k) * xA(i+1, j, k))
            gradtk_z(i, j, k) = gradtk_z(i, j, k) + &
                               (tk_face * xnz(i+1, j, k) * xA(i+1, j, k))

            ! y_face(i, j, k)
            tk_face = (y_tk_left(i,j,k)+ y_tk_right(i,j,k))
            gradtk_x(i, j, k) = gradtk_x(i, j, k) - &
                               (tk_face * ynx(i, j, k) * yA(i, j, k))
            gradtk_y(i, j, k) = gradtk_y(i, j, k) - &
                               (tk_face * yny(i, j, k) * yA(i, j, k))
            gradtk_z(i, j, k) = gradtk_z(i, j, k) - &
                               (tk_face * ynz(i, j, k) * yA(i, j, k))

            ! y_face(i, j+1, k)
            tk_face = (y_tk_left(i,j+1,k)+ y_tk_right(i,j+1,k))
            gradtk_x(i, j, k) = gradtk_x(i, j, k) + &
                               (tk_face * ynx(i, j+1, k) * yA(i, j+1, k))
            gradtk_y(i, j, k) = gradtk_y(i, j, k) + &
                               (tk_face * yny(i, j+1, k) * yA(i, j+1, k))
            gradtk_z(i, j, k) = gradtk_z(i, j, k) + &
                               (tk_face * ynz(i, j+1, k) * yA(i, j+1, k))

            ! z_face(i, j, k)
            tk_face = (z_tk_left(i,j,k)+ z_tk_right(i,j,k))
            gradtk_x(i, j, k) = gradtk_x(i, j, k) - &
                               (tk_face * znx(i, j, k) * zA(i, j, k))
            gradtk_y(i, j, k) = gradtk_y(i, j, k) - &
                               (tk_face * zny(i, j, k) * zA(i, j, k))
            gradtk_z(i, j, k) = gradtk_z(i, j, k) - &
                               (tk_face * znz(i, j, k) * zA(i, j, k))

            ! z_face(i, j, k+1)
            tk_face = (z_tk_left(i,j,k+1)+ z_tk_right(i,j,k+1))
            gradtk_x(i, j, k) = gradtk_x(i, j, k) + &
                               (tk_face * znx(i, j, k+1) * zA(i, j, k+1))
            gradtk_y(i, j, k) = gradtk_y(i, j, k) + &
                               (tk_face * zny(i, j, k+1) * zA(i, j, k+1))
            gradtk_z(i, j, k) = gradtk_z(i, j, k) + &
                               (tk_face * znz(i, j, k+1) * zA(i, j, k+1))

            ! Factor of volume
            gradtk_x(i, j, k) = gradtk_x(i, j, k) * 0.5 / vol(i, j, k)
            gradtk_y(i, j, k) = gradtk_y(i, j, k) * 0.5 / vol(i, j, k)
            gradtk_z(i, j, k) = gradtk_z(i, j, k) * 0.5 / vol(i, j, k)

!omega gradient

   ! x_face(i, j, k)
            tw_face = (x_tw_left(i,j,k)+ x_tw_right(i,j,k))
            gradtw_x(i, j, k) = gradtw_x(i, j, k) - &
                               (tw_face * xnx(i, j, k) * xA(i, j, k))
            gradtw_y(i, j, k) = gradtw_y(i, j, k) - &
                               (tw_face * xny(i, j, k) * xA(i, j, k))
            gradtw_z(i, j, k) = gradtw_z(i, j, k) - &
                               (tw_face * xnz(i, j, k) * xA(i, j, k))

            ! x_face(i+1, j, k)
            tw_face = (x_tw_left(i+1,j,k)+ x_tw_right(i+1,j,k))
            gradtw_x(i, j, k) = gradtw_x(i, j, k) + &
                               (tw_face * xnx(i+1, j, k) * xA(i+1, j, k))
            gradtw_y(i, j, k) = gradtw_y(i, j, k) + &
                               (tw_face * xny(i+1, j, k) * xA(i+1, j, k))
            gradtw_z(i, j, k) = gradtw_z(i, j, k) + &
                               (tw_face * xnz(i+1, j, k) * xA(i+1, j, k))

            ! y_face(i, j, k)
            tw_face = (y_tw_left(i,j,k)+ y_tw_right(i,j,k))
            gradtw_x(i, j, k) = gradtw_x(i, j, k) - &
                               (tw_face * ynx(i, j, k) * yA(i, j, k))
            gradtw_y(i, j, k) = gradtw_y(i, j, k) - &
                               (tw_face * yny(i, j, k) * yA(i, j, k))
            gradtw_z(i, j, k) = gradtw_z(i, j, k) - &
                               (tw_face * ynz(i, j, k) * yA(i, j, k))

            ! y_face(i, j+1, k)
            tw_face = (y_tw_left(i,j+1,k)+ y_tw_right(i,j+1,k))
            gradtw_x(i, j, k) = gradtw_x(i, j, k) + &
                               (tw_face * ynx(i, j+1, k) * yA(i, j+1, k))
            gradtw_y(i, j, k) = gradtw_y(i, j, k) + &
                               (tw_face * yny(i, j+1, k) * yA(i, j+1, k))
            gradtw_z(i, j, k) = gradtw_z(i, j, k) + &
                               (tw_face * ynz(i, j+1, k) * yA(i, j+1, k))

            ! z_face(i, j, k)
            tw_face = (z_tw_left(i,j,k)+ z_tw_right(i,j,k))
            gradtw_x(i, j, k) = gradtw_x(i, j, k) - &
                               (tw_face * znx(i, j, k) * zA(i, j, k))
            gradtw_y(i, j, k) = gradtw_y(i, j, k) - &
                               (tw_face * zny(i, j, k) * zA(i, j, k))
            gradtw_z(i, j, k) = gradtw_z(i, j, k) - &
                               (tw_face * znz(i, j, k) * zA(i, j, k))

            ! z_face(i, j, k+1)
            tw_face = (z_tw_left(i,j,k+1)+ z_tw_right(i,j,k+1))
            gradtw_x(i, j, k) = gradtw_x(i, j, k) + &
                               (tw_face * znx(i, j, k+1) * zA(i, j, k+1))
            gradtw_y(i, j, k) = gradtw_y(i, j, k) + &
                               (tw_face * zny(i, j, k+1) * zA(i, j, k+1))
            gradtw_z(i, j, k) = gradtw_z(i, j, k) + &
                               (tw_face * znz(i, j, k+1) * zA(i, j, k+1))

            ! Factor of volume
            gradtw_x(i, j, k) = gradtw_x(i, j, k) * 0.5 / vol(i, j, k)
            gradtw_y(i, j, k) = gradtw_y(i, j, k) * 0.5 / vol(i, j, k)
            gradtw_z(i, j, k) = gradtw_z(i, j, k) * 0.5 / vol(i, j, k)
