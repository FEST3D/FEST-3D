!include in boundary condition module.
            

            if ((bc_imn(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
                x_tk_left(1, :, :) = 0.
                x_tk_right(1, :, :) = 0.
                do k = 1, kmx-1
                  do j = 1, jmx-1
                     T_L = x_pressure_left(1, j, k) / (x_density_left(1, j, k) * R_gas)  
                     T_R = x_pressure_right(1, j, k) / (x_density_right(1, j, k) * R_gas)
                     T_face = 0.5 * (T_L + T_R)
                     mu = mu_ref * (T_face / T_ref)**1.5 * ((T_ref + &    
                         Sutherland_temp) / (T_face + Sutherland_temp))
                    rho = 0.5*(x_density_left(1, j, k) + x_density_right(1, j, k))
                    x_tw_left(1, j, k) = 60*mu/(rho*beta1*(dist(1, j, k))**2)
                    x_tw_right(1, j, k) = 60*mu/(rho*beta1*dist(1, j, k)**2)
                  end do
                end do
            else if ((bc_imx(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
                x_tk_left(imx, :, :) = 0.
                x_tk_right(imx, :, :) = 0.
                do k = 1, kmx-1
                  do j = 1, jmx-1
                     T_L = x_pressure_left(imx, j, k) / (x_density_left(imx, j, k) * R_gas)  
                     T_R = x_pressure_right(imx, j, k) / (x_density_right(imx, j, k) * R_gas)
                     T_face = 0.5 * (T_L + T_R)
                     mu = mu_ref * (T_face / T_ref)**1.5 * ((T_ref + &    
                         Sutherland_temp) / (T_face + Sutherland_temp))
                    rho = 0.5*(x_density_left(imx, j, k) + x_density_right(imx, j, k))
                    x_tw_left(imx, j, k) = 60*mu/(rho*beta1*dist(imx, j, k)**2)
                    x_tw_right(imx, j, k) = 60*mu/(rho*beta1*dist(imx, j, k)**2)
                  end do
                end do
            else if ((bc_jmn(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
                y_tk_left(:, 1, :) = 0.
                y_tk_right(:, 1, :) = 0.
                do k = 1, kmx-1
                  do i = 1, imx-1
                     T_L = y_pressure_left(i, 1, k) / (y_density_left(i, 1, k) * R_gas)  
                     T_R = y_pressure_right(i, 1, k) / (y_density_right(i, 1, k) * R_gas)
                     T_face = 0.5 * (T_L + T_R)
                     mu = mu_ref * (T_face / T_ref)**1.5 * ((T_ref + &    
                         Sutherland_temp) / (T_face + Sutherland_temp))
                    rho = 0.5*(y_density_left(i, 1, k) + y_density_right(i, 1, k))
                    y_tw_left(i, 1, k) = 60*mu/(rho*beta1*dist(i,1,k)**2)
                    y_tw_right(i, 1, k) = 60*mu/(rho*beta1*dist(i,1,k)**2)
                  end do
                end do
            else if ((bc_jmx(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
                y_tk_left(:, jmx, :) = 0.
                y_tk_right(:, jmx, :) = 0.
                do k = 1, kmx-1
                  do i = 1, imx-1
                     T_L = y_pressure_left(i, jmx, k) / (y_density_left(i, jmx, k) * R_gas)  
                     T_R = y_pressure_right(i, jmx, k) / (y_density_right(i, jmx, k) * R_gas)
                     T_face = 0.5 * (T_L + T_R)
                     mu = mu_ref * (T_face / T_ref)**1.5 * ((T_ref + &    
                         Sutherland_temp) / (T_face + Sutherland_temp))
                    rho = 0.5*(y_density_left(i,jmx,k) + y_density_right(1,jmx,k))
                    y_tw_left(i, jmx, k) = 60*mu/(rho*beta1*dist(i,jmx,k)**2)
                    y_tw_right(i, jmx, k) = 60*mu/(rho*beta1*dist(i,jmx,k)**2)
                  end do
                end do
            else if ((bc_kmn(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
                z_tk_left(:, :, 1) = 0.
                z_tk_right(:, :, 1) = 0.
                do j = 1, jmx-1
                  do i = 1, imx-1
                     T_L = z_pressure_left(i, j, 1) / (z_density_left(i, j, 1) * R_gas)  
                     T_R = z_pressure_right(i, j, 1) / (z_density_right(i, j, 1) * R_gas)
                     T_face = 0.5 * (T_L + T_R)
                     mu = mu_ref * (T_face / T_ref)**1.5 * ((T_ref + &    
                         Sutherland_temp) / (T_face + Sutherland_temp))
                    rho = 0.5*(z_density_left(i,j,1) + z_density_right(i,j,1))
                    z_tw_left(i, j, 1) = 60*mu/(rho*beta1*dist(i,j,1)**2)
                    z_tw_right(i, j, 1) = 60*mu/(rho*beta1*dist(i,j,1)**2)
                  end do
                end do
            else if ((bc_kmx(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
                z_tk_left(:, :, kmx) = 0.
                z_tk_right(:, :, kmx) = 0.
                do j = 1, jmx-1
                  do i = 1, imx-1
                     T_L = z_pressure_left(i, j, kmx) / (z_density_left(i, j, kmx) * R_gas)  
                     T_R = z_pressure_right(i, j, kmx) / (z_density_right(i, j, kmx) * R_gas)
                     T_face = 0.5 * (T_L + T_R)
                     mu = mu_ref * (T_face / T_ref)**1.5 * ((T_ref + &    
                         Sutherland_temp) / (T_face + Sutherland_temp))
                    rho = 0.5*(z_density_left(i,j,kmx) + z_density_right(i,j,kmx))
                    z_tw_left(i, j, kmx) = 60*mu/(rho*beta1*dist(i,j,kmx)**2)
                    z_tw_right(i, j, kmx) = 60*mu/(rho*beta1*dist(i,j,kmx)**2)
                  end do
                end do
            end if
