!should be included in subroutine - no_slip 
!    in module boundary conditons

  if (cell_ind == -1) then
      if (face == "imin") then
          tk(0, 1:jmx-1, 1:kmx-1) = - tk(1, 1:jmx-1, 1:kmx-1)
          mu_t(0, 1:jmx-1, 1:kmx-1) = - mu_t(1, 1:jmx-1, 1:kmx-1)
          do k = 1,kmx-1
            do j = 1,jmx-1
              T_l = pressure(0, j, k)/(density(0, j, k)*R_gas)
              T_r = pressure(1, j, k)/(density(1, j, k)*R_gas)
              T_face =  0.5 * ( T_l + T_r )
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(0, j, k) + density(1, j, k))
              tw(0, j, k) = 120*mu/(rho*beta1*(2*dist(1, j, k))**2) - tw(1, j, k)
            end do
          end do
      else if (face == "imax") then
          tk(imx, 1:jmx-1, 1:kmx-1) = - tk(imx-1, 1:jmx-1, 1:kmx-1)
          mu_t(imx, 1:jmx-1, 1:kmx-1) = - mu_t(imx-1, 1:jmx-1, 1:kmx-1)
          do k = 1,kmx-1
            do j = 1,jmx-1
              T_l = pressure(imx-1, j, k)/(density(imx-1, j, k)*R_gas)
              T_r = pressure(imx, j, k)/(density(imx, j, k)*R_gas)
              T_face =  0.5 * ( T_l + T_r )
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(imx-1, j, k) + density(imx, j, k))
              tw(imx, j, k) = 120*mu/(rho*beta1*(2*dist(imx-1, j, k))**2) - tw(imx-1, j, k)
            end do
          end do
      else if (face == "jmin") then
          tk(1:imx-1, 0, 1:kmx-1) = - tk(1:imx-1, 1, 1:kmx-1)
          mu_t(1:imx-1, 0, 1:kmx-1) = - mu_t(1:imx-1, 1, 1:kmx-1)
          do k = 1,kmx-1
            do i = 1,imx-1
              T_l = pressure(i, 0, k)/(density(i, 0, k)*R_gas)
              T_r = pressure(i, 1, k)/(density(i, 1, k)*R_gas)
              T_face =  0.5 * ( T_l + T_r )
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(i, 0, k) + density(i, 1, k))
              tw(i, 0, k) = 120*mu/(rho*beta1*(2*dist(i, 1, k))**2) - tw(i, 1, k)
            end do
          end do
      else if (face == "jmax") then
          tk(1:imx-1, jmx, 1:kmx-1) = - tk(1:imx-1, jmx-1, 1:kmx-1)
          mu_t(1:imx-1, jmx, 1:kmx-1) = - mu_t(1:imx-1, jmx-1, 1:kmx-1)
          do k = 1,kmx-1
            do i = 1,imx-1
              T_l = pressure(i, jmx-1, k)/(density(i, jmx-1, k)*R_gas)
              T_r = pressure(i, jmx, k)/(density(i, jmx, k)*R_gas)
              T_face =  0.5 * ( T_l + T_r )
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(i, jmx-1, k) + density(i, jmx, k))
              tw(i, jmx, k) = 120*mu/(rho*beta1*(2*dist(i, jmx-1, k))**2) - tw(i, jmx-1, k)
            end do
          end do
      else if (face == "kmin") then
          tk(1:imx-1, 1:jmx-1, 0) = - tk(1:imx-1, 1:jmx-1, 1)
          mu_t(1:imx-1, 1:jmx-1, 0) = - mu_t(1:imx-1, 1:jmx-1, 1)
          do j = 1,jmx-1
            do i = 1,imx-1
              T_l = pressure(i, j, 0)/(density(i, j, 0)*R_gas)
              T_r = pressure(i, j, 1)/(density(i, j, 1)*R_gas)
              T_face =  0.5 * ( T_l + T_r )
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(i, j, 0) + density(i, j, 1))
              tw(i, j, 0) = 120*mu/(rho*beta1*(2*dist(i, j, 1))**2) - tw(i, j, 1)
            end do
          end do
      else if (face == "kmax") then
          tk(1:imx-1, 1:jmx-1, kmx) = - tk(1:imx-1, 1:jmx-1, kmx-1)
          mu_t(1:imx-1, 1:jmx-1, kmx) = - mu_t(1:imx-1, 1:jmx-1, kmx-1)
          do j = 1,jmx-1
            do i = 1,imx-1
              T_l = pressure(i, j, kmx-1)/(density(i, j, kmx-1)*R_gas)
              T_r = pressure(i, j, kmx)/(density(i, j, kmx)*R_gas)
              T_face =  0.5 * ( T_l + T_r )
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(i, j, kmx-1) + density(i, j, kmx))
              tw(i, j, kmx) = 120*mu/(rho*beta1*(2*dist(i, j, kmx-1))**2) - tw(i, j, kmx-1)
            end do
          end do
      end if
  end if
