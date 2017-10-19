


    do k=1,kmx-1
      do j=1,jmx-1
        do i=1,imx-1
       gradw_z(i,j,k) =(-(z_speed(i-1,j  ,k  )+z_speed(i,j,k))*xnz(i,j,k)*xA(i,j,k) &
                        -(z_speed(i  ,j-1,k  )+z_speed(i,j,k))*ynz(i,j,k)*yA(i,j,k) &
                        -(z_speed(i  ,j  ,k-1)+z_speed(i,j,k))*znz(i,j,k)*zA(i,j,k) &
                        +(z_speed(i+1,j  ,k  )+z_speed(i,j,k))*xnz(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                        +(z_speed(i  ,j+1,k  )+z_speed(i,j,k))*ynz(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                        +(z_speed(i  ,j  ,k+1)+z_speed(i,j,k))*znz(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                       )/(2*volume(i,j,k))
        end do
      end do
    end do

