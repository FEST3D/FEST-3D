
    do k=1,kmx-1
      do j=1,jmx-1
        do i=1,imx-1
       gradv_y(i,j,k) =(-(y_speed(i-1,j  ,k  )+y_speed(i,j,k))*xny(i,j,k)*xA(i,j,k) &
                        -(y_speed(i  ,j-1,k  )+y_speed(i,j,k))*yny(i,j,k)*yA(i,j,k) &
                        -(y_speed(i  ,j  ,k-1)+y_speed(i,j,k))*zny(i,j,k)*zA(i,j,k) &
                        +(y_speed(i+1,j  ,k  )+y_speed(i,j,k))*xny(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                        +(y_speed(i  ,j+1,k  )+y_speed(i,j,k))*yny(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                        +(y_speed(i  ,j  ,k+1)+y_speed(i,j,k))*zny(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                       )/(2*volume(i,j,k))
        end do
      end do
    end do

