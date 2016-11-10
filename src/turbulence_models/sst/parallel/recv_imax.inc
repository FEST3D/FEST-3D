! include in paralle module, send_recv_i subroutine

    do k = 0, kmx
      do j = 0, jmx
        do l = 1, layers
          tk(imx+l-1,j,k) = imax_recv_buf(count)
          count = count+1
        end do
      end do
    end do
    
    do k = 0, kmx
      do j = 0, jmx
        do l = 1, layers
           tw(imx+l-1,j,k) = imax_recv_buf(count)
          count = count+1
        end do
      end do
    end do
