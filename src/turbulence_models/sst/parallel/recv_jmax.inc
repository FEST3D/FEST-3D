! include in paralle module, send_recv_j subroutine

    do k = 0, kmx
      do i = 0, imx
        do l = 1, layers
          tk(i,jmx+l-1,k) = jmax_recv_buf(count) 
          count = count+1
        end do
      end do
    end do
    
    do k = 0, kmx
      do i = 0, imx
        do l = 1, layers
          tw(i,jmx+l-1,k) = jmax_recv_buf(count) 
          count = count+1
        end do
      end do
    end do
    
    do k = 0, kmx
      do i = 0, imx
        do l = 1, layers
          sst_F1(i,jmx+l-1,k) = jmax_recv_buf(count) 
          count = count+1
        end do
      end do
    end do
