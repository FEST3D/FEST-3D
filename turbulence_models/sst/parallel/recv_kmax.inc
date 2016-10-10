! include in paralle module, send_recv_k subroutine

    do j = 0, jmx
      do i = 0, imx
        do l = 1, layers
          tk(i,j,kmx+l-1) = kmax_recv_buf(count) 
          count = count+1
        end do
      end do
    end do
    
    do j = 0, jmx
      do i = 0, imx
        do l = 1, layers
          tw(i,j,kmx+l-1) = kmax_recv_buf(count) 
          count = count+1
        end do
      end do
    end do
