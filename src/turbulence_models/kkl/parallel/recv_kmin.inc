! include in paralle module, send_recv_k subroutine

    do j = 0, jmx
      do i = 0, imx
        do l = 1, layers
          tk(i,j,1-l) = kmin_recv_buf(count) 
          count = count+1
        end do
      end do
    end do
    
    do j = 0, jmx
      do i = 0, imx
        do l = 1, layers
          tkl(i,j,1-l) = kmin_recv_buf(count) 
          count = count+1
        end do
      end do
    end do
    
    do j = 0, jmx
      do i = 0, imx
        do l = 1, layers
          mu_t(i,j,1-l) = kmin_recv_buf(count) 
          count = count+1
        end do
      end do
    end do
