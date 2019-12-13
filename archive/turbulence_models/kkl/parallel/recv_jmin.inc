! include in paralle module, send_recv_j subroutine

    do k = 0, kmx
      do i = 0, imx
        do l = 1, layers
          tk(i,1-l,k) = jmin_recv_buf(count) 
          count = count+1
        end do
      end do
    end do
    
    do k = 0, kmx
      do i = 0, imx
        do l = 1, layers
          tkl(i,1-l,k) = jmin_recv_buf(count) 
          count = count+1
        end do
      end do
    end do
    
    do k = 0, kmx
      do i = 0, imx
        do l = 1, layers
          mu_t(i,1-l,k) = jmin_recv_buf(count) 
          count = count+1
        end do
      end do
    end do
