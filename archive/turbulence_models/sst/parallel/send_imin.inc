! include in paralle module, send_recv_i subroutine

    do k = 0, kmx
      do j = 0, jmx
        do l = 1, layers
          imin_send_buf(count) = tk(l,j,k)
          count = count+1
        end do
      end do
    end do
    
    do k = 0, kmx
      do j = 0, jmx
        do l = 1, layers
          imin_send_buf(count) = tw(l,j,k)
          count = count+1
        end do
      end do
    end do
    
    do k = 0, kmx
      do j = 0, jmx
        do l = 1, layers
          imin_send_buf(count) = sst_F1(l,j,k)
          count = count+1
        end do
      end do
    end do
    
    do k = 0, kmx
      do j = 0, jmx
        do l = 1, layers
          imin_send_buf(count) = mu_t(l,j,k)
          count = count+1
        end do
      end do
    end do
