! include in paralle module, send_recv_j subroutine

    do k = 0, kmx
      do i = 0, imx
        do l = 1, layers
          jmax_send_buf(count) = tk(i,jmx-l,k)
          count = count+1
        end do
      end do
    end do
    
    do k = 0, kmx
      do i = 0, imx
        do l = 1, layers
          jmax_send_buf(count) = tw(i,jmx-l,k)
          count = count+1
        end do
      end do
    end do
    
    do k = 0, kmx
      do i = 0, imx
        do l = 1, layers
          jmax_send_buf(count) = sst_F1(i,jmx-l,k)
          count = count+1
        end do
      end do
    end do
