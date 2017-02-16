! include in paralle module, send_recv_i subroutine

    do k = 0, kmx
      do j = 0, jmx
        do l = 1, layers
          imax_send_buf(count) = tk(imx-l,j,k)
          count = count+1
        end do
      end do
    end do
    
    do k = 0, kmx
      do j = 0, jmx
        do l = 1, layers
          imax_send_buf(count) = tw(imx-l,j,k)
          count = count+1
        end do
      end do
    end do
    
    do k = 0, kmx
      do j = 0, jmx
        do l = 1, layers
          imax_send_buf(count) = sst_F1(imx-l,j,k)
          count = count+1
        end do
      end do
    end do
