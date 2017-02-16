! include in paralle module, send_recv_k subroutine

    do j = 0, jmx
      do i = 0, imx
        do l = 1, layers
          kmin_send_buf(count) = tk(i,j,l)
          count = count+1
        end do
      end do
    end do
    
    do j = 0, jmx
      do i = 0, imx
        do l = 1, layers
          kmin_send_buf(count) = tw(i,j,l)
          count = count+1
        end do
      end do
    end do
    
    do j = 0, jmx
      do i = 0, imx
        do l = 1, layers
          kmin_send_buf(count) = sst_F1(i,j,l)
          count = count+1
        end do
      end do
    end do
