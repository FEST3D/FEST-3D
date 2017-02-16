! include in paralle module, send_recv_i subroutine

    do k = 0, kmx
      do j = 0, jmx
        do l = 1, layers
          tk(1-l,j,k) = imin_recv_buf(count)
          count = count+1
        end do
      end do
    end do
    
    do k = 0, kmx
      do j = 0, jmx
        do l = 1, layers
           tw(1-l,j,k) = imin_recv_buf(count)
          count = count+1
        end do
      end do
    end do

    do k = 0, kmx
      do j = 0, jmx
        do l = 1, layers
           sst_F1(1-l,j,k) = imin_recv_buf(count)
          count = count+1
        end do
      end do
    end do
