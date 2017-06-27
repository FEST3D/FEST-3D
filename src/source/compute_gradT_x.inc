    do k=1,kmx-1
      do j=1,jmx-1
        do i=1,imx-1

          cell_T = (pressure(i,j,k)/density(i,j,k))/R_gas

          T(1)   = (pressure(i-1,j,k)/density(i-1,j,k))/R_gas + cell_T
          T(2)   = (pressure(i,j-1,k)/density(i,j-1,k))/R_gas + cell_T
          T(3)   = (pressure(i,j,k-1)/density(i,j,k-1))/R_gas + cell_T
          T(4)   = (pressure(i+1,j,k)/density(i+1,j,k))/R_gas + cell_T
          T(5)   = (pressure(i,j+1,k)/density(i,j+1,k))/R_gas + cell_T
          T(6)   = (pressure(i,j,k+1)/density(i,j,k+1))/R_gas + cell_T

       gradT_x(i,j,k) =(-T(1)*xnx(i,j,k)*xA(i,j,k) &
                        -T(2)*ynx(i,j,k)*yA(i,j,k) &
                        -T(3)*znx(i,j,k)*zA(i,j,k) &
                        +T(4)*xnx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                        +T(5)*ynx(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                        +T(6)*znx(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                       )/(2*volume(i,j,k))
        end do
      end do
    end do
