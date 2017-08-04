  do k = 1,kmx-1
    do j = 1,jmx-1
      do i = 1,imx-1

        ! get turbulent kinetic energy
        select case (trim(turbulence))
          case('sst', 'kkl')
            k = qp(i,j,k,6)
          case DEFAULT
            k = 0.
        end select

        r = qp(i,j,k,1) !density
        u = qp(i,j,k,2) !x_speed
        v = qp(i,j,k,3) !y_speed
        w = qp(i,j,k,4) !z_speed

        p = (gm-1.)*(qc(i,j,k,5) - (0.5*r*(u**2+v**2+w**2)) - r*k )
        qp(i,j,k,5) = p !pressure

      end do
    end do
  end do

