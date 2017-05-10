!should be included in subroutine - no_slip 
!    in module boundary conditons

  if (cell_ind == -1) then
      if (face == "imin") then
            tk(0, 1:jmx-1, 1:kmx-1) =    tk(1, 1:jmx-1, 1:kmx-1)
            tw(0, 1:jmx-1, 1:kmx-1) =    tw(1, 1:jmx-1, 1:kmx-1)
          mu_t(0, 1:jmx-1, 1:kmx-1) =  mu_t(1, 1:jmx-1, 1:kmx-1)

            tk(-1, 1:jmx-1, 1:kmx-1) =    tk(2, 1:jmx-1, 1:kmx-1)
            tw(-1, 1:jmx-1, 1:kmx-1) =    tw(2, 1:jmx-1, 1:kmx-1)
          mu_t(-1, 1:jmx-1, 1:kmx-1) =  mu_t(2, 1:jmx-1, 1:kmx-1)

            tk(-2, 1:jmx-1, 1:kmx-1) =    tk(3, 1:jmx-1, 1:kmx-1)
            tw(-2, 1:jmx-1, 1:kmx-1) =    tw(3, 1:jmx-1, 1:kmx-1)
          mu_t(-2, 1:jmx-1, 1:kmx-1) =  mu_t(3, 1:jmx-1, 1:kmx-1)

      else if (face == "imax") then
            tk(imx, 1:jmx-1, 1:kmx-1) =    tk(imx-1, 1:jmx-1, 1:kmx-1)
            tw(imx, 1:jmx-1, 1:kmx-1) =    tw(imx-1, 1:jmx-1, 1:kmx-1)
          mu_t(imx, 1:jmx-1, 1:kmx-1) =  mu_t(imx-1, 1:jmx-1, 1:kmx-1)

            tk(imx+1, 1:jmx-1, 1:kmx-1) =    tk(imx-2, 1:jmx-1, 1:kmx-1)
            tw(imx+1, 1:jmx-1, 1:kmx-1) =    tw(imx-2, 1:jmx-1, 1:kmx-1)
          mu_t(imx+1, 1:jmx-1, 1:kmx-1) =  mu_t(imx-2, 1:jmx-1, 1:kmx-1)

            tk(imx+2, 1:jmx-1, 1:kmx-1) =    tk(imx-3, 1:jmx-1, 1:kmx-1)
            tw(imx+2, 1:jmx-1, 1:kmx-1) =    tw(imx-3, 1:jmx-1, 1:kmx-1)
          mu_t(imx+2, 1:jmx-1, 1:kmx-1) =  mu_t(imx-3, 1:jmx-1, 1:kmx-1)

      else if (face == "jmin") then
            tk(1:imx-1, 0, 1:kmx-1) =    tk(1:imx-1, 1, 1:kmx-1)
            tw(1:imx-1, 0, 1:kmx-1) =    tw(1:imx-1, 1, 1:kmx-1)
          mu_t(1:imx-1, 0, 1:kmx-1) =  mu_t(1:imx-1, 1, 1:kmx-1)

            tk(1:imx-1, -1, 1:kmx-1) =    tk(1:imx-1, 2, 1:kmx-1)
            tw(1:imx-1, -1, 1:kmx-1) =    tw(1:imx-1, 2, 1:kmx-1)
          mu_t(1:imx-1, -1, 1:kmx-1) =  mu_t(1:imx-1, 2, 1:kmx-1)

            tk(1:imx-1, -2, 1:kmx-1) =    tk(1:imx-1, 3, 1:kmx-1)
            tw(1:imx-1, -2, 1:kmx-1) =    tw(1:imx-1, 3, 1:kmx-1)
          mu_t(1:imx-1, -2, 1:kmx-1) =  mu_t(1:imx-1, 3, 1:kmx-1)

      else if (face == "jmax") then
            tk(1:imx-1, jmx, 1:kmx-1) =    tk(1:imx-1, jmx-1, 1:kmx-1)
            tw(1:imx-1, jmx, 1:kmx-1) =    tw(1:imx-1, jmx-1, 1:kmx-1)
          mu_t(1:imx-1, jmx, 1:kmx-1) =  mu_t(1:imx-1, jmx-1, 1:kmx-1)

            tk(1:imx-1, jmx+1, 1:kmx-1) =    tk(1:imx-1, jmx-2, 1:kmx-1)
            tw(1:imx-1, jmx+1, 1:kmx-1) =    tw(1:imx-1, jmx-2, 1:kmx-1)
          mu_t(1:imx-1, jmx+1, 1:kmx-1) =  mu_t(1:imx-1, jmx-2, 1:kmx-1)

            tk(1:imx-1, jmx+2, 1:kmx-1) =    tk(1:imx-1, jmx-3, 1:kmx-1)
            tw(1:imx-1, jmx+2, 1:kmx-1) =    tw(1:imx-1, jmx-3, 1:kmx-1)
          mu_t(1:imx-1, jmx+2, 1:kmx-1) =  mu_t(1:imx-1, jmx-3, 1:kmx-1)

      else if (face == "kmin") then
            tk(1:imx-1, 1:jmx-1, 0) =    tk(1:imx-1, 1:jmx-1, 1)
            tw(1:imx-1, 1:jmx-1, 0) =    tw(1:imx-1, 1:jmx-1, 1)
          mu_t(1:imx-1, 1:jmx-1, 0) =  mu_t(1:imx-1, 1:jmx-1, 1)

            tk(1:imx-1, 1:jmx-1, -1) =    tk(1:imx-1, 1:jmx-1, 2)
            tw(1:imx-1, 1:jmx-1, -1) =    tw(1:imx-1, 1:jmx-1, 2)
          mu_t(1:imx-1, 1:jmx-1, -1) =  mu_t(1:imx-1, 1:jmx-1, 2)

            tk(1:imx-1, 1:jmx-1, -2) =    tk(1:imx-1, 1:jmx-1, 3)
            tw(1:imx-1, 1:jmx-1, -2) =    tw(1:imx-1, 1:jmx-1, 3)
          mu_t(1:imx-1, 1:jmx-1, -2) =  mu_t(1:imx-1, 1:jmx-1, 3)

      else if (face == "kmax") then
            tk(1:imx-1, 1:jmx-1, kmx) =    tk(1:imx-1, 1:jmx-1, kmx-1)
            tw(1:imx-1, 1:jmx-1, kmx) =    tw(1:imx-1, 1:jmx-1, kmx-1)
          mu_t(1:imx-1, 1:jmx-1, kmx) =  mu_t(1:imx-1, 1:jmx-1, kmx-1)

            tk(1:imx-1, 1:jmx-1, kmx+1) =    tk(1:imx-1, 1:jmx-1, kmx-2)
            tw(1:imx-1, 1:jmx-1, kmx+1) =    tw(1:imx-1, 1:jmx-1, kmx-2)
          mu_t(1:imx-1, 1:jmx-1, kmx+1) =  mu_t(1:imx-1, 1:jmx-1, kmx-2)

            tk(1:imx-1, 1:jmx-1, kmx+2) =    tk(1:imx-1, 1:jmx-1, kmx-3)
            tw(1:imx-1, 1:jmx-1, kmx+2) =    tw(1:imx-1, 1:jmx-1, kmx-3)
          mu_t(1:imx-1, 1:jmx-1, kmx+2) =  mu_t(1:imx-1, 1:jmx-1, kmx-3)

      end if
  end if
