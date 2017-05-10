!module boundary_conditions

        !---------------------------------------------------------------
        ! "Fix" type conditions
        !
        ! These subroutines are named "fix_<parameter>" where 
        ! <parameter> can be tk, tw. Each 
        ! of the subroutines accept an argument face as described above. 
        ! The subroutine will set the value of the ghost cell at the 
        ! face specified. 
        !
        ! The subroutine also accepts an optional val. If val is 
        ! specified, the parameter is fixed to that value, else the 
        ! infinity values are used. 
        !---------------------------------------------------------------

        subroutine fix_tk(face, cell_ind, optional_fixed_tk)

            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind
            real, intent(in), optional :: optional_fixed_tk
            real :: fixed_tk

            if (present(optional_fixed_tk)) then
                fixed_tk = optional_fixed_tk
            else
                fixed_tk = tk_inf
            end if
            
            if (cell_ind == -1) then
                if (face == "imin") then
                    tk(0, 1:jmx-1, 1:kmx-1) = fixed_tk
                else if (face == "imax") then
                    tk(imx, 1:jmx-1, 1:kmx-1) = fixed_tk
                else if (face == "jmin") then
                    tk(1:imx-1, 0, 1:kmx-1) = fixed_tk
                else if (face == "jmax") then
                    tk(1:imx-1, jmx, 1:kmx-1) = fixed_tk
                else if (face == "kmin") then
                    tk(1:imx-1, 1:jmx-1, 0) = fixed_tk
                else if (face == "kmax") then
                    tk(1:imx-1, 1:jmx-1, kmx) = fixed_tk
                end if
            end if

        end subroutine fix_tk

        subroutine fix_tw(face, cell_ind, optional_fixed_tw)

            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind
            real, intent(in), optional :: optional_fixed_tw
            real :: fixed_tw

            if (present(optional_fixed_tw)) then
                fixed_tw = optional_fixed_tw
            else
                fixed_tw = tw_inf
            end if
            
            if (cell_ind == -1) then
                if (face == "imin") then
                    tw(0, 1:jmx-1, 1:kmx-1) = fixed_tw
                else if (face == "imax") then
                    tw(imx, 1:imx-1, 1:kmx-1) = fixed_tw
                else if (face == "jmin") then
                    tw(1:imx-1, 0, 1:kmx-1) = fixed_tw
                else if (face == "jmax") then
                    tw(1:imx-1, jmx, 1:kmx-1) = fixed_tw
                else if (face == "kmin") then
                    tw(1:imx-1, 1:jmx-1, 0) = fixed_tw
                else if (face == "kmax") then
                    tw(1:imx-1, 1:jmx-1, kmx) = fixed_tw
                end if
            end if


        end subroutine fix_tw

        !---------------------------------------------------------------
        ! "Copy" type conditions
        !
        ! These subroutines are named "copy_<parameter>" where 
        ! <parameter> can be tk, tw. Each 
        ! of the subroutines accept an argument face as described above.
        ! The subroutine will set the value of the ghost cell at the 
        ! face specified. 
        !---------------------------------------------------------------

        subroutine copy_tk(face, cell_ind)

            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind

            if (cell_ind == -1) then
                if (face == "imin") then
                    tk(0, 1:jmx-1, 1:kmx-1) = tk(1, 1:jmx-1, 1:kmx-1)
                else if (face == "imax") then
                    tk(imx, 1:jmx-1, 1:kmx-1) = tk(imx-1, 1:jmx-1, 1:kmx-1)
                else if (face == "jmin") then
                    tk(1:imx-1, 0, 1:kmx-1) = tk(1:imx-1, 1, 1:kmx-1)
                else if (face == "jmax") then
                    tk(1:imx-1, jmx, 1:kmx-1) = tk(1:imx-1, jmx-1, 1:kmx-1)
                else if (face == "kmin") then
                    tk(1:imx-1, 1:jmx-1, 0) = tk(1:imx-1, 1:jmx-1, 1)
                else if (face == "kmax") then
                    tk(1:imx-1, 1:jmx-1, kmx) = tk(1:imx-1, 1:jmx-1, kmx-1)
                end if
            end if

        end subroutine copy_tk


        subroutine copy_tw(face, cell_ind)

            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind

            if (cell_ind == -1) then
                if (face == "imin") then
                    tw(0, 1:jmx-1, 1:kmx-1) = tw(1, 1:jmx-1, 1:kmx-1)
                else if (face == "imax") then
                    tw(imx, 1:jmx-1, 1:kmx-1) = tw(imx-1, 1:jmx-1, 1:kmx-1)
                else if (face == "jmin") then
                    tw(1:imx-1, 0, 1:kmx-1) = tw(1:imx-1, 1, 1:kmx-1)
                else if (face == "jmax") then
                    tw(1:imx-1, jmx, 1:kmx-1) = tw(1:imx-1, jmx-1, 1:kmx-1)
                else if (face == "kmin") then
                    tw(1:imx-1, 1:jmx-1, 0) = tw(1:imx-1, 1:jmx-1, 1)
                else if (face == "kmax") then
                    tw(1:imx-1, 1:jmx-1, kmx) = tw(1:imx-1, 1:jmx-1, kmx-1)
                end if
            end if

        end subroutine copy_tw

        subroutine extrapolate_F1(face)
          implicit none
          character(len=*), intent(in) :: face
          if (face == "imin") then
              sst_F1(0, 1:jmx-1, 1:kmx-1) = sst_F1(1, 1:jmx-1, 1:kmx-1)
          else if (face == "imax") then
              sst_F1(imx, 1:jmx-1, 1:kmx-1) = sst_F1(imx-1, 1:jmx-1, 1:kmx-1)
          else if (face == "jmin") then
              sst_F1(1:imx-1, 0, 1:kmx-1) = sst_F1(1:imx-1, 1, 1:kmx-1)
          else if (face == "jmax") then
              sst_F1(1:imx-1, jmx, 1:kmx-1) = sst_F1(1:imx-1, jmx-1, 1:kmx-1)
          else if (face == "kmin") then
              sst_F1(1:imx-1, 1:jmx-1, 0) = sst_F1(1:imx-1, 1:jmx-1, 1)
          else if (face == "kmax") then
              sst_F1(1:imx-1, 1:jmx-1, kmx) = sst_F1(1:imx-1, 1:jmx-1, kmx-1)
          end if
        end subroutine extrapolate_F1
