
        select case(id(1))

          case(-1)
              density(      0, 1:jmx-1, 1:kmx-1) = fix_density
              density(     -1, 1:jmx-1, 1:kmx-1) = fix_density
              density(     -2, 1:jmx-1, 1:kmx-1) = fix_density
              x_speed(      0, 1:jmx-1, 1:kmx-1) = fix_x_speed
              x_speed(     -1, 1:jmx-1, 1:kmx-1) = fix_x_speed
              x_speed(     -2, 1:jmx-1, 1:kmx-1) = fix_x_speed
              density(      0, 1:jmx-1, 1:kmx-1) = fix_density
              density(     -1, 1:jmx-1, 1:kmx-1) = fix_density
              density(     -2, 1:jmx-1, 1:kmx-1) = fix_density
            call fix(x_speed , fixed_x_speed , face)
            call fix(y_speed , fixed_y_speed , face)
            call fix(z_speed , fixed_z_speed , face)
            call fix(pressure, fixed_pressure, face)
            select case (turbulence)
              case('none')
                !do nothing
                continue
              case('sst')
                call check_if_value_fixed("sst")
                call fix(tk, fixed_tk, face)
                call fix(tw, fixed_tw, face)
              case DEFAULT
                call turbulence_read_error()
            end select

          case(-2)
        call copy3(density , "flat", face)
        call copy3(x_speed , "flat", face)
        call copy3(y_speed , "flat", face)
        call copy3(z_speed , "flat", face)
        call copy3(pressure, "flat", face)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sst')
            call copy3(tk, "flat", face)
            call copy3(tw, "flat", face)
          case DEFAULT
            call turbulence_read_error()
        end select

          case(-3)
        call fix(density , fixed_density , face)
        call fix(x_speed , fixed_x_speed , face)
        call fix(y_speed , fixed_y_speed , face)
        call fix(z_speed , fixed_z_speed , face)
        call copy3(pressure, "flat", face)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sst')
            call check_if_value_fixed("sst")
            call fix(tk, fixed_tk, face)
            call fix(tw, fixed_tw, face)
          case DEFAULT
            call turbulence_read_error()
        end select

          case(-4)
        call copy3(density, "flat", face)
        call copy3(x_speed, "flat", face)
        call copy3(y_speed, "flat", face)
        call copy3(z_speed, "flat", face)
        call fix(pressure, fixed_pressure, face)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sst')
            call copy3(tk, "flat", face)
            call copy3(tw, "flat", face)
          case DEFAULT
            call turbulence_read_error()
        end select

          case(-5)
        call copy3(density , "symm",  face)
        call copy3(pressure, "symm",  face)
        call no_slip(face)

          case(-6)
        call copy3(density , "symm", face)
        call copy3(pressure, "symm", face)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sst')
            call copy3(tk, "symm", face)
            call copy3(tw, "symm", face)
          case DEFAULT
            call turbulence_read_error()
        end select
        call flow_tangency(face)

          case Default
            if(id(i)>=0) then
              continue !interface boundary 
            else
              print*, " boundary condition not recognised -> id is :", id(i)
            end if

          end select


              var(      0, 1:jmx-1, 1:kmx-1) = fix_val(1)
              var(     -1, 1:jmx-1, 1:kmx-1) = fix_val(1)
              var(     -2, 1:jmx-1, 1:kmx-1) = fix_val(1)


      subroutine no_slip(face)
        implicit none
        character(len=*), intent(in) :: face
        call copy3(x_speed, "anti", face)
        call copy3(y_speed, "anti", face)
        call copy3(z_speed, "anti", face)
        select case(turbulence)
          case("none")
            !do nothing
            continue
          case("sst")
            call copy3(tk  , "anti", face)
            call set_omega_at_wall(face)
          case DEFAULT
            call turbulence_read_error()
        end select
      end subroutine no_slip

    

      subroutine set_omega_at_wall(face)
        implicit none
        character(len=*), intent(in) :: face
        real :: T_face
        real :: mu
        real :: rho
        integer :: i,j,k
        
        select case(face)
          case("imin")
          do k = 1,kmx-1
            do j = 1,jmx-1
              T_face = 0.5*((pressure(0, j, k)/density(0, j, k))+(pressure(1, j, k)/density(1, j, k)))/R_gas
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(0, j, k) + density(1, j, k))
              tw(0, j, k) = 120*mu/(rho*beta1*(2*dist(1, j, k))**2) - tw(1, j, k)
            end do
          end do
        case("imax")
          do k = 1,kmx-1
            do j = 1,jmx-1
              T_face = 0.5*((pressure(imx-1, j, k)/density(imx-1, j, k))+(pressure(imx, j, k)/density(imx, j, k)))/R_gas
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(imx-1, j, k) + density(imx, j, k))
              tw(imx, j, k) = 120*mu/(rho*beta1*(2*dist(imx-1, j, k))**2) - tw(imx-1, j, k)
            end do
          end do
        case("jmin")
          do k = 1,kmx-1
            do i = 1,imx-1
              T_face = 0.5*((pressure(i, 0, k)/density(i, 0, k))+(pressure(i, 1, k)/density(i, 1, k)))/R_gas
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(i, 0, k) + density(i, 1, k))
              tw(i, 0, k) = 120*mu/(rho*beta1*(2*dist(i, 1, k))**2) - tw(i, 1, k)
            end do
          end do
        case("jmax")
          do k = 1,kmx-1
            do i = 1,imx-1
              T_face = 0.5*((pressure(i, jmx-1, k)/density(i, jmx-1, k))+(pressure(i, jmx, k)/density(i, jmx, k)))/R_gas
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(i, jmx-1, k) + density(i, jmx, k))
              tw(i, jmx, k) = 120*mu/(rho*beta1*(2*dist(i, jmx-1, k))**2) - tw(i, jmx-1, k)
            end do
          end do
        case("kmin")
          do j = 1,jmx-1
            do i = 1,imx-1
              T_face = 0.5*((pressure(i, j, 0)/density(i, j, 0))+(pressure(i, j, 1)/density(i, j, 1)))/R_gas
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(i, j, 0) + density(i, j, 1))
              tw(i, j, 0) = 120*mu/(rho*beta1*(2*dist(i, j, 1))**2) - tw(i, j, 1)
            end do
          end do
        case("kmax")
          do j = 1,jmx-1
            do i = 1,imx-1
              T_face = 0.5*((pressure(i, j, kmx-1)/density(i, j, kmx-1))+(pressure(i, j, kmx)/density(i, j, kmx)))/R_gas
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(i, j, kmx-1) + density(i, j, kmx))
              tw(i, j, kmx) = 120*mu/(rho*beta1*(2*dist(i, j, kmx-1))**2) - tw(i, j, kmx-1)
            end do
          end do

      end select
    end subroutine set_omega_at_wall

    subroutine check_if_value_fixed(model)
      implicit none
      character(len=*), intent(in) :: model

      select case(model)
        case("none")
          !do nothing
          continue
        case("sst")
          if(fixed_tk(face_num)==0.) fixed_tk(face_num)=tk_inf
          if(fixed_tw(face_num)==0.) fixed_tw(face_num)=tw_inf
        case DEFAULT
          call turbulence_read_error()
      end select
    end subroutine check_if_value_fixed

end module bc_primitive
