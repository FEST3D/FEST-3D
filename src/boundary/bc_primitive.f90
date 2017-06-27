module bc_primitive
  !--------------------------------------------
  ! 170515  Jatinder Pal Singh Sandhu
  ! Aim : applying boundary condition to domain
  !-------------------------------------------
#include "../error.inc"

  use global_vars, only: imin_id
  use global_vars, only: imax_id
  use global_vars, only: jmin_id
  use global_vars, only: jmax_id
  use global_vars, only: kmin_id
  use global_vars, only: kmax_id
  use global_vars, only: fixed_density
  use global_vars, only: fixed_x_speed
  use global_vars, only: fixed_y_speed
  use global_vars, only: fixed_z_speed
  use global_vars, only: fixed_pressure
  use global_vars, only: fixed_tk
  use global_vars, only: fixed_tw
  use global_vars, only: fixed_tkl
  use global_vars, only: mu_ref
  use global_vars, only: R_gas
  use global_vars, only: sutherland_temp
  use global_vars, only: pressure
  use global_vars, only: density
  use global_vars, only: x_speed
  use global_vars, only: y_speed
  use global_vars, only: z_speed
  use global_vars, only: tk
  use global_vars, only: tw
  use global_vars, only: tkl
  use global_vars, only: accur
  use global_vars, only: imx
  use global_vars, only: jmx
  use global_vars, only: kmx
  use global_vars, only: turbulence
  use global_vars, only: xnx
  use global_vars, only: xny
  use global_vars, only: xnz
  use global_vars, only: ynx
  use global_vars, only: yny
  use global_vars, only: ynz
  use global_vars, only: znx
  use global_vars, only: zny
  use global_vars, only: znz
  use global_vars, only: mu_t
  use global_vars, only: T_ref
  use global_vars, only: dist
  use global_vars, only: process_id
  use global_vars, only: tk_inf
  use global_vars, only: tw_inf
  use global_vars, only: te_inf
  use global_vars, only: tv_inf
  use global_vars, only: tkl_inf
  use global_vars, only: face_names
  use global_vars, only: id

  use global_sst , only: beta1
  use utils,       only: turbulence_read_error

  use read_bc   , only : read_fixed_values
  use copy_bc   , only : copy3
  use FT_bc     , only : flow_tangency

  implicit none
  private

  integer                        :: face_num

  public :: populate_ghost_primitive


  contains

    subroutine populate_ghost_primitive()
      implicit none
      integer :: i
      character(len=4) :: face

      
      do i = 1,6
        face_num = i
        face = face_names(face_num)

        select case(id(face_num))

          case(-1)
            call supersonic_inlet(face)

          case(-2)
            call supersonic_outlet(face)

          case(-3)
            call pressure_inlet(face)

          case(-4)
            call pressure_outlet(face)

          case(-5)
            call adiabatic_wall(face)

          case(-6)
            call slip_wall(face)

          case Default
            if(id(i)>=0) then
              continue !interface boundary 
            else
              print*, " boundary condition not recognised -> id is :", id(i)
            end if

          end select
        end do
      end subroutine populate_ghost_primitive


      subroutine supersonic_inlet(face)
        implicit none
        character(len=*), intent(in) :: face
        call fix(density , fixed_density , face)
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
          case('kkl')
            call check_if_value_fixed("kkl")
            call fix(tk, fixed_tk, face)
            call fix(tkl, fixed_tkl, face)
          case DEFAULT
            !call turbulence_read_error()
            Error
        end select
      end subroutine supersonic_inlet


      subroutine supersonic_outlet(face)
        implicit none
        character(len=*), intent(in) :: face
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
          case('kkl')
            call copy3(tk, "flat", face)
            call copy3(tkl, "flat", face)
          case DEFAULT
            !call turbulence_read_error()
            Error
        end select
      end subroutine supersonic_outlet


      subroutine pressure_inlet(face)
        implicit none
        character(len=*), intent(in) :: face
        call fix(density , fixed_density , face)
        call fix(x_speed , fixed_x_speed , face)
        call fix(y_speed , fixed_y_speed , face)
        call fix(z_speed , fixed_z_speed , face)
        call copy3(pressure, "flat", face)
        if(face=='imin') density(0,:,:)=pressure(1,:,:)/(R_gas*300.)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sst')
            call check_if_value_fixed("sst")
            call fix(tk, fixed_tk, face)
            call fix(tw, fixed_tw, face)
          case('kkl')
            call check_if_value_fixed("kkl")
            call fix(tk, fixed_tk, face)
            call fix(tkl, fixed_tw, face)
          case DEFAULT
           ! call turbulence_read_error()
           Error
        end select
      end subroutine pressure_inlet


      subroutine pressure_outlet(face)
        implicit none
        character(len=*), intent(in) :: face
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
          case('kkl')
            call copy3(tk, "flat", face)
            call copy3(tkl, "flat", face)
          case DEFAULT
           ! call turbulence_read_error()
           Error
        end select
      end subroutine pressure_outlet

      subroutine adiabatic_wall(face)
        implicit none
        character(len=*), intent(in) :: face
        call copy3(density , "symm",  face)
        call copy3(pressure, "symm",  face)
        call no_slip(face)
      end subroutine adiabatic_wall


      subroutine slip_wall(face)
        implicit none
        character(len=*), intent(in) :: face
        call copy3(density , "symm", face)
        call copy3(pressure, "symm", face)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sst')
            call copy3(tk, "symm", face)
            call copy3(tw, "symm", face)
          case('kkl')
            call copy3(tk, "symm", face)
            call copy3(tkl, "symm", face)
          case DEFAULT
            !call turbulence_read_error()
            Error
        end select
        call flow_tangency(face)
      end subroutine slip_wall


      subroutine fix(var, fix_val, face)
        implicit none
        real, dimension(-2:imx+2, -2:jmx+2, -2:kmx+2) , intent(out) :: var
        real, dimension(1:6)       , intent(in)  :: fix_val
        character(len=*)         , intent(in)  :: face

        select case(face)
          case("imin")
              var(      0, 1:jmx-1, 1:kmx-1) = fix_val(1)
              var(     -1, 1:jmx-1, 1:kmx-1) = fix_val(1)
              var(     -2, 1:jmx-1, 1:kmx-1) = fix_val(1)
           case("imax")
              var(  imx  , 1:jmx-1, 1:kmx-1) = fix_val(2)
              var(  imx+1, 1:jmx-1, 1:kmx-1) = fix_val(2)
              var(  imx+2, 1:jmx-1, 1:kmx-1) = fix_val(2)
          case("jmin")
              var(1:imx-1,       0, 1:kmx-1) = fix_val(3)
              var(1:imx-1,      -1, 1:kmx-1) = fix_val(3)
              var(1:imx-1,      -2, 1:kmx-1) = fix_val(3)
          case("jmax")
              var(1:imx-1,   jmx  , 1:kmx-1) = fix_val(4)
              var(1:imx-1,   jmx+1, 1:kmx-1) = fix_val(4)
              var(1:imx-1,   jmx+2, 1:kmx-1) = fix_val(4)
          case("kmin")
              var(1:imx-1, 1:jmx-1,       0) = fix_val(5)
              var(1:imx-1, 1:jmx-1,      -1) = fix_val(5)
              var(1:imx-1, 1:jmx-1,      -2) = fix_val(5)
          case("kmax")
              var(1:imx-1, 1:jmx-1,   kmx  ) = fix_val(6)
              var(1:imx-1, 1:jmx-1,   kmx+1) = fix_val(6)
              var(1:imx-1, 1:jmx-1,   kmx+2) = fix_val(6)
          case DEFAULT
            !print*, "ERROR: wrong face for boundary condition"
            Error
        end select
            
      end subroutine fix


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
          case("kkl")
            call copy3(tk  , "anti", face)
            call copy3(tkl , "anti", face)
          case DEFAULT
            !call turbulence_read_error()
            Error
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
        case("kkl")
          if(fixed_tk(face_num)==0.) fixed_tk(face_num)=tk_inf
          if(fixed_tkl(face_num)==0.) fixed_tkl(face_num)=tkl_inf
        case DEFAULT
         ! call turbulence_read_error()
         Error
      end select
    end subroutine check_if_value_fixed

end module bc_primitive
