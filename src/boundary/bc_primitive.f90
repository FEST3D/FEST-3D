  !< Apply boundary condition at every iteration
module bc_primitive
  !< Apply boundary condition at every iteration
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
  use global_vars, only: fixed_tgm
  use global_vars, only: fixed_wall_temperature
  use global_vars, only: fixed_Ttemperature
  use global_vars, only: fixed_Tpressure
  use global_vars, only: fixed_tv
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
  use global_vars, only: tv
  use global_vars, only: tgm
  use global_vars, only: accur
  use global_vars, only: imx
  use global_vars, only: jmx
  use global_vars, only: kmx
  use global_vars, only: turbulence
  use global_vars, only: transition
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
  use global_vars, only: tgm_inf
  use global_vars, only: tkl_inf
  use global_vars, only: face_names
  use global_vars, only: id

  use global_vars, only: gm
  use global_vars, only: x_speed_inf
  use global_vars, only: y_speed_inf
  use global_vars, only: z_speed_inf
  use global_vars, only: density_inf
  use global_vars, only: pressure_inf
  use global_vars, only: vel_mag
  use global_vars, only: qp
  use global_vars, only: current_iter

  use global_sst , only: beta1
  use utils,       only: turbulence_read_error

  use read_bc   , only : read_fixed_values
  use copy_bc   , only : copy3
  use FT_bc     , only : flow_tangency

  implicit none
  private

  integer                        :: face_num
  !< Number of the face : 1:imin, 2:imax, 3:jmin, 4:jmax, 5:kmin, 6:kmax

  public :: populate_ghost_primitive


  contains

    subroutine populate_ghost_primitive()
      !< Populate the state variables in the ghost cell
      !< with particular value based on the boundary conditio 
      !< being applied at that face
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
            call subsonic_inlet(face)

          case(-4)
            call subsonic_outlet(face)

          case(-5)
            call wall(face)

          case(-6)
            call slip_wall(face)

          case(-7)
            call pole(face)

          case(-8)
            call far_field(face)

          case(-9)
            call periodic_bc(face)

          case(-11)
            call total_pressure(face)

          case Default
            if(id(i)>=0 .or. id(i)==-10) then
              continue !interface boundary 
            else
              print*, " boundary condition not recognised -> id is :", id(i)
            end if

          end select
        end do
!        qp(0,0,:,:) = 0.5*(qp(0,1,:,:)+qp(1,0,:,:))
!        qp(0,jmx,:,:) = 0.5*(qp(0,jmx-1,:,:)+qp(1,jmx,:,:))
!        qp(imx,0,:,:) = 0.5*(qp(imx,1,:,:)+qp(imx-1,0,:,:))
!        qp(imx,jmx,:,:) = 0.5*(qp(imx,jmx-1,:,:)+qp(imx-1,jmx,:,:))
!        qp(0,:,0,:) = 0.5*(qp(0,:,1,:)+qp(1,:,0,:))
!        qp(0,:,kmx,:) = 0.5*(qp(0,:,kmx-1,:)+qp(1,:,kmx,:))
!        qp(imx,:,0,:) = 0.5*(qp(imx,:,1,:)+qp(imx-1,:,0,:))
!        qp(imx,:,kmx,:) = 0.5*(qp(imx,:,jmx-1,:)+qp(imx-1,:,jmx,:))
         qp(:,0,0,:) = 0.33*(qp(:,1,1,:)+qp(:,0,1,:)+qp(:,1,0,:))
         qp(:,0,kmx,:) = 0.33*(qp(:,1,kmx-1,:)+qp(:,0,kmx-1,:)+qp(:,1,kmx,:))
         qp(:,jmx,0,:)   = 0.33*(qp(:,jmx-1,1,:)+qp(:,jmx,1,:)+qp(:,jmx-1,0,:))
         qp(:,jmx,kmx,:) = 0.33*(qp(:,jmx-1,kmx-1,:)+qp(:,jmx,kmx-1,:)+qp(:,jmx-1,kmx,:))
         qp(imx,0,:,:) = 0.33*(qp(imx-1,1,:,:)+qp(imx-1,0,:,:)+qp(imx,1,:,:))
         qp(0,0,:,:) = 0.33*(qp(1,1,:,:)+qp(1,0,:,:)+qp(0,1,:,:))
         qp(0,jmx,:,:)   = 0.33*(qp(1,jmx-1,:,:)+qp(1,jmx,:,:)+qp(0,jmx-1,:,:))
         qp(imx,jmx,:,:) = 0.33*(qp(imx-1,jmx-1,:,:)+qp(imx-1,jmx,:,:)+qp(imx,jmx-1,:,:))
         qp(0,0,0,:) = 0.33*(qp(1,0,0,:) + qp(0,1,0,:) + qp(0,0,1,:))
         qp(imx,0,0,:) = 0.33*(qp(imx-1,0,0,:) + qp(imx,1,0,:) + qp(imx,0,1,:))
         qp(0,jmx,0,:) = 0.33*(qp(1,jmx,0,:) + qp(0,jmx-1,0,:) + qp(0,jmx,1,:))
         qp(0,0,kmx,:) = 0.33*(qp(1,0,kmx,:) + qp(0,1,kmx,:) + qp(0,0,kmx-1,:))
         qp(imx,jmx,0,:) = 0.33*(qp(imx-1,jmx,0,:) + qp(imx,jmx-1,0,:) + qp(imx,jmx,1,:))
         qp(imx,0,kmx,:) = 0.33*(qp(imx-1,0,kmx,:) + qp(imx,1,kmx,:) + qp(imx,0,kmx-1,:))
         qp(0,jmx,kmx,:) = 0.33*(qp(1,jmx,kmx,:) + qp(0,jmx-1,kmx,:) + qp(0,jmx,kmx-1,:))
         qp(imx,jmx,kmx,:) = 0.33*(qp(imx-1,jmx,kmx,:) + qp(imx,jmx-1,kmx,:) + qp(imx,jmx,kmx-1,:))

      end subroutine populate_ghost_primitive


      subroutine supersonic_inlet(face)
        !< Supersonic inlet boundary condition
        !< All the values of state variables are fixed
        implicit none
        character(len=*), intent(in) :: face
        !< Name of the face at which boundary condition is called
        if(current_iter<=2)then
        call fix(density , fixed_density , face)
        call fix(x_speed , fixed_x_speed , face)
        call fix(y_speed , fixed_y_speed , face)
        call fix(z_speed , fixed_z_speed , face)
        call fix(pressure, fixed_pressure, face)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sa', 'saBC')
            call fix(tv, fixed_tv, face)
          case('sst', 'sst2003')
            call check_if_value_fixed("sst")
            call fix(tk, fixed_tk, face)
            call fix(tw, fixed_tw, face)
          case('kkl')
            call check_if_value_fixed("kkl")
            call fix(tk, fixed_tk, face)
            call fix(tkl, fixed_tkl, face)
          case DEFAULT
            !call turbulence_read_error()
            Fatal_error
        end select
        select case(trim(transition))
          case('lctm2015')
            call check_if_value_fixed("lctm2015")
            call fix(tgm, fixed_tgm, face)
          case DEFAULT
            continue
        end select
        end if
      end subroutine supersonic_inlet


      subroutine supersonic_outlet(face)
        !< Supersonic outlet boundary condition. 
        !< All the values of state variables are copied 
        !< from inside the domain
        implicit none
        character(len=*), intent(in) :: face
        !< Name of the face at which boundary condition is called
        call copy3(density , "flat", face)
        call copy3(x_speed , "flat", face)
        call copy3(y_speed , "flat", face)
        call copy3(z_speed , "flat", face)
        call copy3(pressure, "flat", face)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sa', 'saBC')
            call copy3(tv, "flat", face)
          case('sst', 'sst2003')
            call copy3(tk, "flat", face)
            call copy3(tw, "flat", face)
          case('kkl')
            call copy3(tk, "flat", face)
            call copy3(tkl, "flat", face)
          case DEFAULT
            !call turbulence_read_error()
            Fatal_error
        end select
        select case(trim(transition))
          case('lctm2015')
            call copy3(tgm, "flat", face)
          case DEFAULT
            continue
        end select
      end subroutine supersonic_outlet


      subroutine subsonic_inlet(face)
        !< Subsonic inlet boundary condition. 
        !< All the state variables's value expect pressure
        !< is fixed and pressure is copied from inside the 
        !< domain
        implicit none
        character(len=*), intent(in) :: face
        !< Name of the face at which boundary condition is called
        if(current_iter<=2)then
        call fix(density , fixed_density , face)
        call fix(x_speed , fixed_x_speed , face)
        call fix(y_speed , fixed_y_speed , face)
        call fix(z_speed , fixed_z_speed , face)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sa', 'saBC')
            call fix(tv, fixed_tv, face)
          case('sst', 'sst2003')
            call check_if_value_fixed("sst")
            call fix(tk, fixed_tk, face)
            call fix(tw, fixed_tw, face)
          case('kkl')
            call check_if_value_fixed("kkl")
            call fix(tk, fixed_tk, face)
            call fix(tkl, fixed_tw, face)
          case DEFAULT
           ! call turbulence_read_error()
           Fatal_error
        end select
        select case(trim(transition))
          case('lctm2015')
            call check_if_value_fixed("lctm2015")
            call fix(tgm, fixed_tgm, face)
          case DEFAULT
            continue
        end select
        end if
        call copy3(pressure, "flat", face)
      end subroutine subsonic_inlet


      subroutine subsonic_outlet(face)
        !< Subsonic outlet boundary condition. 
        !< All the state variables's value expect pressure
        !< is copied from the inside of the domain and pressure 
        !< is fixed
        implicit none
        character(len=*), intent(in) :: face
        !< Name of the face at which boundary condition is called
        call copy3(density, "flat", face)
        call copy3(x_speed, "flat", face)
        call copy3(y_speed, "flat", face)
        call copy3(z_speed, "flat", face)
        if(current_iter<=2)then
        call fix(pressure, fixed_pressure, face)
        end if
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sa', 'saBC')
            call copy3(tv, "flat", face)
          case('sst', 'sst2003')
            call copy3(tk, "flat", face)
            call copy3(tw, "flat", face)
          case('kkl')
            call copy3(tk, "flat", face)
            call copy3(tkl, "flat", face)
          case DEFAULT
           ! call turbulence_read_error()
           Fatal_error
        end select
        select case(trim(transition))
          case('lctm2015')
            call copy3(tgm, "flat", face)
          case DEFAULT
            continue
        end select
      end subroutine subsonic_outlet

      subroutine wall(face)
        !< Adiabatic/Isothermal wall boundary condition
        implicit none
        character(len=*), intent(in) :: face
        !< Name of the face at which boundary condition is called
        call copy3(pressure, "symm",  face)
        call temp_based_density(fixed_wall_temperature, face)
        call no_slip(face)
      end subroutine wall


      subroutine slip_wall(face)
        !< Slip wall boundary condition. 
        !< Maintain flow tangency
        implicit none
        character(len=*), intent(in) :: face
        !< Name of the face at which boundary condition is called
        call copy3(density , "symm", face)
        call copy3(pressure, "symm", face)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sa', 'saBC')
            call copy3(tv, "symm", face)
          case('sst', 'sst2003')
            call copy3(tk, "symm", face)
            call copy3(tw, "symm", face)
          case('kkl')
            call copy3(tk, "symm", face)
            call copy3(tkl, "symm", face)
          case DEFAULT
            !call turbulence_read_error()
            Fatal_error
        end select
        select case(trim(transition))
          case('lctm2015')
            call copy3(tgm, "flat", face)
          case DEFAULT
            continue
        end select
        call flow_tangency(face)
      end subroutine slip_wall


      subroutine pole(face)
        !< Boundary condition for the block face
        !< with zero area; turning into a pole
        implicit none
        character(len=*), intent(in) :: face
        !< Name of the face at which boundary condition is called
        call copy3(density , "flat", face)
        call copy3(x_speed , "flat", face)
        call copy3(y_speed , "flat", face)
        call copy3(z_speed , "flat", face)
        call copy3(pressure, "flat", face)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sa', 'saBC')
            call copy3(tv, "flat", face) 
          case('sst', 'sst2003')
            call copy3(tk, "flat", face)
            call copy3(tw, "flat", face)
          case('kkl')
            call copy3(tk, "flat", face)
            call copy3(tkl, "flat", face)
          case DEFAULT
            call turbulence_read_error()
        end select
        select case(trim(transition))
          case('lctm2015')
            call copy3(tgm, "flat", face)
          case DEFAULT
            continue
        end select
      end subroutine pole



      subroutine fix(var, fix_val, face)
        !< Generalized subroutine to fix particular value
        !< at particular face
        implicit none
        real, dimension(-2:imx+2, -2:jmx+2, -2:kmx+2) , intent(out) :: var
        !< Variable of which values are being fixed in the ghost cell
        real, dimension(1:6)       , intent(in)  :: fix_val
        !< Amount of value that need to be fixed.
        character(len=*)         , intent(in)  :: face
        !< Name of the face at which boundary condition is called

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
            Fatal_error
        end select
            
      end subroutine fix


      subroutine no_slip(face)
        !< No-slip wall boundary condition. All the 
        !< component of velocity throught face is zero
        implicit none
        character(len=*), intent(in) :: face
        !< Name of the face at which boundary condition is called
        call copy3(x_speed, "anti", face)
        call copy3(y_speed, "anti", face)
        call copy3(z_speed, "anti", face)
        select case(turbulence)
          case("none")
            !do nothing
            continue
          case('sa', 'saBC')
            call copy3(tv  , "anti", face)
          case("sst", 'sst2003')
            call copy3(tk  , "anti", face)
            call set_omega_at_wall(face)
          case("kkl")
            call copy3(tk  , "anti", face)
            call copy3(tkl , "anti", face)
          case DEFAULT
            !call turbulence_read_error()
            Fatal_error
        end select
        select case(trim(transition))
          case('lctm2015')
            call copy3(tgm, "flat", face)
          case DEFAULT
            continue
        end select
      end subroutine no_slip

    

      subroutine set_omega_at_wall(face)
        !< Set value of turbulence variable: omega (turbulenct dissipation rate). 
        !< Value fixed is accourding to the SST turbulence model
        implicit none
        character(len=*), intent(in) :: face
        real :: T_face
        real :: mu
        real :: rho
        integer :: i,j,k,l
        
        select case(face)
          case("imin")
            do l=1,3
          do k = 1,kmx-1
            do j = 1,jmx-1
              T_face = 0.5*((pressure(0, j, k)/density(0, j, k))+(pressure(1, j, k)/density(1, j, k)))/R_gas
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(0, j, k) + density(1, j, k))
              tw(1-l, j, k) = 120*mu/(rho*beta1*(2*dist(1, j, k))**2) - tw(l, j, k)
            end do
          end do
        end do
        case("imax")
          do l=1,3
          do k = 1,kmx-1
            do j = 1,jmx-1
              T_face = 0.5*((pressure(imx-1, j, k)/density(imx-1, j, k))+(pressure(imx, j, k)/density(imx, j, k)))/R_gas
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(imx-1, j, k) + density(imx, j, k))
              tw(imx+l-1, j, k) = 120*mu/(rho*beta1*(2*dist(imx-1, j, k))**2) - tw(imx-l, j, k)
            end do
          end do
        end do
        case("jmin")
          do l=1,3
          do k = 1,kmx-1
            do i = 1,imx-1
              T_face = 0.5*((pressure(i, 0, k)/density(i, 0, k))+(pressure(i, 1, k)/density(i, 1, k)))/R_gas
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(i, 0, k) + density(i, 1, k))
              tw(i, 1-l, k) = 120*mu/(rho*beta1*(2*dist(i, 1, k))**2) - tw(i, l, k)
            end do
          end do
        end do
        case("jmax")
          do l=1,3
          do k = 1,kmx-1
            do i = 1,imx-1
              T_face = 0.5*((pressure(i, jmx-1, k)/density(i, jmx-1, k))+(pressure(i, jmx, k)/density(i, jmx, k)))/R_gas
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(i, jmx-1, k) + density(i, jmx, k))
              tw(i, jmx+l-1, k) = 120*mu/(rho*beta1*(2*dist(i, jmx-1, k))**2) - tw(i, jmx-l, k)
            end do
          end do
        end do
        case("kmin")
          do l=1,3
          do j = 1,jmx-1
            do i = 1,imx-1
              T_face = 0.5*((pressure(i, j, 0)/density(i, j, 0))+(pressure(i, j, 1)/density(i, j, 1)))/R_gas
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(i, j, 0) + density(i, j, 1))
              tw(i, j, 1-l) = 120*mu/(rho*beta1*(2*dist(i, j, 1))**2) - tw(i, j, l)
            end do
          end do
        end do
        case("kmax")
          do l=1,3
          do j = 1,jmx-1
            do i = 1,imx-1
              T_face = 0.5*((pressure(i, j, kmx-1)/density(i, j, kmx-1))+(pressure(i, j, kmx)/density(i, j, kmx)))/R_gas
              mu = mu_ref * (T_face/T_ref)**1.5*((T_ref + Sutherland_temp )/(T_face + Sutherland_temp))
              rho = 0.5 * (density(i, j, kmx-1) + density(i, j, kmx))
              tw(i, j, kmx+l-1) = 120*mu/(rho*beta1*(2*dist(i, j, kmx-1))**2) - tw(i, j, kmx-l)
            end do
          end do
        end do

      end select
    end subroutine set_omega_at_wall

    subroutine check_if_value_fixed(model)
      !< A Fail-check subroutine which set the freestream
      !< as the fixed value in case not specified explicitly
      implicit none
      character(len=*), intent(in) :: model

      select case(model)
        case("none", "lctm2015")
          !do nothing
          continue
        case("sa", 'saBC')
          if(fixed_tv(face_num)==0.0)fixed_tv(face_num)=tv_inf
        case("sst", 'sst2003')
          if(fixed_tk(face_num)==0.) fixed_tk(face_num)=tk_inf
          if(fixed_tw(face_num)==0.) fixed_tw(face_num)=tw_inf
        case("kkl")
          if(fixed_tk(face_num)==0.) fixed_tk(face_num)=tk_inf
          if(fixed_tkl(face_num)==0.) fixed_tkl(face_num)=tkl_inf
        case DEFAULT
         ! call turbulence_read_error()
         Fatal_error
      end select

      select case(trim(transition))
        case('lctm2015')
          if(fixed_tgm(face_num)==0.0) fixed_tgm(face_num)=tgm_inf
        Case DEFAULT
          continue
      end select
    end subroutine check_if_value_fixed

    subroutine far_field(face)
      !< Far-field Riemann boundary condition
      implicit none
      character(len=*) :: face
      real :: cinf, cexp   ! speed of sound
      real :: Rinf, Rexp   ! Riemann invarient
      real :: Uninf, Unexp ! face normal speed
      real :: Unb ! normal velocity boundary
      real :: Cb  ! speed of sound boundary
      real :: vel_diff
      real :: u,v,w
      real :: uf, vf, wf
      integer :: i,j,k
      real :: s
      integer, dimension(6) :: face_already_has_fixed_values=0!0=.no.

      face_already_has_fixed_values=0
      select case(face)
        case("imin")
          do k = 1,kmx-1
            do j = 1,jmx-1
              do i = 1,1
                ! interior cell
                u = x_speed(i,j,k)
                v = y_speed(i,j,k)
                w = z_speed(i,j,k)
                ! ghost cell
                uf = x_speed_inf!x_speed(i-1,j,k)
                vf = y_speed_inf!y_speed(i-1,j,k)
                wf = z_speed_inf!z_speed(i-1,j,k)
                cexp = sqrt(gm*pressure(i,j,k)/density(i,j,k))
                !cinf = sqrt(gm*pressure(i-1,j,k)/density(i-1,j,k))
                cinf = sqrt(gm*pressure_inf/density_inf)
                Unexp = u *(-xnx(i,j,k)) + v *(-xny(i,j,k)) + w *(-xnz(i,j,k))
                Uninf = uf*(-xnx(i,j,k)) + vf*(-xny(i,j,k)) + wf*(-xnz(i,j,k))
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i-1,j,k) = x_speed(i,j,k) + vel_diff*(-xnx(i,j,k))
                  y_speed(i-1,j,k) = y_speed(i,j,k) + vel_diff*(-xny(i,j,k))
                  z_speed(i-1,j,k) = z_speed(i,j,k) + vel_diff*(-xnz(i,j,k))
                  s = pressure(i,j,k)/(density(i,j,k)**(gm))
                  density(i-1,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i-1,j,k) = (density(i-1,j,k)*Cb*Cb/gm)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face)
                      call copy3(tw, "flat", face)
                    case('kkl')
                      call copy3(tk, "flat", face)
                      call copy3(tkl, "flat", face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face)
                    case DEFAULT
                      continue
                  end select
                  face_already_has_fixed_values(1)=0
                else
                  vel_diff = Unb - Uninf
                  x_speed(i-1,j,k) = x_speed_inf + vel_diff*(-xnx(i,j,k))
                  y_speed(i-1,j,k) = y_speed_inf + vel_diff*(-xny(i,j,k))
                  z_speed(i-1,j,k) = z_speed_inf + vel_diff*(-xnz(i,j,k))
                  s = pressure_inf/(density_inf**(gm))
                  density(i-1,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i-1,j,k) = (density(i-1,j,k)*Cb*Cb/gm)

                  if(face_already_has_fixed_values(1)==0)then
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, fixed_tv, face)
                    case('sst', 'sst2003')
                      call check_if_value_fixed("sst")
                      call fix(tk, fixed_tk, face)
                      call fix(tw, fixed_tw, face)
                    case('kkl')
                      call check_if_value_fixed("kkl")
                      call fix(tk, fixed_tk, face)
                      call fix(tkl, fixed_tkl, face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call check_if_value_fixed("lctm2015")
                      call fix(tgm, fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                  end if
                  face_already_has_fixed_values(1)=1
                end if
              end do
            end do
          end do
          qp(-1,:,:,:) = qp(0,:,:,:)
          qp(-2,:,:,:) = qp(0,:,:,:)
         case("imax")
          do k = 1,kmx-1
            do j = 1,jmx-1
              do i = imx,imx
                ! interior cell
                u = x_speed(i-1,j,k)
                v = y_speed(i-1,j,k)
                w = z_speed(i-1,j,k)
                ! ghost cell
                uf = x_speed_inf!x_speed(i,j,k)
                vf = y_speed_inf!y_speed(i,j,k)
                wf = z_speed_inf!z_speed(i,j,k)
                cexp = sqrt(gm*pressure(i-1,j,k)/density(i-1,j,k))
                !cinf = sqrt(gm*pressure(i,j,k)/density(i,j,k))
                cinf = sqrt(gm*pressure_inf/density_inf)
                Unexp = u *(xnx(i,j,k)) + v *(xny(i,j,k)) + w *(xnz(i,j,k))
                Uninf = uf*(xnx(i,j,k)) + vf*(xny(i,j,k)) + wf*(xnz(i,j,k))
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k) = x_speed(i-1,j,k) + vel_diff*(xnx(i,j,k))
                  y_speed(i,j,k) = y_speed(i-1,j,k) + vel_diff*(xny(i,j,k))
                  z_speed(i,j,k) = z_speed(i-1,j,k) + vel_diff*(xnz(i,j,k))
                  s = pressure(i-1,j,k)/(density(i-1,j,k)**(gm))
                  density(i,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k) = (density(i,j,k)*Cb*Cb/gm)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face)
                      call copy3(tw, "flat", face)
                    case('kkl')
                      call copy3(tk, "flat", face)
                      call copy3(tkl, "flat", face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face)
                    case DEFAULT
                      continue
                  end select
                  face_already_has_fixed_values(2)=0
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k) = x_speed_inf + vel_diff*(xnx(i,j,k))
                  y_speed(i,j,k) = y_speed_inf + vel_diff*(xny(i,j,k))
                  z_speed(i,j,k) = z_speed_inf + vel_diff*(xnz(i,j,k))
                  s = pressure_inf/(density_inf**(gm))
                  density(i,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k) = (density(i,j,k)*Cb*Cb/gm)

                  if(face_already_has_fixed_values(2)==0)then
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, fixed_tv, face)
                    case('sst', 'sst2003')
                      call check_if_value_fixed("sst")
                      call fix(tk, fixed_tk, face)
                      call fix(tw, fixed_tw, face)
                    case('kkl')
                      call check_if_value_fixed("kkl")
                      call fix(tk, fixed_tk, face)
                      call fix(tkl, fixed_tkl, face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call check_if_value_fixed("lctm2015")
                      call fix(tgm, fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                  end if
                  face_already_has_fixed_values(2)=1
                end if
              end do
            end do
          end do
          qp(imx+1,:,:,:) = qp(imx,:,:,:)
          qp(imx+2,:,:,:) = qp(imx,:,:,:)
        case("jmin")
          do k = 1,kmx-1
            do j = 1,1
              do i = 1,imx-1
                ! interior cell
                u = x_speed(i,j,k)
                v = y_speed(i,j,k)
                w = z_speed(i,j,k)
                ! ghost cell
                uf = x_speed_inf!x_speed(i,j-1,k)
                vf = y_speed_inf!y_speed(i,j-1,k)
                wf = z_speed_inf!z_speed(i,j-1,k)
                cexp = sqrt(gm*pressure(i,j,k)/density(i,j,k))
                !cinf = sqrt(gm*pressure(i,j-1,k)/density(i,j-1,k))
                cinf = sqrt(gm*pressure_inf/density_inf)
                Unexp = u *(-ynx(i,j,k)) + v *(-yny(i,j,k)) + w *(-ynz(i,j,k))
                Uninf = uf*(-ynx(i,j,k)) + vf*(-yny(i,j,k)) + wf*(-ynz(i,j,k))
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j-1,k) = x_speed(i,j,k) + vel_diff*(-ynx(i,j,k))
                  y_speed(i,j-1,k) = y_speed(i,j,k) + vel_diff*(-yny(i,j,k))
                  z_speed(i,j-1,k) = z_speed(i,j,k) + vel_diff*(-ynz(i,j,k))
                  s = pressure(i,j,k)/(density(i,j,k)**(gm))
                  density(i,j-1,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j-1,k) = (density(i,j-1,k)*Cb*Cb/gm)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face)
                      call copy3(tw, "flat", face)
                    case('kkl')
                      call copy3(tk, "flat", face)
                      call copy3(tkl, "flat", face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face)
                    case DEFAULT
                      continue
                  end select
                  face_already_has_fixed_values(3)=0
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j-1,k) = x_speed_inf + vel_diff*(-ynx(i,j,k))
                  y_speed(i,j-1,k) = y_speed_inf + vel_diff*(-yny(i,j,k))
                  z_speed(i,j-1,k) = z_speed_inf + vel_diff*(-ynz(i,j,k))
                  s = pressure_inf/(density_inf**(gm))
                  density(i,j-1,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j-1,k) = (density(i,j-1,k)*Cb*Cb/gm)

                  if(face_already_has_fixed_values(3)==0)then
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, fixed_tv, face)
                    case('sst', 'sst2003')
                      call check_if_value_fixed("sst")
                      call fix(tk, fixed_tk, face)
                      call fix(tw, fixed_tw, face)
                    case('kkl')
                      call check_if_value_fixed("kkl")
                      call fix(tk, fixed_tk, face)
                      call fix(tkl, fixed_tkl, face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call check_if_value_fixed("lctm2015")
                      call fix(tgm, fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                  end if
                  face_already_has_fixed_values(3)=1
                end if
              end do
            end do
          end do
          qp(:,-1,:,:) = qp(:,0,:,:)
          qp(:,-2,:,:) = qp(:,0,:,:)
        case("jmax")
          do k = 1,kmx-1
            do j = jmx,jmx
              do i = 1,imx-1
                ! interior cell
                u = x_speed(i,j-1,k)
                v = y_speed(i,j-1,k)
                w = z_speed(i,j-1,k)
                ! ghost cell
                uf = x_speed_inf!x_speed(i,j,k)
                vf = y_speed_inf!y_speed(i,j,k)
                wf = z_speed_inf!z_speed(i,j,k)
                cexp = sqrt(gm*pressure(i,j-1,k)/density(i,j-1,k))
                !cinf = sqrt(gm*pressure(i,j,k)/density(i,j,k))
                cinf = sqrt(gm*pressure_inf/density_inf)
                Unexp = u *(ynx(i,j,k)) + v *(yny(i,j,k)) + w *(ynz(i,j,k))
                Uninf = uf*(ynx(i,j,k)) + vf*(yny(i,j,k)) + wf*(ynz(i,j,k))
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k) = x_speed(i,j-1,k) + vel_diff*(ynx(i,j,k))
                  y_speed(i,j,k) = y_speed(i,j-1,k) + vel_diff*(yny(i,j,k))
                  z_speed(i,j,k) = z_speed(i,j-1,k) + vel_diff*(ynz(i,j,k))
                  s = pressure(i,j-1,k)/(density(i,j-1,k)**(gm))
                  density(i,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k) = (density(i,j,k)*Cb*Cb/gm)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face)
                      call copy3(tw, "flat", face)
                    case('kkl')
                      call copy3(tk, "flat", face)
                      call copy3(tkl, "flat", face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face)
                    case DEFAULT
                      continue
                  end select
                  face_already_has_fixed_values(4)=0
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k) = x_speed_inf + vel_diff*(ynx(i,j,k))
                  y_speed(i,j,k) = y_speed_inf + vel_diff*(yny(i,j,k))
                  z_speed(i,j,k) = z_speed_inf + vel_diff*(ynz(i,j,k))
                  s = pressure_inf/(density_inf**(gm))
                  density(i,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k) = (density(i,j,k)*Cb*Cb/gm)

                  if(face_already_has_fixed_values(4)==0)then
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, fixed_tv, face)
                    case('sst', 'sst2003')
                      call check_if_value_fixed("sst")
                      call fix(tk, fixed_tk, face)
                      call fix(tw, fixed_tw, face)
                    case('kkl')
                      call check_if_value_fixed("kkl")
                      call fix(tk, fixed_tk, face)
                      call fix(tkl, fixed_tkl, face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call check_if_value_fixed("lctm2015")
                      call fix(tgm, fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                  end if
                  face_already_has_fixed_values(4)=1
                end if
              end do
            end do
          end do
          qp(:,jmx+1,:,:) = qp(:,jmx,:,:)
          qp(:,jmx+2,:,:) = qp(:,jmx,:,:)
        case("kmin")
          do k = 1,1
            do j = 1,jmx-1
              do i = 1,imx-1
                ! interior cell
                u = x_speed(i,j,k)
                v = y_speed(i,j,k)
                w = z_speed(i,j,k)
                ! ghost cell
                uf = x_speed_inf!x_speed(i,j,k-1)
                vf = y_speed_inf!y_speed(i,j,k-1)
                wf = z_speed_inf!z_speed(i,j,k-1)
                cexp = sqrt(gm*pressure(i,j,k)/density(i,j,k))
                !cinf = sqrt(gm*pressure(i,j,k-1)/density(i,j,k-1))
                cinf = sqrt(gm*pressure_inf/density_inf)
                Unexp = u *(-znx(i,j,k)) + v *(-zny(i,j,k)) + w *(-znz(i,j,k))
                Uninf = uf*(-znx(i,j,k)) + vf*(-zny(i,j,k)) + wf*(-znz(i,j,k))
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k-1) = x_speed(i,j,k) + vel_diff*(-znx(i,j,k))
                  y_speed(i,j,k-1) = y_speed(i,j,k) + vel_diff*(-zny(i,j,k))
                  z_speed(i,j,k-1) = z_speed(i,j,k) + vel_diff*(-znz(i,j,k))
                  s = pressure(i,j,k)/(density(i,j,k)**(gm))
                  density(i,j,k-1) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k-1) = (density(i,j,k-1)*Cb*Cb/gm)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face)
                      call copy3(tw, "flat", face)
                    case('kkl')
                      call copy3(tk, "flat", face)
                      call copy3(tkl, "flat", face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face)
                    case DEFAULT
                      continue
                  end select
                  face_already_has_fixed_values(5)=0
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k-1) = x_speed_inf + vel_diff*(-znx(i,j,k))
                  y_speed(i,j,k-1) = y_speed_inf + vel_diff*(-zny(i,j,k))
                  z_speed(i,j,k-1) = z_speed_inf + vel_diff*(-znz(i,j,k))
                  s = pressure_inf/(density_inf**(gm))
                  density(i,j,k-1) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k-1) = (density(i,j,k-1)*Cb*Cb/gm)

                  if(face_already_has_fixed_values(5)==0)then
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, fixed_tv, face)
                    case('sst', 'sst2003')
                      call check_if_value_fixed("sst")
                      call fix(tk, fixed_tk, face)
                      call fix(tw, fixed_tw, face)
                    case('kkl')
                      call check_if_value_fixed("kkl")
                      call fix(tk, fixed_tk, face)
                      call fix(tkl, fixed_tkl, face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call check_if_value_fixed("lctm2015")
                      call fix(tgm, fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                  end if
                  face_already_has_fixed_values(5)=1
                end if
              end do
            end do
          end do
          qp(:,:,-1,:) = qp(:,:,0,:)
          qp(:,:,-2,:) = qp(:,:,0,:)
        case("kmax")
          do k = kmx,kmx
            do j = 1,jmx-1
              do i = 1,imx-1
                ! interior cell
                u = x_speed(i,j,k-1)
                v = y_speed(i,j,k-1)
                w = z_speed(i,j,k-1)
                ! ghost cell
                uf = x_speed_inf!x_speed(i,j,k)
                vf = y_speed_inf!y_speed(i,j,k)
                wf = z_speed_inf!z_speed(i,j,k)
                cexp = sqrt(gm*pressure(i,j,k-1)/density(i,j,k-1))
                !cinf = sqrt(gm*pressure(i,j,k)/density(i,j,k))
                cinf = sqrt(gm*pressure_inf/density_inf)
                Unexp = u *(znx(i,j,k)) + v *(zny(i,j,k)) + w *(znz(i,j,k))
                Uninf = uf*(znx(i,j,k)) + vf*(zny(i,j,k)) + wf*(znz(i,j,k))
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k) = x_speed(i,j,k-1) + vel_diff*(znx(i,j,k))
                  y_speed(i,j,k) = y_speed(i,j,k-1) + vel_diff*(zny(i,j,k))
                  z_speed(i,j,k) = z_speed(i,j,k-1) + vel_diff*(znz(i,j,k))
                  s = pressure(i,j,k-1)/(density(i,j,k-1)**(gm))
                  density(i,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k) = (density(i,j,k)*Cb*Cb/gm)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face)
                      call copy3(tw, "flat", face)
                    case('kkl')
                      call copy3(tk, "flat", face)
                      call copy3(tkl, "flat", face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face)
                    case DEFAULT
                      continue
                  end select
                  face_already_has_fixed_values(6)=0
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k) = x_speed_inf + vel_diff*(znx(i,j,k))
                  y_speed(i,j,k) = y_speed_inf + vel_diff*(zny(i,j,k))
                  z_speed(i,j,k) = z_speed_inf + vel_diff*(znz(i,j,k))
                  s = pressure_inf/(density_inf**(gm))
                  density(i,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k) = (density(i,j,k)*Cb*Cb/gm)

                  if(face_already_has_fixed_values(6)==0)then
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, fixed_tv, face)
                    case('sst', 'sst2003')
                      call check_if_value_fixed("sst")
                      call fix(tk, fixed_tk, face)
                      call fix(tw, fixed_tw, face)
                    case('kkl')
                      call check_if_value_fixed("kkl")
                      call fix(tk, fixed_tk, face)
                      call fix(tkl, fixed_tkl, face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call check_if_value_fixed("lctm2015")
                      call fix(tgm, fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                  end if
                  face_already_has_fixed_values(6)=1
                end if
              end do
            end do
          end do
          qp(:,:,kmx+1,:) = qp(:,:,kmx,:)
          qp(:,:,kmx+2,:) = qp(:,:,kmx,:)
        case DEFAULT
          !print*, "ERROR: wrong face for boundary condition"
          Fatal_error
      end select

    end subroutine far_field


    subroutine total_pressure(face)
      !< Total Pressure Riemann boundary condition
      implicit none
      character(len=*) :: face
      real :: cinf, cexp   ! speed of sound
      real :: Rinf, Rexp   ! Riemann invarient
      real :: Uninf, Unexp ! face normal speed
      real :: Unb ! normal velocity boundary
      real :: Cb  ! speed of sound boundary
      real :: vel_diff
      real :: u,v,w
      real :: uf, vf, wf
      real :: Mb
      integer :: i,j,k

      select case(face)
        case("imin")
          do k = 1,kmx-1
            do j = 1,jmx-1
              do i = 1,1
                ! interior cell
                u = x_speed(i,j,k)
                v = y_speed(i,j,k)
                w = z_speed(i,j,k)
                ! ghost cell
                uf = x_speed_inf!x_speed(i-1,j,k)
                vf = y_speed_inf!y_speed(i-1,j,k)
                wf = z_speed_inf!z_speed(i-1,j,k)
                cexp = sqrt(gm*pressure(i,j,k)/density(i,j,k))
                !cinf = sqrt(gm*pressure(i-1,j,k)/density(i-1,j,k))
                cinf = sqrt(gm*pressure_inf/density_inf)
                Unexp = u *(-xnx(i,j,k)) + v *(-xny(i,j,k)) + w *(-xnz(i,j,k))
                Uninf = uf*(-xnx(i,j,k)) + vf*(-xny(i,j,k)) + wf*(-xnz(i,j,k))
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i-1,j,k) = x_speed(i,j,k) + vel_diff*(-xnx(i,j,k))
                  y_speed(i-1,j,k) = y_speed(i,j,k) + vel_diff*(-xny(i,j,k))
                  z_speed(i-1,j,k) = z_speed(i,j,k) + vel_diff*(-xnz(i,j,k))
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face)
                      call copy3(tw, "flat", face)
                    case('kkl')
                      call copy3(tk, "flat", face)
                      call copy3(tkl, "flat", face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face)
                    case DEFAULT
                      continue
                  end select
                else
                  vel_diff = Unb - Uninf
                  x_speed(i-1,j,k) = x_speed_inf + vel_diff*(-xnx(i,j,k))
                  y_speed(i-1,j,k) = y_speed_inf + vel_diff*(-xny(i,j,k))
                  z_speed(i-1,j,k) = z_speed_inf + vel_diff*(-xnz(i,j,k))
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, fixed_tv, face)
                    case('sst', 'sst2003')
                      call check_if_value_fixed("sst")
                      call fix(tk, fixed_tk, face)
                      call fix(tw, fixed_tw, face)
                    case('kkl')
                      call check_if_value_fixed("kkl")
                      call fix(tk, fixed_tk, face)
                      call fix(tkl, fixed_tkl, face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call check_if_value_fixed("lctm2015")
                      call fix(tgm, fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                end if
                Mb = sqrt(x_speed(i-1,j,k)**2+y_speed(i-1,j,k)**2+z_speed(i-1,j,k)**2)/Cb
                pressure(i-1,j,k) = fixed_Tpressure(1)/(((1+0.5*(gm-1.)*Mb*Mb))**(gm/(gm-1.)))
                density(i-1,j,k) = gm*pressure(i-1,j,k)/(Cb*Cb)
              end do
            end do
          end do
          qp(-1,:,:,:) = qp(0,:,:,:)
          qp(-2,:,:,:) = qp(0,:,:,:)
         case("imax")
          do k = 1,kmx-1
            do j = 1,jmx-1
              do i = imx,imx
                ! interior cell
                u = x_speed(i-1,j,k)
                v = y_speed(i-1,j,k)
                w = z_speed(i-1,j,k)
                ! ghost cell
                uf = x_speed_inf!x_speed(i,j,k)
                vf = y_speed_inf!y_speed(i,j,k)
                wf = z_speed_inf!z_speed(i,j,k)
                cexp = sqrt(gm*pressure(i-1,j,k)/density(i-1,j,k))
                !cinf = sqrt(gm*pressure(i,j,k)/density(i,j,k))
                cinf = sqrt(gm*pressure_inf/density_inf)
                Unexp = u *(xnx(i,j,k)) + v *(xny(i,j,k)) + w *(xnz(i,j,k))
                Uninf = uf*(xnx(i,j,k)) + vf*(xny(i,j,k)) + wf*(xnz(i,j,k))
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k) = x_speed(i-1,j,k) + vel_diff*(xnx(i,j,k))
                  y_speed(i,j,k) = y_speed(i-1,j,k) + vel_diff*(xny(i,j,k))
                  z_speed(i,j,k) = z_speed(i-1,j,k) + vel_diff*(xnz(i,j,k))
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face)
                      call copy3(tw, "flat", face)
                    case('kkl')
                      call copy3(tk, "flat", face)
                      call copy3(tkl, "flat", face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face)
                    case DEFAULT
                      continue
                  end select
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k) = x_speed_inf + vel_diff*(xnx(i,j,k))
                  y_speed(i,j,k) = y_speed_inf + vel_diff*(xny(i,j,k))
                  z_speed(i,j,k) = z_speed_inf + vel_diff*(xnz(i,j,k))
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, fixed_tv, face)
                    case('sst', 'sst2003')
                      call check_if_value_fixed("sst")
                      call fix(tk, fixed_tk, face)
                      call fix(tw, fixed_tw, face)
                    case('kkl')
                      call check_if_value_fixed("kkl")
                      call fix(tk, fixed_tk, face)
                      call fix(tkl, fixed_tkl, face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call check_if_value_fixed("lctm2015")
                      call fix(tgm, fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                end if
                Mb = sqrt(x_speed(i,j,k)**2+y_speed(i,j,k)**2+z_speed(i,j,k)**2)/Cb
                pressure(i,j,k) = fixed_Tpressure(2)/(((1+0.5*(gm-1.)*Mb*Mb))**(gm/(gm-1.)))
                density(i,j,k) = gm*pressure(i,j,k)/(Cb*Cb)
              end do
            end do
          end do
          qp(imx+1,:,:,:) = qp(imx,:,:,:)
          qp(imx+2,:,:,:) = qp(imx,:,:,:)
        case("jmin")
          do k = 1,kmx-1
            do j = 1,1
              do i = 1,imx-1
                ! interior cell
                u = x_speed(i,j,k)
                v = y_speed(i,j,k)
                w = z_speed(i,j,k)
                ! ghost cell
                uf = x_speed_inf!x_speed(i,j-1,k)
                vf = y_speed_inf!y_speed(i,j-1,k)
                wf = z_speed_inf!z_speed(i,j-1,k)
                cexp = sqrt(gm*pressure(i,j,k)/density(i,j,k))
                !cinf = sqrt(gm*pressure(i,j-1,k)/density(i,j-1,k))
                cinf = sqrt(gm*pressure_inf/density_inf)
                Unexp = u *(-ynx(i,j,k)) + v *(-yny(i,j,k)) + w *(-ynz(i,j,k))
                Uninf = uf*(-ynx(i,j,k)) + vf*(-yny(i,j,k)) + wf*(-ynz(i,j,k))
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j-1,k) = x_speed(i,j,k) + vel_diff*(-ynx(i,j,k))
                  y_speed(i,j-1,k) = y_speed(i,j,k) + vel_diff*(-yny(i,j,k))
                  z_speed(i,j-1,k) = z_speed(i,j,k) + vel_diff*(-ynz(i,j,k))
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face)
                      call copy3(tw, "flat", face)
                    case('kkl')
                      call copy3(tk, "flat", face)
                      call copy3(tkl, "flat", face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face)
                    case DEFAULT
                      continue
                  end select
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j-1,k) = x_speed_inf + vel_diff*(-ynx(i,j,k))
                  y_speed(i,j-1,k) = y_speed_inf + vel_diff*(-yny(i,j,k))
                  z_speed(i,j-1,k) = z_speed_inf + vel_diff*(-ynz(i,j,k))
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, fixed_tv, face)
                    case('sst', 'sst2003')
                      call check_if_value_fixed("sst")
                      call fix(tk, fixed_tk, face)
                      call fix(tw, fixed_tw, face)
                    case('kkl')
                      call check_if_value_fixed("kkl")
                      call fix(tk, fixed_tk, face)
                      call fix(tkl, fixed_tkl, face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call check_if_value_fixed("lctm2015")
                      call fix(tgm, fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                end if
                Mb = sqrt(x_speed(i,j-1,k)**2+y_speed(i,j-1,k)**2+z_speed(i,j-1,k)**2)/Cb
                pressure(i,j-1,k) = fixed_Tpressure(3)/(((1+0.5*(gm-1.)*Mb*Mb))**(gm/(gm-1.)))
                density(i,j-1,k) = gm*pressure(i,j-1,k)/(Cb*Cb)
              end do
            end do
          end do
          qp(:,-1,:,:) = qp(:,0,:,:)
          qp(:,-2,:,:) = qp(:,0,:,:)
        case("jmax")
          do k = 1,kmx-1
            do j = jmx,jmx
              do i = 1,imx-1
                ! interior cell
                u = x_speed(i,j-1,k)
                v = y_speed(i,j-1,k)
                w = z_speed(i,j-1,k)
                ! ghost cell
                uf = x_speed_inf!x_speed(i,j,k)
                vf = y_speed_inf!y_speed(i,j,k)
                wf = z_speed_inf!z_speed(i,j,k)
                cexp = sqrt(gm*pressure(i,j-1,k)/density(i,j-1,k))
                !cinf = sqrt(gm*pressure(i,j,k)/density(i,j,k))
                cinf = sqrt(gm*pressure_inf/density_inf)
                Unexp = u *(ynx(i,j,k)) + v *(yny(i,j,k)) + w *(ynz(i,j,k))
                Uninf = uf*(ynx(i,j,k)) + vf*(yny(i,j,k)) + wf*(ynz(i,j,k))
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k) = x_speed(i,j-1,k) + vel_diff*(ynx(i,j,k))
                  y_speed(i,j,k) = y_speed(i,j-1,k) + vel_diff*(yny(i,j,k))
                  z_speed(i,j,k) = z_speed(i,j-1,k) + vel_diff*(ynz(i,j,k))
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face)
                      call copy3(tw, "flat", face)
                    case('kkl')
                      call copy3(tk, "flat", face)
                      call copy3(tkl, "flat", face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face)
                    case DEFAULT
                      continue
                  end select
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k) = x_speed_inf + vel_diff*(ynx(i,j,k))
                  y_speed(i,j,k) = y_speed_inf + vel_diff*(yny(i,j,k))
                  z_speed(i,j,k) = z_speed_inf + vel_diff*(ynz(i,j,k))
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, fixed_tv, face)
                    case('sst', 'sst2003')
                      call check_if_value_fixed("sst")
                      call fix(tk, fixed_tk, face)
                      call fix(tw, fixed_tw, face)
                    case('kkl')
                      call check_if_value_fixed("kkl")
                      call fix(tk, fixed_tk, face)
                      call fix(tkl, fixed_tkl, face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call check_if_value_fixed("lctm2015")
                      call fix(tgm, fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                end if
                Mb = sqrt(x_speed(i,j,k)**2+y_speed(i,j,k)**2+z_speed(i,j,k)**2)/Cb
                pressure(i,j,k) = fixed_Tpressure(4)/(((1+0.5*(gm-1.)*Mb*Mb))**(gm/(gm-1.)))
                density(i,j,k) = gm*pressure(i,j,k)/(Cb*Cb)
              end do
            end do
          end do
          qp(:,jmx+1,:,:) = qp(:,jmx,:,:)
          qp(:,jmx+2,:,:) = qp(:,jmx,:,:)
        case("kmin")
          do k = 1,1
            do j = 1,jmx-1
              do i = 1,imx-1
                ! interior cell
                u = x_speed(i,j,k)
                v = y_speed(i,j,k)
                w = z_speed(i,j,k)
                ! ghost cell
                uf = x_speed_inf!x_speed(i,j,k-1)
                vf = y_speed_inf!y_speed(i,j,k-1)
                wf = z_speed_inf!z_speed(i,j,k-1)
                cexp = sqrt(gm*pressure(i,j,k)/density(i,j,k))
                !cinf = sqrt(gm*pressure(i,j,k-1)/density(i,j,k-1))
                cinf = sqrt(gm*pressure_inf/density_inf)
                Unexp = u *(-znx(i,j,k)) + v *(-zny(i,j,k)) + w *(-znz(i,j,k))
                Uninf = uf*(-znx(i,j,k)) + vf*(-zny(i,j,k)) + wf*(-znz(i,j,k))
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k-1) = x_speed(i,j,k) + vel_diff*(-znx(i,j,k))
                  y_speed(i,j,k-1) = y_speed(i,j,k) + vel_diff*(-zny(i,j,k))
                  z_speed(i,j,k-1) = z_speed(i,j,k) + vel_diff*(-znz(i,j,k))
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face)
                      call copy3(tw, "flat", face)
                    case('kkl')
                      call copy3(tk, "flat", face)
                      call copy3(tkl, "flat", face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face)
                    case DEFAULT
                      continue
                  end select
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k-1) = x_speed_inf + vel_diff*(-znx(i,j,k))
                  y_speed(i,j,k-1) = y_speed_inf + vel_diff*(-zny(i,j,k))
                  z_speed(i,j,k-1) = z_speed_inf + vel_diff*(-znz(i,j,k))
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, fixed_tv, face)
                    case('sst', 'sst2003')
                      call check_if_value_fixed("sst")
                      call fix(tk, fixed_tk, face)
                      call fix(tw, fixed_tw, face)
                    case('kkl')
                      call check_if_value_fixed("kkl")
                      call fix(tk, fixed_tk, face)
                      call fix(tkl, fixed_tkl, face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call check_if_value_fixed("lctm2015")
                      call fix(tgm, fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                end if
                Mb = sqrt(x_speed(i,j,k)**2+y_speed(i,j,k)**2+z_speed(i,j,k)**2)/Cb
                pressure(i,j,k-1) = fixed_Tpressure(5)/(((1+0.5*(gm-1.)*Mb*Mb))**(gm/(gm-1.)))
                density(i,j,k-1) = gm*pressure(i,j,k-1)/(Cb*Cb)
              end do
            end do
          end do
          qp(:,:,-1,:) = qp(:,:,0,:)
          qp(:,:,-2,:) = qp(:,:,0,:)
        case("kmax")
          do k = kmx,kmx
            do j = 1,jmx-1
              do i = 1,imx-1
                ! interior cell
                u = x_speed(i,j,k-1)
                v = y_speed(i,j,k-1)
                w = z_speed(i,j,k-1)
                ! ghost cell
                uf = x_speed_inf!x_speed(i,j,k)
                vf = y_speed_inf!y_speed(i,j,k)
                wf = z_speed_inf!z_speed(i,j,k)
                cexp = sqrt(gm*pressure(i,j,k-1)/density(i,j,k-1))
                !cinf = sqrt(gm*pressure(i,j,k)/density(i,j,k))
                cinf = sqrt(gm*pressure_inf/density_inf)
                Unexp = u *(znx(i,j,k)) + v *(zny(i,j,k)) + w *(znz(i,j,k))
                Uninf = uf*(znx(i,j,k)) + vf*(zny(i,j,k)) + wf*(znz(i,j,k))
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k) = x_speed(i,j,k-1) + vel_diff*(znx(i,j,k))
                  y_speed(i,j,k) = y_speed(i,j,k-1) + vel_diff*(zny(i,j,k))
                  z_speed(i,j,k) = z_speed(i,j,k-1) + vel_diff*(znz(i,j,k))
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face)
                      call copy3(tw, "flat", face)
                    case('kkl')
                      call copy3(tk, "flat", face)
                      call copy3(tkl, "flat", face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face)
                    case DEFAULT
                      continue
                  end select
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k) = x_speed_inf + vel_diff*(znx(i,j,k))
                  y_speed(i,j,k) = y_speed_inf + vel_diff*(zny(i,j,k))
                  z_speed(i,j,k) = z_speed_inf + vel_diff*(znz(i,j,k))
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, fixed_tv, face)
                    case('sst', 'sst2003')
                      call check_if_value_fixed("sst")
                      call fix(tk, fixed_tk, face)
                      call fix(tw, fixed_tw, face)
                    case('kkl')
                      call check_if_value_fixed("kkl")
                      call fix(tk, fixed_tk, face)
                      call fix(tkl, fixed_tkl, face)
                    case DEFAULT
                      !call turbulence_read_error()
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call check_if_value_fixed("lctm2015")
                      call fix(tgm, fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                end if
                Mb = sqrt(x_speed(i,j,k)**2+y_speed(i,j,k)**2+z_speed(i,j,k)**2)/Cb
                pressure(i,j,k) = fixed_Tpressure(6)/(((1+0.5*(gm-1.)*Mb*Mb))**(gm/(gm-1.)))
                density(i,j,k) = gm*pressure(i,j,k)/(Cb*Cb)
              end do
            end do
          end do
          qp(:,:,kmx+1,:) = qp(:,:,kmx,:)
          qp(:,:,kmx+2,:) = qp(:,:,kmx,:)
        case DEFAULT
          !print*, "ERROR: wrong face for boundary condition"
          Fatal_error
      end select

    end subroutine total_pressure

    subroutine temp_based_density(temperature, face)
      !< Specify the density in the ghost cell based on the
      !< temperature on the wall. Isothermal or adiabatic
      implicit none
      real, dimension(1:6)     , intent(in)  :: temperature
      character(len=*)         , intent(in)  :: face
      real :: stag_temp
      integer :: i,j,k

      select case(face)
        case("imin")
          if(temperature(1)<0.0)then
            do k = 1,kmx-1
              do j = 1,jmx-1
                do i = 1,1
                  stag_temp = (pressure(i,j,k)/(R_gas*density(i,j,k)))*(1 + (0.5*(gm-1.)*gm*pressure(i,j,k)/density(i,j,k)))
                  density(i-1,j,k) = pressure(i-1,j,k)/(R_gas*stag_temp)
                  density(i-2,j,k) = pressure(i-2,j,k)/(R_gas*stag_temp)
                  density(i-3,j,k) = pressure(i-3,j,k)/(R_gas*stag_temp)
                end do
              end do
            end do
          elseif(temperature(1)>1.0)then
            do k = 1,kmx-1
              do j = 1,jmx-1
                do i = 1,1
                  density(i-1,j,k) = pressure(i-1,j,k)/(R_gas*(2*temperature(1)-(pressure(i+0,j,k)/(R_gas*density(i+0,j,k)))))
                  density(i-2,j,k) = pressure(i-2,j,k)/(R_gas*(2*temperature(1)-(pressure(i+1,j,k)/(R_gas*density(i+1,j,k)))))
                  density(i-3,j,k) = pressure(i-3,j,k)/(R_gas*(2*temperature(1)-(pressure(i+2,j,k)/(R_gas*density(i+2,j,k)))))
                end do
              end do
            end do
          else
            call copy3(density , "symm",  face)
          end if
         case("imax")
          if(temperature(2)<0.0)then
            do k = 1,kmx-1
              do j = 1,jmx-1
                do i = imx-1,imx-1
                  stag_temp = (pressure(i,j,k)/(R_gas*density(i,j,k)))*(1 + (0.5*(gm-1.)*gm*pressure(i,j,k)/density(i,j,k)))
                  density(i+1,j,k) = pressure(i+1,j,k)/(R_gas*stag_temp)
                  density(i+2,j,k) = pressure(i+2,j,k)/(R_gas*stag_temp)
                  density(i+3,j,k) = pressure(i+3,j,k)/(R_gas*stag_temp)
                end do
              end do
            end do
          elseif(temperature(2)>1.0)then
            do k = 1,kmx-1
              do j = 1,jmx-1
                do i = imx-1,imx-1
                  density(i+1,j,k) = pressure(i+1,j,k)/(R_gas*(2*temperature(2)-(pressure(i-0,j,k)/(R_gas*density(i-0,j,k)))))
                  density(i+2,j,k) = pressure(i+2,j,k)/(R_gas*(2*temperature(2)-(pressure(i-1,j,k)/(R_gas*density(i-1,j,k)))))
                  density(i+3,j,k) = pressure(i+3,j,k)/(R_gas*(2*temperature(2)-(pressure(i-2,j,k)/(R_gas*density(i-2,j,k)))))
                end do
              end do
            end do
          else
            call copy3(density , "symm",  face)
          end if
        case("jmin")
          if(temperature(3)<0.0)then
            do k = 1,kmx-1
              do j = 1,1
                do i = 1,imx-1
                  stag_temp = (pressure(i,j,k)/(R_gas*density(i,j,k)))*(1 + (0.5*(gm-1.)*gm*pressure(i,j,k)/density(i,j,k)))
                  density(i,j-1,k) = pressure(i,j-1,k)/(R_gas*stag_temp)
                  density(i,j-2,k) = pressure(i,j-2,k)/(R_gas*stag_temp)
                  density(i,j-3,k) = pressure(i,j-3,k)/(R_gas*stag_temp)
                end do
              end do
            end do
          elseif(temperature(3)>1.0)then
            do k = 1,kmx-1
              do j = 1,1
                do i = 1,imx-1
                  density(i,j-1,k) = pressure(i,j-1,k)/(R_gas*(2*temperature(3)-(pressure(i,j+0,k)/(R_gas*density(i,j+0,k)))))
                  density(i,j-2,k) = pressure(i,j-2,k)/(R_gas*(2*temperature(3)-(pressure(i,j+1,k)/(R_gas*density(i,j+1,k)))))
                  density(i,j-3,k) = pressure(i,j-3,k)/(R_gas*(2*temperature(3)-(pressure(i,j+2,k)/(R_gas*density(i,j+2,k)))))
                end do
              end do
            end do
          else
            call copy3(density , "symm",  face)
          end if
        case("jmax")
          if(temperature(4)<0.0)then
            do k = 1,kmx-1
              do j = jmx-1,jmx-1
                do i = 1,imx-1
                  stag_temp = (pressure(i,j,k)/(R_gas*density(i,j,k)))*(1 + (0.5*(gm-1.)*gm*pressure(i,j,k)/density(i,j,k)))
                  density(i,j+1,k) = pressure(i,j+1,k)/(R_gas*stag_temp)
                  density(i,j+2,k) = pressure(i,j+2,k)/(R_gas*stag_temp)
                  density(i,j+3,k) = pressure(i,j+3,k)/(R_gas*stag_temp)
                end do
              end do
            end do
          elseif(temperature(4)>1.0)then
            do k = 1,kmx-1
              do j = jmx-1,jmx-1
                do i = 1,imx-1
                  density(i,j+1,k) = pressure(i,j+1,k)/(R_gas*(2*temperature(4)-(pressure(i,j-0,k)/(R_gas*density(i,j-0,k)))))
                  density(i,j+2,k) = pressure(i,j+2,k)/(R_gas*(2*temperature(4)-(pressure(i,j-1,k)/(R_gas*density(i,j-1,k)))))
                  density(i,j+3,k) = pressure(i,j+3,k)/(R_gas*(2*temperature(4)-(pressure(i,j-2,k)/(R_gas*density(i,j-2,k)))))
                end do
              end do
            end do
          else
            call copy3(density , "symm",  face)
          end if
        case("kmin")
          if(temperature(5)<0.0)then
            do k = 1,1
              do j = 1,jmx-1
                do i = 1,imx-1
                  stag_temp = (pressure(i,j,k)/(R_gas*density(i,j,k)))*(1 + (0.5*(gm-1.)*gm*pressure(i,j,k)/density(i,j,k)))
                  density(i,j,k-1) = pressure(i,j,k-1)/(R_gas*stag_temp)
                  density(i,j,k-2) = pressure(i,j,k-2)/(R_gas*stag_temp)
                  density(i,j,k-3) = pressure(i,j,k-3)/(R_gas*stag_temp)
                end do
              end do
            end do
          elseif(temperature(5)>1.0)then
            do k = 1,1
              do j = 1,jmx-1
                do i = 1,imx-1
                  density(i,j,k-1) = pressure(i,j,k-1)/(R_gas*(2*temperature(5)-(pressure(i,j,k+0)/(R_gas*density(i,j,k+0)))))
                  density(i,j,k-2) = pressure(i,j,k-2)/(R_gas*(2*temperature(5)-(pressure(i,j,k+1)/(R_gas*density(i,j,k+1)))))
                  density(i,j,k-3) = pressure(i,j,k-3)/(R_gas*(2*temperature(5)-(pressure(i,j,k+2)/(R_gas*density(i,j,k+2)))))
                end do
              end do
            end do
          else
            call copy3(density , "symm",  face)
          end if
        case("kmax")
          if(temperature(6)<0.0)then
            do k = kmx-1,kmx-1
              do j = 1,jmx-1
                do i = 1,imx-1
                  stag_temp = (pressure(i,j,k)/(R_gas*density(i,j,k)))*(1 + (0.5*(gm-1.)*gm*pressure(i,j,k)/density(i,j,k)))
                  density(i,j,k+1) = pressure(i,j,k+1)/(R_gas*stag_temp)
                  density(i,j,k+2) = pressure(i,j,k+2)/(R_gas*stag_temp)
                  density(i,j,k+3) = pressure(i,j,k+3)/(R_gas*stag_temp)
                end do
              end do
            end do
          elseif(temperature(6)>1.0)then
            do k = kmx-1,kmx-1
              do j = 1,jmx-1
                do i = 1,imx-1
                  density(i,j,k+1) = pressure(i,j,k+1)/(R_gas*(2*temperature(6)-(pressure(i,j,k-0)/(R_gas*density(i,j,k-0)))))
                  density(i,j,k+2) = pressure(i,j,k+2)/(R_gas*(2*temperature(6)-(pressure(i,j,k-1)/(R_gas*density(i,j,k-1)))))
                  density(i,j,k+3) = pressure(i,j,k+3)/(R_gas*(2*temperature(6)-(pressure(i,j,k-2)/(R_gas*density(i,j,k-2)))))
                end do
              end do
            end do
          else
            call copy3(density , "symm",  face)
          end if
        case DEFAULT
          !print*, "ERROR: wrong face for boundary condition"
          Fatal_error
      end select

    end subroutine temp_based_density


    subroutine periodic_bc(face)
      !< Single block periodic boundary condition.
      !< Not to be used for multiblock boundary condition
      implicit none
      character(len=*), intent(in) :: face

      select case(trim(face))

        case('imin')
          qp(-2:0,:,:,:) = qp(imx-3:imx-1,:,:,:)

        case('imax')
          qp(imx:imx+2,:,:,:) = qp(1:3,:,:,:)

        case('jmin')
          qp(:,-2:0,:,:) =   qp(:,jmx-3:jmx-1,:,:)

        case('jmax')
          qp(:,jmx:jmx+2,:,:) = qp(:,1:3,:,:)

        case('kmin')
          qp(:,:,-2:0,:) = qp(:,:,kmx-3:kmx-1,:)

        case('kmax')
          qp(:,:,kmx:kmx+2,:) = qp(:,:,1:3,:)

        case Default
          Fatal_error

      end select

    end subroutine periodic_bc
end module bc_primitive
