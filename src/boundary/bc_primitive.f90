  !< Apply boundary condition at every iteration
module bc_primitive
  !< Apply boundary condition at every iteration
  !-------------------------------------------
#include "../error.h"
  use vartypes
  use wall_dist, only: dist
  use global_sst , only: beta1
  use copy_bc   , only : copy3
  use FT_bc     , only : flow_tangency

  implicit none
  private

  integer                        :: face_num
  integer                        :: current_iter, imx, jmx, kmx, n_var
  !< Number of the face : 1:imin, 2:imax, 3:jmin, 4:jmax, 5:kmin, 6:kmax
  character(len=32) :: turbulence, transition
  real(wp) :: gm, R_gas, mu_ref,  T_ref, Sutherland_temp
  real(wp) :: x_speed_inf
  real(wp) :: y_speed_inf
  real(wp) :: z_speed_inf
  real(wp) :: density_inf
  real(wp) :: pressure_inf
  real(wp) :: tk_inf
  real(wp) :: tw_inf
  real(wp) :: te_inf
  real(wp) :: tv_inf
  real(wp) :: tgm_inf
  real(wp) :: tkl_inf
  real(wp), dimension(:, :, :, :), pointer :: qp
  real(wp), dimension(:, :, :), pointer :: density      
   !< Rho pointer, point to slice of qp (:,:,:,1)
  real(wp), dimension(:, :, :), pointer :: x_speed      
   !< U pointer, point to slice of qp (:,:,:,2) 
  real(wp), dimension(:, :, :), pointer :: y_speed      
   !< V pointer, point to slice of qp (:,:,:,3) 
  real(wp), dimension(:, :, :), pointer :: z_speed      
   !< W pointer, point to slice of qp (:,:,:,4)
  real(wp), dimension(:, :, :), pointer :: pressure     
   !< P pointer, point to slice of qp (:,:,:,5)
  ! state variable turbulent
  real(wp), dimension(:, :, :), pointer :: tk        !< TKE/mass
  real(wp), dimension(:, :, :), pointer :: tw        !< Omega
  real(wp), dimension(:, :, :), pointer :: te        !< Dissipation
  real(wp), dimension(:, :, :), pointer :: tv        !< SA visocity
  real(wp), dimension(:, :, :), pointer :: tkl       !< KL K-KL method
  real(wp), dimension(:, :, :), pointer :: tgm       !< Intermittency of LCTM2015

  public :: populate_ghost_primitive


  contains

    subroutine populate_ghost_primitive(state, Ifaces, Jfaces, Kfaces, control, scheme, flow, bc, dims)
      !< Populate the state variables in the ghost cell
      !< with particular value based on the boundary conditio 
      !< being applied at that face
      implicit none
      type(extent), intent(in) :: dims
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(inout), target :: state
      !< state variables
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      type(controltype), intent(in) :: control
      type(schemetype), intent(in) :: scheme
      type(flowtype), intent(in) :: flow
      type(boundarytype), intent(in) :: bc
      integer :: i
      character(len=4) :: face

      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx
      n_var = dims%n_var
      current_iter = control%current_iter
      turbulence = trim(scheme%turbulence)
      transition = trim(scheme%transition)
      mu_ref = flow%mu_ref
      gm = flow%gm
      R_gas = flow%R_gas
      T_ref = flow%T_ref
      sutherland_temp = flow%sutherland_temp
      x_speed_inf  =  flow%x_speed_inf 
      y_speed_inf  =  flow%y_speed_inf 
      z_speed_inf  =  flow%z_speed_inf 
      density_inf  =  flow%density_inf 
      pressure_inf =  flow%pressure_inf
      tk_inf       =  flow%tk_inf      
      tw_inf       =  flow%tw_inf      
      te_inf       =  flow%te_inf      
      tv_inf       =  flow%tv_inf      
      tgm_inf      =  flow%tgm_inf     
      tkl_inf      =  flow%tkl_inf     
     
      qp(-2:imx+2, -2:jmx+2, -2:kmx+2, 1:n_var) => state(:, :, :, :)
      density(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 1)
      x_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 2)
      y_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 3)
      z_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 4)
      pressure(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 5)


      select case (trim(scheme%turbulence))

          case ("none")
              !include nothing
              continue
          
          case ("sst", "sst2003", "bsl", "des-sst", "kw")
              tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 6)
              tw(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 7)

          case ("kkl")
              tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 6)
              tkl(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 7)

          case ("sa", "saBC")
              tv(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 6)

          case ("ke")
              tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 6)
              te(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 7)

          case ("les")
            continue
            ! todo

          case DEFAULT
            Fatal_error

      end select


      ! Transition modeling
      select case(trim(scheme%transition))

        case('lctm2015')
          tgm(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, n_var)
!          tgm_inf => qp_inf(n_var)

        case('bc', 'none')
          !do nothing
          continue

        case DEFAULT
          Fatal_error

      end Select
     
      
      do i = 1,6
        face_num = i
        face = bc%face_names(face_num)

        select case(bc%id(face_num))

          case(-1)
            call supersonic_inlet(face, bc)

          case(-2)
            call supersonic_outlet(face, bc, dims)

          case(-3)
            call subsonic_inlet(face, bc, dims)

          case(-4)
            call subsonic_outlet(face, bc, dims)

          case(-5)
            call wall(face, bc, dims)

          case(-6)
            call slip_wall(face, Ifaces, Jfaces, Kfaces, bc, dims)

          case(-7)
            call pole(face, bc, dims)

          case(-8)
            call far_field(face, Ifaces, Jfaces, Kfaces, bc, dims)

          case(-9)
            call periodic_bc(face)

          case(-11)
            call total_pressure(face, Ifaces, Jfaces, Kfaces, bc, dims)

          case Default
            if(bc%id(i)>=0 .or. bc%id(i)==-10) then
              continue !interface boundary 
            else
              print*, " boundary condition not recognised -> id is :", bc%id(i)
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


      subroutine supersonic_inlet(face, bc)
        !< Supersonic inlet boundary condition
        !< All the values of state variables are fixed
        implicit none
        character(len=*), intent(in) :: face
        type(boundarytype), intent(in) :: bc
        !< Name of the face at which boundary condition is called
        if(current_iter<=2)then
        call fix(density , bc%fixed_density , face)
        call fix(x_speed , bc%fixed_x_speed , face)
        call fix(y_speed , bc%fixed_y_speed , face)
        call fix(z_speed , bc%fixed_z_speed , face)
        call fix(pressure, bc%fixed_pressure, face)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sa', 'saBC')
            call fix(tv, bc%fixed_tv, face)
          case('sst', 'sst2003')
            !call check_if_value_fixed(bc, "sst")
            call fix(tk, bc%fixed_tk, face)
            call fix(tw, bc%fixed_tw, face)
          case('kkl')
            !call check_if_value_fixed(bc, "kkl")
            call fix(tk, bc%fixed_tk, face)
            call fix(tkl, bc%fixed_tkl, face)
          case DEFAULT
            Fatal_error
        end select
        select case(trim(transition))
          case('lctm2015')
            !call check_if_value_fixed(bc, "lctm2015")
            call fix(tgm, bc%fixed_tgm, face)
          case DEFAULT
            continue
        end select
        end if
      end subroutine supersonic_inlet


      subroutine supersonic_outlet(face, bc, dims)
        !< Supersonic outlet boundary condition. 
        !< All the values of state variables are copied 
        !< from inside the domain
        implicit none
        type(extent), intent(in) :: dims
        type(boundarytype), intent(in) :: bc
        character(len=*), intent(in) :: face
        !< Name of the face at which boundary condition is called
        call copy3(density , "flat", face, bc, dims)
        call copy3(x_speed , "flat", face, bc, dims)
        call copy3(y_speed , "flat", face, bc, dims)
        call copy3(z_speed , "flat", face, bc, dims)
        call copy3(pressure, "flat", face, bc, dims)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sa', 'saBC')
            call copy3(tv, "flat", face, bc, dims)
          case('sst', 'sst2003')
            call copy3(tk, "flat", face, bc, dims)
            call copy3(tw, "flat", face, bc, dims)
          case('kkl')
            call copy3(tk, "flat", face, bc, dims)
            call copy3(tkl, "flat", face, bc, dims)
          case DEFAULT
            Fatal_error
        end select
        select case(trim(transition))
          case('lctm2015')
            call copy3(tgm, "flat", face, bc, dims)
          case DEFAULT
            continue
        end select
      end subroutine supersonic_outlet


      subroutine subsonic_inlet(face, bc, dims)
        !< Subsonic inlet boundary condition. 
        !< All the state variables's value expect pressure
        !< is fixed and pressure is copied from inside the 
        !< domain
        implicit none
        type(extent), intent(in) :: dims
        type(boundarytype), intent(in) :: bc
        character(len=*), intent(in) :: face
        !< Name of the face at which boundary condition is called
        if(current_iter<=2)then
        call fix(density , bc%fixed_density , face)
        call fix(x_speed , bc%fixed_x_speed , face)
        call fix(y_speed , bc%fixed_y_speed , face)
        call fix(z_speed , bc%fixed_z_speed , face)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sa', 'saBC')
            call fix(tv, bc%fixed_tv, face)
          case('sst', 'sst2003')
            !call check_if_value_fixed(bc, "sst")
            call fix(tk, bc%fixed_tk, face)
            call fix(tw, bc%fixed_tw, face)
          case('kkl')
            !call check_if_value_fixed(bc, "kkl")
            call fix(tk, bc%fixed_tk, face)
            call fix(tkl, bc%fixed_tw, face)
          case DEFAULT
           Fatal_error
        end select
        select case(trim(transition))
          case('lctm2015')
            !call check_if_value_fixed(bc, "lctm2015")
            call fix(tgm, bc%fixed_tgm, face)
          case DEFAULT
            continue
        end select
        end if
        call copy3(pressure, "flat", face, bc, dims)
      end subroutine subsonic_inlet


      subroutine subsonic_outlet(face, bc, dims)
        !< Subsonic outlet boundary condition. 
        !< All the state variables's value expect pressure
        !< is copied from the inside of the domain and pressure 
        !< is fixed
        implicit none
        type(extent), intent(in) :: dims
        type(boundarytype), intent(in) :: bc
        character(len=*), intent(in) :: face
        !< Name of the face at which boundary condition is called
        call copy3(density, "flat", face, bc, dims)
        call copy3(x_speed, "flat", face, bc, dims)
        call copy3(y_speed, "flat", face, bc, dims)
        call copy3(z_speed, "flat", face, bc, dims)
        if(current_iter<=2)then
        call fix(pressure, bc%fixed_pressure, face)
        end if
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sa', 'saBC')
            call copy3(tv, "flat", face, bc, dims)
          case('sst', 'sst2003')
            call copy3(tk, "flat", face, bc, dims)
            call copy3(tw, "flat", face, bc, dims)
          case('kkl')
            call copy3(tk, "flat", face, bc, dims)
            call copy3(tkl, "flat", face, bc, dims)
          case DEFAULT
           Fatal_error
        end select
        select case(trim(transition))
          case('lctm2015')
            call copy3(tgm, "flat", face, bc, dims)
          case DEFAULT
            continue
        end select
      end subroutine subsonic_outlet

      subroutine wall(face, bc, dims)
        !< Adiabatic/Isothermal wall boundary condition
        implicit none
        type(extent), intent(in) :: dims
        character(len=*), intent(in) :: face
        type(boundarytype), intent(in) :: bc
        !< Name of the face at which boundary condition is called
        call copy3(pressure, "symm",  face, bc, dims)
        call temp_based_density(bc%fixed_wall_temperature, face, bc, dims)
        call no_slip(face, bc, dims)
      end subroutine wall


      subroutine slip_wall(face, Ifaces, Jfaces, Kfaces, bc, dims)
        !< Slip wall boundary condition. 
        !< Maintain flow tangency
        implicit none
        type(extent), intent(in) :: dims
        type(boundarytype), intent(in) :: bc
        character(len=*), intent(in) :: face
        type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
        !< Input varaible which stores I faces' area and unit normal
        type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
        !< Input varaible which stores J faces' area and unit normal
        type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
        !< Input varaible which stores K faces' area and unit normal
        !< Name of the face at which boundary condition is called
        call copy3(density , "symm", face, bc, dims)
        call copy3(pressure, "symm", face, bc, dims)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sa', 'saBC')
            call copy3(tv, "symm", face, bc, dims)
          case('sst', 'sst2003')
            call copy3(tk, "symm", face, bc, dims)
            call copy3(tw, "symm", face, bc, dims)
          case('kkl')
            call copy3(tk, "symm", face, bc, dims)
            call copy3(tkl, "symm", face, bc, dims)
          case DEFAULT
            Fatal_error
        end select
        select case(trim(transition))
          case('lctm2015')
            call copy3(tgm, "flat", face, bc, dims)
          case DEFAULT
            continue
        end select
        call flow_tangency(qp, face, Ifaces, Jfaces, Kfaces, dims)
      end subroutine slip_wall


      subroutine pole(face, bc, dims)
        !< Boundary condition for the block face
        !< with zero area; turning into a pole
        implicit none
        type(extent), intent(in) :: dims
        type(boundarytype), intent(in) :: bc
        character(len=*), intent(in) :: face
        !< Name of the face at which boundary condition is called
        call copy3(density , "flat", face, bc, dims)
        call copy3(x_speed , "flat", face, bc, dims)
        call copy3(y_speed , "flat", face, bc, dims)
        call copy3(z_speed , "flat", face, bc, dims)
        call copy3(pressure, "flat", face, bc, dims)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sa', 'saBC')
            call copy3(tv, "flat", face, bc, dims) 
          case('sst', 'sst2003')
            call copy3(tk, "flat", face, bc, dims)
            call copy3(tw, "flat", face, bc, dims)
          case('kkl')
            call copy3(tk, "flat", face, bc, dims)
            call copy3(tkl, "flat", face, bc, dims)
          case DEFAULT
            Fatal_error
        end select
        select case(trim(transition))
          case('lctm2015')
            call copy3(tgm, "flat", face, bc, dims)
          case DEFAULT
            continue
        end select
      end subroutine pole



      subroutine fix(var, fix_val, face)
        !< Generalized subroutine to fix particular value
        !< at particular face
        implicit none
        real(wp), dimension(-2:imx+2, -2:jmx+2, -2:kmx+2) , intent(out) :: var
        !< Variable of which values are being fixed in the ghost cell
        real(wp), dimension(1:6)       , intent(in)  :: fix_val
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


      subroutine no_slip(face, bc, dims)
        !< No-slip wall boundary condition. All the 
        !< component of velocity throught face is zero
        implicit none
        type(extent), intent(in) :: dims
        type(boundarytype), intent(in) :: bc
        character(len=*), intent(in) :: face
        !< Name of the face at which boundary condition is called
        call copy3(x_speed, "anti", face, bc, dims)
        call copy3(y_speed, "anti", face, bc, dims)
        call copy3(z_speed, "anti", face, bc, dims)
        select case(turbulence)
          case("none")
            !do nothing
            continue
          case('sa', 'saBC')
            call copy3(tv  , "anti", face, bc, dims)
          case("sst", 'sst2003')
            call copy3(tk  , "anti", face, bc, dims)
            call set_omega_at_wall(face)
          case("kkl")
            call copy3(tk  , "anti", face, bc, dims)
            call copy3(tkl , "anti", face, bc, dims)
          case DEFAULT
            Fatal_error
        end select
        select case(trim(transition))
          case('lctm2015')
            call copy3(tgm, "flat", face, bc, dims)
          case DEFAULT
            continue
        end select
      end subroutine no_slip

    

      subroutine set_omega_at_wall(face)
        !< Set value of turbulence variable: omega (turbulenct dissipation rate). 
        !< Value fixed is accourding to the SST turbulence model
        implicit none
        character(len=*), intent(in) :: face
        real(wp) :: T_face
        real(wp) :: mu
        real(wp) :: rho
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

    subroutine far_field(face, Ifaces, Jfaces, Kfaces, bc, dims)
      !< Far-field Riemann boundary condition
      implicit none
      type(extent), intent(in) :: dims
      character(len=*), intent(in) :: face
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      type(boundarytype), intent(in) :: bc
      real(wp) :: cinf, cexp   ! speed of sound
      real(wp) :: Rinf, Rexp   ! Riemann invarient
      real(wp) :: Uninf, Unexp ! face normal speed
      real(wp) :: Unb ! normal velocity boundary
      real(wp) :: Cb  ! speed of sound boundary
      real(wp) :: vel_diff
      real(wp) :: u,v,w
      real(wp) :: uf, vf, wf
      integer :: i,j,k
      real(wp) :: s
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
                Unexp = u *(-Ifaces(i,j,k)%nx) + v *(-Ifaces(i,j,k)%ny) + w *(-Ifaces(i,j,k)%nz)
                Uninf = uf*(-Ifaces(i,j,k)%nx) + vf*(-Ifaces(i,j,k)%ny) + wf*(-Ifaces(i,j,k)%nz)
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i-1,j,k) = x_speed(i,j,k) + vel_diff*(-Ifaces(i,j,k)%nx)
                  y_speed(i-1,j,k) = y_speed(i,j,k) + vel_diff*(-Ifaces(i,j,k)%ny)
                  z_speed(i-1,j,k) = z_speed(i,j,k) + vel_diff*(-Ifaces(i,j,k)%nz)
                  s = pressure(i,j,k)/(density(i,j,k)**(gm))
                  density(i-1,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i-1,j,k) = (density(i-1,j,k)*Cb*Cb/gm)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face, bc, dims)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tw, "flat", face, bc, dims)
                    case('kkl')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tkl, "flat", face, bc, dims)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face, bc, dims)
                    case DEFAULT
                      continue
                  end select
                  face_already_has_fixed_values(1)=0
                else
                  vel_diff = Unb - Uninf
                  x_speed(i-1,j,k) = x_speed_inf + vel_diff*(-Ifaces(i,j,k)%nx)
                  y_speed(i-1,j,k) = y_speed_inf + vel_diff*(-Ifaces(i,j,k)%ny)
                  z_speed(i-1,j,k) = z_speed_inf + vel_diff*(-Ifaces(i,j,k)%nz)
                  s = pressure_inf/(density_inf**(gm))
                  density(i-1,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i-1,j,k) = (density(i-1,j,k)*Cb*Cb/gm)

                  if(face_already_has_fixed_values(1)==0)then
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, bc%fixed_tv, face)
                    case('sst', 'sst2003')
                      !call check_if_value_fixed(bc, "sst")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tw, bc%fixed_tw, face)
                    case('kkl')
                      !call check_if_value_fixed(bc, "kkl")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tkl, bc%fixed_tkl, face)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      !call check_if_value_fixed(bc, "lctm2015")
                      call fix(tgm, bc%fixed_tgm, face)
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
                Unexp = u *(Ifaces(i,j,k)%nx) + v *(Ifaces(i,j,k)%ny) + w *(Ifaces(i,j,k)%nz)
                Uninf = uf*(Ifaces(i,j,k)%nx) + vf*(Ifaces(i,j,k)%ny) + wf*(Ifaces(i,j,k)%nz)
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k) = x_speed(i-1,j,k) + vel_diff*(Ifaces(i,j,k)%nx)
                  y_speed(i,j,k) = y_speed(i-1,j,k) + vel_diff*(Ifaces(i,j,k)%ny)
                  z_speed(i,j,k) = z_speed(i-1,j,k) + vel_diff*(Ifaces(i,j,k)%nz)
                  s = pressure(i-1,j,k)/(density(i-1,j,k)**(gm))
                  density(i,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k) = (density(i,j,k)*Cb*Cb/gm)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face, bc, dims)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tw, "flat", face, bc, dims)
                    case('kkl')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tkl, "flat", face, bc, dims)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face, bc, dims)
                    case DEFAULT
                      continue
                  end select
                  face_already_has_fixed_values(2)=0
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k) = x_speed_inf + vel_diff*(Ifaces(i,j,k)%nx)
                  y_speed(i,j,k) = y_speed_inf + vel_diff*(Ifaces(i,j,k)%ny)
                  z_speed(i,j,k) = z_speed_inf + vel_diff*(Ifaces(i,j,k)%nz)
                  s = pressure_inf/(density_inf**(gm))
                  density(i,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k) = (density(i,j,k)*Cb*Cb/gm)

                  if(face_already_has_fixed_values(2)==0)then
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, bc%fixed_tv, face)
                    case('sst', 'sst2003')
                      !call check_if_value_fixed(bc, "sst")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tw, bc%fixed_tw, face)
                    case('kkl')
                      !call check_if_value_fixed(bc, "kkl")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tkl, bc%fixed_tkl, face)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      !call check_if_value_fixed(bc, "lctm2015")
                      call fix(tgm, bc%fixed_tgm, face)
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
                Unexp = u *(-Jfaces(i,j,k)%nx) + v *(-Jfaces(i,j,k)%ny) + w *(-Jfaces(i,j,k)%nz)
                Uninf = uf*(-Jfaces(i,j,k)%nx) + vf*(-Jfaces(i,j,k)%ny) + wf*(-Jfaces(i,j,k)%nz)
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j-1,k) = x_speed(i,j,k) + vel_diff*(-Jfaces(i,j,k)%nx)
                  y_speed(i,j-1,k) = y_speed(i,j,k) + vel_diff*(-Jfaces(i,j,k)%ny)
                  z_speed(i,j-1,k) = z_speed(i,j,k) + vel_diff*(-Jfaces(i,j,k)%nz)
                  s = pressure(i,j,k)/(density(i,j,k)**(gm))
                  density(i,j-1,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j-1,k) = (density(i,j-1,k)*Cb*Cb/gm)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face, bc, dims)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tw, "flat", face, bc, dims)
                    case('kkl')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tkl, "flat", face, bc, dims)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face, bc, dims)
                    case DEFAULT
                      continue
                  end select
                  face_already_has_fixed_values(3)=0
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j-1,k) = x_speed_inf + vel_diff*(-Jfaces(i,j,k)%nx)
                  y_speed(i,j-1,k) = y_speed_inf + vel_diff*(-Jfaces(i,j,k)%ny)
                  z_speed(i,j-1,k) = z_speed_inf + vel_diff*(-Jfaces(i,j,k)%nz)
                  s = pressure_inf/(density_inf**(gm))
                  density(i,j-1,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j-1,k) = (density(i,j-1,k)*Cb*Cb/gm)

                  if(face_already_has_fixed_values(3)==0)then
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, bc%fixed_tv, face)
                    case('sst', 'sst2003')
                      !call check_if_value_fixed(bc, "sst")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tw, bc%fixed_tw, face)
                    case('kkl')
                      !call check_if_value_fixed(bc, "kkl")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tkl, bc%fixed_tkl, face)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      !call check_if_value_fixed(bc, "lctm2015")
                      call fix(tgm, bc%fixed_tgm, face)
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
                Unexp = u *(Jfaces(i,j,k)%nx) + v *(Jfaces(i,j,k)%ny) + w *(Jfaces(i,j,k)%nz)
                Uninf = uf*(Jfaces(i,j,k)%nx) + vf*(Jfaces(i,j,k)%ny) + wf*(Jfaces(i,j,k)%nz)
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k) = x_speed(i,j-1,k) + vel_diff*(Jfaces(i,j,k)%nx)
                  y_speed(i,j,k) = y_speed(i,j-1,k) + vel_diff*(Jfaces(i,j,k)%ny)
                  z_speed(i,j,k) = z_speed(i,j-1,k) + vel_diff*(Jfaces(i,j,k)%nz)
                  s = pressure(i,j-1,k)/(density(i,j-1,k)**(gm))
                  density(i,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k) = (density(i,j,k)*Cb*Cb/gm)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face, bc, dims)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tw, "flat", face, bc, dims)
                    case('kkl')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tkl, "flat", face, bc, dims)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face, bc, dims)
                    case DEFAULT
                      continue
                  end select
                  face_already_has_fixed_values(4)=0
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k) = x_speed_inf + vel_diff*(Jfaces(i,j,k)%nx)
                  y_speed(i,j,k) = y_speed_inf + vel_diff*(Jfaces(i,j,k)%ny)
                  z_speed(i,j,k) = z_speed_inf + vel_diff*(Jfaces(i,j,k)%nz)
                  s = pressure_inf/(density_inf**(gm))
                  density(i,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k) = (density(i,j,k)*Cb*Cb/gm)

                  if(face_already_has_fixed_values(4)==0)then
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, bc%fixed_tv, face)
                    case('sst', 'sst2003')
                      !call check_if_value_fixed(bc, "sst")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tw, bc%fixed_tw, face)
                    case('kkl')
                      !call check_if_value_fixed(bc, "kkl")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tkl, bc%fixed_tkl, face)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      !call check_if_value_fixed(bc, "lctm2015")
                      call fix(tgm, bc%fixed_tgm, face)
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
                Unexp = u *(-Kfaces(i,j,k)%nx) + v *(-Kfaces(i,j,k)%ny) + w *(-Kfaces(i,j,k)%nz)
                Uninf = uf*(-Kfaces(i,j,k)%nx) + vf*(-Kfaces(i,j,k)%ny) + wf*(-Kfaces(i,j,k)%nz)
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k-1) = x_speed(i,j,k) + vel_diff*(-Kfaces(i,j,k)%nx)
                  y_speed(i,j,k-1) = y_speed(i,j,k) + vel_diff*(-Kfaces(i,j,k)%ny)
                  z_speed(i,j,k-1) = z_speed(i,j,k) + vel_diff*(-Kfaces(i,j,k)%nz)
                  s = pressure(i,j,k)/(density(i,j,k)**(gm))
                  density(i,j,k-1) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k-1) = (density(i,j,k-1)*Cb*Cb/gm)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face, bc, dims)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tw, "flat", face, bc, dims)
                    case('kkl')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tkl, "flat", face, bc, dims)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face, bc, dims)
                    case DEFAULT
                      continue
                  end select
                  face_already_has_fixed_values(5)=0
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k-1) = x_speed_inf + vel_diff*(-Kfaces(i,j,k)%nx)
                  y_speed(i,j,k-1) = y_speed_inf + vel_diff*(-Kfaces(i,j,k)%ny)
                  z_speed(i,j,k-1) = z_speed_inf + vel_diff*(-Kfaces(i,j,k)%nz)
                  s = pressure_inf/(density_inf**(gm))
                  density(i,j,k-1) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k-1) = (density(i,j,k-1)*Cb*Cb/gm)

                  if(face_already_has_fixed_values(5)==0)then
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, bc%fixed_tv, face)
                    case('sst', 'sst2003')
                      !call check_if_value_fixed(bc, "sst")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tw, bc%fixed_tw, face)
                    case('kkl')
                      !call check_if_value_fixed(bc, "kkl")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tkl, bc%fixed_tkl, face)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      !call check_if_value_fixed(bc, "lctm2015")
                      call fix(tgm, bc%fixed_tgm, face)
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
                Unexp = u *(Kfaces(i,j,k)%nx) + v *(Kfaces(i,j,k)%ny) + w *(Kfaces(i,j,k)%nz)
                Uninf = uf*(Kfaces(i,j,k)%nx) + vf*(Kfaces(i,j,k)%ny) + wf*(Kfaces(i,j,k)%nz)
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k) = x_speed(i,j,k-1) + vel_diff*(Kfaces(i,j,k)%nx)
                  y_speed(i,j,k) = y_speed(i,j,k-1) + vel_diff*(Kfaces(i,j,k)%ny)
                  z_speed(i,j,k) = z_speed(i,j,k-1) + vel_diff*(Kfaces(i,j,k)%nz)
                  s = pressure(i,j,k-1)/(density(i,j,k-1)**(gm))
                  density(i,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k) = (density(i,j,k)*Cb*Cb/gm)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face, bc, dims)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tw, "flat", face, bc, dims)
                    case('kkl')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tkl, "flat", face, bc, dims)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face, bc, dims)
                    case DEFAULT
                      continue
                  end select
                  face_already_has_fixed_values(6)=0
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k) = x_speed_inf + vel_diff*(Kfaces(i,j,k)%nx)
                  y_speed(i,j,k) = y_speed_inf + vel_diff*(Kfaces(i,j,k)%ny)
                  z_speed(i,j,k) = z_speed_inf + vel_diff*(Kfaces(i,j,k)%nz)
                  s = pressure_inf/(density_inf**(gm))
                  density(i,j,k) = (Cb*Cb/(gm*s))**(1./(gm-1.))
                  pressure(i,j,k) = (density(i,j,k)*Cb*Cb/gm)

                  if(face_already_has_fixed_values(6)==0)then
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, bc%fixed_tv, face)
                    case('sst', 'sst2003')
                      !call check_if_value_fixed(bc, "sst")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tw, bc%fixed_tw, face)
                    case('kkl')
                      !call check_if_value_fixed(bc, "kkl")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tkl, bc%fixed_tkl, face)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      !call check_if_value_fixed(bc, "lctm2015")
                      call fix(tgm, bc%fixed_tgm, face)
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


    subroutine total_pressure(face, Ifaces, Jfaces, Kfaces, bc, dims)
      !< Total Pressure Riemann boundary condition
      implicit none
      type(extent), intent(in) :: dims
      character(len=*), intent(in) :: face
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      type(boundarytype), intent(in) :: bc
      real(wp) :: cinf, cexp   ! speed of sound
      real(wp) :: Rinf, Rexp   ! Riemann invarient
      real(wp) :: Uninf, Unexp ! face normal speed
      real(wp) :: Unb ! normal velocity boundary
      real(wp) :: Cb  ! speed of sound boundary
      real(wp) :: vel_diff
      real(wp) :: u,v,w
      real(wp) :: uf, vf, wf
      real(wp) :: Mb
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
                Unexp = u *(-Ifaces(i,j,k)%nx) + v *(-Ifaces(i,j,k)%ny) + w *(-Ifaces(i,j,k)%nz)
                Uninf = uf*(-Ifaces(i,j,k)%nx) + vf*(-Ifaces(i,j,k)%ny) + wf*(-Ifaces(i,j,k)%nz)
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i-1,j,k) = x_speed(i,j,k) + vel_diff*(-Ifaces(i,j,k)%nx)
                  y_speed(i-1,j,k) = y_speed(i,j,k) + vel_diff*(-Ifaces(i,j,k)%ny)
                  z_speed(i-1,j,k) = z_speed(i,j,k) + vel_diff*(-Ifaces(i,j,k)%nz)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face, bc, dims)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tw, "flat", face, bc, dims)
                    case('kkl')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tkl, "flat", face, bc, dims)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face, bc, dims)
                    case DEFAULT
                      continue
                  end select
                else
                  vel_diff = Unb - Uninf
                  x_speed(i-1,j,k) = x_speed_inf + vel_diff*(-Ifaces(i,j,k)%nx)
                  y_speed(i-1,j,k) = y_speed_inf + vel_diff*(-Ifaces(i,j,k)%ny)
                  z_speed(i-1,j,k) = z_speed_inf + vel_diff*(-Ifaces(i,j,k)%nz)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, bc%fixed_tv, face)
                    case('sst', 'sst2003')
                      !call check_if_value_fixed(bc, "sst")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tw, bc%fixed_tw, face)
                    case('kkl')
                      !call check_if_value_fixed(bc, "kkl")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tkl, bc%fixed_tkl, face)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      !call check_if_value_fixed(bc, "lctm2015")
                      call fix(tgm, bc%fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                end if
                Mb = sqrt(x_speed(i-1,j,k)**2+y_speed(i-1,j,k)**2+z_speed(i-1,j,k)**2)/Cb
                pressure(i-1,j,k) = bc%fixed_Tpressure(1)/(((1+0.5*(gm-1.)*Mb*Mb))**(gm/(gm-1.)))
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
                Unexp = u *(Ifaces(i,j,k)%nx) + v *(Ifaces(i,j,k)%ny) + w *(Ifaces(i,j,k)%nz)
                Uninf = uf*(Ifaces(i,j,k)%nx) + vf*(Ifaces(i,j,k)%ny) + wf*(Ifaces(i,j,k)%nz)
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k) = x_speed(i-1,j,k) + vel_diff*(Ifaces(i,j,k)%nx)
                  y_speed(i,j,k) = y_speed(i-1,j,k) + vel_diff*(Ifaces(i,j,k)%ny)
                  z_speed(i,j,k) = z_speed(i-1,j,k) + vel_diff*(Ifaces(i,j,k)%nz)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face, bc, dims)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tw, "flat", face, bc, dims)
                    case('kkl')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tkl, "flat", face, bc, dims)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face, bc, dims)
                    case DEFAULT
                      continue
                  end select
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k) = x_speed_inf + vel_diff*(Ifaces(i,j,k)%nx)
                  y_speed(i,j,k) = y_speed_inf + vel_diff*(Ifaces(i,j,k)%ny)
                  z_speed(i,j,k) = z_speed_inf + vel_diff*(Ifaces(i,j,k)%nz)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, bc%fixed_tv, face)
                    case('sst', 'sst2003')
                      !call check_if_value_fixed(bc, "sst")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tw, bc%fixed_tw, face)
                    case('kkl')
                      !call check_if_value_fixed(bc, "kkl")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tkl, bc%fixed_tkl, face)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      !call check_if_value_fixed(bc, "lctm2015")
                      call fix(tgm, bc%fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                end if
                Mb = sqrt(x_speed(i,j,k)**2+y_speed(i,j,k)**2+z_speed(i,j,k)**2)/Cb
                pressure(i,j,k) = bc%fixed_Tpressure(2)/(((1+0.5*(gm-1.)*Mb*Mb))**(gm/(gm-1.)))
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
                Unexp = u *(-Jfaces(i,j,k)%nx) + v *(-Jfaces(i,j,k)%ny) + w *(-Jfaces(i,j,k)%nz)
                Uninf = uf*(-Jfaces(i,j,k)%nx) + vf*(-Jfaces(i,j,k)%ny) + wf*(-Jfaces(i,j,k)%nz)
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j-1,k) = x_speed(i,j,k) + vel_diff*(-Jfaces(i,j,k)%nx)
                  y_speed(i,j-1,k) = y_speed(i,j,k) + vel_diff*(-Jfaces(i,j,k)%ny)
                  z_speed(i,j-1,k) = z_speed(i,j,k) + vel_diff*(-Jfaces(i,j,k)%nz)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face, bc, dims)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tw, "flat", face, bc, dims)
                    case('kkl')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tkl, "flat", face, bc, dims)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face, bc, dims)
                    case DEFAULT
                      continue
                  end select
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j-1,k) = x_speed_inf + vel_diff*(-Jfaces(i,j,k)%nx)
                  y_speed(i,j-1,k) = y_speed_inf + vel_diff*(-Jfaces(i,j,k)%ny)
                  z_speed(i,j-1,k) = z_speed_inf + vel_diff*(-Jfaces(i,j,k)%nz)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, bc%fixed_tv, face)
                    case('sst', 'sst2003')
                      !call check_if_value_fixed(bc, "sst")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tw, bc%fixed_tw, face)
                    case('kkl')
                      !call check_if_value_fixed(bc, "kkl")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tkl, bc%fixed_tkl, face)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      !call check_if_value_fixed(bc, "lctm2015")
                      call fix(tgm, bc%fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                end if
                Mb = sqrt(x_speed(i,j-1,k)**2+y_speed(i,j-1,k)**2+z_speed(i,j-1,k)**2)/Cb
                pressure(i,j-1,k) = bc%fixed_Tpressure(3)/(((1+0.5*(gm-1.)*Mb*Mb))**(gm/(gm-1.)))
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
                Unexp = u *(Jfaces(i,j,k)%nx) + v *(Jfaces(i,j,k)%ny) + w *(Jfaces(i,j,k)%nz)
                Uninf = uf*(Jfaces(i,j,k)%nx) + vf*(Jfaces(i,j,k)%ny) + wf*(Jfaces(i,j,k)%nz)
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k) = x_speed(i,j-1,k) + vel_diff*(Jfaces(i,j,k)%nx)
                  y_speed(i,j,k) = y_speed(i,j-1,k) + vel_diff*(Jfaces(i,j,k)%ny)
                  z_speed(i,j,k) = z_speed(i,j-1,k) + vel_diff*(Jfaces(i,j,k)%nz)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face, bc, dims)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tw, "flat", face, bc, dims)
                    case('kkl')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tkl, "flat", face, bc, dims)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face, bc, dims)
                    case DEFAULT
                      continue
                  end select
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k) = x_speed_inf + vel_diff*(Jfaces(i,j,k)%nx)
                  y_speed(i,j,k) = y_speed_inf + vel_diff*(Jfaces(i,j,k)%ny)
                  z_speed(i,j,k) = z_speed_inf + vel_diff*(Jfaces(i,j,k)%nz)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, bc%fixed_tv, face)
                    case('sst', 'sst2003')
                      !call check_if_value_fixed(bc, "sst")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tw, bc%fixed_tw, face)
                    case('kkl')
                      !call check_if_value_fixed(bc, "kkl")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tkl, bc%fixed_tkl, face)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      !call check_if_value_fixed(bc, "lctm2015")
                      call fix(tgm, bc%fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                end if
                Mb = sqrt(x_speed(i,j,k)**2+y_speed(i,j,k)**2+z_speed(i,j,k)**2)/Cb
                pressure(i,j,k) = bc%fixed_Tpressure(4)/(((1+0.5*(gm-1.)*Mb*Mb))**(gm/(gm-1.)))
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
                Unexp = u *(-Kfaces(i,j,k)%nx) + v *(-Kfaces(i,j,k)%ny) + w *(-Kfaces(i,j,k)%nz)
                Uninf = uf*(-Kfaces(i,j,k)%nx) + vf*(-Kfaces(i,j,k)%ny) + wf*(-Kfaces(i,j,k)%nz)
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k-1) = x_speed(i,j,k) + vel_diff*(-Kfaces(i,j,k)%nx)
                  y_speed(i,j,k-1) = y_speed(i,j,k) + vel_diff*(-Kfaces(i,j,k)%ny)
                  z_speed(i,j,k-1) = z_speed(i,j,k) + vel_diff*(-Kfaces(i,j,k)%nz)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face, bc, dims)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tw, "flat", face, bc, dims)
                    case('kkl')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tkl, "flat", face, bc, dims)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face, bc, dims)
                    case DEFAULT
                      continue
                  end select
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k-1) = x_speed_inf + vel_diff*(-Kfaces(i,j,k)%nx)
                  y_speed(i,j,k-1) = y_speed_inf + vel_diff*(-Kfaces(i,j,k)%ny)
                  z_speed(i,j,k-1) = z_speed_inf + vel_diff*(-Kfaces(i,j,k)%nz)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, bc%fixed_tv, face)
                    case('sst', 'sst2003')
                      !call check_if_value_fixed(bc, "sst")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tw, bc%fixed_tw, face)
                    case('kkl')
                      !call check_if_value_fixed(bc, "kkl")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tkl, bc%fixed_tkl, face)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      !call check_if_value_fixed(bc, "lctm2015")
                      call fix(tgm, bc%fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                end if
                Mb = sqrt(x_speed(i,j,k)**2+y_speed(i,j,k)**2+z_speed(i,j,k)**2)/Cb
                pressure(i,j,k-1) = bc%fixed_Tpressure(5)/(((1+0.5*(gm-1.)*Mb*Mb))**(gm/(gm-1.)))
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
                Unexp = u *(Kfaces(i,j,k)%nx) + v *(Kfaces(i,j,k)%ny) + w *(Kfaces(i,j,k)%nz)
                Uninf = uf*(Kfaces(i,j,k)%nx) + vf*(Kfaces(i,j,k)%ny) + wf*(Kfaces(i,j,k)%nz)
                Rinf  = Uninf - 2*cinf/(gm-1.)
                Rexp  = Unexp + 2*cexp/(gm-1.)
                Unb   = 0.5*(Rexp + Rinf)
                Cb    = 0.25*(gm-1.)*(Rexp - Rinf)
                if(Unb > 0.)then
                  vel_diff = Unb - Unexp
                  x_speed(i,j,k) = x_speed(i,j,k-1) + vel_diff*(Kfaces(i,j,k)%nx)
                  y_speed(i,j,k) = y_speed(i,j,k-1) + vel_diff*(Kfaces(i,j,k)%ny)
                  z_speed(i,j,k) = z_speed(i,j,k-1) + vel_diff*(Kfaces(i,j,k)%nz)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call copy3(tv, "flat", face, bc, dims)
                    case('sst', 'sst2003')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tw, "flat", face, bc, dims)
                    case('kkl')
                      call copy3(tk, "flat", face, bc, dims)
                      call copy3(tkl, "flat", face, bc, dims)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      call copy3(tgm, "flat", face, bc, dims)
                    case DEFAULT
                      continue
                  end select
                else
                  vel_diff = Unb - Uninf
                  x_speed(i,j,k) = x_speed_inf + vel_diff*(Kfaces(i,j,k)%nx)
                  y_speed(i,j,k) = y_speed_inf + vel_diff*(Kfaces(i,j,k)%ny)
                  z_speed(i,j,k) = z_speed_inf + vel_diff*(Kfaces(i,j,k)%nz)
                  select case (turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sa', 'saBC')
                      call fix(tv, bc%fixed_tv, face)
                    case('sst', 'sst2003')
                      !call check_if_value_fixed(bc, "sst")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tw, bc%fixed_tw, face)
                    case('kkl')
                      !call check_if_value_fixed(bc, "kkl")
                      call fix(tk, bc%fixed_tk, face)
                      call fix(tkl, bc%fixed_tkl, face)
                    case DEFAULT
                      Fatal_error
                  end select
                  select case(trim(transition))
                    case('lctm2015')
                      !call check_if_value_fixed(bc, "lctm2015")
                      call fix(tgm, bc%fixed_tgm, face)
                    case DEFAULT
                      continue
                  end select
                end if
                Mb = sqrt(x_speed(i,j,k)**2+y_speed(i,j,k)**2+z_speed(i,j,k)**2)/Cb
                pressure(i,j,k) = bc%fixed_Tpressure(6)/(((1+0.5*(gm-1.)*Mb*Mb))**(gm/(gm-1.)))
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

    subroutine temp_based_density(temperature, face, bc, dims)
      !< Specify the density in the ghost cell based on the
      !< temperature on the wall. Isothermal or adiabatic
      implicit none
      type(extent), intent(in) :: dims
      type(boundarytype), intent(in) :: bc
      real(wp), dimension(1:6)     , intent(in)  :: temperature
      character(len=*)         , intent(in)  :: face
      real(wp) :: stag_temp
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
            call copy3(density , "symm",  face, bc, dims)
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
            call copy3(density , "symm",  face, bc, dims)
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
            call copy3(density , "symm",  face, bc, dims)
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
            call copy3(density , "symm",  face, bc, dims)
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
            call copy3(density , "symm",  face, bc, dims)
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
            call copy3(density , "symm",  face, bc, dims)
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
