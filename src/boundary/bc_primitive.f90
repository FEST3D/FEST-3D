module bc_primitive
  !--------------------------------------------
  ! 170515  Jatinder Pal Singh Sandhu
  ! Aim : applying boundary condition to domain
  !-------------------------------------------
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
  use global_vars, only: face_names
  use global_vars, only: id
  use global_vars, only: c1
  use global_vars, only: c2
  use global_vars, only: c3

  use global_sst , only: beta1
  use utils,       only: turbulence_read_error

  use read_bc   , only: read_fixed_values

  implicit none
  private

  integer                        :: face_num

  public :: populate_ghost_primitive
  public :: set_wall_flux_zero
!  public :: setup_bc


  contains

!    subroutine setup_bc()
!      implicit none
!      face_names(1) = "imin"
!      face_names(2) = "imax"
!      face_names(3) = "jmin"
!      face_names(4) = "jmax"
!      face_names(5) = "kmin"
!      face_names(6) = "kmax"
!      
!      id(1) =  imin_id
!      id(2) =  imax_id
!      id(3) =  jmin_id
!      id(4) =  jmax_id
!      id(5) =  kmin_id
!      id(6) =  kmax_id
!
!      c2 = 1 + accur
!      c3 = 0.5*accur
!      c1 = c2-c3
!      call read_fixed_values()
!
!    end subroutine setup_bc


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
          case DEFAULT
            call turbulence_read_error()
        end select
      end subroutine supersonic_inlet


      subroutine supersonic_outlet(face)
        implicit none
        character(len=*), intent(in) :: face
        call copy(density , "flat", face)
        call copy(x_speed , "flat", face)
        call copy(y_speed , "flat", face)
        call copy(z_speed , "flat", face)
        call copy(pressure, "flat", face)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sst')
            call copy(tk, "flat", face)
            call copy(tw, "flat", face)
          case DEFAULT
            call turbulence_read_error()
        end select
      end subroutine supersonic_outlet


      subroutine pressure_inlet(face)
        implicit none
        character(len=*), intent(in) :: face
        call fix(density , fixed_density , face)
        call fix(x_speed , fixed_x_speed , face)
        call fix(y_speed , fixed_y_speed , face)
        call fix(z_speed , fixed_z_speed , face)
        call copy(pressure, "flat", face)
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
      end subroutine pressure_inlet


      subroutine pressure_outlet(face)
        implicit none
        character(len=*), intent(in) :: face
        call copy(density, "flat", face)
        call copy(x_speed, "flat", face)
        call copy(y_speed, "flat", face)
        call copy(z_speed, "flat", face)
        call fix(pressure, fixed_pressure, face)
        select case (turbulence)
          case('none')
            !do nothing
            continue
          case('sst')
            call copy(tk, "flat", face)
            call copy(tw, "flat", face)
          case DEFAULT
            call turbulence_read_error()
        end select
      end subroutine pressure_outlet

      subroutine adiabatic_wall(face)
        implicit none
        character(len=*), intent(in) :: face
        call copy(density , "symm",  face)
        call copy(pressure, "symm",  face)
        call no_slip(face)
      end subroutine adiabatic_wall


      subroutine slip_wall(face)
        implicit none
        character(len=*), intent(in) :: face
        call copy(density , "symm", face)
        call copy(pressure, "symm", face)
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
            print*, "ERROR: wrong face for boundary condition"
        end select
            
      end subroutine fix


      subroutine copy(var, type, face)
        implicit none
        character(len=*), intent(in) :: face
        character(len=*), intent(in) :: type
        real, dimension(-2:imx+2, -2:jmx+2, -2:kmx+2), intent(inout) :: var

        real :: a1=1
        real :: a2=1
        real :: a3=0

        integer :: i1=1
        integer :: i2=2
        integer :: i3=3

        select case(type)
          case("anti")
            a1 =  1.
            a2 = -1.
            a3 =  0.
          case("flat")
            i1 =  1
            i2 =  1
            i3 =  1
          case("symm")
            a1 =c1
            a2 =c2
            a3 =c3
            ! do nothing
            ! use default value
            continue
          case DEFAULT
            print*, "ERROR: Wrong boundary condition type"
        end select

        select case(face)
          case("imin")
              var(      0, 1:jmx-1, 1:kmx-1) = (a2*var(     i1, 1:jmx-1, 1:kmx-1)-a3*var(    i1+1, 1:jmx-1, 1:kmx-1))/a1
              var(     -1, 1:jmx-1, 1:kmx-1) = (a2*var(     i2, 1:jmx-1, 1:kmx-1)-a3*var(    i2+1, 1:jmx-1, 1:kmx-1))/a1
              var(     -2, 1:jmx-1, 1:kmx-1) = (a2*var(     i3, 1:jmx-1, 1:kmx-1)-a3*var(    i3+1, 1:jmx-1, 1:kmx-1))/a1
          case("imax")
              var(  imx  , 1:jmx-1, 1:kmx-1) = (a2*var( imx-i1, 1:jmx-1, 1:kmx-1)-a3*var(imx-i1-1, 1:jmx-1, 1:kmx-1))/a1
              var(  imx+1, 1:jmx-1, 1:kmx-1) = (a2*var( imx-i2, 1:jmx-1, 1:kmx-1)-a3*var(imx-i2-1, 1:jmx-1, 1:kmx-1))/a1
              var(  imx+2, 1:jmx-1, 1:kmx-1) = (a2*var( imx-i3, 1:jmx-1, 1:kmx-1)-a3*var(imx-i3-1, 1:jmx-1, 1:kmx-1))/a1
          case("jmin")
              var(1:imx-1,       0, 1:kmx-1) = (a2*var(1:imx-1,      i1, 1:kmx-1)-a3*var(1:imx-1 ,    i1+1, 1:kmx-1))/a1
              var(1:imx-1,      -1, 1:kmx-1) = (a2*var(1:imx-1,      i2, 1:kmx-1)-a3*var(1:imx-1 ,    i2+1, 1:kmx-1))/a1
              var(1:imx-1,      -2, 1:kmx-1) = (a2*var(1:imx-1,      i3, 1:kmx-1)-a3*var(1:imx-1 ,    i3+1, 1:kmx-1))/a1
              print*, var(1, i1, 1), var(1,  0, 1)
              print*, var(1, i2, 1), var(1, -1, 1)
              print*, var(1, i3, 1), var(1, -2, 1)
          case("jmax")
              var(1:imx-1,   jmx  , 1:kmx-1) = (a2*var(1:imx-1,  jmx-i1, 1:kmx-1)-a3*var(1:imx-1 ,jmx-i1-1, 1:kmx-1))/a1
              var(1:imx-1,   jmx+1, 1:kmx-1) = (a2*var(1:imx-1,  jmx-i2, 1:kmx-1)-a3*var(1:imx-1 ,jmx-i2-1, 1:kmx-1))/a1
              var(1:imx-1,   jmx+2, 1:kmx-1) = (a2*var(1:imx-1,  jmx-i3, 1:kmx-1)-a3*var(1:imx-1 ,jmx-i3-1, 1:kmx-1))/a1
          case("kmin")
              var(1:imx-1, 1:jmx-1,       0) = (a2*var(1:imx-1, 1:jmx-1,      i1)-a3*var(1:imx-1 , 1:jmx-1,    i1+1))/a1
              var(1:imx-1, 1:jmx-1,      -1) = (a2*var(1:imx-1, 1:jmx-1,      i2)-a3*var(1:imx-1 , 1:jmx-1,    i2+1))/a1
              var(1:imx-1, 1:jmx-1,      -2) = (a2*var(1:imx-1, 1:jmx-1,      i3)-a3*var(1:imx-1 , 1:jmx-1,    i3+1))/a1
          case("kmax")
              var(1:imx-1, 1:jmx-1,   kmx  ) = (a2*var(1:imx-1, 1:jmx-1,  kmx-i1)-a3*var(1:imx-1 , 1:jmx-1,kmx-i1-1))/a1
              var(1:imx-1, 1:jmx-1,   kmx+1) = (a2*var(1:imx-1, 1:jmx-1,  kmx-i2)-a3*var(1:imx-1 , 1:jmx-1,kmx-i2-1))/a1
              var(1:imx-1, 1:jmx-1,   kmx+2) = (a2*var(1:imx-1, 1:jmx-1,  kmx-i3)-a3*var(1:imx-1 , 1:jmx-1,kmx-i3-1))/a1
          case DEFAULT
            print*, "ERROR: wrong face for boundary condition"
        end select
      end subroutine copy


      subroutine flow_tangency(face)
        implicit none
        character(len=*), intent(in) :: face
        real, dimension(:,:), allocatable :: dot


        select case(face)
          case("imin")
            allocate(dot(1:jmx-1, 1:kmx-1))
            call find_dot(dot,1,face)
            x_speed(      0, 1:jmx-1, 1:kmx-1) = x_speed(     1, 1:jmx-1, 1:kmx-1) - (2*dot*xnx(      1, 1:jmx-1, 1:kmx-1))
            y_speed(      0, 1:jmx-1, 1:kmx-1) = y_speed(     1, 1:jmx-1, 1:kmx-1) - (2*dot*xny(      1, 1:jmx-1, 1:kmx-1))
            z_speed(      0, 1:jmx-1, 1:kmx-1) = z_speed(     1, 1:jmx-1, 1:kmx-1) - (2*dot*xnz(      1, 1:jmx-1, 1:kmx-1))
            call find_dot(dot,2,face)
            x_speed(     -1, 1:jmx-1, 1:kmx-1) = x_speed(     2, 1:jmx-1, 1:kmx-1) - (2*dot*xnx(      1, 1:jmx-1, 1:kmx-1))
            y_speed(     -1, 1:jmx-1, 1:kmx-1) = y_speed(     2, 1:jmx-1, 1:kmx-1) - (2*dot*xny(      1, 1:jmx-1, 1:kmx-1))
            z_speed(     -1, 1:jmx-1, 1:kmx-1) = z_speed(     2, 1:jmx-1, 1:kmx-1) - (2*dot*xnz(      1, 1:jmx-1, 1:kmx-1))
            call find_dot(dot,3,face)
            x_speed(     -2, 1:jmx-1, 1:kmx-1) = x_speed(     3, 1:jmx-1, 1:kmx-1) - (2*dot*xnx(      1, 1:jmx-1, 1:kmx-1))
            y_speed(     -2, 1:jmx-1, 1:kmx-1) = y_speed(     3, 1:jmx-1, 1:kmx-1) - (2*dot*xny(      1, 1:jmx-1, 1:kmx-1))
            z_speed(     -2, 1:jmx-1, 1:kmx-1) = z_speed(     3, 1:jmx-1, 1:kmx-1) - (2*dot*xnz(      1, 1:jmx-1, 1:kmx-1))
            deallocate(dot)
          case("imax")
            allocate(dot(1:jmx-1, 1:kmx-1))
            call find_dot(dot,1,face)
            x_speed(  imx+0, 1:jmx-1, 1:kmx-1) = x_speed( imx-1, 1:jmx-1, 1:kmx-1) - (2*dot*xnx(    imx, 1:jmx-1, 1:kmx-1))
            y_speed(  imx+0, 1:jmx-1, 1:kmx-1) = y_speed( imx-1, 1:jmx-1, 1:kmx-1) - (2*dot*xny(    imx, 1:jmx-1, 1:kmx-1))
            z_speed(  imx+0, 1:jmx-1, 1:kmx-1) = z_speed( imx-1, 1:jmx-1, 1:kmx-1) - (2*dot*xnz(    imx, 1:jmx-1, 1:kmx-1))
            call find_dot(dot,2,face)
            x_speed(  imx+1, 1:jmx-1, 1:kmx-1) = x_speed( imx-2, 1:jmx-1, 1:kmx-1) - (2*dot*xnx(    imx, 1:jmx-1, 1:kmx-1))
            y_speed(  imx+1, 1:jmx-1, 1:kmx-1) = y_speed( imx-2, 1:jmx-1, 1:kmx-1) - (2*dot*xny(    imx, 1:jmx-1, 1:kmx-1))
            z_speed(  imx+1, 1:jmx-1, 1:kmx-1) = z_speed( imx-2, 1:jmx-1, 1:kmx-1) - (2*dot*xnz(    imx, 1:jmx-1, 1:kmx-1))
            call find_dot(dot,3,face)
            x_speed(  imx+2, 1:jmx-1, 1:kmx-1) = x_speed( imx-3, 1:jmx-1, 1:kmx-1) - (2*dot*xnx(    imx, 1:jmx-1, 1:kmx-1))
            y_speed(  imx+2, 1:jmx-1, 1:kmx-1) = y_speed( imx-3, 1:jmx-1, 1:kmx-1) - (2*dot*xny(    imx, 1:jmx-1, 1:kmx-1))
            z_speed(  imx+2, 1:jmx-1, 1:kmx-1) = z_speed( imx-3, 1:jmx-1, 1:kmx-1) - (2*dot*xnz(    imx, 1:jmx-1, 1:kmx-1))
            deallocate(dot)
          case ("jmin")
            allocate(dot(1:imx-1, 1:kmx-1))
            call find_dot(dot,1,face)
            x_speed(1:imx-1,       0, 1:kmx-1) = x_speed(1:imx-1,      1, 1:kmx-1) - (2*dot*ynx(1:imx-1,       1, 1:kmx-1))
            y_speed(1:imx-1,       0, 1:kmx-1) = y_speed(1:imx-1,      1, 1:kmx-1) - (2*dot*yny(1:imx-1,       1, 1:kmx-1))
            z_speed(1:imx-1,       0, 1:kmx-1) = z_speed(1:imx-1,      1, 1:kmx-1) - (2*dot*ynz(1:imx-1,       1, 1:kmx-1))
            call find_dot(dot,2,face)
            x_speed(1:imx-1,      -1, 1:kmx-1) = x_speed(1:imx-1,      2, 1:kmx-1) - (2*dot*ynx(1:imx-1,       1, 1:kmx-1))
            y_speed(1:imx-1,      -1, 1:kmx-1) = y_speed(1:imx-1,      2, 1:kmx-1) - (2*dot*yny(1:imx-1,       1, 1:kmx-1))
            z_speed(1:imx-1,      -1, 1:kmx-1) = z_speed(1:imx-1,      2, 1:kmx-1) - (2*dot*ynz(1:imx-1,       1, 1:kmx-1))
            call find_dot(dot,3,face)
            x_speed(1:imx-1,      -2, 1:kmx-1) = x_speed(1:imx-1,      3, 1:kmx-1) - (2*dot*ynx(1:imx-1,       1, 1:kmx-1))
            y_speed(1:imx-1,      -2, 1:kmx-1) = y_speed(1:imx-1,      3, 1:kmx-1) - (2*dot*yny(1:imx-1,       1, 1:kmx-1))
            z_speed(1:imx-1,      -2, 1:kmx-1) = z_speed(1:imx-1,      3, 1:kmx-1) - (2*dot*ynz(1:imx-1,       1, 1:kmx-1))
            deallocate(dot)
          case ("jmax")
            allocate(dot(1:imx-1, 1:kmx-1))
            call find_dot(dot,1,face)
            x_speed(1:imx-1,     jmx, 1:kmx-1) = x_speed(1:imx-1,  jmx-1, 1:kmx-1) - (2*dot*ynx(1:imx-1,     jmx, 1:kmx-1))
            y_speed(1:imx-1,     jmx, 1:kmx-1) = y_speed(1:imx-1,  jmx-1, 1:kmx-1) - (2*dot*yny(1:imx-1,     jmx, 1:kmx-1))
            z_speed(1:imx-1,     jmx, 1:kmx-1) = z_speed(1:imx-1,  jmx-1, 1:kmx-1) - (2*dot*ynz(1:imx-1,     jmx, 1:kmx-1))
            call find_dot(dot,2,face)
            x_speed(1:imx-1,   jmx+1, 1:kmx-1) = x_speed(1:imx-1,  jmx-2, 1:kmx-1) - (2*dot*ynx(1:imx-1,     jmx, 1:kmx-1))
            y_speed(1:imx-1,   jmx+1, 1:kmx-1) = y_speed(1:imx-1,  jmx-2, 1:kmx-1) - (2*dot*yny(1:imx-1,     jmx, 1:kmx-1))
            z_speed(1:imx-1,   jmx+1, 1:kmx-1) = z_speed(1:imx-1,  jmx-2, 1:kmx-1) - (2*dot*ynz(1:imx-1,     jmx, 1:kmx-1))
            call find_dot(dot,3,face)
            x_speed(1:imx-1,   jmx+2, 1:kmx-1) = x_speed(1:imx-1,  jmx-3, 1:kmx-1) - (2*dot*ynx(1:imx-1,     jmx, 1:kmx-1))
            y_speed(1:imx-1,   jmx+2, 1:kmx-1) = y_speed(1:imx-1,  jmx-3, 1:kmx-1) - (2*dot*yny(1:imx-1,     jmx, 1:kmx-1))
            z_speed(1:imx-1,   jmx+2, 1:kmx-1) = z_speed(1:imx-1,  jmx-3, 1:kmx-1) - (2*dot*ynz(1:imx-1,     jmx, 1:kmx-1))
            deallocate(dot)
          case("kmin")
            allocate(dot(1:imx-1, 1:jmx-1))
            call find_dot(dot,1,face)
            x_speed(1:imx-1, 1:jmx-1,       0) = x_speed(1:imx-1, 1:jmx-1,      1) - (2*dot*znx(1:imx-1, 1:jmx-1,       1))
            y_speed(1:imx-1, 1:jmx-1,       0) = y_speed(1:imx-1, 1:jmx-1,      1) - (2*dot*zny(1:imx-1, 1:jmx-1,       1))
            z_speed(1:imx-1, 1:jmx-1,       0) = z_speed(1:imx-1, 1:jmx-1,      1) - (2*dot*znz(1:imx-1, 1:jmx-1,       1))
            call find_dot(dot,2,face)
            x_speed(1:imx-1, 1:jmx-1,      -1) = x_speed(1:imx-1, 1:jmx-1,      2) - (2*dot*znx(1:imx-1, 1:jmx-1,       1))
            y_speed(1:imx-1, 1:jmx-1,      -1) = y_speed(1:imx-1, 1:jmx-1,      2) - (2*dot*zny(1:imx-1, 1:jmx-1,       1))
            z_speed(1:imx-1, 1:jmx-1,      -1) = z_speed(1:imx-1, 1:jmx-1,      2) - (2*dot*znz(1:imx-1, 1:jmx-1,       1))
            call find_dot(dot,3,face)
            x_speed(1:imx-1, 1:jmx-1,      -2) = x_speed(1:imx-1, 1:jmx-1,      3) - (2*dot*znx(1:imx-1, 1:jmx-1,       1))
            y_speed(1:imx-1, 1:jmx-1,      -2) = y_speed(1:imx-1, 1:jmx-1,      3) - (2*dot*zny(1:imx-1, 1:jmx-1,       1))
            z_speed(1:imx-1, 1:jmx-1,      -2) = z_speed(1:imx-1, 1:jmx-1,      3) - (2*dot*znz(1:imx-1, 1:jmx-1,       1))
            deallocate(dot)
          case("kmax")
            allocate(dot(1:imx-1, 1:jmx-1))
            call find_dot(dot,1,face)
            x_speed(1:imx-1, 1:jmx-1,   kmx  ) = x_speed(1:imx-1, 1:jmx-1,  kmx-1) - (2*dot*znx(1:imx-1, 1:jmx-1,     kmx))
            y_speed(1:imx-1, 1:jmx-1,   kmx  ) = y_speed(1:imx-1, 1:jmx-1,  kmx-1) - (2*dot*zny(1:imx-1, 1:jmx-1,     kmx))
            z_speed(1:imx-1, 1:jmx-1,   kmx  ) = z_speed(1:imx-1, 1:jmx-1,  kmx-1) - (2*dot*znz(1:imx-1, 1:jmx-1,     kmx))
            call find_dot(dot,2,face)
            x_speed(1:imx-1, 1:jmx-1,   kmx+1) = x_speed(1:imx-1, 1:jmx-1,  kmx-2) - (2*dot*znx(1:imx-1, 1:jmx-1,     kmx))
            y_speed(1:imx-1, 1:jmx-1,   kmx+1) = y_speed(1:imx-1, 1:jmx-1,  kmx-2) - (2*dot*zny(1:imx-1, 1:jmx-1,     kmx))
            z_speed(1:imx-1, 1:jmx-1,   kmx+1) = z_speed(1:imx-1, 1:jmx-1,  kmx-2) - (2*dot*znz(1:imx-1, 1:jmx-1,     kmx))
            call find_dot(dot,3,face)
            x_speed(1:imx-1, 1:jmx-1,   kmx+2) = x_speed(1:imx-1, 1:jmx-1,  kmx-3) - (2*dot*znx(1:imx-1, 1:jmx-1,     kmx))
            y_speed(1:imx-1, 1:jmx-1,   kmx+2) = y_speed(1:imx-1, 1:jmx-1,  kmx-3) - (2*dot*zny(1:imx-1, 1:jmx-1,     kmx))
            z_speed(1:imx-1, 1:jmx-1,   kmx+2) = z_speed(1:imx-1, 1:jmx-1,  kmx-3) - (2*dot*znz(1:imx-1, 1:jmx-1,     kmx))
            deallocate(dot)
        end select
        select case(turbulence)
          case("none")
            !do nothing
            continue
          case("sst")
            call copy(tk  , "symm", face)
            call copy(tw  , "symm", face)
          case DEFAULT
            call turbulence_read_error()
        end select
      end subroutine flow_tangency
      
      subroutine find_dot(dot,i, face)
        implicit none
        integer, intent(in) :: i
        character(len=*), intent(in) :: face
        real, dimension(:,:), intent(out) :: dot

        select case(face)
          case("imin")
            dot = (x_speed(i, 1:jmx-1, 1:kmx-1) * xnx(1, 1:jmx-1, 1:kmx-1)) + &
                  (y_speed(i, 1:jmx-1, 1:kmx-1) * xny(1, 1:jmx-1, 1:kmx-1)) + &
                  (z_speed(i, 1:jmx-1, 1:kmx-1) * xnz(1, 1:jmx-1, 1:kmx-1))
          case("imax")
            dot = (x_speed(imx-i, 1:jmx-1, 1:kmx-1) * xnx(imx, 1:jmx-1, 1:kmx-1)) + &
                  (y_speed(imx-i, 1:jmx-1, 1:kmx-1) * xny(imx, 1:jmx-1, 1:kmx-1)) + &
                  (z_speed(imx-i, 1:jmx-1, 1:kmx-1) * xnz(imx, 1:jmx-1, 1:kmx-1))
          case ("jmin")
            dot = (x_speed(1:imx-1, i, 1:kmx-1) * ynx(1:imx-1, 1, 1:kmx-1)) + &
                  (y_speed(1:imx-1, i, 1:kmx-1) * yny(1:imx-1, 1, 1:kmx-1)) + &
                  (z_speed(1:imx-1, i, 1:kmx-1) * ynz(1:imx-1, 1, 1:kmx-1))
          case ("jmax")
            dot = (x_speed(1:imx-1, jmx-i, 1:kmx-1) * ynx(1:imx-1, jmx, 1:kmx-1)) + &
                  (y_speed(1:imx-1, jmx-i, 1:kmx-1) * yny(1:imx-1, jmx, 1:kmx-1)) + &
                  (z_speed(1:imx-1, jmx-i, 1:kmx-1) * ynz(1:imx-1, jmx, 1:kmx-1))
          case("kmin")
            dot = (x_speed(1:imx-1, 1:jmx-1, i) * znx(1:imx-1, 1:jmx-1, 1)) + &
                  (y_speed(1:imx-1, 1:jmx-1, i) * zny(1:imx-1, 1:jmx-1, 1)) + &
                  (z_speed(1:imx-1, 1:jmx-1, i) * znz(1:imx-1, 1:jmx-1, 1))
          case("kmax")
            dot = (x_speed(1:imx-1, 1:jmx-1, kmx-i) * znx(1:imx-1, 1:jmx-1, kmx)) + &
                  (y_speed(1:imx-1, 1:jmx-1, kmx-i) * zny(1:imx-1, 1:jmx-1, kmx)) + &
                  (z_speed(1:imx-1, 1:jmx-1, kmx-i) * znz(1:imx-1, 1:jmx-1, kmx))
        end select
      end subroutine find_dot

      subroutine no_slip(face)
        implicit none
        character(len=*), intent(in) :: face
        call copy(x_speed, "anti", face)
        call copy(y_speed, "anti", face)
        call copy(z_speed, "anti", face)
        select case(turbulence)
          case("none")
            !do nothing
            continue
          case("sst")
            call copy(tk  , "anti", face)
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

      
    subroutine set_wall_flux_zero()
      implicit none
      integer :: i
      character(len=4) :: face
      
      do i = 1,6

        face = face_names(i)

        select case(id(i))

          case(-5)
            !call making_flux_zero(face)

          case(-6)
            !call making_flux_zero(face)

        end select
      end do

    end subroutine set_wall_flux_zero

!    subroutine making_flux_zero(face)
!      implicit none
!      character(len=*), intent(in) :: face
!
!      select case(face)
!        case("imin")
!          F(1,:,:)   = 0.
!        case("imax")
!          F(imx,:,:) = 0.
!        case("jmin")
!          G(:,1,:)   = 0.
!        case("jmax")
!          G(:,jmx,:) = 0.
!        case("kmin")
!          H(:,:,1)   = 0.
!        case("kmax")
!          H(:,:,kmx) = 0.
!      end select
!
!    end subroutine making_flux_zero

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
