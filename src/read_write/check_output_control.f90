  !< Check for the input from the output_control.md file
module check_output_control
  !< Check for the input from the output_control.md file
  use vartypes
  implicit none
  private
  public :: verify_write_control
  public :: verify_read_control


  contains

    subroutine verify_write_control(control, scheme, flow)
      !< Verify all the variable being asked to write in the output file. 
      !< This is a fail-safe subroutine which do not allow to write the incorrect input variable
      implicit none
      type(controltype), intent(inout) :: control
      type(schemetype) , intent(in) :: scheme
      type(flowtype)   , intent(in) :: flow
      integer :: n
      character(len=*), parameter :: err="Control Error: can't write variable - "

      do n = 1,control%w_count

        select case (trim(lcase(control%w_list(n))))
        
          case('velocity','vel','speed','u','v')
            control%w_list(n) = "Velocity"

          case('density','rho')
            control%w_list(n) = "Density"
          
          case('pressure','presssure','p')
            control%w_list(n) = "Pressure"

          case('mu','viscosity','mu_l','laminar_viscosity','muv','mu_v')
            if (flow%mu_ref/=0.0) then
              control%w_list(n) = "Mu"
            else
              print*, err//trim(control%w_list(n))//" to file"
              control%w_list(n) = "do not write"
            end if
            
          case('mu_t','turbulent_viscosity','mut')
            if (scheme%turbulence/='none') then
              control%w_list(n) = "Mu_t"
            else
              print*, err//trim(control%w_list(n))//" to file"
              control%w_list(n) = "do not write"
            end if
            
          case('tke','tk','turbulent_kinetic_enrgy','k')
            select case (trim(scheme%turbulence))
              case('sst', 'sst2003','kw','bsl','kkl','ke','des-sst')
                control%w_list(n) = "TKE"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
             end select

          case('omega','tw')
            select case (trim(scheme%turbulence))
              case('sst', 'sst2003','kw','bsl','des-sst')
                control%w_list(n) = "Omega"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
             end select

          case('dissipation','te','teps','eps')
            select case (trim(scheme%turbulence))
              case('ke')
                control%w_list(n) = "Dissipation"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
             end select

          case('kl')
            select case (trim(scheme%turbulence))
              case('kkl')
                control%w_list(n) = "Kl"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
             end select

          case('tv')
            select case (trim(scheme%turbulence))
              case('sa', 'saBC')
                control%w_list(n) = "tv"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
             end select

          case('wall_distance', 'dist', 'wall_dist', 'wdist')
            if(scheme%turbulence/="none") then
              control%w_list(n) = "Wall_distance"
            else
              print*, err//trim(control%w_list(n))//" to file"
              control%w_list(n) = "do not write"
            end if

          case('resnorm')
            control%w_list(n) = "Resnorm"

          case('tke_residue')
            control%w_list(n) = "TKE_residue"

          case('omega_residue')
            control%w_list(n) = "Omega_residue"

          case('tv_residue')
            control%w_list(n) = "Tv_residue"

          case('mass_residue')
            control%w_list(n) = "Mass_residue"

          case('x_mom_residue')
            control%w_list(n) = "X_mom_residue"

          case('y_mom_residue')
            control%w_list(n) = "Y_mom_residue"

          case('z_mom_residue')
            control%w_list(n) = "Z_mom_residue"

          case('energy_residue')
            control%w_list(n) = "Energy_residue"

          case('f1')
            if(trim(scheme%turbulence)=='sst' .or. trim(scheme%turbulence)=='sst2003')then
              control%w_list(n) = "F1"
            else
              control%w_list(n) = 'do not write'
            end if

          case('dudx')
            control%w_list(n) = "Dudx"

          case('dudy')
            control%w_list(n) = "Dudy"

          case('dudz')
            control%w_list(n) = "Dudz"

          case('dvdx')
            control%w_list(n) = "Dvdx"

          case('dvdy')
            control%w_list(n) = "Dvdy"

          case('dvdz')
            control%w_list(n) = "Dvdz"

          case('dwdx')
            control%w_list(n) = "Dwdx"

          case('dwdy')
            control%w_list(n) = "Dwdy"

          case('dwdz')
            control%w_list(n) = "Dwdz"

          case('dTdx')
            control%w_list(n) = "DTdx"

          case('dTdy')
            control%w_list(n) = "DTdy"

          case('dTdz')
            control%w_list(n) = "DTdz"

          case('dtkdx')
            select case (trim(scheme%turbulence))
              case('sst', 'sst2003','kw','bsl','kkl','ke','des-sst')
                control%w_list(n) = "Dtkdx"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
            end select

          case('dtkdy')
            select case (trim(scheme%turbulence))
              case('sst', 'sst2003','kw','bsl','kkl','ke','des-sst')
                control%w_list(n) = "Dtkdy"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
            end select

          case('dtkdz')
            select case (trim(scheme%turbulence))
              case('sst','sst2003','kw','bsl','kkl','ke','des-sst')
                control%w_list(n) = "Dtkdz"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
            end select

          case('dtwdx')
            select case (trim(scheme%turbulence))
              case('sst','sst2003','kw','bsl','des-sst')
                control%w_list(n) = "Dtwdx"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
            end select

          case('dtwdy')
            select case (trim(scheme%turbulence))
              case('sst','sst2003','kw','bsl','des-sst')
                control%w_list(n) = "Dtwdy"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
            end select

          case('dtwdz')
            select case (trim(scheme%turbulence))
              case('sst','sst2003','kw','bsl','des-sst')
                control%w_list(n) = "Dtwdz"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
            end select

          case('dtedx')
            select case (trim(scheme%turbulence))
              case('ke')
                control%w_list(n) = "Dtedx"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
            end select

          case('dtedy')
            select case (trim(scheme%turbulence))
              case('ke')
                control%w_list(n) = "Dtedy"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
            end select

          case('dtedz')
            select case (trim(scheme%turbulence))
              case('ke')
                control%w_list(n) = "Dtedz"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
            end select

          case('dtkldx')
            select case (trim(scheme%turbulence))
              case('kkl')
                control%w_list(n) = "Dtkldx"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
             end select

          case('dtkldy')
            select case (trim(scheme%turbulence))
              case('kkl')
                control%w_list(n) = "Dtkldy"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
             end select

          case('dtkldz')
            select case (trim(scheme%turbulence))
              case('kkl')
                control%w_list(n) = "Dtkldz"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
             end select

          case('dtvdx')
            select case (trim(scheme%turbulence))
              case('sa', 'saBC')
                control%w_list(n) = "Dtvdx"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
             end select

          case('dtvdy')
            select case (trim(scheme%turbulence))
              case('sa', 'saBC')
                control%w_list(n) = "Dtvdy"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
             end select

          case('dtvdz')
            select case (trim(scheme%turbulence))
              case('sa', 'saBC')
                control%w_list(n) = "Dtvdz"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
             end select

          case('intermittency')
            select case (trim(scheme%turbulence))
              case('saBC')
                control%w_list(n) = "Intermittency"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
             end select

           case('tgm')
            select case (trim(scheme%transition))
              case('lctm2015')
                control%w_list(n) = "tgm"
              case DEFAULT
                print*, err//trim(control%w_list(n))//" to file"
                control%w_list(n) = "do not write"
             end select

          case('extravar1','extravar2', 'extravar3', 'extravar4', 'extravar5')
            control%w_list(n) = trim(lcase(control%w_list(n)))

          case Default
            print*, err//trim(control%w_list(n))//" to file"
            control%w_list(n) = "do not write"

        end select
      end do

    end subroutine verify_write_control

    subroutine verify_read_control(control, scheme)
      !< Verify all the variable being asked to read in the output file. 
      !< This is a fail-safe subroutine which do not allow to read the incorrect input variable. 
      !< Based on previous flow type some varible might be skipped
      implicit none
      type(controltype), intent(inout) :: control
      type(schemetype) , intent(in) :: scheme
      integer :: n
      character(len=*), parameter :: err="Control Error: can't read variable - "

      do n = 1,control%r_count

        select case (trim(lcase(control%r_list(n))))
        
          case('velocity','vel','speed','u','v')
            control%r_list(n) = "Velocity"

          case('density','rho')
            control%r_list(n) = "Density"
          
          case('pressure','presssure','p')
            control%r_list(n) = "Pressure"

          case('mu','viscosity','mu_l','laminar_viscosity','muv','mu_v')
            control%r_list(n) = "do not read"
           ! if (flow%mu_ref/=0.0) then
           !   control%r_list(n) = "Mu"
           ! else
           !   print*, err//trim(control%r_list(n))//" from file"
           !   control%r_list(n) = "do not read"
           ! end if
            
          case('mu_t','turbulent_viscosity','mut')
            control%r_list(n) = "do not read"
            !if (scheme%turbulence/='none') then
            !  control%r_list(n) = "Mu_t"
            !else
            !  print*, err//trim(control%r_list(n))//" from file"
            !  control%r_list(n) = "do not read"
            !end if
            
          case('tke','tk','turbulent_kinetic_enrgy','k')
            select case (trim(scheme%turbulence))
              case('sst','sst2003','kw','bsl','kkl','ke','des-sst')
                select case (trim(control%previous_flow_type))
                  case('sst','sst2003','kw','bsl','kkl','ke','des-sst')
                    control%r_list(n) = "TKE"
                end select
              case DEFAULT
                print*, err//trim(control%w_list(n))//" from file"
                control%r_list(n) = "do not read"
            end select

          case('omega','tw')
            select case (trim(scheme%turbulence))
              case('sst','sst2003','kw','bsl','des-sst')
                select case (trim(control%previous_flow_type))
                  case('sst','sst2003','kw','bsl','des-sst')
                    control%r_list(n) = "Omega"
                  case DEFAULT
                    print*, err//trim(control%w_list(n))//" from file"
                    control%r_list(n) = "do not read"
                end select
              case DEFAULT
                print*, err//trim(control%w_list(n))//" from file"
                control%r_list(n) = "do not read"
            end select

          case('dissipation','te','teps','eps')
            select case (trim(scheme%turbulence))
              case('ke')
                select case (trim(control%previous_flow_type))
                  case('ke')
                    control%r_list(n) = "Dissipation"
                  case DEFAULT
                    print*, err//trim(control%w_list(n))//" to file"
                    control%r_list(n) = "do not write"
                end select
              case DEFAULT
                print*, err//trim(control%w_list(n))//" from file"
                control%r_list(n) = "do not read"
            end select

          case('kl')
            select case (trim(scheme%turbulence))
              case('kkl')
                select case (trim(control%previous_flow_type))
                  case('kkl')
                    control%r_list(n) = "Kl"
                  case DEFAULT
                    print*, err//trim(control%w_list(n))//" to file"
                    control%r_list(n) = "do not write"
                end select
              case DEFAULT
                print*, err//trim(control%w_list(n))//" from file"
                control%r_list(n) = "do not read"
            end select

          case('tv')
            select case (trim(scheme%turbulence))
              case('sa', 'saBC')
                select case (trim(control%previous_flow_type))
                  case('sa', 'saBC')
                    control%r_list(n) = "tv"
                  case DEFAULT
                    print*, err//trim(control%w_list(n))//" to file"
                    control%r_list(n) = "do not write"
                end select
              case DEFAULT
                print*, err//trim(control%w_list(n))//" from file"
                control%r_list(n) = "do not read"
            end select

          case('wall_distance', 'dist', 'wall_dist', 'wdist')
            control%r_list(n) = "do not read"
            !if(scheme%turbulence/="none") then
            !  control%r_list(n) = "Wall_distance"
            !else
            !  print*, err//trim(control%r_list(n))//" from file"
            !  control%r_list(n) = "do not read"
            !end if

          case('intermittency')
            select case (trim(scheme%turbulence))
              case('saBC')
                select case(trim(control%previous_flow_type))
                  case('saBC')
                    control%r_list(n) = "Intermittency"
                  case DEFAULT
                    print*, err//trim(control%r_list(n))//" to file"
                    control%r_list(n) = "do not read"
                end select
              case DEFAULT
                print*, err//trim(control%r_list(n))//" to file"
                control%r_list(n) = "do not read"
             end select

          case('tgm')
            select case (trim(scheme%transition))
              case('lctm2015')
                select case(trim(control%previous_flow_type))
                  case('sst','sst2003')
                    control%r_list(n) = "tgm"
                  case DEFAULT
                    print*, err//trim(control%r_list(n))//" to file"
                    control%r_list(n) = "do not read"
                end select
              case DEFAULT
                print*, err//trim(control%r_list(n))//" to file"
                control%r_list(n) = "do not read"
             end select

          case('resnorm')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Resnorm"

          case('tke_residue')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "TKE_residue"

          case('omega_residue')
            control%r_list(n) = "do not read"
            !control%w_list(n) = "Omega_residue"

          case('f1')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "F1"

          case('dudx')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dudx"

          case('dudy')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dudy"

          case('dudz')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dudz"

          case('dvdx')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dvdx"

          case('dvdy')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dvdy"

          case('dvdz')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dvdz"

          case('dwdx')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dwdx"

          case('dwdy')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dwdy"

          case('dwdz')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dwdz"

          case('dTdx')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "DTdx"

          case('dTdy')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "DTdy"

          case('dTdz')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "DTdz"

          case('dtkdx')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dtkdx"

          case('dtkdy')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dtkdy"

          case('dtkdz')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dtkdz"

          case('dtwdx')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dtwdx"

          case('dtwdy')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dtwdy"

          case('dtwdz')
            control%r_list(n) = "do not read"
            !control%r_list(n) = "Dtwdz"

          case('extravar1','extravar2', 'extravar3', 'extravar4', 'extravar5')
            control%r_list(n) = "do not read"
            !control%r_list(n) = trim(lcase(control%w_list(n)))

          case Default
            print*, err//trim(control%r_list(n))//" from file"
            control%r_list(n) = "do not read"

        end select
      end do

    end subroutine verify_read_control

    function lcase(text) result(res)
      !< Make the whole string to lower case
      CHARACTER(len=*), intent(in)         :: text
      !< Input string of any case
      character(len=STRING_BUFFER_LENGTH) :: res
      !< Output string of lower case
      integer ::  I,C
  
      res=text
      DO I = 1,LEN(TEXT)
        C = INDEX("ABCDEFGHIJKLMNOPQRSTUVWXYZ",TEXT(I:I))
        IF (C.GT.0) res(I:I) = "abcdefghijklmnopqrstuvwxyz"(C:C)
      END DO
  
    end function lcase

end module check_output_control
