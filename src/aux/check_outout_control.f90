  !< Check for the input from the output_control.md file
module check_output_control
  !< Check for the input from the output_control.md file
! ---------------------------------------------
! 170619 -jatinder pal singh sandhu
! Aim: to check wheter input are correct or not
!----------------------------------------------
  use global_vars, only: r_list
  use global_vars, only: w_list
  use global_vars, only: previous_flow_type
  use global_vars, only: turbulence
  use global_vars, only: transition
  use global_vars, only: mu_ref
  use global_vars, only: r_count
  use global_vars, only: w_count
  use str_case   , only: lcase
  implicit none
  private
  public :: verify_write_control
  public :: verify_read_control

  contains

    subroutine verify_write_control()
      !< Verify all the variable being asked to write in the output file. 
      !< This is a fail-safe subroutine which do not allow to write the incorrect input variable
      implicit none
      integer :: n
      character(len=*), parameter :: err="Control Error: can't write variable - "

      do n = 1,w_count

        select case (trim(lcase(w_list(n))))
        
          case('velocity','vel','speed','u','v')
            w_list(n) = "Velocity"

          case('density','rho')
            w_list(n) = "Density"
          
          case('pressure','presssure','p')
            w_list(n) = "Pressure"

          case('mu','viscosity','mu_l','laminar_viscosity','muv','mu_v')
            if (mu_ref/=0.0) then
              w_list(n) = "Mu"
            else
              print*, err//trim(w_list(n))//" to file"
              w_list(n) = "do not write"
            end if
            
          case('mu_t','turbulent_viscosity','mut')
            if (turbulence/='none') then
              w_list(n) = "Mu_t"
            else
              print*, err//trim(w_list(n))//" to file"
              w_list(n) = "do not write"
            end if
            
          case('tke','tk','turbulent_kinetic_enrgy','k')
            select case (trim(turbulence))
              case('sst', 'sst2003','kw','bsl','kkl','ke','des-sst')
                w_list(n) = "TKE"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
             end select

          case('omega','tw')
            select case (trim(turbulence))
              case('sst', 'sst2003','kw','bsl','des-sst')
                w_list(n) = "Omega"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
             end select

          case('dissipation','te','teps','eps')
            select case (trim(turbulence))
              case('ke')
                w_list(n) = "Dissipation"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
             end select

          case('kl')
            select case (trim(turbulence))
              case('kkl')
                w_list(n) = "Kl"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
             end select

          case('tv')
            select case (trim(turbulence))
              case('sa', 'saBC')
                w_list(n) = "tv"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
             end select

          case('wall_distance', 'dist', 'wall_dist', 'wdist')
            if(turbulence/="none") then
              w_list(n) = "Wall_distance"
            else
              print*, err//trim(w_list(n))//" to file"
              w_list(n) = "do not write"
            end if

          case('resnorm')
            w_list(n) = "Resnorm"

          case('tke_residue')
            w_list(n) = "TKE_residue"

          case('omega_residue')
            w_list(n) = "Omega_residue"

          case('tv_residue')
            w_list(n) = "Tv_residue"

          case('mass_residue')
            w_list(n) = "Mass_residue"

          case('x_mom_residue')
            w_list(n) = "X_mom_residue"

          case('y_mom_residue')
            w_list(n) = "Y_mom_residue"

          case('z_mom_residue')
            w_list(n) = "Z_mom_residue"

          case('energy_residue')
            w_list(n) = "Energy_residue"

          case('f1')
            if(trim(turbulence)=='sst' .or. trim(turbulence)=='sst2003')then
              w_list(n) = "F1"
            else
              w_list(n) = 'do not write'
            end if

          case('dudx')
            w_list(n) = "Dudx"

          case('dudy')
            w_list(n) = "Dudy"

          case('dudz')
            w_list(n) = "Dudz"

          case('dvdx')
            w_list(n) = "Dvdx"

          case('dvdy')
            w_list(n) = "Dvdy"

          case('dvdz')
            w_list(n) = "Dvdz"

          case('dwdx')
            w_list(n) = "Dwdx"

          case('dwdy')
            w_list(n) = "Dwdy"

          case('dwdz')
            w_list(n) = "Dwdz"

          case('dTdx')
            w_list(n) = "DTdx"

          case('dTdy')
            w_list(n) = "DTdy"

          case('dTdz')
            w_list(n) = "DTdz"

          case('dtkdx')
            select case (trim(turbulence))
              case('sst', 'sst2003','kw','bsl','kkl','ke','des-sst')
                w_list(n) = "Dtkdx"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
            end select

          case('dtkdy')
            select case (trim(turbulence))
              case('sst', 'sst2003','kw','bsl','kkl','ke','des-sst')
                w_list(n) = "Dtkdy"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
            end select

          case('dtkdz')
            select case (trim(turbulence))
              case('sst','sst2003','kw','bsl','kkl','ke','des-sst')
                w_list(n) = "Dtkdz"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
            end select

          case('dtwdx')
            select case (trim(turbulence))
              case('sst','sst2003','kw','bsl','des-sst')
                w_list(n) = "Dtwdx"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
            end select

          case('dtwdy')
            select case (trim(turbulence))
              case('sst','sst2003','kw','bsl','des-sst')
                w_list(n) = "Dtwdy"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
            end select

          case('dtwdz')
            select case (trim(turbulence))
              case('sst','sst2003','kw','bsl','des-sst')
                w_list(n) = "Dtwdz"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
            end select

          case('dtedx')
            select case (trim(turbulence))
              case('ke')
                w_list(n) = "Dtedx"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
            end select

          case('dtedy')
            select case (trim(turbulence))
              case('ke')
                w_list(n) = "Dtedy"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
            end select

          case('dtedz')
            select case (trim(turbulence))
              case('ke')
                w_list(n) = "Dtedz"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
            end select

          case('dtkldx')
            select case (trim(turbulence))
              case('kkl')
                w_list(n) = "Dtkldx"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
             end select

          case('dtkldy')
            select case (trim(turbulence))
              case('kkl')
                w_list(n) = "Dtkldy"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
             end select

          case('dtkldz')
            select case (trim(turbulence))
              case('kkl')
                w_list(n) = "Dtkldz"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
             end select

          case('dtvdx')
            select case (trim(turbulence))
              case('sa', 'saBC')
                w_list(n) = "Dtvdx"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
             end select

          case('dtvdy')
            select case (trim(turbulence))
              case('sa', 'saBC')
                w_list(n) = "Dtvdy"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
             end select

          case('dtvdz')
            select case (trim(turbulence))
              case('sa', 'saBC')
                w_list(n) = "Dtvdz"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
             end select

          case('intermittency')
            select case (trim(turbulence))
              case('saBC')
                w_list(n) = "Intermittency"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
             end select

           case('tgm')
            select case (trim(transition))
              case('lctm2015')
                w_list(n) = "tgm"
              case DEFAULT
                print*, err//trim(w_list(n))//" to file"
                w_list(n) = "do not write"
             end select

          case('extravar1','extravar2', 'extravar3', 'extravar4', 'extravar5')
            w_list(n) = trim(lcase(w_list(n)))

          case Default
            print*, err//trim(w_list(n))//" to file"
            w_list(n) = "do not write"

        end select
      end do

    end subroutine verify_write_control

    subroutine verify_read_control()
      !< Verify all the variable being asked to read in the output file. 
      !< This is a fail-safe subroutine which do not allow to read the incorrect input variable. 
      !< Based on previous flow type some varible might be skipped
      implicit none
      integer :: n
      character(len=*), parameter :: err="Control Error: can't read variable - "

      do n = 1,r_count

        select case (trim(lcase(r_list(n))))
        
          case('velocity','vel','speed','u','v')
            r_list(n) = "Velocity"

          case('density','rho')
            r_list(n) = "Density"
          
          case('pressure','presssure','p')
            r_list(n) = "Pressure"

          case('mu','viscosity','mu_l','laminar_viscosity','muv','mu_v')
            r_list(n) = "do not read"
           ! if (mu_ref/=0.0) then
           !   r_list(n) = "Mu"
           ! else
           !   print*, err//trim(r_list(n))//" from file"
           !   r_list(n) = "do not read"
           ! end if
            
          case('mu_t','turbulent_viscosity','mut')
            r_list(n) = "do not read"
            !if (turbulence/='none') then
            !  r_list(n) = "Mu_t"
            !else
            !  print*, err//trim(r_list(n))//" from file"
            !  r_list(n) = "do not read"
            !end if
            
          case('tke','tk','turbulent_kinetic_enrgy','k')
            select case (trim(turbulence))
              case('sst','sst2003','kw','bsl','kkl','ke','des-sst')
                select case (trim(previous_flow_type))
                  case('sst','sst2003','kw','bsl','kkl','ke','des-sst')
                    r_list(n) = "TKE"
                end select
              case DEFAULT
                print*, err//trim(w_list(n))//" from file"
                r_list(n) = "do not read"
            end select

          case('omega','tw')
            select case (trim(turbulence))
              case('sst','sst2003','kw','bsl','des-sst')
                select case (trim(previous_flow_type))
                  case('sst','sst2003','kw','bsl','des-sst')
                    r_list(n) = "Omega"
                  case DEFAULT
                    print*, err//trim(w_list(n))//" from file"
                    r_list(n) = "do not read"
                end select
              case DEFAULT
                print*, err//trim(w_list(n))//" from file"
                r_list(n) = "do not read"
            end select

          case('dissipation','te','teps','eps')
            select case (trim(turbulence))
              case('ke')
                select case (trim(previous_flow_type))
                  case('ke')
                    r_list(n) = "Dissipation"
                  case DEFAULT
                    print*, err//trim(w_list(n))//" to file"
                    r_list(n) = "do not write"
                end select
              case DEFAULT
                print*, err//trim(w_list(n))//" from file"
                r_list(n) = "do not read"
            end select

          case('kl')
            select case (trim(turbulence))
              case('kkl')
                select case (trim(previous_flow_type))
                  case('kkl')
                    r_list(n) = "Kl"
                  case DEFAULT
                    print*, err//trim(w_list(n))//" to file"
                    r_list(n) = "do not write"
                end select
              case DEFAULT
                print*, err//trim(w_list(n))//" from file"
                r_list(n) = "do not read"
            end select

          case('tv')
            select case (trim(turbulence))
              case('sa', 'saBC')
                select case (trim(previous_flow_type))
                  case('sa', 'saBC')
                    r_list(n) = "tv"
                  case DEFAULT
                    print*, err//trim(w_list(n))//" to file"
                    r_list(n) = "do not write"
                end select
              case DEFAULT
                print*, err//trim(w_list(n))//" from file"
                r_list(n) = "do not read"
            end select

          case('wall_distance', 'dist', 'wall_dist', 'wdist')
            r_list(n) = "do not read"
            !if(turbulence/="none") then
            !  r_list(n) = "Wall_distance"
            !else
            !  print*, err//trim(r_list(n))//" from file"
            !  r_list(n) = "do not read"
            !end if

          case('intermittency')
            select case (trim(turbulence))
              case('saBC')
                select case(trim(previous_flow_type))
                  case('saBC')
                    r_list(n) = "Intermittency"
                  case DEFAULT
                    print*, err//trim(r_list(n))//" to file"
                    r_list(n) = "do not read"
                end select
              case DEFAULT
                print*, err//trim(r_list(n))//" to file"
                r_list(n) = "do not read"
             end select

          case('tgm')
            select case (trim(transition))
              case('lctm2015')
                select case(trim(previous_flow_type))
                  case('sst','sst2003')
                    r_list(n) = "tgm"
                  case DEFAULT
                    print*, err//trim(r_list(n))//" to file"
                    r_list(n) = "do not read"
                end select
              case DEFAULT
                print*, err//trim(r_list(n))//" to file"
                r_list(n) = "do not read"
             end select

          case('resnorm')
            r_list(n) = "do not read"
            !r_list(n) = "Resnorm"

          case('tke_residue')
            r_list(n) = "do not read"
            !r_list(n) = "TKE_residue"

          case('omega_residue')
            r_list(n) = "do not read"
            !w_list(n) = "Omega_residue"

          case('f1')
            r_list(n) = "do not read"
            !r_list(n) = "F1"

          case('dudx')
            r_list(n) = "do not read"
            !r_list(n) = "Dudx"

          case('dudy')
            r_list(n) = "do not read"
            !r_list(n) = "Dudy"

          case('dudz')
            r_list(n) = "do not read"
            !r_list(n) = "Dudz"

          case('dvdx')
            r_list(n) = "do not read"
            !r_list(n) = "Dvdx"

          case('dvdy')
            r_list(n) = "do not read"
            !r_list(n) = "Dvdy"

          case('dvdz')
            r_list(n) = "do not read"
            !r_list(n) = "Dvdz"

          case('dwdx')
            r_list(n) = "do not read"
            !r_list(n) = "Dwdx"

          case('dwdy')
            r_list(n) = "do not read"
            !r_list(n) = "Dwdy"

          case('dwdz')
            r_list(n) = "do not read"
            !r_list(n) = "Dwdz"

          case('dTdx')
            r_list(n) = "do not read"
            !r_list(n) = "DTdx"

          case('dTdy')
            r_list(n) = "do not read"
            !r_list(n) = "DTdy"

          case('dTdz')
            r_list(n) = "do not read"
            !r_list(n) = "DTdz"

          case('dtkdx')
            r_list(n) = "do not read"
            !r_list(n) = "Dtkdx"

          case('dtkdy')
            r_list(n) = "do not read"
            !r_list(n) = "Dtkdy"

          case('dtkdz')
            r_list(n) = "do not read"
            !r_list(n) = "Dtkdz"

          case('dtwdx')
            r_list(n) = "do not read"
            !r_list(n) = "Dtwdx"

          case('dtwdy')
            r_list(n) = "do not read"
            !r_list(n) = "Dtwdy"

          case('dtwdz')
            r_list(n) = "do not read"
            !r_list(n) = "Dtwdz"

          case('extravar1','extravar2', 'extravar3', 'extravar4', 'extravar5')
            r_list(n) = "do not read"
            !r_list(n) = trim(lcase(w_list(n)))

          case Default
            print*, err//trim(r_list(n))//" from file"
            r_list(n) = "do not read"

        end select
      end do

    end subroutine verify_read_control

end module check_output_control
