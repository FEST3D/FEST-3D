  !< Allocate memory to the state variables and initialize them
module state
  !< Allocate memory to the state variables and initialize them
  !< The state of the system is defined using the density, velocity and
  !< pressure (primitive variables qp), and trubulent and transition
  !< variables at the cell-center points.
    !-------------------------------------------------------------------
    ! The state module contains the state variables and the methods that
    ! act on them. 
    !-------------------------------------------------------------------
    
#include "debug.h"
#include "error.h"

    use vartypes
    use global, only: FILE_NAME_LENGTH, STATE_FILE_UNIT, OUT_FILE_UNIT, &
            DESCRIPTION_STRING_LENGTH, STRING_BUFFER_LENGTH
    use global_vars, only : qp
    use global_vars, only : density
    use global_vars, only : x_speed
    use global_vars, only : y_speed
    use global_vars, only : z_speed
    use global_vars, only : pressure
    use global_vars, only : tk
    use global_vars, only : tw
    use global_vars, only : te
    use global_vars, only : tv
    use global_vars, only : tkl
    use global_vars, only : tgm
    use global_vars, only : infile
    use utils,       only: alloc, dealloc
    use read_output, only: read_file

    use check_output_control, only : verify_write_control

    implicit none
    private

    integer :: n_var
    integer :: imx, jmx, kmx

    ! Public methods
    public :: setup_state

    contains


        subroutine link_aliases(scheme)
          !< Setup state variable pointers

            implicit none
            type(schemetype), intent(in) :: scheme

            DebugCall("link_aliases")

            density(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 1)
            x_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 2)
            y_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 3)
            z_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 4)
            pressure(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 5)

!            density_inf => qp_inf(1)
!            x_speed_inf => qp_inf(2)
!            y_speed_inf => qp_inf(3)
!            z_speed_inf => qp_inf(4)
!            pressure_inf => qp_inf(5)

            select case (trim(scheme%turbulence))

                case ("none")
                    !include nothing
                    continue
                
                case ("sst", "sst2003", "bsl", "des-sst", "kw")
                    tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 6)
                    tw(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 7)
!                    tk_inf => qp_inf(6)
!                    tw_inf => qp_inf(7)

                case ("kkl")
                    tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 6)
                    tkl(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 7)
!                    tk_inf => qp_inf(6)
!                    tkl_inf => qp_inf(7)

                case ("sa", "saBC")
                    tv(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 6)
!                    tv_inf => qp_inf(6)

                case ("ke")
                    tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 6)
                    te(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 7)
!                    tk_inf => qp_inf(6)
!                    te_inf => qp_inf(7)

                case ("les")
                  continue
                  ! todo

                case DEFAULT
                  Fatal_error

            end select


            ! Transition modeling
            select case(trim(scheme%transition))

              case('lctm2015')
                tgm(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, n_var)
!                tgm_inf => qp_inf(n_var)

              case('bc', 'none')
                !do nothing
                continue

              case DEFAULT
                Fatal_error

            end Select

        end subroutine link_aliases



!        subroutine unlink_aliases()
!          !< Nullify the pointer link
!
!            implicit none
!
!            DebugCall("unlink_aliases")
!
!            nullify(density)
!            nullify(x_speed)
!            nullify(y_speed)
!            nullify(z_speed)
!            nullify(pressure)
!
!            nullify(density_inf)
!            nullify(x_speed_inf)
!            nullify(y_speed_inf)
!            nullify(z_speed_inf)
!            nullify(pressure_inf)
!
!            select case (trim(turbulence))
!
!                case ("none")
!                    continue
!
!                  case ("sst", "sst2003", "bsl", "kw", "des-sst")
!                    nullify(tk)
!                    nullify(tw)
!                    nullify(tk_inf)
!                    nullify(tw_inf)
!
!                case ("kkl")
!                    nullify(tk)
!                    nullify(tkl)
!                    nullify(tk_inf)
!                    nullify(tkl_inf)
!
!                case ("sa", "saBC")
!                    nullify(tv)
!                    nullify(tv_inf)
!
!                case ("ke")
!                    nullify(tk)
!                    nullify(te)
!                    nullify(tk_inf)
!                    nullify(te_inf)
!
!                case ("les")
!                    continue
!                    ! todo
!
!                case DEFAULT
!                  Fatal_error
!
!            end select
!
!
!            !Transition modeling
!            select case(trim(transition))
!
!              case('lctm2015')
!                nullify(tgm)
!                nullify(tgm_inf)
!
!              case('bc', 'none')
!                !do nothing
!                continue
!
!              case DEFAULT
!                Fatal_error
!
!            end Select
!
!        end subroutine unlink_aliases




        subroutine allocate_memory()
            !< Allocate memory to the state variables
            !-----------------------------------------------------------
            implicit none

            DebugCall("allocate_memory")

            ! The state of the system is defined by the primitive 
            ! variables (density, velocity and pressure) at the grid
            ! cell centers. 
            call alloc(qp, -2, imx+2, -2, jmx+2, -2, kmx+2, 1, n_var, AErrMsg("qp"))
!            call alloc(qp_inf, 1, n_var, AErrMsg("qp_inf"))

        end subroutine allocate_memory



!        subroutine deallocate_memory()
!          !< Deallocate memory from the state variable
!
!            implicit none
!
!            DebugCall("allocate_memory")
!
!            call dealloc(qp)
!
!        end subroutine deallocate_memory



        subroutine setup_state(control, scheme, flow, dims)
            !< Setup the state module.
            !< This subroutine should be run before the state variables
            !< are initilized. This subroutine allocates the memory for 
            !< state variables and sets up the aliases to refer to the 
            !< components of the state
            !-----------------------------------------------------------

            implicit none
            type(controltype), intent(inout) :: control
            type(schemetype), intent(in) :: scheme
            type(flowtype), intent(inout) :: flow
            type(extent), intent(in) :: dims

            DebugCall("setup_state")

            n_var = control%n_var

            imx = dims%imx
            jmx = dims%jmx
            kmx = dims%kmx

            call set_n_var_value(control, scheme)
            call allocate_memory()
            call link_aliases(scheme)
            call init_infinity_values(scheme, flow)
            call initstate(control, scheme, flow, dims)

        end subroutine setup_state



!        subroutine destroy_state()
!            !< Destroy the state module.
!            !< This subroutine destroys the state module which includes
!            !< unlinking the aliases for the state components and 
!            !< deallocating the memory held by the state variables
!            !-----------------------------------------------------------
!
!            implicit none
!            
!            DebugCall("destroy_state")
!
!            call unlink_aliases()
!            call deallocate_memory()
!
!        end subroutine destroy_state



        subroutine init_infinity_values(scheme, flow)
            !< Set the values of the infinity variables "qp_inf"
            !-----------------------------------------------------------

            implicit none
            type(schemetype), intent(in) :: scheme
            type(flowtype), intent(inout) :: flow
            
            DebugCall("init_infinity_values")

            !density_inf = free_stream_density
            !x_speed_inf = free_stream_x_speed
            !y_speed_inf = free_stream_y_speed
            !z_speed_inf = free_stream_z_speed
            !pressure_inf = free_stream_pressure
            flow%vel_mag = sqrt(flow%x_speed_inf**2 + flow%y_speed_inf**2 + flow%z_speed_inf**2)
            flow%MInf    = flow%vel_mag/sqrt(flow%gm*flow%pressure_inf/flow%density_inf)
            flow%Reynolds_number = flow%density_inf*flow%vel_mag*1.0/flow%mu_ref
            flow%Turb_intensity_inf = flow%tu_inf/100

            select case (trim(scheme%turbulence))
                
                case ("none")
                    continue

                case ("sst", "sst2003", "bsl")
                    flow%tk_inf = 1.5*((flow%Vel_mag*flow%Turb_Intensity_inf)**2)
                    flow%tw_inf = flow%density_inf*flow%tk_inf/(flow%mu_ref*flow%mu_ratio_inf)

                case ("kkl")
                    flow%tk_inf = 9*(1e-9)*(sound_speed_inf(flow)**2)
                    flow%tkl_inf = 1.5589*(1e-6)*(flow%mu_ref*sound_speed_inf(flow))/flow%density_inf

                case ("sa")
                     flow%tv_inf = flow%mu_ratio_inf*flow%mu_ref/flow%density_inf

                case ("saBC")
                    flow%tv_inf = 0.005*flow%mu_ratio_inf*flow%mu_ref/flow%density_inf

                case ("kw")
                    flow%tk_inf = 1.5*((flow%Vel_mag*flow%Turb_Intensity_inf)**2)
                    flow%tw_inf = flow%density_inf*flow%tk_inf/(flow%mu_ref*flow%mu_ratio_inf)

                case ("ke")
                    flow%tk_inf = 1.5*((flow%Vel_mag*flow%Turb_Intensity_inf)**2)
                    flow%tw_inf = 0.09*flow%density_inf*flow%tk_inf*flow%tk_inf/(flow%mu_ref*flow%mu_ratio_inf)

                case ("des-sst")
                    flow%tk_inf = 1.5*((flow%Vel_mag*flow%Turb_Intensity_inf)**2)
                    flow%tw_inf = flow%density_inf*flow%tk_inf/(flow%mu_ref*flow%mu_ratio_inf)

                case ("les")
                  continue
                  ! todo

                case DEFAULT
                  Fatal_error

            end select


           ! !Transition modeling
           ! select case(trim(scheme%transition))

           !   case('lctm2015')
           !     tgm_inf = free_stream_tgm

           !   case('bc', 'none')
           !     !do nothing
           !     continue

           !   case DEFAULT
           !     Fatal_error

           ! end Select

        end subroutine init_infinity_values


        
        function sound_speed_inf(flow) result(a)
            !< Return the free stream speed of sound.
            !-----------------------------------------------------------

            implicit none
            type(flowtype), intent(in) :: flow
            real :: a

            a = sqrt(flow%gm * flow%pressure_inf / flow%density_inf)

        end function sound_speed_inf



        subroutine initstate(control, scheme, flow, dims)
            !< Initialize the state.
            !< If load file(start_from) is 0, then the state should be 
            !< set to the infinity values. Otherwise, read the state_file
            !< to get the state values
            !-----------------------------------------------------------

            implicit none
            type(extent), intent(in) :: dims
            type(controltype), intent(inout) :: control
            type(schemetype), intent(in) :: scheme
            type(flowtype), intent(in) :: flow
            
            DebugCall("initstate")

            call  verify_write_control(control, scheme, flow)

            if (control%start_from .eq. 0) then
                ! Set the state to the infinity values
                call init_state_with_infinity_values(scheme, flow)
            else
                write(infile,'(a,i4.4,a,i2.2)') &
                  "time_directories/",control%start_from,"/process_",process_id
                ! Set the state to the infinity values so if some
                ! variable are not restart variable they get free_stream value
                call init_state_with_infinity_values(scheme, flow)
                call read_file(control, scheme, dims)

            end if

        end subroutine initstate



        subroutine init_state_with_infinity_values(scheme, flow)
            !< Initialize the state based on the infinity values
            !-----------------------------------------------------------
            
            implicit none
            type(schemetype), intent(in) :: scheme
            type(flowtype), intent(in) :: flow
            
            DebugCall("init_state_with_infinity_values")

            density(-2:imx+2, -2:jmx+2, -2:kmx+2) = flow%density_inf
            x_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) = flow%x_speed_inf
            y_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) = flow%y_speed_inf
            z_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) = flow%z_speed_inf
            pressure(-2:imx+2, -2:jmx+2, -2:kmx+2)= flow%pressure_inf

            select case (trim(scheme%turbulence))

                case ("none")
                    !include nothing
                    continue
                
                case ("sst", "sst2003", "bsl", "des-sst", "kw")
                    tk(-2:imx+2, -2:jmx+2, -2:kmx+2) = flow%tk_inf
                    tw(-2:imx+2, -2:jmx+2, -2:kmx+2) = flow%tw_inf

                case ("kkl")
                    tk(-2:imx+2, -2:jmx+2, -2:kmx+2) = flow%tk_inf
                    tkl(-2:imx+2, -2:jmx+2, -2:kmx+2) = flow%tkl_inf

                case ("sa", "saBC")
                    tv(-2:imx+2, -2:jmx+2, -2:kmx+2) = flow%tv_inf

                case ("ke")
                    tk(-2:imx+2, -2:jmx+2, -2:kmx+2) = flow%tk_inf
                    te(-2:imx+2, -2:jmx+2, -2:kmx+2) = flow%te_inf

                case ("les")
                  continue
                  ! todo

                case DEFAULT
                  Fatal_error

            end select


            ! Transition modeling
            select case(trim(scheme%transition))

              case('lctm2015')
                tgm(-2:imx+2, -2:jmx+2, -2:kmx+2) = flow%tgm_inf

              case('bc', 'none')
                !do nothing
                continue

              case DEFAULT
                Fatal_error

            end Select
            
        end subroutine init_state_with_infinity_values



        subroutine set_n_var_value(control, scheme)
          !< Set number of variable to solver for based on
          !< the tubulence and transition model being used
          implicit none
          type(controltype), intent(inout) :: control
          type(schemetype),  intent(in) ::scheme

          DebugCall("set_n_var_value")

          select case (trim(scheme%turbulence))
            case('none')
              n_var=5

            case('sa', 'saBC')
              n_var=6
              
            case('sst', "sst2003", 'bsl', 'kw', 'ke', 'kkl', 'Des-kw')
              n_var=7

            case DEFAULT
              n_var=5

          end select


          !Transition modeling
          select case(trim(scheme%transition))
            case('lctm2015')
              n_var = n_var + 1

            case('bc', 'none')
              n_var = n_var + 0

            case DEFAULT
              Fatal_error

          end Select

          control%n_var = n_var

        end subroutine set_n_var_value

end module state
