  !< Calculate the time step for the current iteration
module time
  !< Calculate the time step for the current iteration
  use vartypes
  use mpi
  use viscosity, only : mu
  use viscosity, only : mu_t
  use utils, only: alloc
  use face_interpolant, only: &
          x_qp_left, x_qp_right, &
          y_qp_left, y_qp_right, &
          z_qp_left, z_qp_right
  use read, only : read_input_and_controls

#include "debug.h"
#include "error.h"

    private
    integer :: &
    nb_ticks_initial, & !< Initial value of the clock tick counter
    nb_ticks_final,   & !< Final value of the clock tick counter
    nb_ticks_max,     & !< Maximum value of the clock counter
    nb_ticks_sec,     & !< Number of clock ticks per second
    nb_ticks           !< Number of clock ticks of the code
    real(wp) :: elapsed_time  !< real(wp) time in seconds
    real(wp) :: t1         !< Start clock time
    real(wp) :: t2         !< Finish clock time
    real(wp) :: cpu_time_elapsed
    real(wp) :: sim_clock=0.0

    integer :: imx, jmx, kmx, n_var
    ! Public methods
    public :: setup_time
    public :: destroy_time
    public :: compute_time_step
    public :: update_simulation_clock

    contains

        subroutine setup_time(delta_t, control, dims)
          !< Allocate memeroy and setup initial clock
            implicit none
            type(controltype), intent(in) :: control
            !< Control parameters
            type(extent), intent(in) :: dims
            !< Extent of the domain:imx,jmx,kmx
            real(wp), dimension(:,:,:), allocatable, intent(out) :: delta_t
            !< Local time increment value at each cell center

            DebugCall('initmisc')

            imx = dims%imx
            jmx = dims%jmx
            kmx = dims%kmx
            n_var = control%n_var
            call alloc(delta_t, 1, imx-1, 1, jmx-1, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for delta_t.')
            CALL SYSTEM_CLOCK(COUNT_RATE=nb_ticks_sec, COUNT_MAX=nb_ticks_max)
            CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
            CALL CPU_TIME(t1)

        end subroutine setup_time

        subroutine destroy_time(control)
          !< Deallocate memory and find simulation time.
            implicit none
            type(controltype), intent(in) :: control
            !< Control parameters
            real(wp), dimension(:), allocatable :: total_time 
            !< Total time of executation for each block
            integer :: ierr
            !< error variable for mpi communication
            
            DebugCall('deallocate_misc')

            !simlulation clock data
            if(control%process_id==0) write(*, '(A)') '>> TIME <<'
            if(control%process_id==0) write(*, '(A)') "Simulation Clock : "//trim(write_time(sim_clock))
            call alloc(total_time, 1, control%total_process)
            CALL CPU_TIME(t2)
            CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)

            nb_ticks = nb_ticks_final - nb_ticks_initial
            IF (nb_ticks_final < nb_ticks_initial) &
            nb_ticks = nb_ticks + nb_ticks_max
            elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
            cpu_time_elapsed = t2-t1 
            write(*,'(A,I0,A)') 'process: ',control%process_id,&
                                " > SYSTEM clock <: "//trim(write_time(elapsed_time))//&
                                " /-\ CPU time <: "//trim(write_time(cpu_time_elapsed))
            
            !total time including all blocks
            call MPI_GATHER(elapsed_time, 1, MPI_DOUBLE_PRECISION, &
            total_time, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
            if(control%process_id==0) print*, "Total SYSTEM clock: ", trim(write_time(sum(total_time)))
            call MPI_GATHER(cpu_time_elapsed, 1, MPI_DOUBLE_PRECISION, &
            total_time, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
            if(control%process_id==0) print*, "Total CPU time    : ", trim(write_time(sum(total_time)))

        end subroutine destroy_time

        function write_time(time_in_seconds) result(string)
          !< Particular format to write time in output log file
          implicit none
          real(wp), intent(in) :: time_in_seconds
          !< Time to output
          character(len=64):: string
          !< Time as string in particlar format
          if(time_in_seconds>86400) then
            write(string,'(f0.16,2x,A)') time_in_seconds/86400.,"days"
          elseif(time_in_seconds>3600) then
            write(string,'(f0.16,2x,A)') time_in_seconds/3600.,"Hr."
          elseif(time_in_seconds>60) then
            write(string,'(f0.16,2x,A)') time_in_seconds/60.,"Min."
          elseif(time_in_seconds>0) then
            write(string,'(f0.16,2x,A)') time_in_seconds,"Sec."
          else
            write(string,'(A)') "Not Valid"
          end if
        end function write_time

        subroutine compute_local_time_step(qp, delta_t, cells, Ifaces, Jfaces, Kfaces, CFL, scheme, flow, dims)
            !< Compute the time step to be used at each cell center
            !<
            !< Local time stepping can be used to get the solution 
            !< advance towards steady state faster. If only the steady
            !< state solution is required, i.e., transients are 
            !< irrelevant, use local time stepping. 
            !-----------------------------------------------------------

            implicit none
            real(wp), intent(in) :: CFL
            !< CFL number
            type(schemetype), intent(in) :: scheme
            !< finite-volume Schemes
            type(flowtype), intent(in) :: flow
            !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
            type(extent), intent(in) :: dims
            !< Extent of the domain:imx,jmx,kmx
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in), target :: qp
            !< Store primitive variable at cell center
            real(wp) , dimension(1:dims%imx-1, 1:dims%jmx-1, 1:dims%kmx-1), intent(inout) :: delta_t
            !< Local time increment value at each cell center
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Input cell quantities: volume
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Store face quantites for I faces 
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Store face quantites for J faces 
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
            !< Store face quantites for K faces 

            real(wp) :: lmx1, lmx2, lmx3, lmx4, lmx5, lmx6, lmxsum
            real(wp) :: x_sound_speed_avg, y_sound_speed_avg, z_sound_speed_avg
            integer :: i, j, k

            DebugCall('compute_local_time_step')

            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
               ! For orientation, refer to the report. The standard i,j,k 
               ! direction are marked. All orientation notations are w.r.t 
               ! to the perspective shown in the image.

               ! Faces with lower index
               x_sound_speed_avg = 0.5 * (sqrt(flow%gm * x_qp_left(i, j, k, 5) / &
                                                    x_qp_left(i, j, k, 1)) + &
                                          sqrt(flow%gm * x_qp_right(i, j, k, 5) / &
                                                    x_qp_right(i, j, k, 1)) )
               y_sound_speed_avg = 0.5 * (sqrt(flow%gm * y_qp_left(i, j, k, 5) / &
                                                    y_qp_left(i, j, k, 1)) + &
                                          sqrt(flow%gm * y_qp_right(i, j, k, 5) / &
                                                    y_qp_right(i, j, k, 1)) )
               z_sound_speed_avg = 0.5 * (sqrt(flow%gm * z_qp_left(i, j, k, 5) / &
                                                    z_qp_left(i, j, k, 1)) + &
                                          sqrt(flow%gm * z_qp_right(i, j, k, 5) / &
                                                    z_qp_right(i, j, k, 1)) )
               
               ! For left face: i.e., lower index face along xi direction
               lmx1 = abs( &
                    (qp(i, j, k,2) * Ifaces(i, j, k)%nx) + &
                    (qp(i, j, k,3) * Ifaces(i, j, k)%ny) + &
                    (qp(i, j, k,4) * Ifaces(i, j, k)%nz)) + &
                    x_sound_speed_avg
               ! For front face, i.e., lower index face along eta direction
               lmx2 = abs( &
                    (qp(i, j, k,2) * Jfaces(i, j, k)%nx) + &
                    (qp(i, j, k,3) * Jfaces(i, j, k)%ny) + &
                    (qp(i, j, k,4) * Jfaces(i, j, k)%nz)) + &
                    y_sound_speed_avg
               ! For bottom face, i.e., lower index face along zeta direction
               lmx3 = abs( &
                    (qp(i, j, k,2) * Kfaces(i, j, k)%nx) + &
                    (qp(i, j, k,3) * Kfaces(i, j, k)%ny) + &
                    (qp(i, j, k,4) * Kfaces(i, j, k)%nz)) + &
                    z_sound_speed_avg

               ! Faces with higher index
               x_sound_speed_avg = 0.5 * (sqrt(flow%gm * x_qp_left(i+1,j,k,5) / x_qp_left(i+1,j,k,1)) + &
                                          sqrt(flow%gm * x_qp_right(i+1,j,k,5) / x_qp_right(i+1,j,k,1)) )
               y_sound_speed_avg = 0.5 * (sqrt(flow%gm * y_qp_left(i,j+1,k,5) / y_qp_left(i,j+1,k,1)) + &
                                          sqrt(flow%gm * y_qp_right(i,j+1,k,5) / y_qp_right(i,j+1,k,1)) )
               z_sound_speed_avg = 0.5 * (sqrt(flow%gm * z_qp_left(i,j,k+1,5) / z_qp_left(i,j,k+1,1)) + &
                                          sqrt(flow%gm * z_qp_right(i,j,k+1,5) / z_qp_right(i,j,k+1,1)) )
               
               ! For right face, i.e., higher index face along xi direction
               lmx4 = abs( &
                    (qp(i+1, j, k,2) * Ifaces(i+1, j, k)%nx) + &  !x_speed*xnx
                    (qp(i+1, j, k,3) * Ifaces(i+1, j, k)%ny) + &  !y_speed*xny
                    (qp(i+1, j, k,4) * Ifaces(i+1, j, k)%nz)) + & !z_speed*xnz
                    x_sound_speed_avg
               ! For back face, i.e., higher index face along eta direction
               lmx5 = abs( &
                    (qp(i, j+1, k,2) * Jfaces(i, j+1, k)%nx) + &
                    (qp(i, j+1, k,3) * Jfaces(i, j+1, k)%ny) + &
                    (qp(i, j+1, k,4) * Jfaces(i, j+1, k)%nz)) + &
                    y_sound_speed_avg
               ! For top face, i.e., higher index face along zeta direction
               lmx6 = abs( &
                    (qp(i, j, k+1,2) * Kfaces(i, j, k+1)%nx) + &
                    (qp(i, j, k+1,3) * Kfaces(i, j, k+1)%ny) + &
                    (qp(i, j, k+1,4) * Kfaces(i, j, k+1)%nz)) + &
                    z_sound_speed_avg

               lmxsum = (Ifaces(i, j, k)%A * lmx1) + &
                        (Jfaces(i, j, k)%A * lmx2) + &
                        (Kfaces(i, j, k)%A * lmx3) + &
                        (Ifaces(i+1, j, k)%A * lmx4) + &
                        (Jfaces(i, j+1, k)%A * lmx5) + &
                        (Kfaces(i, j, k+1)%A * lmx6)
            
               delta_t(i, j, k) = 1. / lmxsum
               delta_t(i, j, k) = delta_t(i, j, k) * cells(i, j, k)%volume * CFL
              end do
             end do
            end do

            if(flow%mu_ref/=0.0) then
              call add_viscous_time(qp, delta_t, cells, Ifaces, Jfaces, Kfaces, CFL, flow, dims)
            end if
            if(flow%mu_ref/=0 .and. trim(scheme%turbulence)/='none')then
              call add_turbulent_time(qp, delta_t, cells, Ifaces, Jfaces, Kfaces, CFL, flow, dims)
            end if

        end subroutine compute_local_time_step

        subroutine compute_global_time_step(qp, delta_t, cells, Ifaces, Jfaces, Kfaces, CFL, scheme, flow, dims)
            !< Compute a common time step to be used at all cell centers
            !<
            !< Global time stepping is generally used to get time 
            !< accurate solutions; transients can be studied by 
            !< employing this strategy.
            !<-----------------------------------------------------------

            implicit none
            real(wp), intent(in) :: CFL
            !< CFL number
            type(schemetype), intent(in) :: scheme
            !< finite-volume Schemes
            type(flowtype), intent(in) :: flow
            !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
            type(extent), intent(in) :: dims
            !< Extent of the domain:imx,jmx,kmx
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in), target :: qp
            !< Store primitive variable at cell center
            real(wp) , dimension(1:dims%imx-1, 1:dims%jmx-1, 1:dims%kmx-1), intent(inout) :: delta_t
            !< Local time increment value at each cell center
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Input cell quantities: volume
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Store face quantites for I faces 
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Store face quantites for J faces 
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
            !< Store face quantites for K faces 
            
            DebugCall('compute_global_time_step')

            if (scheme%global_time_step > 0) then
                delta_t = scheme%global_time_step
            else
              call compute_local_time_step(qp, delta_t, cells, Ifaces, Jfaces, Kfaces, CFL, scheme, flow, dims)
                ! The global time step is the minimum of all the local time
                ! steps.
                delta_t = minval(delta_t)
            end if

        end subroutine compute_global_time_step

        subroutine compute_time_step(qp, delta_t, CFL, cells, Ifaces, Jfaces, Kfaces, scheme, flow, dims)
            !< Compute the time step to be used
            !<
            !< This calls either compute_global_time_step() or 
            !< compute_local_time_step() based on what 
            !< time_stepping_method is set to.
            !-----------------------------------------------------------

            implicit none
            real(wp), intent(in) :: CFL
            !< CFL number
            type(schemetype), intent(in) :: scheme
            !< finite-volume Schemes
            type(flowtype), intent(in) :: flow
            !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
            type(extent), intent(in) :: dims
            !< Extent of the domain:imx,jmx,kmx
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in), target :: qp
            !< Store primitive variable at cell center
            real(wp) , dimension(1:dims%imx-1, 1:dims%jmx-1, 1:dims%kmx-1), intent(inout) :: delta_t
            !< Local time increment value at each cell center
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Input cell quantities: volume
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Store face quantites for I faces 
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Store face quantites for J faces 
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
            !< Store face quantites for K faces 
            
            DebugCall('compute_time_step')

            if (scheme%time_stepping_method .eq. 'g') then
                call compute_global_time_step(qp, delta_t, cells, Ifaces, Jfaces, Kfaces, CFL, scheme, flow, dims)
            else if (scheme%time_stepping_method .eq. 'l') then
                call compute_local_time_step(qp, delta_t, cells, Ifaces, Jfaces, Kfaces, CFL, scheme, flow, dims)
            else
                print*,'In compute_time_step: value for time_stepping_method (' //scheme%time_stepping_method // ') not recognized.'
                Fatal_error
            end if
            !update_simulation clock
            call update_simulation_clock(delta_t, scheme, dims)

        end subroutine compute_time_step


      subroutine update_simulation_clock(delta_t, scheme, dims)
          !<  Update the simulation clock
          !< 
          !<  It is sometimes useful to know what the simulation time is
          !<  at every iteration so that a comparison with an analytical
          !<  solution is possible. Since, the global timesteps used may
          !<  not be uniform, we need to track this explicitly.
          !< 
          !<  Of course, it makes sense to track this only if the time 
          !<  stepping is global and not local. If the time stepping is
          !<  local, the simulation clock is set to -1. If it is global
          !<  it is incremented according to the time step found.
          !-----------------------------------------------------------

          implicit none
          type(extent), intent(in) :: dims
          !< Extent of the domain:imx,jmx,kmx
          type(schemetype), intent(in) :: scheme
          !< finite-volume Schemes: time stepping methods
          real(wp) , dimension(1:dims%imx-1, 1:dims%jmx-1, 1:dims%kmx-1), intent(in) :: delta_t
          !< Local time increment value at each cell center
          if (scheme%time_stepping_method .eq. 'g' .and. sim_clock >= 0.) then
              sim_clock = sim_clock + minval(delta_t)
          else if (scheme%time_stepping_method .eq. 'l') then
              sim_clock = -1
          end if

      end subroutine update_simulation_clock

      subroutine add_viscous_time(qp, delta_t, cells, Ifaces, Jfaces, Kfaces, CFL, flow, dims)
        !< Addition to local time step due to viscous effects
        implicit none

        real(wp), intent(in) :: CFL
        !< CFL number
        type(flowtype), intent(in) :: flow
        !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
        type(extent), intent(in) :: dims
        !< Extent of the domain:imx,jmx,kmx
        real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in), target :: qp
        !< Store primitive variable at cell center
        real(wp) , dimension(1:dims%imx-1, 1:dims%jmx-1, 1:dims%kmx-1), intent(inout) :: delta_t
        !< Local time increment value at each cell center
        type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
        !< Input cell quantities: volume
        type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
        !< Store face quantites for I faces 
        type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
        !< Store face quantites for J faces 
        type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
        !< Store face quantites for K faces 
        real(wp) :: lmx1, lmx2, lmx3, lmx4, lmx5, lmx6, lmxsum
        integer :: i, j, k

        DebugCall('add_viscous_time_step')

        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1

           ! Faces with lower index

           
           ! For left face: i.e., lower index face along xi direction
           lmx1 = mu(i,j,k)/(qp(i,j,k,1)*abs( &
                ((cells(i-1,j,k)%centerx - cells(i,j,k)%centerx) * Ifaces(i, j, k)%nx) + &
                ((cells(i-1,j,k)%centery - cells(i,j,k)%centery) * Ifaces(i, j, k)%ny) + &
                ((cells(i-1,j,k)%centerz - cells(i,j,k)%centerz) * Ifaces(i, j, k)%nz)))
           ! For front face, i.e., lower index face along eta direction
           lmx2 = mu(i,j,k)/(qp(i,j,k,1)*abs( &
                ((cells(i,j-1,k)%centerx - cells(i,j,k)%centerx) * Jfaces(i, j, k)%nx) + &
                ((cells(i,j-1,k)%centery - cells(i,j,k)%centery) * Jfaces(i, j, k)%ny) + &
                ((cells(i,j-1,k)%centerz - cells(i,j,k)%centerz) * Jfaces(i, j, k)%nz)))
           ! For bottom face, i.e., lower index face along zeta direction
           lmx3 = mu(i,j,k)/(qp(i,j,k,1)*abs( &
                ((cells(i,j,k-1)%centerx - cells(i,j,k)%centerx) * Kfaces(i, j, k)%nx) + &
                ((cells(i,j,k-1)%centery - cells(i,j,k)%centery) * Kfaces(i, j, k)%ny) + &
                ((cells(i,j,k-1)%centerz - cells(i,j,k)%centerz) * Kfaces(i, j, k)%nz)))

           
           ! For right face, i.e., higher index face along xi direction
           lmx4 = mu(i+1,j,k)/(qp(i+1,j,k,1)*abs( &
                ((cells(i,j,k)%centerx - cells(i+1,j,k)%centerx) * Ifaces(i+1, j, k)%nx) + &
                ((cells(i,j,k)%centery - cells(i+1,j,k)%centery) * Ifaces(i+1, j, k)%ny) + &
                ((cells(i,j,k)%centerz - cells(i+1,j,k)%centerz) * Ifaces(i+1, j, k)%nz)))
           ! For back face, i.e., higher index face along eta direction
           lmx5 = mu(i,j+1,k)/(qp(i,j+1,k,1)*abs( &
                ((cells(i,j,k)%centerx - cells(i,j+1,k)%centerx) * Jfaces(i, j+1, k)%nx) + &
                ((cells(i,j,k)%centery - cells(i,j+1,k)%centery) * Jfaces(i, j+1, k)%ny) + &
                ((cells(i,j,k)%centerz - cells(i,j+1,k)%centerz) * Jfaces(i, j+1, k)%nz)))
           ! For top face, i.e., higher index face along zeta direction
           lmx6 = mu(i,j,k+1)/(qp(i,j,k+1,1)*abs( &
                ((cells(i,j,k)%centerx - cells(i,j,k+1)%centerx) * Kfaces(i, j, k+1)%nx) + &
                ((cells(i,j,k)%centery - cells(i,j,k+1)%centery) * Kfaces(i, j, k+1)%ny) + &
                ((cells(i,j,k)%centerz - cells(i,j,k+1)%centerz) * Kfaces(i, j, k+1)%nz)))

               lmxsum = (Ifaces(i, j, k)%A * lmx1) + &
                        (Jfaces(i, j, k)%A * lmx2) + &
                        (Kfaces(i, j, k)%A * lmx3) + &
                        (Ifaces(i+1, j, k)%A * lmx4) + &
                        (Jfaces(i, j+1, k)%A * lmx5) + &
                        (Kfaces(i, j, k+1)%A * lmx6)

           lmxsum = flow%gm*lmxsum/flow%Pr

           lmxsum = 2./(lmxsum + (2.*CFL*cells(i,j,k)%volume/delta_t(i,j,k)))
        
           delta_t(i, j, k) = CFL*( lmxsum * cells(i, j, k)%volume)
          end do
         end do
        end do
      end subroutine add_viscous_time

      subroutine add_turbulent_time(qp,delta_t,cells,Ifaces,Jfaces,Kfaces,CFL,flow,dims)
        !< Addition to local time step due to turbulence 
        implicit none
        real(wp), intent(in) :: CFL
        !< CFL number
        type(flowtype), intent(in) :: flow
        !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
        type(extent), intent(in) :: dims
        !< Extent of the domain:imx,jmx,kmx
        real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in), target :: qp
        !< Store primitive variable at cell center
        real(wp) , dimension(1:dims%imx-1, 1:dims%jmx-1, 1:dims%kmx-1), intent(inout) :: delta_t
        !< Local time increment value at each cell center
        type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
        !< Input cell quantities: volume
        type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
        !< Store face quantites for I faces 
        type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
        !< Store face quantites for J faces 
        type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
        !< Store face quantites for K faces 
        real(wp) :: lmx1, lmx2, lmx3, lmx4, lmx5, lmx6, lmxsum
        integer :: i, j, k

        DebugCall('add_viscous_time_step')

        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1

           ! Faces with lower index

           
           ! For left face: i.e., lower index face along xi direction
           lmx1 = mu_t(i,j,k)/(qp(i,j,k,1)*abs( &
                ((cells(i-1,j,k)%centerx - cells(i,j,k)%centerx) * Ifaces(i, j, k)%nx) + &
                ((cells(i-1,j,k)%centery - cells(i,j,k)%centery) * Ifaces(i, j, k)%ny) + &
                ((cells(i-1,j,k)%centerz - cells(i,j,k)%centerz) * Ifaces(i, j, k)%nz)))
           ! For front face, i.e., lower index face along eta direction
           lmx2 = mu_t(i,j,k)/(qp(i,j,k,1)*abs( &
                ((cells(i,j-1,k)%centerx - cells(i,j,k)%centerx) * Jfaces(i, j, k)%nx) + &
                ((cells(i,j-1,k)%centery - cells(i,j,k)%centery) * Jfaces(i, j, k)%ny) + &
                ((cells(i,j-1,k)%centerz - cells(i,j,k)%centerz) * Jfaces(i, j, k)%nz)))
           ! For bottom face, i.e., lower index face along zeta direction
           lmx3 = mu_t(i,j,k)/(qp(i,j,k,1)*abs( &
                ((cells(i,j,k-1)%centerx - cells(i,j,k)%centerx) * Kfaces(i, j, k)%nx) + &
                ((cells(i,j,k-1)%centery - cells(i,j,k)%centery) * Kfaces(i, j, k)%ny) + &
                ((cells(i,j,k-1)%centerz - cells(i,j,k)%centerz) * Kfaces(i, j, k)%nz)))

           
           ! For right face, i.e., higher index face along xi direction
           lmx4 = mu_t(i+1,j,k)/(qp(i+1,j,k,1)*abs( &
                ((cells(i,j,k)%centerx - cells(i+1,j,k)%centerx) * Ifaces(i+1, j, k)%nx) + &
                ((cells(i,j,k)%centery - cells(i+1,j,k)%centery) * Ifaces(i+1, j, k)%ny) + &
                ((cells(i,j,k)%centerz - cells(i+1,j,k)%centerz) * Ifaces(i+1, j, k)%nz)))
           ! For back face, i.e., higher index face along eta direction
           lmx5 = mu_t(i,j+1,k)/(qp(i,j+1,k,1)*abs( &
                ((cells(i,j,k)%centerx - cells(i,j+1,k)%centerx) * Jfaces(i, j+1, k)%nx) + &
                ((cells(i,j,k)%centery - cells(i,j+1,k)%centery) * Jfaces(i, j+1, k)%ny) + &
                ((cells(i,j,k)%centerz - cells(i,j+1,k)%centerz) * Jfaces(i, j+1, k)%nz)))
           ! For top face, i.e., higher index face along zeta direction
           lmx6 = mu_t(i,j,k+1)/(qp(i,j,k+1,1)*abs( &
                ((cells(i,j,k)%centerx - cells(i,j,k+1)%centerx) * Kfaces(i, j, k+1)%nx) + &
                ((cells(i,j,k)%centery - cells(i,j,k+1)%centery) * Kfaces(i, j, k+1)%ny) + &
                ((cells(i,j,k)%centerz - cells(i,j,k+1)%centerz) * Kfaces(i, j, k+1)%nz)))

               lmxsum = (Ifaces(i, j, k)%A * lmx1) + &
                        (Jfaces(i, j, k)%A * lmx2) + &
                        (Kfaces(i, j, k)%A * lmx3) + &
                        (Ifaces(i+1, j, k)%A * lmx4) + &
                        (Jfaces(i, j+1, k)%A * lmx5) + &
                        (Kfaces(i, j, k+1)%A * lmx6)

           lmxsum = flow%gm*lmxsum/flow%tPr

           lmxsum = 2./(lmxsum + (2.*CFL*cells(i,j,k)%volume/delta_t(i,j,k)))
        
           delta_t(i, j, k) = CFL*( lmxsum * cells(i, j, k)%volume)
          end do
         end do
        end do
      end subroutine add_turbulent_time
end module time
