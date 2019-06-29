  !< Calculate the time step for the current iteration
module time
  !< Calculate the time step for the current iteration

  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx

  use global_vars, only : xnx, xny, xnz !face unit normal x
  use global_vars, only : ynx, yny, ynz !face unit normal y
  use global_vars, only : znx, zny, znz !face unit normal z
  use global_vars, only : xA, yA, zA    !face area
  use global_vars, only : volume
    
  use global_vars, only : n_var
  use global_vars, only : qp
  use global_vars, only : qp_inf
  use global_vars, only : density
  use global_vars, only : x_speed
  use global_vars, only : y_speed
  use global_vars, only : z_speed
  use global_vars, only : pressure
  use global_vars, only : gm
  use global_vars, only : Pr
  use global_vars, only : tPr
  use global_vars, only : R_gas
  use global_vars, only : mu_ref
  use global_vars, only : mu
  use global_vars, only : mu_t

  use global_vars, only : CFL
  use global_vars, only : total_process
  use global_vars, only : process_id
  use global_vars, only : time_stepping_method
  use global_vars, only : time_step_accuracy
  use global_vars, only : global_time_step
  use global_vars, only : delta_t
  use global_vars, only : sim_clock
  use global_vars, only : turbulence

  use utils, only: alloc
  use utils, only:  dealloc 
  use utils, only:  dmsg
  use utils, only:  DEBUG_LEVEL
  use face_interpolant, only: interpolant, &
          x_qp_left, x_qp_right, &
          y_qp_left, y_qp_right, &
          z_qp_left, z_qp_right, compute_face_interpolant, &
          extrapolate_cell_averages_to_faces

  use string
  use read, only : read_input_and_controls
  use geometry, only : CellCenter

#include "mpi.inc"

    private
    INTEGER :: &
    nb_ticks_initial, & !< Initial value of the clock tick counter
    nb_ticks_final,   & !< Final value of the clock tick counter
    nb_ticks_max,     & !< Maximum value of the clock counter
    nb_ticks_sec,     & !< Number of clock ticks per second
    nb_ticks           !< Number of clock ticks of the code
    REAL :: elapsed_time  !< Real time in seconds
    REAL :: t1         !< Start clock time
    REAL :: t2         !< Finish clock time
    real :: cpu_time_elapsed

    ! Public methods
    public :: setup_time
    public :: destroy_time
    public :: compute_time_step
    public :: update_simulation_clock

    contains

        subroutine setup_time()
          !< Allocate memeroy and setup initial clock
            implicit none
            
            call dmsg(1, 'time', 'initmisc')
            call alloc(delta_t, 1, imx-1, 1, jmx-1, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for delta_t.')
            CALL SYSTEM_CLOCK(COUNT_RATE=nb_ticks_sec, COUNT_MAX=nb_ticks_max)
            CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
            CALL CPU_TIME(t1)

        end subroutine setup_time

        subroutine destroy_time()
          !< Deallocate memory and find simulation time.
            implicit none
            real, dimension(:), allocatable :: total_time 
            integer :: ierr
            
            call dmsg(1, 'solver', 'deallocate_misc')

            !simlulation clock data
            if(process_id==0) write(*, '(A)') '>> TIME <<'
            if(process_id==0) write(*, '(A)') "Simulation Clock : "//trim(write_time(sim_clock))
            call alloc(total_time, 1, total_process)
            CALL CPU_TIME(t2)
            CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
            call dealloc(delta_t)

            nb_ticks = nb_ticks_final - nb_ticks_initial
            IF (nb_ticks_final < nb_ticks_initial) &
            nb_ticks = nb_ticks + nb_ticks_max
            elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
            cpu_time_elapsed = t2-t1 
            write(*,'(A,I0,A)') 'process: ',process_id,&
                                " > SYSTEM clock <: "//trim(write_time(elapsed_time))//&
                                " /-\ CPU time <: "//trim(write_time(cpu_time_elapsed))
            
            !total time including all blocks
            call MPI_GATHER(elapsed_time, 1, MPI_DOUBLE_PRECISION, &
            total_time, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
            if(process_id==0) print*, "Total SYSTEM clock: ", trim(write_time(sum(total_time)))
            call MPI_GATHER(cpu_time_elapsed, 1, MPI_DOUBLE_PRECISION, &
            total_time, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
            if(process_id==0) print*, "Total CPU time    : ", trim(write_time(sum(total_time)))
            call dealloc(total_time)

        end subroutine destroy_time

        function write_time(time_in_seconds) result(string)
          !< Particular format to write time in output log file
          implicit none
          real, intent(in) :: time_in_seconds
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

        subroutine compute_local_time_step()
            !< Compute the time step to be used at each cell center
            !<
            !< Local time stepping can be used to get the solution 
            !< advance towards steady state faster. If only the steady
            !< state solution is required, i.e., transients are 
            !< irrelevant, use local time stepping. 
            !-----------------------------------------------------------

            implicit none

            real :: lmx1, lmx2, lmx3, lmx4, lmx5, lmx6, lmxsum
            real :: x_sound_speed_avg, y_sound_speed_avg, z_sound_speed_avg
            integer :: i, j, k

            call dmsg(1, 'solver', 'compute_local_time_step')

            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
               ! For orientation, refer to the report. The standard i,j,k 
               ! direction are marked. All orientation notations are w.r.t 
               ! to the perspective shown in the image.

               ! Faces with lower index
               x_sound_speed_avg = 0.5 * (sqrt(gm * x_qp_left(i, j, k, 5) / &
                                                    x_qp_left(i, j, k, 1)) + &
                                          sqrt(gm * x_qp_right(i, j, k, 5) / &
                                                    x_qp_right(i, j, k, 1)) )
               y_sound_speed_avg = 0.5 * (sqrt(gm * y_qp_left(i, j, k, 5) / &
                                                    y_qp_left(i, j, k, 1)) + &
                                          sqrt(gm * y_qp_right(i, j, k, 5) / &
                                                    y_qp_right(i, j, k, 1)) )
               z_sound_speed_avg = 0.5 * (sqrt(gm * z_qp_left(i, j, k, 5) / &
                                                    z_qp_left(i, j, k, 1)) + &
                                          sqrt(gm * z_qp_right(i, j, k, 5) / &
                                                    z_qp_right(i, j, k, 1)) )
               
               ! For left face: i.e., lower index face along xi direction
               lmx1 = abs( &
                    (x_speed(i, j, k) * xnx(i, j, k)) + &
                    (y_speed(i, j, k) * xny(i, j, k)) + &
                    (z_speed(i, j, k) * xnz(i, j, k))) + &
                    x_sound_speed_avg
               ! For front face, i.e., lower index face along eta direction
               lmx2 = abs( &
                    (x_speed(i, j, k) * ynx(i, j, k)) + &
                    (y_speed(i, j, k) * yny(i, j, k)) + &
                    (z_speed(i, j, k) * ynz(i, j, k))) + &
                    y_sound_speed_avg
               ! For bottom face, i.e., lower index face along zeta direction
               lmx3 = abs( &
                    (x_speed(i, j, k) * znx(i, j, k)) + &
                    (y_speed(i, j, k) * zny(i, j, k)) + &
                    (z_speed(i, j, k) * znz(i, j, k))) + &
                    z_sound_speed_avg

               ! Faces with higher index
               x_sound_speed_avg = 0.5 * (sqrt(gm * x_qp_left(i+1,j,k,5) / x_qp_left(i+1,j,k,1)) + &
                                          sqrt(gm * x_qp_right(i+1,j,k,5) / x_qp_right(i+1,j,k,1)) )
               y_sound_speed_avg = 0.5 * (sqrt(gm * y_qp_left(i,j+1,k,5) / y_qp_left(i,j+1,k,1)) + &
                                          sqrt(gm * y_qp_right(i,j+1,k,5) / y_qp_right(i,j+1,k,1)) )
               z_sound_speed_avg = 0.5 * (sqrt(gm * z_qp_left(i,j,k+1,5) / z_qp_left(i,j,k+1,1)) + &
                                          sqrt(gm * z_qp_right(i,j,k+1,5) / z_qp_right(i,j,k+1,1)) )
               
               ! For right face, i.e., higher index face along xi direction
               lmx4 = abs( &
                    (x_speed(i+1, j, k) * xnx(i+1, j, k)) + &
                    (y_speed(i+1, j, k) * xny(i+1, j, k)) + &
                    (z_speed(i+1, j, k) * xnz(i+1, j, k))) + &
                    x_sound_speed_avg
               ! For back face, i.e., higher index face along eta direction
               lmx5 = abs( &
                    (x_speed(i, j+1, k) * ynx(i, j+1, k)) + &
                    (y_speed(i, j+1, k) * yny(i, j+1, k)) + &
                    (z_speed(i, j+1, k) * ynz(i, j+1, k))) + &
                    y_sound_speed_avg
               ! For top face, i.e., higher index face along zeta direction
               lmx6 = abs( &
                    (x_speed(i, j, k+1) * znx(i, j, k+1)) + &
                    (y_speed(i, j, k+1) * zny(i, j, k+1)) + &
                    (z_speed(i, j, k+1) * znz(i, j, k+1))) + &
                    z_sound_speed_avg

               lmxsum = (xA(i, j, k) * lmx1) + &
                        (yA(i, j, k) * lmx2) + &
                        (zA(i, j, k) * lmx3) + &
                        (xA(i+1, j, k) * lmx4) + &
                        (yA(i, j+1, k) * lmx5) + &
                        (zA(i, j, k+1) * lmx6)
            
               delta_t(i, j, k) = 1. / lmxsum
               delta_t(i, j, k) = delta_t(i, j, k) * volume(i, j, k) * CFL
              end do
             end do
            end do

            if(mu_ref/=0.0) then
              call add_viscous_time()
            end if
            if(mu_ref/=0 .and. trim(turbulence)/='none')then
              call add_turbulent_time()
            end if

        end subroutine compute_local_time_step

        subroutine compute_global_time_step()
            !< Compute a common time step to be used at all cell centers
            !<
            !< Global time stepping is generally used to get time 
            !< accurate solutions; transients can be studied by 
            !< employing this strategy.
            !<-----------------------------------------------------------

            implicit none
            
            call dmsg(1, 'solver', 'compute_global_time_step')

            if (global_time_step > 0) then
                delta_t = global_time_step
            else
                call compute_local_time_step()
                ! The global time step is the minimum of all the local time
                ! steps.
                delta_t = minval(delta_t)
            end if

        end subroutine compute_global_time_step

        subroutine compute_time_step()
            !< Compute the time step to be used
            !<
            !< This calls either compute_global_time_step() or 
            !< compute_local_time_step() based on what 
            !< time_stepping_method is set to.
            !-----------------------------------------------------------

            implicit none
            
            call dmsg(1, 'solver', 'compute_time_step')

            if (time_stepping_method .eq. 'g') then
                call compute_global_time_step()
            else if (time_stepping_method .eq. 'l') then
                call compute_local_time_step()
            else
                call dmsg(5, 'solver', 'compute_time_step', &
                        msg='Value for time_stepping_method (' // &
                            time_stepping_method // ') not recognized.')
                stop
            end if
            !update_simulation clock
            call update_simulation_clock()

        end subroutine compute_time_step


      subroutine update_simulation_clock
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
          if (time_stepping_method .eq. 'g' .and. sim_clock >= 0.) then
              sim_clock = sim_clock + minval(delta_t)
          else if (time_stepping_method .eq. 'l') then
              sim_clock = -1
          end if

      end subroutine update_simulation_clock

      subroutine add_viscous_time()
        !< Addition to local time step due to viscous effects
        implicit none

        real :: lmx1, lmx2, lmx3, lmx4, lmx5, lmx6, lmxsum
        integer :: i, j, k

        call dmsg(1, 'time', 'add_viscous_time_step')

        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1

           ! Faces with lower index

           
           ! For left face: i.e., lower index face along xi direction
           lmx1 = mu(i,j,k)/(density(i,j,k)*abs( &
                ((CellCenter(i-1,j,k,1) - CellCenter(i,j,k,1)) * xnx(i, j, k)) + &
                ((CellCenter(i-1,j,k,2) - CellCenter(i,j,k,2)) * xny(i, j, k)) + &
                ((CellCenter(i-1,j,k,3) - CellCenter(i,j,k,3)) * xnz(i, j, k))))
           ! For front face, i.e., lower index face along eta direction
           lmx2 = mu(i,j,k)/(density(i,j,k)*abs( &
                ((CellCenter(i,j-1,k,1) - CellCenter(i,j,k,1)) * ynx(i, j, k)) + &
                ((CellCenter(i,j-1,k,2) - CellCenter(i,j,k,2)) * yny(i, j, k)) + &
                ((CellCenter(i,j-1,k,3) - CellCenter(i,j,k,3)) * ynz(i, j, k))))
           ! For bottom face, i.e., lower index face along zeta direction
           lmx3 = mu(i,j,k)/(density(i,j,k)*abs( &
                ((CellCenter(i,j,k-1,1) - CellCenter(i,j,k,1)) * znx(i, j, k)) + &
                ((CellCenter(i,j,k-1,2) - CellCenter(i,j,k,2)) * zny(i, j, k)) + &
                ((CellCenter(i,j,k-1,3) - CellCenter(i,j,k,3)) * znz(i, j, k))))

           
           ! For right face, i.e., higher index face along xi direction
           lmx4 = mu(i+1,j,k)/(density(i+1,j,k)*abs( &
                ((CellCenter(i,j,k,1) - CellCenter(i+1,j,k,1)) * xnx(i+1, j, k)) + &
                ((CellCenter(i,j,k,2) - CellCenter(i+1,j,k,2)) * xny(i+1, j, k)) + &
                ((CellCenter(i,j,k,3) - CellCenter(i+1,j,k,3)) * xnz(i+1, j, k))))
           ! For back face, i.e., higher index face along eta direction
           lmx5 = mu(i,j+1,k)/(density(i,j+1,k)*abs( &
                ((CellCenter(i,j,k,1) - CellCenter(i,j+1,k,1)) * ynx(i, j+1, k)) + &
                ((CellCenter(i,j,k,2) - CellCenter(i,j+1,k,2)) * yny(i, j+1, k)) + &
                ((CellCenter(i,j,k,3) - CellCenter(i,j+1,k,3)) * ynz(i, j+1, k))))
           ! For top face, i.e., higher index face along zeta direction
           lmx6 = mu(i,j,k+1)/(density(i,j,k+1)*abs( &
                ((CellCenter(i,j,k,1) - CellCenter(i,j,k+1,1)) * znx(i, j, k+1)) + &
                ((CellCenter(i,j,k,2) - CellCenter(i,j,k+1,2)) * zny(i, j, k+1)) + &
                ((CellCenter(i,j,k,3) - CellCenter(i,j,k+1,3)) * znz(i, j, k+1))))

           lmxsum = (xA(i, j, k) * lmx1) + &
                    (yA(i, j, k) * lmx2) + &
                    (zA(i, j, k) * lmx3) + &
                    (xA(i+1, j, k) * lmx4) + &
                    (yA(i, j+1, k) * lmx5) + &
                    (zA(i, j, k+1) * lmx6)

           lmxsum = gm*lmxsum/Pr

           lmxsum = 2./(lmxsum + (2.*CFL*volume(i,j,k)/delta_t(i,j,k)))
        
           delta_t(i, j, k) = CFL*( lmxsum * volume(i, j, k))
          end do
         end do
        end do
      end subroutine add_viscous_time

      subroutine add_turbulent_time()
        !< Addition to local time step due to turbulence 
        implicit none

        real :: lmx1, lmx2, lmx3, lmx4, lmx5, lmx6, lmxsum
        integer :: i, j, k

        call dmsg(1, 'time', 'add_viscous_time_step')

        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1

           ! Faces with lower index

           
           ! For left face: i.e., lower index face along xi direction
           lmx1 = mu_t(i,j,k)/(density(i,j,k)*abs( &
                ((CellCenter(i-1,j,k,1) - CellCenter(i,j,k,1)) * xnx(i, j, k)) + &
                ((CellCenter(i-1,j,k,2) - CellCenter(i,j,k,2)) * xny(i, j, k)) + &
                ((CellCenter(i-1,j,k,3) - CellCenter(i,j,k,3)) * xnz(i, j, k))))
           ! For front face, i.e., lower index face along eta direction
           lmx2 = mu_t(i,j,k)/(density(i,j,k)*abs( &
                ((CellCenter(i,j-1,k,1) - CellCenter(i,j,k,1)) * ynx(i, j, k)) + &
                ((CellCenter(i,j-1,k,2) - CellCenter(i,j,k,2)) * yny(i, j, k)) + &
                ((CellCenter(i,j-1,k,3) - CellCenter(i,j,k,3)) * ynz(i, j, k))))
           ! For bottom face, i.e., lower index face along zeta direction
           lmx3 = mu_t(i,j,k)/(density(i,j,k)*abs( &
                ((CellCenter(i,j,k-1,1) - CellCenter(i,j,k,1)) * znx(i, j, k)) + &
                ((CellCenter(i,j,k-1,2) - CellCenter(i,j,k,2)) * zny(i, j, k)) + &
                ((CellCenter(i,j,k-1,3) - CellCenter(i,j,k,3)) * znz(i, j, k))))

           
           ! For right face, i.e., higher index face along xi direction
           lmx4 = mu_t(i+1,j,k)/(density(i+1,j,k)*abs( &
                ((CellCenter(i,j,k,1) - CellCenter(i+1,j,k,1)) * xnx(i+1, j, k)) + &
                ((CellCenter(i,j,k,2) - CellCenter(i+1,j,k,2)) * xny(i+1, j, k)) + &
                ((CellCenter(i,j,k,3) - CellCenter(i+1,j,k,3)) * xnz(i+1, j, k))))
           ! For back face, i.e., higher index face along eta direction
           lmx5 = mu_t(i,j+1,k)/(density(i,j+1,k)*abs( &
                ((CellCenter(i,j,k,1) - CellCenter(i,j+1,k,1)) * ynx(i, j+1, k)) + &
                ((CellCenter(i,j,k,2) - CellCenter(i,j+1,k,2)) * yny(i, j+1, k)) + &
                ((CellCenter(i,j,k,3) - CellCenter(i,j+1,k,3)) * ynz(i, j+1, k))))
           ! For top face, i.e., higher index face along zeta direction
           lmx6 = mu_t(i,j,k+1)/(density(i,j,k+1)*abs( &
                ((CellCenter(i,j,k,1) - CellCenter(i,j,k+1,1)) * znx(i, j, k+1)) + &
                ((CellCenter(i,j,k,2) - CellCenter(i,j,k+1,2)) * zny(i, j, k+1)) + &
                ((CellCenter(i,j,k,3) - CellCenter(i,j,k+1,3)) * znz(i, j, k+1))))

           lmxsum = (xA(i, j, k) * lmx1) + &
                    (yA(i, j, k) * lmx2) + &
                    (zA(i, j, k) * lmx3) + &
                    (xA(i+1, j, k) * lmx4) + &
                    (yA(i, j+1, k) * lmx5) + &
                    (zA(i, j, k+1) * lmx6)

           lmxsum = gm*lmxsum/tPr

           lmxsum = 2./(lmxsum + (2.*CFL*volume(i,j,k)/delta_t(i,j,k)))
        
           delta_t(i, j, k) = CFL*( lmxsum * volume(i, j, k))
          end do
         end do
        end do
      end subroutine add_turbulent_time
end module time
