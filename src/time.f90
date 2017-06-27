module time

  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx

  use global_vars, only : xnx, xny, xnz !face unit normal x
  use global_vars, only : ynx, yny, ynz !face unit normal y
  use global_vars, only : znx, zny, znz !face unit normal z
  use global_vars, only : xA, yA, zA    !face area
  use global_vars, only : volume
    
  use global_vars, only : n_var
  use global_vars, only : sst_n_var
  use global_vars, only : qp
  use global_vars, only : qp_inf
  use global_vars, only : density
  use global_vars, only : x_speed
  use global_vars, only : y_speed
  use global_vars, only : z_speed
  use global_vars, only : pressure
  use global_vars, only : tk
  use global_vars, only : tw
  use global_vars, only : tk_inf
  use global_vars, only : tw_inf
  use global_vars, only : gm
  use global_vars, only : R_gas
  use global_vars, only : mu_ref
  use global_vars, only : T_ref
  use global_vars, only : Sutherland_temp
  use global_vars, only : Pr

  use global_vars, only : CFL
  use global_vars, only : tolerance
  use global_vars, only : min_iter
  use global_vars, only : max_iters
  use global_vars, only : current_iter
  use global_vars, only : checkpoint_iter
  use global_vars, only : checkpoint_iter_count
  use global_vars, only : time_stepping_method
  use global_vars, only : time_step_accuracy
  use global_vars, only : global_time_step
  use global_vars, only : delta_t
  use global_vars, only : sim_clock
  use global_vars, only : turbulence
  use global_vars, only : supersonic_flag

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

    private
    INTEGER :: &
    nb_ticks_initial, & ! initial value of the clock tick counter
    nb_ticks_final,   & ! final value of the clock tick counter
    nb_ticks_max,     & ! maximum value of the clock counter
    nb_ticks_sec,     & ! number of clock ticks per second
    nb_ticks           ! number of clock ticks of the code
    REAL :: elapsed_time  ! real time in seconds
    REAL :: t1         ! start clock time
    REAL :: t2         ! finish clock time
    real :: cpu_time_elapsed

    ! Public methods
    public :: setup_time
    public :: destroy_time
    public :: compute_time_step

    contains

        subroutine setup_time()
            implicit none
            
            call dmsg(1, 'solver', 'initmisc')
            call alloc(delta_t, 1, imx-1, 1, jmx-1, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for delta_t.')
            CALL SYSTEM_CLOCK(COUNT_RATE=nb_ticks_sec, COUNT_MAX=nb_ticks_max)
            CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
            CALL CPU_TIME(t1)

        end subroutine setup_time

        subroutine destroy_time()
            implicit none
            
            call dmsg(1, 'solver', 'deallocate_misc')
            CALL CPU_TIME(t2)
            CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
            call dealloc(delta_t)

            nb_ticks = nb_ticks_final - nb_ticks_initial
            IF (nb_ticks_final < nb_ticks_initial) &
            nb_ticks = nb_ticks + nb_ticks_max
            elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
            cpu_time_elapsed = t2-t1 
            print*, "SYSTEM clock: ", elapsed_time
            print*, "CPU time    : ", cpu_time_elapsed

        end subroutine destroy_time

        subroutine compute_local_time_step()
            !-----------------------------------------------------------
            ! Compute the time step to be used at each cell center
            !
            ! Local time stepping can be used to get the solution 
            ! advance towards steady state faster. If only the steady
            ! state solution is required, i.e., transients are 
            ! irrelevant, use local time stepping. 
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

        end subroutine compute_local_time_step

        subroutine compute_global_time_step()
            !-----------------------------------------------------------
            ! Compute a common time step to be used at all cell centers
            !
            ! Global time stepping is generally used to get time 
            ! accurate solutions; transients can be studied by 
            ! employing this strategy.
            !-----------------------------------------------------------

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
            !-----------------------------------------------------------
            ! Compute the time step to be used
            !
            ! This calls either compute_global_time_step() or 
            ! compute_local_time_step() based on what 
            ! time_stepping_method is set to.
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

        end subroutine compute_time_step

end module time
