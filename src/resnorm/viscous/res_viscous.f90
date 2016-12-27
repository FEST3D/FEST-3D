module res_viscous
  !----------------------------------------------------
  ! this module contains subroutine that 
  ! 1. check if time for resnorm dump is arrived
  ! 2. calculate resnorm
  ! 3. send those resnorm to processor number 0
  ! 4. Recalulate resnorm based on information 
  !    availble from all processors
  ! 5. Append the data to resnorm file
  !----------------------------------------------------

  use global,      only: RESNORM_FILE_UNIT

  use global_vars, only: gm
  use global_vars, only: density_inf
  use global_vars, only: x_speed_inf
  use global_vars, only: y_speed_inf
  use global_vars, only: z_speed_inf
  use global_vars, only: pressure_inf
  use global_vars, only: mass_residue
  use global_vars, only: x_mom_residue
  use global_vars, only: y_mom_residue
  use global_vars, only: z_mom_residue
  use global_vars, only: energy_residue
  use global_vars, only: vis_resnorm
  use global_vars, only: vis_resnorm_0
  use global_vars, only: vis_resnorm_abs
  use global_vars, only: cont_resnorm
  use global_vars, only: cont_resnorm_0
  use global_vars, only: x_mom_resnorm
  use global_vars, only: x_mom_resnorm_0
  use global_vars, only: y_mom_resnorm
  use global_vars, only: y_mom_resnorm_0
  use global_vars, only: z_mom_resnorm
  use global_vars, only: z_mom_resnorm_0
  use global_vars, only: energy_resnorm
  use global_vars, only: energy_resnorm_0
  use global_vars, only: current_iter
  use global_vars, only: res_write_interval

  use utils,      only: dmsg
  use utils,      only: dealloc
  use utils,      only: alloc
  use layout,     only: process_id
  use layout,     only: total_process

  use mpi

  implicit none
  private

  integer, parameter                :: n = 7 ! total viscous variable + 2
  real                              :: speed_inf
  real, dimension(7)                :: res_send_buf
  real, dimension(:),   allocatable :: root_res_recv_buf
  real, dimension(:,:), allocatable :: global_resnorm

  public :: compute_viscous_resnorm

  contains

    subroutine compute_viscous_resnorm()
      implicit none

      speed_inf = sqrt(x_speed_inf ** 2. + y_speed_inf ** 2. + &
                        z_speed_inf ** 2.)

!      if ( mod(current_iter,res_write_interval) == 0 ) then
        call compute_block_resnorm()
        if (current_iter <= 5) then
          call store_intial_resnorm()
        end if
        if (process_id == 0) then
          call alloc(root_res_recv_buf, 1,total_process*n, &
              errmsg='Error: Unable to allocate memory to root_res_recv_buf')
          root_res_recv_buf = 0.
        end if
        call send_resnorm_to_process_0()
        if (process_id == 0) then
          call alloc(global_resnorm, 1,total_process, 1,n, &
              errmsg='Error: Unable to allocate memory to root_res_recv_buf')
          call recv_resnorm_to_process_0()
          call recalculate_collective_resnorm()
          call dealloc(global_resnorm)
          call dealloc(root_res_recv_buf)
        end if
!      end if

    end subroutine compute_viscous_resnorm

    subroutine compute_block_resnorm()

        implicit none
        
        call dmsg(1, 'solver', 'compute_residue_norm')

        call R_cont()
        call R_velocity()
        call R_energy()

        vis_resnorm = sqrt(                       &
                            cont_resnorm   **2  + &
                            x_mom_resnorm  **2  + &
                            y_mom_resnorm  **2  + &
                            z_mom_resnorm  **2  + &
                            energy_resnorm ** 2   &
                          )

    end subroutine compute_block_resnorm


    subroutine R_energy()
      implicit none

      energy_resnorm = sqrt(                                        &
                          sum(                                      &
                              (                                     &
                               energy_residue(:, :, :) /             &
                              (density_inf * speed_inf *            &
                              ((0.5 * speed_inf * speed_inf) +      &
                              (gm/(gm-1)*pressure_inf/density_inf)))&
                              ) ** 2                                &
                             )                                      &
                           )

    end subroutine R_energy

    subroutine R_velocity()
      implicit none

        x_mom_resnorm = sqrt(                                      &
                          sum(                                     &
                              (x_mom_residue(:, :, :) /            &
                              (density_inf * speed_inf ** 2)) ** 2 &
                             )                                     &
                            )
            
        y_mom_resnorm = sqrt(                                      &
                          sum(                                     &
                              (y_mom_residue(:, :, :) /            &
                              (density_inf * speed_inf ** 2)) ** 2 &
                             )                                     &
                            )
            
        z_mom_resnorm = sqrt(                                      &
                          sum(                                     &
                              (z_mom_residue(:, :, :) /            &
                              (density_inf * speed_inf ** 2)) ** 2 &
                             )                                     &
                            )
        
    end subroutine R_velocity

    subroutine R_cont()
      implicit none

        cont_resnorm = sqrt(                                       &
                          sum(                                     &
                              (mass_residue(:, :, :) /             &
                              (density_inf * speed_inf)) ** 2      &
                             )                                     &
                           )
                    
    end subroutine R_cont

    subroutine store_intial_resnorm()
      implicit none

        vis_resnorm_0     = vis_resnorm
        cont_resnorm_0    = cont_resnorm
        x_mom_resnorm_0   = x_mom_resnorm
        y_mom_resnorm_0   = y_mom_resnorm
        z_mom_resnorm_0   = z_mom_resnorm
        energy_resnorm_0  = energy_resnorm

    end subroutine store_intial_resnorm

    subroutine send_resnorm_to_process_0()
      implicit none
      integer :: ierr


      res_send_buf(1) = vis_resnorm
      res_send_buf(2) = vis_resnorm/vis_resnorm_0
      res_send_buf(3) = cont_resnorm/cont_resnorm_0
      res_send_buf(4) = x_mom_resnorm/x_mom_resnorm_0
      res_send_buf(5) = y_mom_resnorm/y_mom_resnorm_0
      res_send_buf(6) = z_mom_resnorm/z_mom_resnorm_0
      res_send_buf(7) = energy_resnorm/energy_resnorm_0

      call MPI_Gather(res_send_buf, 7, MPI_DOUBLE_PRECISION, &
        root_res_recv_buf, 7, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    end subroutine send_resnorm_to_process_0

    subroutine recv_resnorm_to_process_0()
      implicit none
      integer :: id
      do id = 0,total_process-1
        global_resnorm(id+1, 1) = root_res_recv_buf(1+n*id)
        global_resnorm(id+1, 2) = root_res_recv_buf(2+n*id)
        global_resnorm(id+1, 3) = root_res_recv_buf(3+n*id)
        global_resnorm(id+1, 4) = root_res_recv_buf(4+n*id)
        global_resnorm(id+1, 5) = root_res_recv_buf(5+n*id)
        global_resnorm(id+1, 6) = root_res_recv_buf(6+n*id)
        global_resnorm(id+1, 7) = root_res_recv_buf(7+n*id)
      end do
    end  subroutine recv_resnorm_to_process_0

    subroutine recalculate_collective_resnorm()
      implicit none
      integer :: id

        !all resnorm except first are all normalized
        vis_resnorm_abs = 0.
        vis_resnorm     = 0.
        cont_resnorm    = 0.
        x_mom_resnorm   = 0.
        y_mom_resnorm   = 0.
        z_mom_resnorm   = 0.
        energy_resnorm  = 0.

        do id = 1,total_process
          vis_resnorm_abs = vis_resnorm_abs + (global_resnorm(id,1)**2)
          vis_resnorm     = vis_resnorm     + (global_resnorm(id,2)**2)
          cont_resnorm    = cont_resnorm    + (global_resnorm(id,3)**2)
          x_mom_resnorm   = x_mom_resnorm   + (global_resnorm(id,4)**2)
          y_mom_resnorm   = y_mom_resnorm   + (global_resnorm(id,5)**2)
          z_mom_resnorm   = z_mom_resnorm   + (global_resnorm(id,6)**2)
          energy_resnorm  = energy_resnorm  + (global_resnorm(id,7)**2)
        end do

        vis_resnorm_abs = sqrt(vis_resnorm_abs)       
        vis_resnorm     = sqrt(vis_resnorm)      
        cont_resnorm    = sqrt(cont_resnorm) 
        x_mom_resnorm   = sqrt(x_mom_resnorm) 
        y_mom_resnorm   = sqrt(y_mom_resnorm) 
        z_mom_resnorm   = sqrt(z_mom_resnorm) 
        energy_resnorm  = sqrt(energy_resnorm)

    end subroutine recalculate_collective_resnorm

end module res_viscous

