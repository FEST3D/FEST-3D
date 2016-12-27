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
  use global_vars, only:    vis_resnorm
  use global_vars, only:   cont_resnorm
  use global_vars, only:  x_mom_resnorm
  use global_vars, only:  y_mom_resnorm
  use global_vars, only:  z_mom_resnorm
  use global_vars, only: energy_resnorm
  use global_vars, only:    vis_resnorm_0
  use global_vars, only:   cont_resnorm_0
  use global_vars, only:  x_mom_resnorm_0
  use global_vars, only:  y_mom_resnorm_0
  use global_vars, only:  z_mom_resnorm_0
  use global_vars, only: energy_resnorm_0
  use global_vars, only:    vis_resnorm_d1
  use global_vars, only:   cont_resnorm_d1
  use global_vars, only:  x_mom_resnorm_d1
  use global_vars, only:  y_mom_resnorm_d1
  use global_vars, only:  z_mom_resnorm_d1
  use global_vars, only: energy_resnorm_d1
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

  integer, parameter                :: n = 12 ! total viscous variable + 2
  real                              :: speed_inf
  real, dimension(n)                :: res_send_buf
  real, dimension(:),   allocatable :: root_res_recv_buf
  real, dimension(:,:), allocatable :: global_resnorm

  public :: compute_viscous_resnorm

  contains

    subroutine compute_viscous_resnorm()
      implicit none
      call dmsg(1, 'res_viscous', 'compute_viscous_resnorm')

      speed_inf = sqrt(x_speed_inf ** 2. + y_speed_inf ** 2. + &
                        z_speed_inf ** 2.)

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

    end subroutine compute_viscous_resnorm

    subroutine compute_block_resnorm()

        implicit none
        
        call dmsg(1, 'solver', 'compute_residue_norm')

        call R_cont()
        call R_velocity()
        call R_energy()

        vis_resnorm =    (                    &
                            cont_resnorm    + &
                            x_mom_resnorm   + &
                            y_mom_resnorm   + &
                            z_mom_resnorm   + &
                            energy_resnorm    &
                          )

    end subroutine compute_block_resnorm


    subroutine R_energy()
      implicit none

      call dmsg(1, 'res_viscous', 'R_energy')

      ! Important sqrt omitted to be taken care in the 

      energy_resnorm = sum(                                       &
                            (                                     &
                             energy_residue(:, :, :) /            &
                            (density_inf * speed_inf *            &
                            ((0.5 * speed_inf * speed_inf) +      &
                            (gm/(gm-1)*pressure_inf/density_inf)))&
                            ) ** 2                                &
                          )                                       

    end subroutine R_energy

    subroutine R_velocity()
      implicit none

      call dmsg(1, 'res_viscous', 'R_velocity')
        x_mom_resnorm = sum(                                     &
                            (x_mom_residue(:, :, :) /            &
                            (density_inf * speed_inf ** 2)) ** 2 &
                           )                                     
            
        y_mom_resnorm = sum(                                     &
                            (y_mom_residue(:, :, :) /            &
                            (density_inf * speed_inf ** 2)) ** 2 &
                           )                                     
            
        z_mom_resnorm = sum(                                     &
                            (z_mom_residue(:, :, :) /            &
                            (density_inf * speed_inf ** 2)) ** 2 &
                           )                                     
        
    end subroutine R_velocity

    subroutine R_cont()
      implicit none

      call dmsg(1, 'res_viscous', 'R_cont')
        cont_resnorm = sum(                                     &
                           (mass_residue(:, :, :) /             &
                           (density_inf * speed_inf)) ** 2      &
                          )                                     
                    
    end subroutine R_cont

    subroutine store_intial_resnorm()
      implicit none

      call dmsg(1, 'res_viscous', 'store_initial_resnorm')
           vis_resnorm_0  =    vis_resnorm
          cont_resnorm_0  =   cont_resnorm
         x_mom_resnorm_0  =  x_mom_resnorm
         y_mom_resnorm_0  =  y_mom_resnorm
         z_mom_resnorm_0  =  z_mom_resnorm
        energy_resnorm_0  = energy_resnorm

    end subroutine store_intial_resnorm

    subroutine send_resnorm_to_process_0()
      implicit none
      integer :: ierr

      call dmsg(1, 'res_viscous', 'send_resnorm_to_process_0')

      res_send_buf(1)  =    vis_resnorm
      res_send_buf(2)  =   cont_resnorm
      res_send_buf(3)  =  x_mom_resnorm
      res_send_buf(4)  =  y_mom_resnorm
      res_send_buf(5)  =  z_mom_resnorm
      res_send_buf(6)  = energy_resnorm
      res_send_buf(7)  =    vis_resnorm_0
      res_send_buf(8)  =   cont_resnorm_0
      res_send_buf(9)  =  x_mom_resnorm_0
      res_send_buf(10) =  y_mom_resnorm_0
      res_send_buf(11) =  z_mom_resnorm_0
      res_send_buf(12) = energy_resnorm_0

      call MPI_Gather(res_send_buf, n, MPI_DOUBLE_PRECISION, &
        root_res_recv_buf, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    end subroutine send_resnorm_to_process_0

    subroutine recv_resnorm_to_process_0()
      implicit none
      integer :: id
      call dmsg(1, 'res_viscous', 'recv_resnorm_to_process_0')
      do id = 0,total_process-1
        global_resnorm(id+1, 1)  = root_res_recv_buf(1  + n*id)
        global_resnorm(id+1, 2)  = root_res_recv_buf(2  + n*id)
        global_resnorm(id+1, 3)  = root_res_recv_buf(3  + n*id)
        global_resnorm(id+1, 4)  = root_res_recv_buf(4  + n*id)
        global_resnorm(id+1, 5)  = root_res_recv_buf(5  + n*id)
        global_resnorm(id+1, 6)  = root_res_recv_buf(6  + n*id)
        global_resnorm(id+1, 7)  = root_res_recv_buf(7  + n*id)
        global_resnorm(id+1, 8)  = root_res_recv_buf(8  + n*id)
        global_resnorm(id+1, 9)  = root_res_recv_buf(9  + n*id)
        global_resnorm(id+1,10)  = root_res_recv_buf(10 + n*id)
        global_resnorm(id+1,11)  = root_res_recv_buf(11 + n*id)
        global_resnorm(id+1,12)  = root_res_recv_buf(12 + n*id)
      end do
    end  subroutine recv_resnorm_to_process_0

    subroutine recalculate_collective_resnorm()
      implicit none
      integer :: id

      call dmsg(1, 'res_viscous', 'recalculate_collective_resnorm')
        !all resnorm except first are all normalized
!           vis_resnorm     = 0.
!          cont_resnorm     = 0.
!         x_mom_resnorm     = 0.
!         y_mom_resnorm     = 0.
!         z_mom_resnorm     = 0.
!        energy_resnorm     = 0.
!           vis_resnorm_0   = 0.
!          cont_resnorm_0   = 0.
!         x_mom_resnorm_0   = 0.
!         y_mom_resnorm_0   = 0.
!         z_mom_resnorm_0   = 0.
!        energy_resnorm_0   = 0.

!        do id = 1,total_process
!             vis_resnorm   =    vis_resnorm    + (global_resnorm(id,1 ))
!            cont_resnorm   =   cont_resnorm    + (global_resnorm(id,2 ))
!           x_mom_resnorm   =  x_mom_resnorm    + (global_resnorm(id,3 ))
!           y_mom_resnorm   =  y_mom_resnorm    + (global_resnorm(id,4 ))
!           z_mom_resnorm   =  z_mom_resnorm    + (global_resnorm(id,5 ))
!          energy_resnorm   = energy_resnorm    + (global_resnorm(id,6 ))
!             vis_resnorm_0 =    vis_resnorm_0  + (global_resnorm(id,7 ))
!            cont_resnorm_0 =   cont_resnorm_0  + (global_resnorm(id,8 ))
!           x_mom_resnorm_0 =  x_mom_resnorm_0  + (global_resnorm(id,9 ))
!           y_mom_resnorm_0 =  y_mom_resnorm_0  + (global_resnorm(id,10))
!           z_mom_resnorm_0 =  z_mom_resnorm_0  + (global_resnorm(id,11))
!          energy_resnorm_0 = energy_resnorm_0  + (global_resnorm(id,12))
!        end do
           vis_resnorm    = SUM(global_resnorm(: ,1 ))
          cont_resnorm    = SUM(global_resnorm(: ,2 ))
         x_mom_resnorm    = SUM(global_resnorm(: ,3 ))
         y_mom_resnorm    = SUM(global_resnorm(: ,4 ))
         z_mom_resnorm    = SUM(global_resnorm(: ,5 ))
        energy_resnorm    = SUM(global_resnorm(: ,6 ))
           vis_resnorm_0  = SUM(global_resnorm(: ,7 ))
          cont_resnorm_0  = SUM(global_resnorm(: ,8 ))
         x_mom_resnorm_0  = SUM(global_resnorm(: ,9 ))
         y_mom_resnorm_0  = SUM(global_resnorm(: ,10))
         z_mom_resnorm_0  = SUM(global_resnorm(: ,11))
        energy_resnorm_0  = SUM(global_resnorm(: ,12))

           vis_resnorm    = SQRT(   vis_resnorm  )
          cont_resnorm    = SQRT(  cont_resnorm  )
         x_mom_resnorm    = SQRT( x_mom_resnorm  )
         y_mom_resnorm    = SQRT( y_mom_resnorm  )
         z_mom_resnorm    = SQRT( z_mom_resnorm  )
        energy_resnorm    = SQRT(energy_resnorm  )
           vis_resnorm_0  = SQRT(   vis_resnorm_0)
          cont_resnorm_0  = SQRT(  cont_resnorm_0)
         x_mom_resnorm_0  = SQRT( x_mom_resnorm_0)
         y_mom_resnorm_0  = SQRT( y_mom_resnorm_0)
         z_mom_resnorm_0  = SQRT( z_mom_resnorm_0)
        energy_resnorm_0  = SQRT(energy_resnorm_0)

           vis_resnorm_d1 =    vis_resnorm /    vis_resnorm_0
          cont_resnorm_d1 =   cont_resnorm /   cont_resnorm_0
         x_mom_resnorm_d1 =  x_mom_resnorm /  x_mom_resnorm_0
         y_mom_resnorm_d1 =  y_mom_resnorm /  y_mom_resnorm_0
         z_mom_resnorm_d1 =  z_mom_resnorm /  z_mom_resnorm_0
        energy_resnorm_d1 = energy_resnorm / energy_resnorm_0

    end subroutine recalculate_collective_resnorm

end module res_viscous

