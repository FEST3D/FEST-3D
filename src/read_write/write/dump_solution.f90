!<  Check, create, and purge folder in the time_directory folder
module dump_solution
  !< This module contians subroutine that
  !<  1. check if point of dumping condition is arrived.
  !<  2. create particular folder for dump.
  !<  3. dump data in that folder.
  !<  4. purge folders if required.
  !------------------------------------------
#include "../../debug.h"
#include "../../error.h"
  use vartypes
  use global,      only : FILE_NAME_LENGTH
  use global,      only : RESTART_FILE_UNIT
!  use global_vars, only : current_iter
!  use global_vars, only : max_iters
!  use global_vars, only : last_iter
!  use global_vars, only : checkpoint_iter
!  use global_vars, only : checkpoint_iter_count
!  use global_vars, only : purge_write
  use global_vars, only : sim_clock
  use global_vars, only :     outfile
  use global_vars, only : restartfile
  use global_vars, only :        resnorm_0
  use global_vars, only :    vis_resnorm_0
  use global_vars, only :   turb_resnorm_0
  use global_vars, only :   cont_resnorm_0
  use global_vars, only :  x_mom_resnorm_0
  use global_vars, only :  y_mom_resnorm_0
  use global_vars, only :  z_mom_resnorm_0
  use global_vars, only : energy_resnorm_0
  use global_vars, only :    TKE_resnorm_0
  use global_vars, only :  omega_resnorm_0
  use global_vars, only :  turbulence
  use utils
!  use string
  use write_output, only : write_file
  use layout,      only : process_id

  implicit none
  private
  character(len=FILE_NAME_LENGTH) :: dump_dirname
  !< Name(check point number) of the directory to create
  character(len=FILE_NAME_LENGTH) :: purge_dirname
  !< Name(check point number) of the directory to remove

  public :: checkpoint

  contains

    subroutine checkpoint(nodes, control, dims)
      !< Create a checkpoint dump file if the time has come
      !-----------------------------------------------------------

      implicit none
      type(extent), intent(in) :: dims
      type(controltype), intent(inout) :: control
      type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(in) :: nodes


      DebugCall('checkpoint')

      if (control%checkpoint_iter .ne. 0) then
          if (mod(control%current_iter, control%checkpoint_iter) == 0 &
             .or. control%current_iter == control%max_iters) then
              call make_dump_dir(control)
              call dump_data(nodes, control, dims)
              print*, "writing data at: ", control%current_iter, control%checkpoint_iter_count
              call purge_dump_dir(control)
              control%checkpoint_iter_count = control%checkpoint_iter_count + 1
          end if
      end if

    end subroutine checkpoint

    subroutine create_directory(dirname)
      !< Create a directory to keep the solution files from all the processor
      implicit none
      character(len=*), intent(in)    :: dirname
      character(len=FILE_NAME_LENGTH) :: mkdircmd

      mkdircmd = 'mkdir -p '//trim(dirname)
      call system(mkdircmd)

    end subroutine create_directory

    subroutine remove_directory(dirname)
      !< Remove a directory 
      implicit none
      character(len=*), intent(in)    :: dirname
      character(len=FILE_NAME_LENGTH) :: rmdircmd

      rmdircmd = 'rm -rf '//trim(dirname)
      call system(rmdircmd)

    end subroutine remove_directory

    subroutine purge_dump_dir(control)
      !< Purge the directory based on the input
      implicit none
      type(controltype), intent(in) :: control
      integer                         :: purge_num

      purge_num = control%checkpoint_iter_count-control%purge_write
      if (control%purge_write /=0 .and. purge_num > 0) then
        write(purge_dirname,'(A,I4.4)') 'time_directories/', purge_num
        call remove_directory(purge_dirname)
      end if

    end subroutine purge_dump_dir

    subroutine make_dump_dir(control)
      !< Solution directory and sub-directory in created with particular number 
      implicit none
      type(controltype), intent(in) :: control

      write(dump_dirname,'(A,I4.4)') 'time_directories/',control%checkpoint_iter_count
      call create_directory(dump_dirname)
      call create_directory(trim(dump_dirname)//'/restart')

    end subroutine make_dump_dir

    subroutine dump_data(nodes, control, dims)
      !< Call to write save files in the directory
      implicit none
      type(extent), intent(in) :: dims
      type(controltype), intent(in) :: control
      type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(in) :: nodes
!      character(len=FILE_NAME_LENGTH) :: filename

      DebugCall('dump_solution: dump_data')
      write(restartfile, '(A,I2.2)') trim(dump_dirname)//'/restart/process_',process_id
      write(    outfile, '(A,I2.2)') trim(dump_dirname)//'/process_',process_id
      call write_restart_log(control)
      call write_file(nodes, control, dims)

    end subroutine dump_data

    subroutine write_restart_log(control)
      !< Call to write log file in the subdirectory "restart". 
      !< It is useful information while restarting the solver
      implicit none
      type(controltype), intent(in) :: control
      open(RESTART_FILE_UNIT, file=restartfile)
      select case (turbulence)
          
        case ('none')
          write(RESTART_FILE_UNIT, '(A)') 'viscous'
        case('sst','sst2003', 'kkl', 'ke', 'kw', 'sa', 'saBC', 'des-sst')
          write(RESTART_FILE_UNIT, '(A)') trim(turbulence)
        case DEFAULT
           Fatal_error
      end select
      call write_initial_resnorm(control%current_iter, control%last_iter)
      close(RESTART_FILE_UNIT)

    end subroutine write_restart_log

    subroutine write_initial_resnorm(current_iter, last_iter)
      !< Writing Initial resnorom in the log file to 
      !< maintian continuity of resnorm while restrarting
      implicit none
      integer, intent(in) :: current_iter, last_iter
      write(RESTART_FILE_UNIT, '(I0)')    current_iter+last_iter
      write(RESTART_FILE_UNIT, '(f0.16)')        resnorm_0
      write(RESTART_FILE_UNIT, '(f0.16)')    vis_resnorm_0
      write(RESTART_FILE_UNIT, '(f0.16)')   turb_resnorm_0
      write(RESTART_FILE_UNIT, '(f0.16)')   cont_resnorm_0
      write(RESTART_FILE_UNIT, '(f0.16)')  x_mom_resnorm_0
      write(RESTART_FILE_UNIT, '(f0.16)')  y_mom_resnorm_0
      write(RESTART_FILE_UNIT, '(f0.16)')  z_mom_resnorm_0
      write(RESTART_FILE_UNIT, '(f0.16)') energy_resnorm_0
      write(RESTART_FILE_UNIT, '(f0.16)')    TKE_resnorm_0
      write(RESTART_FILE_UNIT, '(f0.16)')  omega_resnorm_0
    end subroutine write_initial_resnorm

end module dump_solution
