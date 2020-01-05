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
  use utils
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

    subroutine checkpoint(files, qp, nodes, control, scheme, dims)
      !< Create a checkpoint dump file if the time has come
      !-----------------------------------------------------------

      implicit none
      type(filetype), intent(inout) :: files
      type(extent), intent(in) :: dims
      type(controltype), intent(inout) :: control
      type(schemetype), intent(in) :: scheme
      type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(in) :: nodes
      real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp


      DebugCall('checkpoint')

      if (control%checkpoint_iter .ne. 0) then
          if (mod(control%current_iter, control%checkpoint_iter) == 0 &
             .or. control%current_iter == control%max_iters) then
              call make_dump_dir(control)
              call dump_data(files, qp, nodes, control, scheme, dims)
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

    subroutine dump_data(files, qp, nodes, control, scheme, dims)
      !< Call to write save files in the directory
      implicit none
      type(filetype), intent(inout) :: files
      type(extent), intent(in) :: dims
      type(controltype), intent(in) :: control
      type(schemetype), intent(in) :: scheme
      type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(in) :: nodes
      real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp

      DebugCall('dump_solution: dump_data')
      write(files%restartfile, '(A,I2.2)') trim(dump_dirname)//'/restart/process_',process_id
      write(files%outfile, '(A,I2.2)') trim(dump_dirname)//'/process_',process_id
      call write_restart_log(files, scheme, control)
      call write_file(files, qp, nodes, control, scheme, dims)

    end subroutine dump_data

    subroutine write_restart_log(files, scheme, control)
      !< Call to write log file in the subdirectory "restart". 
      !< It is useful information while restarting the solver
      implicit none
      type(filetype), intent(in) :: files
      type(controltype), intent(in) :: control
      type(schemetype), intent(in) :: scheme
      open(files%RESTART_FILE_UNIT, file=files%restartfile)
      select case (scheme%turbulence)
          
        case ('none')
          write(files%RESTART_FILE_UNIT, '(A)') 'viscous'
        case('sst','sst2003', 'kkl', 'ke', 'kw', 'sa', 'saBC', 'des-sst')
          write(files%RESTART_FILE_UNIT, '(A)') trim(scheme%turbulence)
        case DEFAULT
           Fatal_error
      end select
      call write_initial_resnorm(files, control)
      close(files%RESTART_FILE_UNIT)

    end subroutine write_restart_log

    subroutine write_initial_resnorm(files, control)
      !< Writing Initial resnorom in the log file to 
      !< maintian continuity of resnorm while restrarting
      implicit none
      type(filetype), intent(in) :: files
      type(controltype), intent(in) :: control
      integer :: i
      !integer, intent(in) :: current_iter, last_iter
      write(files%RESTART_FILE_UNIT, '(I0)')    control%current_iter+control%last_iter
      do i = 1,control%n_var+1
        write(files%RESTART_FILE_UNIT, '(f0.16)')  control%previous_res(i)
      end do
      !write(files%RESTART_FILE_UNIT, '(f0.16)')    vis_resnorm_0
      !write(files%RESTART_FILE_UNIT, '(f0.16)')   turb_resnorm_0
      !write(files%RESTART_FILE_UNIT, '(f0.16)')   cont_resnorm_0
      !write(files%RESTART_FILE_UNIT, '(f0.16)')  x_mom_resnorm_0
      !write(files%RESTART_FILE_UNIT, '(f0.16)')  y_mom_resnorm_0
      !write(files%RESTART_FILE_UNIT, '(f0.16)')  z_mom_resnorm_0
      !write(files%RESTART_FILE_UNIT, '(f0.16)') energy_resnorm_0
      !write(files%RESTART_FILE_UNIT, '(f0.16)')    TKE_resnorm_0
      !write(files%RESTART_FILE_UNIT, '(f0.16)')  omega_resnorm_0
    end subroutine write_initial_resnorm

end module dump_solution
