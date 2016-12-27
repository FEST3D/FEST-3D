module dump_solution
  !-------------------------------------------
  ! This module contians subroutine that
  !  1. check if point of dump is arrived.
  !  2. create particular folder for dump.
  !  3. dump data in that folder.
  !  4. purge folders if required.
  !------------------------------------------

  use global,      only : FILE_NAME_LENGTH
  use global_vars, only : current_iter
  use global_vars, only : checkpoint_iter
  use global_vars, only : checkpoint_iter_count
  use global_vars, only : purge_write
  use global_vars, only : sim_clock
  use utils
  use string
  use state,       only : writestate_vtk
  use layout,      only : process_id

  implicit none
  private
  character(len=FILE_NAME_LENGTH) :: dump_dirname
  character(len=FILE_NAME_LENGTH) :: purge_dirname

  public :: checkpoint

  contains

    subroutine checkpoint()
      !-----------------------------------------------------------
      ! Create a checkpoint dump file if the time has come
      !-----------------------------------------------------------

      implicit none


      if (checkpoint_iter .ne. 0) then
          if (mod(current_iter, checkpoint_iter) == 0) then
              call make_dump_dir()
              call dump_data()
              call purge_dump_dir()
              checkpoint_iter_count = checkpoint_iter_count + 1
              call dmsg(3, 'dump_solution', 'checkpoint', &
                      'Checkpoint created at iteration: ' + current_iter)

          end if
      end if

    end subroutine checkpoint

    subroutine create_directory(dirname)
      implicit none
      character(len=*), intent(in)    :: dirname
      character(len=FILE_NAME_LENGTH) :: mkdircmd

      mkdircmd = 'mkdir -p '//trim(dirname)
      call system(mkdircmd)

    end subroutine create_directory

    subroutine remove_directory(dirname)
      implicit none
      character(len=*), intent(in)    :: dirname
      character(len=FILE_NAME_LENGTH) :: rmdircmd

      rmdircmd = 'rm -rf '//trim(dirname)
      call system(rmdircmd)

    end subroutine remove_directory

    subroutine purge_dump_dir()
      implicit none
      integer                         :: purge_num

      purge_num = checkpoint_iter_count-purge_write
      if (purge_num > 0) then
        write(purge_dirname,'(A,I4.4)') 'time_directories/', purge_num
        call remove_directory(purge_dirname)
      end if

    end subroutine purge_dump_dir

    subroutine make_dump_dir()
      implicit none

      write(dump_dirname,'(A,I4.4)') 'time_directories/',checkpoint_iter_count
      call create_directory(dump_dirname)

    end subroutine make_dump_dir

    subroutine dump_data()
      implicit none
      character(len=FILE_NAME_LENGTH) :: filename

      write(filename, '(A,I2.2,A)') trim(dump_dirname)//'/process_',process_id,'.vtk'
      call writestate_vtk(filename, 'Simulaion clock: ' + sim_clock)

    end subroutine dump_data

end module dump_solution
