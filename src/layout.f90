  !< Contains routine to load layout file and sets the layout variables
  !< and gets process id and total process 
module layout
  !< Contains routine to load layout file and sets the layout variables
  !< and gets process id and total process 
  !------------------------------
  use vartypes
#include "error.h"
#include "debug.h"
#include "mpi.inc"
  
  public :: read_layout_file
  public :: get_process_data

contains


  subroutine get_process_data(control)
    !<Get Processor Id and total number of processors
  implicit none
  type(controltype), intent(inout) :: control
  !< Control parameters
    ! Finds and sets process data
    integer :: ierr
    call MPI_COMM_RANK(MPI_COMM_WORLD,control%process_id,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,control%total_process,ierr)
    process_id = control%process_id

  end subroutine get_process_data

  subroutine get_next_token_parallel(handler, buf)
    !< Extract the next token from the layout file
    !<
    !< Each token is on a separate line.
    !< There may be multiple comments (lines beginning with #) 
    !< and blank lines in between.
    !< The purpose of this subroutine is to ignore all these 
    !< lines and return the next "useful" line.
    !-----------------------------------------------------------

    implicit none
    integer, intent(in) :: handler
    character(len=STRING_BUFFER_LENGTH), intent(out) :: buf
    integer :: ios

    do
       !read(CONFIG_FILE_UNIT, '(A)', iostat=ios) buf
       read(handler, '(A)', iostat=ios) buf
       if (ios /= 0) then
          print *, 'Error while reading config file.'
          print *, 'Current buffer length is set to: ', &
               STRING_BUFFER_LENGTH
          stop
       end if
       if (index(buf, '#') == 1) then
          ! The current line begins with a hash
          ! Ignore it
          continue
       else if (len_trim(buf) == 0) then
          ! The current line is empty
          ! Ignore it
          continue
       else
          ! A new token has been found
          ! Break out
          exit
       end if
    end do

  end subroutine get_next_token_parallel


  subroutine read_layout_file(files,control, bc)
    !< Read the layout file for particular processor
    implicit none
    type(filetype), intent(inout) :: files
    !< Files' name and handler
    type(boundarytype), intent(inout) :: bc
    !< boundary conditions and fixed values
    character(len=STRING_BUFFER_LENGTH) :: buf
    !< read buffer
    character(len=128) :: grid_file_buf
    !< Name of the gridfile to load
    character(len=128) :: bc_file
    !< Name of the boundary condition file to load.
    type(controltype) ,intent(in)::control
    !< Processor id for current block
    integer :: total_entries      
    !< Total enteries in layout.md for each processorS
    integer :: i,buf_id 
    DebugCall('read_layout_file')

    open(files%CONFIG_FILE_UNIT, file=files%layout_file)

    ! Read the parameters from the file
    call get_next_token_parallel(files%CONFIG_FILE_UNIT, buf)
    read(buf,*)!control%total_process
    call get_next_token_parallel(files%CONFIG_FILE_UNIT, buf)
    read(buf,*)total_entries
    i = 0
    !print *, process_id
    call get_next_token_parallel(files%CONFIG_FILE_UNIT, buf)
    do while(i < control%process_id)
          call get_next_token_parallel(files%CONFIG_FILE_UNIT, buf)
       i = i+1
    end do
    read(buf,*) buf_id, grid_file_buf, bc_file, bc%imin_id, bc%imax_id, bc%jmin_id,bc%jmax_id,bc%kmin_id,bc%kmax_id
    write(files%gridfile, '(A)') 'system/mesh/gridfiles/'//trim(grid_file_buf)
    write(files%bcfile, '(A)') 'system/mesh/bc/'//trim(bc_file)
  end subroutine read_layout_file


end module layout
