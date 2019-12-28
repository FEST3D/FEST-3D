 !< Read output files from the restart folder
module read_output
 !< Read output files from the restart folder
  
  !---------------------------------------------------------
  ! This module read state + other variable in output file
  !---------------------------------------------------------
#include "../../debug.h"
#include "../../error.h"
  use vartypes
  use read_output_vtk, only : read_file_vtk => read_file
  use read_output_tec, only : read_file_tec => read_file
  use check_output_control, only: verify_read_control

  use utils

  implicit none
  private
  !< Free-stream velocity magnitude
  character(len=8) :: file_format
  !< Read file format
  character(len=16) :: data_format
  !< Read file data type

  public :: read_file

  contains

    subroutine read_file(files, qp, control, scheme, dims)
      !< Read restart file
      implicit none
      type(filetype), intent(inout) :: files
      type(extent), intent(in) :: dims
      type(controltype), intent(inout) :: control
      type(schemetype) , intent(in) :: scheme
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(inout), target :: qp
      call setup_file(control)
      call open_file(files,  control)
      call read_restart_file(files%RESTART_FILE_UNIT, control)
      call verify_read_control(control, scheme)
        
      select case (control%read_file_format)
        
        case ('vtk')
          call read_file_vtk(files%IN_FILE_UNIT, qp, control, scheme, dims)
        
        case ('tecplot')
          call read_file_tec(files%IN_FILE_UNIT, qp, control, scheme, dims)
        
        case DEFAULT
          Fatal_error
      end select

      call close_file(files)
    end subroutine read_file


    subroutine setup_file(control)
      !< Steup the file to read the restart state.
      implicit none
      type(controltype), intent(in) :: control
      DebugCall('setup_file')
      if (control%read_file_format == "vtk") then
        file_format = ".vtk"
      elseif (control%read_file_format == "tecplot") then
        file_format = ".dat"
      else
        print*, "File format not recoganised. Accepted formats are"
        print*, "'vtk' and 'tecplot' "
      end if

      if (control%read_data_format == "ASCII") then
        data_format = "formatted"
      elseif (control%read_data_format == "BINARY") then
        data_format = "unformatted"
      else
        print*, "Data format not recoganised. Accepted formats are"
        print*, "'ASCII' and 'BINARY' "
      end if

      !write(infile,'(a,i4.4,a,i2.2)') 'time_directories/',start_from,'process_',process_id

    end subroutine setup_file

    subroutine open_file(files, control)
      !< Open file from the restart folder 
      implicit none
      type(filetype), intent(inout) :: files
      type(controltype), intent(in) :: control
      DebugCall('open_file')

      write(files%restartfile, '(A,I4.4,A,I2.2)') 'time_directories/',control%start_from,&
                          '/restart/process_', process_id
      open(files%IN_FILE_UNIT, file=trim(files%infile)//trim(file_format))!, form=trim(data_format))
      open(files%RESTART_FILE_UNIT, file=files%restartfile, status='old')

    end subroutine open_file

    subroutine close_file(files)
      !< Close the file after reading 
      implicit none
      type(filetype), intent(in) :: files

      DebugCall('close_files')
      close(files%IN_FILE_UNIT)
      close(files%RESTART_FILE_UNIT)

    end subroutine close_file

    subroutine read_restart_file(RESTART_FILE_UNIT, control)
      !< Read the sub-directory log file in the restart folder
      implicit none
      integer, intent(in) :: RESTART_FILE_UNIT
      type(controltype), intent(inout) :: control
      integer :: i
      read(RESTART_FILE_UNIT, *) control%previous_flow_type

      read(RESTART_FILE_UNIT, *)        control%last_iter
      do i = 1,control%n_var+1
        read(RESTART_FILE_UNIT, *)  control%previous_res(i)
      end do
      !read(RESTART_FILE_UNIT, *)    vis_resnorm_0
      !read(RESTART_FILE_UNIT, *)   turb_resnorm_0
      !read(RESTART_FILE_UNIT, *)   cont_resnorm_0
      !read(RESTART_FILE_UNIT, *)  x_mom_resnorm_0
      !read(RESTART_FILE_UNIT, *)  y_mom_resnorm_0
      !read(RESTART_FILE_UNIT, *)  z_mom_resnorm_0
      !read(RESTART_FILE_UNIT, *) energy_resnorm_0
      !read(RESTART_FILE_UNIT, *)    TKE_resnorm_0
      !read(RESTART_FILE_UNIT, *)  omega_resnorm_0
    end subroutine read_restart_file

end module read_output
