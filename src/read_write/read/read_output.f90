 !< Read output files from the restart folder
module read_output
 !< Read output files from the restart folder
  
  !---------------------------------------------------------
  ! This module read state + other variable in output file
  !---------------------------------------------------------
  use global     , only :      IN_FILE_UNIT
  use global     , only : RESTART_FILE_UNIT
  use global_vars, only :      infile
  use global_vars, only : restartfile

  use global_vars, only : read_data_format
  use global_vars, only : read_file_format
  use global_vars, only : start_from
  use global_vars, only : process_id
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
  use global_vars, only : previous_flow_type
  use global_vars, only : last_iter
  use global_vars, only : mu_ref

  use read_output_vtk, only : read_file_vtk => read_file
  use read_output_tec, only : read_file_tec => read_file
  use check_output_control, only: verify_read_control

  use utils
  use string

  implicit none
  private
  integer :: i,j,k
  real    :: speed_inf
  !< Free-stream velocity magnitude
  character(len=8) :: file_format
  !< Read file format
  character(len=16) :: data_format
  !< Read file data type
  character(len=16) :: read_flow_type
  !< Previous flow type 

  public :: read_file

  contains

    subroutine read_file()
      !< Read restart file
      implicit none
      call setup_file
      call open_file(infile)
      call read_restart_file()
      call verify_read_control()
        
      select case (read_file_format)
        
        case ('vtk')
          call read_file_vtk()
        
        case ('tecplot')
          call read_file_tec()
        
        case DEFAULT
        call dmsg(5, 'read_output', 'read_file',&
          'ERROR: read file format not recognised. READ format -> '//read_file_format)
      end select

      call close_file()
    end subroutine read_file


    subroutine setup_file()
      !< Steup the file to read the restart state.
      implicit none
      call dmsg(1, 'read_output_vtk', 'setup_file')
      if (read_file_format == "vtk") then
        file_format = ".vtk"
      elseif (read_file_format == "tecplot") then
        file_format = ".dat"
      else
        print*, "File format not recoganised. Accepted formats are"
        print*, "'vtk' and 'tecplot' "
      end if

      if (read_data_format == "ASCII") then
        data_format = "formatted"
      elseif (read_data_format == "BINARY") then
        data_format = "unformatted"
      else
        print*, "Data format not recoganised. Accepted formats are"
        print*, "'ASCII' and 'BINARY' "
      end if

      !write(infile,'(a,i4.4,a,i2.2)') 'time_directories/',start_from,'process_',process_id

    end subroutine setup_file

    subroutine open_file(filename)
      !< Open file from the restart folder 
      implicit none
      character(len=*), intent(in) :: filename 
      call dmsg(1, 'read_output_vtk', 'open_file')

      write(restartfile, '(A,I4.4,A,I2.2)') 'time_directories/',start_from,&
                          '/restart/process_', process_id
      open(IN_FILE_UNIT, file=trim(filename)//trim(file_format))!, form=trim(data_format))
      open(RESTART_FILE_UNIT, file=restartfile, status='old')

    end subroutine open_file

    subroutine close_file()
      !< Close the file after reading 
      implicit none

      call dmsg(1, 'read_output_vtk', 'close_files')
      close(IN_FILE_UNIT)
      close(RESTART_FILE_UNIT)

    end subroutine close_file

    subroutine read_restart_file()
      !< Read the sub-directory log file in the restart folder
      implicit none
      read(RESTART_FILE_UNIT, *) previous_flow_type

      read(RESTART_FILE_UNIT, *)        last_iter
      read(RESTART_FILE_UNIT, *)        resnorm_0
      read(RESTART_FILE_UNIT, *)    vis_resnorm_0
      read(RESTART_FILE_UNIT, *)   turb_resnorm_0
      read(RESTART_FILE_UNIT, *)   cont_resnorm_0
      read(RESTART_FILE_UNIT, *)  x_mom_resnorm_0
      read(RESTART_FILE_UNIT, *)  y_mom_resnorm_0
      read(RESTART_FILE_UNIT, *)  z_mom_resnorm_0
      read(RESTART_FILE_UNIT, *) energy_resnorm_0
      read(RESTART_FILE_UNIT, *)    TKE_resnorm_0
      read(RESTART_FILE_UNIT, *)  omega_resnorm_0
    end subroutine read_restart_file

end module read_output
