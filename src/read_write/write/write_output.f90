  !< Open/close and call other modules for writing solution
module write_output
  !< Open/close and call other modules for writing solution
  !< based on the input: type of file, either vtk or tecplot
  !< modules are called.
  use global                 ,only : OUT_FILE_UNIT
  use global_vars            ,only : outfile
  use global_vars            ,only : outfile
  use global_vars            ,only : write_data_format
  use global_vars            ,only : write_file_format
  use utils
  use string
  use write_output_vtk       ,only : write_file_vtk     => write_file
  use write_output_tec       ,only : write_file_tec => write_file
  use write_output_tec_node  ,only : write_file_tec_nodal => write_file

  implicit none
  private

  character(len=16) :: data_format
  character(len=16) :: file_format
  public write_file

  contains

    subroutine setup_file()
      !< Setup the file type based on the input
      implicit none
      call dmsg(1, 'write_output_vtk', 'setup_file')
      if (write_file_format == "vtk") then
        file_format = ".vtk"
      elseif (write_file_format == "tecplot" .or. write_file_format == "tecplot_nodal") then
        file_format = ".dat"
      else
        print*, "File format not recoganised. Accepted formats are"
        print*, "'vtk', 'tecplot' and 'tecplot_nodal' "
      end if

      if (write_data_format == "ASCII") then
        data_format = "formatted"
      elseif (write_data_format == "BINARY") then
        data_format = "unformatted"
      else
        print*, "Data format not recoganised. Accepted formats are"
        print*, "'ASCII' and 'BINARY' "
      end if

    end subroutine setup_file

    subroutine open_file(filename)
      !< open the file to write the solution
      implicit none
      character(len=*), intent(in) :: filename 
      call dmsg(1, 'write_output_vtk', 'open_file')

      open(OUT_FILE_UNIT, file=trim(filename)//trim(file_format) + '.part', form=trim(data_format))

    end subroutine open_file

    subroutine close_file(filename)
      !< close the file after writing solution.
      implicit none

      character(len=*), intent(in) :: filename 
      call dmsg(1, 'write_output_vtk', 'close_file')
      call rename(trim(filename)//trim(file_format) + '.part', trim(filename)//trim(file_format))
      close(OUT_FILE_UNIT)
    end subroutine close_file

    subroutine write_file()
      !< Writing output in the file according to the input file type
      implicit none

      call setup_file()
      call open_file(outfile)

      select case (write_file_format)

        case ('vtk')
          call write_file_vtk()

        case ('tecplot')
          call write_file_tec()

        case ('tecplot_nodal')
          call write_file_tec_nodal()

        case DEFAULT
          call dmsg(5, 'write_output', 'write_file',&
            'ERROR: write file format nor recognised. READ format -> '//write_file_format)

      end select

      call close_file(outfile)

    end subroutine write_file

end module write_output
