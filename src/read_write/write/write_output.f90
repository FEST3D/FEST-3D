  !< Open/close and call other modules for writing solution
module write_output
  !< Open/close and call other modules for writing solution
  !< based on the input: type of file, either vtk or tecplot
  !< modules are called
#include "../../debug.h"
#include "../../error.h"
  use vartypes
  use utils
  use write_output_vtk       ,only : write_file_vtk     => write_file
  use write_output_tec       ,only : write_file_tec => write_file
  use write_output_tec_node  ,only : write_file_tec_nodal => write_file

  implicit none
  private

  character(len=16) :: data_format
  character(len=16) :: file_format
  public write_file

  contains

    subroutine setup_file(control)
      !< Setup the file type based on the input
      implicit none
      type(controltype), intent(in) :: control
      DebugCall('write_output_vtk: setup_file')
      if (control%write_file_format == "vtk") then
        file_format = ".vtk"
      elseif (control%write_file_format == "tecplot" .or. control%write_file_format == "tecplot_nodal") then
        file_format = ".dat"
      else
        print*, "File format not recoganised. Accepted formats are"
        print*, "'vtk', 'tecplot' and 'tecplot_nodal' "
      end if

      if (control%write_data_format == "ASCII") then
        data_format = "formatted"
      elseif (control%write_data_format == "BINARY") then
        data_format = "unformatted"
      else
        print*, "Data format not recoganised. Accepted formats are"
        print*, "'ASCII' and 'BINARY' "
      end if

    end subroutine setup_file

    subroutine open_file(file_handler, filename)
      !< Open the file to write the solution
      implicit none
      integer, intent(in) :: file_handler
      character(len=*), intent(in) :: filename 
      DebugCall('write_output: open_file')

      open(file_handler, file=trim(filename)//trim(file_format)//'.part', form=trim(data_format))

    end subroutine open_file

    subroutine close_file(file_handler, filename)
      !< Close the file after writing solution.
      implicit none
      integer, intent(in) :: file_handler
      character(len=*), intent(in) :: filename 
      DebugCall('write_output_vtk: close_file')
      call rename(trim(filename)//trim(file_format)//'.part', trim(filename)//trim(file_format))
      close(file_handler)
    end subroutine close_file

    subroutine write_file(files, qp, nodes, control, scheme,  dims)
      !< Writing output in the file according to the input file type
      implicit none
      type(filetype), intent(in) :: files
      type(extent), intent(in) :: dims
      type(controltype), intent(in) :: control
      type(schemetype), intent(in) :: scheme
      type(nodetype), dimension(-2:dims%imx+3, -2:dims%jmx+3,-2:dims%kmx+3), intent(in) :: nodes
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2,-2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      integer:: file_handler

      file_handler = files%OUT_FILE_UNIT

      call setup_file(control)
      call open_file(file_handler, files%outfile)

      select case (control%write_file_format)

        case ('vtk')
          call write_file_vtk(file_handler, qp, nodes, control, scheme, dims)

        case ('tecplot')
          call write_file_tec(file_handler, qp, nodes, control, scheme, dims)

        case ('tecplot_nodal')
          call write_file_tec_nodal(file_handler, qp, nodes, control, scheme, dims)

        case DEFAULT
          Fatal_error

      end select

      call close_file(file_handler, files%outfile)

    end subroutine write_file

end module write_output
