  !< Read the restart file in the vtk format
module read_output_vtk
  !< Read the restart file in the vtk format
  !---------------------------------------------------------
  ! This module read state + other variable in output file
  !---------------------------------------------------------
#include "../../debug.h"
#include "../../error.h"
  use vartypes
  use utils

  implicit none
  private
  integer :: IN_FILE_UNIT
  integer :: imx, jmx, kmx
  integer :: i,j,k
  real(wp), dimension(:, :, :), pointer :: density      
   !< Rho pointer, point to slice of qp (:,:,:,1)
  real(wp), dimension(:, :, :), pointer :: x_speed      
   !< U pointer, point to slice of qp (:,:,:,2) 
  real(wp), dimension(:, :, :), pointer :: y_speed      
   !< V pointer, point to slice of qp (:,:,:,3) 
  real(wp), dimension(:, :, :), pointer :: z_speed      
   !< W pointer, point to slice of qp (:,:,:,4)
  real(wp), dimension(:, :, :), pointer :: pressure     
   !< P pointer, point to slice of qp (:,:,:,5)
  real(wp), dimension(:, :, :), pointer :: tk        
  !< TKE, point to slice of qp (:,:,:,6)
  real(wp), dimension(:, :, :), pointer :: tw        
  !< Omega, point to slice of qp (:,:,:,7)
  real(wp), dimension(:, :, :), pointer :: te        
  !< Dissipation, point to slice of qp (:,:,:,7)
  real(wp), dimension(:, :, :), pointer :: tv        
  !< SA visocity, point to slice of qp (:,:,:,6)
  real(wp), dimension(:, :, :), pointer :: tkl       
  !< KL K-KL method, point to slice of qp (:,:,:,7)
  real(wp), dimension(:, :, :), pointer :: tgm       
  !< Intermittency of LCTM2015, point to slice of qp (:,:,:,8)

  public :: read_file

  contains

    subroutine read_file(file_handler, state, control, scheme, dims)
      !< Read all the variable for the vtk restart file
      implicit none
      integer, intent(in) :: file_handler
      type(controltype), intent(in) :: control
      type(schemetype), intent(in) :: scheme
      type(extent), intent(in) :: dims
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(inout), target :: state
      integer :: n

      IN_FILE_UNIT = file_handler

      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx

      density(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 1)
      x_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 2)
      y_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 3)
      z_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 4)
      pressure(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 5)

      select case (trim(scheme%turbulence))
          case ("none")
              !include nothing
              continue
          
          case ("sst", "sst2003", "bsl", "des-sst", "kw")
              tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 6)
              tw(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 7)

          case ("kkl")
              tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 6)
              tkl(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 7)

          case ("sa", "saBC")
              tv(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 6)

          case ("ke")
              tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 6)
              te(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 7)

          case DEFAULT
            Fatal_error
      end select

      ! Transition modeling
      select case(trim(scheme%transition))
        case('lctm2015')
          tgm(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 8)

        case('bc', 'none')
          !do nothing
          continue

        case DEFAULT
          Fatal_error
      end Select

      call read_header()
      call read_grid()
      do n = 1,control%r_count

        select case (trim(control%r_list(n)))
        
          case('Velocity')
            call read_velocity()

          case('Density')
            call read_scalar(density, 'Density', -2)
          
          case('Pressure')
            call read_scalar(pressure, 'Pressure', -2)
            
          case('TKE')
            call read_scalar(tk, 'TKE', -2)

          case('Omega')
            call read_scalar(tw, 'Omega', -2)

          case('Kl')
            call read_scalar(tkl, 'Kl', -2)

          case('tv')
            call read_scalar(tv, 'tv', -2)

          case('tgm')
            call read_scalar(tgm, 'tgm', -2)

          case('do not read')
            call skip_scalar()

          case Default
            Print*, "read error: list var : "//trim(control%r_list(n))

        end select
      end do

    end subroutine read_file

    subroutine read_header()
      !< Skip read the header in the vtk file
      implicit none

      DebugCall('read_output_vtk: read_header')
      read(IN_FILE_UNIT, *) !'# vtk DataFile Version 3.1'
      read(IN_FILE_UNIT, *) !'cfd-iitm output'   ! comment line
      read(IN_FILE_UNIT, *) !trim(read_data_format)
      read(IN_FILE_UNIT, *) !'DATASET STRUCTURED_GRID'
      !read(IN_FILE_UNIT, *)


    end subroutine read_header

    subroutine read_grid()
      !< Skip the grid read in the restart file
      implicit none

      ! read grid point coordinates
      DebugCall('read_output_vtk: read_grid')
      read(IN_FILE_UNIT, * ) !'DIMENSIONS ', imx, ' ', jmx, ' ', kmx
      read(IN_FILE_UNIT, * ) !'POINTS ', imx*jmx*kmx, ' DOUBLE'
      do k = 1, kmx
       do j = 1, jmx
        do i = 1, imx
          read(IN_FILE_UNIT, *) !point(i, j, k)%x, ' ', point(i, j, k)%y, ' ', point(i, j, k)%z
        end do
       end do
      end do
      read(IN_FILE_UNIT, *) 

    end subroutine read_grid

    subroutine read_velocity()
      !< Read velocity vector from the vtk file
      implicit none

      DebugCall('read_output_vtk: read_velocity')
      read(IN_FILE_UNIT, *) !'CELL_DATA ', (imx-1)*(jmx-1)*(kmx-1)
      read(IN_FILE_UNIT, *) !'VECTORS Velocity FLOAT'
      do k = 1, kmx - 1
       do j = 1, jmx - 1
        do i = 1, imx - 1
          read(IN_FILE_UNIT, *) x_speed(i, j, k), y_speed(i, j, k), z_speed(i, j, k)
        end do
       end do
      end do
      read(IN_FILE_UNIT, *) 

    end subroutine read_velocity

    subroutine read_scalar(var, name,  index)
      !< Read scalar from the vtk file
      implicit none
      integer, intent(in) :: index
      real(wp), dimension(index:imx-index,index:jmx-index,index:kmx-index), intent(out) :: var
      character(len=*), intent(in) :: name

      DebugCall('read_output_vtk'//trim(name))
      read(IN_FILE_UNIT, *) !'SCALARS '//trim(name)//' FLOAT'
      read(IN_FILE_UNIT, *) !'LOOKUP_TABLE default'
      do k = 1, kmx - 1
       do j = 1, jmx - 1
        do i = 1, imx - 1
          read(IN_FILE_UNIT, *) var(i, j, k)
        end do
       end do
      end do
      read(IN_FILE_UNIT, *)

    end subroutine read_scalar

    subroutine skip_scalar()
      !< Skip read scalar from the vtk file
      implicit none

      DebugCall('read_output_vtk: skip_scalar')
      read(IN_FILE_UNIT, *) !'SCALARS '//trim(name)//' FLOAT'
      read(IN_FILE_UNIT, *) !'LOOKUP_TABLE default'
      do k = 1, kmx - 1
       do j = 1, jmx - 1
        do i = 1, imx - 1
          read(IN_FILE_UNIT, *)
        end do
       end do
      end do
      read(IN_FILE_UNIT, *)

    end subroutine skip_scalar

end module read_output_vtk
