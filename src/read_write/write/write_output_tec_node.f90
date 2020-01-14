  !< Writing solution in the output file in tecplot format with
  !< node data instead of cell-center data.
module write_output_tec_node
  !< Writing solution in the output file in tecplot format with
  !< node data instead of cell-center data.
  !---------------------------------------------------------
  ! This module write state + other variable in output file
  !---------------------------------------------------------
#include "../../debug.h"
#include "../../error.h"
  use vartypes
  use viscosity, only : mu 
  use viscosity, only : mu_t 
  use wall_dist, only : dist
  use global_sst , only : sst_F1
  use gradients, only : gradu_x
  use gradients, only : gradu_y
  use gradients, only : gradu_z
  use gradients, only : gradv_x
  use gradients, only : gradv_y
  use gradients, only : gradv_z
  use gradients, only : gradw_x
  use gradients, only : gradw_y
  use gradients, only : gradw_z
  use gradients, only : gradT_x
  use gradients, only : gradT_y
  use gradients, only : gradT_z
  use gradients, only : gradtk_x
  use gradients, only : gradtk_y
  use gradients, only : gradtk_z
  use gradients, only : gradtw_x
  use gradients, only : gradtw_y
  use gradients, only : gradtw_z

  use utils

  implicit none
  private
  integer :: OUT_FILE_UNIT
  integer :: i,j,k
  character(len=*), parameter :: format="(35e25.15)"
  integer :: imx, jmx, kmx
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
  public :: write_file

  contains

    subroutine write_file(file_handler, state, nodes, control, scheme, dims)
      !< Write output file in the tecplot format with node data
      implicit none
      integer, intent(in) :: file_handler
      type(controltype), intent(in) :: control
      type(schemetype), intent(in) :: scheme
      type(extent), intent(in) :: dims
      type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(in) :: nodes 
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in), target :: state
      integer :: n
      character(len=*), parameter :: err="Write error: Asked to write non-existing variable- "

      DebugCall("write_file")
      
      OUT_FILE_UNIT = file_handler

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

      call write_header(control)
      call write_grid(nodes)

      do n = 1,control%w_count

        select case (trim(control%w_list(n)))
        
          case('Velocity')
            call write_scalar(x_speed, "u", -2)
            call write_scalar(y_speed, "v", -2)
            call write_scalar(z_speed, "w", -2)

          case('Density')
            call write_scalar(density, "Density", -2)
          
          case('Pressure')
            call write_scalar(pressure, "Pressure", -2)
            
          case('Mu')
            call write_scalar(mu, "Mu", -2)
            
          case('Mu_t')
            call write_scalar(mu_t, "Mu_t", -2)
            
          case('TKE')
            call write_scalar(tk, "TKE",  -2)

          case('Omega')
            call write_scalar(tw, "Omega", -2)

          case('Kl')
            call write_scalar(tkl, "Kl", -2)

          case('tv')
            call write_scalar(tv, "tv", -2)

         ! case('Wall_distance')
         !   call write_scalar(dist, "Wall_dist", 1)

          case('F1')
            call write_scalar(sst_F1, "F1",  -2)

          case('Dudx')
            call write_scalar(gradu_x ,"dudx ", 0)
                                               
          case('Dudy')                         
            call write_scalar(gradu_y ,"dudy ", 0)
                                               
          case('Dudz')                         
            call write_scalar(gradu_z ,"dudz ", 0)
                                               
          case('Dvdx')                         
            call write_scalar(gradv_x ,"dvdx ", 0)
                                               
          case('Dvdy')                         
            call write_scalar(gradv_y ,"dvdy ", 0)
                                               
          case('Dvdz')                         
            call write_scalar(gradv_z ,"dvdz ", 0)
                                               
          case('Dwdx')                         
            call write_scalar(gradw_x ,"dwdx ", 0)
                                               
          case('Dwdy')                         
            call write_scalar(gradw_y ,"dwdy ", 0)
                                               
          case('Dwdz')                         
            call write_scalar(gradw_z ,"dwdz ", 0)
                                               
          case('DTdx')                         
            call write_scalar(gradT_x ,"dTdx ", 0)
                                               
          case('DTdy')                         
            call write_scalar(gradT_y ,"dTdy ", 0)
                                               
          case('DTdz')                         
            call write_scalar(gradT_z ,"dTdz ", 0)
                                               
          case('Dtkdx')                        
            call write_scalar(gradtk_x,"dtkdx", 0)
                                               
          case('Dtkdy')                        
            call write_scalar(gradtk_y,"dtkdy", 0)
                                               
          case('Dtkdz')                        
            call write_scalar(gradtk_z,"dtkdz", 0)
                                               
          case('Dtwdx')                        
            call write_scalar(gradtw_x,"dtwdx", 0)
                                               
          case('Dtwdy')                        
            call write_scalar(gradtw_y,"dtwdy", 0)
                                               
          case('Dtwdz')                        
            call write_scalar(gradtw_z,"dtwdz", 0)

          case('do not write')
            ! do not write
            continue

          case Default
            print*, err//trim(control%w_list(n))//" to file"

        end select
      end do


    end subroutine write_file


    subroutine write_header(control)
      !< Write the header in the output file
      implicit none
      type(controltype), intent(in) :: control
      integer :: n
      integer :: total

      DebugCall('write_output_vtk: write_header')
      write(OUT_FILE_UNIT,'(a)') "variables = x y z "

      total=3
      do n = 1,control%w_count

        select case (trim(control%w_list(n)))
        
          case('Velocity')
            write(OUT_FILE_UNIT, '(a)') " u v w "
            total = total+3

          case('do not write')
            !skip 
            continue

          case Default
            write(OUT_FILE_UNIT, '(a)') trim(control%w_list(n))//" "
            total = total+1

        end select
      end do

      write(OUT_FILE_UNIT, '(a,i4.4,3(a,i5.5),a)') "zone T=block",process_id,"  i=",imx," j=",jmx, " k=",kmx-1, " Datapacking=Block"

      write(OUT_FILE_UNIT,*) "Varlocation=([1-3]=Nodal)"
      write(OUT_FILE_UNIT,'(a,i2.2,a)') "Varlocation=([4-",total,"]=Nodal)"
      write(OUT_FILE_UNIT,"(a,i4.4)") "STRANDID=",1
      write(OUT_FILE_UNIT,"(a,i4.4)") "SOLUTIONTIME=",control%checkpoint_iter_count


    end subroutine write_header

    subroutine write_grid(nodes)
      !< Write grid information in the output file
      implicit none
      type(nodetype), dimension(-2:imx+3,-2:jmx+3,-2:kmx+3), intent(in) :: nodes 

      ! write grid point coordinates
      DebugCall('write_output_tec_node: write_grid')
      write(OUT_FILE_UNIT, format) (((nodes(i, j, k)%x,i=1,imx), j=1,jmx), k=1,kmx-1)
      write(OUT_FILE_UNIT, format) (((nodes(i, j, k)%y,i=1,imx), j=1,jmx), k=1,kmx-1)
      write(OUT_FILE_UNIT, format) (((nodes(i, j, k)%z,i=1,imx), j=1,jmx), k=1,kmx-1)

    end subroutine write_grid

    subroutine write_scalar(var, name, index)
      !< Write scalar variable in the output file
      implicit none
      integer, intent(in) :: index
      real(wp), dimension(index:imx-index,index:jmx-index,index:kmx-index), intent(in) :: var
      character(len=*),       intent(in):: name

      DebugCall('write_scalar'//trim(name))

      write(OUT_FILE_UNIT, format) (((0.25*(var(i, j, k) + var(i-1,j,k)&
                                          + var(i, j-1, k) + var(i-1,j-1,k))&
                                          ,i=1,imx), j=1,jmx), k=1,kmx-1)

    end subroutine write_scalar


end module write_output_tec_node
