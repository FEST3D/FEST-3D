    !< The grid module read grid file and allocate memory to storing variables
module grid
    !< The grid module contains the grid definition (locations of the 
    !< grid points) as well as procedures to load these from a file.
    !-------------------------------------------------------------------
    use vartypes
    use mpi
    use mapping, only : read_interface_map
    
#include "error.h"
#include "debug.h"
    private

    ! Public methods
    public :: setup_grid

    contains

        subroutine setup_grid(files, nodes, control, bc, dims)
            !< Read the grid file and initialize the grid
            !-----------------------------------------------------------

            implicit none
            type(filetype), intent(in) :: files
            !< Files' name and handler
            type(controltype), intent(in) :: control
            !< Control parameters
            type(boundarytype), intent(inout) :: bc
            !< boundary conditions and fixed values
            type(nodetype), dimension(:,:,:), allocatable, intent(out) :: nodes
            !< Grid points 
            type(extent), intent(out) :: dims
            !< Extent of the domain:imx,jmx,kmx
            
            DebugCall('setup_grid')

            open(files%GRID_FILE_UNIT, file=files%gridfile)

            call extract_grid_size(files%GRID_FILE_UNIT, dims)

            ! allocate memory for storing grid points
            allocate(nodes(-2:dims%imx+3, -2:dims%jmx+3, -2:dims%kmx+3))

            !read interface mapping
            call read_interface_map(files, control, bc, dims)

            ! ghost grid exchange
            call populate_grid_points(files%GRID_FILE_UNIT, nodes, dims)

            close(files%GRID_FILE_UNIT)

            ! populate ghost grid points
            call ghost_grid(nodes, dims)

        
        end subroutine setup_grid
        
        subroutine extract_grid_size(file_handler, dims)
            !< Extract the grid size from the grid file header
            !
            ! We assume that the grid could be in 1 or 2 dimensions. If
            ! the grid is in 1 dimension, jmx will be set to 1.
            ! We assume that at least one number is specified in the 
            ! header, i.e., the grid has atleast one dimension.
            !-----------------------------------------------------------

            implicit none
            integer, intent(in) :: file_handler
            !< (input)file handling unit
            character(len=STRING_BUFFER_LENGTH) :: header
            !< store header
            type(extent), intent(out) :: dims
            !< Extent of the domain:imx,jmx,kmx
            integer :: ios  ! io operation status

            DebugCall('extract_grid_size')

            read(file_handler, '(A)', iostat=ios) header
            if (ios /= 0) then
                print *, 'Error while reading grid file header.'
                print *, 'Current buffer length is set to: ', &
                        STRING_BUFFER_LENGTH
                !stop
            end if

            ! Try to read constants corresponding to two dimensions.
            read(header, *, iostat=ios) dims%imx, dims%jmx, dims%kmx
            if (ios /= 0) then
              print*, "Not able to read dimension from the grid file"
              print*, "Make sure you provdie 3D grid"
              Fatal_error
            end if

        end subroutine extract_grid_size

        subroutine populate_grid_points(file_handler, nodes, dims)
            !< Use the grid file to populate the grid points.
            !-----------------------------------------------------------

            implicit none
            integer, intent(in) :: file_handler
            !< (input)file handling unit
            type(extent), intent(in) :: dims
            !< Extent of the domain:imx,jmx,kmx
            type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(out) :: nodes
            !< Grid points
            character(len=STRING_BUFFER_LENGTH) :: line
            !< store read line
            integer :: i, j, k
            integer :: ios  
            !< input/output  status

            DebugCall('populate_grid_point')
         !  print *, imx, jmx, kmx

            ! Read grid points from the grid file
            do k = 1, dims%kmx
                do j = 1, dims%jmx
                    do i = 1, dims%imx
                        read(file_handler, '(A)', iostat=ios) line
                        if (ios /= 0) then
                            print *, 'Error while reading grid line.'
                            print *, 'Current grid point: ', i, j, k
                            print *, 'Current buffer length is set to: ', &
                                     STRING_BUFFER_LENGTH
                            print *, 'Exiting program.'
                            !stop
                        end if
                        !call extract_grid_point(line, i, j, k)
                        read(line, *) nodes(i, j, k)%x, nodes(i, j, k)%y, nodes(i, j, k)%z
                    end do
                end do
            end do

        end subroutine populate_grid_points

        subroutine ghost_grid(nodes, dims)
          !< generate ghost grid for the various operations later.
          implicit none
          type(extent), intent(in) :: dims
          !< Extent of the domain:imx,jmx,kmx
          type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(inout) :: nodes
          !< grid points

          DebugCall('ghost_grid')
          !-------------------------------------------------------------------
          !getting ghost cell for all faces even if it is a interface cell
          ! <algorithm>
          ! Point_ghost = 2*Point_first_inner_cell - Point_second_inner_cell
          ! </algorithm>
          !-------------------------------------------------------------------

          !--- I faces ---!
          !imin face -> 0 grid point
          !imin face -> -1 grid point
          !imin face -> -2 grid point
          nodes( 0,:,:)%x=2*nodes( 1,:,:)%x-nodes(2,:,:)%x
          nodes( 0,:,:)%y=2*nodes( 1,:,:)%y-nodes(2,:,:)%y
          nodes( 0,:,:)%z=2*nodes( 1,:,:)%z-nodes(2,:,:)%z
          nodes(-1,:,:)%x=2*nodes( 0,:,:)%x-nodes(1,:,:)%x
          nodes(-1,:,:)%y=2*nodes( 0,:,:)%y-nodes(1,:,:)%y
          nodes(-1,:,:)%z=2*nodes( 0,:,:)%z-nodes(1,:,:)%z
          nodes(-2,:,:)%x=2*nodes(-1,:,:)%x-nodes(0,:,:)%x
          nodes(-2,:,:)%y=2*nodes(-1,:,:)%y-nodes(0,:,:)%y
          nodes(-2,:,:)%z=2*nodes(-1,:,:)%z-nodes(0,:,:)%z

          !imax face -> imx+1 grid point
          !imax face -> imx+2 grid point
          !imax face -> imx+3 grid point
          nodes(dims%imx+1,:,:)%x=2*nodes(dims%imx+0,:,:)%x-nodes(dims%imx-1,:,:)%x
          nodes(dims%imx+1,:,:)%y=2*nodes(dims%imx+0,:,:)%y-nodes(dims%imx-1,:,:)%y
          nodes(dims%imx+1,:,:)%z=2*nodes(dims%imx+0,:,:)%z-nodes(dims%imx-1,:,:)%z
          nodes(dims%imx+2,:,:)%x=2*nodes(dims%imx+1,:,:)%x-nodes(dims%imx-0,:,:)%x
          nodes(dims%imx+2,:,:)%y=2*nodes(dims%imx+1,:,:)%y-nodes(dims%imx-0,:,:)%y
          nodes(dims%imx+2,:,:)%z=2*nodes(dims%imx+1,:,:)%z-nodes(dims%imx-0,:,:)%z
          nodes(dims%imx+3,:,:)%x=2*nodes(dims%imx+2,:,:)%x-nodes(dims%imx+1,:,:)%x
          nodes(dims%imx+3,:,:)%y=2*nodes(dims%imx+2,:,:)%y-nodes(dims%imx+1,:,:)%y
          nodes(dims%imx+3,:,:)%z=2*nodes(dims%imx+2,:,:)%z-nodes(dims%imx+1,:,:)%z


          !--- Jmin faces ---!
          !jmin faces -> 0 grid point
          !jmin face -> -1 grid point
          !jmin face -> -2 grid point
          nodes(:, 0,:)%x=2*nodes(:, 1,:)%x-nodes(:,2,:)%x
          nodes(:, 0,:)%y=2*nodes(:, 1,:)%y-nodes(:,2,:)%y
          nodes(:, 0,:)%z=2*nodes(:, 1,:)%z-nodes(:,2,:)%z
          nodes(:,-1,:)%x=2*nodes(:, 0,:)%x-nodes(:,1,:)%x
          nodes(:,-1,:)%y=2*nodes(:, 0,:)%y-nodes(:,1,:)%y
          nodes(:,-1,:)%z=2*nodes(:, 0,:)%z-nodes(:,1,:)%z
          nodes(:,-2,:)%x=2*nodes(:,-1,:)%x-nodes(:,0,:)%x
          nodes(:,-2,:)%y=2*nodes(:,-1,:)%y-nodes(:,0,:)%y
          nodes(:,-2,:)%z=2*nodes(:,-1,:)%z-nodes(:,0,:)%z

          !jmax face -> jmx+1 grid point
          !jmax face -> jmx+3 grid point
          !jmax face -> jmx+2 grid point
          nodes(:,dims%jmx+1,:)%x=2*nodes(:,dims%jmx+0,:)%x-nodes(:,dims%jmx-1,:)%x
          nodes(:,dims%jmx+1,:)%y=2*nodes(:,dims%jmx+0,:)%y-nodes(:,dims%jmx-1,:)%y
          nodes(:,dims%jmx+1,:)%z=2*nodes(:,dims%jmx+0,:)%z-nodes(:,dims%jmx-1,:)%z
          nodes(:,dims%jmx+2,:)%x=2*nodes(:,dims%jmx+1,:)%x-nodes(:,dims%jmx-0,:)%x
          nodes(:,dims%jmx+2,:)%y=2*nodes(:,dims%jmx+1,:)%y-nodes(:,dims%jmx-0,:)%y
          nodes(:,dims%jmx+2,:)%z=2*nodes(:,dims%jmx+1,:)%z-nodes(:,dims%jmx-0,:)%z
          nodes(:,dims%jmx+3,:)%x=2*nodes(:,dims%jmx+2,:)%x-nodes(:,dims%jmx+1,:)%x
          nodes(:,dims%jmx+3,:)%y=2*nodes(:,dims%jmx+2,:)%y-nodes(:,dims%jmx+1,:)%y
          nodes(:,dims%jmx+3,:)%z=2*nodes(:,dims%jmx+2,:)%z-nodes(:,dims%jmx+1,:)%z


          !--- Kmax faces ---!
          !kmin faces -> 0 grid point
          !kmin face -> -1 grid point
          !kmin face -> -2 grid point
          nodes(:,:, 0)%x=2*nodes(:,:, 1)%x-nodes(:,:,2)%x
          nodes(:,:, 0)%y=2*nodes(:,:, 1)%y-nodes(:,:,2)%y
          nodes(:,:, 0)%z=2*nodes(:,:, 1)%z-nodes(:,:,2)%z
          nodes(:,:,-1)%x=2*nodes(:,:, 0)%x-nodes(:,:,1)%x
          nodes(:,:,-1)%y=2*nodes(:,:, 0)%y-nodes(:,:,1)%y
          nodes(:,:,-1)%z=2*nodes(:,:, 0)%z-nodes(:,:,1)%z
          nodes(:,:,-2)%x=2*nodes(:,:,-1)%x-nodes(:,:,0)%x
          nodes(:,:,-2)%y=2*nodes(:,:,-1)%y-nodes(:,:,0)%y
          nodes(:,:,-2)%z=2*nodes(:,:,-1)%z-nodes(:,:,0)%z

          !kmax face -> kmx+1 grid point
          !kmax face -> kmx+2 grid point
          !kmax face -> kmx+3 grid point
          nodes(:,:,dims%kmx+1)%x=2*nodes(:,:,dims%kmx+0)%x-nodes(:,:,dims%kmx-1)%x
          nodes(:,:,dims%kmx+1)%y=2*nodes(:,:,dims%kmx+0)%y-nodes(:,:,dims%kmx-1)%y
          nodes(:,:,dims%kmx+1)%z=2*nodes(:,:,dims%kmx+0)%z-nodes(:,:,dims%kmx-1)%z
          nodes(:,:,dims%kmx+2)%x=2*nodes(:,:,dims%kmx+1)%x-nodes(:,:,dims%kmx-0)%x
          nodes(:,:,dims%kmx+2)%y=2*nodes(:,:,dims%kmx+1)%y-nodes(:,:,dims%kmx-0)%y
          nodes(:,:,dims%kmx+2)%z=2*nodes(:,:,dims%kmx+1)%z-nodes(:,:,dims%kmx-0)%z
          nodes(:,:,dims%kmx+3)%x=2*nodes(:,:,dims%kmx+2)%x-nodes(:,:,dims%kmx+1)%x
          nodes(:,:,dims%kmx+3)%y=2*nodes(:,:,dims%kmx+2)%y-nodes(:,:,dims%kmx+1)%y
          nodes(:,:,dims%kmx+3)%z=2*nodes(:,:,dims%kmx+2)%z-nodes(:,:,dims%kmx+1)%z

        end subroutine ghost_grid

end module grid
