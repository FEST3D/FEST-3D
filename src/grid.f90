    !< The grid module read grid file and allocate memory to storing variables
module grid
    !< The grid module contains the grid definition (locations of the 
    !< grid points) as well as procedures to load these from a file.
    !-------------------------------------------------------------------
    
    use vartypes
    use global, only: STRING_BUFFER_LENGTH, GRID_FILE_UNIT
!    use global_vars, only : imx
!    use global_vars, only : jmx
!    use global_vars, only : kmx
    use global_vars, only : imin_id
    use global_vars, only : jmin_id
    use global_vars, only : kmin_id
    use global_vars, only : imax_id
    use global_vars, only : jmax_id
    use global_vars, only : kmax_id
    use mapping, only : Gilo
    use mapping, only : Gjlo
    use mapping, only : Gklo
    use mapping, only : Gihi
    use mapping, only : Gjhi
    use mapping, only : Gkhi
    use mapping, only : mpi_class
    use mapping, only : read_interface_map
    use global_vars, only : dir_switch
    use global_vars, only : layers
    use global_vars, only : process_id
    use utils, only: alloc, dealloc
    
#include "error.inc"
#include "mpi.inc"
#include "debug.h"
    private

!    type, public :: location
!      real :: x
!      real :: y
!      real :: z
!    end type location
!
!    type, public :: extent
!      integer :: imx
!      integer :: jmx
!      integer :: kmx
!    end type extent
!
!    public operator(+)
!    interface operator(+)
!      procedure pointadd
!    end interface operator(+)
!
!    public operator(-)
!    interface operator(-)
!      procedure pointsub
!    end interface operator(-)
!
!    public operator(*)
!    interface operator(*)
!      procedure pointIntscalar,&
!                 pointRealscalar
!    end interface operator(*)
!
!    public assignment(=)
!    interface assignment(=)
!      procedure pointassign
!    end interface assignment(=)
!
!    type(location), dimension(:,:,:), allocatable ::  point
!    
    ! Public methods
    public :: setup_grid
    public :: destroy_grid
    !public :: point

    contains

!        subroutine pointassign(p1,p2)
!          implicit none
!          type(nodetype), intent(out) :: p1
!          type(nodetype), intent(in) :: p2
!          p1%x = p2%x
!          p1%y = p2%y
!          p1%z = p2%z
!        end subroutine pointassign
!
!        function pointadd(p2, p3) result(p1)
!          implicit none
!          type(nodetype) :: p1
!          type(nodetype), intent(in) :: p2,p3
!          p1%x = p2%x + p3%x
!          p1%y = p2%y + p3%y
!          p1%z = p2%z + p3%z
!        end function pointadd
!    
!        function pointsub(p2, p3) result(p1)
!          implicit none
!          type(nodetype) :: p1
!          type(nodetype), intent(in) :: p2,p3
!          p1%x = p2%x - p3%x
!          p1%y = p2%y - p3%y
!          p1%z = p2%z - p3%z
!        end function pointsub
!    
!        function pointIntscalar(p2,p3) result(p1)
!          implicit none
!          integer, intent(in)    :: p2
!          type(nodetype), intent(in) :: p3
!          type(nodetype) :: p1
!          p1%x = p2*p3%x
!          p1%y = p2*p3%y
!          p1%z = p2*p3%z
!        end function pointIntscalar
!
!        function pointRealscalar(p2,p3) result(p1)
!          implicit none
!          type(nodetype)  :: p1
!          real , intent(in) :: p2
!          type(nodetype), intent(in) :: p3
!          p1%x = p2*p3%x
!          p1%y = p2*p3%y
!          p1%z = p2*p3%z
!        end function pointRealscalar
    
    
        subroutine allocate_memory(nodes, dims)
            !< Allocate memory to store the grid
            !-----------------------------------------------------------

            implicit none
            type(nodetype), dimension(:,:,:), allocatable, intent(out) :: nodes
            type(extent), intent(in) :: dims

            DebugCall("allocate memory to grid")
            allocate(nodes(-2:dims%imx+3, -2:dims%jmx+3, -2:dims%kmx+3))
            !allocate(point(-2:dims%imx+3, -2:dims%jmx+3, -2:dims%kmx+3))
            
        end subroutine allocate_memory

        subroutine destroy_grid()
            !< Deallocate the memory allocated for the grid.
            !-----------------------------------------------------------

            implicit none

            DebugCall('destroy_memory')

            !call dealloc(grid_x)
            !call dealloc(grid_y)
            !call dealloc(grid_z)

!           call dealloc(sphere_indices)

        end subroutine destroy_grid

        subroutine setup_grid(gridfile, mapfile, periodicfile, nodes, dims)
            !< Read the grid file and initialize the grid
            !-----------------------------------------------------------

            implicit none
            character(len=*), intent(in) :: gridfile
            character(len=*), intent(in) :: mapfile
            character(len=*), intent(in) :: periodicfile
            type(nodetype), dimension(:,:,:), allocatable, intent(out) :: nodes
            type(extent), intent(out) :: dims
            
            DebugCall('setup_grid')

            open(GRID_FILE_UNIT, file=gridfile)
            call extract_grid_size(dims)
            call allocate_memory(nodes, dims)
            !read interface mapping
            call read_interface_map(mapfile, periodicfile, dims)

            ! ghost grid exchange
            call populate_grid_points(nodes, dims)

            close(GRID_FILE_UNIT)

            ! populate ghost grid points
            call ghost_grid(nodes, dims)

            !point%x = nodes%x
            !point%y = nodes%y
            !point%z = nodes%z
        
        end subroutine setup_grid
        
        subroutine extract_grid_size(dims)
            !< Extract the grid size from the grid file header
            !
            ! We assume that the grid could be in 1 or 2 dimensions. If
            ! the grid is in 1 dimension, jmx will be set to 1.
            ! We assume that at least one number is specified in the 
            ! header, i.e., the grid has atleast one dimension.
            !-----------------------------------------------------------

            implicit none
            character(len=STRING_BUFFER_LENGTH) :: header
            type(extent), intent(out) :: dims
            integer :: ios  ! io operation status

            DebugCall('extract_grid_size')

            read(GRID_FILE_UNIT, '(A)', iostat=ios) header
            if (ios /= 0) then
                print *, 'Error while reading grid file header.'
                print *, 'Current buffer length is set to: ', &
                        STRING_BUFFER_LENGTH
                !stop
            end if

            ! Try to read constants corresponding to two dimensions.
            read(header, *, iostat=ios) dims%imx, dims%jmx, dims%kmx
           ! if (ios /= 0) then
           !     ! An io error means it was not possible to read kmx
           !     ! This means the file does not have kmx  and so, set
           !     ! the extent of this direction to 1. Read the remaining
           !     ! dimension from the header.
           !     read(header, *, iostat=ios) dims%imx, dims%jmx
           !     if (ios /= 0) then
           !         ! This means that jmx does not exist. Repeat again
           !         read(header, *, iostat=ios) dims%imx
           !         if (ios /= 0) then
           !             ! Error while reading
           !             print *, 'Unable to read grid extent.'
           !             stop
           !         end if
           !         dims%jmx = 1
           !     end if
           !     dims%kmx = 1
           ! end if
           ! print*, dims%imx, dims%jmx, dims%kmx
           ! imx = dims%imx
           ! jmx = dims%jmx
           ! kmx = dims%kmx

        end subroutine extract_grid_size
!
!
!        subroutine extract_grid_point(line, i, j, k)
!            !< Extract a grid point from a line of the grid file. 
!            !-----------------------------------------------------------
!
!            implicit none
!            character(len=STRING_BUFFER_LENGTH), intent(in) :: line
!            integer, intent(in) :: i, j, k
!
!            DebugCall('extract_grid_point')
!
!            if (kmx > 1) then
!                read(line, *) nodes(i, j, k)%x, nodes(i, j, k)%y, nodes(i, j, k)%z
!            else    
!                if (jmx > 1) then
!                    read(line, *) nodes(i, j, k)%x, nodes(i, j, k)%y
!                else
!                    read(line, *) nodes(i, j, k)%x
!                    nodes(i, j, k)%y = 0.
!                end if
!                nodes(i, j, k)%z = 0.
!            end if
!        end subroutine extract_grid_point

        subroutine populate_grid_points(nodes, dims)
            !< Use the grid file to populate the grid points.
            !-----------------------------------------------------------

            implicit none
            type(extent), intent(in) :: dims
            type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(out) :: nodes
            character(len=STRING_BUFFER_LENGTH) :: line
            integer :: i, j, k
            integer :: ios  ! io status

            DebugCall('populate_grid_point')
         !  print *, imx, jmx, kmx

            ! Read grid points from the grid file
            do k = 1, dims%kmx
                do j = 1, dims%jmx
                    do i = 1, dims%imx
                        read(GRID_FILE_UNIT, '(A)', iostat=ios) line
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
          type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(inout) :: nodes

          DebugCall('ghost_grid')
!          !-------------------------------------------------------------------
!          !getting ghost cell for all faces even if it is a interface cell
!          ! <algorithm>
!          ! Point_ghost = 2*Point_first_inner_cell - Point_second_inner_cell
!          ! </algorithm>
!          !-------------------------------------------------------------------
!
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
