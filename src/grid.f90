    !< The grid module read grid file and allocate memory to storing variables
module grid
    !< The grid module contains the grid definition (locations of the 
    !< grid points) as well as procedures to load these from a file.
    !-------------------------------------------------------------------
    
    use global, only: STRING_BUFFER_LENGTH, GRID_FILE_UNIT
    use global_vars, only : imx
    use global_vars, only : jmx
    use global_vars, only : kmx
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
    use utils, only: alloc, dealloc, dmsg
    
#include "error.inc"
#include "mpi.inc"
#include "debug.h"
    private

    type :: position
      real :: x
      real :: y
      real :: z
    end type position

    public operator(+)
    interface operator(+)
      procedure pointadd
    end interface operator(+)

    public operator(-)
    interface operator(-)
      procedure pointsub
    end interface operator(-)

    public operator(*)
    interface operator(*)
      procedure pointIntscalar,&
                 pointRealscalar
    end interface operator(*)

    public assignment(=)
    interface assignment(=)
      procedure pointassign
    end interface assignment(=)

    type(position), dimension(:,:,:), allocatable ::  point
    
    ! Public methods
    public :: setup_grid
    public :: destroy_grid
    public :: point

    contains

        subroutine pointassign(p1,p2)
          implicit none
          type(position), intent(out) :: p1
          type(position), intent(in) :: p2
          p1%x = p2%x
          p1%y = p2%y
          p1%z = p2%z
        end subroutine pointassign

        function pointadd(p2, p3) result(p1)
          implicit none
          type(position) :: p1
          type(position), intent(in) :: p2,p3
          p1%x = p2%x + p3%x
          p1%y = p2%y + p3%y
          p1%z = p2%z + p3%z
        end function pointadd
    
        function pointsub(p2, p3) result(p1)
          implicit none
          type(position) :: p1
          type(position), intent(in) :: p2,p3
          p1%x = p2%x - p3%x
          p1%y = p2%y - p3%y
          p1%z = p2%z - p3%z
        end function pointsub
    
        function pointIntscalar(p2,p3) result(p1)
          implicit none
          integer, intent(in)    :: p2
          type(position), intent(in) :: p3
          type(position) :: p1
          p1%x = p2*p3%x
          p1%y = p2*p3%y
          p1%z = p2*p3%z
        end function pointIntscalar

        function pointRealscalar(p2,p3) result(p1)
          implicit none
          type(position)  :: p1
          real , intent(in) :: p2
          type(position), intent(in) :: p3
          p1%x = p2*p3%x
          p1%y = p2*p3%y
          p1%z = p2*p3%z
        end function pointRealscalar
    
    
        subroutine allocate_memory()
            !< Allocate memory to store the grid
            !-----------------------------------------------------------

            implicit none

            DebugCall("allocate memory to grid")
            allocate(point(-2:imx+3, -2:jmx+3, -2:kmx+3))
            !call dmsg(1, 'grid', 'allocate_memory')
            
        end subroutine allocate_memory

        subroutine destroy_grid()
            !< Deallocate the memory allocated for the grid.
            !-----------------------------------------------------------

            implicit none

            call dmsg(1, 'grid', 'destroy_memory')

            !call dealloc(grid_x)
            !call dealloc(grid_y)
            !call dealloc(grid_z)

!           call dealloc(sphere_indices)

        end subroutine destroy_grid

        subroutine setup_grid(gridfile)
            !< Read the grid file and initialize the grid
            !-----------------------------------------------------------

            implicit none
            character(len=*), intent(in) :: gridfile
            
            DebugCall('setup_grid')

            open(GRID_FILE_UNIT, file=gridfile)
            call extract_grid_size()
            call allocate_memory()
            !read interface mapping
            call read_interface_map()

            ! ghost grid exchange
            call populate_grid_points()

            close(GRID_FILE_UNIT)

            ! populate ghost grid points
            call ghost_grid()
        
        end subroutine setup_grid
        
        subroutine extract_grid_size()
            !< Extract the grid size from the grid file header
            !
            ! We assume that the grid could be in 1 or 2 dimensions. If
            ! the grid is in 1 dimension, jmx will be set to 1.
            ! We assume that at least one number is specified in the 
            ! header, i.e., the grid has atleast one dimension.
            !-----------------------------------------------------------

            implicit none
            character(len=STRING_BUFFER_LENGTH) :: header
            integer :: ios  ! io operation status

            call dmsg(1, 'grid', 'extract_grid_size')

            read(GRID_FILE_UNIT, '(A)', iostat=ios) header
            if (ios /= 0) then
                print *, 'Error while reading grid file header.'
                print *, 'Current buffer length is set to: ', &
                        STRING_BUFFER_LENGTH
                stop
            end if

            ! Try to read constants corresponding to two dimensions.
            read(header, *, iostat=ios) imx, jmx, kmx
            if (ios /= 0) then
                ! An io error means it was not possible to read kmx
                ! This means the file does not have kmx  and so, set
                ! the extent of this direction to 1. Read the remaining
                ! dimension from the header.
                read(header, *, iostat=ios) imx, jmx
                if (ios /= 0) then
                    ! This means that jmx does not exist. Repeat again
                    read(header, *, iostat=ios) imx
                    if (ios /= 0) then
                        ! Error while reading
                        print *, 'Unable to read grid extent.'
                        stop
                    end if
                    jmx = 1
                end if
                kmx = 1
            end if

        end subroutine extract_grid_size


        subroutine extract_grid_point(line, i, j, k)
            !< Extract a grid point from a line of the grid file. 
            !-----------------------------------------------------------

            implicit none
            character(len=STRING_BUFFER_LENGTH), intent(in) :: line
            integer, intent(in) :: i, j, k

            call dmsg(0, 'grid', 'extract_grid_point')

            if (kmx > 1) then
                read(line, *) point(i, j, k)%x, point(i, j, k)%y, point(i, j, k)%z
            else    
                if (jmx > 1) then
                    read(line, *) point(i, j, k)%x, point(i, j, k)%y
                else
                    read(line, *) point(i, j, k)%x
                    point(i, j, k)%y = 0.
                end if
                point(i, j, k)%z = 0.
            end if
        end subroutine extract_grid_point

        subroutine populate_grid_points()
            !< Use the grid file to populate the grid points.
            !-----------------------------------------------------------

            implicit none
            character(len=STRING_BUFFER_LENGTH) :: line
            integer :: i, j, k
            integer :: ios  ! io status

            call dmsg(1, 'grid', 'populate_grid_point')
         !  print *, imx, jmx, kmx

            ! Read grid points from the grid file
            do k = 1, kmx
                do j = 1, jmx
                    do i = 1, imx
                        read(GRID_FILE_UNIT, '(A)', iostat=ios) line
                        if (ios /= 0) then
                            print *, 'Error while reading grid line.'
                            print *, 'Current grid point: ', i, j, k
                            print *, 'Current buffer length is set to: ', &
                                     STRING_BUFFER_LENGTH
                            print *, 'Exiting program.'
                            stop
                        end if
                        call extract_grid_point(line, i, j, k)
                    end do
                end do
            end do

        end subroutine populate_grid_points

        subroutine ghost_grid()
          !< generate ghost grid for the various operations later.
          DebugCall('ghost_grid')
          implicit none
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
          point( 0,:,:)%x=2*point( 1,:,:)%x-point(2,:,:)%x
          point( 0,:,:)%y=2*point( 1,:,:)%y-point(2,:,:)%y
          point( 0,:,:)%z=2*point( 1,:,:)%z-point(2,:,:)%z
          point(-1,:,:)%x=2*point( 0,:,:)%x-point(1,:,:)%x
          point(-1,:,:)%y=2*point( 0,:,:)%y-point(1,:,:)%y
          point(-1,:,:)%z=2*point( 0,:,:)%z-point(1,:,:)%z
          point(-2,:,:)%x=2*point(-1,:,:)%x-point(0,:,:)%x
          point(-2,:,:)%y=2*point(-1,:,:)%y-point(0,:,:)%y
          point(-2,:,:)%z=2*point(-1,:,:)%z-point(0,:,:)%z

          !imax face -> imx+1 grid point
          !imax face -> imx+2 grid point
          !imax face -> imx+3 grid point
          point(imx+1,:,:)%x=2*point(imx+0,:,:)%x-point(imx-1,:,:)%x
          point(imx+1,:,:)%y=2*point(imx+0,:,:)%y-point(imx-1,:,:)%y
          point(imx+1,:,:)%z=2*point(imx+0,:,:)%z-point(imx-1,:,:)%z
          point(imx+2,:,:)%x=2*point(imx+1,:,:)%x-point(imx-0,:,:)%x
          point(imx+2,:,:)%y=2*point(imx+1,:,:)%y-point(imx-0,:,:)%y
          point(imx+2,:,:)%z=2*point(imx+1,:,:)%z-point(imx-0,:,:)%z
          point(imx+3,:,:)%x=2*point(imx+2,:,:)%x-point(imx+1,:,:)%x
          point(imx+3,:,:)%y=2*point(imx+2,:,:)%y-point(imx+1,:,:)%y
          point(imx+3,:,:)%z=2*point(imx+2,:,:)%z-point(imx+1,:,:)%z


          !--- Jmin faces ---!
          !jmin faces -> 0 grid point
          !jmin face -> -1 grid point
          !jmin face -> -2 grid point
          point(:, 0,:)%x=2*point(:, 1,:)%x-point(:,2,:)%x
          point(:, 0,:)%y=2*point(:, 1,:)%y-point(:,2,:)%y
          point(:, 0,:)%z=2*point(:, 1,:)%z-point(:,2,:)%z
          point(:,-1,:)%x=2*point(:, 0,:)%x-point(:,1,:)%x
          point(:,-1,:)%y=2*point(:, 0,:)%y-point(:,1,:)%y
          point(:,-1,:)%z=2*point(:, 0,:)%z-point(:,1,:)%z
          point(:,-2,:)%x=2*point(:,-1,:)%x-point(:,0,:)%x
          point(:,-2,:)%y=2*point(:,-1,:)%y-point(:,0,:)%y
          point(:,-2,:)%z=2*point(:,-1,:)%z-point(:,0,:)%z

          !jmax face -> jmx+1 grid point
          !jmax face -> jmx+3 grid point
          !jmax face -> jmx+2 grid point
          point(:,jmx+1,:)%x=2*point(:,jmx+0,:)%x-point(:,jmx-1,:)%x
          point(:,jmx+1,:)%y=2*point(:,jmx+0,:)%y-point(:,jmx-1,:)%y
          point(:,jmx+1,:)%z=2*point(:,jmx+0,:)%z-point(:,jmx-1,:)%z
          point(:,jmx+2,:)%x=2*point(:,jmx+1,:)%x-point(:,jmx-0,:)%x
          point(:,jmx+2,:)%y=2*point(:,jmx+1,:)%y-point(:,jmx-0,:)%y
          point(:,jmx+2,:)%z=2*point(:,jmx+1,:)%z-point(:,jmx-0,:)%z
          point(:,jmx+3,:)%x=2*point(:,jmx+2,:)%x-point(:,jmx+1,:)%x
          point(:,jmx+3,:)%y=2*point(:,jmx+2,:)%y-point(:,jmx+1,:)%y
          point(:,jmx+3,:)%z=2*point(:,jmx+2,:)%z-point(:,jmx+1,:)%z


          !--- Kmax faces ---!
          !kmin faces -> 0 grid point
          !kmin face -> -1 grid point
          !kmin face -> -2 grid point
          point(:,:, 0)%x=2*point(:,:, 1)%x-point(:,:,2)%x
          point(:,:, 0)%y=2*point(:,:, 1)%y-point(:,:,2)%y
          point(:,:, 0)%z=2*point(:,:, 1)%z-point(:,:,2)%z
          point(:,:,-1)%x=2*point(:,:, 0)%x-point(:,:,1)%x
          point(:,:,-1)%y=2*point(:,:, 0)%y-point(:,:,1)%y
          point(:,:,-1)%z=2*point(:,:, 0)%z-point(:,:,1)%z
          point(:,:,-2)%x=2*point(:,:,-1)%x-point(:,:,0)%x
          point(:,:,-2)%y=2*point(:,:,-1)%y-point(:,:,0)%y
          point(:,:,-2)%z=2*point(:,:,-1)%z-point(:,:,0)%z

          !kmax face -> kmx+1 grid point
          !kmax face -> kmx+2 grid point
          !kmax face -> kmx+3 grid point
          point(:,:,kmx+1)%x=2*point(:,:,kmx+0)%x-point(:,:,kmx-1)%x
          point(:,:,kmx+1)%y=2*point(:,:,kmx+0)%y-point(:,:,kmx-1)%y
          point(:,:,kmx+1)%z=2*point(:,:,kmx+0)%z-point(:,:,kmx-1)%z
          point(:,:,kmx+2)%x=2*point(:,:,kmx+1)%x-point(:,:,kmx-0)%x
          point(:,:,kmx+2)%y=2*point(:,:,kmx+1)%y-point(:,:,kmx-0)%y
          point(:,:,kmx+2)%z=2*point(:,:,kmx+1)%z-point(:,:,kmx-0)%z
          point(:,:,kmx+3)%x=2*point(:,:,kmx+2)%x-point(:,:,kmx+1)%x
          point(:,:,kmx+3)%y=2*point(:,:,kmx+2)%y-point(:,:,kmx+1)%y
          point(:,:,kmx+3)%z=2*point(:,:,kmx+2)%z-point(:,:,kmx+1)%z

        end subroutine ghost_grid

end module grid
