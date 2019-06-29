  !< Calculate the distance from the wall 
  !< for each cell-center in the domain
module wall_dist
  !< Calculate the distance from the wall 
  !< for each cell-center in the domain
  use global,  only: NODESURF_FILE_UNIT
  use global,  only: WALL_DIST_FILE_UNIT
  use global,  only: wall_dist_file
  use global,  only: surface_node_points

  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx
  use global_vars, only : grid_x
  use global_vars, only : grid_y
  use global_vars, only : grid_z
  use global_vars, only : dist

  use utils, only: alloc, dealloc, dmsg
!  use grid, only: imx,jmx, kmx, grid_x, grid_y, grid_z

  implicit none
  private

  integer                                     :: n_surfnodes
  !< Number of surfce node points
!  real, public, dimension(:,:,:), allocatable :: dist
  real, private,dimension(:)    , allocatable :: wall_x
  !< X component of wall surface node point
  real, private,dimension(:)    , allocatable :: wall_y
  !< Y component of wall surface node point
  real, private,dimension(:)    , allocatable :: wall_z
  !< Z component of wall surface node point

  public :: setup_wall_dist
  public :: destroy_wall_dist
  public :: find_wall_dist

  contains

    subroutine setup_wall_dist()
      !< Allocate memory to the wall_distance variables
      !< and read the surface node file

      implicit none

      call dmsg(1, 'wall_dist', 'setup_wall_dist')
      call setup_nodefile()
      call alloc(wall_x, 1, n_surfnodes,&
                "ERROR: unale to allocate memory to 'Dist' variable " )
      call alloc(wall_y, 1, n_surfnodes,&
                "ERROR: unale to allocate memory to 'Dist' variable " )
      call alloc(wall_z, 1, n_surfnodes,&
                "ERROR: unale to allocate memory to 'Dist' variable " )
      call alloc(dist, -2, imx+2, -2, jmx+2, -2, kmx+2, &
                "ERROR: unale to allocate memory to 'Dist' variable " )
      call read_destroy_nodefile

    end subroutine setup_wall_dist



    subroutine destroy_wall_dist()
      !< Deallocate the memory of wall_distance variable,
      !< wall_x, wall_y, and wall_z

      implicit none

      call dmsg(1, 'wall_dist', 'destroy_wall_dist')
      call dealloc(wall_x)
      call dealloc(wall_y)
      call dealloc(wall_z)
      call dealloc(dist)

    end subroutine destroy_wall_dist


    subroutine setup_nodefile()
      !< Open and read first line of surface_node_point file
      implicit none
      integer :: ios
      open(NODESURF_FILE_UNIT, file=surface_node_points, status='old', IOSTAT=ios)
      if(ios/=0) then
        call dmsg(5, 'wall_dist', 'setup_nodefile', &
          "!!! -->file containg surface nodepoints not found" )
      end if
      read(NODESURF_FILE_UNIT, *) n_surfnodes

    end subroutine setup_nodefile


    subroutine read_destroy_nodefile()
      !< Read, and close surface_node_point file
      implicit none
      integer :: i
      do i = 1, n_surfnodes
        read(NODESURF_FILE_UNIT, *) wall_x(i), wall_y(i), wall_z(i)
      end do
      close(NODESURF_FILE_UNIT)
    end subroutine read_destroy_nodefile

    subroutine find_wall_dist()
      !< Determine the minimum wall distance from the wall surface node points

      implicit none

      integer :: i,j,k,n
      real :: current_dist
      real, dimension(:,:,:), allocatable :: node_dist
      call dmsg(1, 'wall_dist', 'find_wall_dist')
      call alloc(node_dist,-2,imx+3,-2,jmx+3,-2,kmx+3)

      do k = -2,kmx+3
        do j = -2,jmx+3
          do i = -2,imx+3
            node_dist(i,j,k) = 1.e+20
            do n = 1,n_surfnodes

            current_dist = sqrt((wall_x(n)-grid_x(i,j,k))**2&
                               +(wall_y(n)-grid_y(i,j,k))**2&
                               +(wall_z(n)-grid_z(i,j,k))**2&
                               ) 
            node_dist(i,j,k) = min(node_dist(i,j,k),current_dist)

            end do
          end do
        end do
      end do
      do k=-2,kmx+2
        do j=-2,jmx+2
          do i=-2,imx+2
            dist(i,j,k) = 0.125*(node_dist(i  ,j  ,k  )&
                                +node_dist(i  ,j+1,k  )&
                                +node_dist(i  ,j+1,k+1)&
                                +node_dist(i  ,j  ,k+1)&
                                +node_dist(i+1,j  ,k+1)&
                                +node_dist(i+1,j  ,k  )&
                                +node_dist(i+1,j+1,k  )&
                                +node_dist(i+1,j+1,k+1)&
                                )
          end do
        end do
      end do
      call dealloc(node_dist)
      call dmsg(1, 'wall_dist', 'find_wall_dist-> complete')

    end subroutine find_wall_dist


end module wall_dist
