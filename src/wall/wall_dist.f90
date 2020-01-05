  !< Calculate the distance from the wall 
  !< for each cell-center in the domain
module wall_dist
  !< Calculate the distance from the wall 
  !< for each cell-center in the domain
  use vartypes
!  use global_vars, only : dist
  use utils, only: alloc
#include "../debug.h"
#include "../error.h"

  implicit none
  private

  integer                                     :: n_surfnodes
  !< Number of surfce node points
!  real(wp), public, dimension(:,:,:), allocatable :: dist
  real(wp), private,dimension(:)    , allocatable :: wall_x
  !< X component of wall surface node point
  real(wp), private,dimension(:)    , allocatable :: wall_y
  !< Y component of wall surface node point
  real(wp), private,dimension(:)    , allocatable :: wall_z
  !< Z component of wall surface node point
  real(wp), dimension(:, :, :), allocatable             :: dist 

  integer :: imx, jmx, kmx

  public :: setup_wall_dist
  public :: find_wall_dist
  public :: dist

  contains

    subroutine setup_wall_dist(files, dims)
      !< Allocate memory to the wall_distance variables
      !< and read the surface node file

      implicit none
      type(filetype), intent(in) :: files
      type(extent), intent(in) :: dims
      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx

      DebugCall('setup_wall_dist')
      call setup_nodefile(files)
      call alloc(wall_x, 1, n_surfnodes,&
                "ERROR: unale to allocate memory to 'Dist' variable " )
      call alloc(wall_y, 1, n_surfnodes,&
                "ERROR: unale to allocate memory to 'Dist' variable " )
      call alloc(wall_z, 1, n_surfnodes,&
                "ERROR: unale to allocate memory to 'Dist' variable " )
      call alloc(dist, -2, imx+2, -2, jmx+2, -2, kmx+2, &
                "ERROR: unale to allocate memory to 'Dist' variable " )
      call read_destroy_nodefile(files)

    end subroutine setup_wall_dist


!
!    subroutine destroy_wall_dist()
!      !< Deallocate the memory of wall_distance variable,
!      !< wall_x, wall_y, and wall_z
!
!      implicit none
!
!      DebugCall('destroy_wall_dist')
!      call dealloc(wall_x)
!      call dealloc(wall_y)
!      call dealloc(wall_z)
!      call dealloc(dist)
!
!    end subroutine destroy_wall_dist
!

    subroutine setup_nodefile(files)
      !< Open and read first line of surface_node_point file
      implicit none
      type(filetype), intent(in) :: files
      integer :: ios
      open(files%NODESURF_FILE_UNIT, file=files%surface_node_points, status='old', IOSTAT=ios)
      if(ios/=0) then
        print*, "!!! -->file containg surface nodepoints not found"
        Fatal_error
      end if
      read(files%NODESURF_FILE_UNIT, *) n_surfnodes

    end subroutine setup_nodefile


    subroutine read_destroy_nodefile(files)
      !< Read, and close surface_node_point file
      implicit none
      type(filetype), intent(in) :: files
      integer :: i
      do i = 1, n_surfnodes
        read(files%NODESURF_FILE_UNIT, *) wall_x(i), wall_y(i), wall_z(i)
      end do
      close(files%NODESURF_FILE_UNIT)
    end subroutine read_destroy_nodefile

    subroutine find_wall_dist(nodes, dims)
      !< Determine the minimum wall distance from the wall surface node points

      implicit none
      type(extent), intent(in) :: dims
      type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(in) :: nodes

      integer :: i,j,k,n
      real(wp) :: current_dist
      real(wp), dimension(:,:,:), allocatable :: node_dist
      DebugCall('find_wall_dist')
      call alloc(node_dist,-2,imx+3,-2,jmx+3,-2,kmx+3)

      do k = -2,dims%kmx+3
        do j = -2,dims%jmx+3
          do i = -2,dims%imx+3
            node_dist(i,j,k) = 1.e+20
            do n = 1,n_surfnodes

            current_dist = sqrt((wall_x(n)-nodes(i,j,k)%x)**2&
                               +(wall_y(n)-nodes(i,j,k)%y)**2&
                               +(wall_z(n)-nodes(i,j,k)%z)**2&
                               ) 
            node_dist(i,j,k) = min(node_dist(i,j,k),current_dist)

            end do
          end do
        end do
      end do
      do k=-2,dims%kmx+2
        do j=-2,dims%jmx+2
          do i=-2,dims%imx+2
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
      deallocate(node_dist)
      DebugCall('find_wall_dist-> complete')

    end subroutine find_wall_dist


end module wall_dist
