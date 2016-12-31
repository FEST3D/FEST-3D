module wall_dist
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
!  real, public, dimension(:,:,:), allocatable :: dist
  real, private,dimension(:)    , allocatable :: wall_x, wall_y, wall_z

  public :: setup_wall_dist
  public :: destroy_wall_dist
  public :: find_wall_dist
  public :: write_wall_dist

  contains

    subroutine setup_wall_dist()

      implicit none

      call dmsg(1, 'wall_dist', 'setup_wall_dist')
      call setup_nodefile()
      call alloc(wall_x, 1, n_surfnodes,&
                "ERROR: unale to allocate memory to 'Dist' variable " )
      call alloc(wall_y, 1, n_surfnodes,&
                "ERROR: unale to allocate memory to 'Dist' variable " )
      call alloc(wall_z, 1, n_surfnodes,&
                "ERROR: unale to allocate memory to 'Dist' variable " )
      call alloc(dist, 1, imx-1, 1, jmx-1, 1, kmx-1, &
                "ERROR: unale to allocate memory to 'Dist' variable " )
      call read_destroy_nodefile

    end subroutine setup_wall_dist



    subroutine destroy_wall_dist()

      implicit none

      call dmsg(1, 'wall_dist', 'destroy_wall_dist')
      call dealloc(wall_x)
      call dealloc(wall_y)
      call dealloc(wall_z)
      call dealloc(dist)

    end subroutine destroy_wall_dist


    subroutine setup_nodefile()
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
      implicit none
      integer :: i
      do i = 1, n_surfnodes
        read(NODESURF_FILE_UNIT, *) wall_x(i), wall_y(i), wall_z(i)
      end do
      close(NODESURF_FILE_UNIT)
    end subroutine read_destroy_nodefile

    subroutine find_wall_dist()

      implicit none

      integer :: i,j,k
      integer :: id,jd,kd
      real :: cc_x, cc_y, cc_z  !cell center
      call dmsg(1, 'wall_dist', 'find_wall_dist')

      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1

            id = 1
            jd = 1
            kd = 1

            cc_x = 0.125 *(                                   &  
                                  grid_x(i, j, k )            & 
                                  + grid_x(i+id, j, k)        & 
                                  + grid_x(i, j+jd, k)        &  
                                  + grid_x(i, j, k+kd)        &    
                                  + grid_x(i, j+jd, k+kd)     &    
                                  + grid_x(i+id, j, k+kd)     &  
                                  + grid_x(i+id, j+jd, k)     &  
                                  + grid_x(i+id, j+jd, k+kd)  &  
                          )
            cc_y = 0.125 *(                                   &  
                                    grid_y(i, j, k )          & 
                                  + grid_y(i+id, j, k)        & 
                                  + grid_y(i, j+jd, k)        &  
                                  + grid_y(i, j, k+kd)        &    
                                  + grid_y(i, j+jd, k+kd)     &    
                                  + grid_y(i+id, j, k+kd)     &  
                                  + grid_y(i+id, j+jd, k)     &  
                                  + grid_y(i+id, j+jd, k+kd)  &  
                          )
            cc_z = 0.125 *(                                   &  
                                    grid_z(i, j, k )          & 
                                  + grid_z(i+id, j, k)        & 
                                  + grid_z(i, j+jd, k)        &  
                                  + grid_z(i, j, k+kd)        &    
                                  + grid_z(i, j+jd, k+kd)     &    
                                  + grid_z(i+id, j, k+kd)     &  
                                  + grid_z(i+id, j+jd, k)     &  
                                  + grid_z(i+id, j+jd, k+kd)  &  
                          )
            dist(i, j, k) =  minval(                                   & 
                                   sqrt(                            &
                                         (wall_x(:) - cc_x)**2      &
                                       + (wall_y(:) - cc_y)**2      &
                                       + (wall_z(:) - cc_z)**2      &
                              )                                     &
                      )
         
          end do
        end do
      end do

    end subroutine find_wall_dist

    subroutine write_wall_dist()
        implicit none
      integer :: i,j,k
      call dmsg(1, 'wall_dist', 'write_wall_dist')

      open(WALL_DIST_FILE_UNIT, file=wall_dist_file)

      write(WALL_DIST_FILE_UNIT,'(a)') '# vtk DataFile Version 3.1'
      write(WALL_DIST_FILE_UNIT,'(a)')  'wall dist'
      write(WALL_DIST_FILE_UNIT,'(a)')  'ASCII'
      write(WALL_DIST_FILE_UNIT,'(a)')  'DATASET STRUCTURED_GRID'
      write(WALL_DIST_FILE_UNIT,*)
      write(WALL_DIST_FILE_UNIT,'(a, i0, a, i0, a, i0)')  'DIMENSIONS ', imx, ' ', jmx, ' ', kmx
      write(WALL_DIST_FILE_UNIT, '(a, i0, a)')  'POINTS ', imx*jmx*kmx, ' DOUBLE'


      do k = 1,kmx
        do j = 1,jmx
          do i = 1,imx
            write(WALL_DIST_FILE_UNIT, '(f0.16, a, f0.16, a, f0.16)') grid_x(i,j,k), ' ', &
                      grid_y(i,j,k), ' ', grid_z(i,j,k)
            end do
          end do
        end do

        write(WALL_DIST_FILE_UNIT,*) 


        write(WALL_DIST_FILE_UNIT, '(a, i0)') 'CELL_DATA ', (imx-1)*(jmx-1)*(kmx-1)

        write(WALL_DIST_FILE_UNIT, '(a)') 'SCALARS WALL_DIST FLOAT'
        write(WALL_DIST_FILE_UNIT, '(a)') 'LOOKUP_TABLE default'

      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1
            write(WALL_DIST_FILE_UNIT, '(f0.16)') dist(i,j,k)
          end do
        end do
      end do
    end subroutine write_wall_dist


end module wall_dist
