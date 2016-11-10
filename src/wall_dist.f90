module wall_dist
  use utils, only: alloc, dealloc, dmsg
  use grid, only: imx,jmx, kmx, grid_x, grid_y, grid_z
  use surface, only: setup_surface, destroy_surface, surface_points, wallc,&
                    wall_x, wall_y, wall_z

  implicit none
  private

  real, public, dimension(:,:,:), allocatable :: dist

  public :: setup_wall_dist
  public :: destroy_wall_dist
  public :: find_wall_dist

  contains

    subroutine setup_wall_dist()

      implicit none

      call dmsg(1, 'wall_dist', 'setup_wall_dist')
      call setup_surface
      call alloc(dist, 1, imx-1, 1, jmx-1, 1, kmx-1, &
                "ERROR: unale to allocate memory to 'Dist' variable " )

    end subroutine setup_wall_dist



    subroutine destroy_wall_dist()

      implicit none

      call dmsg(1, 'wall_dist', 'destroy_wall_dist')
      call destroy_surface
      call dealloc(dist)

    end subroutine destroy_wall_dist



    subroutine find_wall_dist()

      implicit none

      integer :: i,j,k
      integer :: id,jd,kd
      real :: cc_x, cc_y, cc_z  !cell center
      call dmsg(1, 'wall_dist', 'find_wall_dist')

      call surface_points()
      open(607, file='distance.vtk')

      write(607,'(a)') '# vtk DataFile Version 3.1'
      write(607,'(a)')  'wall dist'
      write(607,'(a)')  'ASCII'
      write(607,'(a)')  'DATASET STRUCTURED_GRID'
      write(607,*)
      write(607,'(a, i0, a, i0, a, i0)')  'DIMENSIONS ', imx, ' ', jmx, ' ', kmx
      write(607, '(a, i0, a)')  'POINTS ', imx*jmx*kmx, ' DOUBLE'


      do k = 1,kmx
        do j = 1,jmx
          do i = 1,imx
            write(607, '(f0.16, a, f0.16, a, f0.16)') grid_x(i,j,k), ' ', &
                      grid_y(i,j,k), ' ', grid_z(i,j,k)
            end do
          end do
        end do

        write(607,*) 


        write(607, '(a, i0)') 'CELL_DATA ', (imx-1)*(jmx-1)*(kmx-1)

        write(607, '(a)') 'SCALARS WALL_DIST FLOAT'
        write(607, '(a)') 'LOOKUP_TABLE default'

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
         
                    write(607, '(f0.16)') dist(i,j,k)
          end do
        end do
      end do

      end subroutine find_wall_dist

end module wall_dist
