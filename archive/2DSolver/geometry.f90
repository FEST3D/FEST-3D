module geometry
    !-------------------------------------------------------------------
    ! The geometry module contains various geometrical quantities like 
    ! normals, areas and volumes to be used in computations. 
    !
    ! This version of the geometry module assumes the grid is atmost
    ! 2-dimensional. 
    !-------------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx, grid_x, grid_y

    implicit none
    private

    ! Grid face normals
    real, public, dimension(:, :), allocatable, target :: xnx, xny
    real, public, dimension(:, :), allocatable, target :: ynx, yny
    ! Grid face areas
    real, public, dimension(:, :), allocatable, target :: xA, yA
    ! Grid cell volumes
    real, public, dimension(:, :), allocatable :: volume
    ! Ghost cell centroid
    real, public, dimension(:, :), allocatable, target :: left_ghost_centroid, &
             right_ghost_centroid, top_ghost_centroid, bottom_ghost_centroid

    ! Public methods
    public :: setup_geometry
    public :: destroy_geometry

    contains

        subroutine allocate_memory_volumes()
            !-----------------------------------------------------------
            ! Allocate memory for the volume variables.
            !-----------------------------------------------------------
            
            implicit none

            call alloc(volume, 1, imx-1, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for volume.')

        end subroutine allocate_memory_volumes

        subroutine allocate_memory_areas()
            !-----------------------------------------------------------
            ! Allocate memory for the area variables.
            !-----------------------------------------------------------
            
            implicit none

            call alloc(xA, 1, imx, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for xA.')
            call alloc(yA, 1, imx-1, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for yA.')

        end subroutine allocate_memory_areas

        subroutine allocate_memory_normals()
            !-----------------------------------------------------------
            ! Allocate memory for the normal variables.
            !-----------------------------------------------------------
                        
            implicit none

            call alloc(xnx, 1, imx, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for xnx.')
            call alloc(xny, 1, imx, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for xny.')
            call alloc(ynx, 1, imx-1, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for ynx.')
            call alloc(yny, 1, imx-1, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for yny.')

        end subroutine allocate_memory_normals

        subroutine allocate_memory_ghost_centroids()

            implicit none

            call alloc(left_ghost_centroid, 1, jmx-1, 1, 2, &
                    errmsg='Error: Unable to allocate memory for left_ghost_centroid')
            call alloc(right_ghost_centroid, 1, jmx-1, 1, 2, &
                    errmsg='Error: Unable to allocate memory for right_ghost_centroid')
            call alloc(top_ghost_centroid, 1, imx-1, 1, 2, &
                    errmsg='Error: Unable to allocate memory for top_ghost_centroid')
            call alloc(bottom_ghost_centroid, 1, imx-1, 1, 2, &
                    errmsg='Error: Unable to allocate memory for bottom_ghost_centroid')

        end subroutine allocate_memory_ghost_centroids

        subroutine allocate_memory()
            !-----------------------------------------------------------
            ! Allocate memory for the required variables.
            !-----------------------------------------------------------
            
            implicit none

            call dmsg(1, 'geometry', 'allocate_memory')

            call allocate_memory_normals()
            call allocate_memory_areas()
            call allocate_memory_volumes()
            call allocate_memory_ghost_centroids()

        end subroutine allocate_memory

        subroutine deallocate_memory()

            implicit none

            call dmsg(1, 'geometry', 'deallocate_memory')

            call dealloc(xnx)
            call dealloc(xny)
            call dealloc(ynx)
            call dealloc(yny)
            call dealloc(xA)
            call dealloc(yA)
            call dealloc(volume)
            call dealloc(left_ghost_centroid)
            call dealloc(right_ghost_centroid)
            call dealloc(top_ghost_centroid)
            call dealloc(bottom_ghost_centroid)

        end subroutine deallocate_memory

        subroutine normalize_face_normals()
            !-----------------------------------------------------------
            ! Normalize the face area vectors computed to get normals.
            !
            ! This should be used after the face area vectors have been 
            ! computed and placed in the face normal variables. 
            ! The face areas should also have been computed and placed 
            ! in the face area variables. 
            !-----------------------------------------------------------
            
            implicit none

            xnx(:, :) = xnx(:, :) / xA(:, :)
            xny(:, :) = xny(:, :) / xA(:, :)

            ynx(:, :) = ynx(:, :) / yA(:, :)
            yny(:, :) = yny(:, :) / yA(:, :)
        
        end subroutine normalize_face_normals

        subroutine compute_face_areas()
            !-----------------------------------------------------------
            ! Compute face areas based on area vectors
            !
            ! The face areas are computed using the face area vectors. 
            ! Prior to using this subroutine, the face area vectors must
            ! computed and placed in the face normal variables. 
            !-----------------------------------------------------------
            
            implicit none

            xA(:, :) = sqrt((xnx(:, :)) ** 2. + (xny(:, :)) ** 2.)

            yA(:, :) = sqrt((ynx(:, :)) ** 2. + (yny(:, :)) ** 2.)

        end subroutine compute_face_areas

        subroutine compute_face_area_vectors()
            !-----------------------------------------------------------
            ! Compute face area vectors
            !
            ! The face area vectors denote the face area both in 
            ! magnitude and direction. They are placed in the face 
            ! normal variables for further calculation.
            !
            ! This is the 2-dimensional version: in this case, the face
            ! areas default to edge lengths.
            !-----------------------------------------------------------
            
            implicit none

            xnx(1:imx, 1:jmx-1) = grid_y(1:imx, 2:jmx) - grid_y(1:imx, 1:jmx-1)
            xny(1:imx, 1:jmx-1) = - (grid_x(1:imx, 2:jmx) - &
                    grid_x(1:imx, 1:jmx-1))
            
            ynx(1:imx-1, 1:jmx) = - (grid_y(2:imx, 1:jmx) - &
                    grid_y(1:imx-1, 1:jmx))
            yny(1:imx-1, 1:jmx) = grid_x(2:imx, 1:jmx) - grid_x(1:imx-1, 1:jmx)

        end subroutine compute_face_area_vectors

        subroutine compute_face_areas_and_normals()
            !-----------------------------------------------------------
            ! Compute the face areas and normals
            !
            ! This is the 2-dimensional version. In this case, the face 
            ! areas default to edge lengths.
            !-----------------------------------------------------------

            implicit none

            call compute_face_area_vectors()
            call compute_face_areas()
            call normalize_face_normals()
        
        end subroutine compute_face_areas_and_normals

        subroutine compute_volumes()
            !-----------------------------------------------------------
            ! Compute the grid cell volumes using the grid points.
            !-----------------------------------------------------------

            implicit none

            volume(1:imx-1, 1:jmx-1) = 0.5 * abs( &
                    (grid_x(1:imx-1, 1:jmx-1) * grid_y(1:imx-1, 2:jmx)) - &
                    (grid_y(1:imx-1, 1:jmx-1) * grid_x(1:imx-1, 2:jmx)) + &
                    (grid_x(1:imx-1, 2:jmx) * grid_y(2:imx, 2:jmx)) - &
                    (grid_y(1:imx-1, 2:jmx) * grid_x(2:imx, 2:jmx)) + &
                    (grid_x(2:imx, 2:jmx) * grid_y(2:imx, 1:jmx-1)) - &
                    (grid_y(2:imx, 2:jmx) * grid_x(2:imx, 1:jmx-1)) + &
                    (grid_x(2:imx, 1:jmx-1) * grid_y(1:imx-1, 1:jmx-1)) - &
                    (grid_y(2:imx, 1:jmx-1) * grid_x(1:imx-1, 1:jmx-1)) &
                    )
            
        end subroutine compute_volumes

        subroutine compute_geometric_parameters()
            !-----------------------------------------------------------
            ! Compute the geometric parameters based on the grid points
            !
            ! The geometric parameters include the face normals and 
            ! areas and the cell volumes.
            !-----------------------------------------------------------
            
            implicit none

            call dmsg(1, 'geometry', 'compute_geometric_parameters')

            call compute_face_areas_and_normals()
            call compute_volumes()

        end subroutine compute_geometric_parameters

        subroutine compute_ghost_cell_centroid()
            !-----------------------------------------------------------
            ! Computes the centroid of the ghost cell. To be used in
            ! viscous module
            !
            ! The ghost cell centroid is found as follows:
            ! The face plane is idealised as a plane passing through the
            ! centroid of the face (average of 4 points of a face) and 
            ! the face normal. 
            !
            ! The required centroid is then found as the mirror of the
            ! centroid of the entire element with respect to the above
            ! defined plane
            !-----------------------------------------------------------

            ! a is the vector from the face centroid to the centroid of
            ! the element
            ! The vector r_ghost_centroid - r_face_centroid = r (say)
            ! is given by the equation:
            ! r = a - 2(a.n)n
            ! Hence, r_ghost_centroid = r_face_centoid + a - 2(a.n)n
            ! n is the face normals, which are calculated in previous
            ! subroutines
            ! Note that the formula is invariant of the direction of n
            real, dimension(2) :: a
            real, dimension(2) :: face_centroid, centroid
            integer :: i, j
            
            ! left face ghost cell centroids. i = 1
            do j = 1, jmx - 1
                centroid(1) = (grid_x(1, j) + grid_x(2, j) + &
                               grid_x(2, j+1) + grid_x(1, j+1)) * 0.25
                centroid(2) = (grid_y(1, j) + grid_y(2, j) + &
                               grid_y(2, j+1) + grid_y(1, j+1)) * 0.25
                
                face_centroid(1) = (grid_x(1, j) + grid_x(1, j+1)) * 0.5 
                face_centroid(2) = (grid_y(1, j) + grid_y(1, j+1)) * 0.5 

                a(1) = centroid(1) - face_centroid(1)
                a(2) = centroid(2) - face_centroid(2)

                left_ghost_centroid(j, 1) = face_centroid(1) + a(1) - &
                                    (2*(a(1)*xnx(1, j) + a(2)*xny(1, j)) * &
                                    xnx(1, j) )
                left_ghost_centroid(j, 2) = face_centroid(2) + a(2) - &
                                    (2*(a(1)*xnx(1, j) + a(2)*xny(1, j)) * &
                                    xny(1, j) )
            end do

            ! right face ghost cell centroids. i = imx
            do j = 1, jmx - 1
                centroid(1) = (grid_x(imx-1, j) + grid_x(imx, j) + &
                               grid_x(imx, j+1) + grid_x(imx-1, j+1)) * 0.25
                centroid(2) = (grid_y(imx-1, j) + grid_y(imx, j) + &
                               grid_y(imx, j+1) + grid_y(imx-1, j+1)) * 0.25
                
                face_centroid(1) = (grid_x(imx, j) + grid_x(imx, j+1)) * 0.5 
                face_centroid(2) = (grid_y(imx, j) + grid_y(imx, j+1)) * 0.5 

                a(1) = centroid(1) - face_centroid(1)
                a(2) = centroid(2) - face_centroid(2)

                right_ghost_centroid(j, 1) = face_centroid(1) + a(1) - &
                                    (2*(a(1)*xnx(imx, j) + a(2)*xny(imx, j)) * &
                                    xnx(imx, j) )
                right_ghost_centroid(j, 2) = face_centroid(2) + a(2) - &
                                    (2*(a(1)*xnx(imx, j) + a(2)*xny(imx, j)) * &
                                    xny(imx, j) )
            end do

            ! bottom face ghost cell centroids. j = 1
            do i = 1, imx - 1
                centroid(1) = (grid_x(i, 1) + grid_x(i+1, 1) + &
                               grid_x(i+1, 2) + grid_x(i, 2)) * 0.25
                centroid(2) = (grid_y(i, 1) + grid_y(i+1, 1) + &
                               grid_y(i+1, 2) + grid_y(i, 2)) * 0.25
                
                face_centroid(1) = (grid_x(i, 1) + grid_x(i+1, 1)) * 0.5
                face_centroid(2) = (grid_y(i, 1) + grid_y(i+1, 1)) * 0.5
                a(1) = centroid(1) - face_centroid(1)
                a(2) = centroid(2) - face_centroid(2)

                bottom_ghost_centroid(i, 1) = face_centroid(1) + a(1) - &
                                    (2*(a(1)*ynx(i, 1) + a(2)*yny(i, 1))*ynx(i, 1) )
                bottom_ghost_centroid(i, 2) = face_centroid(2) + a(2) - &
                                    (2*(a(1)*ynx(i, 1) + a(2)*yny(i, 1))*yny(i, 1) )
            end do

            ! top face ghost cell centroids. j = jmx
            do i = 1, imx - 1
                centroid(1) = (grid_x(i, jmx-1) + grid_x(i+1, jmx-1) + &
                               grid_x(i+1, jmx) + grid_x(i, jmx)) * 0.25
                centroid(2) = (grid_y(i, jmx-1) + grid_y(i+1, jmx-1) + &
                               grid_y(i+1, jmx) + grid_y(i, jmx)) * 0.25
                
                face_centroid(1) = (grid_x(i, jmx) + grid_x(i+1, jmx)) * 0.5
                face_centroid(2) = (grid_y(i, jmx) + grid_y(i+1, jmx)) * 0.5
                a(1) = centroid(1) - face_centroid(1)
                a(2) = centroid(2) - face_centroid(2)

                top_ghost_centroid(i, 1) = face_centroid(1) + a(1) - &
                                    (2*(a(1)*ynx(i, jmx) + a(2)*yny(i, jmx))*ynx(i, jmx) )
                top_ghost_centroid(i, 2) = face_centroid(2) + a(2) - &
                                    (2*(a(1)*ynx(i, jmx) + a(2)*yny(i, jmx))*yny(i, jmx) )
            end do

          ! open(24, file='ghost_centroids.txt')
          ! do j = 1, jmx - 1
          !     write (24, *) left_ghost_centroid(j, 1), left_ghost_centroid(j, 2)
          !     write (24, *) right_ghost_centroid(j, 1), right_ghost_centroid(j, 2)
          ! end do
          ! do i = 1, imx - 1
          !     write (24, *) bottom_ghost_centroid(i, 1), bottom_ghost_centroid(i, 2)
          !     write (24, *) top_ghost_centroid(i, 1), top_ghost_centroid(i, 2)
          ! end do
          ! close(24)

        end subroutine compute_ghost_cell_centroid

        subroutine setup_geometry()
            !-----------------------------------------------------------
            ! Make the geometry module useful
            !
            ! Allocates memory to the variables and initializes them.
            !-----------------------------------------------------------

            implicit none

            call dmsg(1, 'geometry', 'setup_geometry')

            call allocate_memory()
            call compute_geometric_parameters()
            call compute_ghost_cell_centroid()

        end subroutine setup_geometry

        subroutine destroy_geometry()

            implicit none
            
            call dmsg(1, 'geometry', 'destroy_geometry')

            call deallocate_memory()

        end subroutine destroy_geometry

end module geometry
