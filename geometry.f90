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
    real, public, dimension(:, :), allocatable :: xnx, xny
    real, public, dimension(:, :), allocatable :: ynx, yny
    ! Grid face areas
    real, public, dimension(:, :), allocatable :: xA, yA
    ! Grid cell volumes
    real, public, dimension(:, :), allocatable :: volume

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

        subroutine allocate_memory()
            !-----------------------------------------------------------
            ! Allocate memory for the required variables.
            !-----------------------------------------------------------
            
            implicit none

            call dmsg(1, 'geometry', 'allocate_memory')

            call allocate_memory_normals()
            call allocate_memory_areas()
            call allocate_memory_volumes()

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

        end subroutine setup_geometry

        subroutine destroy_geometry()

            implicit none
            
            call dmsg(1, 'geometry', 'destroy_geometry')

            call deallocate_memory()

        end subroutine destroy_geometry

end module geometry
