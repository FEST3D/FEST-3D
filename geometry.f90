module geometry
    !-------------------------------------------------------------------
    ! The geometry module contains various geometrical quantities like 
    ! normals, areas and volumes to be used in computations. 
    !
    ! This version of the geometry module assumes the grid is atmost
    ! 2-dimensional. 
    !-------------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx, kmx, grid_x, grid_y, grid_z

    implicit none
    private

    ! Grid face normals
    real, public, dimension(:, :, :), allocatable :: xnx, xny, xnz
    real, public, dimension(:, :, :), allocatable :: ynx, yny, ynz
    real, public, dimension(:, :, :), allocatable :: znx, zny, znz
    ! Grid face areas
    real, public, dimension(:, :, :), allocatable :: xA, yA, zA
    ! Grid cell volumes
    real, public, dimension(:, :, :), allocatable :: volume
    ! Ghost cell centroid
    real, public, dimension(:, :, :), allocatable :: left_ghost_centroid, &
        right_ghost_centroid, front_ghost_centroid, back_ghost_centroid, &
        top_ghost_centroid, bottom_ghost_centroid

    ! Public methods
    public :: setup_geometry
    public :: destroy_geometry

    contains

        subroutine allocate_memory_volumes()
            !-----------------------------------------------------------
            ! Allocate memory for the volume variables.
            !-----------------------------------------------------------
            
            implicit none

            call alloc(volume, 1, imx-1, 1, jmx-1, 1, kmx -1, &
                    errmsg='Error: Unable to allocate memory for volume.')

        end subroutine allocate_memory_volumes

        subroutine allocate_memory_areas()
            !-----------------------------------------------------------
            ! Allocate memory for the area variables.
            !-----------------------------------------------------------
            
            implicit none

            call alloc(xA, 1, imx, 1, jmx-1, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for xA.')
            call alloc(yA, 1, imx-1, 1, jmx, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for yA.')
            call alloc(zA, 1, imx-1, 1, jmx-1, 1, kmx, &
                    errmsg='Error: Unable to allocate memory for yA.')

        end subroutine allocate_memory_areas

        subroutine allocate_memory_normals()
            !-----------------------------------------------------------
            ! Allocate memory for the normal variables.
            !-----------------------------------------------------------
                        
            implicit none

            call alloc(xnx, 1, imx, 1, jmx-1, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for xnx.')
            call alloc(xny, 1, imx, 1, jmx-1, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for xny.')
            call alloc(xnz, 1, imx, 1, jmx-1, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for xny.')
            call alloc(ynx, 1, imx-1, 1, jmx, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for ynx.')
            call alloc(yny, 1, imx-1, 1, jmx, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for yny.')
            call alloc(ynz, 1, imx-1, 1, jmx, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for ynx.')
            call alloc(znx, 1, imx-1, 1, jmx-1, 1, kmx, &
                    errmsg='Error: Unable to allocate memory for ynx.')
            call alloc(zny, 1, imx-1, 1, jmx-1, 1, kmx, &
                    errmsg='Error: Unable to allocate memory for yny.')
            call alloc(znz, 1, imx-1, 1, jmx-1, 1, kmx, &
                    errmsg='Error: Unable to allocate memory for ynx.')

        end subroutine allocate_memory_normals

        subroutine allocate_memory_ghost_centroids()
            !-----------------------------------------------------------
            ! Allocate memory for centroids of ghost cells
            !-----------------------------------------------------------
            
            implicit none

            call alloc(left_ghost_centroid, 1, jmx-1, 1, kmx-1, 1, 3, &
                    errmsg='Error: Unable to allocate memory for left_ghost_centroid')
            call alloc(right_ghost_centroid, 1, jmx-1, 1, kmx-1, 1, 3, &
                    errmsg='Error: Unable to allocate memory for right_ghost_centroid')
            call alloc(front_ghost_centroid, 1, imx-1, 1, kmx-1, 1, 3, &
                    errmsg='Error: Unable to allocate memory for front_ghost_centroid')
            call alloc(back_ghost_centroid, 1, imx-1, 1, kmx-1, 1, 3, &
                    errmsg='Error: Unable to allocate memory for back_ghost_centroid')
            call alloc(top_ghost_centroid, 1, imx-1, 1, jmx-1, 1, 3, &
                    errmsg='Error: Unable to allocate memory for top_ghost_centroid')
            call alloc(bottom_ghost_centroid, 1, imx-1, 1, jmx-1, 1, 3, &
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
            call dealloc(front_ghost_centroid)
            call dealloc(back_ghost_centroid)
            call dealloc(top_ghost_centroid)
            call dealloc(bottom_ghost_centroid)
    
        end subroutine deallocate_memory

        subroutine normalize_face_normals()
            !-----------------------------------------------------------
            ! Normalize the face normal vectors computed to get unit
            ! vectors
            !-----------------------------------------------------------
            
            implicit none

            xnx(:, :, :) = xnx(:, :, :) / xA(:, :, :)
            xny(:, :, :) = xny(:, :, :) / xA(:, :, :) 
            xnz(:, :, :) = xnz(:, :, :) / xA(:, :, :)

            ynx(:, :, :) = ynx(:, :, :) / yA(:, :, :)
            yny(:, :, :) = yny(:, :, :) / yA(:, :, :)
            ynz(:, :, :) = ynz(:, :, :) / yA(:, :, :)

            znx(:, :, :) = znx(:, :, :) / zA(:, :, :)
            zny(:, :, :) = zny(:, :, :) / zA(:, :, :)
            znz(:, :, :) = znz(:, :, :) / zA(:, :, :)
            
        end subroutine normalize_face_normals

        subroutine compute_face_areas()
            !-----------------------------------------------------------
            ! Compute face areas based on area vectors
            !
            ! The face areas are computed using the face area vectors. 
            ! Prior to using this subroutine, the face area vectors must
            ! computed and placed in the face normal variables. 
            !
            ! Since the area is given by abs(d1 x d2), the areas are
            ! calculated using the normal vectors calculated in 
            ! compute_face_area_vectors, but before normalising them
            !-----------------------------------------------------------
            
            implicit none

            xA(:, :, :) = sqrt((xnx(:, :, :)) ** 2. + (xny(:, :, :)) ** 2. + &
                          (xnz(:, :, :)) ** 2.)

            yA(:, :, :) = sqrt((ynx(:, :, :)) ** 2. + (yny(:, :, :)) ** 2. + &
                          (ynz(:, :, :)) ** 2.)
            
            zA(:, :, :) = sqrt((znx(:, :, :)) ** 2. + (zny(:, :, :)) ** 2. + &
                          (znz(:, :, :)) ** 2.)

        end subroutine compute_face_areas

        subroutine compute_face_area_vectors()
            !-----------------------------------------------------------
            ! Compute face area vectors
            !
            ! The face area vectors denote the face area both in 
            ! magnitude and direction. They are placed in the face 
            ! normal variables for further calculation.
            !
            ! The face normal is given by d1 x d2, where d1 and d2 are
            ! the diagonals of a face
            !-----------------------------------------------------------
            
            implicit none

!           real, dimension(1:imx, 1:jmx-1, 1:kmx-1) :: xd1x, xd1y, xd1z, &
!                                                       xd2x, xd2y, xd2z 
!           real, dimension(1:imx-1, 1:jmx, 1:kmx-1) :: yd1x, yd1y, yd1z, &
!                                                       yd2x, yd2y, yd2z 
!           real, dimension(1:imx-1, 1:jmx-1, 1:kmx) :: zd1x, zd1y, zd1z, &
!                                                       zd2x, zd2y, zd2z 
    
            real :: d1x, d2x, d1y, d2y, d1z, d2z
            integer :: i, j, k

            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx
                d1x = grid_x(i, j+1, k+1) - grid_x(i, j, k)
                d1y = grid_y(i, j+1, k+1) - grid_y(i, j, k)
                d1z = grid_z(i, j+1, k+1) - grid_z(i, j, k)
                d2x = grid_x(i, j, k+1) - grid_x(i, j+1, k)
                d2y = grid_y(i, j, k+1) - grid_y(i, j+1, k)
                d2z = grid_z(i, j, k+1) - grid_z(i, j+1, k)
                xnx(i, j, k) = 0.5 * (d1y*d2z - d1z*d2y)
                xny(i, j, k) = 0.5 * (d1z*d2x - d1x*d2z)
                xnz(i, j, k) = 0.5 * (d1x*d2y - d1y*d2x)
               end do
              end do
             end do

          ! xd1x(:, :, :) = grid_x(1:imx, 2:jmx, 2:kmx) - grid_x(1:imx, 1:jmx-1, 1:kmx-1)
          ! xd1y(:, :, :) = grid_y(1:imx, 2:jmx, 2:kmx) - grid_y(1:imx, 1:jmx-1, 1:kmx-1)
          ! xd1z(:, :, :) = grid_z(1:imx, 2:jmx, 2:kmx) - grid_z(1:imx, 1:jmx-1, 1:kmx-1)
          ! xd2x(:, :, :) = grid_x(1:imx, 1:jmx-1, 2:kmx) - grid_x(1:imx, 2:jmx, 1:kmx-1)
          ! xd2y(:, :, :) = grid_y(1:imx, 1:jmx-1, 2:kmx) - grid_y(1:imx, 2:jmx, 1:kmx-1)
          ! xd2z(:, :, :) = grid_z(1:imx, 1:jmx-1, 2:kmx) - grid_z(1:imx, 2:jmx, 1:kmx-1)
          ! 
          ! xnx(:, :, :) = xd1y(:,:,:)*xd2z(:,:,:) - xd1z(:,:,:)*xd2y(:,:,:)
          ! xny(:, :, :) = xd1z(:,:,:)*xd2x(:,:,:) - xd1x(:,:,:)*xd2z(:,:,:)
          ! xnz(:, :, :) = xd1x(:,:,:)*xd2y(:,:,:) - xd1y(:,:,:)*xd2x(:,:,:)

            do k = 1, kmx - 1
             do j = 1, jmx
              do i = 1, imx - 1
               d1x = grid_x(i+1, j, k+1) - grid_x(i, j, k)
               d1y = grid_y(i+1, j, k+1) - grid_y(i, j, k)
               d1z = grid_z(i+1, j, k+1) - grid_z(i, j, k)
               d2x = grid_x(i+1, j, k) - grid_x(i, j, k+1)
               d2y = grid_y(i+1, j, k) - grid_y(i, j, k+1)
               d2z = grid_z(i+1, j, k) - grid_z(i, j, k+1)
            
               ynx(i, j, k) = 0.5 * (d1y*d2z - d1z*d2y)
               yny(i, j, k) = 0.5 * (d1z*d2x - d1x*d2z)
               ynz(i, j, k) = 0.5 * (d1x*d2y - d1y*d2x)
              end do
             end do
            end do

          ! yd1x(:, :, :) = grid_x(2:imx, 1:jmx, 2:kmx) - grid_x(1:imx-1, 1:jmx, 1:kmx-1)
          ! yd1y(:, :, :) = grid_y(2:imx, 1:jmx, 2:kmx) - grid_y(1:imx-1, 1:jmx, 1:kmx-1)
          ! yd1z(:, :, :) = grid_z(2:imx, 1:jmx, 2:kmx) - grid_z(1:imx-1, 1:jmx, 1:kmx-1)
          ! yd2x(:, :, :) = grid_x(2:imx, 1:jmx, 1:kmx-1) - grid_x(1:imx-1, 1:jmx, 2:kmx)
          ! yd2y(:, :, :) = grid_y(2:imx, 1:jmx, 1:kmx-1) - grid_y(1:imx-1, 1:jmx, 2:kmx)
          ! yd2z(:, :, :) = grid_z(2:imx, 1:jmx, 1:kmx-1) - grid_z(1:imx-1, 1:jmx, 2:kmx)
          ! 
          ! ynx(:, :, :) = yd1y(:,:,:)*yd2z(:,:,:) - yd1z(:,:,:)*yd2y(:,:,:)
          ! yny(:, :, :) = yd1z(:,:,:)*yd2x(:,:,:) - yd1x(:,:,:)*yd2z(:,:,:)
          ! ynz(:, :, :) = yd1x(:,:,:)*yd2y(:,:,:) - yd1y(:,:,:)*yd2x(:,:,:)
            
            do k = 1, kmx
             do j = 1, jmx - 1
              do i = 1, imx - 1
               d1x = grid_x(i+1, j+1, k) - grid_x(i, j, k)
               d1y = grid_y(i+1, j+1, k) - grid_y(i, j, k)
               d1z = grid_z(i+1, j+1, k) - grid_z(i, j, k)
               d2x = grid_x(i, j+1, k) - grid_x(i+1, j, k)
               d2y = grid_y(i, j+1, k) - grid_y(i+1, j, k)
               d2z = grid_z(i, j+1, k) - grid_z(i+1, j, k)
            
               znx(i, j, k) = 0.5 * (d1y*d2z - d1z*d2y)
               zny(i, j, k) = 0.5 * (d1z*d2x - d1x*d2z)
               znz(i, j, k) = 0.5 * (d1x*d2y - d1y*d2x)
              end do
             end do
            end do

          ! zd1x(:, :, :) = grid_x(2:imx, 2:jmx, 1:kmx) - grid_x(1:imx-1, 1:jmx-1, 1:kmx)
          ! zd1y(:, :, :) = grid_y(2:imx, 2:jmx, 1:kmx) - grid_y(1:imx-1, 1:jmx-1, 1:kmx)
          ! zd1z(:, :, :) = grid_z(2:imx, 2:jmx, 1:kmx) - grid_z(1:imx-1, 1:jmx-1, 1:kmx)
          ! zd2x(:, :, :) = grid_x(1:imx-1, 2:jmx, 1:kmx) - grid_x(2:imx, 1:jmx-1, 1:kmx)
          ! zd2y(:, :, :) = grid_y(1:imx-1, 2:jmx, 1:kmx) - grid_y(2:imx, 1:jmx-1, 1:kmx)
          ! zd2z(:, :, :) = grid_z(1:imx-1, 2:jmx, 1:kmx) - grid_z(2:imx, 1:jmx-1, 1:kmx)
          ! 
          ! znx(:, :, :) = zd1y(:,:,:)*zd2z(:,:,:) - zd1z(:,:,:)*zd2y(:,:,:)
          ! zny(:, :, :) = zd1z(:,:,:)*zd2x(:,:,:) - zd1x(:,:,:)*zd2z(:,:,:)
          ! znz(:, :, :) = zd1x(:,:,:)*zd2y(:,:,:) - zd1y(:,:,:)*zd2x(:,:,:)

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
       
        function vol_tetrahedron(p1, p2, p3, p4)
            !-----------------------------------------------------------
            ! Compute the volume of a tetrahedron, given 4 points which
            ! are 1-D arrays
            ! Since we know that the determinant is to be evaluated of 
            ! a 3x3 matrix, we write the expression itself
            !-----------------------------------------------------------

            implicit none
            real, dimension(:), intent(in):: p1, p2, p3, p4
            real, dimension(1:3,1:3) :: A
            real :: vol_tetrahedron

            A(:, 1) = p1 - p4
            A(:, 2) = p2 - p4
            A(:, 3) = p3 - p4

            vol_tetrahedron = A(1,1) * (A(2,2)*A(3,3) - A(2,3)*A(3,2)) + &
                              A(1,2) * (A(2,3)*A(3,1) - A(2,1)*A(3,3)) + &
                              A(1,3) * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
            vol_tetrahedron = vol_tetrahedron / 6.                  
        
        end function vol_tetrahedron
        
        function vol_hexahedron(p_list)
            !-----------------------------------------------------------
            ! Compute the volume of a hexahedron, given a list of points
            ! The points are arranged in a specific order. For the 
            ! element i,j,k, the order of nodes required are:
            ! i,j,k
            ! i+1, j, k
            ! i+1, j+1, k
            ! i, j+1, k
            ! i, j, k+1
            ! i+1, j, k+1
            ! i+1, j+1, k+1
            ! i, j+1, k+1
            !
            ! The hexahedron is to be split into 5 tetrahedrons. 
            ! Source: Split hex into 5 tetrahedron:
            ! No assumptions about planarity seem to be made. All cuts 
            ! were made with a plane containing only 3 vertices at a time.
            ! https://ieeexplore.ieee.org/ieee_pilot/articles/06/ttg2009061587/assets/img/article_1/fig_6/large.gif
            !
            ! The indices of the 5 split tetrahedra can be visualised from
            ! the above link. But since the volume of each tetrahedron 
            ! depends on the determinant calculated, it is IMPERATIVE to 
            ! ensure that a "correct" order is followed for the 4 points. 
            !            
            ! The logic to get the "correct" order is explained as 
            ! follows (Refer wiki article on parallelepiped):
            ! The determinant is taken of a matrix of pi - p4, i = 1, 2, 3. 
            ! Graphically it denotes the sides with p4 as common vertex, 
            ! with direction outward from p4, i.e., directed from 
            ! p4 to pi, i = 1, 2, 3
            ! Hence, if you ensure that  cross(p1-p4, p2-p4) is along 
            ! p3-p4, then the determinant will be positive.
            !
            ! From the above link, a set of 5 tetrahedra was obtained. 
            ! Each tetrahedra has 4 points, and in the function calls 
            ! below, care was taken to ensure that the order is observed
            ! while passing parameters into the vol_tetrahedron function
            !-----------------------------------------------------------

            implicit none
            real, dimension(1:3, 1:8), intent(in) :: p_list
            real :: vol_hexahedron
            
            vol_hexahedron = 0.
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,1), p_list(:,5), &
                                             p_list(:,8), p_list(:,6))
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,7), p_list(:,8), &
                                             p_list(:,6), p_list(:,3))
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,8), p_list(:,4), &
                                             p_list(:,1), p_list(:,3))
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,6), p_list(:,1), &
                                             p_list(:,3), p_list(:,8))
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,1), p_list(:,2), &
                                             p_list(:,6), p_list(:,3))
            
        end function vol_hexahedron
        
        subroutine compute_volumes()
            !-----------------------------------------------------------
            ! Compute the grid cell volumes
            ! Each grid is a hexahedron, whose volume is calculated by
            ! splitting it into 5 tetrahedrons, whose volume is known
            !-----------------------------------------------------------

            implicit none
            integer :: i,j,k
            real, dimension(1:3, 1:8) :: p_list

            do k = 1, kmx - 1
                do j = 1, jmx - 1
                    do i = 1, imx -1
                        p_list(:, :) = 0.
                        p_list(:, 1) = (/ grid_x(i,j,k), grid_y(i,j,k), grid_z(i,j,k) /)
                        p_list(:, 2) = (/ grid_x(i+1,j,k), grid_y(i+1,j,k), grid_z(i+1,j,k) /)
                        p_list(:, 3) = (/ grid_x(i+1,j+1,k), grid_y(i+1,j+1,k), grid_z(i+1,j+1,k) /)
                        p_list(:, 4) = (/ grid_x(i,j+1,k), grid_y(i,j+1,k), grid_z(i,j+1,k) /)
                        p_list(:, 5) = (/ grid_x(i,j,k+1), grid_y(i,j,k+1), grid_z(i,j,k+1) /)
                        p_list(:, 6) = (/ grid_x(i+1,j,k+1), grid_y(i+1,j,k+1), grid_z(i+1,j,k+1) /)
                        p_list(:, 7) = (/ grid_x(i+1,j+1,k+1), grid_y(i+1,j+1,k+1), grid_z(i+1,j+1,k+1) /)
                        p_list(:, 8) = (/ grid_x(i,j+1,k+1), grid_y(i,j+1,k+1), grid_z(i,j+1,k+1) /)
                        volume(i, j, k) = vol_hexahedron(p_list)
                    end do
                end do
            end do
            
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
            real, dimension(3) :: a
            real, dimension(3) :: face_centroid, centroid
            integer :: i, j, k
            
            ! left face ghost cell centroids. i = 1
            do k = 1, kmx - 1
             do j = 1, jmx - 1
                centroid(1) = (grid_x(1, j, k) + grid_x(2, j, k) + &
                        grid_x(2, j+1, k) + grid_x(1, j+1, k) + &
                        grid_x(1, j, k+1) + grid_x(2, j, k+1) + &
                        grid_x(2, j+1, k+1) + grid_x(1, j+1, k+1) &
                        ) * 0.125
                centroid(2) = (grid_y(1, j, k) + grid_y(2, j, k) + &
                        grid_y(2, j+1, k) + grid_y(1, j+1, k) + &
                        grid_y(1, j, k+1) + grid_y(2, j, k+1) + &
                        grid_y(2, j+1, k+1) + grid_y(1, j+1, k+1) &
                        ) * 0.125
                centroid(3) = (grid_z(1, j, k) + grid_z(2, j, k) + &
                        grid_z(2, j+1, k) + grid_z(1, j+1, k) + &
                        grid_z(1, j, k+1) + grid_z(2, j, k+1) + &
                        grid_z(2, j+1, k+1) + grid_z(1, j+1, k+1) &
                        ) * 0.125

                face_centroid(1) = (grid_x(1, j, k) + grid_x(1, j+1, k) + &
                        grid_x(1, j, k+1) + grid_x(1, j+1, k+1) &
                        ) * 0.25
                face_centroid(2) = (grid_y(1, j, k) + grid_y(1, j+1, k) + &
                        grid_y(1, j, k+1) + grid_y(1, j+1, k+1) &
                        ) * 0.25
                face_centroid(3) = (grid_z(1, j, k) + grid_z(1, j+1, k) + &
                        grid_z(1, j, k+1) + grid_z(1, j+1, k+1) &
                        ) * 0.25

                a(1) = centroid(1) - face_centroid(1)
                a(2) = centroid(2) - face_centroid(2)
                a(3) = centroid(3) - face_centroid(3)

                left_ghost_centroid(j, k, 1) = face_centroid(1) + a(1) - &
                                    (2*(a(1)*xnx(1, j, k) + a(2)*xny(1, j, k) + &
                                        a(3)*xnz(1, j, k))*(xnx(1, j, k)) )
                left_ghost_centroid(j, k, 2) = face_centroid(2) + a(2) - &
                                    (2*(a(1)*xnx(1, j, k) + a(2)*xny(1, j, k) + &
                                        a(3)*xnz(1, j, k))*(xny(1, j, k)) )
                left_ghost_centroid(j, k, 3) = face_centroid(3) + a(3) - &
                                    (2*(a(1)*xnx(1, j, k) + a(2)*xny(1, j, k) + &
                                        a(3)*xnz(1, j, k))*(xnz(1, j, k)) )               
             end do
            end do

            ! right face ghost cell centroids. i = imx
            do k = 1, kmx - 1
             do j = 1, jmx - 1
                centroid(1) = (grid_x(imx-1, j, k) + grid_x(imx, j, k) + &
                        grid_x(imx, j+1, k) + grid_x(imx-1, j+1, k) + &
                        grid_x(imx-1, j, k+1) + grid_x(imx, j, k+1) + &
                        grid_x(imx, j+1, k+1) + grid_x(imx-1, j+1, k+1) &
                        ) * 0.125
                centroid(2) = (grid_y(imx-1, j, k) + grid_y(imx, j, k) + &
                        grid_y(imx, j+1, k) + grid_y(imx-1, j+1, k) + &
                        grid_y(imx-1, j, k+1) + grid_y(imx, j, k+1) + &
                        grid_y(imx, j+1, k+1) + grid_y(imx-1, j+1, k+1) &
                        ) * 0.125
                centroid(3) = (grid_z(imx-1, j, k) + grid_z(imx, j, k) + &
                        grid_z(imx, j+1, k) + grid_z(imx-1, j+1, k) + &
                        grid_z(imx-1, j, k+1) + grid_z(imx, j, k+1) + &
                        grid_z(imx, j+1, k+1) + grid_z(imx-1, j+1, k+1) &
                        ) * 0.125

                face_centroid(1) = (grid_x(imx, j, k) + grid_x(imx, j+1, k) + &
                        grid_x(imx, j, k+1) + grid_x(imx, j+1, k+1) &
                        ) * 0.25
                face_centroid(2) = (grid_y(imx, j, k) + grid_y(imx, j+1, k) + &
                        grid_y(imx, j, k+1) + grid_y(imx, j+1, k+1) &
                        ) * 0.25
                face_centroid(3) = (grid_z(imx, j, k) + grid_z(imx, j+1, k) + &
                        grid_z(imx, j, k+1) + grid_z(imx, j+1, k+1) &
                        ) * 0.25

                a(1) = centroid(1) - face_centroid(1)
                a(2) = centroid(2) - face_centroid(2)
                a(3) = centroid(3) - face_centroid(3)

                right_ghost_centroid(j, k, 1) = face_centroid(1) + a(1) - &
                                    (2*(a(1)*xnx(imx, j, k) + a(2)*xny(imx, j, k) + &
                                        a(3)*xnz(imx, j, k))*(xnx(imx, j, k)) )
                right_ghost_centroid(j, k, 2) = face_centroid(2) + a(2) - &
                                    (2*(a(1)*xnx(imx, j, k) + a(2)*xny(imx, j, k) + &
                                        a(3)*xnz(imx, j, k))*(xny(imx, j, k)) )
                right_ghost_centroid(j, k, 3) = face_centroid(3) + a(3) - &
                                    (2*(a(1)*xnx(imx, j, k) + a(2)*xny(imx, j, k) + &
                                        a(3)*xnz(imx, j, k))*(xnz(imx, j, k)) )                
             end do
            end do

            ! front face ghost cell centroids. j = 1
            do k = 1, kmx - 1
             do i = 1, imx - 1
                centroid(1) = (grid_x(i, 1, k) + grid_x(i+1, 1, k) + &
                        grid_x(i+1, 2, k) + grid_x(i, 2, k) + &
                        grid_x(i, 1, k+1) + grid_x(i+1, 1, k+1) + &
                        grid_x(i+1, 2, k+1) + grid_x(i, 2, k+1) &
                        ) * 0.125
                centroid(2) = (grid_y(i, 1, k) + grid_y(i+1, 1, k) + &
                        grid_y(i+1, 2, k) + grid_y(i, 2, k) + &
                        grid_y(i, 1, k+1) + grid_y(i+1, 1, k+1) + &
                        grid_y(i+1, 2, k+1) + grid_y(i, 2, k+1) &
                        ) * 0.125
                centroid(3) = (grid_z(i, 1, k) + grid_z(i+1, 1, k) + &
                        grid_z(i+1, 2, k) + grid_z(i, 2, k) + &
                        grid_z(i, 1, k+1) + grid_z(i+1, 1, k+1) + &
                        grid_z(i+1, 2, k+1) + grid_z(i, 2, k+1) &
                        ) * 0.125

                face_centroid(1) = (grid_x(i, 1, k) + grid_x(i+1, 1, k) + &
                        grid_x(i+1, 1, k+1) + grid_x(i, 1, k+1) &
                        ) * 0.25
                face_centroid(2) = (grid_y(i, 1, k) + grid_y(i+1, 1, k) + &
                        grid_y(i+1, 1, k+1) + grid_y(i, 1, k+1) &
                        ) * 0.25
                face_centroid(3) = (grid_z(i, 1, k) + grid_z(i+1, 1, k) + &
                        grid_z(i+1, 1, k+1) + grid_z(i, 1, k+1) &
                        ) * 0.25

                a(1) = centroid(1) - face_centroid(1)
                a(2) = centroid(2) - face_centroid(2)
                a(3) = centroid(3) - face_centroid(3)

                front_ghost_centroid(i, k, 1) = face_centroid(1) + a(1) - &
                                    (2*(a(1)*ynx(i, 1, k) + a(2)*yny(i, 1, k) + &
                                        a(3)*ynz(i, 1, k))*(ynx(i, 1, k)) )
                front_ghost_centroid(i, k, 2) = face_centroid(2) + a(2) - &
                                    (2*(a(1)*ynx(i, 1, k) + a(2)*yny(i, 1, k) + &
                                        a(3)*ynz(i, 1, k))*(yny(i, 1, k)) )
                front_ghost_centroid(i, k, 3) = face_centroid(3) + a(3) - &
                                    (2*(a(1)*ynx(i, 1, k) + a(2)*yny(i, 1, k) + &
                                        a(3)*ynz(i, 1, k))*(ynz(i, 1, k)) )                
             end do
            end do

            ! back face ghost cell centroids. j = jmx
            do k = 1, kmx - 1
             do i = 1, imx - 1
                centroid(1) = (grid_x(i, jmx-1, k) + grid_x(i+1, jmx-1, k) + &
                        grid_x(i+1, jmx, k) + grid_x(i, jmx, k) + &
                        grid_x(i, jmx-1, k+1) + grid_x(i+1, jmx-1, k+1) + &
                        grid_x(i+1, jmx, k+1) + grid_x(i, jmx, k+1) &
                        ) * 0.125
                centroid(2) = (grid_y(i, jmx-1, k) + grid_y(i+1, jmx-1, k) + &
                        grid_y(i+1, jmx, k) + grid_y(i, jmx, k) + &
                        grid_y(i, jmx-1, k+1) + grid_y(i+1, jmx-1, k+1) + &
                        grid_y(i+1, jmx, k+1) + grid_y(i, jmx, k+1) &
                        ) * 0.125
                centroid(3) = (grid_z(i, jmx-1, k) + grid_z(i+1, jmx-1, k) + &
                        grid_z(i+1, jmx, k) + grid_z(i, jmx, k) + &
                        grid_z(i, jmx-1, k+1) + grid_z(i+1, jmx-1, k+1) + &
                        grid_z(i+1, jmx, k+1) + grid_z(i, jmx, k+1) &
                        ) * 0.125

                face_centroid(1) = (grid_x(i, jmx, k) + grid_x(i+1, jmx, k) + &
                        grid_x(i+1, jmx, k+1) + grid_x(i, jmx, k+1) &
                        ) * 0.25
                face_centroid(2) = (grid_y(i, jmx, k) + grid_y(i+1, jmx, k) + &
                        grid_y(i+1, jmx, k+1) + grid_y(i, jmx, k+1) &
                        ) * 0.25
                face_centroid(3) = (grid_z(i, jmx, k) + grid_z(i+1, jmx, k) + &
                        grid_z(i+1, jmx, k+1) + grid_z(i, jmx, k+1) &
                        ) * 0.25

                a(1) = centroid(1) - face_centroid(1)
                a(2) = centroid(2) - face_centroid(2)
                a(3) = centroid(3) - face_centroid(3)

                back_ghost_centroid(i, k, 1) = face_centroid(1) + a(1) - &
                                    (2*(a(1)*ynx(i, jmx, k) + a(2)*yny(i, jmx, k) + &
                                        a(3)*ynz(i, jmx, k))*(ynx(i, jmx, k)) )
                back_ghost_centroid(i, k, 2) = face_centroid(2) + a(2) - &
                                    (2*(a(1)*ynx(i, jmx, k) + a(2)*yny(i, jmx, k) + &
                                        a(3)*ynz(i, jmx, k))*(yny(i, jmx, k)) )
                back_ghost_centroid(i, k, 3) = face_centroid(3) + a(3) - &
                                    (2*(a(1)*ynx(i, jmx, k) + a(2)*yny(i, jmx, k) + &
                                        a(3)*ynz(i, jmx, k))*(ynz(i, jmx, k)) )                
             end do
            end do

            ! bottom face ghost cell centroids. k = 1
            do j = 1, jmx - 1
             do i = 1, imx - 1
                centroid(1) = (grid_x(i, j, 1) + grid_x(i+1, j, 1) + &
                        grid_x(i+1, j+1, 1) + grid_x(i, j, 1) + &
                        grid_x(i, j, 2) + grid_x(i+1, j, 2) + &
                        grid_x(i+1, j+1, 2) + grid_x(i, j+1, 2) &
                        ) * 0.125
                centroid(2) = (grid_y(i, j, 1) + grid_y(i+1, j, 1) + &
                        grid_y(i+1, j+1, 1) + grid_y(i, j, 1) + &
                        grid_y(i, j, 2) + grid_y(i+1, j, 2) + &
                        grid_y(i+1, j+1, 2) + grid_y(i, j+1, 2) &
                        ) * 0.125
                centroid(3) = (grid_z(i, j, 1) + grid_z(i+1, j, 1) + &
                        grid_z(i+1, j+1, 1) + grid_z(i, j, 1) + &
                        grid_z(i, j, 2) + grid_z(i+1, j, 2) + &
                        grid_z(i+1, j+1, 2) + grid_z(i, j+1, 2) &
                        ) * 0.125

                face_centroid(1) = (grid_x(i, j, 1) + grid_x(i+1, j, 1) + &
                        grid_x(i+1, j+1, 1) + grid_x(i, j+1, 1) &
                        ) * 0.25
                face_centroid(2) = (grid_y(i, j, 1) + grid_y(i+1, j, 1) + &
                        grid_y(i+1, j+1, 1) + grid_y(i, j+1, 1) &
                        ) * 0.25
                face_centroid(3) = (grid_z(i, j, 1) + grid_z(i+1, j, 1) + &
                        grid_z(i+1, j+1, 1) + grid_z(i, j+1, 1) &
                        ) * 0.25

                a(1) = centroid(1) - face_centroid(1)
                a(2) = centroid(2) - face_centroid(2)
                a(3) = centroid(3) - face_centroid(3)

                bottom_ghost_centroid(i, j, 1) = face_centroid(1) + a(1) - &
                                    (2*(a(1)*znx(i, j, 1) + a(2)*zny(i, j, 1) + &
                                        a(3)*znz(i, j, 1))*(znx(i, j, 1)) )
                bottom_ghost_centroid(i, j, 2) = face_centroid(2) + a(2) - &
                                    (2*(a(1)*znx(i, j, 1) + a(2)*zny(i, j, 1) + &
                                        a(3)*znz(i, j, 1))*(zny(i, j, 1)) )
                bottom_ghost_centroid(i, j, 3) = face_centroid(3) + a(3) - &
                                    (2*(a(1)*znx(i, j, 1) + a(2)*zny(i, j, 1) + &
                                        a(3)*znz(i, j, 1))*(znz(i, j, 1)) )          
             end do
            end do

            ! top face ghost cell centroids. k = kmx
            do j = 1, jmx - 1
             do i = 1, imx - 1
                centroid(1) = (grid_x(i, j, kmx-1) + grid_x(i+1, j, kmx-1) + &
                        grid_x(i+1, j+1, kmx-1) + grid_x(i, j, kmx-1) + &
                        grid_x(i, j, kmx) + grid_x(i+1, j, kmx) + &
                        grid_x(i+1, j+1, kmx) + grid_x(i, j+1, kmx) &
                        ) * 0.125
                centroid(2) = (grid_y(i, j, kmx-1) + grid_y(i+1, j, kmx-1) + &
                        grid_y(i+1, j+1, kmx-1) + grid_y(i, j, kmx-1) + &
                        grid_y(i, j, kmx) + grid_y(i+1, j, kmx) + &
                        grid_y(i+1, j+1, kmx) + grid_y(i, j+1, kmx) &
                        ) * 0.125
                centroid(3) = (grid_z(i, j, kmx-1) + grid_z(i+1, j, kmx-1) + &
                        grid_z(i+1, j+1, kmx-1) + grid_z(i, j, kmx-1) + &
                        grid_z(i, j, kmx) + grid_z(i+1, j, kmx) + &
                        grid_z(i+1, j+1, kmx) + grid_z(i, j+1, kmx) &
                        ) * 0.125

                face_centroid(1) = (grid_x(i, j, kmx) + grid_x(i+1, j, kmx) + &
                        grid_x(i+1, j+1, kmx) + grid_x(i, j+1, kmx) &
                        ) * 0.25
                face_centroid(2) = (grid_y(i, j, kmx) + grid_y(i+1, j, kmx) + &
                        grid_y(i+1, j+1, kmx) + grid_y(i, j+1, kmx) &
                        ) * 0.25
                face_centroid(3) = (grid_z(i, j, kmx) + grid_z(i+1, j, kmx) + &
                        grid_z(i+1, j+1, kmx) + grid_z(i, j+1, kmx) &
                        ) * 0.25

                a(1) = centroid(1) - face_centroid(1)
                a(2) = centroid(2) - face_centroid(2)
                a(3) = centroid(3) - face_centroid(3)

                top_ghost_centroid(i, j, 1) = face_centroid(1) + a(1) - &
                                    (2*(a(1)*znx(i, j, kmx) + a(2)*zny(i, j, kmx) + &
                                        a(3)*znz(i, j, kmx))*(znx(i, j, kmx)) )
                top_ghost_centroid(i, j, 2) = face_centroid(2) + a(2) - &
                                    (2*(a(1)*znx(i, j, kmx) + a(2)*zny(i, j, kmx) + &
                                        a(3)*znz(i, j, kmx))*(zny(i, j, kmx)) )
                top_ghost_centroid(i, j, 3) = face_centroid(3) + a(3) - &
                                    (2*(a(1)*znx(i, j, kmx) + a(2)*zny(i, j, kmx) + &
                                        a(3)*znz(i, j, kmx))*(znz(i, j, kmx)) )              
             end do
            end do

!           open(24, file='left_ghost.txt')

!           do k = 1, kmx-1
!            do j = 1, jmx-1
!               write (24, *) left_ghost_centroid(j, k, 1), left_ghost_centroid(j, k, 2), &
!                             left_ghost_centroid(j, k, 3)
!               write (24, *) right_ghost_centroid(j, k, 1), right_ghost_centroid(j, k, 2), &
!                             right_ghost_centroid(j, k, 3)
!            end do
!           end do
!           
!           do k = 1, kmx-1
!            do i = 1, imx-1
!               write (24, *) front_ghost_centroid(i, k, 1), front_ghost_centroid(i, k, 2), &
!                             front_ghost_centroid(i, k, 3)
!               write (24, *) back_ghost_centroid(i, k, 1), back_ghost_centroid(i, k, 2), &
!                             back_ghost_centroid(i, k, 3)
!            end do
!           end do

!           do j = 1, jmx-1
!            do i = 1, imx-1
!               write (24, *) top_ghost_centroid(i, j, 1), top_ghost_centroid(i, j, 2), &
!                             top_ghost_centroid(i, j, 3)
!               write (24, *) bottom_ghost_centroid(i, j, 1), bottom_ghost_centroid(i, j, 2), &
!                             bottom_ghost_centroid(i, j, 3)
!            end do
!           end do

!           close(24)

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
