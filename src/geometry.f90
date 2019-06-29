    !< The geometry module calculates various geometrical quantities like 
    !< face-normals, face-areas and cell-volumes to be used in computations. 
module geometry
    !< The geometry module calculates various geometrical quantities like 
    !< face-normals, face-areas and cell-volumes to be used in computations. 
    !-------------------------------------------------------------------
#include "error.inc"
    use global_vars, only : imx
    use global_vars, only : jmx
    use global_vars, only : kmx
    use global_vars, only : grid_x
    use global_vars, only : grid_y
    use global_vars, only : grid_z

    use global_vars, only : xn           !face unit norm x
    use global_vars, only : yn           !face unit norm y
    use global_vars, only : zn           !face unit norm z
    use global_vars, only : xnx, xny, xnz !face unit normal x
    use global_vars, only : ynx, yny, ynz !face unit normal y
    use global_vars, only : znx, zny, znz !face unit normal z
    use global_vars, only : xA, yA, zA    !face area
    use global_vars, only : volume
    use global_vars, only :   left_ghost_centroid
    use global_vars, only :  right_ghost_centroid
    use global_vars, only :  front_ghost_centroid
    use global_vars, only :   back_ghost_centroid
    use global_vars, only :    top_ghost_centroid
    use global_vars, only : bottom_ghost_centroid
    use global_vars, only : imin_id
    use global_vars, only : imax_id
    use global_vars, only : jmin_id
    use global_vars, only : jmax_id
    use global_vars, only : kmin_id
    use global_vars, only : kmax_id
    use global_vars, only : process_id
    
    use utils, only: alloc, dealloc, dmsg

    implicit none
    private

    real, dimension(:,:,:,:), allocatable, public:: CellCenter
    !< Store Cell-center location 

    ! Public methods
    public :: setup_geometry
    public :: destroy_geometry

    contains

        subroutine allocate_memory_volumes()
            !< Allocate memory for the volume variables.
            !-----------------------------------------------------------
            
            implicit none

            call alloc(volume, -2, imx+2, -2, jmx+2, -2, kmx+2, &
                    errmsg='Error: Unable to allocate memory for volume.')

        end subroutine allocate_memory_volumes

        subroutine allocate_memory_areas()
            !< Allocate memory for the area variables.
            !-----------------------------------------------------------
            
            implicit none

            call alloc(xA, -2, imx+3, -2, jmx+2, -2, kmx+2, &
                    errmsg='Error: Unable to allocate memory for xA.')
            call alloc(yA, -2, imx+2, -2, jmx+3, -2, kmx+2, &
                    errmsg='Error: Unable to allocate memory for yA.')
            call alloc(zA, -2, imx+2, -2, jmx+2, -2, kmx+3, &
                    errmsg='Error: Unable to allocate memory for yA.')

        end subroutine allocate_memory_areas

        subroutine allocate_memory_normals()
            !< Allocate memory for the normal variables.
            !-----------------------------------------------------------
                        
            implicit none

            call alloc(xn, -2, imx+3, -2, jmx+2, -2, kmx+2, 1,3, &
                    errmsg='Error: Unable to allocate memory for xnx.')
            call alloc(yn, -2, imx+2, -2, jmx+3, -2, kmx+2, 1,3, &
                    errmsg='Error: Unable to allocate memory for ynx.')
            call alloc(zn, -2, imx+2, -2, jmx+2, -2, kmx+3, 1,3, &
                    errmsg='Error: Unable to allocate memory for ynx.')

            xnx(-2:imx+3,-2:jmx+2,-2:kmx+2) => xn(:,:,:,1)
            xny(-2:imx+3,-2:jmx+2,-2:kmx+2) => xn(:,:,:,2)
            xnz(-2:imx+3,-2:jmx+2,-2:kmx+2) => xn(:,:,:,3)

            ynx(-2:imx+2,-2:jmx+3,-2:kmx+2) => yn(:,:,:,1)
            yny(-2:imx+2,-2:jmx+3,-2:kmx+2) => yn(:,:,:,2)
            ynz(-2:imx+2,-2:jmx+3,-2:kmx+2) => yn(:,:,:,3)

            znx(-2:imx+2,-2:jmx+2,-2:kmx+3) => zn(:,:,:,1)
            zny(-2:imx+2,-2:jmx+2,-2:kmx+3) => zn(:,:,:,2)
            znz(-2:imx+2,-2:jmx+2,-2:kmx+3) => zn(:,:,:,3)

        end subroutine allocate_memory_normals

        subroutine allocate_memory_ghost_centroids()
            !< Allocate memory for centroids of ghost cells
            !-----------------------------------------------------------
            
            implicit none

            call alloc(CellCenter, -2, imx+2, -2, jmx+2, -2, kmx+2, 1, 3, &
                    errmsg='Error: Unable to allocate memory for volume.')

!
!            call alloc(left_ghost_centroid, 1, jmx-1, 1, kmx-1, 1, 3, &
!                    errmsg='Error: Unable to allocate memory for left_ghost_centroid')
!            call alloc(right_ghost_centroid, 1, jmx-1, 1, kmx-1, 1, 3, &
!                    errmsg='Error: Unable to allocate memory for right_ghost_centroid')
!            call alloc(front_ghost_centroid, 1, imx-1, 1, kmx-1, 1, 3, &
!                    errmsg='Error: Unable to allocate memory for front_ghost_centroid')
!            call alloc(back_ghost_centroid, 1, imx-1, 1, kmx-1, 1, 3, &
!                    errmsg='Error: Unable to allocate memory for back_ghost_centroid')
!            call alloc(top_ghost_centroid, 1, imx-1, 1, jmx-1, 1, 3, &
!                    errmsg='Error: Unable to allocate memory for top_ghost_centroid')
!            call alloc(bottom_ghost_centroid, 1, imx-1, 1, jmx-1, 1, 3, &
!                    errmsg='Error: Unable to allocate memory for bottom_ghost_centroid')

        end subroutine allocate_memory_ghost_centroids

        subroutine allocate_memory()
            !< Allocate memory for the required variables.
            !-----------------------------------------------------------
            
            implicit none

            call dmsg(1, 'geometry', 'allocate_memory')

            call allocate_memory_normals()
            call allocate_memory_areas()
            call allocate_memory_volumes()
            call allocate_memory_ghost_centroids()

        end subroutine allocate_memory

        subroutine deallocate_memory()
          !< Deallocate the memoery used by the geometry variables

            implicit none

            call dmsg(1, 'geometry', 'deallocate_memory')

            call dealloc(xn)
            call dealloc(yn)
            call dealloc(zn)
            call dealloc(xA)
            call dealloc(yA)
            call dealloc(zA)
            call dealloc(volume)
            call dealloc(CellCenter)
!            call dealloc(left_ghost_centroid)
!            call dealloc(right_ghost_centroid)
!            call dealloc(front_ghost_centroid)
!            call dealloc(back_ghost_centroid)
!            call dealloc(top_ghost_centroid)
!            call dealloc(bottom_ghost_centroid)
    
        end subroutine deallocate_memory

        subroutine normalize_face_normals()
            !< Normalize the face normal vectors computed to get unit
            !< vectors
            !-----------------------------------------------------------
            
            implicit none
            integer :: i,j,k

            do k = -2,kmx+2
              do j = -2,jmx+2
                do i = -2,imx+3
                  if(xA(i,j,k)/=0.) then
                    xnx(i,j,k) = xnx(i,j,k) / xA(i,j,k)
                    xny(i,j,k) = xny(i,j,k) / xA(i,j,k)
                    xnz(i,j,k) = xnz(i,j,k) / xA(i,j,k)
                  end if
                end do
              end do
            end do

            do k = -2,kmx+2
              do j = -2,jmx+3
                do i = -2,imx+2
                  if(yA(i,j,k)/=0.) then
                    ynx(i,j,k) = ynx(i,j,k) / yA(i,j,k)
                    yny(i,j,k) = yny(i,j,k) / yA(i,j,k)
                    ynz(i,j,k) = ynz(i,j,k) / yA(i,j,k)
                  end if
                end do
              end do
            end do

            do k = -2,kmx+3
              do j = -2,jmx+2
                do i = -2,imx+2
                  if(zA(i,j,k)/=0.) then
                    znx(i,j,k) = znx(i,j,k) / zA(i,j,k)
                    zny(i,j,k) = zny(i,j,k) / zA(i,j,k)
                    znz(i,j,k) = znz(i,j,k) / zA(i,j,k)
                  end if
                end do
              end do
            end do

            ! pole boundary condition
            if(imin_id==-7) then
              xn( 1,:,:,:)=xn(2,:,:,:)
              xn( 0,:,:,:)=xn(2,:,:,:)
              xn(-1,:,:,:)=xn(2,:,:,:)
              xn(-2,:,:,:)=xn(2,:,:,:)
            end if

            if(imax_id==-7) then
              xn(imx+0,:,:,:)=xn(imx-1,:,:,:)
              xn(imx+1,:,:,:)=xn(imx-1,:,:,:)
              xn(imx+2,:,:,:)=xn(imx-1,:,:,:)
              xn(imx+3,:,:,:)=xn(imx-1,:,:,:)
            end if

            if(jmin_id==-7) then
              yn(:, 1,:,:)=yn(:,2,:,:)
              yn(:, 0,:,:)=yn(:,2,:,:)
              yn(:,-1,:,:)=yn(:,2,:,:)
              yn(:,-2,:,:)=yn(:,2,:,:)
            end if

            if(jmax_id==-7) then
              yn(:,jmx+0,:,:)=yn(:,jmx-1,:,:)
              yn(:,jmx+1,:,:)=yn(:,jmx-1,:,:)
              yn(:,jmx+2,:,:)=yn(:,jmx-1,:,:)
              yn(:,jmx+3,:,:)=yn(:,jmx-1,:,:)
            end if

            if(kmin_id==-7) then
              zn(:,:, 1,:)=zn(:,:,2,:)
              zn(:,:, 0,:)=zn(:,:,2,:)
              zn(:,:,-1,:)=zn(:,:,2,:)
              zn(:,:,-2,:)=zn(:,:,2,:)
            end if

            if(kmax_id==-7) then
              zn(:,:,kmx+0,:)=zn(:,:,kmx-1,:)
              zn(:,:,kmx+1,:)=zn(:,:,kmx-1,:)
              zn(:,:,kmx+2,:)=zn(:,:,kmx-1,:)
              zn(:,:,kmx+3,:)=zn(:,:,kmx-1,:)
            end if

            
        end subroutine normalize_face_normals

        subroutine compute_face_areas()
            !< Compute face areas based on area vectors
            !<
            !< The face areas are computed using the face area vectors. 
            !< Prior to using this subroutine, the face area vectors must
            !< computed and placed in the face normal variables. 
            !<
            !< Since the area is given by abs(d1 x d2), the areas are
            !< calculated using the normal vectors calculated in 
            !< compute_face_area_vectors, but before normalising them
            !-----------------------------------------------------------
            
            implicit none

            xA(:, :, :) = sqrt((xnx(:, :, :)) ** 2. + (xny(:, :, :)) ** 2. + &
                          (xnz(:, :, :)) ** 2.)

            yA(:, :, :) = sqrt((ynx(:, :, :)) ** 2. + (yny(:, :, :)) ** 2. + &
                          (ynz(:, :, :)) ** 2.)
            
            zA(:, :, :) = sqrt((znx(:, :, :)) ** 2. + (zny(:, :, :)) ** 2. + &
                          (znz(:, :, :)) ** 2.)

            ! Pole boundary conditions
            ! making sure face area is exactly equal zero
            if(imin_id==-7) xA(-2:1     ,  :,  :)=0.
            if(imax_id==-7) xA(imx:imx+3,  :,  :)=0.
            if(jmin_id==-7) yA(  :,     -2:1,  :)=0.
            if(jmax_id==-7) yA(  :,jmx:jmx+3,  :)=0.
            if(kmin_id==-7) zA(  :,  :,     -2:1)=0.
            if(kmax_id==-7) zA(  :,  :,kmx:kmx+3)=0.

        end subroutine compute_face_areas

        subroutine compute_face_area_vectors()
            !< Compute face area vectors
            !<
            !< The face area vectors denote the face area both in 
            !< magnitude and direction. They are placed in the face 
            !< normal variables for further calculation.
            !<
            !< The face normal is given by d1 x d2, where d1 and d2 are
            !< the diagonals of a face
            !-----------------------------------------------------------
            
            implicit none

    
            real :: d1x, d2x, d1y, d2y, d1z, d2z
            integer :: i, j, k

            do k = -2, kmx+2
             do j = -2, jmx+2
              do i = -2, imx+3
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


            do k = -2, kmx+2
             do j = -2, jmx+3
              do i = -2, imx+2
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

            
            do k = -2, kmx+3
             do j = -2, jmx+2
              do i = -2, imx+2
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


        end subroutine compute_face_area_vectors

        subroutine compute_face_areas_and_normals()
            !< Compute the face areas and normals
            !<
            !< This is the 2-dimensional version. In this case, the face 
            !< areas default to edge lengths.
            !-----------------------------------------------------------

            implicit none

            call compute_face_area_vectors()
            call compute_face_areas()
            call normalize_face_normals()
        
        end subroutine compute_face_areas_and_normals
       
        function vol_tetrahedron(p1, p2, p3, p4)
            !< Compute the volume of a tetrahedron, given 4 points which
            !< are 1-D arrays
            !< Since we know that the determinant is to be evaluated of 
            !< a 3x3 matrix, we write the expression itself
            !-----------------------------------------------------------

            implicit none
            real, dimension(:), intent(in):: p1, p2, p3, p4
            real, dimension(1:3,1:3) :: A
            real :: vol_tetrahedron

            A(:, 1) = p1 - p4
            A(:, 2) = p2 - p4
            A(:, 3) = p3 - p4

            !vol_tetrahedron = A(1,1) * (A(2,2)*A(3,3) - A(2,3)*A(3,2)) + &
            !                  A(1,2) * (A(2,3)*A(3,1) - A(2,1)*A(3,3)) + &
            !                  A(1,3) * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
            vol_tetrahedron = (p4(1) - p1(1))*((p2(2)-p1(2))*(p3(3)-p1(3)) - (p2(3)-p1(3))*(p3(2)-p1(2))) &
                            + (p4(2) - p1(2))*((p2(3)-p1(3))*(p3(1)-p1(1)) - (p2(1)-p1(1))*(p3(3)-p1(3))) &
                            + (p4(3) - p1(3))*((p2(1)-p1(1))*(p3(2)-p1(2)) - (p2(2)-p1(2))*(p3(1)-p1(1)))
            vol_tetrahedron = -vol_tetrahedron / 6.                  
        
        end function vol_tetrahedron
        
        function vol_hexahedron(p_list)
            !< Compute the volume of a hexahedron, given a list of points
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
            real :: vol_hexahedron1
            
            vol_hexahedron1 = 0.
            vol_hexahedron1 = vol_hexahedron1 + &
                             vol_tetrahedron(p_list(:,1), p_list(:,5), &
                                             p_list(:,8), p_list(:,6))
            vol_hexahedron1 = vol_hexahedron1 + &
                             vol_tetrahedron(p_list(:,7), p_list(:,8), &
                                             p_list(:,6), p_list(:,3))
            vol_hexahedron1 = vol_hexahedron1 + &
                             vol_tetrahedron(p_list(:,8), p_list(:,4), &
                                             p_list(:,1), p_list(:,3))
            vol_hexahedron1 = vol_hexahedron1 + &
                             vol_tetrahedron(p_list(:,6), p_list(:,1), &
                                             p_list(:,3), p_list(:,8))
            vol_hexahedron1 = vol_hexahedron1 + &
                             vol_tetrahedron(p_list(:,1), p_list(:,2), &
                                             p_list(:,6), p_list(:,3))
            vol_hexahedron = 0.
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,2), p_list(:,6), &
                                             p_list(:,5), p_list(:,7))
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,8), p_list(:,5), &
                                             p_list(:,7), p_list(:,4))
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,5), p_list(:,1), &
                                             p_list(:,2), p_list(:,4))
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,7), p_list(:,2), &
                                             p_list(:,4), p_list(:,5))
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,2), p_list(:,3), &
                                             p_list(:,7), p_list(:,4))

            vol_hexahedron = max(vol_hexahedron,vol_hexahedron1)
            
        end function vol_hexahedron
        
        subroutine compute_volumes()
            !< Compute the grid cell volumes
            !< Each grid is a hexahedron, whose volume is calculated by
            !< splitting it into 5 tetrahedrons, whose volume is known
            !-----------------------------------------------------------

            implicit none
            integer :: i,j,k
            real, dimension(1:3, 1:8) :: p_list

            volume=1.
            do k = 0, kmx+0
                do j = 0, jmx+0
                    do i = 0, imx+0
                        p_list(:, :) = 0.
                        p_list(:, 1) = (/ grid_x(i,j,k), grid_y(i,j,k), grid_z(i,j,k) /)
                        p_list(:, 2) = (/ grid_x(i+1,j,k), grid_y(i+1,j,k), grid_z(i+1,j,k) /)
                        p_list(:, 3) = (/ grid_x(i+1,j+1,k), grid_y(i+1,j+1,k), grid_z(i+1,j+1,k) /)
                        p_list(:, 4) = (/ grid_x(i,j+1,k), grid_y(i,j+1,k), grid_z(i,j+1,k) /)
                        p_list(:, 5) = (/ grid_x(i,j,k+1), grid_y(i,j,k+1), grid_z(i,j,k+1) /)
                        p_list(:, 6) = (/ grid_x(i+1,j,k+1), grid_y(i+1,j,k+1), grid_z(i+1,j,k+1) /)
                        p_list(:, 7) = (/ grid_x(i+1,j+1,k+1), grid_y(i+1,j+1,k+1), grid_z(i+1,j+1,k+1) /)
                        p_list(:, 8) = (/ grid_x(i,j+1,k+1), grid_y(i,j+1,k+1), grid_z(i,j+1,k+1) /)
                        volume(i, j, k) = (vol_hexahedron(p_list))
                        if(volume(i,j,k)<=0.0) then
                          if(i==0 .or. i==imx .or. j==0 .or. j==jmx .or. k==0 .or. k==kmx) then
                            !print*, "Ghost Cell volume negative"
                            volume(i, j, k) = abs(vol_hexahedron(p_list))
                          else
                            print*, process_id, i,j,k
                            print*, "negative volume :", (vol_hexahedron(p_list))
                            STOP
                          end if
                        end if
                    end do
                end do
            end do
            if(any(volume==0.0))then
              Fatal_error
            end if
            if(any((volume)<0.0))then
              Fatal_error
            end if
            
        end subroutine compute_volumes

        subroutine compute_geometric_parameters()
            !< Compute the geometric parameters based on the grid points
            !<
            !< The geometric parameters include the face normals and 
            !< areas and the cell volumes.
            !-----------------------------------------------------------
            
            implicit none

            call dmsg(1, 'geometry', 'compute_geometric_parameters')

            call compute_face_areas_and_normals()
            call compute_volumes()

        end subroutine compute_geometric_parameters

        subroutine compute_ghost_cell_centroid()
          !< Compute cell center of all cell including ghost cells
          implicit none
          integer :: i,j,k

          do k = -2, kmx+2
            do j = -2, jmx+2
              do i = -2, imx+2
                CellCenter(i,j,k,1) = 0.125 * ( grid_x(i  ,j  ,k  ) &
                                              + grid_x(i+1,j  ,k  ) &
                                              + grid_x(i+1,j+1,k  ) &
                                              + grid_x(i+1,j+1,k+1) &
                                              + grid_x(i+1,j  ,k+1) &
                                              + grid_x(i  ,j+1,k  ) &
                                              + grid_x(i  ,j+1,k+1) &
                                              + grid_x(i  ,j  ,k+1) &
                                              )

                CellCenter(i,j,k,2) = 0.125 * ( grid_y(i  ,j  ,k  ) &
                                              + grid_y(i+1,j  ,k  ) &
                                              + grid_y(i+1,j+1,k  ) &
                                              + grid_y(i+1,j+1,k+1) &
                                              + grid_y(i+1,j  ,k+1) &
                                              + grid_y(i  ,j+1,k  ) &
                                              + grid_y(i  ,j+1,k+1) &
                                              + grid_y(i  ,j  ,k+1) &
                                              )

                CellCenter(i,j,k,3) = 0.125 * ( grid_z(i  ,j  ,k  ) &
                                              + grid_z(i+1,j  ,k  ) &
                                              + grid_z(i+1,j+1,k  ) &
                                              + grid_z(i+1,j+1,k+1) &
                                              + grid_z(i+1,j  ,k+1) &
                                              + grid_z(i  ,j+1,k  ) &
                                              + grid_z(i  ,j+1,k+1) &
                                              + grid_z(i  ,j  ,k+1) &
                                              )
              end do
            end do
          end do

        end subroutine compute_ghost_cell_centroid

        subroutine setup_geometry()
            !< Make the geometry module useful
            !<
            !< Allocates memory to the variables and initializes them.
            !-----------------------------------------------------------

            implicit none

            call dmsg(1, 'geometry', 'setup_geometry')

            call allocate_memory()
            call compute_geometric_parameters()
            call compute_ghost_cell_centroid()

        end subroutine setup_geometry

        subroutine destroy_geometry()
          !< Nullify all the face normal pionter 

            implicit none
            
            call dmsg(1, 'geometry', 'destroy_geometry')

            nullify(xnx)
            nullify(xny)
            nullify(xnz)
            nullify(ynx)
            nullify(yny)
            nullify(ynz)
            nullify(znx)
            nullify(zny)
            nullify(znz)
            call deallocate_memory()

        end subroutine destroy_geometry

end module geometry
