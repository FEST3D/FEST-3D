    !< The geometry module calculates various geometrical quantities like 
    !< face-normals, face-areas and cell-volumes to be used in computations. 
module geometry
    !< The geometry module calculates various geometrical quantities like 
    !< face-normals, face-areas and cell-volumes to be used in computations. 
    !-------------------------------------------------------------------
#include "error.inc"
#include "debug.h"
!    use global_vars, only : imx
!    use global_vars, only : jmx
!    use global_vars, only : kmx
    !use grid, only : nodetype, extent
    use vartypes

    use global_vars, only : xn           !face unit norm x
    use global_vars, only : yn           !face unit norm y
    use global_vars, only : zn           !face unit norm z
    use global_vars, only : xnx, xny, xnz !face unit normal x
    use global_vars, only : ynx, yny, ynz !face unit normal y
    use global_vars, only : znx, zny, znz !face unit normal z
    use global_vars, only : xA, yA, zA    !face area
    use global_vars, only : volume
!    use global_vars, only :   left_ghost_centroid
!    use global_vars, only :  right_ghost_centroid
!    use global_vars, only :  front_ghost_centroid
!    use global_vars, only :   back_ghost_centroid
!    use global_vars, only :    top_ghost_centroid
!    use global_vars, only : bottom_ghost_centroid
    use global_vars, only : imin_id
    use global_vars, only : imax_id
    use global_vars, only : jmin_id
    use global_vars, only : jmax_id
    use global_vars, only : kmin_id
    use global_vars, only : kmax_id
    use global_vars, only : process_id
    
    use utils, only: alloc

    implicit none
    private
      real, dimension(:,:,:,:), allocatable, public:: CellCenter
      !< Store Cell-center location 
      integer :: imx, jmx, kmx

    ! Public methods
    public :: setup_geometry
!    public :: destroy_geometry

    contains

!        subroutine allocate_memory_volumes(volume)
!            !< Allocate memory for the volume variables.
!            !-----------------------------------------------------------
!            
!            implicit none
!            real, dimension(:,:,:), allocatable, intent(out) :: volume
!
!            call alloc(volume, -2, imx+2, -2, jmx+2, -2, kmx+2, &
!                    errmsg='Error: Unable to allocate memory for volume.')
!
!        end subroutine allocate_memory_volumes
!
!        subroutine allocate_memory_areas()
!            !< Allocate memory for the area variables.
!            !-----------------------------------------------------------
!            
!            implicit none
!
!            call alloc(xA, -2, imx+3, -2, jmx+2, -2, kmx+2, &
!                    errmsg='Error: Unable to allocate memory for xA.')
!            call alloc(yA, -2, imx+2, -2, jmx+3, -2, kmx+2, &
!                    errmsg='Error: Unable to allocate memory for yA.')
!            call alloc(zA, -2, imx+2, -2, jmx+2, -2, kmx+3, &
!                    errmsg='Error: Unable to allocate memory for yA.')
!
!        end subroutine allocate_memory_areas
!
!        subroutine allocate_memory_normals()
!            !< Allocate memory for the normal variables.
!            !-----------------------------------------------------------
!                        
!            implicit none
!
!            call alloc(xn, -2, imx+3, -2, jmx+2, -2, kmx+2, 1,3, &
!                    errmsg='Error: Unable to allocate memory for xnx.')
!            call alloc(yn, -2, imx+2, -2, jmx+3, -2, kmx+2, 1,3, &
!                    errmsg='Error: Unable to allocate memory for ynx.')
!            call alloc(zn, -2, imx+2, -2, jmx+2, -2, kmx+3, 1,3, &
!                    errmsg='Error: Unable to allocate memory for ynx.')
!
!            xnx(-2:imx+3,-2:jmx+2,-2:kmx+2) => xn(:,:,:,1)
!            xny(-2:imx+3,-2:jmx+2,-2:kmx+2) => xn(:,:,:,2)
!            xnz(-2:imx+3,-2:jmx+2,-2:kmx+2) => xn(:,:,:,3)
!
!            ynx(-2:imx+2,-2:jmx+3,-2:kmx+2) => yn(:,:,:,1)
!            yny(-2:imx+2,-2:jmx+3,-2:kmx+2) => yn(:,:,:,2)
!            ynz(-2:imx+2,-2:jmx+3,-2:kmx+2) => yn(:,:,:,3)
!
!            znx(-2:imx+2,-2:jmx+2,-2:kmx+3) => zn(:,:,:,1)
!            zny(-2:imx+2,-2:jmx+2,-2:kmx+3) => zn(:,:,:,2)
!            znz(-2:imx+2,-2:jmx+2,-2:kmx+3) => zn(:,:,:,3)
!
!        end subroutine allocate_memory_normals
!
!        subroutine allocate_memory_ghost_centroids()
!            !< Allocate memory for centroids of ghost cells
!            !-----------------------------------------------------------
!            
!            implicit none
!
!            call alloc(CellCenter, -2, imx+2, -2, jmx+2, -2, kmx+2, 1, 3, &
!                    errmsg='Error: Unable to allocate memory for volume.')
!
!!
!!            call alloc(left_ghost_centroid, 1, jmx-1, 1, kmx-1, 1, 3, &
!!                    errmsg='Error: Unable to allocate memory for left_ghost_centroid')
!!            call alloc(right_ghost_centroid, 1, jmx-1, 1, kmx-1, 1, 3, &
!!                    errmsg='Error: Unable to allocate memory for right_ghost_centroid')
!!            call alloc(front_ghost_centroid, 1, imx-1, 1, kmx-1, 1, 3, &
!!                    errmsg='Error: Unable to allocate memory for front_ghost_centroid')
!!            call alloc(back_ghost_centroid, 1, imx-1, 1, kmx-1, 1, 3, &
!!                    errmsg='Error: Unable to allocate memory for back_ghost_centroid')
!!            call alloc(top_ghost_centroid, 1, imx-1, 1, jmx-1, 1, 3, &
!!                    errmsg='Error: Unable to allocate memory for top_ghost_centroid')
!!            call alloc(bottom_ghost_centroid, 1, imx-1, 1, jmx-1, 1, 3, &
!!                    errmsg='Error: Unable to allocate memory for bottom_ghost_centroid')
!
!        end subroutine allocate_memory_ghost_centroids

        subroutine allocate_memory(cells, Ifaces, Jfaces, Kfaces)
            !< Allocate memory for the required variables.
            !-----------------------------------------------------------
            implicit none
            type(celltype), dimension(:,:,:), allocatable, intent(out) :: cells
            type(facetype), dimension(:,:,:), allocatable, intent(out) :: Ifaces
            type(facetype), dimension(:,:,:), allocatable, intent(out) :: Jfaces
            type(facetype), dimension(:,:,:), allocatable, intent(out) :: Kfaces

            DebugCall('allocate_memory')

            allocate(cells(-2:imx+2, -2:jmx+2, -2:kmx+2))
            !< Allocate memory for cells
            !-----------------------------------------------------------

            allocate(Ifaces(-2:imx+3, -2:jmx+2, -2:kmx+2))
            allocate(Jfaces(-2:imx+2, -2:jmx+3, -2:kmx+2))
            allocate(Kfaces(-2:imx+2, -2:jmx+2, -2:kmx+3))
            !< Allocate memory for the face variables.
            !-----------------------------------------------------------


            call alloc(volume, -2, imx+2, -2, jmx+2, -2, kmx+2, &
                    errmsg='Error: Unable to allocate memory for volume.')

            call alloc(xA, -2, imx+3, -2, jmx+2, -2, kmx+2, &
                    errmsg='Error: Unable to allocate memory for xA.')
            call alloc(yA, -2, imx+2, -2, jmx+3, -2, kmx+2, &
                    errmsg='Error: Unable to allocate memory for yA.')
            call alloc(zA, -2, imx+2, -2, jmx+2, -2, kmx+3, &
                    errmsg='Error: Unable to allocate memory for yA.')

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

            call alloc(CellCenter, -2, imx+2, -2, jmx+2, -2, kmx+2, 1, 3, &
                    errmsg='Error: Unable to allocate memory for volume.')

        end subroutine allocate_memory

!        subroutine deallocate_memory()
!          !< Deallocate the memoery used by the geometry variables
!
!            implicit none
!
!            DebugCall('deallocate_memory')
!
!            call dealloc(xn)
!            call dealloc(yn)
!            call dealloc(zn)
!            call dealloc(xA)
!            call dealloc(yA)
!            call dealloc(zA)
!            call dealloc(volume)
!            call dealloc(CellCenter)
!!            call dealloc(left_ghost_centroid)
!!            call dealloc(right_ghost_centroid)
!!            call dealloc(front_ghost_centroid)
!!            call dealloc(back_ghost_centroid)
!!            call dealloc(top_ghost_centroid)
!!            call dealloc(bottom_ghost_centroid)
!    
!        end subroutine deallocate_memory

        subroutine normalize_face_normals(Ifaces, Jfaces, Kfaces)
            !< Normalize the face normal vectors computed to get unit
            !< vectors
            !-----------------------------------------------------------
            
            implicit none
            type(facetype), dimension(-2:imx+3,-2:jmx+2,-2:kmx+2), intent(inout) :: Ifaces
            type(facetype), dimension(-2:imx+2,-2:jmx+3,-2:kmx+2), intent(inout) :: Jfaces
            type(facetype), dimension(-2:imx+2,-2:jmx+2,-2:kmx+3), intent(inout) :: Kfaces
            integer :: i,j,k

            do k = -2,kmx+2
              do j = -2,jmx+2
                do i = -2,imx+3
                  if(Ifaces(i,j,k)%A/=0.) then
                    Ifaces(i,j,k)%nx = Ifaces(i,j,k)%nx/Ifaces(i,j,k)%A
                    Ifaces(i,j,k)%ny = Ifaces(i,j,k)%ny/Ifaces(i,j,k)%A
                    Ifaces(i,j,k)%nz = Ifaces(i,j,k)%nz/Ifaces(i,j,k)%A
                  end if
                end do
              end do
            end do

            do k = -2,kmx+2
              do j = -2,jmx+3
                do i = -2,imx+2
                  if(Jfaces(i,j,k)%A/=0.) then
                    Jfaces(i,j,k)%nx = Jfaces(i,j,k)%nx/Jfaces(i,j,k)%A
                    Jfaces(i,j,k)%ny = Jfaces(i,j,k)%ny/Jfaces(i,j,k)%A
                    Jfaces(i,j,k)%nz = Jfaces(i,j,k)%nz/Jfaces(i,j,k)%A
                  end if
                end do
              end do
            end do

            do k = -2,kmx+3
              do j = -2,jmx+2
                do i = -2,imx+2
                  if(Kfaces(i,j,k)%A/=0.) then
                    Kfaces(i,j,k)%nx = Kfaces(i,j,k)%nx/Kfaces(i,j,k)%A
                    Kfaces(i,j,k)%ny = Kfaces(i,j,k)%ny/Kfaces(i,j,k)%A
                    Kfaces(i,j,k)%nz = Kfaces(i,j,k)%nz/Kfaces(i,j,k)%A
                  end if
                end do
              end do
            end do

            ! pole boundary condition
            if(imin_id==-7) then
              Ifaces( 1,:,:)%nx=Ifaces(2,:,:)%nx
              Ifaces( 0,:,:)%nx=Ifaces(2,:,:)%nx
              Ifaces(-1,:,:)%nx=Ifaces(2,:,:)%nx
              Ifaces(-2,:,:)%nx=Ifaces(2,:,:)%nx
              Ifaces( 1,:,:)%ny=Ifaces(2,:,:)%ny
              Ifaces( 0,:,:)%ny=Ifaces(2,:,:)%ny
              Ifaces(-1,:,:)%ny=Ifaces(2,:,:)%ny
              Ifaces(-2,:,:)%ny=Ifaces(2,:,:)%ny
              Ifaces( 1,:,:)%nz=Ifaces(2,:,:)%nz
              Ifaces( 0,:,:)%nz=Ifaces(2,:,:)%nz
              Ifaces(-1,:,:)%nz=Ifaces(2,:,:)%nz
              Ifaces(-2,:,:)%nz=Ifaces(2,:,:)%nz
            end if

            if(imax_id==-7) then
              Ifaces(imx+0,:,:)%nx=Ifaces(imx-1,:,:)%nx
              Ifaces(imx+1,:,:)%nx=Ifaces(imx-1,:,:)%nx
              Ifaces(imx+2,:,:)%nx=Ifaces(imx-1,:,:)%nx
              Ifaces(imx+3,:,:)%nx=Ifaces(imx-1,:,:)%nx
              Ifaces(imx+0,:,:)%ny=Ifaces(imx-1,:,:)%ny
              Ifaces(imx+1,:,:)%ny=Ifaces(imx-1,:,:)%ny
              Ifaces(imx+2,:,:)%ny=Ifaces(imx-1,:,:)%ny
              Ifaces(imx+3,:,:)%ny=Ifaces(imx-1,:,:)%ny
              Ifaces(imx+0,:,:)%nz=Ifaces(imx-1,:,:)%nz
              Ifaces(imx+1,:,:)%nz=Ifaces(imx-1,:,:)%nz
              Ifaces(imx+2,:,:)%nz=Ifaces(imx-1,:,:)%nz
              Ifaces(imx+3,:,:)%nz=Ifaces(imx-1,:,:)%nz
            end if

            if(jmin_id==-7) then
              Jfaces(:, 1,:)%nx=Jfaces(:,2,:)%nx
              Jfaces(:, 0,:)%nx=Jfaces(:,2,:)%nx
              Jfaces(:,-1,:)%nx=Jfaces(:,2,:)%nx
              Jfaces(:,-2,:)%nx=Jfaces(:,2,:)%nx
              Jfaces(:, 1,:)%ny=Jfaces(:,2,:)%ny
              Jfaces(:, 0,:)%ny=Jfaces(:,2,:)%ny
              Jfaces(:,-1,:)%ny=Jfaces(:,2,:)%ny
              Jfaces(:,-2,:)%ny=Jfaces(:,2,:)%ny
              Jfaces(:, 1,:)%nz=Jfaces(:,2,:)%nz
              Jfaces(:, 0,:)%nz=Jfaces(:,2,:)%nz
              Jfaces(:,-1,:)%nz=Jfaces(:,2,:)%nz
              Jfaces(:,-2,:)%nz=Jfaces(:,2,:)%nz
            end if

            if(jmax_id==-7) then
              Jfaces(:,jmx+0,:)%nx=Jfaces(:,jmx-1,:)%nx
              Jfaces(:,jmx+1,:)%nx=Jfaces(:,jmx-1,:)%nx
              Jfaces(:,jmx+2,:)%nx=Jfaces(:,jmx-1,:)%nx
              Jfaces(:,jmx+3,:)%nx=Jfaces(:,jmx-1,:)%nx
              Jfaces(:,jmx+0,:)%ny=Jfaces(:,jmx-1,:)%ny
              Jfaces(:,jmx+1,:)%ny=Jfaces(:,jmx-1,:)%ny
              Jfaces(:,jmx+2,:)%ny=Jfaces(:,jmx-1,:)%ny
              Jfaces(:,jmx+3,:)%ny=Jfaces(:,jmx-1,:)%ny
              Jfaces(:,jmx+0,:)%nz=Jfaces(:,jmx-1,:)%nz
              Jfaces(:,jmx+1,:)%nz=Jfaces(:,jmx-1,:)%nz
              Jfaces(:,jmx+2,:)%nz=Jfaces(:,jmx-1,:)%nz
              Jfaces(:,jmx+3,:)%nz=Jfaces(:,jmx-1,:)%nz
            end if

            if(kmin_id==-7) then
              Kfaces(:,:, 1)%nx=Kfaces(:,:,2)%nx
              Kfaces(:,:, 0)%nx=Kfaces(:,:,2)%nx
              Kfaces(:,:,-1)%nx=Kfaces(:,:,2)%nx
              Kfaces(:,:,-2)%nx=Kfaces(:,:,2)%nx
              Kfaces(:,:, 1)%ny=Kfaces(:,:,2)%ny
              Kfaces(:,:, 0)%ny=Kfaces(:,:,2)%ny
              Kfaces(:,:,-1)%ny=Kfaces(:,:,2)%ny
              Kfaces(:,:,-2)%ny=Kfaces(:,:,2)%ny
              Kfaces(:,:, 1)%nz=Kfaces(:,:,2)%nz
              Kfaces(:,:, 0)%nz=Kfaces(:,:,2)%nz
              Kfaces(:,:,-1)%nz=Kfaces(:,:,2)%nz
              Kfaces(:,:,-2)%nz=Kfaces(:,:,2)%nz
            end if

            if(kmax_id==-7) then
              Kfaces(:,:,kmx+0)%nx=Kfaces(:,:,kmx-1)%nx
              Kfaces(:,:,kmx+1)%nx=Kfaces(:,:,kmx-1)%nx
              Kfaces(:,:,kmx+2)%nx=Kfaces(:,:,kmx-1)%nx
              Kfaces(:,:,kmx+3)%nx=Kfaces(:,:,kmx-1)%nx
              Kfaces(:,:,kmx+0)%ny=Kfaces(:,:,kmx-1)%ny
              Kfaces(:,:,kmx+1)%ny=Kfaces(:,:,kmx-1)%ny
              Kfaces(:,:,kmx+2)%ny=Kfaces(:,:,kmx-1)%ny
              Kfaces(:,:,kmx+3)%ny=Kfaces(:,:,kmx-1)%ny
              Kfaces(:,:,kmx+0)%nz=Kfaces(:,:,kmx-1)%nz
              Kfaces(:,:,kmx+1)%nz=Kfaces(:,:,kmx-1)%nz
              Kfaces(:,:,kmx+2)%nz=Kfaces(:,:,kmx-1)%nz
              Kfaces(:,:,kmx+3)%nz=Kfaces(:,:,kmx-1)%nz
            end if

            
        end subroutine normalize_face_normals

        subroutine compute_face_areas(Ifaces, Jfaces, Kfaces)
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
            type(facetype), dimension(-2:imx+3,-2:jmx+2,-2:kmx+2), intent(inout) :: Ifaces
            type(facetype), dimension(-2:imx+2,-2:jmx+3,-2:kmx+2), intent(inout) :: Jfaces
            type(facetype), dimension(-2:imx+2,-2:jmx+2,-2:kmx+3), intent(inout) :: Kfaces

            Ifaces(:,:,:)%A = sqrt((Ifaces(:,:,:)%nx)**2 + (Ifaces(:,:,:)%ny)**2 + &
                                  (Ifaces(:,:,:)%nz)**2)

            Jfaces(:,:,:)%A = sqrt((Jfaces(:,:,:)%nx)**2 + (Jfaces(:,:,:)%ny)**2 + &
                                  (Jfaces(:,:,:)%nz)**2)

            Kfaces(:,:,:)%A = sqrt((Kfaces(:,:,:)%nx)**2 + (Kfaces(:,:,:)%ny)**2 + &
                                  (Kfaces(:,:,:)%nz)**2)

            ! Pole boundary conditions
            ! making sure face area is exactly equal zero
            if(imin_id==-7) Ifaces(-2:1     ,  :,  :)%A=0.
            if(imax_id==-7) Ifaces(imx:imx+3,  :,  :)%A=0.
            if(jmin_id==-7) Jfaces(  :,     -2:1,  :)%A=0.
            if(jmax_id==-7) Jfaces(  :,jmx:jmx+3,  :)%A=0.
            if(kmin_id==-7) Kfaces(  :,  :,     -2:1)%A=0.
            if(kmax_id==-7) Kfaces(  :,  :,kmx:kmx+3)%A=0.

        end subroutine compute_face_areas

        subroutine compute_face_area_vectors(Ifaces, Jfaces, Kfaces, nodes)
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
            type(facetype), dimension(-2:imx+3,-2:jmx+2,-2:kmx+2), intent(inout) :: Ifaces
            type(facetype), dimension(-2:imx+2,-2:jmx+3,-2:kmx+2), intent(inout) :: Jfaces
            type(facetype), dimension(-2:imx+2,-2:jmx+2,-2:kmx+3), intent(inout) :: Kfaces
            type(nodetype), dimension(-2:imx+3,-2:jmx+3,-2:kmx+3), intent(in)  :: nodes

    
            real :: d1x, d2x, d1y, d2y, d1z, d2z
            integer :: i, j, k

            do k = -2, kmx+2
             do j = -2, jmx+2
              do i = -2, imx+3
                d1x = nodes(i, j+1, k+1)%x - nodes(i, j, k)%x
                d1y = nodes(i, j+1, k+1)%y - nodes(i, j, k)%y
                d1z = nodes(i, j+1, k+1)%z - nodes(i, j, k)%z
                d2x = nodes(i, j, k+1)%x - nodes(i, j+1, k)%x
                d2y = nodes(i, j, k+1)%y - nodes(i, j+1, k)%y
                d2z = nodes(i, j, k+1)%z - nodes(i, j+1, k)%z
                Ifaces(i, j, k)%nx = 0.5 * (d1y*d2z - d1z*d2y)
                Ifaces(i, j, k)%ny = 0.5 * (d1z*d2x - d1x*d2z)
                Ifaces(i, j, k)%nz = 0.5 * (d1x*d2y - d1y*d2x)
               end do
              end do
             end do


            do k = -2, kmx+2
             do j = -2, jmx+3
              do i = -2, imx+2
               d1x = nodes(i+1, j, k+1)%x - nodes(i, j, k)%x
               d1y = nodes(i+1, j, k+1)%y - nodes(i, j, k)%y
               d1z = nodes(i+1, j, k+1)%z - nodes(i, j, k)%z
               d2x = nodes(i+1, j, k)%x - nodes(i, j, k+1)%x
               d2y = nodes(i+1, j, k)%y - nodes(i, j, k+1)%y
               d2z = nodes(i+1, j, k)%z - nodes(i, j, k+1)%z
            
               Jfaces(i, j, k)%nx = 0.5 * (d1y*d2z - d1z*d2y)
               Jfaces(i, j, k)%ny = 0.5 * (d1z*d2x - d1x*d2z)
               Jfaces(i, j, k)%nz = 0.5 * (d1x*d2y - d1y*d2x)
              end do
             end do
            end do

            
            do k = -2, kmx+3
             do j = -2, jmx+2
              do i = -2, imx+2
               d1x = nodes(i+1, j+1, k)%x - nodes(i, j, k)%x
               d1y = nodes(i+1, j+1, k)%y - nodes(i, j, k)%y
               d1z = nodes(i+1, j+1, k)%z - nodes(i, j, k)%z
               d2x = nodes(i, j+1, k)%x - nodes(i+1, j, k)%x
               d2y = nodes(i, j+1, k)%y - nodes(i+1, j, k)%y
               d2z = nodes(i, j+1, k)%z - nodes(i+1, j, k)%z
            
               Kfaces(i, j, k)%nx = 0.5 * (d1y*d2z - d1z*d2y)
               Kfaces(i, j, k)%ny = 0.5 * (d1z*d2x - d1x*d2z)
               Kfaces(i, j, k)%nz = 0.5 * (d1x*d2y - d1y*d2x)
              end do
             end do
            end do


        end subroutine compute_face_area_vectors

        subroutine compute_face_areas_and_normals(Ifaces,Jfaces,Kfaces, nodes)
            !< Compute the face areas and normals
            !<
            !< This is the 2-dimensional version. In this case, the face 
            !< areas default to edge lengths.
            !-----------------------------------------------------------

            implicit none
            type(facetype), dimension(-2:imx+3,-2:jmx+2,-2:kmx+2), intent(inout) :: Ifaces
            type(facetype), dimension(-2:imx+2,-2:jmx+3,-2:kmx+2), intent(inout) :: Jfaces
            type(facetype), dimension(-2:imx+2,-2:jmx+2,-2:kmx+3), intent(inout) :: Kfaces
            type(nodetype), dimension(-2:imx+3,-2:jmx+3,-2:kmx+3), intent(in)  :: nodes

            call compute_face_area_vectors(Ifaces,Jfaces,Kfaces, nodes)
            call compute_face_areas(Ifaces,Jfaces,Kfaces)
            call normalize_face_normals(Ifaces,Jfaces,Kfaces)
        
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
        
        subroutine compute_volumes(cells, nodes)
            !< Compute the grid cell volumes
            !< Each grid is a hexahedron, whose volume is calculated by
            !< splitting it into 5 tetrahedrons, whose volume is known
            !-----------------------------------------------------------

            implicit none
            type(celltype), dimension(-2:imx+2,-2:jmx+2,-2:kmx+2), intent(out) :: cells
            type(nodetype), dimension(-2:imx+3,-2:jmx+3,-2:kmx+3), intent(in)  :: nodes
            integer :: i,j,k
            real, dimension(1:3, 1:8) :: p_list

            cells(:,:,:)%volume=1.
            do k = 0, kmx+0
                do j = 0, jmx+0
                    do i = 0, imx+0
                        p_list(:, :) = 0.
                        p_list(:, 1) = (/ nodes(i,j,k)%x, nodes(i,j,k)%y, nodes(i,j,k)%z /)
                        p_list(:, 2) = (/ nodes(i+1,j,k)%x, nodes(i+1,j,k)%y, nodes(i+1,j,k)%z /)
                        p_list(:, 3) = (/ nodes(i+1,j+1,k)%x, nodes(i+1,j+1,k)%y, nodes(i+1,j+1,k)%z /)
                        p_list(:, 4) = (/ nodes(i,j+1,k)%x, nodes(i,j+1,k)%y, nodes(i,j+1,k)%z /)
                        p_list(:, 5) = (/ nodes(i,j,k+1)%x, nodes(i,j,k+1)%y, nodes(i,j,k+1)%z /)
                        p_list(:, 6) = (/ nodes(i+1,j,k+1)%x, nodes(i+1,j,k+1)%y, nodes(i+1,j,k+1)%z /)
                        p_list(:, 7) = (/ nodes(i+1,j+1,k+1)%x, nodes(i+1,j+1,k+1)%y, nodes(i+1,j+1,k+1)%z /)
                        p_list(:, 8) = (/ nodes(i,j+1,k+1)%x, nodes(i,j+1,k+1)%y, nodes(i,j+1,k+1)%z /)
                        cells(i, j, k)%volume = (vol_hexahedron(p_list))
                        if(cells(i,j,k)%volume<=0.0) then
                          if(i==0 .or. i==imx .or. j==0 .or. j==jmx .or. k==0 .or. k==kmx) then
                            !print*, "Ghost Cell volume negative"
                            cells(i, j, k)%volume = (vol_hexahedron(p_list))
                          else
                            print*, process_id, i,j,k
                            print*, "negative volume :", (vol_hexahedron(p_list))
                            STOP
                          end if
                        end if
                    end do
                end do
            end do
            if(any(cells(:,:,:)%volume==0.0))then
              Fatal_error
            end if
            if(any((cells(:,:,:)%volume)<0.0))then
              Fatal_error
            end if
            
        end subroutine compute_volumes

!        subroutine compute_geometric_parameters()
!            !< Compute the geometric parameters based on the grid points
!            !<
!            !< The geometric parameters include the face normals and 
!            !< areas and the cell volumes.
!            !-----------------------------------------------------------
!            
!            implicit none
!
!            DebugCall('compute_geometric_parameters')
!
!            call compute_face_areas_and_normals()
!            call compute_volumes()
!
!        end subroutine compute_geometric_parameters


        subroutine compute_ghost_cell_centroid(cells, nodes)
          !< Compute cell center of all cell including ghost cells
          implicit none
          type(celltype), dimension(-2:imx+2,-2:jmx+2,-2:kmx+2), intent(out) :: cells
          type(nodetype), dimension(-2:imx+3,-2:jmx+3,-2:kmx+3), intent(in)  :: nodes
          integer :: i,j,k

          do k = -2, kmx+2
            do j = -2, jmx+2
              do i = -2, imx+2
                cells(i,j,k)%centerx = 0.125 *( nodes(i  ,j  ,k  )%x &
                                              + nodes(i+1,j  ,k  )%x &
                                              + nodes(i+1,j+1,k  )%x &
                                              + nodes(i+1,j+1,k+1)%x &
                                              + nodes(i+1,j  ,k+1)%x &
                                              + nodes(i  ,j+1,k  )%x &
                                              + nodes(i  ,j+1,k+1)%x &
                                              + nodes(i  ,j  ,k+1)%x &
                                              )

                cells(i,j,k)%centery = 0.125 *( nodes(i  ,j  ,k  )%y &
                                              + nodes(i+1,j  ,k  )%y &
                                              + nodes(i+1,j+1,k  )%y &
                                              + nodes(i+1,j+1,k+1)%y &
                                              + nodes(i+1,j  ,k+1)%y &
                                              + nodes(i  ,j+1,k  )%y &
                                              + nodes(i  ,j+1,k+1)%y &
                                              + nodes(i  ,j  ,k+1)%y &
                                              )

                cells(i,j,k)%centerz = 0.125 *( nodes(i  ,j  ,k  )%z  &
                                              + nodes(i+1,j  ,k  )%z  &
                                              + nodes(i+1,j+1,k  )%z  &
                                              + nodes(i+1,j+1,k+1)%z  &
                                              + nodes(i+1,j  ,k+1)%z  &
                                              + nodes(i  ,j+1,k  )%z  &
                                              + nodes(i  ,j+1,k+1)%z  &
                                              + nodes(i  ,j  ,k+1)%z  &
                                              )
              end do
            end do
          end do

        end subroutine compute_ghost_cell_centroid


        subroutine setup_geometry(cells, Ifaces, Jfaces, Kfaces, nodes, dims)
            !< Make the geometry module useful
            !<
            !< Allocates memory to the variables and initializes them.
            !-----------------------------------------------------------

            implicit none
            type(extent), intent(in) :: dims
            type(celltype), dimension(:,:,:), allocatable, intent(inout) :: cells
            type(facetype), dimension(:,:,:), allocatable, intent(inout) :: Ifaces
            type(facetype), dimension(:,:,:), allocatable, intent(inout) :: Jfaces
            type(facetype), dimension(:,:,:), allocatable, intent(inout) :: Kfaces
            type(nodetype), dimension(-2:imx+3,-2:jmx+3,-2:kmx+3), intent(in)  :: nodes

            DebugCall('setup_geometry')

            imx = dims%imx
            jmx = dims%jmx
            kmx = dims%kmx

            call allocate_memory(cells, Ifaces, Jfaces, Kfaces)
            call compute_face_areas_and_normals(Ifaces,Jfaces,Kfaces, nodes)
            call compute_volumes(cells, nodes)
            call compute_ghost_cell_centroid(cells, nodes)

            volume = cells(:,:,:)%volume
            cellcenter(:,:,:,1) = cells(-2:imx+2,-2:jmx+2,-2:kmx+2)%centerx
            cellcenter(:,:,:,2) = cells(-2:imx+2,-2:jmx+2,-2:kmx+2)%centery
            cellcenter(:,:,:,3) = cells(-2:imx+2,-2:jmx+2,-2:kmx+2)%centerz
            xA = Ifaces(:,:,:)%A
            yA = Jfaces(:,:,:)%A
            zA = Kfaces(:,:,:)%A
            xn(:,:,:,1) = Ifaces(:,:,:)%nx
            xn(:,:,:,2) = Ifaces(:,:,:)%ny
            xn(:,:,:,3) = Ifaces(:,:,:)%nz
            yn(:,:,:,1) = Jfaces(:,:,:)%nx
            yn(:,:,:,2) = Jfaces(:,:,:)%ny
            yn(:,:,:,3) = Jfaces(:,:,:)%nz
            zn(:,:,:,1) = Kfaces(:,:,:)%nx
            zn(:,:,:,2) = Kfaces(:,:,:)%ny
            zn(:,:,:,3) = Kfaces(:,:,:)%nz

        end subroutine setup_geometry

!        subroutine destroy_geometry()
!          !< Nullify all the face normal pionter 
!
!            implicit none
!            
!            DebugCall('destroy_geometry')
!
!            nullify(xnx)
!            nullify(xny)
!            nullify(xnz)
!            nullify(ynx)
!            nullify(yny)
!            nullify(ynz)
!            nullify(znx)
!            nullify(zny)
!            nullify(znz)
!            call deallocate_memory()
!
!        end subroutine destroy_geometry

end module geometry
