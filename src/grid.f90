    !< The grid module read grid file and allocate memory to storing variables
module grid
    !< The grid module contains the grid definition (locations of the 
    !< grid points) as well as procedures to load these from a file.
    !-------------------------------------------------------------------
    
    use global, only: STRING_BUFFER_LENGTH, GRID_FILE_UNIT
    use global_vars, only : imx
    use global_vars, only : jmx
    use global_vars, only : kmx
    use global_vars, only : grid_x
    use global_vars, only : grid_y
    use global_vars, only : grid_z
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
    private
    
    ! Public methods
    public :: setup_grid
    public :: destroy_grid

    contains
    
        subroutine allocate_memory()
            !< Allocate memory to store the grid
            !-----------------------------------------------------------

            implicit none

            call dmsg(1, 'grid', 'allocate_memory')

            call alloc(grid_x, -2, imx+3, -2, jmx+3, -2, kmx+3, &
                    errmsg='Error: Unable to allocate memory for grid_x.')
            call alloc(grid_y, -2, imx+3, -2, jmx+3, -2, kmx+3, &
                    errmsg='Error: Unable to allocate memory for grid_y.')
            call alloc(grid_z, -2, imx+3, -2, jmx+3, -2, kmx+3, &
                    errmsg='Error: Unable to allocate memory for grid_z.')
            
            ! The alloc function earlier allocates only real arrays.
            ! A new function has been written to allocate for integers
            ! as well

!           call alloc(sphere_indices, 1, 3, 1, n_sph_ind, &
!                   errmsg='Error: Unable to allocate memory for sphere_indices.')
            
        end subroutine allocate_memory

        subroutine destroy_grid()
            !< Deallocate the memory allocated for the grid.
            !-----------------------------------------------------------

            implicit none

            call dmsg(1, 'grid', 'destroy_memory')

            call dealloc(grid_x)
            call dealloc(grid_y)
            call dealloc(grid_z)

!           call dealloc(sphere_indices)

        end subroutine destroy_grid

        subroutine setup_grid(gridfile)
            !< Read the grid file and initialize the grid
            !-----------------------------------------------------------

            implicit none
            character(len=64), intent(in) :: gridfile
   !        character(len=32) :: sphindfile
            
            call dmsg(1, 'grid', 'setup_grid')
!           sphindfile = 'sphere-indices.txt'

            open(GRID_FILE_UNIT, file=gridfile)

!           open(SPHERE_INDICES_FILE_UNIT, file=sphindfile)

            call extract_grid_size()
!           call extract_sphere_indices_size()

            call allocate_memory()

            !read interface mapping
            call read_interface_map()

            ! ghost grid exchange
            call populate_grid_points()
!           call populate_sphere_indices()

            close(GRID_FILE_UNIT)
!           close(SPHERE_INDICES_FILE_UNIT)

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
                read(line, *) grid_x(i, j, k), grid_y(i, j, k), grid_z(i, j, k)
            else    
                if (jmx > 1) then
                    read(line, *) grid_x(i, j, k), grid_y(i, j, k)
                else
                    read(line, *) grid_x(i, j, k)
                    grid_y(i, j, k) = 0.
                end if
                grid_z(i, j, k) = 0.
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
          implicit none
          integer :: count=0
          integer :: i,j,k,l
          integer :: ierr
          integer :: status(MPI_STATUS_SIZE)
          real, dimension(:),allocatable :: imin_send_buffer
          real, dimension(:),allocatable :: jmin_send_buffer
          real, dimension(:),allocatable :: kmin_send_buffer
          real, dimension(:),allocatable :: imin_recv_buffer
          real, dimension(:),allocatable :: jmin_recv_buffer
          real, dimension(:),allocatable :: kmin_recv_buffer
          real, dimension(:),allocatable :: imax_send_buffer
          real, dimension(:),allocatable :: jmax_send_buffer
          real, dimension(:),allocatable :: kmax_send_buffer
          real, dimension(:),allocatable :: imax_recv_buffer
          real, dimension(:),allocatable :: jmax_recv_buffer
          real, dimension(:),allocatable :: kmax_recv_buffer
!          real, dimension(3) :: delI
!          real, dimension(3) :: N
!          real  :: dot
!          real  :: magnitude
!          real :: d1x, d2x, d1y, d2y, d1z, d2z
!          integer :: ii,jj,kk
!
!          call dmsg(1, 'grid', 'ghost_grid_extrapolate')
!          !---IMIN---!
!          do k = 1, kmx
!           do j = 1, jmx
!             ! finding face normal
!              ii=1
!              jj=j
!              kk=k
!              if(j==jmx) jj=j-1
!              if(k==kmx) kk=k-1
!              d1x = grid_x(ii, jj+1, kk+1) - grid_x(ii, jj  , kk)
!              d1y = grid_y(ii, jj+1, kk+1) - grid_y(ii, jj  , kk)
!              d1z = grid_z(ii, jj+1, kk+1) - grid_z(ii, jj  , kk)
!              d2x = grid_x(ii, jj  , kk+1) - grid_x(ii, jj+1, kk)
!              d2y = grid_y(ii, jj  , kk+1) - grid_y(ii, jj+1, kk)
!              d2z = grid_z(ii, jj  , kk+1) - grid_z(ii, jj+1, kk)
!              N(1) = -0.5 * (d1y*d2z - d1z*d2y)
!              N(2) = -0.5 * (d1z*d2x - d1x*d2z)
!              N(3) = -0.5 * (d1x*d2y - d1y*d2x)
!              magnitude = SUM(N**2)
!              if(magnitude==0)then
!                d1x = grid_x(ii+1, jj+1, kk+1) - grid_x(ii+1, jj  , kk)
!                d1y = grid_y(ii+1, jj+1, kk+1) - grid_y(ii+1, jj  , kk)
!                d1z = grid_z(ii+1, jj+1, kk+1) - grid_z(ii+1, jj  , kk)
!                d2x = grid_x(ii+1, jj  , kk+1) - grid_x(ii+1, jj+1, kk)
!                d2y = grid_y(ii+1, jj  , kk+1) - grid_y(ii+1, jj+1, kk)
!                d2z = grid_z(ii+1, jj  , kk+1) - grid_z(ii+1, jj+1, kk)
!                N(1) = -0.5 * (d1y*d2z - d1z*d2y)
!                N(2) = -0.5 * (d1z*d2x - d1x*d2z)
!                N(3) = -0.5 * (d1x*d2y - d1y*d2x)
!                magnitude = SUM(N**2)
!                if(magnitude==0)then
!                  Fatal_error
!                end if
!              end if
!              do l=1,layers
!                delI(1) = grid_x(ii+1,j,k)-grid_x(ii,j,k)
!                delI(2) = grid_y(ii+1,j,k)-grid_y(ii,j,k)
!                delI(3) = grid_z(ii+1,j,k)-grid_z(ii,j,k)
!                dot=SUM(N*delI)/magnitude
!                grid_x(ii-l,j,k)= grid_x(ii,j,k) + delI(1) + (l+1)*dot*N(1) 
!                grid_y(ii-l,j,k)= grid_y(ii,j,k) + delI(2) + (l+1)*dot*N(2) 
!                grid_z(ii-l,j,k)= grid_z(ii,j,k) + delI(3) + (l+1)*dot*N(3) 
!             end do
!            end do
!           end do
!
!
!           !--- IMAX ---!
!          do k = 1, kmx
!           do j = 1, jmx
!             ! finding face normal
!              ii=imx
!              jj=j
!              kk=k
!              if(j==jmx) jj=j-1
!              if(k==kmx) kk=k-1
!              d1x = grid_x(ii, jj+1, kk+1) - grid_x(ii, jj  , kk)
!              d1y = grid_y(ii, jj+1, kk+1) - grid_y(ii, jj  , kk)
!              d1z = grid_z(ii, jj+1, kk+1) - grid_z(ii, jj  , kk)
!              d2x = grid_x(ii, jj  , kk+1) - grid_x(ii, jj+1, kk)
!              d2y = grid_y(ii, jj  , kk+1) - grid_y(ii, jj+1, kk)
!              d2z = grid_z(ii, jj  , kk+1) - grid_z(ii, jj+1, kk)
!              N(1) = 0.5 * (d1y*d2z - d1z*d2y)
!              N(2) = 0.5 * (d1z*d2x - d1x*d2z)
!              N(3) = 0.5 * (d1x*d2y - d1y*d2x)
!              magnitude = SUM(N**2)
!              if(magnitude==0)then
!                d1x = grid_x(ii-1, jj+1, kk+1) - grid_x(ii-1, jj  , kk)
!                d1y = grid_y(ii-1, jj+1, kk+1) - grid_y(ii-1, jj  , kk)
!                d1z = grid_z(ii-1, jj+1, kk+1) - grid_z(ii-1, jj  , kk)
!                d2x = grid_x(ii-1, jj  , kk+1) - grid_x(ii-1, jj+1, kk)
!                d2y = grid_y(ii-1, jj  , kk+1) - grid_y(ii-1, jj+1, kk)
!                d2z = grid_z(ii-1, jj  , kk+1) - grid_z(ii-1, jj+1, kk)
!                N(1) = 0.5 * (d1y*d2z - d1z*d2y)
!                N(2) = 0.5 * (d1z*d2x - d1x*d2z)
!                N(3) = 0.5 * (d1x*d2y - d1y*d2x)
!                magnitude = SUM(N**2)
!                if(magnitude==0)then
!                  Fatal_error
!                end if
!              end if
!              do l=1,layers
!                delI(1) = grid_x(ii-1,j,k)-grid_x(ii,j,k)
!                delI(2) = grid_y(ii-1,j,k)-grid_y(ii,j,k)
!                delI(3) = grid_z(ii-1,j,k)-grid_z(ii,j,k)
!                dot=SUM(N*delI)/magnitude
!                grid_x(ii+l,j,k)= grid_x(ii,j,k) + delI(1) + (l+1)*dot*N(1) 
!                grid_y(ii+l,j,k)= grid_y(ii,j,k) + delI(2) + (l+1)*dot*N(2) 
!                grid_z(ii+l,j,k)= grid_z(ii,j,k) + delI(3) + (l+1)*dot*N(3) 
!            end do
!           end do
!          end do
!
!          !--- JMIN ---!
!          do k = 1, kmx
!           do i = -2, imx+3
!             ! finding face normal
!              ii=i
!              jj=1
!              kk=k
!              if(i==imx+3) ii=i-1
!              if(k==kmx) kk=k-1
!              d1x = grid_x(ii+1, jj, kk+1) - grid_x(ii, jj, kk)
!              d1y = grid_y(ii+1, jj, kk+1) - grid_y(ii, jj, kk)
!              d1z = grid_z(ii+1, jj, kk+1) - grid_z(ii, jj, kk)
!              d2x = grid_x(ii+1, jj, kk)   - grid_x(ii, jj, kk+1)
!              d2y = grid_y(ii+1, jj, kk)   - grid_y(ii, jj, kk+1)
!              d2z = grid_z(ii+1, jj, kk)   - grid_z(ii, jj, kk+1)
!              N(1) = 0.5 * (d1y*d2z - d1z*d2y)
!              N(2) = 0.5 * (d1z*d2x - d1x*d2z)
!              N(3) = 0.5 * (d1x*d2y - d1y*d2x)
!              magnitude = SUM(N**2)
!              if(magnitude==0)then
!                d1x = grid_x(ii+1, jj+1, kk+1) - grid_x(ii, jj+1, kk)
!                d1y = grid_y(ii+1, jj+1, kk+1) - grid_y(ii, jj+1, kk)
!                d1z = grid_z(ii+1, jj+1, kk+1) - grid_z(ii, jj+1, kk)
!                d2x = grid_x(ii+1, jj+1, kk)   - grid_x(ii, jj+1, kk+1)
!                d2y = grid_y(ii+1, jj+1, kk)   - grid_y(ii, jj+1, kk+1)
!                d2z = grid_z(ii+1, jj+1, kk)   - grid_z(ii, jj+1, kk+1)
!                N(1) = 0.5 * (d1y*d2z - d1z*d2y)
!                N(2) = 0.5 * (d1z*d2x - d1x*d2z)
!                N(3) = 0.5 * (d1x*d2y - d1y*d2x)
!                magnitude = SUM(N**2)
!                if(magnitude==0)then
!                  Fatal_error
!                end if
!              end if
!              do l=1,layers
!                delI(1) = grid_x(i,jj+1,k)-grid_x(i,jj,k)
!                delI(2) = grid_y(i,jj+1,k)-grid_y(i,jj,k)
!                delI(3) = grid_z(i,jj+1,k)-grid_z(i,jj,k)
!                dot=SUM(N*delI)/magnitude
!                grid_x(i,jj-l,k)= grid_x(i,jj,k) + delI(1) - (l+1)*dot*N(1) 
!                grid_y(i,jj-l,k)= grid_y(i,jj,k) + delI(2) - (l+1)*dot*N(2) 
!                grid_z(i,jj-l,k)= grid_z(i,jj,k) + delI(3) - (l+1)*dot*N(3) 
!              end do
!            end do
!          end do
!
!          
!          !--- JMAX ---!
!          do k = 1, kmx
!           do i = -2, imx+3
!             ! finding face normal
!              ii=i
!              jj=jmx
!              kk=k
!              if(i==imx+3) ii=i-1
!              if(k==kmx) kk=k-1
!              d1x = grid_x(ii+1, jj, kk+1) - grid_x(ii, jj, kk)
!              d1y = grid_y(ii+1, jj, kk+1) - grid_y(ii, jj, kk)
!              d1z = grid_z(ii+1, jj, kk+1) - grid_z(ii, jj, kk)
!              d2x = grid_x(ii+1, jj, kk  ) - grid_x(ii, jj, kk+1)
!              d2y = grid_y(ii+1, jj, kk  ) - grid_y(ii, jj, kk+1)
!              d2z = grid_z(ii+1, jj, kk  ) - grid_z(ii, jj, kk+1)
!              N(1) = 0.5 * (d1y*d2z - d1z*d2y)
!              N(2) = 0.5 * (d1z*d2x - d1x*d2z)
!              N(3) = 0.5 * (d1x*d2y - d1y*d2x)
!              magnitude = SUM(N**2)
!              if(magnitude==0)then
!                d1x = grid_x(ii+1, jj-1, kk+1) - grid_x(ii, jj-1, kk)
!                d1y = grid_y(ii+1, jj-1, kk+1) - grid_y(ii, jj-1, kk)
!                d1z = grid_z(ii+1, jj-1, kk+1) - grid_z(ii, jj-1, kk)
!                d2x = grid_x(ii+1, jj-1, kk  ) - grid_x(ii, jj-1, kk+1)
!                d2y = grid_y(ii+1, jj-1, kk  ) - grid_y(ii, jj-1, kk+1)
!                d2z = grid_z(ii+1, jj-1, kk  ) - grid_z(ii, jj-1, kk+1)
!                N(1) = 0.5 * (d1y*d2z - d1z*d2y)
!                N(2) = 0.5 * (d1z*d2x - d1x*d2z)
!                N(3) = 0.5 * (d1x*d2y - d1y*d2x)
!                magnitude = SUM(N**2)
!                if(magnitude==0)then
!                  Fatal_error
!                end if
!              end if
!              do l=1,layers
!                delI(1) = grid_x(i,jj-1,k)-grid_x(i,jj,k)
!                delI(2) = grid_y(i,jj-1,k)-grid_y(i,jj,k)
!                delI(3) = grid_z(i,jj-1,k)-grid_z(i,jj,k)
!                dot=SUM(N*delI)/magnitude
!                grid_x(i,jj+l,k)= grid_x(i,jj,k) + delI(1) - (l+1)*dot*N(1) 
!                grid_y(i,jj+l,k)= grid_y(i,jj,k) + delI(2) - (l+1)*dot*N(2) 
!                grid_z(i,jj+l,k)= grid_z(i,jj,k) + delI(3) - (l+1)*dot*N(3) 
!              end do
!            end do
!          end do
!
!          !--- KMIN ---!
!          do j = -2, jmx+3
!           do i = -2, imx+3
!              ii=i
!              jj=j
!              kk=1
!              if(i==imx+3) ii=i-1
!              if(j==jmx+3) jj=j-1
!              d1x = grid_x(ii+1, jj+1, kk) - grid_x(ii  , jj, kk)
!              d1y = grid_y(ii+1, jj+1, kk) - grid_y(ii  , jj, kk)
!              d1z = grid_z(ii+1, jj+1, kk) - grid_z(ii  , jj, kk)
!              d2x = grid_x(ii  , jj+1, kk) - grid_x(ii+1, jj, kk)
!              d2y = grid_y(ii  , jj+1, kk) - grid_y(ii+1, jj, kk)
!              d2z = grid_z(ii  , jj+1, kk) - grid_z(ii+1, jj, kk)
!              N(1) = 0.5 * (d1y*d2z - d1z*d2y)
!              N(2) = 0.5 * (d1z*d2x - d1x*d2z)
!              N(3) = 0.5 * (d1x*d2y - d1y*d2x)
!              magnitude = SUM(N**2)
!              if(magnitude==0)then
!                d1x = grid_x(ii+1, jj+1, kk+1) - grid_x(ii  , jj, kk+1)
!                d1y = grid_y(ii+1, jj+1, kk+1) - grid_y(ii  , jj, kk+1)
!                d1z = grid_z(ii+1, jj+1, kk+1) - grid_z(ii  , jj, kk+1)
!                d2x = grid_x(ii  , jj+1, kk+1) - grid_x(ii+1, jj, kk+1)
!                d2y = grid_y(ii  , jj+1, kk+1) - grid_y(ii+1, jj, kk+1)
!                d2z = grid_z(ii  , jj+1, kk+1) - grid_z(ii+1, jj, kk+1)
!                N(1) = 0.5 * (d1y*d2z - d1z*d2y)
!                N(2) = 0.5 * (d1z*d2x - d1x*d2z)
!                N(3) = 0.5 * (d1x*d2y - d1y*d2x)
!                magnitude = SUM(N**2)
!                if(magnitude==0)then
!                  Fatal_error
!                end if
!              end if
!              do l=1,layers
!                delI(1) = grid_x(i,j,kk+1)-grid_x(i,j,kk)
!                delI(2) = grid_y(i,j,kk+1)-grid_y(i,j,kk)
!                delI(3) = grid_z(i,j,kk+1)-grid_z(i,j,kk)
!                dot=SUM(N*delI)/magnitude
!                grid_x(i,j,kk-l)= grid_x(i,j,kk) + delI(1) - (l+1)*dot*N(1) 
!                grid_y(i,j,kk-l)= grid_y(i,j,kk) + delI(2) - (l+1)*dot*N(2) 
!                grid_z(i,j,kk-l)= grid_z(i,j,kk) + delI(3) - (l+1)*dot*N(3) 
!            end do
!           end do
!          end do
!
!          !--- KMAX ---!
!          do j = -2, jmx+3
!           do i = -2, imx+3
!              ii=i
!              jj=j
!              kk=kmx
!              if(i==imx+3) ii=i-1
!              if(j==jmx+3) jj=j-1
!              d1x = grid_x(ii+1, jj+1, kk) - grid_x(ii  , jj, kk)
!              d1y = grid_y(ii+1, jj+1, kk) - grid_y(ii  , jj, kk)
!              d1z = grid_z(ii+1, jj+1, kk) - grid_z(ii  , jj, kk)
!              d2x = grid_x(ii  , jj+1, kk) - grid_x(ii+1, jj, kk)
!              d2y = grid_y(ii  , jj+1, kk) - grid_y(ii+1, jj, kk)
!              d2z = grid_z(ii  , jj+1, kk) - grid_z(ii+1, jj, kk)
!              N(1) = 0.5 * (d1y*d2z - d1z*d2y)
!              N(2) = 0.5 * (d1z*d2x - d1x*d2z)
!              N(3) = 0.5 * (d1x*d2y - d1y*d2x)
!              magnitude = SUM(N**2)
!              if(magnitude==0)then
!                d1x = grid_x(ii+1, jj+1, kk-1) - grid_x(ii  , jj, kk-1)
!                d1y = grid_y(ii+1, jj+1, kk-1) - grid_y(ii  , jj, kk-1)
!                d1z = grid_z(ii+1, jj+1, kk-1) - grid_z(ii  , jj, kk-1)
!                d2x = grid_x(ii  , jj+1, kk-1) - grid_x(ii+1, jj, kk-1)
!                d2y = grid_y(ii  , jj+1, kk-1) - grid_y(ii+1, jj, kk-1)
!                d2z = grid_z(ii  , jj+1, kk-1) - grid_z(ii+1, jj, kk-1)
!                N(1) = 0.5 * (d1y*d2z - d1z*d2y)
!                N(2) = 0.5 * (d1z*d2x - d1x*d2z)
!                N(3) = 0.5 * (d1x*d2y - d1y*d2x)
!                magnitude = SUM(N**2)
!                if(magnitude==0)then
!                  Fatal_error
!                end if
!              end if
!              do l=1,layers
!                delI(1) = grid_x(i,j,kk-1)-grid_x(i,j,kk)
!                delI(2) = grid_y(i,j,kk-1)-grid_y(i,j,kk)
!                delI(3) = grid_z(i,j,kk-1)-grid_z(i,j,kk)
!                dot=SUM(N*delI)/magnitude
!                grid_x(i,j,kk+l)= grid_x(i,j,kk) + delI(1) - (l+1)*dot*N(1) 
!                grid_y(i,j,kk+l)= grid_y(i,j,kk) + delI(2) - (l+1)*dot*N(2) 
!                grid_z(i,j,kk+l)= grid_z(i,j,kk) + delI(3) - (l+1)*dot*N(3) 
!                print*, grid_x(ii+l,j,k), grid_y(ii+l,j,k), grid_z(ii+l,j,k)
!            end do
!           end do
!          end do


!          !-------------------------------------------------------------------
!          !getting ghost cell for all faces even if it is a interface cell
!          ! <algorithm>
!          ! Point_ghost = 2*Point_first_inner_cell - Point_second_inner_cell
!          ! </algorithm>
!          !-------------------------------------------------------------------
!
          !--- I faces ---!
          !imin face -> 0 grid point
          grid_x( 0,:,:)=2*grid_x( 1,:,:)-grid_x(2,:,:)
          grid_y( 0,:,:)=2*grid_y( 1,:,:)-grid_y(2,:,:)
          grid_z( 0,:,:)=2*grid_z( 1,:,:)-grid_z(2,:,:)
          !imin face -> -1 grid point
          grid_x(-1,:,:)=2*grid_x( 0,:,:)-grid_x(1,:,:)
          grid_y(-1,:,:)=2*grid_y( 0,:,:)-grid_y(1,:,:)
          grid_z(-1,:,:)=2*grid_z( 0,:,:)-grid_z(1,:,:)
          !imin face -> -2 grid point
          grid_x(-2,:,:)=2*grid_x(-1,:,:)-grid_x(0,:,:)
          grid_y(-2,:,:)=2*grid_y(-1,:,:)-grid_y(0,:,:)
          grid_z(-2,:,:)=2*grid_z(-1,:,:)-grid_z(0,:,:)

          !imax face -> imx+1 grid point
          grid_x(imx+1,:,:)=2*grid_x(imx+0,:,:)-grid_x(imx-1,:,:)
          grid_y(imx+1,:,:)=2*grid_y(imx+0,:,:)-grid_y(imx-1,:,:)
          grid_z(imx+1,:,:)=2*grid_z(imx+0,:,:)-grid_z(imx-1,:,:)
          !imax face -> imx+2 grid point
          grid_x(imx+2,:,:)=2*grid_x(imx+1,:,:)-grid_x(imx-0,:,:)
          grid_y(imx+2,:,:)=2*grid_y(imx+1,:,:)-grid_y(imx-0,:,:)
          grid_z(imx+2,:,:)=2*grid_z(imx+1,:,:)-grid_z(imx-0,:,:)
          !imax face -> imx+3 grid point
          grid_x(imx+3,:,:)=2*grid_x(imx+2,:,:)-grid_x(imx+1,:,:)
          grid_y(imx+3,:,:)=2*grid_y(imx+2,:,:)-grid_y(imx+1,:,:)
          grid_z(imx+3,:,:)=2*grid_z(imx+2,:,:)-grid_z(imx+1,:,:)


          !--- Jmin faces ---!
          !jmin faces -> 0 grid point
          grid_x(:, 0,:)=2*grid_x(:, 1,:)-grid_x(:,2,:)
          grid_y(:, 0,:)=2*grid_y(:, 1,:)-grid_y(:,2,:)
          grid_z(:, 0,:)=2*grid_z(:, 1,:)-grid_z(:,2,:)
          !jmin face -> -1 grid point
          grid_x(:,-1,:)=2*grid_x(:, 0,:)-grid_x(:,1,:)
          grid_y(:,-1,:)=2*grid_y(:, 0,:)-grid_y(:,1,:)
          grid_z(:,-1,:)=2*grid_z(:, 0,:)-grid_z(:,1,:)
          !jmin face -> -2 grid point
          grid_x(:,-2,:)=2*grid_x(:,-1,:)-grid_x(:,0,:)
          grid_y(:,-2,:)=2*grid_y(:,-1,:)-grid_y(:,0,:)
          grid_z(:,-2,:)=2*grid_z(:,-1,:)-grid_z(:,0,:)

          !jmax face -> imx+1 grid point
          grid_x(:,jmx+1,:)=2*grid_x(:,jmx+0,:)-grid_x(:,jmx-1,:)
          grid_y(:,jmx+1,:)=2*grid_y(:,jmx+0,:)-grid_y(:,jmx-1,:)
          grid_z(:,jmx+1,:)=2*grid_z(:,jmx+0,:)-grid_z(:,jmx-1,:)
          !jmax face -> imx+2 grid point
          grid_x(:,jmx+2,:)=2*grid_x(:,jmx+1,:)-grid_x(:,jmx-0,:)
          grid_y(:,jmx+2,:)=2*grid_y(:,jmx+1,:)-grid_y(:,jmx-0,:)
          grid_z(:,jmx+2,:)=2*grid_z(:,jmx+1,:)-grid_z(:,jmx-0,:)
          !jmax face -> imx+3 grid point
          grid_x(:,jmx+3,:)=2*grid_x(:,jmx+2,:)-grid_x(:,jmx+1,:)
          grid_y(:,jmx+3,:)=2*grid_y(:,jmx+2,:)-grid_y(:,jmx+1,:)
          grid_z(:,jmx+3,:)=2*grid_z(:,jmx+2,:)-grid_z(:,jmx+1,:)


          !--- Kmax faces ---!
          !kmin faces -> 0 grid point
          grid_x(:,:, 0)=2*grid_x(:,:, 1)-grid_x(:,:,2)
          grid_y(:,:, 0)=2*grid_y(:,:, 1)-grid_y(:,:,2)
          grid_z(:,:, 0)=2*grid_z(:,:, 1)-grid_z(:,:,2)
          !kmin face -> -1 grid point
          grid_x(:,:,-1)=2*grid_x(:,:, 0)-grid_x(:,:,1)
          grid_y(:,:,-1)=2*grid_y(:,:, 0)-grid_y(:,:,1)
          grid_z(:,:,-1)=2*grid_z(:,:, 0)-grid_z(:,:,1)
          !kmin face -> -2 grid point
          grid_x(:,:,-2)=2*grid_x(:,:,-1)-grid_x(:,:,0)
          grid_y(:,:,-2)=2*grid_y(:,:,-1)-grid_y(:,:,0)
          grid_z(:,:,-2)=2*grid_z(:,:,-1)-grid_z(:,:,0)

          !kmax face -> imx+1 grid point
          grid_x(:,:,kmx+1)=2*grid_x(:,:,kmx+0)-grid_x(:,:,kmx-1)
          grid_y(:,:,kmx+1)=2*grid_y(:,:,kmx+0)-grid_y(:,:,kmx-1)
          grid_z(:,:,kmx+1)=2*grid_z(:,:,kmx+0)-grid_z(:,:,kmx-1)
          !kmax face -> imx+2 grid point
          grid_x(:,:,kmx+2)=2*grid_x(:,:,kmx+1)-grid_x(:,:,kmx-0)
          grid_y(:,:,kmx+2)=2*grid_y(:,:,kmx+1)-grid_y(:,:,kmx-0)
          grid_z(:,:,kmx+2)=2*grid_z(:,:,kmx+1)-grid_z(:,:,kmx-0)
          !kmax face -> imx+3 grid point
          grid_x(:,:,kmx+3)=2*grid_x(:,:,kmx+2)-grid_x(:,:,kmx+1)
          grid_y(:,:,kmx+3)=2*grid_y(:,:,kmx+2)-grid_y(:,:,kmx+1)
          grid_z(:,:,kmx+3)=2*grid_z(:,:,kmx+2)-grid_z(:,:,kmx+1)

          !print*, grid_x(:,:,kmx)
          !print*, grid_y(:,:,kmx)
          !print*, grid_z(:,:,kmx)
          call dmsg(1, 'grid', 'ghost_grid_interface')
          !---  MPI transfer of grid point across interface  ---!
          !--- imin face ---!
          allocate(imin_send_buffer(3*layers*(jmx+6)*(kmx+6)))
          allocate(jmin_send_buffer(3*layers*(imx+6)*(kmx+6)))
          allocate(kmin_send_buffer(3*layers*(imx+6)*(jmx+6)))
          allocate(imin_recv_buffer(3*layers*(jmx+6)*(kmx+6)))
          allocate(jmin_recv_buffer(3*layers*(imx+6)*(kmx+6)))
          allocate(kmin_recv_buffer(3*layers*(imx+6)*(jmx+6)))
          allocate(imax_send_buffer(3*layers*(jmx+6)*(kmx+6)))
          allocate(jmax_send_buffer(3*layers*(imx+6)*(kmx+6)))
          allocate(kmax_send_buffer(3*layers*(imx+6)*(jmx+6)))
          allocate(imax_recv_buffer(3*layers*(jmx+6)*(kmx+6)))
          allocate(jmax_recv_buffer(3*layers*(imx+6)*(kmx+6)))
          allocate(kmax_recv_buffer(3*layers*(imx+6)*(jmx+6)))
!          if(imin_id>=0)then
!            !collect grid point in 1d array
!            count=0
!            do l=1,layers
!              do k=-2,kmx+3
!                do j=-2,jmx+3
!                  count=count+1
!                  imin_send_buffer(count) = grid_x(l+1,j,k)
!                end do
!              end do
!            end do
!
!            do l=1,layers
!              do k=-2,kmx+3
!                do j=-2,jmx+3
!                  count=count+1
!                  imin_send_buffer(count) = grid_y(l+1,j,k)
!                end do
!              end do
!            end do
!
!            do l=1,layers
!              do k=-2,kmx+3
!                do j=-2,jmx+3
!                  count=count+1
!                  imin_send_buffer(count) = grid_z(l+1,j,k)
!                end do
!              end do
!            end do
!
!        call MPI_SENDRECV(imin_send_buffer,count, MPI_DOUBLE_PRECISION, imin_id,1,&
!                          imin_recv_buffer,count, MPI_DOUBLE_PRECISION, imin_id,1,&
!                          MPI_COMM_WORLD,status,ierr)
!        !    if(mpi_class(1)==0)then
!        !      print*, Process_id, "imin master"
!        !      call MPI_SEND(imin_send_buffer, count,MPI_DOUBLE_PRECISION,imin_id,1,MPI_COMM_WORLD, ierr)
!        !      call MPI_RECV(imin_recv_buffer, count,MPI_DOUBLE_PRECISION,imin_id,1,MPI_COMM_WORLD,status,ierr)
!        !    else
!        !      print*, Process_id, "imin slave"
!        !      call MPI_RECV(imin_recv_buffer, count,MPI_DOUBLE_PRECISION,imin_id,1,MPI_COMM_WORLD,status,ierr)
!        !      call MPI_SEND(imin_send_buffer, count,MPI_DOUBLE_PRECISION,imin_id,1,MPI_COMM_WORLD, ierr)
!        !    end if
!             ! distribute grid points
!            if(dir_switch(1)==0)then
!              count=0
!              do l=1,layers
!                do k=Gklo(1),Gkhi(1)
!                  do j=Gjlo(1),Gjhi(1)
!                    count=count+1
!                    grid_x(1-l,j,k) = imin_recv_buffer(count)
!                  end do
!                end do
!              end do
!
!              do l=1,layers
!                do k=Gklo(1),Gkhi(1)
!                  do j=Gjlo(1),Gjhi(1)
!                    count=count+1
!                    grid_y(1-l,j,k) = imin_recv_buffer(count)
!                  end do
!                end do
!              end do
!
!              do l=1,layers
!                do k=Gklo(1),Gkhi(1)
!                  do j=Gjlo(1),Gjhi(1)
!                    count=count+1
!                    grid_z(1-l,j,k) = imin_recv_buffer(count)
!                  end do
!                end do
!              end do
!            else
!              count=0
!              do l=1,layers
!                do j=Gjlo(1),Gjhi(1)
!                  do k=Gklo(1),Gkhi(1)
!                    count=count+1
!                    grid_x(1-l,j,k) = imin_recv_buffer(count)
!                  end do
!                end do
!              end do
!
!              do l=1,layers
!                do j=Gjlo(1),Gjhi(1)
!                  do k=Gklo(1),Gkhi(1)
!                    count=count+1
!                    grid_y(1-l,j,k) = imin_recv_buffer(count)
!                  end do
!                end do
!              end do
!
!              do l=1,layers
!                do j=Gjlo(1),Gjhi(1)
!                  do k=Gklo(1),Gkhi(1)
!                    count=count+1
!                    grid_z(1-l,j,k) = imin_recv_buffer(count)
!                  end do
!                end do
!              end do
!            end if
!          end if
!
!          !--- IMAX ---!
!          if(imax_id>=0)then
!            !collect grid point in 1d array
!            count=0
!            do l=1,layers
!              do k=-2,kmx+3
!                do j=-2,jmx+3
!                  count=count+1
!                  imax_send_buffer(count) = grid_x(imx-l,j,k)
!                end do
!              end do
!            end do
!
!            do l=1,layers
!              do k=-2,kmx+3
!                do j=-2,jmx+3
!                  count=count+1
!                  imax_send_buffer(count) = grid_y(imx-l,j,k)
!                end do
!              end do
!            end do
!
!            do l=1,layers
!              do k=-2,kmx+3
!                do j=-2,jmx+3
!                  count=count+1
!                  imax_send_buffer(count) = grid_z(imx-l,j,k)
!                end do
!              end do
!            end do
!
!        call MPI_SENDRECV(imax_send_buffer,count, MPI_DOUBLE_PRECISION, imax_id,1,&
!                          imax_recv_buffer,count, MPI_DOUBLE_PRECISION, imax_id,1,&
!                          MPI_COMM_WORLD,status,ierr)
!        !    if(mpi_class(2)==0)then
!        !      print*, Process_id, "imax master"
!        !      call MPI_SEND(imax_send_buffer, count,MPI_DOUBLE_PRECISION,imax_id,1,MPI_COMM_WORLD, ierr)
!        !      call MPI_RECV(imax_recv_buffer, count,MPI_DOUBLE_PRECISION,imax_id,1,MPI_COMM_WORLD,status,ierr)
!        !    else
!        !      print*, Process_id, "imax slave"
!        !      call MPI_RECV(imax_recv_buffer, count,MPI_DOUBLE_PRECISION,imax_id,1,MPI_COMM_WORLD,status,ierr)
!        !      call MPI_SEND(imax_send_buffer, count,MPI_DOUBLE_PRECISION,imax_id,1,MPI_COMM_WORLD, ierr)
!        !    end if
!             ! distribute grid points
!            if(dir_switch(2)==0)then
!              count=0
!              do l=1,layers
!                do k=Gklo(2),Gkhi(2)
!                  do j=Gjlo(2),Gjhi(2)
!                    count=count+1
!                    grid_x(imx+l,j,k) = imax_recv_buffer(count)
!                  end do
!                end do
!              end do
!
!              do l=1,layers
!                do k=Gklo(2),Gkhi(2)
!                  do j=Gjlo(2),Gjhi(2)
!                    count=count+1
!                    grid_y(imx+l,j,k) = imax_recv_buffer(count)
!                  end do
!                end do
!              end do
!
!              do l=1,layers
!                do k=Gklo(2),Gkhi(2)
!                  do j=Gjlo(2),Gjhi(2)
!                    count=count+1
!                    grid_z(imx+l,j,k) = imax_recv_buffer(count)
!                  end do
!                end do
!              end do
!            else
!              count=0
!              do l=1,layers
!                do j=Gjlo(2),Gjhi(2)
!                  do k=Gklo(2),Gkhi(2)
!                    count=count+1
!                    grid_x(imx+l,j,k) = imax_recv_buffer(count)
!                  end do
!                end do
!              end do
!
!              do l=1,layers
!                do j=Gjlo(2),Gjhi(2)
!                  do k=Gklo(2),Gkhi(2)
!                    count=count+1
!                    grid_y(imx+l,j,k) = imax_recv_buffer(count)
!                  end do
!                end do
!              end do
!
!              do l=1,layers
!                do j=Gjlo(2),Gjhi(2)
!                  do k=Gklo(2),Gkhi(2)
!                    count=count+1
!                    grid_z(imx+l,j,k) = imax_recv_buffer(count)
!                  end do
!                end do
!              end do
!            end if
!          end if

          !--- JMIN ---!
          if(jmin_id>=0)then
            !collect grid point in 1d array
            count=0
            do l=1,layers
              do k=-2,kmx+3
                do i=-2,imx+3
                  count=count+1
                  jmin_send_buffer(count) = grid_x(i,l+1,k)
                end do
              end do
            end do

            do l=1,layers
              do k=-2,kmx+3
                do i=-2,imx+3
                  count=count+1
                  jmin_send_buffer(count) = grid_y(i,l+1,k)
                end do
              end do
            end do

            do l=1,layers
              do k=-2,kmx+3
                do i=-2,imx+3
                  count=count+1
                  jmin_send_buffer(count) = grid_z(i,l+1,k)
                end do
              end do
            end do

            if(mpi_class(3)==0)then
              call MPI_SEND(jmin_send_buffer, count,MPI_DOUBLE_PRECISION,jmin_id,1,MPI_COMM_WORLD, ierr)
              call MPI_RECV(jmin_recv_buffer, count,MPI_DOUBLE_PRECISION,jmin_id,1,MPI_COMM_WORLD,status,ierr)
            else
              call MPI_RECV(jmin_recv_buffer, count,MPI_DOUBLE_PRECISION,jmin_id,1,MPI_COMM_WORLD,status,ierr)
              call MPI_SEND(jmin_send_buffer, count,MPI_DOUBLE_PRECISION,jmin_id,1,MPI_COMM_WORLD, ierr)
            end if
             ! distribute grid points
            if(dir_switch(3)==0)then
              count=0
              do l=1,layers
                do k=Gklo(3),Gkhi(3)
                  do i=Gilo(3),Gihi(3)
                    count=count+1
                    grid_x(i,1-l,k) = jmin_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do k=Gklo(3),Gkhi(3)
                  do i=Gilo(3),Gihi(3)
                    count=count+1
                    grid_y(i,1-l,k) = jmin_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do k=Gklo(3),Gkhi(3)
                  do i=Gilo(3),Gihi(3)
                    count=count+1
                    grid_z(i,1-l,k) = jmin_recv_buffer(count)
                  end do
                end do
              end do
            else
              count=0
              do l=1,layers
                do i=Gilo(3),Gihi(3)
                  do k=Gklo(3),Gkhi(3)
                    count=count+1
                    grid_x(i,1-l,k) = jmin_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do i=Gilo(3),Gihi(3)
                  do k=Gklo(3),Gkhi(3)
                    count=count+1
                    grid_y(i,1-l,k) = jmin_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do i=Gilo(3),Gihi(3)
                  do k=Gklo(3),Gkhi(3)
                    count=count+1
                    grid_z(i,1-l,k) = jmin_recv_buffer(count)
                  end do
                end do
              end do
            end if
          end if

          !--- JMAX ---!
          if(jmax_id>=0)then
            !collect grid point in 1d array
            count=0
            do l=1,layers
              do k=-2,kmx+3
                do i=-2,imx+3
                  count=count+1
                  jmax_send_buffer(count) = grid_x(i,jmx-l,k)
                end do
              end do
            end do

            do l=1,layers
              do k=-2,kmx+3
                do i=-2,imx+3
                  count=count+1
                  jmax_send_buffer(count) = grid_y(i,jmx-l,k)
                end do
              end do
            end do

            do l=1,layers
              do k=-2,kmx+3
                do i=-2,imx+3
                  count=count+1
                  jmax_send_buffer(count) = grid_z(i,jmx-l,k)
                end do
              end do
            end do

            if(mpi_class(4)==0)then
              call MPI_SEND(jmax_send_buffer, count,MPI_DOUBLE_PRECISION,jmax_id,1,MPI_COMM_WORLD, ierr)
              call MPI_RECV(jmax_recv_buffer, count,MPI_DOUBLE_PRECISION,jmax_id,1,MPI_COMM_WORLD,status,ierr)
            else
              call MPI_RECV(jmax_recv_buffer, count,MPI_DOUBLE_PRECISION,jmax_id,1,MPI_COMM_WORLD,status,ierr)
              call MPI_SEND(jmax_send_buffer, count,MPI_DOUBLE_PRECISION,jmax_id,1,MPI_COMM_WORLD, ierr)
            end if
             ! distribute grid points
            if(dir_switch(4)==0)then
              count=0
              do l=1,layers
                do k=Gklo(4),Gkhi(4)
                  do i=Gilo(4),Gihi(4)
                    count=count+1
                    grid_x(i,jmx+l,k) = jmax_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do k=Gklo(4),Gkhi(4)
                  do i=Gilo(4),Gihi(4)
                    count=count+1
                    grid_y(i,jmx+l,k) = jmax_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do k=Gklo(4),Gkhi(4)
                  do i=Gilo(4),Gihi(4)
                    count=count+1
                    grid_z(i,jmx+l,k) = jmax_recv_buffer(count)
                  end do
                end do
              end do
            else
              count=0
              do l=1,layers
                do i=Gilo(4),Gihi(4)
                  do k=Gklo(4),Gkhi(4)
                    count=count+1
                    grid_x(i,jmx+l,k) = jmax_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do i=Gilo(4),Gihi(4)
                  do k=Gklo(4),Gkhi(4)
                    count=count+1
                    grid_y(i,jmx+l,k) = jmax_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do i=Gilo(4),Gihi(4)
                  do k=Gklo(4),Gkhi(4)
                    count=count+1
                    grid_z(i,jmx+l,k) = jmax_recv_buffer(count)
                  end do
                end do
              end do
            end if
          end if
            
          !--- KMIN ---!
          if(kmin_id>=0)then
            !collect grid point in 1d array
            count=0
            do l=1,layers
              do j=-2,jmx+3
                do i=-2,imx+3
                  count=count+1
                  kmin_send_buffer(count) = grid_x(i,j,1+l)
                end do
              end do
            end do

            do l=1,layers
              do j=-2,jmx+3
                do i=-2,imx+3
                  count=count+1
                  kmin_send_buffer(count) = grid_y(i,j,1+l)
                end do
              end do
            end do

            do l=1,layers
              do j=-2,jmx+3
                do i=-2,imx+3
                  count=count+1
                  kmin_send_buffer(count) = grid_z(i,j,1+l)
                end do
              end do
            end do

            if(mpi_class(5)==0)then
              call MPI_SEND(kmin_send_buffer, count,MPI_DOUBLE_PRECISION,kmin_id,1,MPI_COMM_WORLD, ierr)
              call MPI_RECV(kmin_recv_buffer, count,MPI_DOUBLE_PRECISION,kmin_id,1,MPI_COMM_WORLD,status,ierr)
            else
              call MPI_RECV(kmin_recv_buffer, count,MPI_DOUBLE_PRECISION,kmin_id,1,MPI_COMM_WORLD,status,ierr)
              call MPI_SEND(kmin_send_buffer, count,MPI_DOUBLE_PRECISION,kmin_id,1,MPI_COMM_WORLD, ierr)
            end if
             ! distribute grid points
            if(dir_switch(5)==0)then
              count=0
              do l=1,layers
                do j=Gjlo(5),Gjhi(5)
                  do i=Gilo(5),Gihi(5)
                    count=count+1
                    grid_x(i,j,1-l) = kmin_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do j=Gjlo(5),Gjhi(5)
                  do i=Gilo(5),Gihi(5)
                    count=count+1
                    grid_y(i,j,1-l) = kmin_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do j=Gjlo(5),Gjhi(5)
                  do i=Gilo(5),Gihi(5)
                    count=count+1
                    grid_z(i,j,1-l) = kmin_recv_buffer(count)
                  end do
                end do
              end do
            else
              count=0
              do l=1,layers
                do i=Gilo(5),Gihi(5)
                  do j=Gjlo(5),Gjhi(5)
                    count=count+1
                    grid_x(i,j,1-l) = kmin_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do i=Gilo(5),Gihi(5)
                  do j=Gjlo(5),Gjhi(5)
                    count=count+1
                    grid_y(i,j,1-l) = kmin_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do i=Gilo(5),Gihi(5)
                  do j=Gjlo(5),Gjhi(5)
                    count=count+1
                    grid_z(i,j,1-l) = kmin_recv_buffer(count)
                  end do
                end do
              end do
            end if
          end if

          !--- KMAX ---!
          if(kmax_id>=0)then
            !collect grid point in 1d array
            count=0
            do l=1,layers
              do j=-2,jmx+3
                do i=-2,imx+3
                  count=count+1
                  kmax_send_buffer(count) = grid_x(i,j,kmx-l)
                end do
              end do
            end do

            do l=1,layers
              do j=-2,jmx+3
                do i=-2,imx+3
                  count=count+1
                  kmax_send_buffer(count) = grid_y(i,j,kmx-l)
                end do
              end do
            end do

            do l=1,layers
              do j=-2,jmx+3
                do i=-2,imx+3
                  count=count+1
                  kmax_send_buffer(count) = grid_z(i,j,kmx-l)
                end do
              end do
            end do

            if(mpi_class(6)==0)then
              call MPI_SEND(kmax_send_buffer, count,MPI_DOUBLE_PRECISION,kmax_id,1,MPI_COMM_WORLD, ierr)
              call MPI_RECV(kmax_recv_buffer, count,MPI_DOUBLE_PRECISION,kmax_id,1,MPI_COMM_WORLD,status,ierr)
            else
              call MPI_RECV(kmax_recv_buffer, count,MPI_DOUBLE_PRECISION,kmax_id,1,MPI_COMM_WORLD,status,ierr)
              call MPI_SEND(kmax_send_buffer, count,MPI_DOUBLE_PRECISION,kmax_id,1,MPI_COMM_WORLD, ierr)
            end if
             ! distribute grid points
            if(dir_switch(6)==0)then
              count=0
              do l=1,layers
                do j=Gjlo(6),Gjhi(6)
                  do i=Gilo(6),Gihi(6)
                    count=count+1
                    grid_x(i,j,kmx+l) = kmax_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do j=Gjlo(6),Gjhi(6)
                  do i=Gilo(6),Gihi(6)
                    count=count+1
                    grid_y(i,j,kmx+l) = kmax_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do j=Gjlo(6),Gjhi(6)
                  do i=Gilo(6),Gihi(6)
                    count=count+1
                    grid_z(i,j,kmx+l) = kmax_recv_buffer(count)
                  end do
                end do
              end do
            else
              count=0
              do l=1,layers
                do i=Gilo(6),Gihi(6)
                  do j=Gjlo(6),Gjhi(6)
                    count=count+1
                    grid_x(i,j,kmx+l) = kmax_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do i=Gilo(6),Gihi(6)
                  do j=Gjlo(6),Gjhi(6)
                    count=count+1
                    grid_y(i,j,kmx+l) = kmax_recv_buffer(count)
                  end do
                end do
              end do

              do l=1,layers
                do i=Gilo(6),Gihi(6)
                  do j=Gjlo(6),Gjhi(6)
                    count=count+1
                    grid_z(i,j,kmx+l) = kmax_recv_buffer(count)
                  end do
                end do
              end do
            end if
          end if
            
          deallocate(imin_send_buffer)
          deallocate(jmin_send_buffer)
          deallocate(kmin_send_buffer)
          deallocate(imin_recv_buffer)
          deallocate(jmin_recv_buffer)
          deallocate(kmin_recv_buffer)
          deallocate(imax_send_buffer)
          deallocate(jmax_send_buffer)
          deallocate(kmax_send_buffer)
          deallocate(imax_recv_buffer)
          deallocate(jmax_recv_buffer)
          deallocate(kmax_recv_buffer)
          call mpi_barrier(MPI_COMM_WORLD,ierr)
          call dmsg(1, 'grid', 'done with ghost_grid')
        end subroutine ghost_grid
end module grid
