module grid
    !-------------------------------------------------------------------
    ! The grid module contains the grid definition (locations of the 
    ! grid points) as well as procedures to load these from a file.
    !
    ! This version of the grid module only allows for a 1-dimensional
    ! or 2-dimensional grid.
    !-------------------------------------------------------------------
    
    use global, only: STRING_BUFFER_LENGTH, GRID_FILE_UNIT
!                     SPHERE_INDICES_FILE_UNIT
    use utils, only: alloc, dealloc, dmsg
    
    implicit none
    private
    
    ! Grid point coordinates
    real, public, dimension(:, :, :), allocatable :: grid_x, grid_y, grid_z
    ! Grid size
    integer, public :: imx, jmx, kmx
!   integer, public, dimension(:, :), allocatable :: sphere_indices
!   integer, public :: n_sph_ind

    ! Public methods
    public :: setup_grid
    public :: destroy_grid

    contains
    
        subroutine allocate_memory()
            !-----------------------------------------------------------
            ! Allocate memory to store the grid
            !-----------------------------------------------------------

            implicit none

            call dmsg(1, 'grid', 'allocate_memory')

            call alloc(grid_x, 1, imx, 1, jmx, 1, kmx, &
                    errmsg='Error: Unable to allocate memory for grid_x.')
            call alloc(grid_y, 1, imx, 1, jmx, 1, kmx, &
                    errmsg='Error: Unable to allocate memory for grid_y.')
            call alloc(grid_z, 1, imx, 1, jmx, 1, kmx, &
                    errmsg='Error: Unable to allocate memory for grid_z.')
            
            ! The alloc function earlier allocates only real arrays.
            ! A new function has been written to allocate for integers
            ! as well

!           call alloc(sphere_indices, 1, 3, 1, n_sph_ind, &
!                   errmsg='Error: Unable to allocate memory for sphere_indices.')
            
        end subroutine allocate_memory

        subroutine destroy_grid()
            !-----------------------------------------------------------
            ! Deallocate the memory allocated for the grid.
            !-----------------------------------------------------------

            implicit none

            call dmsg(1, 'grid', 'destroy_memory')

            call dealloc(grid_x)
            call dealloc(grid_y)
            call dealloc(grid_z)

!           call dealloc(sphere_indices)

        end subroutine destroy_grid

        subroutine setup_grid(gridfile)
            !-----------------------------------------------------------
            ! Read the grid file and initialize the grid
            !-----------------------------------------------------------

            implicit none
            character(len=32), intent(in) :: gridfile
   !        character(len=32) :: sphindfile
            
            call dmsg(1, 'grid', 'setup_grid')
!           sphindfile = 'sphere-indices.txt'

            open(GRID_FILE_UNIT, file=gridfile)

!           open(SPHERE_INDICES_FILE_UNIT, file=sphindfile)

            call extract_grid_size()
!           call extract_sphere_indices_size()

            call allocate_memory()

            call populate_grid_points()
!           call populate_sphere_indices()

            close(GRID_FILE_UNIT)
!           close(SPHERE_INDICES_FILE_UNIT)
        
        end subroutine setup_grid
        
        subroutine extract_grid_size()
            !-----------------------------------------------------------
            ! Extract the grid size from the grid file header
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

!       subroutine extract_sphere_indices_size
!       
!           implicit none
!           integer :: ios
!           character(len=STRING_BUFFER_LENGTH) :: header
!           
!           call dmsg(1, 'grid', 'extract_sphere_indices_size')
!           
!     
!           read(SPHERE_INDICES_FILE_UNIT, '(A)', iostat=ios) header
!           if (ios /= 0) then
!               ! We have an error
!               print *, 'Error reading sphere indices file. Aborting..'
!               stop
!           end if

!           read(header, *, iostat=ios) n_sph_ind
!           if (ios /= 0) then
!               ! We have an error
!               print *, 'Error reading number of sphere indices. Aborting..'
!               stop
!           end if

!       end subroutine extract_sphere_indices_size

        subroutine extract_grid_point(line, i, j, k)
            !-----------------------------------------------------------
            ! Extract a grid point from a line of the grid file. 
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
            !-----------------------------------------------------------
            ! Use the grid file to populate the grid points.
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

!       subroutine populate_sphere_indices()
!       
!           implicit none
!           integer :: i
!           integer :: ios
!           character(len=STRING_BUFFER_LENGTH) :: line
!           
!           call dmsg(1, 'grid', 'populate_sphere_indices')

!        !  print *, 'Number sphere_indices: ', n_sph_ind

!           do i = 1, n_sph_ind
!               read(SPHERE_INDICES_FILE_UNIT, '(A)', iostat=ios) line
!               read(line, *) sphere_indices(1, i), sphere_indices(2, i), &
!                             sphere_indices(3, i)
!               if (ios /= 0) then
!                   print *, 'Error while reading line'
!                   print *, 'Line number: ', i
!                   stop
!               end if
!           end do
!       
!       end subroutine populate_sphere_indices

end module grid
