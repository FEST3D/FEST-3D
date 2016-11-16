module grid
    !-------------------------------------------------------------------
    ! The grid module contains the grid definition (locations of the 
    ! grid points) as well as procedures to load these from a file.
    !
    ! This version of the grid module only allows for a 1-dimensional
    ! or 2-dimensional grid.
    !-------------------------------------------------------------------
    
    use global, only: STRING_BUFFER_LENGTH, GRID_FILE_UNIT
    use utils, only: alloc, dealloc, dmsg
    
    implicit none
    private
    
    ! Grid point coordinates
    real, public, dimension(:, :), allocatable :: grid_x, grid_y
    ! Grid size
    integer, public :: imx, jmx

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

            call alloc(grid_x, 1, imx, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for grid_x.')
            call alloc(grid_y, 1, imx, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for grid_y.')

        end subroutine allocate_memory

        subroutine destroy_grid()
            !-----------------------------------------------------------
            ! Deallocate the memory allocated for the grid.
            !-----------------------------------------------------------

            implicit none

            call dmsg(1, 'grid', 'destroy_memory')

            call dealloc(grid_x)
            call dealloc(grid_y)

        end subroutine destroy_grid

        subroutine setup_grid(gridfile)
            !-----------------------------------------------------------
            ! Read the grid file and initialize the grid
            !-----------------------------------------------------------

            implicit none
            character(len=32), intent(in) :: gridfile
            
            call dmsg(1, 'grid', 'setup_grid')

            open(GRID_FILE_UNIT, file=gridfile)

            call extract_grid_size()

            call allocate_memory()

            call populate_grid_points()

            close(GRID_FILE_UNIT)
        
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
            read(header, *, iostat=ios) imx, jmx
            if (ios /= 0) then
                ! An io error means it was not possible to read jmx.
                ! This means the file does not have a jmx and so, set
                ! the extent of this direction to 1. Read the remaining
                ! dimension from the header.
                read(header, *, iostat=ios) imx
                if (ios /= 0) then
                    ! There was an error reading the extent.
                    ! This is an error.
                    print *, 'Unable to read grid extent.'
                    stop
                end if
                jmx = 1
            end if
        end subroutine extract_grid_size

        subroutine extract_grid_point(line, i, j)
            !-----------------------------------------------------------
            ! Extract a grid point from a line of the grid file. 
            !-----------------------------------------------------------

            implicit none
            character(len=STRING_BUFFER_LENGTH), intent(in) :: line
            integer, intent(in) :: i, j

            call dmsg(0, 'grid', 'extract_grid_point')

            if (jmx > 1) then
                read(line, *) grid_x(i, j), grid_y(i, j)
            else
                read(line, *) grid_x(i, j)
                grid_y(i, j) = 0.
            end if
        end subroutine extract_grid_point

        subroutine populate_grid_points()
            !-----------------------------------------------------------
            ! Use the grid file to populate the grid points.
            !-----------------------------------------------------------

            implicit none
            character(len=STRING_BUFFER_LENGTH) :: line
            integer :: i, j
            integer :: ios  ! io status

            call dmsg(1, 'grid', 'populate_grid_point')
            !print *, imx, jmx

            ! Read grid points from the grid file
            do j = 1, jmx
                do i = 1, imx
                    read(GRID_FILE_UNIT, '(A)', iostat=ios) line
                    if (ios /= 0) then
                        print *, 'Error while reading grid line.'
                        print *, 'Current grid point: ', i, j
                        print *, 'Current buffer length is set to: ', &
                                STRING_BUFFER_LENGTH
                        print *, 'Exiting program.'
                        stop
                    end if
                    call extract_grid_point(line, i, j)
                end do
            end do

        end subroutine populate_grid_points

end module grid
