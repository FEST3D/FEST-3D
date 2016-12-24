module surfnode
  use global,     only: STRING_BUFFER_LENGTH
  use global,     only: FILE_NAME_LENGTH
  use global,     only: LAYOUT_FILE_UNIT
  use global,     only: NODESURF_FILE_UNIT
  use global,     only: TEMP_NODE_FILE_UNIT
  use global,     only: nodefile_temp
  use global,     only: surface_node_points
  use grid,       only: setup_grid
  use grid,       only: destroy_grid                        
  use wall_find,  only: setup_surface 
  use wall_find,  only: destroy_surface
  use wall_find,  only: surface_points
  use wall_find,  only: n_wall

  implicit none
  private

  integer                                                    :: n_total_surfnodes
  integer                                                    :: Nblocks
  integer, dimension(:), allocatable                         :: n_surfnodes
  real, dimension(:), allocatable                            :: surface_grid_all
  character(len=FILE_NAME_LENGTH), dimension(:), allocatable :: gridfiles
  character(len=FILE_NAME_LENGTH), dimension(:), allocatable :: bcfiles
  public                                                     :: extract_surfnodes
  contains


    subroutine setup_tempfile()
      !---------------------------------------------------------------------
      ! This subroutine creats a temp file which will store all the
      ! grid points which lies on no slip boundary conditoion. Files
      ! is just opened and closed here to just create a file in directory
      ! which will be futher used by wall_find module to append data.
      ! This file will be destroy on complition of extraciton process.
      !--------------------------------------------------------------------
      open(TEMP_NODE_FILE_UNIT, file=nodefile_temp)
      close(TEMP_NODE_FILE_UNIT)
    end subroutine setup_tempfile

    subroutine destroy_tempfile()
      open(TEMP_NODE_FILE_UNIT, file=nodefile_temp, status='old')
      close(TEMP_NODE_FILE_UNIT, status='delete')
    end subroutine destroy_tempfile



    subroutine setup_nodefile()
      open(NODESURF_FILE_UNIT, file=surface_node_points)
      allocate(surface_grid_all(1:n_total_surfnodes))
    end subroutine setup_nodefile

    subroutine destroy_nodefile()
      close(NODESURF_FILE_UNIT)
      deallocate(surface_grid_all)
    end subroutine destroy_nodefile


    subroutine get_next_token_parallel(buf)
      !-----------------------------------------------------------
      ! Extract the next token from the layout file
      !
      ! Each token is on a separate line.
      ! There may be multiple comments (lines beginning with #) 
      ! and blank lines in between.
      ! The purpose of this subroutine is to ignore all these 
      ! lines and return the next "useful" line.
      !-----------------------------------------------------------
  
      implicit none
      character(len=STRING_BUFFER_LENGTH), intent(out) :: buf
      integer :: ios
  
      do
         read(LAYOUT_FILE_UNIT, '(A)', iostat=ios) buf
         if (ios /= 0) then
            print *, 'Error while reading config file.'
            print *, 'Current buffer length is set to: ', &
                 STRING_BUFFER_LENGTH
            stop
         end if
         if (index(buf, '#') == 1) then
            ! The current line begins with a hash
            ! Ignore it
            continue
         else if (len_trim(buf) == 0) then
            ! The current line is empty
            ! Ignore it
            continue
         else
            ! A new token has been found
            ! Break out
            exit
         end if
      end do
!      call dmsg(0, 'solver', 'get_next_token', 'Returning: ' // trim(buf))
  
    end subroutine get_next_token_parallel


    subroutine read_layout_file()
      implicit none
      integer                             :: i
      character(len=STRING_BUFFER_LENGTH) :: buf
      character(len=STRING_BUFFER_LENGTH) :: dump

      open(LAYOUT_FILE_UNIT, file = "layout/layout.md")

      ! Read the parameters from the file
      call get_next_token_parallel(buf)
      read(buf,*) Nblocks !total number of blocks
      call get_next_token_parallel(buf)
      read(buf,*) dump    !no of entry per block, not required

      allocate(gridfiles(Nblocks))
      allocate(bcfiles(Nblocks))
      !reading block data
      do i = 1,Nblocks
      call get_next_token_parallel(buf)
      read(buf,*) dump, gridfiles(i), bcfiles(i), dump
      write(gridfiles(i), '(A)') 'gridfiles/'//trim(gridfiles(i))
      write(bcfiles(i),       '(A)') 'bc/'//trim(bcfiles(i))
      end do

    end subroutine read_layout_file


    subroutine setup_destroy_surfnodes()
      implicit none
      integer :: i

      do i = 1,Nblocks
        call setup_grid(gridfiles(i))
        call setup_surface(bcfiles(i))
        n_surfnodes(i) = n_wall
        call surface_points()
        call destroy_surface()
        call destroy_grid()
      end do

    end subroutine setup_destroy_surfnodes

    subroutine write_surfnodes()
      implicit none
      integer   :: i
      real      :: x,y,z

      close(TEMP_NODE_FILE_UNIT)
      open(TEMP_NODE_FILE_UNIT, file=nodefile_temp, status='old')
      write(NODESURF_FILE_UNIT,'(I0)') n_total_surfnodes
      do i = 1,n_total_surfnodes
        read(TEMP_NODE_FILE_UNIT, *) x, y, z
        write(NODESURF_FILE_UNIT, '(3(f0.16, A))') x,' ',y,' ',z,' '
      end do
      close(TEMP_NODE_FILE_UNIT)

    end subroutine write_surfnodes

    subroutine extract_surfnodes()
      call read_layout_file()
      call setup_tempfile()
      allocate(n_surfnodes(Nblocks))
      call setup_destroy_surfnodes()
      n_total_surfnodes = sum(n_surfnodes)
      call setup_nodefile()
      call write_surfnodes()
      call destroy_tempfile()
      call destroy_nodefile()
    end subroutine extract_surfnodes

end module surfnode




