 !< Detect all the grid points on the wall boundary condition
module wall
 !< Detect all the grid points on the wall boundary condition
 !< and store them in a single file
  use vartypes
  use mpi
  use utils, only: alloc

#include "../debug.h"
#include "../error.h"

  private
  integer :: ierr
  !< Integer to store error 
  integer :: new_type
  !< Create new type for MPI
  integer :: thisfile
  !< File hadler
  integer, parameter :: maxlen=70
  !< Maximum length for string

  real(wp), private, dimension(:, :), allocatable, target :: wallc 
  !< Centre of wall surface
  real(wp), private, dimension(:), pointer :: wall_x 
  !< X coordiante of center of wall surface
  real(wp), private, dimension(:), pointer :: wall_y 
  !< Y coordiante of center of wall surface
  real(wp), private, dimension(:), pointer :: wall_z 
  !< Z coordiante of center of wall surface
  integer, dimension(6) :: no_slip_flag=0 
  !< Flag to detect wall
  integer, public :: n_wall
  !< Number of points on the wall
  integer, public :: total_n_wall
  !< Total number of points on the block across all processes
  character(len=maxlen), dimension(:), allocatable :: str
  !< Store all wall corridnate of current process in a string vector
  character(len=maxlen) :: line
  !< Line to write in output file
  character , parameter :: lf=Achar(10)
  !< End of line character
  

  ! For gather all the data to process 0
  integer, dimension(:), allocatable :: n_wall_buf
  !< Store n_wall points of all processors in a array form
  integer, dimension(:), allocatable :: write_flag
  !< Check if current processor has any wall points to write

  integer :: imx, jmx, kmx
  public :: write_surfnode

  contains 

    subroutine write_surfnode(files, nodes, control, bc, dims)
      !< Extract and write the wall surface node points
      !< in a file shared by all the MPI processes
      implicit none
      type(filetype), intent(in) :: files
      type(controltype), intent(in) :: control
      type(extent), intent(in) :: dims
      type(boundarytype), intent(in) :: bc
      type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(in) :: nodes
      integer :: count
      integer :: i

      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx

      call setup_surface(files, control, bc)
      call surface_points(nodes)
      call MPI_GATHER(n_wall, 1, MPI_Integer, n_wall_buf, 1, &
                      MPI_integer,0, MPI_COMM_WORLD, ierr)
      total_n_wall = sum(n_wall_buf(:))
      call MPI_Bcast(total_n_wall,1, MPI_Integer, 0, &
                       MPI_COMM_WORLD, ierr)
      call MPI_Bcast(n_wall_buf, control%total_process, MPI_Integer, 0, &
                       MPI_COMM_WORLD, ierr)

      write_flag=0
      count=0
      do i=1,control%total_process
        if(n_wall_buf(i)>0) then
          write_flag(i)=count
          count = count+1
        end if
      end do
      call MPI_TYPE_CONTIGUOUS(maxlen,MPI_Character, new_type, ierr)
      call MPI_TYPE_COMMIT(new_type, ierr)
      if(process_id==0)then
      write(line, '(I0)') total_n_wall
      line(len(line):len(line))=lf
      call MPI_FILE_WRITE_shared(thisfile, line, 1, &
                              new_type, &
                              MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(n_wall>0)then
      do i=1,n_wall
        write(line, '(3(ES18.10E3,4x))') wall_x(i), wall_y(i), wall_z(i)
        line(len(line):len(line))=lf
        call MPI_FILE_WRITE_shared(thisfile, line, 1, &
                                new_type, &
                                MPI_STATUS_IGNORE, ierr)
      end do
      end if
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call MPI_FILE_CLOSE(thisfile, ierr)

    end subroutine write_surfnode


    subroutine allocate_memory(control)
      !< Allocate memory to str and wallc variable array
        implicit none
        type(controltype), intent(in) :: control

        DebugCall('setup_surface')

        n_wall = -1

        n_wall =  ((jmx)*(kmx)*(NO_SLIP_flag(1)) &
                  +(jmx)*(kmx)*(NO_SLIP_flag(2)) &
                  +(kmx)*(imx)*(NO_SLIP_flag(3)) &
                  +(kmx)*(imx)*(NO_SLIP_flag(4)) &
                  +(imx)*(jmx)*(NO_SLIP_flag(5)) &
                  +(imx)*(jmx)*(NO_SLIP_flag(6)) &
                  )

        allocate(str(1:n_wall))
        call alloc(wallc, 1, n_wall, 1, 3 ,&
                  errmsg='Error: Unable to allocate memory for wallc')

        allocate(n_wall_buf(1:control%total_process))
        allocate(write_flag(1:control%total_process))
    end subroutine allocate_memory



    subroutine link_aliases()
      !< Link pointers wall_x, wall_y, wall_z to wallc

      implicit none

      DebugCall('link_aliases')
      wall_x(1:n_wall) => wallc(1:n_wall,1)
      wall_y(1:n_wall) => wallc(1:n_wall,2)
      wall_z(1:n_wall) => wallc(1:n_wall,3)

    end subroutine link_aliases


    

    subroutine setup_surface(files, control, bc)
      !< Open MPI_shared write file, allocate memory and
      !< setup pointers

      implicit none
      type(boundarytype), intent(in) :: bc
      type(filetype), intent(in) :: files
      type(controltype), intent(in) :: control
      integer :: stat

      DebugCall('setup_surface')
      if(process_id==0)then
      open(files%NODESURF_FILE_UNIT, file=files%surface_node_points, iostat=stat)
      if(stat==0) close(files%NODESURF_FILE_UNIT, status='delete')
      end if
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call MPI_FILE_OPEN(MPI_COMM_WORLD, files%surface_node_points, &
                        MPI_MODE_WRONLY + MPI_MODE_CREATE + MPI_MODE_EXCL, &
                        MPI_INFO_NULL, thisfile, ierr)
      call find_wall(bc)
      call allocate_memory(control)
      call link_aliases()

    end subroutine setup_surface


    subroutine find_wall(bc)
      !< Setup wall flag for all six boundary of the block

      implicit none
      type(boundarytype), intent(in) :: bc
      integer :: i

      NO_slip_flag=0
      do i = 1,6
        if(bc%id(i)==-5) NO_SLIP_FLAG(i)=1
      end do

    end subroutine find_wall



    subroutine surface_points(nodes)
      !< Extract surface points and store them
      !< in a string vector str(ind)


      implicit none
      type(nodetype), dimension(-2:imx+3,-2:jmx+3,-2:kmx+3), intent(in) :: nodes
      integer :: OL
      integer :: i, j, k, ind
      integer :: im=1, ix=1
      integer :: jm=1, jx=1
      integer :: km=1, kx=1

      DebugCall('surface_points')


      ind = 0

      do OL = 1,6

        if (NO_SLIP_flag(OL) == 1 )then
          select case (OL)
            case (1)
              km = 1
              jm = 1
              im = 1
              kx = kmx
              jx = jmx
              ix = 1
            case (2)
              km = 1
              jm = 1
              im = imx
              kx = kmx
              jx = jmx
              ix = imx
            case (3)
              km = 1
              jm = 1
              im = 1
              kx = kmx
              jx = 1
              ix = imx
            case (4)
              km = 1
              jm = jmx
              im = 1
              kx = kmx
              jx = jmx
              ix = imx
            case (5)
              km = 1
              jm = 1
              im = 1
              kx = 1
              jx = jmx
              ix = imx
            case (6)
              km = kmx
              jm = 1
              im = 1
              kx = kmx
              jx = jmx
              ix = imx
            case DEFAULT
              Fatal_error
              km = 1
              jm = 1
              im = 1
              kx = -1
              jx = -1
              ix = -1
          end select

        do k = km,kx
          do j = jm,jx
            do i = im,ix 
              ind = ind + 1
              wall_x(ind) = nodes(i, j, k )%x
              wall_y(ind) = nodes(i, j, k )%y
              wall_z(ind) = nodes(i, j, k )%z
              write(str(ind),'(3(f0.16, 4x))') wall_x(ind), wall_y(ind), wall_z(ind)
            end do
          end do
        end do

        end if

    end do

    end subroutine surface_points

end module wall
