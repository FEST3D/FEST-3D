module wall
  use global     , only : surface_node_points
  use global     , only : NODESURF_FILE_UNIT
  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx
  use global_vars, only : grid_x
  use global_vars, only : grid_y
  use global_vars, only : grid_z
  use global_vars, only : process_id
  use global_vars, only : total_process

  use string
  use bitwise
  use utils, only: alloc, dealloc, dmsg, DEBUG_LEVEL
  use boundary_conditions, only : bc_imn
  use boundary_conditions, only : bc_imx
  use boundary_conditions, only : bc_jmn
  use boundary_conditions, only : bc_jmx
  use boundary_conditions, only : bc_kmn
  use boundary_conditions, only : bc_kmx
  use boundary_conditions, only : BC_NO_SLIP

#ifdef __GFORTRAN__
  use mpi
#endif
  implicit none
#ifdef __INTEL_COMPILER
  include "mpif.h"
#endif

  private
  integer :: ierr
  real :: buf
  integer :: BUFSIZE
  integer :: new_type
  integer :: thisfile
  integer, parameter :: maxlen=70

  real, private, dimension(:, :), allocatable, target :: wallc ! centre of wall surface
  real, private, dimension(:), pointer :: wall_x
  real, private, dimension(:), pointer :: wall_y
  real, private, dimension(:), pointer :: wall_z
  integer, dimension(6) :: no_slip_flag=0
  integer, public :: n_wall
  integer, public :: total_n_wall
  character(len=maxlen), dimension(:), allocatable :: str
  character(len=maxlen) :: line
  character , parameter :: lf=Achar(10)
  

  ! for gather all the data to process 0
  integer, dimension(:), allocatable :: n_wall_buf
  integer, dimension(:), allocatable :: write_flag
  real, dimension(:), allocatable :: wall_buf

  public :: write_surfnode

  contains 

    subroutine write_surfnode()
      implicit none
      integer :: count
      integer :: i
      call setup_surface()
      call surface_points()
      call MPI_GATHER(n_wall, 1, MPI_Integer, n_wall_buf, 1, &
                      MPI_integer,0, MPI_COMM_WORLD, ierr)
      total_n_wall = sum(n_wall_buf(:))
      call MPI_Bcast(total_n_wall,1, MPI_Integer, 0, &
                       MPI_COMM_WORLD, ierr)
      call MPI_Bcast(n_wall_buf, total_process, MPI_Integer, 0, &
                       MPI_COMM_WORLD, ierr)

      write_flag=0
      count=0
      do i=1,total_process
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
      call destroy_surface()

    end subroutine write_surfnode


    subroutine allocate_memory()
        implicit none

        call dmsg(1, 'wall_find', 'setup_surface')

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

        allocate(n_wall_buf(1:total_process))
        allocate(write_flag(1:total_process))
    end subroutine allocate_memory



    subroutine link_aliases()

      implicit none

      call dmsg(1, 'wall_find', 'link_aliases')
      wall_x(1:n_wall) => wallc(1:n_wall,1)
      wall_y(1:n_wall) => wallc(1:n_wall,2)
      wall_z(1:n_wall) => wallc(1:n_wall,3)

    end subroutine link_aliases



    subroutine unlink_aliases()

      implicit none

      call dmsg(1, 'wall_find', 'unlink_aliases')
      nullify(wall_x)
      nullify(wall_y)
      nullify(wall_z)

    end subroutine unlink_aliases



    subroutine deallocate_memory()

      implicit none

      call dmsg(1, 'wall_find', 'dealloate_memory')
      call dealloc(wallc)

    end subroutine deallocate_memory

    

    subroutine setup_surface()

      implicit none
      integer :: stat

      call dmsg(1, 'wall_find', 'setup_surface')
      if(process_id==0)then
      open(NODESURF_FILE_UNIT, file=surface_node_points, iostat=stat)
      if(stat==0) close(NODESURF_FILE_UNIT, status='delete')
      end if
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call MPI_FILE_OPEN(MPI_COMM_WORLD, surface_node_points, &
                        MPI_MODE_WRONLY + MPI_MODE_CREATE + MPI_MODE_EXCL, &
                        MPI_INFO_NULL, thisfile, ierr)
      call find_wall()
      call allocate_memory()
      call link_aliases()

    end subroutine setup_surface



    subroutine destroy_surface()

      implicit none

      call dmsg(1, 'wall_find', 'destroy_surface')
      call deallocate_memory()
      call unlink_aliases()
      call MPI_FILE_CLOSE(thisfile, ierr)

    end subroutine destroy_surface


    subroutine find_wall()

      implicit none


      if    ((bc_imn(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
        NO_slip_flag(1)=1
      elseif((bc_imx(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
        NO_slip_flag(2)=1
      elseif((bc_jmn(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
        NO_slip_flag(3)=1
      elseif((bc_jmx(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
        NO_slip_flag(4)=1
      elseif((bc_kmn(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
        NO_slip_flag(5)=1
      elseif((bc_kmx(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
        NO_slip_flag(6)=1
      end if

    end subroutine find_wall



    subroutine surface_points()

      implicit none
      integer :: OL
      integer :: i, j, k, ind
      integer :: im=1, ix=1
      integer :: jm=1, jx=1
      integer :: km=1, kx=1

      call dmsg(1, 'wall_find', 'surface_points')


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
              call dmsg(5, "wall_find", 'Surface_points', 'FATAL  ERROR: select case')
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
              wall_x(ind) = grid_x(i, j, k )        
              wall_y(ind) = grid_y(i, j, k )        
              wall_z(ind) = grid_z(i, j, k )        
              write(str(ind),'(3(f0.16, 4x))') wallc(i,:)
            end do
          end do
        end do

        end if

    end do

    end subroutine surface_points

end module wall
