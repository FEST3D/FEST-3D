module surface
  use global, only: CONFIG_FILE_UNIT, FILE_NAME_LENGTH, STRING_BUFFER_LENGTH,&
    BOUNDARY_CONDITIONS_FILE_UNIT
  use bitwise
  use utils, only: alloc, dealloc, dmsg, DEBUG_LEVEL
  use grid, only: imx, jmx, kmx, grid_x, grid_y, grid_z

  implicit none
  private

  real, public, dimension(:, :), allocatable, target :: wallc ! centre of wall surface
  real, public, dimension(:), pointer :: wall_x
  real, public, dimension(:), pointer :: wall_y
  real, public, dimension(:), pointer :: wall_z
 ! real , dimension(:,:,:), allocatable, pointer :: face_points()
  integer, dimension(6) :: NO_SLIP_flag
  integer :: n_wall


  public :: setup_surface
  public :: destroy_surface
  public :: surface_points   

  contains


        subroutine allocate_memory()
            implicit none

            call dmsg(1, 'surface', 'setup_surface')

            n_wall = 0

            n_wall =  ((jmx-1)*(kmx-1)*(NO_SLIP_flag(1)) &
                      +(jmx-1)*(kmx-1)*(NO_SLIP_flag(2)) &
                      +(kmx-1)*(imx-1)*(NO_SLIP_flag(3)) &
                      +(kmx-1)*(imx-1)*(NO_SLIP_flag(4)) &
                      +(imx-1)*(jmx-1)*(NO_SLIP_flag(5)) &
                      +(imx-1)*(jmx-1)*(NO_SLIP_flag(6)) &
                      )

            call alloc(wallc, 1, n_wall, 1, 3 ,&
                      errmsg='Error: Unable to allocate memory for wallc')

        end subroutine allocate_memory



        subroutine link_aliases()

          implicit none

          call dmsg(1, 'surface', 'link_aliases')
          wall_x(1:n_wall) => wallc(1:n_wall,1)
          wall_y(1:n_wall) => wallc(1:n_wall,2)
          wall_z(1:n_wall) => wallc(1:n_wall,3)

        end subroutine link_aliases



        subroutine unlink_aliases()

          implicit none

          call dmsg(1, 'surface', 'unlink_aliases')
          nullify(wall_x)
          nullify(wall_y)
          nullify(wall_z)

        end subroutine unlink_aliases
  
  
  
        subroutine deallocate_memory()
  
          implicit none
  
          call dmsg(1, 'surface', 'dealloate_memory')
          call dealloc(wallc)
  
        end subroutine deallocate_memory
  
        

        subroutine setup_surface()
  
          implicit none
  
          call dmsg(1, 'surface', 'setup_surface')
          call find_wall()
          call allocate_memory()
          call link_aliases()
  
        end subroutine setup_surface
  

  
        subroutine destroy_surface()
  
          implicit none
  
          call dmsg(1, 'surface', 'destroy_surface')
          call deallocate_memory()
          call unlink_aliases()
  
        end subroutine destroy_surface



        subroutine update_flag(i)

          implicit none


          integer, intent(in) :: i !index of array
          integer :: ios
          character(len=STRING_BUFFER_LENGTH) :: buf  

          call dmsg(1, 'surface', 'update_flag')
          buf = "none_zero_length"
          do while ( len_trim(buf) /= 0 )
            read(BOUNDARY_CONDITIONS_FILE_UNIT, '(A)', iostat =ios)  buf
            if(ios /=0 )then
              print*, "  /!\"
              print*, " /_|_\ ERROR: Boundary conditon read unsuccesful "
              print*, " checkpoint : module_surface - subroutine_ update_flag"
              STOP
            end if
  
            print*, buf(3:index(buf(3:), ' ')+1)
            if( buf(3:index(buf(3:), ' ')+1)  == "NO_SLIP" )then
              NO_SLIP_flag(i) = 1
            end if
          end do

        end subroutine update_flag



        subroutine find_wall()

          implicit none


          integer :: ios
          character(len=STRING_BUFFER_LENGTH) :: buf 

          call dmsg(1, 'surface', 'find_wall')
          
          NO_SLIP_flag = 0

          ! opening the boundary condtion file (bc.config.d)
          open(BOUNDARY_CONDITIONS_FILE_UNIT, file="bc.config.md")
          !Ignore the file header
          read(BOUNDARY_CONDITIONS_FILE_UNIT, *)
          read(BOUNDARY_CONDITIONS_FILE_UNIT, *)
          read(BOUNDARY_CONDITIONS_FILE_UNIT, *)


          buf = "NOT_FIN"
          do while (buf /= "FIN")

            read(BOUNDARY_CONDITIONS_FILE_UNIT, '(A)', iostat = ios) buf
            if(ios /=0 )then
              print*, "  /!\"
              print*, " /_|_\ ERROR: Boundary conditon read unsuccesful "
              print*, " checkpoint : module_surface - subroutine_ find_wall"
              STOP
            end if

            if (len_trim(buf) /= 0 ) then

            select case ( buf )

              case ( "# imn" )
                call dmsg(1, 'surface', 'find_wall', '# imn')
                call update_flag(1)
              case ( "# imx" )
                call dmsg(1, 'surface', 'find_wall', '# imx')
                call update_flag(2)
              case ( "# jmn" )
                call dmsg(1, 'surface', 'find_wall', '# jmn')
                call update_flag(3)
              case ( "# jmx" )
                call dmsg(1, 'surface', 'find_wall', '# jmx')
                call update_flag(4)
              case ( "# kmn" )
                call dmsg(1, 'surface', 'find_wall', '# kmn')
                call update_flag(5)
              case ( "# kmx" )
                call dmsg(1, 'surface', 'find_wall', '# kmx')
                call update_flag(6)
              case DEFAULT
                print*, "------"
                print*, "|  0  |" ," NO wall found " 
                print*, "------"

            end select

            end if

          
          end do
          close(BOUNDARY_CONDITIONS_FILE_UNIT)
        end subroutine find_wall

 

        subroutine surface_points()

          implicit none
          integer :: OL
          integer :: i, j, k, ind
          integer :: im=1, ix=1, id=1
          integer :: jm=1, jx=1, jd=1
          integer :: km=1, kx=1, kd=1

          call dmsg(1, 'surface', 'surface_points')


          ind = 0

          do OL = 1,6

            if (NO_SLIP_flag(OL) == 1 )then
              select case (OL)
                case (1)
                  km = 1
                  jm = 1
                  im = 1
                  kx = kmx-1
                  jx = jmx-1
                  ix = 1
                case (2)
                  km = 1
                  jm = 1
                  im = imx
                  kx = kmx-1
                  jx = jmx-1
                  ix = imx
                case (3)
                  km = 1
                  jm = 1
                  im = 1
                  kx = kmx-1
                  jx = 1
                  ix = imx-1
                case (4)
                  km = 1
                  jm = jmx
                  im = 1
                  kx = kmx-1
                  jx = jmx
                  ix = imx-1
                case (5)
                  km = 1
                  jm = 1
                  im = 1
                  kx = 1
                  jx = jmx-1
                  ix = imx-1
                case (6)
                  km = kmx
                  jm = 1
                  im = 1
                  kx = kmx
                  jx = jmx-1
                  ix = imx-1
                case DEFAULT
                  call dmsg(5, "Surface", 'Surface_points', 'FATAL  ERROR: select case')
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
                  id = min(ix-im, 1)
                  jd = min(jx-jm, 1)
                  kd = min(kx-km, 1)
                  wall_x(ind) = 0.5 * 0.25 * (                        &
                                              grid_x(i, j, k )        &
                                            + grid_x(i+id, j, k)      &
                                            + grid_x(i, j+jd, k)      &
                                            + grid_x(i, j, k+kd)      &
                                            + grid_x(i, j+jd, k+kd)   &
                                            + grid_x(i+id, j, k+kd)   &
                                            + grid_x(i+id, j+jd, k)   &
                                            + grid_x(i+id, j+jd, k+kd)  &
                                            )


                  wall_y(ind) = 0.5 * 0.25 * (                        &
                                              grid_y(i, j, k )        &
                                            + grid_y(i+id, j, k)      &
                                            + grid_y(i, j+jd, k)      &
                                            + grid_y(i, j, k+kd)      &
                                            + grid_y(i, j+jd, k+kd)   &
                                            + grid_y(i+id, j, k+kd)   &
                                            + grid_y(i+id, j+jd, k)   &
                                            + grid_y(i+id, j+jd, k+kd)  &
                                             )


                  wall_z(ind) = 0.5 * 0.25 * (                        &
                                              grid_z(i, j, k )        &
                                            + grid_z(i+id, j, k)      &
                                            + grid_z(i, j+jd, k)      &
                                            + grid_z(i, j, k+kd)      &
                                            + grid_z(i, j+jd, k+kd)   &
                                            + grid_z(i+id, j, k+kd)   &
                                            + grid_z(i+id, j+jd, k)   &
                                            + grid_z(i+id, j+jd, k+kd)  &
                                             )

                end do
              end do
            end do

            end if

        end do

        end subroutine surface_points


          
end module surface
