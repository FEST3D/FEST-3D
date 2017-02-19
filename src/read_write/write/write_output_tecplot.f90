module write_output_tecplot
  !---------------------------------------------------------
  ! This module write state + other variable in output.dat file
  !  data written at nodes
  ! BUG:
  !    node temp and node mach require node velocity, pressure
  !    and node density to be switched on
  !---------------------------------------------------------
  use global     , only : OUT_FILE_UNIT
  use global     , only : OUTIN_FILE_UNIT
  use global     , only : LONG_BUFFER_LENGTH
  use global     , only : outin_file

  use global_vars, only : write_data_format
  use global_vars, only : write_file_format
  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx
  use global_vars, only : grid_x
  use global_vars, only : grid_y
  use global_vars, only : grid_z
  use global_vars, only : density 
  use global_vars, only : x_speed 
  use global_vars, only : y_speed 
  use global_vars, only : z_speed 
  use global_vars, only : pressure 
  use global_vars, only : tk 
  use global_vars, only : tw 
  use global_vars, only : density_inf
  use global_vars, only : x_speed_inf
  use global_vars, only : y_speed_inf
  use global_vars, only : z_speed_inf
  use global_vars, only : pressure_inf 
  use global_vars, only : gm
  use global_vars, only : R_gas
  use global_vars, only : dist
  use global_vars, only : mu
  use global_vars, only : mu_t

  use global_vars, only : turbulence

  use utils, only : alloc, dmsg, dealloc
  use string

!  use source   , only :sst_mu

  implicit none
  private
  integer :: i,j,k
  real    :: speed_inf
  character(len=8) :: file_format
  character(len=16) :: data_format
  character(len=128) ::write_format

  real                                  :: node_density 
  real                                  :: node_x_speed 
  real                                  :: node_y_speed 
  real                                  :: node_z_speed 
  real                                  :: node_pressure 
  real                                  :: node_tk 
  real                                  :: node_tw
  real                                  :: node_mach
  real                                  :: node_temp
  real                                  :: node_mu
  real                                  :: node_mu_t
  real                                  :: node_dist

  character(len=64), dimension(:), allocatable  :: data_list
  ! string1 and string2 are for storage when string0 is full
  character(len=LONG_BUFFER_LENGTH)   :: string0
  character(len=LONG_BUFFER_LENGTH)   :: string1
  character(len=LONG_BUFFER_LENGTH)   :: string2
  integer :: kk, kj, ki
  integer           :: counter

  public :: write_file

  contains

    subroutine write_file()

      call read_input_file()
      call write_header()
      call write_data
      deallocate(data_list)

    end subroutine write_file

    subroutine write_data()

      implicit none
      integer :: n

      do n = 1,counter

        select case (data_list(n))
          
          case ('Density')
            call patch_8_verticies(density)

          case ('U')
            call patch_8_verticies(x_speed)

          case ('V')
            call patch_8_verticies(y_speed)

          case ('W')
            call patch_8_verticies(z_speed)

          case ('Pressure')
            call patch_8_verticies(pressure)

          case ('TKE')
            call patch_8_verticies(tk)

          case ('Omega')
            call patch_8_verticies(tw)

          case DEFAULT
            call dmsg(1,'write_output_tecplot', 'write_data', &
              'not apply patch to'//data_list(n))

        end select

      end do


      do kk = 1,kmx
        do  kj = 1,jmx
          do ki = 1,imx
            call get_write_string()
            if (write_data_format == 'ASCII') then
              write(OUT_FILE_UNIT, '(A)') trim(string2)//trim(string1)//trim(string0)
            else if (write_data_format == 'BINARY') then
              write(OUT_FILE_UNIT) trim(string2)//trim(string1)//trim(string0)
            end if
          end do
        end do
      end do

    end subroutine write_data
    
    subroutine read_input_file()

      implicit none
      character(len=64) :: buf
!      integer           :: dump
      open(OUTIN_FILE_UNIT, file=outin_file, status='old', action='read')
      ! reading only counter first for dimension
      read(OUTIN_FILE_UNIT, *)
      read(OUTIN_FILE_UNIT, *)
      buf=" not } "
      counter = 0
      do while (.true.)
        read(OUTIN_FILE_UNIT, *) buf
        if (trim(buf)=='}') EXIT
        counter = counter + 1
      end do
      allocate(data_list(1:counter))
      ! reading data types only, dumping counter
      rewind(OUTIN_FILE_UNIT)
      read(OUTIN_FILE_UNIT, *)
      read(OUTIN_FILE_UNIT, *)
      buf=" "
      do i = 1,counter
        read(OUTIN_FILE_UNIT, *) buf
        read(buf,*) data_list(i)
      end do
      close(OUTIN_FILE_UNIT)

    end subroutine read_input_file

    subroutine write_header()
      implicit none
      integer :: n
      character(len=LONG_BUFFER_LENGTH) :: header
      character(len=16) :: simx, sjmx, skmx

      write(simx, '(I0)') imx
      write(sjmx, '(I0)') jmx
      write(skmx, '(I0)') kmx


      if (write_data_format == 'ASCII') then

          write(header, '(A)') "VARIABLES = X Y Z "
        do n =  1, counter
          write(header, '(A)') trim(header)//" "//trim(data_list(n))
        end do
          write(OUT_FILE_UNIT, '(A)') trim(header)
          write(OUT_FILE_UNIT, '(A)') "zone i= "//trim(simx)//&
                                          " j= "//trim(sjmx)//&
                                          " k= "//trim(skmx)

      else if (write_data_format == 'BINARY') then

          write(header, '(A)') "VARIABLES = X Y Z "
        do n =  1, counter
          write(header, '(A)') trim(header)//" "//trim(data_list(n))
        end do
          write(OUT_FILE_UNIT) trim(header)
          write(OUT_FILE_UNIT) "zone i= "//trim(simx)//&
                                   " j= "//trim(sjmx)//&
                                   " k= "//trim(skmx)

      end if


    end subroutine write_header


    subroutine get_write_string()
      implicit none
      integer :: n

      string0 = ""
      string1 = ""
      string2 = ""

      write(string0, '(3(f0.16, 1x))') grid_x(ki,kj,kk), &
                                      grid_y(ki,kj,kk), &
                                      grid_z(ki,kj,kk)

      do n = 1,counter
        call check_string_limit()
        call add_to_string(data_list(n))
      end do

    end subroutine get_write_string

    subroutine check_string_limit()
      implicit none
      if (len(trim(string0))>200) then
        if (len(trim(string1))>200) then 
          if (len(trim(string2))>200) then
            call dmsg(5, 'write_output_tecplot', 'check_string_limit',&
                      'ERROR: No string found empty to write output data' )
            STOP
          else
            write(string2,'(A)')  trim(string0)
            string0  = ""
          end if
        else
          write(string1,'(A)') trim(string0)
          string0  = ""
        end if
      end if
    end subroutine check_string_limit


    subroutine add_to_string(arg)

      implicit none
      character(len=*) , intent(in)   :: arg
      character(len=32)               :: to_add
      character(len=*), parameter     :: add_format='(f0.16, 1x)'

      select case (arg)
        
        case ('Density')
          call calculate_node_data(density, node_density)
          write(to_add, add_format) node_density

        case ('U')
          call calculate_node_data(x_speed, node_x_speed)
          write(to_add, add_format) node_x_speed

        case ('V')
          call calculate_node_data(y_speed, node_y_speed)
          write(to_add, add_format) node_y_speed

        case ('W')
          call calculate_node_data(z_speed, node_z_speed)
          write(to_add, add_format) node_z_speed

        case ('Pressure')
          call calculate_node_data(pressure, node_pressure)
          write(to_add, add_format) node_pressure

        case ('Mach')
          call calculate_from_node_data('Mach')
          write(to_add, add_format) node_mach

        case ('Temperature')
          call calculate_from_node_data('Temperature')
          write(to_add, add_format) node_temp

        case ('TKE')
          call calculate_node_data(tk, node_tk)
          write(to_add, add_format) node_tk

        case ('Omega')
          call calculate_node_data(tw, node_tw)
          write(to_add, add_format) node_tw

        case ('Mu')
          call calculate_node_data(mu, node_mu)
          write(to_add, add_format) node_mu

        case ('Mu_t')
          call calculate_node_data(mu_t, node_mu_t)
          write(to_add, add_format) node_mu_t

        case ('Wall_dist')
          call calculate_node_data(dist, node_dist)
          write(to_add, add_format) node_dist

      end select
      string0 = trim(string0)//" "//trim(to_add)

    end subroutine add_to_string


    subroutine patch_8_verticies(arg)
      implicit none
      real, dimension(-2:imx, -2:jmx, -2:kmx), intent(inout) :: arg
      !copy cell center to vertex ghost cell
      arg(0  ,  0,  0)    = arg(    1,    1,    1)
      arg(imx,  0,  0)    = arg(imx-1,    1,    1)
      arg(  0,jmx,  0)    = arg(    1,jmx-1,    1)
      arg(  0,  0,kmx)    = arg(    1,    1,kmx-1)
      arg(imx,jmx,  0)    = arg(imx-1,jmx-1,    1)
      arg(  0,jmx,kmx)    = arg(    1,jmx-1,kmx-1)
      arg(imx,  0,kmx)    = arg(imx-1,    1,kmx-1)
      arg(imx,jmx,kmx)    = arg(imx-1,jmx-1,kmx-1)

    end subroutine patch_8_verticies

    subroutine calculate_node_data(arg, nodearg)

      implicit none
      real, dimension(-2:imx,-2:jmx,-2:kmx), intent(in)    :: arg
      real,                   intent(inout) :: nodearg

      nodearg = 0.125*(                                &
                         arg(ki-1, kj-1, kk-1)&
                       + arg(ki  , kj-1, kk-1)&
                       + arg(ki-1, kj  , kk-1)&
                       + arg(ki-1, kj-1, kk  )&
                       + arg(ki  , kj  , kk-1)&
                       + arg(ki-1, kj  , kk  )&
                       + arg(ki  , kj-1, kk  )&
                       + arg(ki  , kj  , kk  )&
                      )


    end subroutine calculate_node_data

    subroutine calculate_from_node_data(arg)
      implicit none
      character(len=*), intent(in) :: arg
      real :: speed2
      real :: sound
      
      select case (arg)
        
        case ('Mach')
          speed2 = node_x_speed**2 + node_y_speed**2 + node_z_speed**2
          sound  = gm*node_pressure/node_density
          node_mach = sqrt(speed2/sound) 

        case ('Temperature')
          node_temp = node_pressure/(node_density*R_gas)

        case DEFAULT
          call dmsg(5, 'write_output_tecplot', 'calculate_from_node_data',&
            'ERROR: calculate variables not recognised. -> ' + arg)

      end select
    end subroutine calculate_from_node_data

end module write_output_tecplot
