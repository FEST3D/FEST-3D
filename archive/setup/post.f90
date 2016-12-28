program main
  implicit none

  integer, parameter                         :: N_blocks = 16
  integer, parameter                         :: imx = 49
  integer, parameter                         :: jmx = 13
  integer, parameter                         :: kmx = 5
  integer, parameter                         :: level = 10       !state file number

  integer, parameter                         :: I_blocks = 4
  integer, parameter                         :: j_blocks = 4
  integer, parameter                         :: k_blocks = 1

  integer, dimension(16,2)                   :: im
  integer, dimension(16,2)                   :: jm
  integer, dimension(16,2)                   :: km

  real*8, dimension(:,:,:,:), allocatable  :: field
  real*8, dimension(:,:,:)    , allocatable  :: grid_x
  real*8, dimension(:,:,:)    , allocatable  :: grid_y
  real*8, dimension(:,:,:)    , allocatable  :: grid_z

  integer                                    :: i, j, k
  integer                                    :: gia=0, gja=0, gka=0
  integer                                    :: via=0, vja=0, vka=0
  integer                                    :: dia=0, dja=0, dka=0
  integer                                    :: pia=0, pja=0, pka=0
  integer                                    :: gip=0, gjp=0, gkp=0
  integer                                    :: vip=0, vjp=0, vkp=0
  integer                                    :: dip=0, djp=0, dkp=0
  integer                                    :: pip=0, pjp=0, pkp=0
  integer                                    :: fu
  character(len=2)                           :: num1   ! both
  character(len=5)                           :: num2  ! used for reading file
  character(len=32)                          :: ofile, dummy
  integer, parameter, dimension(N_blocks)    :: infile = (/(i+10 , i=1,N_blocks)/)
  integer, parameter                         :: outfile = 100
  integer                                    :: ios   ! for file opening status         
  real*8, dimension(:,:,:), allocatable      :: volume
  real*8, dimension(imx-1,jmx-1,kmx-1,3) :: cell

  allocate(field(1:imx-1, 1:jmx-1, 1:kmx-1, 1:8))
  allocate(grid_x(1:imx, 1:jmx, 1:kmx))
  allocate(grid_y(1:imx, 1:jmx, 1:kmx))
  allocate(grid_z(1:imx, 1:jmx, 1:kmx))
  allocate(volume(1:imx-1, 1:jmx-1, 1:kmx -1))

  do fu = 1,N_blocks

    write(num1,'(i2.2)') fu-1
    write(num2,'(i4.4)') level
    write(ofile,'(a)')  trim(num2)//"/process_"//trim(num1)//".vtk"
    print*, "Requested file: ", ofile
    open(unit = infile(fu), file = ofile , iostat=ios, status='old', action='read')
      if(ios/=0)then
        print*, "not able to open file requested"
        print*, "Requested file: ", ofile
      end if


    read(infile(fu), *) ! Skip first line
    read(infile(fu), *) ! Skip comment
    read(infile(fu), *) ! Skip ASCII
    read(infile(fu), *) ! Skip DATASET
    read(infile(fu), *) ! Skip Extra line
  
    im(fu,1) = 1
    jm(fu,1) = 1
    km(fu,1) = 1

    read(infile(fu), *) dummy, im(fu,2), jm(fu,2), km(fu,2)

    write(*,*), im(fu,2), jm(fu,2), km(fu,2)
    read(infile(fu), *) ! Skip POINTS
        if(mod(fu-1,4)==0 .and. fu/=1)then
          gip = 0
          gjp = gjp + jm(fu,2)-1
        end if
    do k = km(fu,1), km(fu,2)
     do j = jm(fu,1), jm(fu,2)
      do i = im(fu,1), im(fu,2)

        gia = gip + i
        gja = gjp + j
        gka = k
        read(infile(fu), *) grid_x(gia,gja,gka), grid_y(gia,gja,gka), grid_z(gia,gja,gka)
      end do
     end do
    end do
    gip = gia-1
    read(infile(fu), *) ! Skip blank space
  
    ! Cell data
    read(infile(fu), *) ! Skip CELL_DATA
    read(infile(fu), *) ! Skip VECTORS Velocity
  
        if(mod(fu-1,4)==0 .and. fu/=1)then
          vip = 0
          vjp = vjp + jm(fu,2)-1
        end if
    do k = km(fu,1), km(fu,2)-1
     do j = jm(fu,1), jm(fu,2)-1
      do i = im(fu,1), im(fu,2)-1
        via = vip + i
        vja = vjp + j
        vka = k
        read(infile(fu), *) field(via, vja, vka, 4), field(via, vja, vka,5), field(via, vja, vka,6)
      end do
     end do
    end do
    vip = via
  
        if(mod(fu-1,4)==0 .and. fu/=1)then
          dip = 0
          djp = djp + jm(fu,2)-1
        end if
    read(infile(fu), *) ! Skip Blank line
    read(infile(fu), *) ! Skip SCALARS DENSITY
    read(infile(fu), *) ! Skip LOOKUP_TABLE
    do k = km(fu,1), km(fu,2)-1
     do j = jm(fu,1), jm(fu,2)-1
      do i = im(fu,1), im(fu,2)-1
        dia = dip + i
        dja = djp + j
        dka = k
        read(infile(fu), *) field(dia, dja, dka, 7)
      end do
     end do
    end do
    dip =dia
  
        if(mod(fu-1,4)==0 .and. fu/=1)then
          pip = 0
          pjp = pjp + jm(fu,2)-1
        end if
    read(infile(fu), *) ! Skip Blank line
    read(infile(fu), *) ! Skip SCALARS Pressure
    read(infile(fu), *) ! Skip LOOKUP_TABLE
    do k = km(fu,1), km(fu,2)-1
     do j = jm(fu,1), jm(fu,2)-1
      do i = im(fu,1), im(fu,2)-1
        pia = pip + i
        pja = pjp + j
        pka = k
        read(infile(fu), *) field(pia, pja, pka, 8)
      end do
     end do
    end do
    pip = pia
    close(infile(fu))
  end do 

  open(unit=outfile, file="output"//trim(num2)//".vtk")



  write(outfile, fmt='(a)') '# vtk DataFile Version 3.1'

  write(outfile, '(a)') 'cfd-iitm output'

  write(outfile, '(a)') 'ASCII'
  write(outfile, '(a)') 'DATASET STRUCTURED_GRID'
  write(outfile, *) 

  write(outfile, fmt='(a, i0, a, i0, a, i0)') &
      'DIMENSIONS ', imx, ' ', jmx, ' ', kmx
  write(outfile, fmt='(a, i0, a)') &
      'POINTS ', imx*jmx*kmx, ' DOUBLE'

  do k = 1, kmx
   do j = 1, jmx
    do i = 1, imx
      write(outfile, fmt='(f0.16, a, f0.16, a, f0.16)') &
          grid_x(i, j, k), ' ', grid_y(i, j, k), ' ', grid_z(i, j, k)
    end do
   end do
  end do
  write(outfile, *) 

  ! Cell data
  write(outfile, fmt='(a, i0)') &
      'CELL_DATA ', (imx-1)*(jmx-1)*(kmx-1)

  ! Writing Velocity
  write(outfile, '(a)') 'VECTORS Velocity FLOAT'
  do k = 1, kmx - 1
   do j = 1, jmx - 1
    do i = 1, imx - 1
      write(outfile, fmt='(f0.16, a, f0.16, a, f0.16)') &
          field(i, j, k, 4), ' ', field(i, j, k, 5), ' ', field(i, j, k, 6)
    end do
   end do
  end do
  write(outfile, *) 

  ! Writing Density
  write(outfile, '(a)') 'SCALARS Density FLOAT'
  write(outfile, '(a)') 'LOOKUP_TABLE default'
  do k = 1, kmx - 1
   do j = 1, jmx - 1
    do i = 1, imx - 1
      write(outfile, fmt='(f0.16)') field(i, j, k, 7)
    end do
   end do
  end do
  write(outfile, *) 

  ! Writing Pressure
  write(outfile, '(a)') 'SCALARS Pressure FLOAT'
  write(outfile, '(a)') 'LOOKUP_TABLE default'
  do k = 1, kmx - 1
   do j = 1, jmx - 1
    do i = 1, imx - 1
      write(outfile, fmt='(f0.16)') field(i, j, k, 8)
    end do
   end do
  end do
  write(outfile, *)
  close(outfile)

  call compute_volumes()
  call fetch_entropy_l2_norm()
  call write_tecplot()
  deallocate(volume)
  deallocate(grid_x)
  deallocate(grid_y)
  deallocate(grid_z)
  deallocate(field)

  contains

    subroutine fetch_entropy_l2_norm()
      implicit none

      integer :: i,j,k
      integer :: outfile
      real*8  :: entropy_error
      real*8  :: entropy_inf
      real*8  :: entropy
      real*8  :: R_gas = 287
      real*8  :: gm = 1.4
      real*8  :: pressure_inf = 101325
      real*8  :: density_inf   = 1.225

      entropy       = 0.
      entropy_error = 0.
      entropy_inf   =  (R_gas/(gm - 1.))*log(pressure_inf/(density_inf**gm))

      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1
            entropy =  (R_gas/(gm - 1.))*log(field(i,j,k,8)/(field(i,j,k,7)**gm))
            entropy_error = entropy_error + (( entropy_inf - entropy )**2)*volume(i,j,k)
          end do
        end do
      end do

      entropy_error = sqrt( entropy_error/sum(volume) )

!      if()
      open(outfile, file = "entropy.txt")
!      open(667,file='entropy_temp', status="old", position="append", action="write")
      write(outfile,'(f0.16)') entropy_error
      close(outfile)

    end subroutine fetch_entropy_l2_norm

   
        function vol_tetrahedron(p1, p2, p3, p4)
            !-----------------------------------------------------------
            ! Compute the volume of a tetrahedron, given 4 points which
            ! are 1-D arrays
            ! Since we know that the determinant is to be evaluated of 
            ! a 3x3 matrix, we write the expression itself
            !-----------------------------------------------------------

            implicit none
            real, dimension(:), intent(in):: p1, p2, p3, p4
            real, dimension(1:3,1:3) :: A
            real :: vol_tetrahedron

            A(:, 1) = p1 - p4
            A(:, 2) = p2 - p4
            A(:, 3) = p3 - p4

            vol_tetrahedron = A(1,1) * (A(2,2)*A(3,3) - A(2,3)*A(3,2)) + &
                              A(1,2) * (A(2,3)*A(3,1) - A(2,1)*A(3,3)) + &
                              A(1,3) * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
            vol_tetrahedron = vol_tetrahedron / 6.                  
        
        end function vol_tetrahedron


        function vol_hexahedron(p_list)
            real, dimension(1:3, 1:8), intent(in) :: p_list
            real :: vol_hexahedron
            
            vol_hexahedron = 0.
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,1), p_list(:,5), &
                                             p_list(:,8), p_list(:,6))
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,7), p_list(:,8), &
                                             p_list(:,6), p_list(:,3))
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,8), p_list(:,4), &
                                             p_list(:,1), p_list(:,3))
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,6), p_list(:,1), &
                                             p_list(:,3), p_list(:,8))
            vol_hexahedron = vol_hexahedron + &
                             vol_tetrahedron(p_list(:,1), p_list(:,2), &
                                             p_list(:,6), p_list(:,3))
            
        end function vol_hexahedron
        
        subroutine compute_volumes()
            !-----------------------------------------------------------
            ! Compute the grid cell volumes
            ! Each grid is a hexahedron, whose volume is calculated by
            ! splitting it into 5 tetrahedrons, whose volume is known
            !-----------------------------------------------------------

            implicit none
            integer :: i,j,k
            real, dimension(1:3, 1:8) :: p_list

            do k = 1, kmx - 1
                do j = 1, jmx - 1
                    do i = 1, imx -1
                        p_list(:, :) = 0.
                        p_list(:, 1) = (/ grid_x(i,j,k), grid_y(i,j,k), grid_z(i,j,k) /)
                        p_list(:, 2) = (/ grid_x(i+1,j,k), grid_y(i+1,j,k), grid_z(i+1,j,k) /)
                        p_list(:, 3) = (/ grid_x(i+1,j+1,k), grid_y(i+1,j+1,k), grid_z(i+1,j+1,k) /)
                        p_list(:, 4) = (/ grid_x(i,j+1,k), grid_y(i,j+1,k), grid_z(i,j+1,k) /)
                        p_list(:, 5) = (/ grid_x(i,j,k+1), grid_y(i,j,k+1), grid_z(i,j,k+1) /)
                        p_list(:, 6) = (/ grid_x(i+1,j,k+1), grid_y(i+1,j,k+1), grid_z(i+1,j,k+1) /)
                        p_list(:, 7) = (/ grid_x(i+1,j+1,k+1), grid_y(i+1,j+1,k+1), grid_z(i+1,j+1,k+1) /)
                        p_list(:, 8) = (/ grid_x(i,j+1,k+1), grid_y(i,j+1,k+1), grid_z(i,j+1,k+1) /)
                        volume(i, j, k) = vol_hexahedron(p_list)
                    end do
                end do
            end do
            
        end subroutine compute_volumes

      subroutine write_tecplot()
        implicit none
        integer :: i,j,k

        call find_cell_centre()

        open(outfile, file="output"//trim(num2)//".dat")
        
        write(outfile,'(a)') "variables = X Y Z U V W Rho P"
        write(outfile, '(a, i0, a, i0, a, i0)') 'zone T="'//trim(num2)//'" I=',imx-1, ' J=',jmx-1, ' K=',kmx-1
        do k = 1,kmx-1
          do j = 1,jmx-1
            do i = 1,imx-1
              write(outfile, '(8(f0.16,2x))') cell(i,j,k,1:3), field(i,j,k,4:8)
            end do
          end do
        end do
        close(outfile)

      end subroutine write_tecplot

      subroutine find_cell_centre()
        implicit none
        integer :: i,j,k

        do k = 1,kmx-1
          do j = 1,jmx-1
            do i = 1,imx-1
              cell(i,j,k,1) = 0.125 * ( &
                                       grid_x(i,j,k)        &
                                      +grid_x(i+1,j,k)      &
                                      +grid_x(i,j+1,k)      &
                                      +grid_x(i,j,k+1)      &
                                      +grid_x(i+1,j+1,k)    &
                                      +grid_x(i,j+1,k+1)    &
                                      +grid_x(i+1,j,k+1)    &
                                      +grid_x(i+1,j+1,k+1)  &
                                    )

              cell(i,j,k,2) = 0.125 * ( &
                                       grid_y(i,j,k)        &
                                      +grid_y(i+1,j,k)      &
                                      +grid_y(i,j+1,k)      &
                                      +grid_y(i,j,k+1)      &
                                      +grid_y(i+1,j+1,k)    &
                                      +grid_y(i,j+1,k+1)    &
                                      +grid_y(i+1,j,k+1)    &
                                      +grid_y(i+1,j+1,k+1)  &
                                    )
              cell(i,j,k,3) = 0.125 * ( &
                                       grid_z(i,j,k)        &
                                      +grid_z(i+1,j,k)      &
                                      +grid_z(i,j+1,k)      &
                                      +grid_z(i,j,k+1)      &
                                      +grid_z(i+1,j+1,k)    &
                                      +grid_z(i,j+1,k+1)    &
                                      +grid_z(i+1,j,k+1)    &
                                      +grid_z(i+1,j+1,k+1)  &
                                    )
                end do
              end do
            end do
          end subroutine find_cell_centre



end program main
