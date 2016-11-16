module ppm
    !-------------------------------------------------------------------
    ! PPM (Piecewise Parabolic Method) is an interpolation technique
    ! which can be used to create higher order extensions of schemes
    ! like the Van-Leer and LDFSS schemes.
    !-------------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx, kmx, grid_x, grid_y, grid_z
    use state, only: qp, n_var, pressure, pressure_inf

    implicit none
    private

    real, dimension(:, :, :, :), allocatable, target :: x_qp_face_estimate
    real, dimension(:, :, :, :), allocatable, target :: y_qp_face_estimate
    real, dimension(:, :, :, :), allocatable, target :: z_qp_face_estimate
    real, dimension(:, :, :, :), allocatable, target :: x_qp_left, x_qp_right
    real, dimension(:, :, :, :), allocatable, target :: y_qp_left, y_qp_right
    real, dimension(:, :, :, :), allocatable, target :: z_qp_left, z_qp_right
    real, dimension(:, :, :, :), pointer :: f_qp_left, f_qp_right
    real, dimension(:, :, :), allocatable :: pdif
    integer :: iter_no

    ! Public members
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_ppm_states
!   public :: output_data
    public :: x_qp_left, x_qp_right
    public :: y_qp_left, y_qp_right
    public :: z_qp_left, z_qp_right

    contains

        subroutine setup_scheme()

            implicit none

            call dmsg(1, 'ppm', 'setup_ppm')
            iter_no = 0

            call alloc(x_qp_face_estimate, 0, imx+1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_face_estimate.')
            call alloc(y_qp_face_estimate, 1, imx-1, 0, jmx+1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_face_estimate.')
            call alloc(z_qp_face_estimate, 1, imx-1, 1, jmx-1, 0, kmx+1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'z_qp_face_estimate.')
            call alloc(x_qp_left, 0, imx+1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_left.')
            call alloc(x_qp_right, 0, imx+1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_right.')
            call alloc(y_qp_left, 1, imx-1, 0, jmx+1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_left.')
            call alloc(y_qp_right, 1, imx-1, 0, jmx+1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_right.')
            call alloc(z_qp_left, 1, imx-1, 1, jmx-1, 0, kmx+1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'z_qp_left.')
            call alloc(z_qp_right, 1, imx-1, 1, jmx-1, 0, kmx+1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'z_qp_right.')
            call alloc(pdif, 0, imx, 0, jmx, 0, kmx, &
                    errmsg='Error: Unable to allocate memory for' // &
                        'pdif')

        end subroutine setup_scheme

        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'ppm', 'destroy_ppm')

            call dealloc(x_qp_face_estimate)
            call dealloc(y_qp_face_estimate)
            call dealloc(z_qp_face_estimate)
            call dealloc(x_qp_left)
            call dealloc(x_qp_right)
            call dealloc(y_qp_left)
            call dealloc(y_qp_right)
            call dealloc(z_qp_left)
            call dealloc(z_qp_right)
            call dealloc(pdif)

        end subroutine destroy_scheme

        subroutine compute_face_estimates(f_dir)

            implicit none
            character, intent(in) :: f_dir
            integer :: i, j, k 
            integer :: i_f, j_f, k_f ! Flags to determine face direction
            real, dimension(:, :, :, :), pointer :: f_qp_estimate

            call dmsg(1, 'ppm', 'compute_face_estimates')
            
            select case (f_dir)
                case ('x')
                    i_f = 1
                    j_f = 0
                    k_f = 0
                    f_qp_estimate => x_qp_face_estimate
                case ('y')
                    i_f = 0
                    j_f = 1
                    k_f = 0
                    f_qp_estimate => y_qp_face_estimate
                case ('z')
                    i_f = 0
                    j_f = 0
                    k_f = 1
                    f_qp_estimate => z_qp_face_estimate
                case default
                    call dmsg(5, 'ppm', 'pressure_based_switching', &
                            'Direction not recognised')
                    stop
            end select

            !TODO: Vectorize this??
            ! Interior faces
            do k = (1 - k_f), kmx - 1 + 2*k_f
             do j = (1 - j_f), jmx - 1 + 2*j_f
              do i = (1 - i_f), imx - 1 + 2*i_f
                f_qp_estimate(i, j, k, :) = (7. * (qp(i, j, k, :) + &
                    qp(i - i_f, j - j_f, k - k_f, :)) - (qp(i + i_f, j + j_f, k + k_f, :) + &
                    qp(i - 2*i_f, j - 2*j_f, k - 2*k_f, :))) / 12.
              end do
             end do
            end do

        !   ! Boundary of physical domain
        !   f_qp_estimate((1-i_f)+i_f*1 : (1-i_f)*(imx-1)+i_f*1, &
        !                 (1-j_f)+j_f*1 : (1-j_f)*(jmx-1)+j_f*1, &
        !                 (1-k_f)+k_f*1 : (1-k_f)*(kmx-1)+k_f*1, :) = &
        !    (qp((1-i_f)+i_f*0 : (1-i_f)*(imx-1)+i_f*0, &
        !        (1-j_f)+j_f*0 : (1-j_f)*(jmx-1)+j_f*0, &
        !        (1-k_f)+k_f*0 : (1-k_f)*(kmx-1)+k_f*0, :) + &
        !     qp((1-i_f)+i_f*1 : (1-i_f)*(imx-1)+i_f*1, &
        !        (1-j_f)+j_f*1 : (1-j_f)*(jmx-1)+j_f*1, &
        !        (1-k_f)+k_f*1 : (1-k_f)*(kmx-1)+k_f*1, :)) /2. 
        !   
        !   f_qp_estimate((1-i_f)+i_f*imx : (1-i_f)*(imx-1)+i_f*imx, &
        !                 (1-j_f)+j_f*jmx : (1-j_f)*(jmx-1)+j_f*jmx, &
        !                 (1-k_f)+k_f*kmx : (1-k_f)*(kmx-1)+k_f*kmx, :) = &
        !    (qp((1-i_f)+i_f*imx : (1-i_f)*(imx-1)+i_f*imx, &
        !        (1-j_f)+j_f*jmx : (1-j_f)*(jmx-1)+j_f*jmx, &
        !        (1-k_f)+k_f*kmx : (1-k_f)*(kmx-1)+k_f*kmx, :) + &
        !     qp((1-i_f)+i_f*(imx-1) : (1-i_f)*(imx-1)+i_f*(imx-1), &
        !        (1-j_f)+j_f*(jmx-1) : (1-j_f)*(jmx-1)+j_f*(jmx-1), &
        !        (1-k_f)+k_f*(kmx-1) : (1-k_f)*(kmx-1)+k_f*(kmx-1), :)) / 2.

        !   ! Boundary of ghost cells
        !   f_qp_estimate((1-i_f)+i_f*0 : (1-i_f)*(imx-1)+i_f*0, &
        !                 (1-j_f)+j_f*0 : (1-j_f)*(jmx-1)+j_f*0, &
        !                 (1-k_f)+k_f*0 : (1-k_f)*(kmx-1)+k_f*0, :) = &
        !    (4.*qp((1-i_f)+i_f*0 : (1-i_f)*(imx-1)+i_f*0, &
        !           (1-j_f)+j_f*0 : (1-j_f)*(jmx-1)+j_f*0, &
        !           (1-k_f)+k_f*0 : (1-k_f)*(kmx-1)+k_f*0, :) - &
        !        qp((1-i_f)+i_f*1 : (1-i_f)*(imx-1)+i_f*1, &
        !           (1-j_f)+j_f*1 : (1-j_f)*(jmx-1)+j_f*1, &
        !           (1-k_f)+k_f*1 : (1-k_f)*(kmx-1)+k_f*1, :)) /3. 
        !   
        !   f_qp_estimate((1-i_f)+i_f*(imx+1) : (1-i_f)*(imx-1)+i_f*(imx+1), &
        !                 (1-j_f)+j_f*(jmx+1) : (1-j_f)*(jmx-1)+j_f*(jmx+1), &
        !                 (1-k_f)+k_f*(kmx+1) : (1-k_f)*(kmx-1)+k_f*(kmx+1), :) = &
        !    (4.*qp((1-i_f)+i_f*imx : (1-i_f)*(imx-1)+i_f*imx, &
        !           (1-j_f)+j_f*jmx : (1-j_f)*(jmx-1)+j_f*jmx, &
        !           (1-k_f)+k_f*kmx : (1-k_f)*(kmx-1)+k_f*kmx, :) - &
        !        qp((1-i_f)+i_f*(imx-1) : (1-i_f)*(imx-1)+i_f*(imx-1), &
        !           (1-j_f)+j_f*(jmx-1) : (1-j_f)*(jmx-1)+j_f*(jmx-1), &
        !           (1-k_f)+k_f*(kmx-1) : (1-k_f)*(kmx-1)+k_f*(kmx-1), :)) / 3.
        
        end subroutine compute_face_estimates

        subroutine remove_extrema(f_dir)

            implicit none
            character, intent(in) :: f_dir
            integer :: i, j, k, l
            integer :: i_f, j_f, k_f ! Flags to determine face direction
            real :: dqrl, dq6

            call dmsg(1, 'ppm', 'remove_extrema')
            
            select case (f_dir)
                case ('x')
                    i_f = 1
                    j_f = 0
                    k_f = 0
                    f_qp_left => x_qp_left
                    f_qp_right => x_qp_right
                case ('y')
                    i_f = 0
                    j_f = 1
                    k_f = 0
                    f_qp_left => y_qp_left
                    f_qp_right => y_qp_right
                case ('z')
                    i_f = 0
                    j_f = 0
                    k_f = 1
                    f_qp_left => z_qp_left
                    f_qp_right => z_qp_right
                case default
                    call dmsg(5, 'ppm', 'remove_extrema', &
                            'Direction not recognised')
                    stop
            end select
            
            !TODO: Vectorize this? Or will it be ugly?
            ! Loop over cells (including ghost cells)
            do l = 1, n_var            
             do k = 1 - k_f, kmx - 1 + k_f
              do j = 1 - j_f, jmx - 1 + j_f
               do i = 1 - i_f, imx - 1 + i_f
                if ((f_qp_left(i+i_f, j+j_f, k+k_f, l) - qp(i, j, k, l)) * &
                    (qp(i, j, k, l) - f_qp_right(i, j, k, l)) <= 0) then
                    f_qp_left(i+i_f, j+j_f, k+k_f, l) = qp(i, j, k, l)
                    f_qp_right(i, j, k, l) = qp(i, j, k, l)
                else      
                    dqrl = f_qp_left(i+i_f, j+j_f, k+k_f, l) - f_qp_right(i, j, k, l)
                    dq6 = 6. * (qp(i, j, k, l) - 0.5*(f_qp_left(i+i_f, j+j_f, k+k_f, l) + &
                                                      f_qp_right(i, j, k, l)))
                    if (dqrl * dq6 > dqrl*dqrl) then
                        f_qp_right(i, j, k, l) = 3.*qp(i, j, k, l) - &
                                                 2.*f_qp_left(i+i_f, j+j_f, k+k_f, l)
                    else if (-dqrl*dqrl > dqrl * dq6) then
                        f_qp_left(i+i_f, j+j_f, k+k_f, l) = 3.*qp(i, j, k, l) - &
                                                 2.*f_qp_right(i, j, k, l)
                    end if
                end if
               end do
              end do 
             end do
            end do

        end subroutine remove_extrema

        subroutine pressure_based_switching(f_dir)

            implicit none
            ! Character can be x or y or z
            character, intent(in) :: f_dir
            integer :: i, j, k, i_end, j_end, k_end
            integer :: i_f, j_f, k_f  ! Flags to determine face direction
            real :: pd2

            call dmsg(1, 'ppm', 'pressure_based_switching')

            select case (f_dir)
                case ('x')
                    f_qp_left => x_qp_left
                    f_qp_right => x_qp_right
                    i_f = 1
                    j_f = 0
                    k_f = 0
                    i_end = imx
                    j_end = jmx - 1
                    k_end = kmx - 1
                case ('y')
                    f_qp_left => y_qp_left
                    f_qp_right => y_qp_right
                    i_f = 0
                    j_f = 1
                    k_f = 0
                    i_end = imx - 1
                    j_end = jmx 
                    k_end = kmx - 1
                case ('z')
                    f_qp_left => z_qp_left
                    f_qp_right => z_qp_right
                    i_f = 0
                    j_f = 0
                    k_f = 1
                    i_end = imx - 1
                    j_end = jmx - 1 
                    k_end = kmx
                case default
                    call dmsg(5, 'ppm', 'pressure_based_switching', &
                            'Direction not recognised')
                    stop
            end select

            ! i_end and j_end denote number of faces
            ! Total number of cells including ghost_cells is
            ! (i_end+1) * j_end for xi faces and i_end*(j_end+1) for
            ! eta faces. 

            ! Loop over cells (physical)
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                pd2 = abs(pressure(i + i_f*1, j + j_f*1, k + k_f*1) - &
                          pressure(i - i_f*1, j - j_f*1, k - k_f*1))
                pdif(i, j, k) = 1 - (pd2/(pd2 + pressure_inf))
              end do
             end do
            end do

            ! Update at ghost cells
            pdif((1-i_f):(1-i_f)*(imx-1), (1-j_f):(1-j_f)*(jmx-1), &
                 (1-k_f):(1-k_f)*(kmx-1)) = &
                pdif(1:imx-1 - i_f*(imx-2), 1:jmx-1 - j_f*(jmx-2), &
                     1:kmx-1 - k_f*(kmx-2))

            ! Loop over faces
            do k = 1, kmx - (1 - k_f)            
             do j = 1, jmx - (1 - j_f)
              do i = 1, imx - (1 - i_f)
                f_qp_left(i, j, k, :) = qp(i - i_f*1, j - j_f*1, k - k_f*1, :) + (&
                    pdif(i - i_f*1, j - j_f*1, k - k_f*1) * ( &
                    f_qp_left(i, j, k, :) - qp(i - i_f*1, j - j_f*1, k - k_f*1, :)))

                f_qp_right(i, j, k, :) = qp(i, j, k, :) - (&
                    pdif(i, j, k) * ( &
                    qp(i, j, k, :) - f_qp_right(i, j, k, :)))
              end do
             end do
            end do

        end subroutine pressure_based_switching
        
        subroutine init_left_and_right_xi_estimates()

            implicit none

            ! x_qp_left and x_qp_right are stored at faces.
            x_qp_left = x_qp_face_estimate
            x_qp_right = x_qp_face_estimate

        end subroutine init_left_and_right_xi_estimates

        subroutine init_left_and_right_eta_estimates()
            
            implicit none

            ! y_qp_left and y_qp_right are stored at faces.
            y_qp_left = y_qp_face_estimate
            y_qp_right = y_qp_face_estimate

        end subroutine init_left_and_right_eta_estimates

        subroutine init_left_and_right_zeta_estimates()
            
            implicit none

            ! y_qp_left and y_qp_right are stored at faces.
            z_qp_left = z_qp_face_estimate
            z_qp_right = z_qp_face_estimate

        end subroutine init_left_and_right_zeta_estimates
     
        subroutine compute_ppm_states()

            implicit none

            iter_no = iter_no + 1

            call compute_face_estimates('x')
            call init_left_and_right_xi_estimates()
            if(ilimiter_switch==1)then
              call remove_extrema('x')
              call pressure_based_switching('x')
            end if

    !       if (iter_no == 5000) then 
    !           call write_line_data_x('qp-before.txt')
    !       end if
    !       if (iter_no == 5000) then 
    !           call write_line_data_x('qp-after.txt')
    !       end if

            call compute_face_estimates('y')
            call init_left_and_right_eta_estimates()
            if(ilimiter_switch==1)then
              call remove_extrema('y')
              call pressure_based_switching('y')
            end if

            call compute_face_estimates('z')
            call init_left_and_right_zeta_estimates()
            if(ilimiter_switch==1)then
              call remove_extrema('z')
              call pressure_based_switching('z')
            end if

        end subroutine compute_ppm_states

        function lin_interp(lambda, x1, x2)

            implicit none
            real, intent(in) :: lambda, x1, x2
            real :: lin_interp
            lin_interp = (1 - lambda)*x1 + (lambda*x2)

        end function lin_interp

 !      subroutine write_line_data_x(filename)

 !          implicit none

 !          integer :: i, j, k
 !          character(len=*), intent(in) :: filename
 !          real :: x1, x2

 !          ! Writing out a row of data
 !          j = 5
 !          
 !          ! Writing out Pressure
 !          k = 4

 !          open(71, file=filename)
 !          ! Writing it as cell based
 !          do i = 1, imx-1
 !              x1 = grid_x(i, j)
 !              x2 = grid_x(i+1, j)
 !              write(71, *) lin_interp(0.02, x1, x2), x_qp_right(i, j, k), &
 !                           lin_interp(0.5, x1, x2), qp(i, j, k), &
 !                           lin_interp(0.98, x1, x2), x_qp_left(i+1, j, k)
 !          end do
 !          close(71)

 !      end subroutine write_line_data_x

     !  subroutine write_line_data()

     !      implicit none

     !      integer :: i, j, k
     !      character(len=20) :: f1, f2, f3

     !      ! Writing out a row of data
     !      j = 1
     !      
     !      ! Writing out Pressure
     !      k = 4
     !      i = 100

     !      f1 = 'qp_left_after'
     !      f2 = 'qp_after'
     !      f3 = 'qp_right_after'

     !      open(71, file=f1)
     !      do i = 1, imx-1
     !          write(71, *) grid_x(i, j), x_qp_right(i, j, k)
     !      end do
     !      close(71)

     !      open(71, file=f2)
     !      do i = 1, imx-1
     !          write(71, *) 0.5 * (grid_x(i, j) + grid_x(i+1, j)), qp(i, j, k)
     !      end do
     !      close(71)

     !      open(71, file=f3)
     !      do i = 2, imx
     !          write(71, *) grid_x(i, j), x_qp_left(i, j, k)
     !      end do
     !      close(71)

     !  end subroutine write_line_data

!       subroutine write_line_data_y()

!           implicit none

!           integer :: i, j, k
!           character(len=20) :: f1, f2, f3

!           ! Writing out a row of data
!           i = 49
!           
!           ! Writing out Pressure
!           k = 4

!           f1 = 'y_qp_left_after'
!           f2 = 'y_qp_after'
!           f3 = 'y_qp_right_after'

!           open(71, file=f1)
!           do j = 1, jmx-1
!               write(71, *) grid_y(i, j), y_qp_right(i, j, k)
!           end do
!           close(71)

!           open(71, file=f2)
!           do j = 1, jmx-1
!               write(71, *) 0.5 * (grid_y(i, j) + grid_y(i, j+1)), qp(i, j, k)
!           end do
!           close(71)

!           open(71, file=f3)
!           do j = 2, jmx
!               write(71, *) grid_y(i, j), y_qp_left(i, j, k)
!           end do
!           close(71)

!       end subroutine write_line_data_y

!       subroutine output_data(filename)
!           implicit none
!           character(len=*), intent(in) :: filename
!           integer :: i, j
!           open(71, file=filename)
!           write(71, *) 'density'
!           write(71, *) 'xi_left'
!           do j = 1, jmx-1
!               do i = 0, imx+1
!                   write(71, *) x_qp_left(i, j, 1)
!               end do
!           end do
!           write(71, *) 'xi_right'
!           do j = 1, jmx-1
!               do i = 0, imx+1
!                   write(71, *) x_qp_right(i, j, 1)
!               end do
!           end do
!           write(71, *) 'eta_left'
!           do j = 0, jmx+1
!               do i = 1, imx-1
!                   write(71, *) y_qp_left(i, j, 1)
!               end do
!           end do
!           write(71, *) 'eta_right'
!           do j = 0, jmx+1
!               do i = 1, imx-1
!                   write(71, *) y_qp_right(i, j, 1)
!               end do
!           end do
!           write(71, *) 'x_speed'
!           write(71, *) 'xi_left'
!           do j = 1, jmx-1
!               do i = 0, imx+1
!                   write(71, *) x_qp_left(i, j, 2)
!               end do
!           end do
!           write(71, *) 'xi_right'
!           do j = 1, jmx-1
!               do i = 0, imx+1
!                   write(71, *) x_qp_right(i, j, 2)
!               end do
!           end do
!           write(71, *) 'eta_left'
!           do j = 0, jmx+1
!               do i = 1, imx-1
!                   write(71, *) y_qp_left(i, j, 2)
!               end do
!           end do
!           write(71, *) 'eta_right'
!           do j = 0, jmx+1
!               do i = 1, imx-1
!                   write(71, *) y_qp_right(i, j, 2)
!               end do
!           end do
!           write(71, *) 'y_speed'
!           write(71, *) 'xi_left'
!           do j = 1, jmx-1
!               do i = 0, imx+1
!                   write(71, *) x_qp_left(i, j, 3)
!               end do
!           end do
!           write(71, *) 'xi_right'
!           do j = 1, jmx-1
!               do i = 0, imx+1
!                   write(71, *) x_qp_right(i, j, 3)
!               end do
!           end do
!           write(71, *) 'eta_left'
!           do j = 0, jmx+1
!               do i = 1, imx-1
!                   write(71, *) y_qp_left(i, j, 3)
!               end do
!           end do
!           write(71, *) 'eta_right'
!           do j = 0, jmx+1
!               do i = 1, imx-1
!                   write(71, *) y_qp_right(i, j, 3)
!               end do
!           end do
!           write(71, *) 'pressure'
!           write(71, *) 'xi_left'
!           do j = 1, jmx-1
!               do i = 0, imx+1
!                   write(71, *) x_qp_left(i, j, 4)
!               end do
!           end do
!           write(71, *) 'xi_right'
!           do j = 1, jmx-1
!               do i = 0, imx+1
!                   write(71, *) x_qp_right(i, j, 4)
!               end do
!           end do
!           write(71, *) 'eta_left'
!           do j = 0, jmx+1
!               do i = 1, imx-1
!                   write(71, *) y_qp_left(i, j, 4)
!               end do
!           end do
!           write(71, *) 'eta_right'
!           do j = 0, jmx+1
!               do i = 1, imx-1
!                   write(71, *) y_qp_right(i, j, 4)
!               end do
!           end do
!           close(71)
!       end subroutine output_data

end module ppm
