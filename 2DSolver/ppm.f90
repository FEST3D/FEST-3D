module ppm
    !-------------------------------------------------------------------
    ! PPM (Piecewise Parabolic Method) is an interpolation technique
    ! which can be used to create higher order extensions of schemes
    ! like the Van-Leer and LDFSS schemes.
    !-------------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx, grid_x, grid_y
    use state, only: qp, pressure, pressure_inf

    implicit none
    private

    real, dimension(:, :, :), allocatable, target :: x_qp_face_estimate
    real, dimension(:, :, :), allocatable, target :: y_qp_face_estimate
    real, dimension(:, :, :), allocatable, target :: x_qp_left, x_qp_right
    real, dimension(:, :, :), allocatable, target :: y_qp_left, y_qp_right
    real, dimension(:, :, :), pointer :: f_qp_left, f_qp_right
    real, dimension(:, :), allocatable :: pdif

    ! Public members
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_ppm_states
    public :: x_qp_left, x_qp_right
    public :: y_qp_left, y_qp_right

    contains

        subroutine setup_scheme()

            implicit none

            call dmsg(1, 'ppm', 'setup_ppm')

            call alloc(x_qp_face_estimate, 0, imx+1, 1, jmx-1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_face_estimate.')
            call alloc(y_qp_face_estimate, 1, imx-1, 0, jmx+1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_face_estimate.')
            call alloc(x_qp_left, 0, imx+1, 1, jmx-1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_left.')
            call alloc(x_qp_right, 0, imx+1, 1, jmx-1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_right.')
            call alloc(y_qp_left, 1, imx-1, 0, jmx+1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_left.')
            call alloc(y_qp_right, 1, imx-1, 0, jmx+1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_right.')
            call alloc(pdif, 0, imx, 0, jmx, &
                    errmsg='Error: Unable to allocate memory for' // &
                        'pdif')

        end subroutine setup_scheme

        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'ppm', 'destroy_ppm')

            call dealloc(x_qp_face_estimate)
            call dealloc(y_qp_face_estimate)
            call dealloc(x_qp_left)
            call dealloc(x_qp_right)
            call dealloc(y_qp_left)
            call dealloc(y_qp_right)
            call dealloc(pdif)

        end subroutine destroy_scheme

        subroutine compute_face_estimates(f_dir)

            implicit none
            character, intent(in) :: f_dir
            integer :: i, j 
            integer :: i_f, j_f ! Flags to determine face direction
            real, dimension(:, :, :), pointer :: f_qp_estimate

            call dmsg(1, 'ppm', 'compute_face_estimates')
            
            select case (f_dir)
                case ('x')
                    i_f = 1
                    j_f = 0
                    f_qp_estimate => x_qp_face_estimate
                case ('y')
                    i_f = 0
                    j_f = 1
                    f_qp_estimate => y_qp_face_estimate
                case default
                    call dmsg(5, 'ppm', 'pressure_based_switching', &
                            'Direction not recognised')
                    stop
            end select

            ! Interior faces
            do j = (1 + j_f), jmx - 1
             do i = (1 + i_f), imx - 1
                f_qp_estimate(i, j, :) = (7. * (qp(i, j, :) + &
                    qp(i - i_f, j - j_f, :)) - (qp(i + i_f, j + j_f, :) + &
                    qp(i - 2*i_f, j - 2*j_f, :))) / 12.
             end do
            end do

            ! Boundary of physical domain
            f_qp_estimate((1-i_f)+i_f*1 : (1-i_f)*(imx-1)+i_f*1, &
                          (1-j_f)+j_f*1 : (1-j_f)*(jmx-1)+j_f*1, :) = &
             (qp((1-i_f)+i_f*0 : (1-i_f)*(imx-1)+i_f*0, &
                 (1-j_f)+j_f*0 : (1-j_f)*(jmx-1)+j_f*0, :) + &
              qp((1-i_f)+i_f*1 : (1-i_f)*(imx-1)+i_f*1, &
                 (1-j_f)+j_f*1 : (1-j_f)*(jmx-1)+j_f*1, :)) /2. 
            
            f_qp_estimate((1-i_f)+i_f*imx : (1-i_f)*(imx-1)+i_f*imx, &
                          (1-j_f)+j_f*jmx : (1-j_f)*(jmx-1)+j_f*jmx, :) = &
             (qp((1-i_f)+i_f*imx : (1-i_f)*(imx-1)+i_f*imx, &
                 (1-j_f)+j_f*jmx : (1-j_f)*(jmx-1)+j_f*jmx, :) + &
              qp((1-i_f)+i_f*(imx-1) : (1-i_f)*(imx-1)+i_f*(imx-1), &
                 (1-j_f)+j_f*(jmx-1) : (1-j_f)*(jmx-1)+j_f*(jmx-1), :)) / 2.

            ! Boundary of ghost cells
            f_qp_estimate((1-i_f)+i_f*0 : (1-i_f)*(imx-1)+i_f*0, &
                          (1-j_f)+j_f*0 : (1-j_f)*(jmx-1)+j_f*0, :) = &
             (4.*qp((1-i_f)+i_f*0 : (1-i_f)*(imx-1)+i_f*0, &
                    (1-j_f)+j_f*0 : (1-j_f)*(jmx-1)+j_f*0, :) - &
                 qp((1-i_f)+i_f*1 : (1-i_f)*(imx-1)+i_f*1, &
                    (1-j_f)+j_f*1 : (1-j_f)*(jmx-1)+j_f*1, :)) /3. 
            
            f_qp_estimate((1-i_f)+i_f*(imx+1) : (1-i_f)*(imx-1)+i_f*(imx+1), &
                          (1-j_f)+j_f*(jmx+1) : (1-j_f)*(jmx-1)+j_f*(jmx+1), :) = &
             (4.*qp((1-i_f)+i_f*imx : (1-i_f)*(imx-1)+i_f*imx, &
                    (1-j_f)+j_f*jmx : (1-j_f)*(jmx-1)+j_f*jmx, :) - &
                 qp((1-i_f)+i_f*(imx-1) : (1-i_f)*(imx-1)+i_f*(imx-1), &
                    (1-j_f)+j_f*(jmx-1) : (1-j_f)*(jmx-1)+j_f*(jmx-1), :)) / 3.
        
        end subroutine compute_face_estimates

        subroutine remove_extrema(f_dir)

            implicit none
            character, intent(in) :: f_dir
            integer :: i, j, k
            integer :: i_f, j_f ! Flags to determine face direction
            real :: dqrl, dq6

            call dmsg(1, 'ppm', 'remove_extrema')
            
            select case (f_dir)
                case ('x')
                    i_f = 1
                    j_f = 0
                    f_qp_left => x_qp_left
                    f_qp_right => x_qp_right
                case ('y')
                    i_f = 0
                    j_f = 1
                    f_qp_left => y_qp_left
                    f_qp_right => y_qp_right
                case default
                    call dmsg(5, 'ppm', 'remove_extrema', &
                            'Direction not recognised')
                    stop
            end select
            
            ! Loop over cells (including ghost cells)
            do k = 1, 4
             do j = 1 - j_f, jmx - 1 + j_f
              do i = 1 - i_f, imx - 1 + i_f
                if ((f_qp_left(i+i_f, j+j_f, k) - qp(i, j, k)) * &
                    (qp(i, j, k) - f_qp_right(i, j, k)) <= 0) then
                    f_qp_left(i+i_f, j+j_f, k) = qp(i, j, k)
                    f_qp_right(i, j, k) = qp(i, j, k)
                else      
                    dqrl = f_qp_left(i+i_f, j+j_f, k) - f_qp_right(i, j, k)
                    dq6 = 6. * (qp(i,j,k) - 0.5*(f_qp_left(i+i_f, j+j_f, k) + &
                                                 f_qp_right(i, j, k)))
                    if (dqrl * dq6 > dqrl*dqrl) then
                        f_qp_right(i, j, k) = 3.*qp(i, j, k) - &
                                              2.*f_qp_left(i+i_f, j+j_f, k)
                    else if (-dqrl*dqrl > dqrl * dq6) then
                        f_qp_left(i+i_f, j+j_f, k) = 3.*qp(i, j, k) - &
                                              2.*f_qp_right(i, j, k)
                    end if
                end if
              end do
             end do 
            end do

        end subroutine remove_extrema

        subroutine pressure_based_switching(f_dir)

            implicit none
            ! Character can be x or y or z
            character, intent(in) :: f_dir
            integer :: i, j, i_end, j_end
            integer :: i_f, j_f  ! Flags to determine face direction
            real :: pd2

            call dmsg(1, 'ppm', 'pressure_based_switching')

            select case (f_dir)
                case ('x')
                    f_qp_left => x_qp_left
                    f_qp_right => x_qp_right
                    i_f = 1
                    j_f = 0
                    i_end = imx
                    j_end = jmx - 1
                case ('y')
                    f_qp_left => y_qp_left
                    f_qp_right => y_qp_right
                    i_f = 0
                    j_f = 1
                    i_end = imx - 1
                    j_end = jmx 
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
            do j = 1, jmx - 1
             do i = 1, imx - 1
                pd2 = abs(pressure(i + i_f*1, j + j_f*1) - &
                          pressure(i - i_f*1, j - j_f*1))
                pdif(i, j) = 1 - (pd2/(pd2 + pressure_inf))
             end do
            end do

            ! Update at ghost cells
            pdif((1-i_f):(1-i_f)*(imx-1), (1-j_f):(1-j_f)*(jmx-1)) = &
                pdif(1:imx-1 - i_f*(imx-2), 1:jmx-1 - j_f*(jmx-2))

            ! Loop over faces
            do j = 1, jmx - (1 - j_f)
             do i = 1, imx - (1 - i_f)
                f_qp_left(i, j, :) = qp(i - i_f*1, j - j_f*1, :) + (&
                    pdif(i - i_f*1, j - j_f*1) * ( &
                    f_qp_left(i, j, :) - qp(i - i_f*1, j - j_f*1, :)))

                f_qp_right(i, j, :) = qp(i, j, :) - (&
                    pdif(i, j) * ( &
                    qp(i, j, :) - f_qp_right(i, j, :)))
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
     
        subroutine compute_ppm_states()

            implicit none

            call compute_face_estimates('x')
            call init_left_and_right_xi_estimates()
            call remove_extrema('x')
            call pressure_based_switching('x')

            call compute_face_estimates('y')
            call init_left_and_right_eta_estimates()
            call remove_extrema('y')
            call pressure_based_switching('y')

        end subroutine compute_ppm_states

        function lin_interp(lambda, x1, x2)

            implicit none
            real, intent(in) :: lambda, x1, x2
            real :: lin_interp
            lin_interp = (1 - lambda)*x1 + (lambda*x2)

        end function lin_interp

end module ppm
