    !< Higher order face state reconstruction method:PPM
module ppm
    !<
    !<Reference: Colella, P. and Woodward, P.R., The piecewise 
    !<parabolic method (PPM) for gas-dynamical simulations, Journal
    !<of computational physics, vol. 54, no. 1, pp.174-201, 1984
    !-------------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg

    use global_vars, only : imx
    use global_vars, only : jmx
    use global_vars, only : kmx
    use global_vars, only : grid_x
    use global_vars, only : grid_y
    use global_vars, only : grid_z

    use global_vars, only : qp
    use global_vars, only : n_var
    use global_vars, only : pressure
    use global_vars, only : pressure_inf
    use global_vars, only : ilimiter_switch
    use global_vars, only : jlimiter_switch
    use global_vars, only : klimiter_switch
    use global_vars, only : iPB_switch
    use global_vars, only : jPB_switch
    use global_vars, only : kPB_switch


    implicit none
    private

    real, dimension(:, :, :, :), allocatable, target :: x_qp_face_estimate
    !< Store the I face estimate from 4th order reconstruction
    real, dimension(:, :, :, :), allocatable, target :: y_qp_face_estimate
    !< Store the J face estimate from 4th order reconstruction
    real, dimension(:, :, :, :), allocatable, target :: z_qp_face_estimate
    !< Store the K face estimate from 4th order reconstruction
    real, dimension(:, :, :, :), allocatable, target :: x_qp_left
      !< Store primitive state at the I-face left side
    real, dimension(:, :, :, :), allocatable, target :: x_qp_right
      !< Store primitive state at the I-face right side
    real, dimension(:, :, :, :), allocatable, target :: y_qp_left
      !< Store primitive state at the J-face left side
    real, dimension(:, :, :, :), allocatable, target :: y_qp_right
      !< Store primitive state at the J-face right side
    real, dimension(:, :, :, :), allocatable, target :: z_qp_left
      !< Store primitive state at the K-face left side
    real, dimension(:, :, :, :), allocatable, target :: z_qp_right
      !< Store primitive state at the K-face right side
    real, dimension(:, :, :, :), pointer :: f_qp_left
    !< Generalized pointer for any I-J-K direction> f_qp_left can 
    !< either point to x_qp_left, y_qp_left or z_qp_left
    real, dimension(:, :, :, :), pointer :: f_qp_right
    !< Generalized pointer for any I-J-K direction> f_qp_right can 
    !< either point to x_qp_right, y_qp_right or z_qp_right
    real, dimension(:, :, :), allocatable :: pdif
    !< Used for pressure based witch

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
          !< Allocate memoery to all array which store state
          !< the face.

            implicit none

            call dmsg(1, 'ppm', 'setup_ppm')

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
          !< Deallocate all the array used 

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
          !< Subroutine to calculate state at the face, generalized for

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

        end subroutine compute_face_estimates

        subroutine remove_extrema(f_dir)
          !< Remove extrema from the state estimated. 
          !< Limiting the value in case of PPM

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
          !< Pressure based switching. 
          !< User x,y, or z for I,J,or K face respectively
          !----------------------------------------------

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
            pdif(((imx-1)*i_f)+1:imx-1+i_f, &
                 ((jmx-1)*j_f)+1:jmx-1+j_f, &
                 ((kmx-1)*k_f)+1:kmx-1+k_f) &
                                        =   &
                pdif(i_f*(imx-2)+1:imx-1, &
                     j_f*(jmx-2)+1:jmx-1, &
                     k_f*(kmx-2)+1:kmx-1)

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
          !< Initialize the left and right state at I direction face

            implicit none

            ! x_qp_left and x_qp_right are stored at faces.
            x_qp_left = x_qp_face_estimate
            x_qp_right = x_qp_face_estimate

        end subroutine init_left_and_right_xi_estimates

        subroutine init_left_and_right_eta_estimates()
          !< Initialize the left and right state at J direction face
            
            implicit none

            ! y_qp_left and y_qp_right are stored at faces.
            y_qp_left = y_qp_face_estimate
            y_qp_right = y_qp_face_estimate

        end subroutine init_left_and_right_eta_estimates

        subroutine init_left_and_right_zeta_estimates()
          !< Initialize the left and right state at K direction face
            
            implicit none

            ! y_qp_left and y_qp_right are stored at faces.
            z_qp_left = z_qp_face_estimate
            z_qp_right = z_qp_face_estimate

        end subroutine init_left_and_right_zeta_estimates
     
        subroutine compute_ppm_states()
          !< Call PPM face-state reconstruction for each face
          !< with optional call for remove extrema based on
          !< input limter switch and call pressure based switching
          !< based on input pressure based switch

            implicit none

            call compute_face_estimates('x')
            call init_left_and_right_xi_estimates()
            if(ilimiter_switch==1)then
              call remove_extrema('x')
            end if
            if (iPB_switch==1)then
              call pressure_based_switching('x')
            end if

            call compute_face_estimates('y')
            call init_left_and_right_eta_estimates()
            if(jlimiter_switch==1)then
              call remove_extrema('y')
            end if
            if (jPB_switch==1)then
              call pressure_based_switching('y')
            end if

            call compute_face_estimates('z')
            call init_left_and_right_zeta_estimates()
            if(klimiter_switch==1)then
              call remove_extrema('z')
            end if
            if (kPB_switch==1)then
              call pressure_based_switching('z')
            end if

        end subroutine compute_ppm_states


end module ppm
