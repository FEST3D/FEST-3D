! this include file calcuolate time step 
  if(time_stepping_method .eq. 'g' .or. time_stepping_method .eq. 'l') then
    if (time_stepping_method .eq. 'g' .and. global_time_step > 0) then
      delta_t = global_time_step
    else
      do k = 1, kmx - 1
       do j = 1, jmx - 1
        do i = 1, imx - 1
         ! For orientation, refer to the report. The standard i,j,k 
         ! direction are marked. All orientation notations are w.r.t 
         ! to the perspective shown in the image.
  
         ! Faces with lower index
         delta_t(i,j,k) = 2*CFL*volume(i,j,k)/&
                (&! For left face: i.e., lower index face along xi direction
                (xA(i,j,k)*(abs( &
              (x_speed(i, j, k) * xnx(i, j, k)) + &
              (y_speed(i, j, k) * xny(i, j, k)) + &
              (z_speed(i, j, k) * xnz(i, j, k))) + &       !U           
              (0.5 * (sqrt(gm * x_qp_left(i, j, k, 5) / &   
                               x_qp_left(i, j, k, 1)) + &  !+
                     sqrt(gm * x_qp_right(i, j, k, 5) / &   
                               x_qp_right(i, j, k, 1)) ))))&  !a
           +&! For front face, i.e., lower index face along eta direction
           (yA(i,j,k)*(abs( &
              (x_speed(i, j, k) * ynx(i, j, k)) + &
              (y_speed(i, j, k) * yny(i, j, k)) + &
              (z_speed(i, j, k) * ynz(i, j, k))) + &
              (0.5 * (sqrt(gm * y_qp_left(i, j, k, 5) / &
                                y_qp_left(i, j, k, 1)) + &
                      sqrt(gm * y_qp_right(i, j, k, 5) / &
                                y_qp_right(i, j, k, 1)) ))))&
         +&! For bottom face, i.e., lower index face along zeta direction
         (zA(i,j,k)*(abs( &
              (x_speed(i, j, k) * znx(i, j, k)) + &
              (y_speed(i, j, k) * zny(i, j, k)) + &
              (z_speed(i, j, k) * znz(i, j, k))) + &
              (0.5 * (sqrt(gm * z_qp_left(i, j, k, 5) / &
                                z_qp_left(i, j, k, 1)) + &
                      sqrt(gm * z_qp_right(i, j, k, 5) / &
                                z_qp_right(i, j, k, 1)) ))))&
         +&! For right face, i.e., higher index face along xi direction
         (xA(i+1,j,k)*(abs( &
              (x_speed(i+1, j, k) * xnx(i+1, j, k)) + &
              (y_speed(i+1, j, k) * xny(i+1, j, k)) + &
              (z_speed(i+1, j, k) * xnz(i+1, j, k))) + &
          (0.5 * (sqrt(gm * x_qp_left(i+1,j,k,5) / x_qp_left(i+1,j,k,1)) + &
                 sqrt(gm * x_qp_right(i+1,j,k,5) / x_qp_right(i+1,j,k,1)) ))))&
         +&! For back face, i.e., higher index face along eta direction
         (yA(i,j+1,k)*(abs( &
              (x_speed(i, j+1, k) * ynx(i, j+1, k)) + &
              (y_speed(i, j+1, k) * yny(i, j+1, k)) + &
              (z_speed(i, j+1, k) * ynz(i, j+1, k))) + &
          (0.5 * (sqrt(gm * y_qp_left(i,j+1,k,5) / y_qp_left(i,j+1,k,1)) + &
                 sqrt(gm * y_qp_right(i,j+1,k,5) / y_qp_right(i,j+1,k,1)) ))))&
         +&! For top face, i.e., higher index face along zeta direction
         (zA(i,j,k+1)*(abs( &
              (x_speed(i, j, k+1) * znx(i, j, k+1)) + &
              (y_speed(i, j, k+1) * zny(i, j, k+1)) + &
              (z_speed(i, j, k+1) * znz(i, j, k+1))) + &
         ( 0.5 * (sqrt(gm * z_qp_left(i,j,k+1,5) / z_qp_left(i,j,k+1,1)) + &
                 sqrt(gm * z_qp_right(i,j,k+1,5) / z_qp_right(i,j,k+1,1)) ))))&
                 )
        end do
       end do
      end do
      if (time_stepping_method .eq. 'g' .and. global_time_step <= 0) delta_t=Minval(delta_t)
    end if
  else
      call dmsg(5, 'solver', 'compute_time_step', &
              msg='Value for time_stepping_method (' // &
                  time_stepping_method // ') not recognized.')
      stop
  end if

