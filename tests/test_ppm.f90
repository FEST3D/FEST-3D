integer function test_ppm() result(r)
   use vartypes
   use ppm
    
   type(schemetype)  :: scheme
   type(flowtype) :: flow
   type(extent) :: dims
   !imx,jmx,kmx, = 5,2,2
   !n_var = 5
   real(wp), dimension(-2:7,-2:4,-2:4,1:5) :: qp=0.0
   real(wp), dimension(0:6,1:1,1:1,5) :: x_qp_left, x_qp_right
   real(wp), dimension(1:4,0:3,1:1,5) :: y_qp_left, y_qp_right
   real(wp), dimension(1:4,1:1,0:3,5) :: z_qp_left, z_qp_right
   real(wp), dimension(0:5,0:2,0:2) :: pdif

   r = 1
   
   dims%imx  = 5
   dims%jmx  = 2
   dims%kmx  = 2
   dims%n_var = 5
   scheme%interpolant = "ppm"
   scheme%ilimiter_switch = 1
   scheme%jlimiter_switch = 1
   scheme%klimiter_switch = 1
   qp(-2,1,1,:) = -2.0
   qp(-1,1,1,:) = -1.0
   qp( 0,1,1,:) =  0.0
   qp( 1,1,1,:) =  1.0
   qp( 2,1,1,:) =  2.0
   qp( 3,1,1,:) =  3.0
   qp( 4,1,1,:) =  5.0
   qp( 5,1,1,:) =  7.0
   qp( 6,1,1,:) =  7.0
   qp( 7,1,1,:) =  7.0
   call compute_ppm_states(qp, x_qp_left, x_qp_right, y_qp_left, y_qp_right, z_qp_left, z_qp_right, pdif, scheme, flow, dims)
   if( dims%imx==5 .and.dims%jmx==2 .and. dims%kmx==2)then
     if( x_qp_left(1,1,1,1)>0.48 .and. x_qp_left(1,1,1,1)<0.52 .and. x_qp_right(1,1,1,1)>0.48 .and. x_qp_right(1,1,1,1)<0.52)then
     if( x_qp_left(2,1,1,1)>1.48 .and. x_qp_left(2,1,1,1)<1.52 .and. x_qp_right(2,1,1,1)>1.48 .and. x_qp_right(2,1,1,1)<1.52)then
     if( x_qp_left(3,1,1,1)>2.3 .and. x_qp_left(3,1,1,1)<2.5 .and. x_qp_right(3,1,1,1)>2.3 .and. x_qp_right(3,1,1,1)<2.5)then
     if( x_qp_left(4,1,1,1)>3.8 .and. x_qp_left(4,1,1,1)<4.0 .and. x_qp_right(4,1,1,1)>3.8 .and. x_qp_right(4,1,1,1)<4.0)then
     if( x_qp_left(5,1,1,1)>6.0 .and. x_qp_left(5,1,1,1)<6.2 .and. x_qp_right(5,1,1,1)>6.5 .and. x_qp_right(5,1,1,1)<6.7)then
       r = 0
     end if
     end if
     end if
     end if
     end if
   end if

end function
