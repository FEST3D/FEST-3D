integer function test_muscl() result(r)
   use vartypes
   use grid,  only : setup_grid
   use muscl
    
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
   scheme%interpolant = "muscl"
   scheme%ilimiter_switch = 0
   scheme%jlimiter_switch = 0
   scheme%klimiter_switch = 0
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
   call compute_muscl_states(qp, x_qp_left, x_qp_right, y_qp_left, y_qp_right, z_qp_left, z_qp_right, pdif, scheme, flow, dims)
   if( dims%imx==5 .and.dims%jmx==2 .and. dims%kmx==2)then
     if( x_qp_left(1,1,1,1)==0.5 .and. x_qp_right(1,1,1,1)==0.5)then
     if( x_qp_left(2,1,1,1)==1.5 .and. x_qp_right(2,1,1,1)==1.5)then
     if( x_qp_left(3,1,1,1)==2.5 .and. x_qp_right(3,1,1,1)>2.3 .and. x_qp_right(3,1,1,1)<2.4)then
     if( x_qp_left(4,1,1,1)>3.82  .and. x_qp_left(4,1,1,1)<3.85 .and. x_qp_right(4,1,1,1)==4.0)then
     if( x_qp_left(5,1,1,1)==6.0 .and. x_qp_right(5,1,1,1)>6.32 .and. x_qp_right(5,1,1,1)<6.35)then
       r = 0
     end if
     end if
     end if
     end if
     end if
   end if

end function
