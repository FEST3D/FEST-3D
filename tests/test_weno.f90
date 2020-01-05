integer function test_weno() result(r)
   use vartypes
   use weno
    
   type(extent) :: dims
   !imx,jmx,kmx, = 5,2,2
   !n_var = 5
   real(wp), dimension(-2:7,-2:4,-2:4,1:5) :: qp=0.0
   real(wp), dimension(0:6,1:1,1:1,5) :: x_qp_left, x_qp_right
   real(wp), dimension(1:4,0:3,1:1,5) :: y_qp_left, y_qp_right
   real(wp), dimension(1:4,1:1,0:3,5) :: z_qp_left, z_qp_right

   r = 1
   
   dims%imx  = 5
   dims%jmx  = 2
   dims%kmx  = 2
   dims%n_var = 5
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
   call compute_weno_states(qp, x_qp_left, x_qp_right, y_qp_left, y_qp_right, z_qp_left, z_qp_right, dims)
   if( dims%imx==5 .and.dims%jmx==2 .and. dims%kmx==2)then
     if( x_qp_left(1,1,1,1)>=0.48 .and. x_qp_left(1,1,1,1)<0.52 .and. x_qp_right(1,1,1,1)>0.48 .and. x_qp_right(1,1,1,1)<0.52)then
     if( x_qp_left(2,1,1,1)>=1.48 .and. x_qp_left(2,1,1,1)<1.52 .and. x_qp_right(2,1,1,1)>1.48 .and. x_qp_right(2,1,1,1)<1.52)then
     if( x_qp_left(3,1,1,1)>=2.4 .and. x_qp_left(3,1,1,1)<2.5 .and. x_qp_right(3,1,1,1)>2.4 .and. x_qp_right(3,1,1,1)<2.5)then
     if( x_qp_left(4,1,1,1)>=3.6 .and. x_qp_left(4,1,1,1)<3.7 .and. x_qp_right(4,1,1,1)>3.8 .and. x_qp_right(4,1,1,1)<4.1)then
     if( x_qp_left(5,1,1,1)>=5.9 .and. x_qp_left(5,1,1,1)<6.1 .and. x_qp_right(5,1,1,1)>6.9 .and. x_qp_right(5,1,1,1)<7.1)then
       r = 0
     end if
     end if
     end if
     end if
     end if
   end if

end function
