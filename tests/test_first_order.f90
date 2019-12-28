integer function test_first_order() result(r)
   use vartypes
   use face_interpolant
    
   type(schemetype)  :: scheme
   type(flowtype) :: flow
   type(extent) :: dims
   !imx,jmx,kmx, = 5,2,2
   !n_var = 5
   real(wp), dimension(-2:7,-2:4,-2:4,1:5) :: qp=0.0
   type(celltype), dimension(-2:7,-2:4,-2:4)  :: cells

   r = 1
   
   dims%imx   = 5
   dims%jmx   = 2
   dims%kmx   = 2
   dims%n_var = 5
   scheme%interpolant = "none"
   scheme%ilimiter_switch = 0
   scheme%jlimiter_switch = 0
   scheme%klimiter_switch = 0
   call setup_interpolant_scheme(dims)
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
   call compute_face_interpolant(qp, cells, scheme, flow, dims)
   if( dims%imx==5 .and.dims%jmx==2 .and. dims%kmx==2)then
     if( x_qp_left(1,1,1,1)==0 .and. x_qp_right(1,1,1,1)==1)then
     if( x_qp_left(2,1,1,1)==1 .and. x_qp_right(2,1,1,1)==2)then
     if( x_qp_left(3,1,1,1)==2 .and. x_qp_right(3,1,1,1)==3)then
     if( x_qp_left(4,1,1,1)==3 .and. x_qp_right(4,1,1,1)==5)then
     if( x_qp_left(5,1,1,1)==5 .and. x_qp_right(5,1,1,1)==7)then
       r = 0
     end if
     end if
     end if
     end if
     end if
   end if

end function
