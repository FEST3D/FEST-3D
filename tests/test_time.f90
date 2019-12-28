integer function test_time() result(r)
   use vartypes
   use time
   use face_interpolant
    
   type(controltype) :: control
   type(schemetype) :: scheme
   type(flowtype) :: flow
   type(extent) :: dims
   !imx,jmx,kmx, = 2,2,2
   !n_var = 5
   type(facetype), dimension(-2:5,-2:4,-2:4):: Ifaces
   type(facetype), dimension(-2:4,-2:5,-2:4):: Jfaces
   type(facetype), dimension(-2:4,-2:4,-2:5):: Kfaces
   real(wp), dimension(-2:4,-2:4,-2:4,5) :: qp
   type(celltype), dimension(-2:4,-2:4,-2:4) :: cells
   real(wp), dimension(:, :, :), allocatable     :: delta_t  
   !< Local time increment value at each cell center

   r = 1
   
   dims%imx  = 2
   dims%jmx  = 2
   dims%kmx  = 2
   dims%n_var = 5
   allocate(x_qp_left(0:3,1:1,1:1,5))
   allocate(y_qp_left(1:1,0:3,1:1,5))
   allocate(z_qp_left(1:1,1:1,0:3,5))
   allocate(x_qp_right(0:3,1:1,1:1,5))
   allocate(y_qp_right(1:1,0:3,1:1,5))
   allocate(z_qp_right(1:1,1:1,0:3,5))
   Ifaces(:,:,:)%nx = 1.0
   Ifaces(:,:,:)%ny = 0.0
   Ifaces(:,:,:)%nz = 0.0
   Ifaces(:,:,:)%A = 1.0
   Jfaces(:,:,:)%nx = 0.0
   Jfaces(:,:,:)%ny = 1.0
   Jfaces(:,:,:)%nz = 0.0
   Jfaces(:,:,:)%A = 1.0
   Kfaces(:,:,:)%nx = 0.0
   Kfaces(:,:,:)%ny = 0.0
   Kfaces(:,:,:)%nz = 1.0
   Kfaces(:,:,:)%A = 1.0
   x_qp_left  = 1.0
   x_qp_right = 1.0
   y_qp_left  = 1.0
   y_qp_right = 1.0
   z_qp_left  = 1.0
   z_qp_right = 1.0
   cells%volume = 1.0
   call setup_time(delta_t, control, dims)
   call compute_time_step(qp, delta_t, 1.0, cells, Ifaces, Jfaces, Kfaces, scheme, flow, dims)
   print*, delta_t(1,1,1)
   if(delta_t(1,1,1)>0.135 .and. delta_t(1,1,1)<0.145)then
       r = 0
   end if

end function
