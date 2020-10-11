integer function test_gradient() result(r)
   use vartypes
   use gradients
    
   type(controltype) :: control
   type(schemetype) :: scheme
   type(boundarytype) :: bc
   type(flowtype) :: flow
   type(celltype), dimension(-2:4,-2:4,-2:4) :: cells
   type(facetype), dimension(-2:5,-2:4,-2:4) :: Ifaces
   type(facetype), dimension(-2:4,-2:5,-2:4) :: Jfaces
   type(facetype), dimension(-2:4,-2:4,-2:5) :: Kfaces
   real(wp), dimension(-2:4,-2:4,-2:4,5) :: qp           
   real(wp), dimension(-2:4,-2:4,-2:4) :: Temp
   type(extent) :: dims

   dims%imx = 2
   dims%jmx = 2
   dims%kmx = 2

   dims%n_var = 5
   allocate(bc%make_F_flux_zero(1:dims%imx))
   allocate(bc%make_G_flux_zero(1:dims%jmx))
   allocate(bc%make_H_flux_zero(1:dims%kmx))
   bc%make_F_flux_zero=1
   bc%make_G_flux_zero=1
   bc%make_H_flux_zero=1
   bc%imin_id = 3
   bc%imax_id = 3
   bc%jmin_id = 3
   bc%jmax_id = 3
   bc%kmin_id = 3
   bc%kmax_id = 3

   cells%volume = 2.0
   Ifaces%nx=1.0
   Ifaces%ny=0.0
   Ifaces%nz=0.0
   Ifaces%A=2.0
   Jfaces%nx=0.0
   Jfaces%ny=1.0
   Jfaces%nz=0.0
   Jfaces%A=1.0
   Kfaces%nx=0.0
   Kfaces%ny=0.0
   Kfaces%nz=1.0
   Kfaces%A=1.0

   flow%mu_ref=1.0
   qp = 1.0
   Temp = 1.0
   qp(0,:,:,:) = 2.0
   qp(1,:,:,:) = 4.0
   qp(2,:,:,:) = 8.0
   qp(3,:,:,:) = 16.0

   r = 1
   call setup_gradients(control, scheme, flow, dims)
   call evaluate_all_gradients(qp, Temp, cells, Ifaces, Jfaces, Kfaces, scheme, bc, dims)
   if(gradqp_x(0,1,1,1)==1.5 .and. gradqp_x(1,1,1,1)==3.0 .and. gradqp_x(2,1,1,1)==6.0)then
     r = 0
   end if

end function
