integer function test_viscosity() result(r)
   use vartypes
   use viscosity
   use wall_dist
   use gradients
    
   type(controltype) :: control
   type(schemetype) :: scheme
   type(boundarytype) :: bc
   type(flowtype) :: flow
   real(wp), dimension(-2:4,-2:4,-2:4,7) :: qp           
   type(extent) :: dims

   dims%imx = 2
   dims%jmx = 2
   dims%kmx = 2

   dims%n_var = 7
   allocate(bc%make_F_flux_zero(1:dims%imx))
   allocate(bc%make_G_flux_zero(1:dims%jmx))
   allocate(bc%make_H_flux_zero(1:dims%kmx))
   bc%make_F_flux_zero=1
   bc%make_G_flux_zero=1
   bc%make_H_flux_zero=1

   allocate(dist(-2:4,-2:4,-2:4))
   scheme%turbulence='sst'

   flow%mu_ref=1.0
   flow%mu_variation="sutherland_law"
   qp =1.0
   qp(1,1,1,:) = (/1.2,1.0,1.0,0.0,103320.0,1.0,1.0/)
   dist = 1.0

   r = 1
   call setup_gradients(control, scheme, flow, dims)
   gradqp_x=0.0
   gradqp_y=0.0
   gradqp_z=0.0
   call setup_viscosity(scheme, flow, dims)
   call calculate_viscosity(qp, scheme, flow, bc, dims)

   if(mu(1,1,1)>0.99 .and. mu(1,1,1)<1.01)then
    r = 0
   end if

end function
