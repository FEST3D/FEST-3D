integer function test_state() result(r)
   use vartypes
   use state
    
   type(filetype) :: files
   type(controltype) :: control
   type(schemetype) scheme
   type(flowtype) :: flow
   type(extent) :: dims
   real(wp), dimension(:,:,:,:), allocatable :: qp

   r = 1
   scheme%turbulence='sst'
   flow%mu_ref=1
   dims%imx = 5
   dims%jmx = 5
   dims%kmx = 5
   call setup_state(files, qp, control, scheme, flow, dims)
   if( dims%n_var==7 .and. flow%Turb_intensity_inf==0.01 .and. flow%Reynolds_number==120 .and. flow%tk_inf == 1.5)then
       r = 0
   end if

end function
