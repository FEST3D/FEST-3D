integer function test_convergence() result(r)
   use vartypes
   use convergence
   use resnorm    , only: Res_abs
    
   type(controltype) :: control
   control%tolerance_type="Viscous_abs"
   control%tolerance = 1.e-8
   control%current_iter = 11

   allocate(Res_abs(0:5))
   Res_abs(1:5) = 1e-9
   r = 1
   if(converged(control))then
    r = 0
   end if

end function
