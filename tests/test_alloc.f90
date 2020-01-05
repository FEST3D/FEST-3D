integer function test_alloc() result(r)
  use iso_fortran_env, only : wp => real64
  use utils
    
   real(wp), dimension(:,:,:), allocatable :: points

   r = 1
   call alloc(points,-2,11,-2,11,-2,11)
   points=0.0
   r = 0

end function
