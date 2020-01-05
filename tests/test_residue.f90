integer function test_residue() result(r)
   use vartypes
   use scheme
    
   type(extent) :: dims
   !imx,jmx,kmx, = 2,2,2
   !n_var = 5
   real(wp), dimension(1:2,1:1,1:1,5) :: F
   real(wp), dimension(1:1,1:2,1:1,5) :: G
   real(wp), dimension(1:1,1:1,1:2,5) :: H
   real(wp), dimension(1:1,1:1,1:1,5) :: residue

   r = 1
   
   dims%imx  = 2
   dims%jmx  = 2
   dims%kmx  = 2
   dims%n_var = 5
   F(1,1,1,:) = (/1.0, 2.0, 3.0, 4.0, 5.0/)
   F(2,1,1,:) = (/2.0, 2.0, 9.0, 0.0, 5.1/)
   G(1,1,1,:) = 1
   G(1,2,1,:) = 1
   H(1,1,1,:) = 1
   H(1,1,2,:) = 1

   call compute_residue(residue, F, G, H, dims)
   if(residue(1,1,1,1)==1.0 .and. residue(1,1,1,2)==0.0)then
   if(residue(1,1,1,3)==6.0 .and. residue(1,1,1,5)>0.09 .and. residue(1,1,1,5)<1.1)then
       r = 0
   end if
   end if

end function
