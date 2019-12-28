integer function test_ausmp() result(r)
   use vartypes
   use ausmP
    
   type(flowtype) :: flow
   type(boundarytype) :: bc
   type(extent) :: dims
   !imx,jmx,kmx, = 2,2,2
   !n_var = 5
   real(wp), dimension(0:3,1:1,1:1,5) :: x_qp_left, x_qp_right
   real(wp), dimension(1:2,0:3,1:1,5) :: y_qp_left, y_qp_right
   real(wp), dimension(1:2,1:1,0:3,5) :: z_qp_left, z_qp_right
   type(facetype), dimension(-2:5,-2:4,-2:4):: Ifaces
   type(facetype), dimension(-2:4,-2:5,-2:4):: Jfaces
   type(facetype), dimension(-2:4,-2:4,-2:5):: Kfaces
   real(wp), dimension(1:2,1:1,1:1,5) :: F
   real(wp), dimension(1:1,1:2,1:1,5) :: G
   real(wp), dimension(1:1,1:1,1:2,5) :: H

   r = 1
   
   dims%imx  = 2
   dims%jmx  = 2
   dims%kmx  = 2
   dims%n_var = 5
   allocate(bc%make_F_flux_zero(1:dims%imx))
   allocate(bc%make_G_flux_zero(1:dims%jmx))
   allocate(bc%make_H_flux_zero(1:dims%kmx))
   bc%make_F_flux_zero=1
   bc%make_G_flux_zero=1
   bc%make_H_flux_zero=1
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
   x_qp_left(1,1,1,:)  = (/1.0,2.0,3.0,4.0,5.0/)
   x_qp_left(2,1,1,:)  = (/1.0,0.0,0.0,0.0,1.0/)
   x_qp_left(2,1,1,:)  = (/1.0,2.0,0.0,0.0,1.0/)
   x_qp_right(1,1,1,:) = (/1.0,2.0,3.0,4.0,5.0/)
   x_qp_right(2,1,1,:) = (/0.125,0.0,0.0,0.0,0.1/)
   y_qp_left  = 1.0
   y_qp_right = 1.0
   z_qp_left  = 1.0
   z_qp_right = 1.0
   call compute_fluxes(F,G,H,x_qp_left,x_qp_right,y_qp_left,y_qp_right,z_qp_left,z_qp_right,Ifaces,Jfaces,Kfaces,flow,bc,dims)
   if(F(1,1,1,2)>8.9 .and. F(1,1,1,2)<9.1)then
   if(F(2,1,1,2)>4.5 .and. F(2,1,1,2)<4.6)then
       r = 0
   end if
   end if

end function
