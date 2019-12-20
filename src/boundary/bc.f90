   !< Setup boundary condition for the domain
module bc
   !< Setup boundary condition for the domain
  !-------------------------------------------
  use vartypes
  use global_vars, only: imin_id
  use global_vars, only: imax_id
  use global_vars, only: jmin_id
  use global_vars, only: jmax_id
  use global_vars, only: kmin_id
  use global_vars, only: kmax_id
  use global_vars, only: c1
  use global_vars, only: c2
  use global_vars, only: c3
  use global_vars, only: id
  use global_vars, only: face_names
  use global_vars, only: make_F_flux_zero
  use global_vars, only: make_G_flux_zero
  use global_vars, only: make_H_flux_zero
  use global_vars, only: PbcId
  use utils, only: alloc

  use read_bc   ,  only: read_fixed_values


  implicit none
  private

  !integer                        :: face_num
  !< Number of the face : 1:imin, 2:imax, 3:jmin, 4:jmax, 5:kmin, 6:kmax

  public :: setup_bc
!  public :: destroy_bc


  contains

    subroutine setup_bc(files, scheme, flow, dims)
      !< Initialization and allocate memory of boundary condition variables
      implicit none
      type(filetype), intent(in) :: files
      type(schemetype), intent(in) :: scheme
      type(flowtype), intent(in) :: flow
      type(extent), intent(in) :: dims
      !check for periodic bc
      if(PbcId(1)>=0) imin_id=-10
      if(PbcId(2)>=0) imax_id=-10
      if(PbcId(3)>=0) jmin_id=-10
      if(PbcId(4)>=0) jmax_id=-10
      if(PbcId(5)>=0) kmin_id=-10
      if(PbcId(6)>=0) kmax_id=-10
      ! assign name to each face
      face_names(1) = "imin"
      face_names(2) = "imax"
      face_names(3) = "jmin"
      face_names(4) = "jmax"
      face_names(5) = "kmin"
      face_names(6) = "kmax"
      
      id(1) =  imin_id
      id(2) =  imax_id
      id(3) =  jmin_id
      id(4) =  jmax_id
      id(5) =  kmin_id
      id(6) =  kmax_id

      c2 = 1 + scheme%accur
      c3 = 0.5*scheme%accur
      c1 = c2-c3
      call read_fixed_values(files, scheme, flow)

      call alloc(make_F_flux_zero, 1,dims%imx)
      call alloc(make_G_flux_zero, 1,dims%jmx)
      call alloc(make_H_flux_zero, 1,dims%kmx)

      make_F_flux_zero=1
      make_G_flux_zero=1
      make_H_flux_zero=1

      if(imin_id==-5 .or. imin_id==-6 .or. imin_id==-7) make_F_flux_zero(1)=0
      if(jmin_id==-5 .or. jmin_id==-6 .or. jmin_id==-7) make_G_flux_zero(1)=0
      if(kmin_id==-5 .or. kmin_id==-6 .or. kmin_id==-7) make_H_flux_zero(1)=0
      if(imax_id==-5 .or. imax_id==-6 .or. imax_id==-7) make_F_flux_zero(dims%imx)=0
      if(jmax_id==-5 .or. jmax_id==-6 .or. jmax_id==-7) make_G_flux_zero(dims%jmx)=0
      if(kmax_id==-5 .or. kmax_id==-6 .or. kmax_id==-7) make_H_flux_zero(dims%kmx)=0

    end subroutine setup_bc

!    subroutine destroy_bc()
!      !< Deallocate memory from boundary condition variables
!      implicit none
!      call dealloc(make_F_flux_zero)
!      call dealloc(make_G_flux_zero)
!      call dealloc(make_H_flux_zero)
!    end subroutine destroy_bc

end module bc
