   !< Setup boundary condition for the domain
module bc
   !< Setup boundary condition for the domain
  !-------------------------------------------
  use vartypes
!  use global_vars, only: imin_id
!  use global_vars, only: imax_id
!  use global_vars, only: jmin_id
!  use global_vars, only: jmax_id
!  use global_vars, only: kmin_id
!  use global_vars, only: kmax_id
!  use global_vars, only: c1
!  use global_vars, only: c2
!  use global_vars, only: c3
!  use global_vars, only: id
!  use global_vars, only: face_names
!  use global_vars, only: make_F_flux_zero
!  use global_vars, only: make_G_flux_zero
!  use global_vars, only: make_H_flux_zero
!  use global_vars, only: PbcId
  use utils, only: alloc

  use read_bc   ,  only: read_fixed_values


  implicit none
  private

  !integer                        :: face_num
  !< Number of the face : 1:imin, 2:imax, 3:jmin, 4:jmax, 5:kmin, 6:kmax

  public :: setup_bc
!  public :: destroy_bc


  contains

    subroutine setup_bc(files, scheme, flow, bc, dims)
      !< Initialization and allocate memory of boundary condition variables
      implicit none
      type(filetype), intent(in) :: files
      type(schemetype), intent(in) :: scheme
      type(flowtype), intent(in) :: flow
      type(boundarytype), intent(inout) :: bc
      type(extent), intent(in) :: dims
      !check for periodic bc
      if(bc%PbcId(1)>=0) bc%imin_id=-10
      if(bc%PbcId(2)>=0) bc%imax_id=-10
      if(bc%PbcId(3)>=0) bc%jmin_id=-10
      if(bc%PbcId(4)>=0) bc%jmax_id=-10
      if(bc%PbcId(5)>=0) bc%kmin_id=-10
      if(bc%PbcId(6)>=0) bc%kmax_id=-10
      ! assign name to each face
      bc%face_names(1) = "imin"
      bc%face_names(2) = "imax"
      bc%face_names(3) = "jmin"
      bc%face_names(4) = "jmax"
      bc%face_names(5) = "kmin"
      bc%face_names(6) = "kmax"
      
      bc%id(1) =  bc%imin_id
      bc%id(2) =  bc%imax_id
      bc%id(3) =  bc%jmin_id
      bc%id(4) =  bc%jmax_id
      bc%id(5) =  bc%kmin_id
      bc%id(6) =  bc%kmax_id

      bc%c2 = 1 + scheme%accur
      bc%c3 = 0.5*scheme%accur
      bc%c1 = bc%c2-bc%c3
      call read_fixed_values(files, scheme, flow, bc)

      call alloc(bc%make_F_flux_zero, 1,dims%imx)
      call alloc(bc%make_G_flux_zero, 1,dims%jmx)
      call alloc(bc%make_H_flux_zero, 1,dims%kmx)

      bc%make_F_flux_zero=1
      bc%make_G_flux_zero=1
      bc%make_H_flux_zero=1

      if(bc%imin_id==-5 .or. bc%imin_id==-6 .or. bc%imin_id==-7) bc%make_F_flux_zero(1)=0
      if(bc%jmin_id==-5 .or. bc%jmin_id==-6 .or. bc%jmin_id==-7) bc%make_G_flux_zero(1)=0
      if(bc%kmin_id==-5 .or. bc%kmin_id==-6 .or. bc%kmin_id==-7) bc%make_H_flux_zero(1)=0
      if(bc%imax_id==-5 .or. bc%imax_id==-6 .or. bc%imax_id==-7) bc%make_F_flux_zero(dims%imx)=0
      if(bc%jmax_id==-5 .or. bc%jmax_id==-6 .or. bc%jmax_id==-7) bc%make_G_flux_zero(dims%jmx)=0
      if(bc%kmax_id==-5 .or. bc%kmax_id==-6 .or. bc%kmax_id==-7) bc%make_H_flux_zero(dims%kmx)=0

    end subroutine setup_bc


end module bc
