module bc
  !--------------------------------------------
  ! 170515  Jatinder Pal Singh Sandhu
  ! Aim : setup boundary condition to domain
  !-------------------------------------------
  use global_vars, only: imin_id
  use global_vars, only: imax_id
  use global_vars, only: jmin_id
  use global_vars, only: jmax_id
  use global_vars, only: kmin_id
  use global_vars, only: kmax_id
  use global_vars, only: accur
  use global_vars, only: c1
  use global_vars, only: c2
  use global_vars, only: c3
  use global_vars, only: id
  use global_vars, only: face_names

  use read_bc   ,  only: read_fixed_values


  implicit none
  private

  integer                        :: face_num

  public :: setup_bc


  contains

    subroutine setup_bc()
      implicit none
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

      c2 = 1 + accur
      c3 = 0.5*accur
      c1 = c2-c3
      call read_fixed_values()

    end subroutine setup_bc


end module bc
