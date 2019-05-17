  !< A module contains generalized subroutine to copy variable in ghost cells
module copy_bc
  !< A module contains generalized subroutine to copy variable in ghost cells
  !--------------------------------------------
  ! 170515  Jatinder Pal Singh Sandhu
  ! Aim : applying boundary condition to domain
  !-------------------------------------------
 
  use global_vars, only: imx
  use global_vars, only: jmx
  use global_vars, only: kmx
  use global_vars, only: c1
  use global_vars, only: c2
  use global_vars, only: c3

  implicit none
  private

  public :: copy1
  public :: copy3

  contains

    subroutine copy1(var, type, face)
      !< Copy 1 layer of interior cell to first ghost cell layer
      implicit none
      character(len=*), intent(in) :: face
      !< Face over which boundary condition is being called.
      character(len=*), intent(in) :: type
      !< type of copy: flat, symmetry, anti-symmetry
      real, dimension(-2:imx+2, -2:jmx+2, -2:kmx+2), intent(inout) :: var
      !< Varible over which these operation has to be performed
      real :: a2=1

      select case(type)
        case("anti")
          a2 = -1
        case("symm")
          a2 =  1
        case DEFAULT
          print*, "ERROR: Wrong boundary condition type"
      end select

      select case(face)
        case("imin")
            var(      0, 1:jmx-1, 1:kmx-1) = a2*var(     1, 1:jmx-1, 1:kmx-1)
        case("imax")
            var(  imx  , 1:jmx-1, 1:kmx-1) = a2*var( imx-1, 1:jmx-1, 1:kmx-1)
        case("jmin")
            var(1:imx-1,       0, 1:kmx-1) = a2*var(1:imx-1,      1, 1:kmx-1)
        case("jmax")
            var(1:imx-1,   jmx  , 1:kmx-1) = a2*var(1:imx-1,  jmx-1, 1:kmx-1)
        case("kmin")
            var(1:imx-1, 1:jmx-1,       0) = a2*var(1:imx-1, 1:jmx-1,      1)
        case("kmax")
            var(1:imx-1, 1:jmx-1,   kmx  ) = a2*var(1:imx-1, 1:jmx-1,  kmx-1)
        case DEFAULT
          print*, "ERROR: wrong face for boundary condition"
      end select
    end subroutine copy1

    
    subroutine copy3(var, type, face)
      !< Copy 3 layer of interior cell to three ghost cell layer
      implicit none
      character(len=*), intent(in) :: face
      !< Face over which boundary condition is being called.
      character(len=*), intent(in) :: type
      !< type of copy: flat, symmetry, anti-symmetry
      real, dimension(-2:imx+2, -2:jmx+2, -2:kmx+2), intent(inout) :: var
      !< Varible over which these operation has to be performed

      real :: a1=1
      real :: a2=1
      real :: a3=0

      integer :: i1=1
      integer :: i2=2
      integer :: i3=3

      select case(type)
        case("anti")
          a1 =  1. ; i1 = 1
          a2 = -1. ; i2 = 2
          a3 =  0. ; i3 = 3
        case("flat")
          a1 =  1. ; i1 = 1
          a2 =  1. ; i2 = 1
          a3 =  0. ; i3 = 1
        case("symm")
          a1 =  c1 ; i1 = 1
          a2 =  c2 ; i2 = 2
          a3 =  c3 ; i3 = 3
          ! do nothing
          ! use default value
          continue
        case DEFAULT
          print*, "ERROR: Wrong boundary condition type"
      end select

      select case(face)
        case("imin")
            var(      0, 1:jmx-1, 1:kmx-1) = (a2*var(     i1, 1:jmx-1, 1:kmx-1)-a3*var(    i1+1, 1:jmx-1, 1:kmx-1))/a1
            var(     -1, 1:jmx-1, 1:kmx-1) = (a2*var(     i2, 1:jmx-1, 1:kmx-1)-a3*var(    i2+1, 1:jmx-1, 1:kmx-1))/a1
            var(     -2, 1:jmx-1, 1:kmx-1) = (a2*var(     i3, 1:jmx-1, 1:kmx-1)-a3*var(    i3+1, 1:jmx-1, 1:kmx-1))/a1
        case("imax")
            var(  imx  , 1:jmx-1, 1:kmx-1) = (a2*var( imx-i1, 1:jmx-1, 1:kmx-1)-a3*var(imx-i1-1, 1:jmx-1, 1:kmx-1))/a1
            var(  imx+1, 1:jmx-1, 1:kmx-1) = (a2*var( imx-i2, 1:jmx-1, 1:kmx-1)-a3*var(imx-i2-1, 1:jmx-1, 1:kmx-1))/a1
            var(  imx+2, 1:jmx-1, 1:kmx-1) = (a2*var( imx-i3, 1:jmx-1, 1:kmx-1)-a3*var(imx-i3-1, 1:jmx-1, 1:kmx-1))/a1
        case("jmin")
            var(1:imx-1,       0, 1:kmx-1) = (a2*var(1:imx-1,      i1, 1:kmx-1)-a3*var(1:imx-1 ,    i1+1, 1:kmx-1))/a1
            var(1:imx-1,      -1, 1:kmx-1) = (a2*var(1:imx-1,      i2, 1:kmx-1)-a3*var(1:imx-1 ,    i2+1, 1:kmx-1))/a1
            var(1:imx-1,      -2, 1:kmx-1) = (a2*var(1:imx-1,      i3, 1:kmx-1)-a3*var(1:imx-1 ,    i3+1, 1:kmx-1))/a1
        case("jmax")
            var(1:imx-1,   jmx  , 1:kmx-1) = (a2*var(1:imx-1,  jmx-i1, 1:kmx-1)-a3*var(1:imx-1 ,jmx-i1-1, 1:kmx-1))/a1
            var(1:imx-1,   jmx+1, 1:kmx-1) = (a2*var(1:imx-1,  jmx-i2, 1:kmx-1)-a3*var(1:imx-1 ,jmx-i2-1, 1:kmx-1))/a1
            var(1:imx-1,   jmx+2, 1:kmx-1) = (a2*var(1:imx-1,  jmx-i3, 1:kmx-1)-a3*var(1:imx-1 ,jmx-i3-1, 1:kmx-1))/a1
        case("kmin")
            var(1:imx-1, 1:jmx-1,       0) = (a2*var(1:imx-1, 1:jmx-1,      i1)-a3*var(1:imx-1 , 1:jmx-1,    i1+1))/a1
            var(1:imx-1, 1:jmx-1,      -1) = (a2*var(1:imx-1, 1:jmx-1,      i2)-a3*var(1:imx-1 , 1:jmx-1,    i2+1))/a1
            var(1:imx-1, 1:jmx-1,      -2) = (a2*var(1:imx-1, 1:jmx-1,      i3)-a3*var(1:imx-1 , 1:jmx-1,    i3+1))/a1
        case("kmax")
            var(1:imx-1, 1:jmx-1,   kmx  ) = (a2*var(1:imx-1, 1:jmx-1,  kmx-i1)-a3*var(1:imx-1 , 1:jmx-1,kmx-i1-1))/a1
            var(1:imx-1, 1:jmx-1,   kmx+1) = (a2*var(1:imx-1, 1:jmx-1,  kmx-i2)-a3*var(1:imx-1 , 1:jmx-1,kmx-i2-1))/a1
            var(1:imx-1, 1:jmx-1,   kmx+2) = (a2*var(1:imx-1, 1:jmx-1,  kmx-i3)-a3*var(1:imx-1 , 1:jmx-1,kmx-i3-1))/a1
        case DEFAULT
          print*, "ERROR: wrong face for boundary condition"
      end select
    end subroutine copy3

end module copy_bc
