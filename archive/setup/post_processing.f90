module post_processing
  use utils,          only: dmsg
  use grid,           only: imx, jmx, kmx
  use geometry,       only: volume
  use state,          only: density
  use state,          only: pressure
  use state,          only: R_gas
  use state,          only: gm
  use state,          only: density_inf
  use state,          only: pressure_inf

  implicit none
  private


  public :: fetch_entropy_l2_norm

  contains

    subroutine fetch_entropy_l2_norm()
      implicit none
      integer :: i,j,k
      real    :: entropy_error
      real    :: entropy_inf
      real    :: entropy

      call dmsg(5, 'post_processing', 'fetch_entropy_l2_norm')
      entropy       = 0.
      entropy_error = 0.
      entropy_inf   =  (R_gas/(gm - 1.))*log(pressure_inf/(density_inf**gm))

      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1
            entropy =  (R_gas/(gm - 1.))*log(pressure(i,j,k)/(density(i,j,k)**gm))
            entropy_error = entropy_error + (( entropy_inf - entropy )**2)*volume(i,j,k)
          end do
        end do
      end do

      entropy_error = sqrt( entropy_error/sum(volume) )

      if()
      open(667,file='entropy_temp', status="old", position="append", action="write")
      write(667,'(f0.16)') entropy_error
      close(667)

    end subroutine fetch_entropy_l2_norm

end module post_processing


