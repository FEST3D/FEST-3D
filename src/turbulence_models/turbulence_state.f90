module turbulence_state
  use utils, only: dmsg 
  use state, only x_speed_inf, y_speed_inf, z_speed_infi, mu_ref
  use sst_state, only: sstsetup_state, sstdestroy_state, sstwritestate_vtk

  implicit none 
  private

  public :: setup_turbulence_state
  public :: destroy_turbulence_state

  contains

    subroutine setup_turbulence_state()
      implicit none
      real :: free_stream_tk, free_stream_tw
      call dmsg(1, 'turbulence/turbulence_state', 'setup_turbulence_state')
      
      select case (turbulence)
        case ('sst')
          !TODO check for units of mu_ref
          free_stream_tk = 1.5*((0.05**2)*(x_speed_inf**2 + &
                                           y_speed_inf**2 + &
                                           z_speed_inf**2)/3)
          free_stream_tw = free_stream_tk/(mu_ref*0.001)
          call sstsetup_state(free_stream_tk, free_stream_tw,state_load_file)
        case default
          call dmsg(5, 'turbulence/turbulence_state', 'setup_turbulece_state',&
                    'Turbulence model not recognised')
          stop
      end select
    end subroutine setup_turbulence_state

    subroutine destroy_turbulence_state()
      implicit none
      call dmsg(1, 'turbulence/turbulence_state', 'destroy_turbulence_state')
      
      select case (turbulence)
        case ('sst')
          call sstdestroy_state()
        case default
          call dmsg(5, 'turbulence/turbulence_state', 'setup_turbulece_state',&
                    'Turbulence model not recognised')
          stop
      end select
    end subroutine destroy_turbulence_state

end module turbulence_state
