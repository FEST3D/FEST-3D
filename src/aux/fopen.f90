  !< Open all files required 
module fopen
  !< Open all files required by the solver. Input and auxillary file
  !---------------------------------------------
  ! 17082  Jatinder Pal Singh Sandhu
  !  Aim : efficient open of files only with handler
  !---------------------------------------------
    ! File unit numbers
  use global, only:      CONFIG_FILE_UNIT
  use global, only:        GRID_FILE_UNIT
  use global, only:       STATE_FILE_UNIT
  use global, only:          IN_FILE_UNIT
  use global, only:         OUT_FILE_UNIT
  use global, only:     RESNORM_FILE_UNIT
  use global, only:   TEMP_NODE_FILE_UNIT
  use global, only:      LAYOUT_FILE_UNIT
  use global, only:    NODESURF_FILE_UNIT
  use global, only:   WALL_DIST_FILE_UNIT
  use global, only: RES_CONTROL_FILE_UNIT
  use global, only:        INFO_FILE_UNIT
  use global, only:     CONTROL_FILE_UNIT
  use global, only:      SCHEME_FILE_UNIT
  use global, only:        FLOW_FILE_UNIT
  use global, only:     RESTART_FILE_UNIT
  use global, only:       OUTIN_FILE_UNIT
  use global, only:        STOP_FILE_UNIT
  use global, only: BOUNDARY_CONDITIONS_FILE_UNIT
  use fclose, only: close_file

  implicit none
  private

!  public :: open_all_files
  public :: open_file

  contains

    subroutine open_file(handler)
      !< Open single file
      implicit none
      integer, intent(in) :: handler
      select case(handler)
        case(1)
          close_file(handler)
          open(handler,)
        case(2)
      end select
    end subroutine open_file

end module fopen
