  !< Close all the opened files
module fclose
  !< Close all the opened files
  !---------------------------------------------
  ! 170513  Jatinder Pal Singh Sandhu
  !  Aim : close all opened file.
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

  implicit none
  private

  public :: close_all_files
  public :: close_file

  contains

    subroutine close_all_files
      !< Call to close all files
      implicit none
      call close_file(      CONFIG_FILE_UNIT)
      call close_file(        GRID_FILE_UNIT)
      call close_file(       STATE_FILE_UNIT)
      call close_file(          IN_FILE_UNIT)
      call close_file(         OUT_FILE_UNIT)
      call close_file(     RESNORM_FILE_UNIT)
      call close_file(   TEMP_NODE_FILE_UNIT)
      call close_file(      LAYOUT_FILE_UNIT)
      call close_file(    NODESURF_FILE_UNIT)
      call close_file(   WALL_DIST_FILE_UNIT)
      call close_file( RES_CONTROL_FILE_UNIT)
      call close_file(        INFO_FILE_UNIT)
      call close_file(     CONTROL_FILE_UNIT)
      call close_file(      SCHEME_FILE_UNIT)
      call close_file(        FLOW_FILE_UNIT)
      call close_file(     RESTART_FILE_UNIT)
      call close_file(       OUTIN_FILE_UNIT)
      call close_file(        STOP_FILE_UNIT)
      call close_file( BOUNDARY_CONDITIONS_FILE_UNIT)
    end subroutine close_all_files

    subroutine close_file(handler)
      !Generalized subroutine to close single file
      implicit none
      integer, intent(in) :: handler
      logical :: ok
      inquire(handler, opened=ok)
      if(ok)  close(handler)
    end subroutine close_file

end module fclose
