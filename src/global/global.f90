module global
    !-------------------------------------------------------------------
    ! The global module holds global variables
    !
    ! Global variables include string buffer lengths, file unit lengths,
    ! etc.
    !-------------------------------------------------------------------

    implicit none
    ! String buffer lengths
    integer, parameter :: FILE_NAME_LENGTH = 64
    integer, parameter :: SCHEME_NAME_LENGTH = 16
    integer, parameter :: INTERPOLANT_NAME_LENGTH = 10
    integer, parameter :: DESCRIPTION_STRING_LENGTH = 64
    integer, parameter :: STRING_BUFFER_LENGTH = 128
    integer, parameter :: ERROR_MESSAGE_LENGTH = 256
    integer, parameter :: LONG_BUFFER_LENGTH = 256
    integer, parameter :: FORMAT_LENGTH = 16
    integer, parameter :: STATE_NAME_LENGTH = 64
    integer, parameter :: FLOW_TYPE_LENGTH = 64
    integer, parameter :: TOLERANCE_LENGTH = 32
    ! File unit numbers
    integer, parameter :: CONFIG_FILE_UNIT = 1
    integer, parameter :: GRID_FILE_UNIT = 2
    integer, parameter :: BOUNDARY_CONDITIONS_FILE_UNIT= 3
    integer, parameter :: STATE_FILE_UNIT = 10
    integer, parameter ::  IN_FILE_UNIT = 19
    integer, parameter :: OUT_FILE_UNIT = 20
    integer, parameter :: RESNORM_FILE_UNIT = 21
    integer, parameter :: SPHERE_INDICES_FILE_UNIT = 22
    integer, parameter ::   TEMP_NODE_FILE_UNIT = 30
    integer, parameter ::      LAYOUT_FILE_UNIT = 31
    integer, parameter ::    NODESURF_FILE_UNIT = 32
    integer, parameter ::   WALL_DIST_FILE_UNIT = 33
    integer, parameter :: RES_CONTROL_FILE_UNIT = 34
    integer, parameter ::        INFO_FILE_UNIT = 35
    integer, parameter ::     CONTROL_FILE_UNIT = 36
    integer, parameter ::      SCHEME_FILE_UNIT = 37
    integer, parameter ::        FLOW_FILE_UNIT = 38
    integer, parameter ::     RESTART_FILE_UNIT = 39
    integer, parameter ::       OUTIN_FILE_UNIT = 40
    integer, parameter ::         MAP_FILE_UNIT = 41
    integer, parameter ::    PERIODIC_FILE_UNIT = 42

    ! stop file
    integer, parameter ::       STOP_FILE_UNIT = 41

    !Fixed file names
    character(len=*), Parameter :: control_file="system/control.md"
    character(len=*), Parameter ::  scheme_file="system/fvscheme.md"
    character(len=*), Parameter ::    flow_file="system/flow.md"
    character(len=*), Parameter ::   outin_file="system/output_control.md"
    character(len=*), parameter :: layout_file='system/mesh/layout/layout.md'
    character(len=*), Parameter :: nodefile_temp="scratch.dat"
    character(len=*), Parameter :: surface_node_points='time_directories/aux/surfnode.dat'
    character(len=*), Parameter :: wall_dist_file='distance.vtk'
    character(len=*), parameter :: res_control_file='system/res_control.md'
    character(len=*), parameter :: resnorm_file='time_directories/aux/resnorm'
    character(len=*), parameter :: stop_file='system/stopfile'
    character(len=*), parameter :: mapfile='system/mesh/layout/mapping.txt'
    character(len=*), parameter :: periodicfile='system/mesh/layout/periodic.txt'

    ! Fixed variable
    integer, parameter :: resnorm_number=10 

    character(len=FILE_NAME_LENGTH), dimension(101:114)   :: &
      files=(/                                               &!no. unit
              'system/mesh/gridfiles/',                      &!1    101
              'system/mesh/bc/',                             &!2    102
              'system/mesh/layout/layout.md',                &!3    103
              'system/contorl.md',                           &!4    104
              'system/fvscheme.md',                          &!5    105
              'system/flow.md',                              &!6    106
              'system/output_control.md',                    &!7    107
              'system/res_control.md',                       &!8    108
              'system/stopfile',                             &!9    109
              'time_directories/aux/surfnode.dat',           &!10   110
              'time_directories/aux/resnorm',                &!11   111
              'time_directories/',                           &!12   112
              'time_directories/',                           &!13   113
              'time_directories/'                            &!14   114
            /)
    integer, parameter :: gridfile        = 101
    integer, parameter :: bcfile          = 102
    integer, parameter :: layoutfile      = 103
    integer, parameter :: controlfile     = 104
    integer, parameter :: fvschemefile    = 105
    integer, parameter :: flowfile        = 106
    integer, parameter :: outcontrolfile  = 107
    integer, parameter :: rescontrolfile  = 108
    integer, parameter :: stopfile        = 109
    integer, parameter :: surfnodefile    = 110
    integer, parameter :: resnorm         = 111
    integer, parameter :: soldumpfile     = 112
    integer, parameter :: solreadfile     = 113
    integer, parameter :: restartfile     = 114
end module global
