    !< Common constant/parameters variables used by most other modules
module global
    !< Common constant/parameters variables used by most other modules

    implicit none
    ! String buffer lengths
    integer, parameter :: FILE_NAME_LENGTH = 64
    !< Length of string used for defining any filename
    integer, parameter :: SCHEME_NAME_LENGTH = 16
    !< Length of string used for storing Scheme
    integer, parameter :: INTERPOLANT_NAME_LENGTH = 10
    !< Length of string used for storing  higher order method
    integer, parameter :: DESCRIPTION_STRING_LENGTH = 64
    !< Length of string used for description in message call
    integer, parameter :: STRING_BUFFER_LENGTH = 128
    !< User to define a string of medium length
    integer, parameter :: ERROR_MESSAGE_LENGTH = 256
    !< Length of string used for passing error message during call
    integer, parameter :: LONG_BUFFER_LENGTH = 256
    !<Used to define a string of large size 
    integer, parameter :: FORMAT_LENGTH = 16
    !< Length of string used for file format: tecplot or vtk
    integer, parameter :: STATE_NAME_LENGTH = 64
    !< Length of string used in array user for sotring and reading Output/input variable list
    integer, parameter :: FLOW_TYPE_LENGTH = 64
    !< Length of string used for storing type of flow: inviscid, laminar, etc.
    integer, parameter :: TOLERANCE_LENGTH = 32
    !< Length of string used for resnorm types: abs or relative
    ! File unit numbers
    integer, parameter :: CONFIG_FILE_UNIT = 1
    !< Handler unit for config.md file
    integer, parameter :: GRID_FILE_UNIT = 2
    !< Handler for input Gridfile; eg: grid_00.txt
    integer, parameter :: BOUNDARY_CONDITIONS_FILE_UNIT= 3
    !< Handler for Boundary condition file; eg: bc_00.md
    integer, parameter :: STATE_FILE_UNIT = 10
    !< __Handler no longer in use__
    integer, parameter ::  IN_FILE_UNIT = 19
    !< Handler for restart file for block: eg: time_drectories/0001/process_00.dat
    integer, parameter :: OUT_FILE_UNIT = 20
    !< Handler for output file for each block
    integer, parameter :: RESNORM_FILE_UNIT = 21
    !< Handler for Residual output file. filename: time_directories/aux/resnorm
    integer, parameter ::   TEMP_NODE_FILE_UNIT = 30
    !< __Handler no longer in use__
    integer, parameter ::      LAYOUT_FILE_UNIT = 31
    !< Handler for input multi-block layout and boundary condition file.
    integer, parameter ::    NODESURF_FILE_UNIT = 32
    !< Handler for storing node point on the wall
    integer, parameter ::   WALL_DIST_FILE_UNIT = 33
    !< __Handler no longer in use__
    integer, parameter :: RES_CONTROL_FILE_UNIT = 34
    !< Handler for residual control file. filename: system/res_control.md
    integer, parameter ::        INFO_FILE_UNIT = 35
    !< __Handler NO longer in user__; info is handled using print*, command
    integer, parameter ::     CONTROL_FILE_UNIT = 36
    !< Handler for input system/control.md file
    integer, parameter ::      SCHEME_FILE_UNIT = 37
    !< Handler for input system/fvscheme.md file
    integer, parameter ::        FLOW_FILE_UNIT = 38
    !< Handler for input system/flow.md  file
    integer, parameter ::     RESTART_FILE_UNIT = 39
    !< Handler for Restart file in Restart folder. eg: time_directories/0001/Restart/process_00
    integer, parameter ::       OUTIN_FILE_UNIT = 40
    !< Handler for file which controls what variables will be read or stored. system/output_control.md
    integer, parameter ::         MAP_FILE_UNIT = 41
    !< Handler for input multi-block mapping file with index and direction.
    integer, parameter ::    PERIODIC_FILE_UNIT = 42
    !< Handler for input periodic boundary condition file

    ! stop file
    integer, parameter ::       STOP_FILE_UNIT = 41
    !< Handler for Stop file

    !Fixed file names
    character(len=*), Parameter :: control_file="system/control.md"
    !< FILENAME string: Control file
    character(len=*), Parameter ::  scheme_file="system/fvscheme.md"
    !< FILENAME string: Scheme file
    character(len=*), Parameter ::    flow_file="system/flow.md"
    !< FILENAME string: FLow condition file
    character(len=*), Parameter ::   outin_file="system/output_control.md"
    !< FILENAME string: Ouput/Input variable control file
    character(len=*), parameter :: layout_file='system/mesh/layout/layout.md'
    !< FILENAME string: Multiple layout/boundary condition file
    character(len=*), Parameter :: nodefile_temp="scratch.dat"
    !< FILENAME string: Temperory file for nodesurface points 
    character(len=*), Parameter :: surface_node_points='time_directories/aux/surfnode.dat'
    !< FILENAME string: Wall surface node points
    character(len=*), Parameter :: wall_dist_file='distance.vtk'
    !< FILENAME string: Wall distance for debug-- not in use anymore
    character(len=*), parameter :: res_control_file='system/res_control.md'
    !< FILENAME string: Residual write control file
    character(len=*), parameter :: resnorm_file='time_directories/aux/resnorm'
    !< FILENAME string: Residual output file
    character(len=*), parameter :: stop_file='system/stopfile'
    !< FILENAME string: Halt/stop file
    character(len=*), parameter :: mapfile='system/mesh/layout/mapping.txt'
    !< FILENAME string: Detailed multiblock mapping file with indicies and direction information at interface
    character(len=*), parameter :: periodicfile='system/mesh/layout/periodic.txt'
    !< FILENAME string: Multiblock periodic boundary condition detials

end module global
