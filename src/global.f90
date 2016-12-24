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
    ! File unit numbers
    integer, parameter :: CONFIG_FILE_UNIT = 1
    integer, parameter :: GRID_FILE_UNIT = 2
    integer, parameter :: BOUNDARY_CONDITIONS_FILE_UNIT= 3
    integer, parameter :: STATE_FILE_UNIT = 10
    integer, parameter :: OUT_FILE_UNIT = 20
    integer, parameter :: RESNORM_FILE_UNIT = 21
    integer, parameter :: SPHERE_INDICES_FILE_UNIT = 22
    integer, parameter :: TEMP_NODE_FILE_UNIT = 30
    integer, parameter :: LAYOUT_FILE_UNIT = 31
    integer, parameter :: NODESURF_FILE_UNIT=32
    integer, parameter :: WALL_DIST_FILE_UNIT=33

    !Fixed file names
    character(len=*), Parameter :: nodefile_temp="scratch.dat"
    character(len=*), Parameter :: surface_node_points='surfnode.dat'
    character(len=*), Parameter :: wall_dist_file='distance.vtk'

end module global
