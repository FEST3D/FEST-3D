module vartypes
  use global
    implicit none

    type, public :: nodetype
        real :: x
        real :: y
        real :: z
    end type nodetype


    type, public :: extent
      integer :: imx
      integer :: jmx
      integer :: kmx
    end type extent


    type, public :: celltype
      real :: volume
      !< Store cell volume
      real :: centerx
      real :: centery
      real :: centerz
      !< Store Cell-center location 
    end type celltype


    type, public :: facetype
      real :: A
       !< Store magnitude of face area vector of direction faces
      real :: nx
      real :: ny
      real :: nz
       !< Store unit face normal vector for all faces 
    end type facetype

    
    type, public :: filetype
        integer :: FILE_NAME_LENGTH = 64
        !< Length of string used for defining any filename
        integer ::      CONFIG_FILE_UNIT = 16
        !< Handler unit for config.md file
        integer ::        GRID_FILE_UNIT = 17
        !< Handler for input Gridfile; eg: grid_00.txt
        integer :: BOUNDARY_CONDITIONS_FILE_UNIT= 18
        !< Handler for Boundary condition file; eg: bc_00.md
        !integer, parameter :: STATE_FILE_UNIT = 10
        !< __Handler no longer in use__
        integer ::          IN_FILE_UNIT = 19
        !< Handler for restart file for block: eg: time_drectories/0001/process_00.dat
        integer ::         OUT_FILE_UNIT = 20
        !< Handler for output file for each block
        integer ::     RESNORM_FILE_UNIT = 21
        !< Handler for Residual output file. filename: time_directories/aux/resnorm
        integer ::      LAYOUT_FILE_UNIT = 31
        !< Handler for input multi-block layout and boundary condition file.
        integer ::    NODESURF_FILE_UNIT = 32
        !< Handler for storing node point on the wall
        !integer, parameter ::   WALL_DIST_FILE_UNIT = 33
        !< __Handler no longer in use__
        integer :: RES_CONTROL_FILE_UNIT = 34
        !< Handler for residual control file. filename: system/res_control.md
        !integer, parameter ::        INFO_FILE_UNIT = 35
        !< __Handler NO longer in user__; info is handled using print*, command
        integer ::     CONTROL_FILE_UNIT = 36
        !< Handler for input system/control.md file
        integer ::      SCHEME_FILE_UNIT = 37
        !< Handler for input system/fvscheme.md file
        integer ::        FLOW_FILE_UNIT = 38
        !< Handler for input system/flow.md  file
        integer ::     RESTART_FILE_UNIT = 39
        !< Handler for Restart file in Restart folder. eg: time_directories/0001/Restart/process_00
        integer ::       OUTIN_FILE_UNIT = 40
        !< Handler for file which controls what variables will be read or stored. system/output_control.md
        integer ::         MAP_FILE_UNIT = 41
        !< Handler for input multi-block mapping file with index and direction.
        integer ::    PERIODIC_FILE_UNIT = 42
        !< Handler for input periodic boundary condition file
        integer ::       STOP_FILE_UNIT   = 43
        !< Handler for Stop file
        !file names
        character(len=FILE_NAME_LENGTH) :: control_file="system/control.md"
        !< FILENAME string: Control file
        character(len=FILE_NAME_LENGTH) ::  scheme_file="system/fvscheme.md"
        !< FILENAME string: Scheme file
        character(len=FILE_NAME_LENGTH) ::    flow_file="system/flow.md"
        !< FILENAME string: FLow condition file
        character(len=FILE_NAME_LENGTH) ::   outin_file="system/output_control.md"
        !< FILENAME string: Ouput/Input variable control file
        character(len=FILE_NAME_LENGTH) :: layout_file='system/mesh/layout/layout.md'
        !< FILENAME string: Multiple layout/boundary condition file
        character(len=FILE_NAME_LENGTH) :: nodefile_temp="scratch.dat"
        !< FILENAME string: Temperory file for nodesurface points 
        character(len=FILE_NAME_LENGTH) :: surface_node_points='time_directories/aux/surfnode.dat'
        !< FILENAME string: Wall surface node points
        character(len=FILE_NAME_LENGTH) :: res_control_file='system/res_control.md'
        !< FILENAME string: Residual write control file
        character(len=FILE_NAME_LENGTH) :: resnorm_file='time_directories/aux/resnorm'
        !< FILENAME string: Residual output file
        character(len=FILE_NAME_LENGTH) :: stop_file='system/stopfile'
        !< FILENAME string: Halt/stop file
        character(len=FILE_NAME_LENGTH) :: mapfile='system/mesh/layout/mapping.txt'
        !< FILENAME string: Detailed multiblock mapping file with indicies and direction information at interface
        character(len=FILE_NAME_LENGTH) :: periodicfile='system/mesh/layout/periodic.txt'
        !< FILENAME string: Multiblock periodic boundary condition detials
    end type filetype


    type, public :: controltype
        real :: CFL=1.0
        !< Courant–Friedrichs–Lewy (CFL) (Read from input)
        integer :: start_from=0 
        !< Number of the folder (in time_directories) to load stored state from to restart computation
        integer :: min_iter=1     
        !< Minimum iteration value, starting iteration value
        integer :: max_iters=1 
        !< Maximum iteration value, stop after these many iteration
        integer :: checkpoint_iter=0  
        !< Write interval for output file. Number of iteration after which solver will dump/store a state in a folder in time_directories
        integer :: checkpoint_iter_count=0 
        !< Counter of folder number to write in time_directories/
        integer :: current_iter=0  
        !< Current iteration number
        integer :: res_write_interval       
        !< Write resnorm after every "res_write_interval" iteration
        integer :: purge_write    
        !< Remove unwanted folder. If Purge_write=2, latest two folder in time_direcotires are kept and 0=no purge
        integer :: last_iter=0    
        !< Last iteration that is stored in the restart file
        integer :: write_percision=6 
        !< Number of place after decimal. Only used for resnorm file
        character(len=FORMAT_LENGTH):: write_data_format  
        !< write data type. Either ASCII or BINARY
        character(len=FORMAT_LENGTH):: write_file_format 
        !< Write file type. Either vtk or tecplot
        character(len=FORMAT_LENGTH)::  read_data_format='ASCII' 
        !< Read data type in file. Either ASCII or BINARY
        character(len=FORMAT_LENGTH)::  read_file_format="vtk" 
        !< Read file type. Either vtk or tecplot
        real :: tolerance
        !< Minimum value of resnorm after which simulation stop
        character(len=TOLERANCE_LENGTH) :: tolerance_type="abs"
        !< Type of tolerance to check:absolute or relative
        integer, public :: DEBUG_LEVEL = 1
        !< Debug level is an input from the control file.
        !< 5-> important calls only, and, 
        !< 1-> all the calls
        character(len=FLOW_TYPE_LENGTH):: previous_flow_type="none" 
        !< Type of flow:inviscid, laminar, etc, stored in the load file 
        integer                                           :: n_var=5
        ! Freestram variable used to read file before inf pointer are linked and allocated
        character(len=STATE_NAME_LENGTH), dimension(:), allocatable ::  r_list 
        !< Read variable list
        character(len=STATE_NAME_LENGTH), dimension(:), allocatable ::  w_list 
        !< Write variable list
        integer :: r_count=0                               
        !< Number of variable to read from the restart file
        integer :: w_count=0                               
        !< Number of variable to write in the output file
        character(len=STATE_NAME_LENGTH), dimension(:), allocatable :: Res_list
        !< Write residual variable list
        integer            :: Res_count       
        !< No of residual variable to save
    end type controltype

    type, public :: schemetype
      character(len=SCHEME_NAME_LENGTH) :: scheme_name    
      !< Flux Scheme to use: ausm, ldfss0, vanleer, ausmup, ausmp, slau
      character(len=INTERPOLANT_NAME_LENGTH) :: interpolant
      !< Face state reconstruction  method to user: muscl, ppm, none, weno, and wenoNM
      real                                      :: global_time_step   
      !< Value of global time step to march the solution with
      character                                 :: time_stepping_method
      !< Either local time stepping or global time stepping
      character(len=INTERPOLANT_NAME_LENGTH)    :: time_step_accuracy 
      !< Type of time_integration scheme: RK4, none(firt order explicit) implicit,
      integer                                           :: ilimiter_switch
       !< Turn on/off application of limiter for MUSCL (higer order face state reconstiion) for I direction faces.
      integer                                           :: jlimiter_switch
       !< Turn on/off application of limiter for MUSCL (higer order face state reconstiion) for J direction faces.
      integer                                           :: klimiter_switch
       !< Turn on/off application of limiter for MUSCL (higer order face state reconstiion) for K direction faces.
      integer                                           :: itlimiter_switch
       !< Turn on/off application of limiter for MUSCL (higer order face turbulent variable state reconstiion) for I direction faces.
      integer                                           :: jtlimiter_switch
       !< Turn on/off application of limiter for MUSCL (higer order face turbulent variable state reconstiion) for J direction faces.
      integer                                           :: ktlimiter_switch
       !< Turn on/off application of limiter for MUSCL (higer order face turbulent variable state reconstiion) for K direction faces.
      integer                                           :: iPB_switch
       !< Turn on/off application of pressure based switching for higher order methods for I direction faces.
      integer                                           :: jPB_switch
       !< Turn on/off application of pressure based switching for higher order methods for J direction faces.
      integer                                           :: kPB_switch
       !< Turn on/off application of pressure based switching for higher order methods for K direction faces.
      character(len=8)                                  :: turbulence
       !< Store Turbulence model name
      character(len=8)                                  :: transition
       !< Store Transition model name
      integer  :: accur=1                          
      !< Switch for higher order boundary condition
   end type schemetype


   type, public ::  flowtype
      real                                              :: density_inf
       !< Read freestream Density from control file
      real                                              :: x_speed_inf
       !< Read freestream U from control file
      real                                              :: y_speed_inf
       !< Read freestream V from control file
      real                                              :: z_speed_inf
       !< Read freestream W from control file
      real                                              :: pressure_inf
       !< Read freestream Pressure from control file
      real                                              :: tk_inf
       !< Read freestream turbulent kinetic energy rate from control file
      real                                              :: tw_inf       
       !< Read freestream turbulent dissipation rate from control file
      real                                              :: te_inf
       !< Read freestream turbulent dissipation from control file
      real                                              :: tv_inf
       !< Read freestream turbulent viscosity(SA) from control file
      real                                              :: tkl_inf
       !< Read freestream kL variable from control file
      real                                              :: tu_inf
       !< Read freestream turbulence intensity (percentage) from control file
      real                                              :: tgm_inf
       !< Read freestream turbulence intermittency from control file
      real                                              :: vel_mag
       !< Calulated freestream Velocity Magnitude from control file
      real                                              :: MInf
       !< Calulated freestream Mach number
      real                                              :: Reynolds_number 
       !< Calculated free_stream Reynolds_number
      real                                              :: mu_ratio_inf 
       !< Read freestream turbulent viscosity to molecular viscosity ratio
      real                                              :: Turb_intensity_inf 
       !< Calculate free_stream turbulence intensity 
      real                                              :: gm=1.4    
      !< Gamma commonly 1.4
      real                                              :: R_gas=287 
      !< Univarsal gas constant
      real                                              :: mu_ref=0.0
      !< Molecular viscoity reference
      character(len=FILE_NAME_LENGTH)                   :: mu_variation="constant"
      !< Type of viscosity variaiton: Sutherland or constant
      real                                              :: T_ref=300
      !< Reference Temperature of flow for viscosity calculation
      real                                              :: Sutherland_temp=110
      !< Sutherland temperature for viscosity calculation
      real                                              :: Pr=0.7 
      !< prandtl number
      real                                              :: tPr=0.9 
      !< turbulent Prandtl number
    end type flowtype


end module vartypes
