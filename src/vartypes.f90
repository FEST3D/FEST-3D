module vartypes
  !< Derived data types
  use iso_fortran_env, only : wp => real64
    implicit none
    integer :: process_id=0
    !< Id no. of each processor assinged by MPICH library
    integer, parameter :: FILE_NAME_LENGTH = 64
    !< Length of string used for defining any filename
    integer, parameter :: STRING_BUFFER_LENGTH = 128
    !< User to define a string of medium length
    integer, parameter :: FORMAT_LENGTH = 16
    !< Length of string used for file format: tecplot or vtk

    type, public :: nodetype
        real(wp) :: x
        real(wp) :: y
        real(wp) :: z
    end type nodetype


    type, public :: extent
      integer :: imx
      integer :: jmx
      integer :: kmx
      integer :: n_var
    end type extent


    type, public :: celltype
      real(wp) :: volume
      !< Store cell volume
      real(wp) :: centerx
      real(wp) :: centery
      real(wp) :: centerz
      !< Store Cell-center location 
    end type celltype


    type, public :: facetype
      real(wp) :: A
       !< Store magnitude of face area vector of direction faces
      real(wp) :: nx
      real(wp) :: ny
      real(wp) :: nz
       !< Store unit face normal vector for all faces 
    end type facetype

    
    type, public :: filetype
!        integer :: FILE_NAME_LENGTH = 64
!        !< Length of string used for defining any filename
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
        character(len=FILE_NAME_LENGTH) :: gridfile
        !< FILENAME string: Grid file
        character(len=FILE_NAME_LENGTH) :: bcfile
        !< FILENAME string: single block boundary condition detials
        character(len=FILE_NAME_LENGTH) :: outfile
        !< FILENAME string: single block solution output file
        character(len=FILE_NAME_LENGTH) :: infile
        !< FILENAME string: single block solution input file
        character(len=FILE_NAME_LENGTH) :: restartfile
        !< FILENAME string: single block restart file
    end type filetype


    type, public :: controltype
        real(wp) :: CFL=1.0
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
        integer :: res_write_interval=10
        !< Write resnorm after every "res_write_interval" iteration
        integer :: purge_write=1
        !< Remove unwanted folder. If Purge_write=2, latest two folder in time_direcotires are kept and 0=no purge
        integer :: last_iter=0    
        !< Last iteration that is stored in the restart file
        integer :: write_percision=6 
        !< Number of place after decimal. Only used for resnorm file
        character(len=FORMAT_LENGTH):: write_data_format='ASCII'
        !< write data type. Either ASCII or BINARY
        character(len=FORMAT_LENGTH):: write_file_format='tecplot'
        !< Write file type. Either vtk or tecplot
        character(len=FORMAT_LENGTH)::  read_data_format='ASCII' 
        !< Read data type in file. Either ASCII or BINARY
        character(len=FORMAT_LENGTH)::  read_file_format='tecplot'
        !< Read file type. Either vtk or tecplot
        real(wp) :: tolerance=1e-14
        !< Minimum value of resnorm after which simulation stop
        character(len=STRING_BUFFER_LENGTH) :: tolerance_type="abs"
        !< Type of tolerance to check:absolute or relative
        integer, public :: DEBUG_LEVEL = 1
        !< Debug level is an input from the control file.
        !< 5-> important calls only, and, 
        !< 1-> all the calls
        character(len=STRING_BUFFER_LENGTH):: previous_flow_type="none" 
        !< Type of flow:inviscid, laminar, etc, stored in the load file 
        integer                                           :: n_var=5
        ! Freestram variable used to read file before inf pointer are linked and allocated
        character(len=STRING_BUFFER_LENGTH), dimension(:), allocatable ::  r_list 
        !< Read variable list
        character(len=STRING_BUFFER_LENGTH), dimension(:), allocatable ::  w_list 
        !< Write variable list
        integer :: r_count=0                               
        !< Number of variable to read from the restart file
        integer :: w_count=0                               
        !< Number of variable to write in the output file
        character(len=STRING_BUFFER_LENGTH), dimension(:), allocatable :: Res_list
        !< Write residual variable list
        integer            :: Res_count       
        !< No of residual variable to save
        integer :: total_process=1
        !< Total number of process to be used for computation
        integer :: process_id=0
        !< Id no. of each processor assinged by MPICH library
        integer :: want_to_stop=0  
        !< 0: continue the solver; 1=Stop the solver
        logical :: Halt = .FALSE.
        !< Logical value used to stop the solver in main program file.
        real(wp), dimension(:), allocatable :: previous_res
        !< starting resnrom of previous run
    end type controltype

    type, public :: schemetype
      character(len=STRING_BUFFER_LENGTH) :: scheme_name='ausm'
      !< Flux Scheme to use: ausm, ldfss0, vanleer, ausmup, ausmp, slau
      character(len=STRING_BUFFER_LENGTH) :: interpolant='muscl'
      !< Face state reconstruction  method to user: muscl, ppm, none, weno, and wenoNM
      real(wp)                                      :: global_time_step=1e-5
      !< Value of global time step to march the solution with
      character                                 :: time_stepping_method='l'
      !< Either local time stepping or global time stepping
      character(len=STRING_BUFFER_LENGTH)    :: time_step_accuracy='implicit'
      !< Type of time_integration scheme: RK4, none(firt order explicit) implicit,
      integer                                           :: ilimiter_switch=1
       !< Turn on/off application of limiter for MUSCL (higer order face state reconstiion) for I direction faces.
      integer                                           :: jlimiter_switch=1
       !< Turn on/off application of limiter for MUSCL (higer order face state reconstiion) for J direction faces.
      integer                                           :: klimiter_switch=1
       !< Turn on/off application of limiter for MUSCL (higer order face state reconstiion) for K direction faces.
      integer                                           :: itlimiter_switch=1
       !< Turn on/off application of limiter for MUSCL (higer order face turbulent variable state reconstiion) for I direction faces.
      integer                                           :: jtlimiter_switch=1
       !< Turn on/off application of limiter for MUSCL (higer order face turbulent variable state reconstiion) for J direction faces.
      integer                                           :: ktlimiter_switch=1
       !< Turn on/off application of limiter for MUSCL (higer order face turbulent variable state reconstiion) for K direction faces.
      integer                                           :: iPB_switch=0
       !< Turn on/off application of pressure based switching for higher order methods for I direction faces.
      integer                                           :: jPB_switch=0
       !< Turn on/off application of pressure based switching for higher order methods for J direction faces.
      integer                                           :: kPB_switch=0
       !< Turn on/off application of pressure based switching for higher order methods for K direction faces.
      character(len=8)                                  :: turbulence='none'
       !< Store Turbulence model name
      character(len=8)                                  :: transition='none'
       !< Store Transition model name
      integer  :: accur=1                          
      !< Switch for higher order boundary condition
   end type schemetype


   type, public ::  flowtype
      real(wp)                                              :: density_inf=1.2
       !< Read freestream Density from control file
      real(wp)                                              :: x_speed_inf=100.0
       !< Read freestream U from control file
      real(wp)                                              :: y_speed_inf=0.0
       !< Read freestream V from control file
      real(wp)                                              :: z_speed_inf=0.0
       !< Read freestream W from control file
      real(wp)                                              :: pressure_inf=101325
       !< Read freestream Pressure from control file
      real(wp)                                              :: tk_inf=0.0
       !< Read freestream turbulent kinetic energy rate from control file
      real(wp)                                              :: tw_inf=0.0
       !< Read freestream turbulent dissipation rate from control file
      real(wp)                                              :: te_inf=0.0
       !< Read freestream turbulent dissipation from control file
      real(wp)                                              :: tv_inf=0.0
       !< Read freestream turbulent viscosity(SA) from control file
      real(wp)                                              :: tkl_inf=0.0
       !< Read freestream kL variable from control file
      real(wp)                                              :: tu_inf=1.0
       !< Read freestream turbulence intensity (percentage) from control file
      real(wp)                                              :: tgm_inf=1.0
       !< Read freestream turbulence intermittency from control file
      real(wp)                                              :: vel_mag=100.0
       !< Calulated freestream Velocity Magnitude from control file
      real(wp)                                              :: MInf=0.0
       !< Calulated freestream Mach number
      real(wp)                                              :: Reynolds_number=0.0
       !< Calculated free_stream Reynolds_number
      real(wp)                                              :: mu_ratio_inf=1.0
       !< Read freestream turbulent viscosity to molecular viscosity ratio
      real(wp)                                              :: Turb_intensity_inf=0.01
       !< Calculate free_stream turbulence intensity 
      real(wp)                                              :: gm=1.4    
      !< Gamma commonly 1.4
      real(wp)                                              :: R_gas=287 
      !< Univarsal gas constant
      real(wp)                                              :: mu_ref=0.0
      !< Molecular viscoity reference
      character(len=FILE_NAME_LENGTH)                   :: mu_variation="constant"
      !< Type of viscosity variaiton: Sutherland or constant
      real(wp)                                              :: T_ref=300
      !< Reference Temperature of flow for viscosity calculation
      real(wp)                                              :: Sutherland_temp=110
      !< Sutherland temperature for viscosity calculation
      real(wp)                                              :: Pr=0.7 
      !< prandtl number
      real(wp)                                              :: tPr=0.9 
      !< turbulent Prandtl number
    end type flowtype


    type :: boundarytype
      integer :: imin_id            
      !< Boundary condition number/ID at imin for particulat processor
      integer :: imax_id            
      !< Boundary condition number/ID at imax for particulat processor
      integer :: jmin_id            
      !< Boundary condition number/ID at jmin for particulat processor
      integer :: jmax_id            
      !< Boundary condition number/ID at jmax for particulat processor
      integer :: kmin_id           
      !< Boundary condition number/ID at kmin for particulat processor
      integer :: kmax_id  
      !< Boundary condition number/ID at kmax for particulat processor
      character(len=4), dimension(6) :: face_names 
      !< Store name of all six boundary faces
      integer,          dimension(6) :: id         
      !< Store the boundary condition ID of all six faces
      real(wp)                           :: c1         
      !< First coefficient user for higher order boundary condition
      real(wp)                           :: c2         
      !< Second coefficient user for higher order boundary condition
      real(wp)                           :: c3         
      !< Third coefficient user for higher order boundary condition

      ! store fix values for 6 faces of domain
      real(wp), dimension(6) :: fixed_density  = 0.
      !< Density value to fix at particular boundary
      real(wp), dimension(6) :: fixed_pressure = 0.
      !< Pressure value to fix at particular boundary
      real(wp), dimension(6) :: fixed_x_speed  = 0.
      !< X component of velocity to fix at particular boundary condition
      real(wp), dimension(6) :: fixed_y_speed  = 0.
      !< Y component of velocity to fix at particular boundary condition
      real(wp), dimension(6) :: fixed_z_speed  = 0.
      !< Z component of velocity to fix at particular boundary condition
      real(wp), dimension(6) :: fixed_tk       = 0.
      !< Turbulent kinetic energy value to fix at particular boundary condition
      real(wp), dimension(6) :: fixed_tw       = 0.
      !< Turbulent kinetic energy dissiaption rate value to fix at particular boundary condition(k-omega and SST model)
      real(wp), dimension(6) :: fixed_te       = 0.
      !< Turbulent kinetic energy dissiaption value to fix at particular boundary condition (K-eplision model)
      real(wp), dimension(6) :: fixed_tv       = 0.
      !< Turbulent viscosity varialble value to fix at particular boundary condition (for SA turbulence model)
      real(wp), dimension(6) :: fixed_tkl       = 0.
      !< (Turbulent kinetic energy x length) varialble value to fix at particular boundary condition (for k-kL turbulence model)
      real(wp), dimension(6) :: fixed_tgm       = 0.
      !<  Fixed intermittency value to apply at particular boundary condition (for SST2003-gamma transition model)
      real(wp), dimension(6) :: fixed_wall_temperature  = 0.
      !<  Fixed wall temperature value to apply at isothermal wall boundary condition.
      real(wp), dimension(6) :: fixed_Tpressure         = 0.
      !<  Fixed Total Pressure value to apply at particular boundary condition
      real(wp), dimension(6) :: fixed_Ttemperature      = 0.
      !<  Fixed Total Temperature value to apply at particular boundary condition


      !interface mapping
      integer, dimension(6) :: ilo, ihi 
       !< Store the lower and upper bound of the indecies of I loop for the interface mapping
      integer, dimension(6) :: jlo, jhi
       !< Store the lower and upper bound of the indecies of J loop for the interface mapping
      integer, dimension(6) :: klo, khi
       !< Store the lower and upper bound of the indecies of K loop for the interface mapping
      integer, dimension(6) :: dir_switch=0
       !< Switch for each boundary face. Activated only if ( for eg I-direction in the mapping is mapped with J-direction)
      integer, dimension(6) :: otherface
       !< Store the face number with which the current interface is connected.

      !zero flux faces
      integer,dimension(:),allocatable::make_F_flux_zero
       !< Store zero to boundary face, which has wall ID, to make F flux zero
      integer,dimension(:),allocatable::make_G_flux_zero
       !< Store zero to boundary face, which has wall ID, to make G flux zero
      integer,dimension(:),allocatable::make_H_flux_zero  
       !< Store zero to boundary face, which has wall ID, to make H flux zero

      !periodic boundary condition
      integer, dimension(6) :: PbcId = -1 !< Block ID for Periodic boundary condition

    end type boundarytype

end module vartypes
