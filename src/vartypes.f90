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
    end type controltype

    type, public :: schemetype
      character(len=SCHEME_NAME_LENGTH) :: scheme_name    
      !< Flux Scheme to use: ausm, ldfss0, vanleer, ausmup, ausmp, slau
      character(len=INTERPOLANT_NAME_LENGTH) :: interpolant
      !< Face state reconstruction  method to user: muscl, ppm, none, weno, and wenoNM
      !solver time secific
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
   end type schemetype


   type, public ::  flowtype
      real                                              :: free_stream_density  
       !< Read freestream Density from control file
      real                                              :: free_stream_x_speed  
       !< Read freestream U from control file
      real                                              :: free_stream_y_speed  
       !< Read freestream V from control file
      real                                              :: free_stream_z_speed  
       !< Read freestream W from control file
      real                                              :: free_stream_pressure 
       !< Read freestream Pressure from control file
      real                                              :: free_stream_tk       
       !< Read freestream turbulent kinetic energy rate from control file
      real                                              :: free_stream_tw       
       !< Read freestream turbulent dissipation rate from control file
      real                                              :: free_stream_te
       !< Read freestream turbulent dissipation from control file
      real                                              :: free_stream_tv
       !< Read freestream turbulent viscosity(SA) from control file
      real                                              :: free_stream_tkl
       !< Read freestream kL variable from control file
      real                                              :: free_stream_tu
       !< Read freestream turbulence intensity (percentage) from control file
      real                                              :: free_stream_tgm
       !< Read freestream turbulence intermittency from control file
      real                                              :: vel_mag
       !< Calulated freestream Velocity Magnitude from control file
      real                                              :: Reynolds_number 
       !< Calculated free_stream Reynolds_number
      real                                              :: mu_ratio_inf 
       !< Read freestream turbulent viscosity to molecular viscosity ratio
      real                                              :: Turb_intensity_inf 
       !< Calculate free_stream turbulence intensity 
      ! thermal properties
      real                                              :: gm    !< Gamma commonly 1.4
      real                                              :: R_gas !< Univarsal gas constant
      
      ! Transport properties
      real                                              :: mu_ref 
      !< Molecular viscoity reference
      character(len=FILE_NAME_LENGTH)                   :: mu_variation
      !< Type of viscosity variaiton: Sutherland or constant

      ! sutherland law variable
      real                                              :: T_ref
      !< Reference Temperature of flow for viscosity calculation
      real                                              :: Sutherland_temp
      !< Sutherland temperature for viscosity calculation

      ! nondimensional numbers
      real                                              :: Pr=0.7 !< prandtl number
      real                                              :: tPr=0.9 !< turbulent Prandtl number

      ! switches
      logical                                           :: supersonic_flag
       !< Switch for boundary condition. No longer in use
    end type flowtype

end module vartypes
