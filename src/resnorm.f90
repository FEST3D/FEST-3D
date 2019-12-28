  !< This module contains subroutine that 
  !< 1. check if time for resnorm dump is arrived
  !< 2. calculate resnorm
  !< 3. send those resnorm to processor number 0
  !< 4. Recalulate resnorm based on information 
  !<    availble from all processors
  !< 5. Append the data to resnorm file
module resnorm
  !< This module contains subroutine that 
  !< 1. check if time for resnorm dump is arrived
  !< 2. calculate resnorm
  !< 3. send those resnorm to processor number 0
  !< 4. Recalulate resnorm based on information 
  !<    availble from all processors
  !< 5. Append the data to resnorm file
  !----------------------------------------------------

  use vartypes
  use utils,      only: alloc

#include "error.inc"
#include "mpi.inc"
  private

  real(wp) :: merror
  !< mass error
  real(wp), dimension(:), allocatable :: buffer
  !< Buffer for mpi communication
  integer, parameter :: Res_itr = 3
  !< Iteration after which Res_save is stores
  real(wp), dimension(:), allocatable :: Res_abs       
  !< Absolute value of residual norm
  real(wp), dimension(:), allocatable :: Res_rel       
  !< Relative value of residual norm
  real(wp), dimension(:), allocatable :: Res_save      
  !< Saved iteration for relative values
  real(wp), dimension(:), allocatable :: Res_scale     
  !< Scaling factor for normalization

  public :: setup_resnorm
  public :: Res_abs, Res_rel
  public :: find_resnorm

  contains

    subroutine setup_resnorm(files, control, scheme, flow)
      !< Allocate memory, setup scale and file to write
      implicit none
      type(filetype), intent(in) :: files
      !< Files' name and handler
      type(controltype), intent(in) :: control
      !< Control parameters
      type(schemetype), intent(in) :: scheme
      !< finite-volume Schemes
      type(flowtype), intent(in) :: flow
      !< Information about fluid flow: freestream-{u,v,rho,p}, etc.
      call allocate_memory(control)
      call setup_scale(scheme, flow)
      call setup_file(files, control)
    end subroutine setup_resnorm

    subroutine find_resnorm(file_handler, residue, F,G,H, control, scheme, dims)
      !< Find the normalized residual for each processor
      implicit none
      integer, intent(in) :: file_handler
      !< residual file handler
      type(controltype), intent(inout) :: control
      !< Control parameters
      type(schemetype) , intent(in) :: scheme
      !< finite-volume Schemes
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      real(wp), dimension(:, :, :, :), intent(in)  :: residue
      !< Store residue at each cell-center
      real(wp), dimension(:, :, :, :), intent(in) :: F
      !< Store fluxes throught the I faces
      real(wp), dimension(:, :, :, :), intent(in) :: G
      !< Store fluxes throught the J faces
      real(wp), dimension(:, :, :, :), intent(in) :: H
      !< Store fluxes throught the K faces
      call get_absolute_resnorm(residue, F,G,H, control, dims)
      call collect_resnorm_from_all_blocks(control)
      call assemble_resnom_at_each_process(control)
      call get_relative_resnorm(control)
      if((mod(control%current_iter,control%res_write_interval)==0 .or. &
              control%current_iter==Res_itr .or.               &
              control%current_iter==1)      .and.              &
              process_id  ==0)      then
        call write_resnorm(file_handler, control, scheme)
      end if
    end subroutine find_resnorm


    subroutine setup_file(files, control)
      !< Open the residual file to write
      implicit none
      type(filetype), intent(in) :: files
      !< Files' name and handler
      type(controltype), intent(in) :: control
      !< Control parameters
      integer :: i
      if(process_id==0)then
        if(control%start_from==0)then
          open(files%RESNORM_FILE_UNIT,file=files%resnorm_file)
        else
          open(files%RESNORM_FILE_UNIT,file=files%resnorm_file, status='old', position='append', action='write')
        end if
        write(files%RESNORM_FILE_UNIT, '(A,2x)', advance='no') "Iteration"
        do i=1,control%Res_count
          write(files%RESNORM_FILE_UNIT, '(A,2x)', advance='no') trim(control%Res_list(i))
        end do
        write(files%RESNORM_FILE_UNIT, *)
      end if
    end subroutine setup_file

    subroutine allocate_memory(control)
      !< Allocate memory to MPI Communication
      implicit none
      type(controltype), intent(in) :: control
      call alloc(Res_abs  , 0,control%n_var)
      call alloc(Res_rel  , 0,control%n_var)
      call alloc(Res_scale, 0,control%n_var)
      call alloc(Res_save , 0,control%n_var)
      call alloc(buffer   , 1,(control%n_var+1)*control%total_process)
    end subroutine allocate_memory


    subroutine setup_scale(scheme, flow)
      !< Setup scale required for relative and absolute
      !< residual for writing in the file.
      implicit none
      type(schemetype), intent(in) :: scheme
      !< finite-volume Schemes
      type(flowtype), intent(in) :: flow
      !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
      Res_scale(0) = 1.
      Res_scale(1) = flow%density_inf*flow%vel_mag
      Res_scale(2) = flow%density_inf*flow%vel_mag*flow%vel_mag
      Res_scale(3) = flow%density_inf*flow%vel_mag*flow%vel_mag
      Res_scale(4) = flow%density_inf*flow%vel_mag*flow%vel_mag
      Res_scale(5) = (0.5*flow%density_inf*flow%vel_mag**3 + &
                     ((flow%gm/(flow%gm-1.))*flow%pressure_inf))

      select case(trim(scheme%turbulence))
        case('none')
          !do nothing
          continue
        case('sst', 'sst2003')
          Res_scale(6) = flow%density_inf*flow%vel_mag*flow%tk_inf
          Res_scale(7) = flow%density_inf*flow%vel_mag*flow%tw_inf
        case('kkl')
          Res_scale(6) = flow%density_inf*flow%vel_mag*flow%tk_inf
          Res_scale(7) = flow%density_inf*flow%vel_mag*flow%tkl_inf
        case('des')
          Res_scale(6) = flow%density_inf*flow%vel_mag*flow%tk_inf
          Res_scale(7) = flow%density_inf*flow%vel_mag*flow%tw_inf
        case('sa', 'saBC')
          Res_scale(6) = flow%density_inf*flow%vel_mag*flow%tv_inf
        case('kw')
          Res_scale(6) = flow%density_inf*flow%vel_mag*flow%tk_inf
          Res_scale(7) = flow%density_inf*flow%vel_mag*flow%tw_inf
        case('ke')
          Res_scale(6) = flow%density_inf*flow%vel_mag*flow%tk_inf
          Res_scale(7) = flow%density_inf*flow%vel_mag*flow%te_inf
        case DEFAULT
          Fatal_error
      end select

    end subroutine setup_scale

    subroutine get_absolute_resnorm(residue, F,G,H, control, dims)
      !< Get absolute residual for current process
      implicit none
      type(controltype), intent(in) :: control
      !< Control parameters: number of variables
      type(extent), intent(in) :: dims
      !< extent of the 3D domain
      real(wp), dimension(:, :, :, :), intent(in)  :: residue
      !< Store residue at each cell-center
      real(wp), dimension(:, :, :, :), intent(in) :: F
      !< Store fluxes throught the I faces
      real(wp), dimension(:, :, :, :), intent(in) :: G
      !< Store fluxes throught the J faces
      real(wp), dimension(:, :, :, :), intent(in) :: H
      !< Store fluxes throught the K faces
      integer :: i
      do i=1,control%n_var
        Res_abs(i) =(sum(Residue(:,:,:,i)**2)/Res_scale(i)**2)
      end do
      merror = (                                     &
               sum(F(  1,1:dims%jmx-1,1:dims%kmx-1,1)) &
              -sum(F(dims%imx,1:dims%jmx-1,1:dims%kmx-1,1)) &
              +sum(G(1:dims%imx-1,  1,1:dims%kmx-1,1)) &
              -sum(G(1:dims%imx-1,dims%jmx,1:dims%kmx-1,1)) &
              +sum(H(1:dims%imx-1,1:dims%jmx-1,  1,1)) &
              -sum(H(1:dims%imx-1,1:dims%jmx-1,dims%kmx,1)) &
              )
      Res_abs(0) = (merror/Res_scale(0))
    end subroutine get_absolute_resnorm

    subroutine collect_resnorm_from_all_blocks(control)
      !< MPI Communication to gather residual from all processes
      implicit none
      type(controltype), intent(in) :: control
      !< Control parameters: number of variables
      integer :: ierr
      call MPI_ALLGATHER(Res_abs, control%n_var+1, MPI_DOUBLE_PRECISION, &
      buffer, control%n_var+1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
    end subroutine collect_resnorm_from_all_blocks

    subroutine assemble_resnom_at_each_process(control)
      !< Sum residual obtained from all the processes after MPI_Communication
      implicit none
      type(controltype), intent(in) :: control
      !< Control parameters: number of variables and total mpi processes
      integer :: i,j
      Res_abs=0.
      do i=0,control%total_process-1
        do j = 0,control%n_var
          Res_abs(j) =  Res_abs(j)+buffer((j+1)+(control%n_var+1)*i)
        end do
      end do
      Res_abs(1:) = sqrt(Res_abs(1:))
      Res_abs(0) =  abs(Res_abs(0))
    end subroutine assemble_resnom_at_each_process

    subroutine get_relative_resnorm(control)
      !< Get relative residual with respect to first iteration residual
      implicit none
      type(controltype), intent(inout) :: control
      !< Control parameters: iterations
      if(control%current_iter<=Res_itr) Res_save=Res_abs
      if(control%start_from/=0) then
        Res_save=control%previous_Res
      else
        control%previous_Res = Res_save
      end if
      Res_rel = Res_abs/Res_save
    end subroutine get_relative_resnorm

    subroutine write_resnorm(RESNORM_FILE_UNIT, control, scheme)
      !< Writing the residual in the file to save.
      implicit none
      integer, intent(in) :: RESNORM_FILE_UNIT
      !<Resnorm file handler unit
      type(controltype), intent(in) :: control
      !< Control parameters
      type(schemetype) , intent(in) :: scheme
      !< turbulenca and transition schemes
      integer :: i
      integer :: n=6
      character(len=20) :: frm

      n=control%write_percision
      write(frm, '(A,I0,A,I0,A)') "(e",n+8,".",n,"E2, 4x)"

      write(RESNORM_FILE_UNIT, '(I0,4x)', advance='no') control%current_iter+control%last_iter
      do i=1,control%Res_count
        select case(trim(control%Res_list(i)))
          !include "resnorm_write_cases.inc"
          case('Mass_abs')
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(0)

          case('Resnorm_abs')
            write(RESNORM_FILE_UNIT, frm, advance='no') sqrt(sum(Res_abs(1:)**2))

          case('Viscous_abs')
            write(RESNORM_FILE_UNIT, frm, advance='no') sqrt(sum(Res_abs(1:5)**2))

          case('Turbulent_abs')
            if(trim(scheme%turbulence)/='none')then
            write(RESNORM_FILE_UNIT, frm, advance='no') sqrt(sum(Res_abs(6:)**2))
            end if

          case('Continuity_abs')
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(1)

          case('X_mom_abs')
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(2)

          case('Y_mom_abs')
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(3)

          case('Z_mom_abs')
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(4)

          case('Energy_abs')
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(5)

          case('Mass_rel')
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(0)

          case('Resnorm_rel')
            write(RESNORM_FILE_UNIT, frm, advance='no') sqrt(sum(Res_rel(1:)**2))

          case('Viscous_rel')
            write(RESNORM_FILE_UNIT, frm, advance='no') sqrt(sum(Res_rel(1:5)**2))

          case('Turbulent_rel')
            if(trim(scheme%turbulence)/='none')then
            write(RESNORM_FILE_UNIT, frm, advance='no') sqrt(sum(Res_rel(6:)**2))
            end if

          case('Continuity_rel')
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(1)

          case('X-mom_rel')
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(2)

          case('Y-mom_rel')
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(3)

          case('Z-mom_rel')
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(4)

          case('Energy_rel')
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(5)

          case('TKE_abs')
            if(trim(scheme%turbulence)=='sst' .or. trim(scheme%turbulence)=='kkl'.or. trim(scheme%turbulence)=='sst2003' )then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(6)
            end if

          case('Tv_abs')
            if(trim(scheme%turbulence)=='sa' .or. trim(scheme%turbulence)=='saBC')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(6)
            end if

          case('Dissipation_abs')
            if(trim(scheme%turbulence)=='ke')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(7)
            end if

          case('Omega_abs')
            if(trim(scheme%turbulence)=='sst'.or. trim(scheme%turbulence)=='sst2003')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(7)
            end if

          case('Kl_abs')
            if(trim(scheme%turbulence)=='kkl')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(7)
            end if

          case('TKE_rel')
            if(trim(scheme%turbulence)=='sst' .or. trim(scheme%turbulence)=='kkl'.or.  trim(scheme%turbulence)=='sst2003')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(6)
            end if

          case('Tv_rel')
            if(trim(scheme%turbulence)=='sa' .or. trim(scheme%turbulence)=='saBC')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(6)
            end if

          case('Dissipation_rel')
            if(trim(scheme%turbulence)=='ke')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(7)
            end if

          case('Omega_rel')
            if(trim(scheme%turbulence)=='sst'.or. trim(scheme%turbulence)=='sst2003')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(7)
            end if

          case('Kl_rel')
            if(trim(scheme%turbulence)=='kkl')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(7)
            end if

          case DEFAULT
            ! making absolute resnorm default
            write(RESNORM_FILE_UNIT, frm, advance='no') sqrt(sum(Res_abs(1:)**2))
            Issue_warning
        end select
      end do
      write(RESNORM_FILE_UNIT, *)

    end subroutine write_resnorm

end module resnorm

