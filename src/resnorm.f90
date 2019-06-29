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

  use global,      only: RESNORM_FILE_UNIT
  use global,      only: resnorm_file

  use global_vars, only: imx
  use global_vars, only: jmx
  use global_vars, only: kmx
  use global_vars, only: gm
  use global_vars, only: n_var
  use global_vars, only: density_inf
  use global_vars, only: vel_mag
  use global_vars, only: pressure_inf
  use global_vars, only: tk_inf
  use global_vars, only: tw_inf
  use global_vars, only: tkl_inf
  use global_vars, only: te_inf
  use global_vars, only: tv_inf
  use global_vars, only: current_iter
  use global_vars, only: res_write_interval
  use global_vars, only: write_percision
  use global_vars, only: Res_abs
  use global_vars, only: Res_rel
  use global_vars, only: Res_save
  use global_vars, only: Res_scale
  use global_vars, only: Res_count
  use global_vars, only: Res_list
  use global_vars, only: Res_itr
  use global_vars, only: turbulence
  use global_vars, only: residue
  use global_vars, only: start_from
  use global_vars, only: last_iter
  use global_vars, only: F_p
  use global_vars, only: G_p
  use global_vars, only: H_p


  use utils,      only: dmsg
  use utils,      only: dealloc
  use utils,      only: alloc
  use layout,     only: process_id
  use layout,     only: total_process
  use string
  use fclose,     only: close_file

#include "error.inc"
#include "mpi.inc"
  private
  real :: merror
  real, dimension(:), allocatable :: buffer

  public :: setup_resnorm
  public :: destroy_resnorm
  public :: find_resnorm

  contains

    subroutine setup_resnorm()
      !< Allocate memory, setup scale and file to write
      implicit none
      call allocate_memory()
      call setup_scale()
      call setup_file()
    end subroutine setup_resnorm

    subroutine find_resnorm()
      !< Find the normalized residual for each processor
      implicit none
      call get_absolute_resnorm()
      call collect_resnorm_from_all_blocks()
      call assemble_resnom_at_each_process()
      call get_relative_resnorm()
      if((mod(current_iter,res_write_interval)==0 .or. &
              current_iter==Res_itr .or.               &
              current_iter==1)      .and.              &
              process_id  ==0)      then
        call write_resnorm()
      end if
    end subroutine find_resnorm

    subroutine destroy_resnorm()
      !< Deallocate memory and close residual file
      implicit none
      call deallocate_memory()
      call close_file(RESNORM_FILE_UNIT)
    end subroutine destroy_resnorm

    subroutine setup_file()
      !< Open the residual file to write
      implicit none
      integer :: i
      if(process_id==0)then
        if(start_from==0)then
          open(RESNORM_FILE_UNIT,file=resnorm_file)
        else
          open(RESNORM_FILE_UNIT,file=resnorm_file, status='old', position='append', action='write')
        end if
        write(RESNORM_FILE_UNIT, '(A,2x)', advance='no') "Iteration"
        do i=1,Res_count
          write(RESNORM_FILE_UNIT, '(A,2x)', advance='no') trim(Res_list(i))
        end do
        write(RESNORM_FILE_UNIT, *)
      end if
    end subroutine setup_file

    subroutine allocate_memory()
      !< Allocate memory to MPI Communication
      implicit none
      call alloc(Res_abs  , 0,n_var)
      call alloc(Res_rel  , 0,n_var)
      call alloc(Res_scale, 0,n_var)
      call alloc(Res_save , 0,n_var)
      call alloc(buffer   , 1,(n_var+1)*total_process)
    end subroutine allocate_memory

    subroutine deallocate_memory()
      !< Deallocate memory required for MPI Communication
      implicit none
      call dealloc(Res_abs)
      call dealloc(Res_rel)
      call dealloc(Res_scale)
      call dealloc(Res_save)
      call dealloc(buffer)
      if(allocated(Res_list)) deallocate(Res_list)
    end subroutine deallocate_memory

    subroutine setup_scale()
      !< Setup scale required for relative and absolute
      !< residual for writing in the file.
      implicit none
      Res_scale(0) = 1.
      Res_scale(1) = density_inf*vel_mag
      Res_scale(2) = density_inf*vel_mag*vel_mag
      Res_scale(3) = density_inf*vel_mag*vel_mag
      Res_scale(4) = density_inf*vel_mag*vel_mag
      Res_scale(5) = (0.5*density_inf*vel_mag**3 + &
                     ((gm/(gm-1.))*pressure_inf))

      select case(trim(turbulence))
        case('none')
          !do nothing
          continue
        case('sst', 'sst2003')
          Res_scale(6) = density_inf*vel_mag*tk_inf
          Res_scale(7) = density_inf*vel_mag*tw_inf
        case('kkl')
          Res_scale(6) = density_inf*vel_mag*tk_inf
          Res_scale(7) = density_inf*vel_mag*tkl_inf
        case('des')
          Res_scale(6) = density_inf*vel_mag*tk_inf
          Res_scale(7) = density_inf*vel_mag*tw_inf
        case('sa', 'saBC')
          Res_scale(6) = density_inf*vel_mag*tv_inf
        case('kw')
          Res_scale(6) = density_inf*vel_mag*tk_inf
          Res_scale(7) = density_inf*vel_mag*tw_inf
        case('ke')
          Res_scale(6) = density_inf*vel_mag*tk_inf
          Res_scale(7) = density_inf*vel_mag*te_inf
        case DEFAULT
          Fatal_error
      end select

    end subroutine setup_scale

    subroutine get_absolute_resnorm()
      !< Get absolute residual for current process
      implicit none
      integer :: i
      do i=1,n_var
        Res_abs(i) =(sum(Residue(:,:,:,i)**2)/Res_scale(i)**2)
      end do
      merror = (                                     &
               sum(F_p(  1,1:jmx-1,1:kmx-1,1)) &
              -sum(F_p(imx,1:jmx-1,1:kmx-1,1)) &
              +sum(G_p(1:imx-1,  1,1:kmx-1,1)) &
              -sum(G_p(1:imx-1,jmx,1:kmx-1,1)) &
              +sum(H_p(1:imx-1,1:jmx-1,  1,1)) &
              -sum(H_p(1:imx-1,1:jmx-1,kmx,1)) &
              )
      Res_abs(0) = (merror/Res_scale(0))
    end subroutine get_absolute_resnorm

    subroutine collect_resnorm_from_all_blocks()
      !< MPI Communication to gather residual from all processes
      implicit none
      integer :: ierr
      call MPI_ALLGATHER(Res_abs, n_var+1, MPI_DOUBLE_PRECISION, &
      buffer, n_var+1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
    end subroutine collect_resnorm_from_all_blocks

    subroutine assemble_resnom_at_each_process()
      !< Sum residual obtained from all the processes after MPI_Communication
      implicit none
      integer :: i,j
      Res_abs=0.
      do i=0,total_process-1
        do j = 0,n_var
          Res_abs(j) =  Res_abs(j)+buffer((j+1)+(n_var+1)*i)
        end do
      end do
      Res_abs(1:) = sqrt(Res_abs(1:))
      Res_abs(0) =  abs(Res_abs(0))
    end subroutine assemble_resnom_at_each_process

    subroutine get_relative_resnorm()
      !< Get relative residual with respect to first iteration residual
      implicit none
      if(current_iter<=Res_itr) Res_save=Res_abs
      Res_rel = Res_abs/Res_save
    end subroutine get_relative_resnorm

    subroutine write_resnorm()
      !< Writing the residual in the file to save.
      implicit none
      integer :: i
      integer :: n=6
      character(len=20) :: frm

      n=write_percision
      write(frm, '(A,I0,A,I0,A)') "(e",n+8,".",n,"E2, 4x)"

      write(RESNORM_FILE_UNIT, '(I0,4x)', advance='no') current_iter+last_iter
      do i=1,Res_count
        select case(trim(Res_list(i)))
          !include "resnorm_write_cases.inc"
          case('Mass_abs')
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(0)

          case('Resnorm_abs')
            write(RESNORM_FILE_UNIT, frm, advance='no') sqrt(sum(Res_abs(1:)**2))

          case('Viscous_abs')
            write(RESNORM_FILE_UNIT, frm, advance='no') sqrt(sum(Res_abs(1:5)**2))

          case('Turbulent_abs')
            if(trim(turbulence)/='none')then
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
            if(trim(turbulence)/='none')then
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
            if(trim(turbulence)=='sst' .or. trim(turbulence)=='kkl'.or. trim(turbulence)=='sst2003' )then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(6)
            end if

          case('Tv_abs')
            if(trim(turbulence)=='sa' .or. trim(turbulence)=='saBC')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(6)
            end if

          case('Dissipation_abs')
            if(trim(turbulence)=='ke')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(7)
            end if

          case('Omega_abs')
            if(trim(turbulence)=='sst'.or. trim(turbulence)=='sst2003')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(7)
            end if

          case('Kl_abs')
            if(trim(turbulence)=='kkl')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(7)
            end if

          case('TKE_rel')
            if(trim(turbulence)=='sst' .or. trim(turbulence)=='kkl'.or.  trim(turbulence)=='sst2003')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(6)
            end if

          case('Tv_rel')
            if(trim(turbulence)=='sa' .or. trim(turbulence)=='saBC')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(6)
            end if

          case('Dissipation_rel')
            if(trim(turbulence)=='ke')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(7)
            end if

          case('Omega_rel')
            if(trim(turbulence)=='sst'.or. trim(turbulence)=='sst2003')then
            write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(7)
            end if

          case('Kl_rel')
            if(trim(turbulence)=='kkl')then
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

