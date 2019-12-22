  !< Preconditioned LU-SGS scheme
  !<  maxtrix-free implicit time-integration method
  !<  for low speed flows
module plusgs
  !<
  !< Reference: Kitamura, K., Shima, E., Fujimoto, K. and Wang, Z.J.,
  !< Performance of low-dissipation Euler fluxes and preconditioned LU-SGS 
  !< at low speeds, Communications in Computational Physics, vol. 10 no. 1, pp.90-119, 2011
  !-----------------------------------------------
  use vartypes
  use global_kkl , only : cphi1
  use global_kkl , only : cphi2
  use global_kkl , only : fphi
  use global_kkl , only : eta
  use global_kkl , only : cd1
  use global_kkl , only : cmu
  use global_sst , only : beta1
  use global_sst , only : beta2
  use global_sst , only : bstar
  use global_sst , only : sst_F1
  use global_sa  , only : sigma_sa
  use global_sa , only : cb1
  use global_sa , only : cb2
  use global_sa , only : cw1
  use global_sa , only : cw2
  use global_sa , only : cw3
  use global_sa , only : cv1
  use global_sa , only : sigma_sa
  use global_sa , only : kappa_sa
  use global_sa , only : cv1_3
  use global_sa , only : cw3_6
  use global_vars, only : DCCVnX
  use global_vars, only : DCCVnY
  use global_vars, only : DCCVnZ
  use global_vars, only : CCnormalX
  use global_vars, only : CCnormalY
  use global_vars, only : CCnormalZ

!  use global_vars, only : volume
!  use global_vars, only : xnx, xny, xnz !face unit normal x
!  use global_vars, only : ynx, yny, ynz !face unit normal y
!  use global_vars, only : znx, zny, znz !face unit normal z
!  use global_vars, only : xA, yA, zA    !face area
  use global_vars, only : dist
  use global_vars, only : mu
  use global_vars, only : mu_t
  use global_vars, only : delta_t
  use global_vars, only : process_id
!  use global_vars, only: residue

  use gradients  , only: gradu_x
  use gradients  , only: gradu_y
  use gradients  , only: gradu_z
  use gradients  , only: gradv_x
  use gradients  , only: gradv_y
  use gradients  , only: gradv_z
  use gradients  , only: gradw_x
  use gradients  , only: gradw_y
  use gradients  , only: gradw_z
!  use geometry   , only: CellCenter

  use utils, only: alloc

  !--- sst implicit update ---!
  use global_sst, only : sst_F1
  use global_sst, only : sigma_k1
  use global_sst, only : sigma_k2
  use global_sst, only : sigma_w1
  use global_sst, only : sigma_w2
  
  use global_kkl, only : sigma_k
  use global_kkl, only : sigma_phi

  !-------mapping --------------!
  use mapping, only : PiDir
  use mapping, only : PjDir
  use mapping, only : PkDir
  use mapping, only : Pilo
  use mapping, only : Pjlo
  use mapping, only : Pklo
  use mapping, only : Pihi
  use mapping, only : Pjhi
  use mapping, only : Pkhi
  use mapping, only : mpi_class
  use global_vars, only: imin_id
  use global_vars, only: jmin_id
  use global_vars, only: kmin_id
  use global_vars, only: imax_id
  use global_vars, only: jmax_id
  use global_vars, only: kmax_id
  use global_vars, only : dir_switch
  use global_vars, only: PbcId


#include "debug.h"
#include "error.h"
#include "mpi.inc"


  real, dimension(:,:,:,:), allocatable :: delQ
  !< Change of state variable (solution) over one time-step
  real, dimension(:,:,:,:), allocatable :: delQstar
  !< Intermediate change of state variable over one time-step
  real, dimension(:,:,:), allocatable, target :: dummy
  !< Dummy variable
  real, dimension(:,:,:), pointer :: tmu
  !< Pointer to turbulent viscosity
  real, dimension(:,:,:), pointer :: mmu
  !< Pointer to molecular viscosity

  !parallel communication
  integer :: ibuf_size
  !< Size of the buffer for I face interface
  integer :: jbuf_size
  !< Size of the buffer for J face interface
  integer :: kbuf_size
  !< Size of the buffer for K face interface
  real, dimension(:), allocatable :: imin_send_buf
  !< Array to store data to send data for Imin face
  real, dimension(:), allocatable :: jmin_send_buf
  !< Array to store data to send data for Jmin face
  real, dimension(:), allocatable :: kmin_send_buf
  !< Array to store data to send data for Kmin face
  real, dimension(:), allocatable :: imin_recv_buf
  !< Array to store data to receive data for Imin face
  real, dimension(:), allocatable :: jmin_recv_buf
  !< Array to store data to receive data for Jmin face
  real, dimension(:), allocatable :: kmin_recv_buf
  !< Array to store data to receive data for Kmin face
  real, dimension(:), allocatable :: imax_send_buf
  !< Array to store data to send data for Imax face
  real, dimension(:), allocatable :: jmax_send_buf
  !< Array to store data to send data for Jmax face
  real, dimension(:), allocatable :: kmax_send_buf
  !< Array to store data to send data for Kmax face
  real, dimension(:), allocatable :: imax_recv_buf
  !< Array to store data to receive data for Imax face
  real, dimension(:), allocatable :: jmax_recv_buf
  !< Array to store data to receive data for Jmax face
  real, dimension(:), allocatable :: kmax_recv_buf
  !< Array to store data to receive data for Kmax face
  integer :: imx, jmx, kmx, n_var
  real :: gm, mu_ref, Reynolds_number, free_stream_tu
  real :: tk_inf
  real :: tkl_inf
  real :: tPr, Pr, R_gas
  real :: MInf

  public :: update_with_plusgs
  public :: setup_plusgs
!  public :: destroy_plusgs

  contains

    subroutine setup_plusgs(control, scheme, flow, dims)
      !< Allocate array memory for data communication
      implicit none
      type(controltype), intent(in) :: control
      type(schemetype), intent(in) :: scheme
      type(flowtype), intent(in) :: flow
      type(extent), intent(in) :: dims
      character(len=*), parameter :: &
        errmsg="module: LUSGS, subrouinte setup"

      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx

      n_var = control%n_var
      gm = flow%gm
      mu_ref = flow%mu_ref
      Reynolds_number = flow%Reynolds_number
      free_stream_tu = flow%tu_inf
      tk_inf = flow%tk_inf
      tkl_inf = flow%tkl_inf
      tpr = flow%tpr
      pr = flow%pr
      MInf = flow%MInf
      R_gas = flow%R_gas

      ibuf_size = (jmx-1)*(kmx-1)*n_var*1
      jbuf_size = (imx-1)*(kmx-1)*n_var*1
      kbuf_size = (imx-1)*(jmx-1)*n_var*1
      call alloc(imin_send_buf,1,ibuf_size, errmsg)
      call alloc(jmin_send_buf,1,jbuf_size, errmsg)
      call alloc(kmin_send_buf,1,kbuf_size, errmsg)
      call alloc(imin_recv_buf,1,ibuf_size, errmsg)
      call alloc(jmin_recv_buf,1,jbuf_size, errmsg)
      call alloc(kmin_recv_buf,1,kbuf_size, errmsg)
      call alloc(imax_send_buf,1,ibuf_size, errmsg)
      call alloc(jmax_send_buf,1,jbuf_size, errmsg)
      call alloc(kmax_send_buf,1,kbuf_size, errmsg)
      call alloc(imax_recv_buf,1,ibuf_size, errmsg)
      call alloc(jmax_recv_buf,1,jbuf_size, errmsg)
      call alloc(kmax_recv_buf,1,kbuf_size, errmsg)

      call alloc(delQ, 0, imx, 0, jmx, 0, kmx, 1, n_var)
      call alloc(delQstar, 0, imx, 0, jmx, 0, kmx, 1, n_var)


      if(mu_ref==0.0 .or. scheme%turbulence=='none') then
        call alloc(dummy, 0, imx, 0, jmx, 0, kmx)
        dummy = 0.0
      end if
      if(mu_ref==0.0)then
        mmu => dummy
      else
        mmu => mu
      end if
      if(trim(scheme%turbulence)=='none')then
        tmu => dummy
      else
        tmu => mu_t
      end if
    end subroutine setup_plusgs

!    subroutine destroy_plusgs()
!      !< Unallocate the memory required by LU-SGS module
!      implicit none
!      call dealloc(imin_send_buf)
!      call dealloc(jmin_send_buf)
!      call dealloc(kmin_send_buf)
!      call dealloc(imin_recv_buf)
!      call dealloc(jmin_recv_buf)
!      call dealloc(kmin_recv_buf)
!      call dealloc(imax_send_buf)
!      call dealloc(jmax_send_buf)
!      call dealloc(kmax_send_buf)
!      call dealloc(imax_recv_buf)
!      call dealloc(jmax_recv_buf)
!      call dealloc(kmax_recv_buf)
!      call dealloc(delQ)
!      call dealloc(delQstar)
!      call dealloc(dummy)
!    end subroutine destroy_plusgs

    subroutine update_with_plusgs(qp, cells, Ifaces,Jfaces,Kfaces,residue, scheme, dims)
      !< Time-integrate with LU_SGS method
      implicit none
      type(schemetype), intent(in) :: scheme
      type(extent), intent(in) :: dims
      real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(inout) :: qp
      real, dimension(:, :, :, :), intent(in)  :: residue
      !< Store residue at each cell-center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal

      select case(trim(scheme%turbulence))
        case('none')
          call update_laminar_variables(qp, residue, cells, Ifaces, Jfaces, Kfaces, dims)

        case('sst', 'sst2003')
          select case(trim(scheme%transition))
            case('none', 'bc')
              call update_SST_variables(qp, residue, cells, Ifaces, Jfaces, Kfaces, dims)
            case('lctm2015')
              call update_lctm2015(qp, residue, cells, Ifaces, Jfaces, Kfaces, dims)
            case DEFAULT
              Fatal_error
          end select

        case('kkl')
!          call update_KKL_variables()

        case('sa', 'saBC')
          call update_SA_variables(qp, residue, cells, Ifaces, Jfaces, Kfaces, dims)

        case Default
          Fatal_error

      end select


    end subroutine update_with_plusgs



    subroutine update_laminar_variables(qp, residue, cells, Ifaces, Jfaces, Kfaces, dims)
      !< Update laminar flow with LU-SGS scheme
      implicit none
      type(extent), intent(in) :: dims
      real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(inout) :: qp
      real, dimension(:, :, :, :), intent(in)  :: residue
      !< Store residue at each cell-center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      integer :: i,j,k
        real, dimension(1:5)     :: deltaU
        real                     :: D
        real, dimension(1:5)     :: conservativeQ
        real, dimension(1:5)     :: OldIminusFlux
        real, dimension(1:5)     :: OldJminusFlux
        real, dimension(1:5)     :: OldKminusFlux
        real, dimension(1:5)     :: NewIminusFlux
        real, dimension(1:5)     :: NewJminusFlux
        real, dimension(1:5)     :: NewKminusFlux
        real, dimension(1:5)     :: DelIminusFlux
        real, dimension(1:5)     :: DelJminusFlux
        real, dimension(1:5)     :: DelKminusFlux
        real, dimension(1:6)     :: LambdaTimesArea
        real, dimension(1:5)     :: Q0 ! state at cell
        real, dimension(1:5)     :: Q1 ! state at neighbours 
        real, dimension(1:5)     :: Q2
        real, dimension(1:5)     :: Q3
        real, dimension(1:5)     :: Q4
        real, dimension(1:5)     :: Q5
        real, dimension(1:5)     :: Q6
        real, dimension(1:5)     :: DQ0! change in state
        real, dimension(1:5)     :: DQ1
        real, dimension(1:5)     :: DQ2
        real, dimension(1:5)     :: DQ3
        real, dimension(1:5)     :: DQ4
        real, dimension(1:5)     :: DQ5
        real, dimension(1:5)     :: DQ6
        real, dimension(1:7)     :: Flist1
        real, dimension(1:7)     :: Flist2
        real, dimension(1:7)     :: Flist3
        real, dimension(1:7)     :: Flist4
        real, dimension(1:7)     :: Flist5
        real, dimension(1:7)     :: Flist6
        real, dimension(1:3)     :: C0
        real, dimension(1:3)     :: C1
        real, dimension(1:3)     :: C2
        real, dimension(1:3)     :: C3
        real, dimension(1:3)     :: C4
        real, dimension(1:3)     :: C5
        real, dimension(1:3)     :: C6
        real                     :: eps
        real                     :: M
        real                     :: VMag
        real                     :: SoundMag
        real                     :: u,v,w,r,p, H
        real                     :: factor
        real, dimension(1:5,1:5) :: PrecondInv



        !intialize delQ
        delQstar = 0.0

        !forward sweep
        do k=1,dims%kmx-1
          do j=1,dims%jmx-1
            do i=1,dims%imx-1
              C0  = (/Cells(i  ,j  ,k  )%Centerx,Cells(i  ,j  ,k  )%Centery,Cells(i  ,j  ,k  )%Centerz/)
              C1  = (/Cells(i-1,j  ,k  )%Centerx,Cells(i-1,j  ,k  )%Centery,Cells(i-1,j  ,k  )%Centerz/)
              C2  = (/Cells(i  ,j-1,k  )%Centerx,Cells(i  ,j-1,k  )%Centery,Cells(i  ,j-1,k  )%Centerz/)
              C3  = (/Cells(i  ,j  ,k-1)%Centerx,Cells(i  ,j  ,k-1)%Centery,Cells(i  ,j  ,k-1)%Centerz/)
              C4  = (/Cells(i+1,j  ,k  )%Centerx,Cells(i+1,j  ,k  )%Centery,Cells(i+1,j  ,k  )%Centerz/)
              C5  = (/Cells(i  ,j+1,k  )%Centerx,Cells(i  ,j+1,k  )%Centery,Cells(i  ,j+1,k  )%Centerz/)
              C6  = (/Cells(i  ,j  ,k+1)%Centerx,Cells(i  ,j  ,k+1)%Centery,Cells(i  ,j  ,k+1)%Centerz/)

              Q0  = qp(i  , j  , k  , 1:5)
              Q1  = qp(i-1, j  , k  , 1:5)
              Q2  = qp(i  , j-1, k  , 1:5)
              Q3  = qp(i  , j  , k-1, 1:5)
              Q4  = qp(i+1, j  , k  , 1:5)
              Q5  = qp(i  , j+1, k  , 1:5)
              Q6  = qp(i  , j  , k+1, 1:5)

              DQ0 = 0.0
              DQ1 = delQstar(i-1, j  , k  , 1:5)
              DQ2 = delQstar(i  , j-1, k  , 1:5)
              DQ3 = delQstar(i  , j  , k-1, 1:5)

              Flist1(1) =  Ifaces(i,j,k)%A
              Flist1(2) = -Ifaces(i,j,k)%nx
              Flist1(3) = -Ifaces(i,j,k)%ny
              Flist1(4) = -Ifaces(i,j,k)%nz
              Flist1(5) = 0.5*(cells(i-1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))

              Flist2(1) =  Jfaces(i,j,k)%A
              Flist2(2) = -Jfaces(i,j,k)%nx
              Flist2(3) = -Jfaces(i,j,k)%ny
              Flist2(4) = -Jfaces(i,j,k)%nz
              Flist2(5) = 0.5*(cells(i  , j-1, k  )%volume + cells(i,j,k)%volume)
              Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))

              Flist3(1) =  Kfaces(i,j,k)%A
              Flist3(2) = -Kfaces(i,j,k)%nx
              Flist3(3) = -Kfaces(i,j,k)%ny
              Flist3(4) = -Kfaces(i,j,k)%nz
              Flist3(5) = 0.5*(cells(i  , j  , k-1)%volume + cells(i,j,k)%volume)
              Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))

              Flist4(1) =  Ifaces(i+1,j,k)%A
              Flist4(2) = +Ifaces(i+1,j,k)%nx
              Flist4(3) = +Ifaces(i+1,j,k)%ny
              Flist4(4) = +Ifaces(i+1,j,k)%nz
              Flist4(5) = 0.5*(cells(i+1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))

              Flist5(1) =  Jfaces(i,j+1,k)%A
              Flist5(2) = +Jfaces(i,j+1,k)%nx
              Flist5(3) = +Jfaces(i,j+1,k)%ny
              Flist5(4) = +Jfaces(i,j+1,k)%nz
              Flist5(5) = 0.5*(cells(i  , j+1, k  )%volume + cells(i,j,k)%volume)
              Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))

              Flist6(1) =  Kfaces(i,j,k+1)%A
              Flist6(2) = +Kfaces(i,j,k+1)%nx
              Flist6(3) = +Kfaces(i,j,k+1)%ny
              Flist6(4) = +Kfaces(i,j,k+1)%nz
              Flist6(5) = 0.5*(cells(i  , j  , k+1)%volume + cells(i,j,k)%volume)
              Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))

              NewIminusFlux     = Flux(Q1, Q0, DQ1, Flist1)
              NewJminusFlux     = Flux(Q2, Q0, DQ2, Flist2)
              NewKminusFlux     = Flux(Q3, Q0, DQ3, Flist3)
              OldIminusFlux     = Flux(Q1, Q0, DQ0, Flist1)
              OldJminusFlux     = Flux(Q2, Q0, DQ0, Flist2)
              OldKminusFlux     = Flux(Q3, Q0, DQ0, Flist3)

              r  = Q0(1)
              u  = Q0(2)
              v  = Q0(3)
              w  = Q0(4)
              p  = Q0(5)
              VMag     = sqrt(u*u + v*v + w*w)
              SoundMag = sqrt(gm*p/r)
              M        = VMag/SoundMag 
              H  = (gm*p/(r*(gm-1.0))) + 0.5*(VMag)
              eps = min(1.0, max(M*M, Minf*Minf))
              factor = (1.0-eps)*(gm-1.0)/(SoundMag*SoundMag)
              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0, eps)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0, eps)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0, eps)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0, eps)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0, eps)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0, eps)


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux


              D = (cells(i,j,k)%volume/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
              !storing D in Iflux array for backward sweep
              !F_p(i,j,k,1) = D
              PrecondInv(1,1) = 1.0 - factor*1*VMag*VMag/2.0
              PrecondInv(2,1) = 0.0 - factor*u*VMag*VMag/2.0
              PrecondInv(3,1) = 0.0 - factor*v*VMag*VMag/2.0
              PrecondInv(4,1) = 0.0 - factor*w*VMag*VMag/2.0
              PrecondInv(5,1) = 0.0 - factor*H*VMag*VMag/2.0
              PrecondInv(1,2) = 0.0 - factor*1*(-u)
              PrecondInv(2,2) = 1.0 - factor*u*(-u)
              PrecondInv(3,2) = 0.0 - factor*v*(-u)
              PrecondInv(4,2) = 0.0 - factor*w*(-u)
              PrecondInv(5,2) = 0.0 - factor*H*(-u)
              PrecondInv(1,3) = 0.0 - factor*1*(-v)
              PrecondInv(2,3) = 0.0 - factor*u*(-v)
              PrecondInv(3,3) = 1.0 - factor*v*(-v)
              PrecondInv(4,3) = 0.0 - factor*w*(-v)
              PrecondInv(5,3) = 0.0 - factor*H*(-v)
              PrecondInv(1,4) = 0.0 - factor*1*(-w)
              PrecondInv(2,4) = 0.0 - factor*u*(-w)
              PrecondInv(3,4) = 0.0 - factor*v*(-w)
              PrecondInv(4,4) = 1.0 - factor*w*(-w)
              PrecondInv(5,4) = 0.0 - factor*H*(-w)
              PrecondInv(1,5) = 0.0 - factor*1*(1.)
              PrecondInv(2,5) = 0.0 - factor*u*(1.)
              PrecondInv(3,5) = 0.0 - factor*v*(1.)
              PrecondInv(4,5) = 0.0 - factor*w*(1.)
              PrecondInv(5,5) = 1.0 - factor*H*(1.)


              !deltaU(1:5) = -residue(i,j,k,1:5) &
              deltaU(1:5) = -matmul(PrecondInv,residue(i,j,k,1:5)) &
                - 0.5*((matmul(PrecondInv,DelIminusFlux) - LambdaTimesArea(1)*delQstar(i-1,j,k,1:5)) &
                     + (matmul(PrecondInv,DelJminusFlux) - LambdaTimesArea(2)*delQstar(i,j-1,k,1:5)) &
                     + (matmul(PrecondInv,DelKminusFlux) - LambdaTimesArea(3)*delQstar(i,j,k-1,1:5)) )

              delQstar(i,j,k,1:5) = deltaU(1:5)/D
            end do
          end do
        end do

        delQ=0.0
        !backward sweep
            do i=dims%imx-1,1,-1
          do j=dims%jmx-1,1,-1
        do k=dims%kmx-1,1,-1
              C0  = (/Cells(i  ,j  ,k  )%Centerx,Cells(i  ,j  ,k  )%Centery,Cells(i  ,j  ,k  )%Centerz/)
              C1  = (/Cells(i-1,j  ,k  )%Centerx,Cells(i-1,j  ,k  )%Centery,Cells(i-1,j  ,k  )%Centerz/)
              C2  = (/Cells(i  ,j-1,k  )%Centerx,Cells(i  ,j-1,k  )%Centery,Cells(i  ,j-1,k  )%Centerz/)
              C3  = (/Cells(i  ,j  ,k-1)%Centerx,Cells(i  ,j  ,k-1)%Centery,Cells(i  ,j  ,k-1)%Centerz/)
              C4  = (/Cells(i+1,j  ,k  )%Centerx,Cells(i+1,j  ,k  )%Centery,Cells(i+1,j  ,k  )%Centerz/)
              C5  = (/Cells(i  ,j+1,k  )%Centerx,Cells(i  ,j+1,k  )%Centery,Cells(i  ,j+1,k  )%Centerz/)
              C6  = (/Cells(i  ,j  ,k+1)%Centerx,Cells(i  ,j  ,k+1)%Centery,Cells(i  ,j  ,k+1)%Centerz/)

              Q0  = qp(i  , j  , k  , 1:5)
              Q1  = qp(i-1, j  , k  , 1:5)
              Q2  = qp(i  , j-1, k  , 1:5)
              Q3  = qp(i  , j  , k-1, 1:5)
              Q4  = qp(i+1, j  , k  , 1:5)
              Q5  = qp(i  , j+1, k  , 1:5)
              Q6  = qp(i  , j  , k+1, 1:5)

              DQ0 = 0.0
              DQ4 = delQ(i+1, j  , k  , 1:5)
              DQ5 = delQ(i  , j+1, k  , 1:5)
              DQ6 = delQ(i  , j  , k+1, 1:5)

              Flist1(1) =  Ifaces(i,j,k)%A
              Flist1(2) = -Ifaces(i,j,k)%nx
              Flist1(3) = -Ifaces(i,j,k)%ny
              Flist1(4) = -Ifaces(i,j,k)%nz
              Flist1(5) = 0.5*(cells(i-1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))

              Flist2(1) =  Jfaces(i,j,k)%A
              Flist2(2) = -Jfaces(i,j,k)%nx
              Flist2(3) = -Jfaces(i,j,k)%ny
              Flist2(4) = -Jfaces(i,j,k)%nz
              Flist2(5) = 0.5*(cells(i  , j-1, k  )%volume + cells(i,j,k)%volume)
              Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))

              Flist3(1) =  Kfaces(i,j,k)%A
              Flist3(2) = -Kfaces(i,j,k)%nx
              Flist3(3) = -Kfaces(i,j,k)%ny
              Flist3(4) = -Kfaces(i,j,k)%nz
              Flist3(5) = 0.5*(cells(i  , j  , k-1)%volume + cells(i,j,k)%volume)
              Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))

              Flist4(1) =  Ifaces(i+1,j,k)%A
              Flist4(2) = +Ifaces(i+1,j,k)%nx
              Flist4(3) = +Ifaces(i+1,j,k)%ny
              Flist4(4) = +Ifaces(i+1,j,k)%nz
              Flist4(5) = 0.5*(cells(i+1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))

              Flist5(1) =  Jfaces(i,j+1,k)%A
              Flist5(2) = +Jfaces(i,j+1,k)%nx
              Flist5(3) = +Jfaces(i,j+1,k)%ny
              Flist5(4) = +Jfaces(i,j+1,k)%nz
              Flist5(5) = 0.5*(cells(i  , j+1, k  )%volume + cells(i,j,k)%volume)
              Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))

              Flist6(1) =  Kfaces(i,j,k+1)%A
              Flist6(2) = +Kfaces(i,j,k+1)%nx
              Flist6(3) = +Kfaces(i,j,k+1)%ny
              Flist6(4) = +Kfaces(i,j,k+1)%nz
              Flist6(5) = 0.5*(cells(i  , j  , k+1)%volume + cells(i,j,k)%volume)
              Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))

              NewIminusFlux     = Flux(Q4, Q0, DQ4, Flist4)
              NewJminusFlux     = Flux(Q5, Q0, DQ5, Flist5)
              NewKminusFlux     = Flux(Q6, Q0, DQ6, Flist6)
              OldIminusFlux     = Flux(Q4, Q0, DQ0, Flist4)
              OldJminusFlux     = Flux(Q5, Q0, DQ0, Flist5)
              OldKminusFlux     = Flux(Q6, Q0, DQ0, Flist6)


              r  = Q0(1)
              u  = Q0(2)
              v  = Q0(3)
              w  = Q0(4)
              p  = Q0(5)
              VMag     = sqrt(u*u + v*v + w*w)
              SoundMag = sqrt(gm*p/r)
              M        = VMag/SoundMag 
              H  = (gm*p/(r*(gm-1.0))) + 0.5*(VMag)
              eps = min(1.0, max(M*M, Minf*Minf))
              factor = (1.0-eps)*(gm-1.0)/(SoundMag*SoundMag)

              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0, eps)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0, eps)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0, eps)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0, eps)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0, eps)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0, eps)


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux

              D = (cells(i,j,k)%volume/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)

              PrecondInv(1,1) = 1.0 - factor*1*VMag*VMag/2.0
              PrecondInv(2,1) = 0.0 - factor*u*VMag*VMag/2.0
              PrecondInv(3,1) = 0.0 - factor*v*VMag*VMag/2.0
              PrecondInv(4,1) = 0.0 - factor*w*VMag*VMag/2.0
              PrecondInv(5,1) = 0.0 - factor*H*VMag*VMag/2.0
              PrecondInv(1,2) = 0.0 - factor*1*(-u)
              PrecondInv(2,2) = 1.0 - factor*u*(-u)
              PrecondInv(3,2) = 0.0 - factor*v*(-u)
              PrecondInv(4,2) = 0.0 - factor*w*(-u)
              PrecondInv(5,2) = 0.0 - factor*H*(-u)
              PrecondInv(1,3) = 0.0 - factor*1*(-v)
              PrecondInv(2,3) = 0.0 - factor*u*(-v)
              PrecondInv(3,3) = 1.0 - factor*v*(-v)
              PrecondInv(4,3) = 0.0 - factor*w*(-v)
              PrecondInv(5,3) = 0.0 - factor*H*(-v)
              PrecondInv(1,4) = 0.0 - factor*1*(-w)
              PrecondInv(2,4) = 0.0 - factor*u*(-w)
              PrecondInv(3,4) = 0.0 - factor*v*(-w)
              PrecondInv(4,4) = 1.0 - factor*w*(-w)
              PrecondInv(5,4) = 0.0 - factor*H*(-w)
              PrecondInv(1,5) = 0.0 - factor*1*(1.)
              PrecondInv(2,5) = 0.0 - factor*u*(1.)
              PrecondInv(3,5) = 0.0 - factor*v*(1.)
              PrecondInv(4,5) = 0.0 - factor*w*(1.)
              PrecondInv(5,5) = 1.0 - factor*H*(1.)


              delQ(i,j,k,1:5) = delQstar(i,j,k,1:5) &
                - 0.5*((matmul(PrecondInv,DelIminusFlux) - LambdaTimesArea(4)*delQ(i+1,j,k,1:5)) &
                     + (matmul(PrecondInv,DelJminusFlux) - LambdaTimesArea(5)*delQ(i,j+1,k,1:5)) &
                     + (matmul(PrecondInv,DelKminusFlux) - LambdaTimesArea(6)*delQ(i,j,k+1,1:5)) )/D

            end do
          end do
        end do
        
        do k=1,dims%kmx-1
          do j = 1,dims%jmx-1
            do i = 1,dims%imx-1
              conservativeQ(1) = qp(i,j,k,1)
              conservativeQ(2) = qp(i,j,k,1) * qp(i,j,k,2)
              conservativeQ(3) = qp(i,j,k,1) * qp(i,j,k,3)
              conservativeQ(4) = qp(i,j,k,1) * qp(i,j,k,4)
              conservativeQ(5) = (qp(i,j,k,5) / (gm-1.0)) + ( 0.5 * qp(i,j,k,1) * sum( qp(i,j,k,2:4)**2) )
              
              ! add new change into conservative solution
              conservativeQ(1:5) = conservativeQ(1:5) + delQ(i,j,k,1:5)

              ! convert back conservative to primitive
              qp(i,j,k,1) = conservativeQ(1)
              qp(i,j,k,2) = conservativeQ(2) / conservativeQ(1)
              qp(i,j,k,3) = conservativeQ(3) / conservativeQ(1)
              qp(i,j,k,4) = conservativeQ(4) / conservativeQ(1)
              qp(i,j,k,5) = (gm-1.0) * ( conservativeQ(5) - (0.5 * sum(conservativeQ(2:4)**2) / conservativeQ(1)) )
            end do
          end do
        end do

    end subroutine update_laminar_variables


    function Flux(ql, qr, du, inputs)
      !< Calculate the total flux through face for laminar flow.
      !--------------------------------------
      ! calculate the total flux through face
      !---------------------------------------
      implicit none
      real, dimension(1:n_var), intent(in) :: ql !left state
      real, dimension(1:n_var), intent(in) :: qr !right state
      !conservative form of updated neighbour
      real, dimension(1:n_var), intent(in) :: du
      real, dimension(1:7)    , intent(in) :: inputs
      real, dimension(1:n_var)             :: Flux
      real, dimension(1:n_var)             :: U ! conservative variables
      real, dimension(1:n_var)             :: W ! new primitive variables
      real, dimension(1:n_var)             :: P ! primitive variables of right cell

      !for extraction of the inputs
      real :: area
      real :: nx
      real :: ny
      real :: nz
      real :: volume
      real :: mmu
      real :: tmu


      real    :: dudx
      real    :: dudy
      real    :: dudz
      real    :: dvdx
      real    :: dvdy
      real    :: dvdz
      real    :: dwdx
      real    :: dwdy
      real    :: dwdz
      real    :: dTdx
      real    :: dTdy
      real    :: dTdz
      real    :: T1, T2
      real    :: uface
      real    :: vface
      real    :: wface
      real    :: trace
      real    :: Tauxx
      real    :: Tauyy
      real    :: Tauzz
      real    :: Tauxy
      real    :: Tauxz
      real    :: Tauyz
      real    :: Qx
      real    :: Qy
      real    :: Qz
      real    :: HalfRhoUsquare
      real    :: RhoHt
      real    :: K_heat
      real    :: FaceNormalVelocity
      real    :: mu

      area   = inputs(1)
      nx     = inputs(2)
      ny     = inputs(3)
      nz     = inputs(4)
      volume = inputs(5)
      mmu    = inputs(6)
      tmu    = inputs(7)


      !save the old stat in P
      P = qr

      ! find conservative variable
      U(1)   =   ql(1)
      U(2)   =   ql(1) * ql(2)
      U(3)   =   ql(1) * ql(3)
      U(4)   =   ql(1) * ql(4)
      U(5)   = ( ql(5) / (gm-1.0) ) + ( 0.5 * ql(1) * sum(ql(2:4)**2) )

      U(1:5) = U(1:5) + du(1:5)
      

      W(1)   =   U(1)
      W(2)   =   U(2) / U(1)
      W(3)   =   U(3) / U(1)
      W(4)   =   U(4) / U(1)
      W(5)   = (gm-1.0) * ( U(5) - ( 0.5 * SUM(U(2:4)**2) / U(1) ) )

      FaceNormalVelocity = (W(2) * nx) + (W(3) * ny) + (W(4) * nz)
      uface = 0.5 * ( W(2) + P(2) )
      vface = 0.5 * ( W(3) + P(3) )
      wface = 0.5 * ( W(4) + P(4) )


      Flux(1) =   W(1) * FaceNormalVelocity
      Flux(2) = ( W(2) * Flux(1) ) + ( W(5) * nx )
      Flux(3) = ( W(3) * Flux(1) ) + ( W(5) * ny )
      Flux(4) = ( W(4) * Flux(1) ) + ( W(5) * nz )

      HalfRhoUsquare = 0.5 * W(1) * ( W(2)*W(2) + W(3)*W(3) + W(4)*W(4) )
      RhoHt          = ( (gm/(gm-1.0)) * W(5) ) + HalfRhoUsquare
      Flux(5)        = RhoHt * FaceNormalVelocity


      ! viscous terms
      mu = mmu + tmu
      T1     =    W(5) / ( W(1) * R_gas )
      T2     =    P(5) / ( P(1) * R_gas )
      dTdx   =  ( T2   - T1   ) * nx * Area / Volume
      dTdy   =  ( T2   - T1   ) * ny * Area / Volume
      dTdz   =  ( T2   - T1   ) * nz * Area / Volume
      dudx   =  ( P(2) - W(2) ) * nx * Area / Volume
      dudy   =  ( P(2) - W(2) ) * ny * Area / Volume
      dudz   =  ( P(2) - W(2) ) * nz * Area / Volume
      dvdx   =  ( P(3) - W(3) ) * nx * Area / Volume
      dvdy   =  ( P(3) - W(3) ) * ny * Area / Volume
      dvdz   =  ( P(3) - W(3) ) * nz * Area / Volume
      dwdx   =  ( P(4) - W(4) ) * nx * Area / Volume
      dwdy   =  ( P(4) - W(4) ) * ny * Area / Volume
      dwdz   =  ( P(4) - W(4) ) * nz * Area / Volume

      trace = dudx + dvdy + dwdz
      Tauxx =  2. * mu * (dudx - trace/3.0)
      Tauyy =  2. * mu * (dvdy - trace/3.0)
      Tauzz =  2. * mu * (dwdz - trace/3.0)
      Tauxy = mu * (dvdx + dudy)
      Tauxz = mu * (dwdx + dudz)
      Tauyz = mu * (dwdy + dvdz)

      K_heat = ( mmu / Pr  + tmu/tpr) * gm * R_gas / ( gm - 1.0 )
      Qx = K_heat*dTdx
      Qy = K_heat*dTdy
      Qz = K_heat*dTdz

      Flux(2) = Flux(2) - ( Tauxx * nx + Tauxy * ny + Tauxz * nz )
      Flux(3) = Flux(3) - ( Tauxy * nx + Tauyy * ny + Tauyz * nz )
      Flux(4) = Flux(4) - ( Tauxz * nx + Tauyz * ny + Tauzz * nz )
      Flux(5) = Flux(5) - ( Tauxx * uface + Tauxy * vface + Tauxz * wface + Qx ) * nx
      Flux(5) = Flux(5) - ( Tauxy * uface + Tauyy * vface + Tauyz * wface + Qy ) * ny
      Flux(5) = Flux(5) - ( Tauxz * uface + Tauyz * vface + Tauzz * wface + Qz ) * nz

      Flux    = Flux * Area

    end function Flux 

    function SpectralRadius(ql, qr, inputs, c1, c2,eps)
      !< Calculated spectral radius
      implicit none
      real, dimension(1:n_var), intent(in) :: ql
      real, dimension(1:n_var), intent(in) :: qr
      real, dimension(1:7)    , intent(in) :: inputs
      real, dimension(1:3)    , intent(in) :: c1
      real, dimension(1:3)    , intent(in) :: c2
      real                    , intent(in) :: eps

      ! local variables
      real                                 :: SpectralRadius
      real                                 :: NormalSpeed
      real                                 :: SpeedOfSound
      real                                 :: vis
      real                                 :: mu
      real                                 :: rho
      real                                 :: distance

      !extract inputs
      real :: Area
      real :: nx
      real :: ny
      real :: nz
      real :: volume
      real :: mm
      real :: tm

      Area = inputs(1)
      nx   = inputs(2)
      ny   = inputs(3)
      nz   = inputs(4)
      volume = inputs(5)
      mm   = inputs(6)
      tm   = inputs(7)
      ! in state vector q (2-4) are the cell center velocity
      NormalSpeed = 0.5 * ( ( ( ql(2) + qr(2) ) * nx ) &
                          + ( ( ql(3) + qr(3) ) * ny ) &
                          + ( ( ql(4) + qr(4) ) * nz ) &
                          )
      NormalSpeed = abs(NormalSpeed)

      SpeedOfSound = 0.5*( sqrt(gm*ql(5)/ql(1)) + sqrt(gm*qr(5)/qr(1)) )

      ! visocus part
      mu  = mm/Pr + tm/tPr
      rho = 0.5*( ql(1) + qr(1) )
      distance = sqrt((c1(1)-c2(1))**2 + (c1(2)-c2(2))**2 +(c1(3)-c2(3))**2) 
      vis = 2.0* gm * (mm/pr + tm/tpr) / ( rho * distance )
      SpectralRadius = ( 0.5*((1.0+eps)*NormalSpeed &
                      + sqrt(((eps-1.0)**2)*(NormalSpeed**2) &
                      + 4*eps*(SpeedOfSound**2))) + vis) * Area

    end function SpectralRadius


    subroutine update_SST_variables(qp, residue, cells, Ifaces, Jfaces, Kfaces, dims)
      !< Update the RANS (SST) equation with LU-SGS
      implicit none
      type(extent), intent(in) :: dims
      real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(inout) :: qp
      real, dimension(:, :, :, :), intent(in)  :: residue
      !< Store residue at each cell-center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      integer :: i,j,k
        real, dimension(1:7)     :: deltaU
        real, dimension(1:7)     :: D
        real, dimension(1:7)     :: conservativeQ
        real, dimension(1:7)     :: OldIminusFlux
        real, dimension(1:7)     :: OldJminusFlux
        real, dimension(1:7)     :: OldKminusFlux
        real, dimension(1:7)     :: NewIminusFlux
        real, dimension(1:7)     :: NewJminusFlux
        real, dimension(1:7)     :: NewKminusFlux
        real, dimension(1:7)     :: DelIminusFlux
        real, dimension(1:7)     :: DelJminusFlux
        real, dimension(1:7)     :: DelKminusFlux
        real, dimension(1:6)     :: LambdaTimesArea
        real, dimension(1:7)     :: Q0 ! state at cell
        real, dimension(1:7)     :: Q1 ! state at neighbours 
        real, dimension(1:7)     :: Q2
        real, dimension(1:7)     :: Q3
        real, dimension(1:7)     :: Q4
        real, dimension(1:7)     :: Q5
        real, dimension(1:7)     :: Q6
        real, dimension(1:7)     :: DQ0! change in state
        real, dimension(1:7)     :: DQ1
        real, dimension(1:7)     :: DQ2
        real, dimension(1:7)     :: DQ3
        real, dimension(1:7)     :: DQ4
        real, dimension(1:7)     :: DQ5
        real, dimension(1:7)     :: DQ6
        real, dimension(1:8)     :: Flist1
        real, dimension(1:8)     :: Flist2
        real, dimension(1:8)     :: Flist3
        real, dimension(1:8)     :: Flist4
        real, dimension(1:8)     :: Flist5
        real, dimension(1:8)     :: Flist6
        real, dimension(1:3)     :: C0
        real, dimension(1:3)     :: C1
        real, dimension(1:3)     :: C2
        real, dimension(1:3)     :: C3
        real, dimension(1:3)     :: C4
        real, dimension(1:3)     :: C5
        real, dimension(1:3)     :: C6
        real                     :: beta
        real                     :: eps
        real                     :: M
        real                     :: VMag
        real                     :: SoundMag
        real                     :: u,v,w,r,p,kk,ww,H
        real                     :: factor
        real, dimension(1:7,1:7) :: PrecondInv

        ! intermittency
        real :: De, Dp

        De = 0.0
        Dp = 0.0


        !intialize delQ
        delQstar = 0.0

        !forward sweep
        do k=1,dims%kmx-1
          do j=1,dims%jmx-1
            do i=1,dims%imx-1
              C0  = (/Cells(i  ,j  ,k  )%Centerx,Cells(i  ,j  ,k  )%Centery,Cells(i  ,j  ,k  )%Centerz/)
              C1  = (/Cells(i-1,j  ,k  )%Centerx,Cells(i-1,j  ,k  )%Centery,Cells(i-1,j  ,k  )%Centerz/)
              C2  = (/Cells(i  ,j-1,k  )%Centerx,Cells(i  ,j-1,k  )%Centery,Cells(i  ,j-1,k  )%Centerz/)
              C3  = (/Cells(i  ,j  ,k-1)%Centerx,Cells(i  ,j  ,k-1)%Centery,Cells(i  ,j  ,k-1)%Centerz/)
              C4  = (/Cells(i+1,j  ,k  )%Centerx,Cells(i+1,j  ,k  )%Centery,Cells(i+1,j  ,k  )%Centerz/)
              C5  = (/Cells(i  ,j+1,k  )%Centerx,Cells(i  ,j+1,k  )%Centery,Cells(i  ,j+1,k  )%Centerz/)
              C6  = (/Cells(i  ,j  ,k+1)%Centerx,Cells(i  ,j  ,k+1)%Centery,Cells(i  ,j  ,k+1)%Centerz/)

              Q0  = qp(i  , j  , k  , 1:7)
              Q1  = qp(i-1, j  , k  , 1:7)
              Q2  = qp(i  , j-1, k  , 1:7)
              Q3  = qp(i  , j  , k-1, 1:7)
              Q4  = qp(i+1, j  , k  , 1:7)
              Q5  = qp(i  , j+1, k  , 1:7)
              Q6  = qp(i  , j  , k+1, 1:7)

              DQ0 = 0.0
              DQ1 = delQstar(i-1, j  , k  , 1:7)
              DQ2 = delQstar(i  , j-1, k  , 1:7)
              DQ3 = delQstar(i  , j  , k-1, 1:7)

              Flist1(1) =  Ifaces(i,j,k)%A
              Flist1(2) = -Ifaces(i,j,k)%nx
              Flist1(3) = -Ifaces(i,j,k)%ny
              Flist1(4) = -Ifaces(i,j,k)%nz
              Flist1(5) = 0.5*(cells(i-1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))
              Flist1(8) = 0.5*(sst_F1(i-1, j  , k  ) + sst_F1(i,j,k))

              Flist2(1) =  Jfaces(i,j,k)%A
              Flist2(2) = -Jfaces(i,j,k)%nx
              Flist2(3) = -Jfaces(i,j,k)%ny
              Flist2(4) = -Jfaces(i,j,k)%nz
              Flist2(5) = 0.5*(cells(i  , j-1, k  )%volume + cells(i,j,k)%volume)
              Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))
              Flist2(8) = 0.5*(sst_F1(i  , j-1, k  ) + sst_F1(i,j,k))

              Flist3(1) =  Kfaces(i,j,k)%A
              Flist3(2) = -Kfaces(i,j,k)%nx
              Flist3(3) = -Kfaces(i,j,k)%ny
              Flist3(4) = -Kfaces(i,j,k)%nz
              Flist3(5) = 0.5*(cells(i  , j  , k-1)%volume + cells(i,j,k)%volume)
              Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))
              Flist3(8) = 0.5*(sst_F1(i  , j  , k-1) + sst_F1(i,j,k))

              Flist4(1) =  Ifaces(i+1,j,k)%A
              Flist4(2) = +Ifaces(i+1,j,k)%nx
              Flist4(3) = +Ifaces(i+1,j,k)%ny
              Flist4(4) = +Ifaces(i+1,j,k)%nz
              Flist4(5) = 0.5*(cells(i+1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))
              Flist4(8) = 0.5*(sst_F1(i+1, j  , k  ) + sst_F1(i,j,k))

              Flist5(1) =  Jfaces(i,j+1,k)%A
              Flist5(2) = +Jfaces(i,j+1,k)%nx
              Flist5(3) = +Jfaces(i,j+1,k)%ny
              Flist5(4) = +Jfaces(i,j+1,k)%nz
              Flist5(5) = 0.5*(cells(i  , j+1, k  )%volume + cells(i,j,k)%volume)
              Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))
              Flist5(8) = 0.5*(sst_F1(i  , j+1, k  ) + sst_F1(i,j,k))

              Flist6(1) =  Kfaces(i,j,k+1)%A
              Flist6(2) = +Kfaces(i,j,k+1)%nx
              Flist6(3) = +Kfaces(i,j,k+1)%ny
              Flist6(4) = +Kfaces(i,j,k+1)%nz
              Flist6(5) = 0.5*(cells(i  , j  , k+1)%volume + cells(i,j,k)%volume)
              Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))
              Flist6(8) = 0.5*(sst_F1(i  , j  , k+1) + sst_F1(i,j,k))

              NewIminusFlux     = SSTFlux(Q1, Q0, DQ1, Flist1)
              NewJminusFlux     = SSTFlux(Q2, Q0, DQ2, Flist2)
              NewKminusFlux     = SSTFlux(Q3, Q0, DQ3, Flist3)
              OldIminusFlux     = SSTFlux(Q1, Q0, DQ0, Flist1)
              OldJminusFlux     = SSTFlux(Q2, Q0, DQ0, Flist2)
              OldKminusFlux     = SSTFlux(Q3, Q0, DQ0, Flist3)

              !---preconditioning---
              r  = Q0(1)
              u  = Q0(2)
              v  = Q0(3)
              w  = Q0(4)
              p  = Q0(5)
              kk = Q0(6)
              ww = Q0(7)
              VMag     = sqrt(u*u + v*v + w*w)
              SoundMag = sqrt(gm*p/r)
              M        = VMag/SoundMag 
              H  = (gm*p/(r*(gm-1.0))) + 0.5*(VMag)
              eps = min(1.0, max(M*M, Minf*Minf))
              factor = (1.0-eps)*(gm-1.0)/(SoundMag*SoundMag)

              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0, eps)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0, eps)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0, eps)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0, eps)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0, eps)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0, eps)


              PrecondInv(1,1) = 1.0 - factor*1*VMag*VMag/2.0
              PrecondInv(2,1) = 0.0 - factor*u*VMag*VMag/2.0
              PrecondInv(3,1) = 0.0 - factor*v*VMag*VMag/2.0
              PrecondInv(4,1) = 0.0 - factor*w*VMag*VMag/2.0
              PrecondInv(5,1) = 0.0 - factor*H*VMag*VMag/2.0
              PrecondInv(6,1) = 0.0 - factor*kk*VMag*VMag/2.0
              PrecondInv(7,1) = 0.0 - factor*ww*VMag*VMag/2.0
              PrecondInv(1,2) = 0.0 - factor*1*(-u)
              PrecondInv(2,2) = 1.0 - factor*u*(-u)
              PrecondInv(3,2) = 0.0 - factor*v*(-u)
              PrecondInv(4,2) = 0.0 - factor*w*(-u)
              PrecondInv(5,2) = 0.0 - factor*H*(-u)
              PrecondInv(6,2) = 0.0 - factor*kk*(-u)
              PrecondInv(7,2) = 0.0 - factor*ww*(-u)
              PrecondInv(1,3) = 0.0 - factor*1*(-v)
              PrecondInv(2,3) = 0.0 - factor*u*(-v)
              PrecondInv(3,3) = 1.0 - factor*v*(-v)
              PrecondInv(4,3) = 0.0 - factor*w*(-v)
              PrecondInv(5,3) = 0.0 - factor*H*(-v)
              PrecondInv(6,3) = 0.0 - factor*kk*(-v)
              PrecondInv(7,3) = 0.0 - factor*ww*(-v)
              PrecondInv(1,4) = 0.0 - factor*1*(-w)
              PrecondInv(2,4) = 0.0 - factor*u*(-w)
              PrecondInv(3,4) = 0.0 - factor*v*(-w)
              PrecondInv(4,4) = 1.0 - factor*w*(-w)
              PrecondInv(5,4) = 0.0 - factor*H*(-w)
              PrecondInv(6,4) = 0.0 - factor*kk*(-w)
              PrecondInv(7,4) = 0.0 - factor*ww*(-w)
              PrecondInv(1,5) = 0.0 - factor*1*(1.)
              PrecondInv(2,5) = 0.0 - factor*u*(1.)
              PrecondInv(3,5) = 0.0 - factor*v*(1.)
              PrecondInv(4,5) = 0.0 - factor*w*(1.)
              PrecondInv(5,5) = 1.0 - factor*H*(1.)
              PrecondInv(6,5) = 0.0 - factor*kk*(1.)
              PrecondInv(7,5) = 0.0 - factor*ww*(1.)
              PrecondInv(1,6) = 0.0 - factor*1*(-1.)
              PrecondInv(2,6) = 0.0 - factor*u*(-1.)
              PrecondInv(3,6) = 0.0 - factor*v*(-1.)
              PrecondInv(4,6) = 0.0 - factor*w*(-1.)
              PrecondInv(5,6) = 0.0 - factor*H*(-1.)
              PrecondInv(6,6) = 1.0 - factor*kk*(-1.)
              PrecondInv(7,6) = 0.0 - factor*ww*(-1.)
              PrecondInv(1,7) = 0.0 - factor*1*(0.)
              PrecondInv(2,7) = 0.0 - factor*u*(0.)
              PrecondInv(3,7) = 0.0 - factor*v*(0.)
              PrecondInv(4,7) = 0.0 - factor*w*(0.)
              PrecondInv(5,7) = 0.0 - factor*H*(0.)
              PrecondInv(6,7) = 0.0 - factor*kk*(0.)
              PrecondInv(7,7) = 1.0 - factor*ww*(0.)
              !---end preconditioning


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux


              D = (cells(i,j,k)%volume/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
              beta = sst_F1(i,j,k)*beta1 + (1.0-sst_F1(i,j,k))*beta2
              D(6) = (D(6) + (bstar*qp(i,j,k,7))*cells(i,j,k)%volume)
              D(7) = (D(7) + 2.0*beta*qp(i,j,k,7)*cells(i,j,k)%volume)
              !storing D in Iflux array for backward sweep
              !F_p(i,j,k,1) = D
              

              deltaU(1:7) = -matmul(PrecondInv,residue(i,j,k,1:7)) &
                - 0.5*((matmul(PrecondInv,DelIminusFlux) - LambdaTimesArea(1)*delQstar(i-1,j,k,1:7)) &
                     + (matmul(PrecondInv,DelJminusFlux) - LambdaTimesArea(2)*delQstar(i,j-1,k,1:7)) &
                     + (matmul(PrecondInv,DelKminusFlux) - LambdaTimesArea(3)*delQstar(i,j,k-1,1:7)) )

              delQstar(i,j,k,1:7) = deltaU(1:7)/D
            end do
          end do
        end do

        delQ=0.0
        !backward sweep
            do i=dims%imx-1,1,-1
          do j=dims%jmx-1,1,-1
        do k=dims%kmx-1,1,-1
              C0  = (/Cells(i  ,j  ,k  )%Centerx,Cells(i  ,j  ,k  )%Centery,Cells(i  ,j  ,k  )%Centerz/)
              C1  = (/Cells(i-1,j  ,k  )%Centerx,Cells(i-1,j  ,k  )%Centery,Cells(i-1,j  ,k  )%Centerz/)
              C2  = (/Cells(i  ,j-1,k  )%Centerx,Cells(i  ,j-1,k  )%Centery,Cells(i  ,j-1,k  )%Centerz/)
              C3  = (/Cells(i  ,j  ,k-1)%Centerx,Cells(i  ,j  ,k-1)%Centery,Cells(i  ,j  ,k-1)%Centerz/)
              C4  = (/Cells(i+1,j  ,k  )%Centerx,Cells(i+1,j  ,k  )%Centery,Cells(i+1,j  ,k  )%Centerz/)
              C5  = (/Cells(i  ,j+1,k  )%Centerx,Cells(i  ,j+1,k  )%Centery,Cells(i  ,j+1,k  )%Centerz/)
              C6  = (/Cells(i  ,j  ,k+1)%Centerx,Cells(i  ,j  ,k+1)%Centery,Cells(i  ,j  ,k+1)%Centerz/)

              Q0  = qp(i  , j  , k  , 1:7)
              Q1  = qp(i-1, j  , k  , 1:7)
              Q2  = qp(i  , j-1, k  , 1:7)
              Q3  = qp(i  , j  , k-1, 1:7)
              Q4  = qp(i+1, j  , k  , 1:7)
              Q5  = qp(i  , j+1, k  , 1:7)
              Q6  = qp(i  , j  , k+1, 1:7)

              DQ0 = 0.0
              DQ4 = delQ(i+1, j  , k  , 1:7)
              DQ5 = delQ(i  , j+1, k  , 1:7)
              DQ6 = delQ(i  , j  , k+1, 1:7)

              Flist1(1) =  Ifaces(i,j,k)%A
              Flist1(2) = -Ifaces(i,j,k)%nx
              Flist1(3) = -Ifaces(i,j,k)%ny
              Flist1(4) = -Ifaces(i,j,k)%nz
              Flist1(5) = 0.5*(cells(i-1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))
              Flist1(8) = 0.5*(sst_F1(i-1, j  , k  ) + sst_F1(i,j,k))

              Flist2(1) =  Jfaces(i,j,k)%A
              Flist2(2) = -Jfaces(i,j,k)%nx
              Flist2(3) = -Jfaces(i,j,k)%ny
              Flist2(4) = -Jfaces(i,j,k)%nz
              Flist2(5) = 0.5*(cells(i  , j-1, k  )%volume + cells(i,j,k)%volume)
              Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))
              Flist2(8) = 0.5*(sst_F1(i  , j-1, k  ) + sst_F1(i,j,k))

              Flist3(1) =  Kfaces(i,j,k)%A
              Flist3(2) = -Kfaces(i,j,k)%nx
              Flist3(3) = -Kfaces(i,j,k)%ny
              Flist3(4) = -Kfaces(i,j,k)%nz
              Flist3(5) = 0.5*(cells(i  , j  , k-1)%volume + cells(i,j,k)%volume)
              Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))
              Flist3(8) = 0.5*(sst_F1(i  , j  , k-1) + sst_F1(i,j,k))

              Flist4(1) =  Ifaces(i+1,j,k)%A
              Flist4(2) = +Ifaces(i+1,j,k)%nx
              Flist4(3) = +Ifaces(i+1,j,k)%ny
              Flist4(4) = +Ifaces(i+1,j,k)%nz
              Flist4(5) = 0.5*(cells(i+1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))
              Flist4(8) = 0.5*(sst_F1(i+1, j  , k  ) + sst_F1(i,j,k))

              Flist5(1) =  Jfaces(i,j+1,k)%A
              Flist5(2) = +Jfaces(i,j+1,k)%nx
              Flist5(3) = +Jfaces(i,j+1,k)%ny
              Flist5(4) = +Jfaces(i,j+1,k)%nz
              Flist5(5) = 0.5*(cells(i  , j+1, k  )%volume + cells(i,j,k)%volume)
              Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))
              Flist5(8) = 0.5*(sst_F1(i  , j+1, k  ) + sst_F1(i,j,k))

              Flist6(1) =  Kfaces(i,j,k+1)%A
              Flist6(2) = +Kfaces(i,j,k+1)%nx
              Flist6(3) = +Kfaces(i,j,k+1)%ny
              Flist6(4) = +Kfaces(i,j,k+1)%nz
              Flist6(5) = 0.5*(cells(i  , j  , k+1)%volume + cells(i,j,k)%volume)
              Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))
              Flist6(8) = 0.5*(sst_F1(i  , j  , k+1) + sst_F1(i,j,k))

              NewIminusFlux     = SSTFlux(Q4, Q0, DQ4, Flist4)
              NewJminusFlux     = SSTFlux(Q5, Q0, DQ5, Flist5)
              NewKminusFlux     = SSTFlux(Q6, Q0, DQ6, Flist6)
              OldIminusFlux     = SSTFlux(Q4, Q0, DQ0, Flist4)
              OldJminusFlux     = SSTFlux(Q5, Q0, DQ0, Flist5)
              OldKminusFlux     = SSTFlux(Q6, Q0, DQ0, Flist6)

              !---preconditioning---
              r  = Q0(1)
              u  = Q0(2)
              v  = Q0(3)
              w  = Q0(4)
              p  = Q0(5)
              kk = Q0(6)
              ww = Q0(7)
              VMag     = sqrt(u*u + v*v + w*w)
              SoundMag = sqrt(gm*p/r)
              M        = VMag/SoundMag 
              H  = (gm*p/(r*(gm-1.0))) + 0.5*(VMag)
              eps = min(1.0, max(M*M, Minf*Minf))
              factor = (1.0-eps)*(gm-1.0)/(SoundMag*SoundMag)

              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0, eps)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0, eps)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0, eps)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0, eps)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0, eps)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0, eps)


              PrecondInv(1,1) = 1.0 - factor*1*VMag*VMag/2.0
              PrecondInv(2,1) = 0.0 - factor*u*VMag*VMag/2.0
              PrecondInv(3,1) = 0.0 - factor*v*VMag*VMag/2.0
              PrecondInv(4,1) = 0.0 - factor*w*VMag*VMag/2.0
              PrecondInv(5,1) = 0.0 - factor*H*VMag*VMag/2.0
              PrecondInv(6,1) = 0.0 - factor*kk*VMag*VMag/2.0
              PrecondInv(7,1) = 0.0 - factor*ww*VMag*VMag/2.0
              PrecondInv(1,2) = 0.0 - factor*1*(-u)
              PrecondInv(2,2) = 1.0 - factor*u*(-u)
              PrecondInv(3,2) = 0.0 - factor*v*(-u)
              PrecondInv(4,2) = 0.0 - factor*w*(-u)
              PrecondInv(5,2) = 0.0 - factor*H*(-u)
              PrecondInv(6,2) = 0.0 - factor*kk*(-u)
              PrecondInv(7,2) = 0.0 - factor*ww*(-u)
              PrecondInv(1,3) = 0.0 - factor*1*(-v)
              PrecondInv(2,3) = 0.0 - factor*u*(-v)
              PrecondInv(3,3) = 1.0 - factor*v*(-v)
              PrecondInv(4,3) = 0.0 - factor*w*(-v)
              PrecondInv(5,3) = 0.0 - factor*H*(-v)
              PrecondInv(6,3) = 0.0 - factor*kk*(-v)
              PrecondInv(7,3) = 0.0 - factor*ww*(-v)
              PrecondInv(1,4) = 0.0 - factor*1*(-w)
              PrecondInv(2,4) = 0.0 - factor*u*(-w)
              PrecondInv(3,4) = 0.0 - factor*v*(-w)
              PrecondInv(4,4) = 1.0 - factor*w*(-w)
              PrecondInv(5,4) = 0.0 - factor*H*(-w)
              PrecondInv(6,4) = 0.0 - factor*kk*(-w)
              PrecondInv(7,4) = 0.0 - factor*ww*(-w)
              PrecondInv(1,5) = 0.0 - factor*1*(1.)
              PrecondInv(2,5) = 0.0 - factor*u*(1.)
              PrecondInv(3,5) = 0.0 - factor*v*(1.)
              PrecondInv(4,5) = 0.0 - factor*w*(1.)
              PrecondInv(5,5) = 1.0 - factor*H*(1.)
              PrecondInv(6,5) = 0.0 - factor*kk*(1.)
              PrecondInv(7,5) = 0.0 - factor*ww*(1.)
              PrecondInv(1,6) = 0.0 - factor*1*(-1.)
              PrecondInv(2,6) = 0.0 - factor*u*(-1.)
              PrecondInv(3,6) = 0.0 - factor*v*(-1.)
              PrecondInv(4,6) = 0.0 - factor*w*(-1.)
              PrecondInv(5,6) = 0.0 - factor*H*(-1.)
              PrecondInv(6,6) = 1.0 - factor*kk*(-1.)
              PrecondInv(7,6) = 0.0 - factor*ww*(-1.)
              PrecondInv(1,7) = 0.0 - factor*1*(0.)
              PrecondInv(2,7) = 0.0 - factor*u*(0.)
              PrecondInv(3,7) = 0.0 - factor*v*(0.)
              PrecondInv(4,7) = 0.0 - factor*w*(0.)
              PrecondInv(5,7) = 0.0 - factor*H*(0.)
              PrecondInv(6,7) = 0.0 - factor*kk*(0.)
              PrecondInv(7,7) = 1.0 - factor*ww*(0.)
              !---end preconditioning


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux

              D = (cells(i,j,k)%volume/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
              beta = sst_F1(i,j,k)*beta1 + (1.0-sst_F1(i,j,k))*beta2
              D(6) = (D(6) + (bstar*qp(i,j,k,7))*cells(i,j,k)%volume)
              D(7) = (D(7) + 2.0*beta*qp(i,j,k,7)*cells(i,j,k)%volume)


              delQ(i,j,k,1:7) = delQstar(i,j,k,1:7) &
                - 0.5*((matmul(PrecondInv,DelIminusFlux) - LambdaTimesArea(4)*delQ(i+1,j,k,1:7)) &
                     + (matmul(PrecondInv,DelJminusFlux) - LambdaTimesArea(5)*delQ(i,j+1,k,1:7)) &
                     + (matmul(PrecondInv,DelKminusFlux) - LambdaTimesArea(6)*delQ(i,j,k+1,1:7)) )/D

            end do
          end do
        end do

        
        do k=1,dims%kmx-1
          do j = 1,dims%jmx-1
            do i = 1,dims%imx-1
              conservativeQ(1) = qp(i,j,k,1)
              conservativeQ(2) = qp(i,j,k,1) * qp(i,j,k,2)
              conservativeQ(3) = qp(i,j,k,1) * qp(i,j,k,3)
              conservativeQ(4) = qp(i,j,k,1) * qp(i,j,k,4)
              conservativeQ(5) = (qp(i,j,k,5) / (gm-1.0)) + ( 0.5 * qp(i,j,k,1) * sum( qp(i,j,k,2:4)**2) )
              conservativeQ(6) = qp(i,j,k,1) * qp(i,j,k,6)
              conservativeQ(7) = qp(i,j,k,1) * qp(i,j,k,7)
              
              ! add new change into conservative solution
              conservativeQ(1:7) = conservativeQ(1:7) + delQ(i,j,k,1:7)

              ! convert back conservative to primitive
              qp(i,j,k,1) = conservativeQ(1)
              qp(i,j,k,2) = conservativeQ(2) / conservativeQ(1)
              qp(i,j,k,3) = conservativeQ(3) / conservativeQ(1)
              qp(i,j,k,4) = conservativeQ(4) / conservativeQ(1)
              qp(i,j,k,5) = (gm-1.0) * ( conservativeQ(5) - (0.5 * sum(conservativeQ(2:4)**2) / conservativeQ(1)) )
              if(conservativeQ(6)>0)then
              qp(i,j,k,6) = conservativeQ(6) / conservativeQ(1)
              end if
              if(conservativeQ(7)>0)then
              qp(i,j,k,7) = conservativeQ(7) / conservativeQ(1)
              end if
            end do
          end do
        end do

    end subroutine update_SST_variables


    function SSTFlux(ql, qr, du, inputs)
      !< Calculate the total flux through face for turbulent flow (SST)
      !--------------------------------------
      ! calculate the total flux through face
      !---------------------------------------
      implicit none
      real, dimension(1:n_var), intent(in) :: ql !left state
      real, dimension(1:n_var), intent(in) :: qr !right state
      !conservative form of updated neighbour
      real, dimension(1:n_var), intent(in) :: du
      real, dimension(1:8)    , intent(in) :: inputs
      real, dimension(1:n_var)             :: Flux
      real, dimension(1:n_var)             :: SSTFlux
      real, dimension(1:n_var)             :: U ! conservative variables
      real, dimension(1:n_var)             :: W ! new primitive variables
      real, dimension(1:n_var)             :: P ! primitive variables of right cell

      !for extraction of the inputs
      real :: area
      real :: nx
      real :: ny
      real :: nz
      real :: volume
      real :: mmu
      real :: tmu


      real    :: dudx
      real    :: dudy
      real    :: dudz
      real    :: dvdx
      real    :: dvdy
      real    :: dvdz
      real    :: dwdx
      real    :: dwdy
      real    :: dwdz
      real    :: dTdx
      real    :: dTdy
      real    :: dTdz
      real    :: dtkdx
      real    :: dtkdy
      real    :: dtkdz
      real    :: dtwdx
      real    :: dtwdy
      real    :: dtwdz
      real    :: T1, T2
      real    :: uface
      real    :: vface
      real    :: wface
      real    :: trace
      real    :: Tauxx
      real    :: Tauyy
      real    :: Tauzz
      real    :: Tauxy
      real    :: Tauxz
      real    :: Tauyz
      real    :: Qx
      real    :: Qy
      real    :: Qz
      real    :: HalfRhoUsquare
      real    :: RhoHt
      real    :: K_heat
      real    :: FaceNormalVelocity
      real    :: mu
      real    :: sigma_k
      real    :: sigma_w
      real    :: F1

      area   = inputs(1)
      nx     = inputs(2)
      ny     = inputs(3)
      nz     = inputs(4)
      volume = inputs(5)
      mmu    = inputs(6)
      tmu    = inputs(7)
      F1     = inputs(8)


      !save the old stat in P
      P = qr

      ! find conservative variable
      U(1)   =   ql(1)
      U(2)   =   ql(1) * ql(2)
      U(3)   =   ql(1) * ql(3)
      U(4)   =   ql(1) * ql(4)
      U(5)   = ( ql(5) / (gm-1.0) ) + ( 0.5 * ql(1) * sum(ql(2:4)**2) )
      U(6)   =   ql(1) * ql(6)
      U(7)   =   ql(1) * ql(7)

      U(1:n_var) = U(1:n_var) + du(1:n_var)
      

      W(1)   =   U(1)
      W(2)   =   U(2) / U(1)
      W(3)   =   U(3) / U(1)
      W(4)   =   U(4) / U(1)
      W(5)   = (gm-1.0) * ( U(5) - ( 0.5 * SUM(U(2:4)**2) / U(1) ) )
      W(6)   =   U(6) / U(1)
      W(7)   =   U(7) / U(1)
      W(6)   = W(6) + 0.5*(1.-sign(1.,W(6)))*(ql(6)-W(6))
      W(7)   = W(7) + 0.5*(1.-sign(1.,W(7)))*(ql(7)-W(7))

      FaceNormalVelocity = (W(2) * nx) + (W(3) * ny) + (W(4) * nz)
      uface = 0.5 * ( W(2) + P(2) )
      vface = 0.5 * ( W(3) + P(3) )
      wface = 0.5 * ( W(4) + P(4) )


      Flux(1) =   W(1) * FaceNormalVelocity
      Flux(2) = ( W(2) * Flux(1) ) + ( W(5) * nx )
      Flux(3) = ( W(3) * Flux(1) ) + ( W(5) * ny )
      Flux(4) = ( W(4) * Flux(1) ) + ( W(5) * nz )

      HalfRhoUsquare = 0.5 * W(1) * ( W(2)*W(2) + W(3)*W(3) + W(4)*W(4) )
      RhoHt          = ( (gm/(gm-1.0)) * W(5) ) + HalfRhoUsquare
      Flux(5)        = RhoHt * FaceNormalVelocity
      Flux(6) = ( W(6) * Flux(1) )   
      Flux(7) = ( W(7) * Flux(1) )   


      ! viscous terms
      mu = mmu + tmu
      T1     =    W(5) / ( W(1) * R_gas )
      T2     =    P(5) / ( P(1) * R_gas )
      dTdx   =  ( T2   - T1   ) * nx * Area / Volume
      dTdy   =  ( T2   - T1   ) * ny * Area / Volume
      dTdz   =  ( T2   - T1   ) * nz * Area / Volume
      dudx   =  ( P(2) - W(2) ) * nx * Area / Volume
      dudy   =  ( P(2) - W(2) ) * ny * Area / Volume
      dudz   =  ( P(2) - W(2) ) * nz * Area / Volume
      dvdx   =  ( P(3) - W(3) ) * nx * Area / Volume
      dvdy   =  ( P(3) - W(3) ) * ny * Area / Volume
      dvdz   =  ( P(3) - W(3) ) * nz * Area / Volume
      dwdx   =  ( P(4) - W(4) ) * nx * Area / Volume
      dwdy   =  ( P(4) - W(4) ) * ny * Area / Volume
      dwdz   =  ( P(4) - W(4) ) * nz * Area / Volume
      dtkdx  =  ( P(6) - W(6) ) * nx * Area / Volume
      dtkdy  =  ( P(6) - W(6) ) * ny * Area / Volume
      dtkdz  =  ( P(6) - W(6) ) * nz * Area / Volume
      dtwdx  =  ( P(7) - W(7) ) * nx * Area / Volume
      dtwdy  =  ( P(7) - W(7) ) * ny * Area / Volume
      dtwdz  =  ( P(7) - W(7) ) * nz * Area / Volume

      trace = dudx + dvdy + dwdz
      Tauxx =  2. * mu * (dudx - trace/3.0)
      Tauyy =  2. * mu * (dvdy - trace/3.0)
      Tauzz =  2. * mu * (dwdz - trace/3.0)
      Tauxy = mu * (dvdx + dudy)
      Tauxz = mu * (dwdx + dudz)
      Tauyz = mu * (dwdy + dvdz)

      K_heat = ( mmu / Pr  + tmu/tpr) * gm * R_gas / ( gm - 1.0 )
      Qx = K_heat*dTdx
      Qy = K_heat*dTdy
      Qz = K_heat*dTdz

      sigma_k = sigma_k1*F1 + sigma_k2*(1.0 - F1)
      sigma_w = sigma_w1*F1 + sigma_w2*(1.0 - F1)

      Flux(2) = Flux(2) - ( Tauxx * nx + Tauxy * ny + Tauxz * nz )
      Flux(3) = Flux(3) - ( Tauxy * nx + Tauyy * ny + Tauyz * nz )
      Flux(4) = Flux(4) - ( Tauxz * nx + Tauyz * ny + Tauzz * nz )
      Flux(5) = Flux(5) - ( Tauxx * uface + Tauxy * vface + Tauxz * wface + Qx ) * nx
      Flux(5) = Flux(5) - ( Tauxy * uface + Tauyy * vface + Tauyz * wface + Qy ) * ny
      Flux(5) = Flux(5) - ( Tauxz * uface + Tauyz * vface + Tauzz * wface + Qz ) * nz
      Flux(6) = Flux(6) + (mmu + sigma_k*tmu)*(dtkdx*nx + dtkdy*ny + dtkdz*nz)
      Flux(7) = Flux(7) + (mmu + sigma_w*tmu)*(dtwdx*nx + dtwdy*ny + dtwdz*nz)

      Flux    = Flux * Area
      SSTFlux = Flux

    end function SSTFlux 

    subroutine update_SA_variables(qp, residue, cells, Ifaces, Jfaces, Kfaces, dims)
      !< Update the RANS (SA) equation with LU-SGS
      implicit none
      type(extent), intent(in) :: dims
      real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(inout) :: qp
      real, dimension(:, :, :, :), intent(in)  :: residue
      !< Store residue at each cell-center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      integer :: i,j,k
        real, dimension(1:6)     :: deltaU
        real, dimension(1:6)     :: D
        real, dimension(1:6)     :: conservativeQ
        real, dimension(1:6)     :: OldIminusFlux
        real, dimension(1:6)     :: OldJminusFlux
        real, dimension(1:6)     :: OldKminusFlux
        real, dimension(1:6)     :: NewIminusFlux
        real, dimension(1:6)     :: NewJminusFlux
        real, dimension(1:6)     :: NewKminusFlux
        real, dimension(1:6)     :: DelIminusFlux
        real, dimension(1:6)     :: DelJminusFlux
        real, dimension(1:6)     :: DelKminusFlux
        real, dimension(1:6)     :: LambdaTimesArea
        real, dimension(1:6)     :: Q0 ! state at cell
        real, dimension(1:6)     :: Q1 ! state at neighbours 
        real, dimension(1:6)     :: Q2
        real, dimension(1:6)     :: Q3
        real, dimension(1:6)     :: Q4
        real, dimension(1:6)     :: Q5
        real, dimension(1:6)     :: Q6
        real, dimension(1:6)     :: DQ0! change in state
        real, dimension(1:6)     :: DQ1
        real, dimension(1:6)     :: DQ2
        real, dimension(1:6)     :: DQ3
        real, dimension(1:6)     :: DQ4
        real, dimension(1:6)     :: DQ5
        real, dimension(1:6)     :: DQ6
        real, dimension(1:7)     :: Flist1
        real, dimension(1:7)     :: Flist2
        real, dimension(1:7)     :: Flist3
        real, dimension(1:7)     :: Flist4
        real, dimension(1:7)     :: Flist5
        real, dimension(1:7)     :: Flist6
        real, dimension(1:3)     :: C0
        real, dimension(1:3)     :: C1
        real, dimension(1:3)     :: C2
        real, dimension(1:3)     :: C3
        real, dimension(1:3)     :: C4
        real, dimension(1:3)     :: C5
        real, dimension(1:3)     :: C6
        real                     :: eps
        real                     :: M
        real                     :: VMag
        real                     :: SoundMag
        real                     :: u,v,w,p,H,tv!r
        real                     :: factor
        real, dimension(1:6,1:6) :: PrecondInv
      real :: fv1
      real :: fv2
      real :: fw
      real :: g
      real :: r
      real :: dist_i
      real :: dist_i_2
      real :: Ji
      real :: Ji_2
      real :: Ji_3
      real :: S
      real :: Omega
      real :: k2
      real :: inv_k2_d2
      real :: Shat
      real :: inv_Shat
      real :: nu
      real :: glim
      real :: g_6
      real :: dfv1
      real :: dfv2
      real :: dfw
      real :: dShat
      real :: dr
      real :: dg
      real :: density



        !intialize delQ
        delQstar = 0.0

        !forward sweep
        do k=1,dims%kmx-1
          do j=1,dims%jmx-1
            do i=1,dims%imx-1
              density = qp(i,j,k,1)
              C0  = (/Cells(i  ,j  ,k  )%Centerx,Cells(i  ,j  ,k  )%Centery,Cells(i  ,j  ,k  )%Centerz/)
              C1  = (/Cells(i-1,j  ,k  )%Centerx,Cells(i-1,j  ,k  )%Centery,Cells(i-1,j  ,k  )%Centerz/)
              C2  = (/Cells(i  ,j-1,k  )%Centerx,Cells(i  ,j-1,k  )%Centery,Cells(i  ,j-1,k  )%Centerz/)
              C3  = (/Cells(i  ,j  ,k-1)%Centerx,Cells(i  ,j  ,k-1)%Centery,Cells(i  ,j  ,k-1)%Centerz/)
              C4  = (/Cells(i+1,j  ,k  )%Centerx,Cells(i+1,j  ,k  )%Centery,Cells(i+1,j  ,k  )%Centerz/)
              C5  = (/Cells(i  ,j+1,k  )%Centerx,Cells(i  ,j+1,k  )%Centery,Cells(i  ,j+1,k  )%Centerz/)
              C6  = (/Cells(i  ,j  ,k+1)%Centerx,Cells(i  ,j  ,k+1)%Centery,Cells(i  ,j  ,k+1)%Centerz/)

              Q0  = qp(i  , j  , k  , 1:6)
              Q1  = qp(i-1, j  , k  , 1:6)
              Q2  = qp(i  , j-1, k  , 1:6)
              Q3  = qp(i  , j  , k-1, 1:6)
              Q4  = qp(i+1, j  , k  , 1:6)
              Q5  = qp(i  , j+1, k  , 1:6)
              Q6  = qp(i  , j  , k+1, 1:6)

              DQ0 = 0.0
              DQ1 = delQstar(i-1, j  , k  , 1:6)
              DQ2 = delQstar(i  , j-1, k  , 1:6)
              DQ3 = delQstar(i  , j  , k-1, 1:6)

              Flist1(1) =  Ifaces(i,j,k)%A
              Flist1(2) = -Ifaces(i,j,k)%nx
              Flist1(3) = -Ifaces(i,j,k)%ny
              Flist1(4) = -Ifaces(i,j,k)%nz
              Flist1(5) = 0.5*(cells(i-1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))

              Flist2(1) =  Jfaces(i,j,k)%A
              Flist2(2) = -Jfaces(i,j,k)%nx
              Flist2(3) = -Jfaces(i,j,k)%ny
              Flist2(4) = -Jfaces(i,j,k)%nz
              Flist2(5) = 0.5*(cells(i  , j-1, k  )%volume + cells(i,j,k)%volume)
              Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))

              Flist3(1) =  Kfaces(i,j,k)%A
              Flist3(2) = -Kfaces(i,j,k)%nx
              Flist3(3) = -Kfaces(i,j,k)%ny
              Flist3(4) = -Kfaces(i,j,k)%nz
              Flist3(5) = 0.5*(cells(i  , j  , k-1)%volume + cells(i,j,k)%volume)
              Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))

              Flist4(1) =  Ifaces(i+1,j,k)%A
              Flist4(2) = +Ifaces(i+1,j,k)%nx
              Flist4(3) = +Ifaces(i+1,j,k)%ny
              Flist4(4) = +Ifaces(i+1,j,k)%nz
              Flist4(5) = 0.5*(cells(i+1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))

              Flist5(1) =  Jfaces(i,j+1,k)%A
              Flist5(2) = +Jfaces(i,j+1,k)%nx
              Flist5(3) = +Jfaces(i,j+1,k)%ny
              Flist5(4) = +Jfaces(i,j+1,k)%nz
              Flist5(5) = 0.5*(cells(i  , j+1, k  )%volume + cells(i,j,k)%volume)
              Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))

              Flist6(1) =  Kfaces(i,j,k+1)%A
              Flist6(2) = +Kfaces(i,j,k+1)%nx
              Flist6(3) = +Kfaces(i,j,k+1)%ny
              Flist6(4) = +Kfaces(i,j,k+1)%nz
              Flist6(5) = 0.5*(cells(i  , j  , k+1)%volume + cells(i,j,k)%volume)
              Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))

              NewIminusFlux     = SAFlux(Q1, Q0, DQ1, Flist1)
              NewJminusFlux     = SAFlux(Q2, Q0, DQ2, Flist2)
              NewKminusFlux     = SAFlux(Q3, Q0, DQ3, Flist3)
              OldIminusFlux     = SAFlux(Q1, Q0, DQ0, Flist1)
              OldJminusFlux     = SAFlux(Q2, Q0, DQ0, Flist2)
              OldKminusFlux     = SAFlux(Q3, Q0, DQ0, Flist3)

              !---preconditioning---
              r  = Q0(1)
              u  = Q0(2)
              v  = Q0(3)
              w  = Q0(4)
              p  = Q0(5)
              tv = Q0(6)
              VMag     = sqrt(u*u + v*v + w*w)
              SoundMag = sqrt(gm*p/r)
              M        = VMag/SoundMag 
              H  = (gm*p/(r*(gm-1.0))) + 0.5*(VMag)
              eps = min(1.0, max(M*M, Minf*Minf))
              factor = (1.0-eps)*(gm-1.0)/(SoundMag*SoundMag)

              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0, eps)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0, eps)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0, eps)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0, eps)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0, eps)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0, eps)


              PrecondInv(1,1) = 1.0 - factor*1*VMag*VMag/2.0
              PrecondInv(2,1) = 0.0 - factor*u*VMag*VMag/2.0
              PrecondInv(3,1) = 0.0 - factor*v*VMag*VMag/2.0
              PrecondInv(4,1) = 0.0 - factor*w*VMag*VMag/2.0
              PrecondInv(5,1) = 0.0 - factor*H*VMag*VMag/2.0
              PrecondInv(6,1) = 0.0 - factor*tv*VMag*VMag/2.0
              PrecondInv(1,2) = 0.0 - factor*1*(-u)
              PrecondInv(2,2) = 1.0 - factor*u*(-u)
              PrecondInv(3,2) = 0.0 - factor*v*(-u)
              PrecondInv(4,2) = 0.0 - factor*w*(-u)
              PrecondInv(5,2) = 0.0 - factor*H*(-u)
              PrecondInv(6,2) = 0.0 - factor*tv*(-u)
              PrecondInv(1,3) = 0.0 - factor*1*(-v)
              PrecondInv(2,3) = 0.0 - factor*u*(-v)
              PrecondInv(3,3) = 1.0 - factor*v*(-v)
              PrecondInv(4,3) = 0.0 - factor*w*(-v)
              PrecondInv(5,3) = 0.0 - factor*H*(-v)
              PrecondInv(6,3) = 0.0 - factor*tv*(-v)
              PrecondInv(1,4) = 0.0 - factor*1*(-w)
              PrecondInv(2,4) = 0.0 - factor*u*(-w)
              PrecondInv(3,4) = 0.0 - factor*v*(-w)
              PrecondInv(4,4) = 1.0 - factor*w*(-w)
              PrecondInv(5,4) = 0.0 - factor*H*(-w)
              PrecondInv(6,4) = 0.0 - factor*tv*(-w)
              PrecondInv(1,5) = 0.0 - factor*1*(1.)
              PrecondInv(2,5) = 0.0 - factor*u*(1.)
              PrecondInv(3,5) = 0.0 - factor*v*(1.)
              PrecondInv(4,5) = 0.0 - factor*w*(1.)
              PrecondInv(5,5) = 1.0 - factor*H*(1.)
              PrecondInv(6,5) = 0.0 - factor*tv*(1.)
              PrecondInv(1,6) = 0.0 - factor*1*(-1.)
              PrecondInv(2,6) = 0.0 - factor*u*(-1.)
              PrecondInv(3,6) = 0.0 - factor*v*(-1.)
              PrecondInv(4,6) = 0.0 - factor*w*(-1.)
              PrecondInv(5,6) = 0.0 - factor*H*(-1.)
              PrecondInv(6,6) = 1.0 - factor*tv*(-1.)
              !---end preconditioning



              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux


              D = (cells(i,j,k)%volume/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
              !storing D in Iflux array for backward sweep
              !F_p(i,j,k,1) = D
              ! -- source term derivatives -- !
              Omega = sqrt( ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
                           + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
                           + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
                            )&
                       )
              dist_i = dist(i,j,k)
              dist_i_2 = dist_i*dist_i
              k2 = kappa_sa*kappa_sa
              nu   = mu(i,j,k)/density
              Ji   = Q0(6)/nu
              Ji_2 = Ji*Ji
              Ji_3 = Ji_2*ji


              ! ___ functions ___
              fv1  = (Ji_3)/((Ji_3) + (cv1_3))
              fv2  = 1.0 - Ji/(1.0 + (Ji*fv1))

              ! ___ Shear stress for production ___
              S = Omega
              inv_k2_d2 = 1.0/(k2*dist_i_2)
              Shat      = S + Q0(6)*fv2*inv_k2_d2
              Shat      = max(Shat, 1.0e-10)
              inv_Shat  = 1.0/Shat
              dfv1 = 3.0*Ji_2*cv1_3/(nu*(Ji_3+cv1_3)**2)
              dfv2 = -((1.0/nu) - Ji_2*dfv1)/((1.0+Ji*fv1)**2)
              dShat = (fv2+Q0(6)*dfv2)*inv_k2_d2

              D = D - cb1*(Q0(6)*dShat+Shat)*cells(i,j,k)%volume

              ! ___ Destruction term___ !
              r    = min(Q0(6)*inv_Shat*inv_k2_d2, 10.0)
              g    = r + cw2*((r**6) - r)
              g_6  = g**6
              glim = ((1.0+cw3_6)/(g_6+cw3_6))**(1.0/6.0)
              fw   = g*glim
              dr = (Shat-Q0(6)*dShat)*inv_Shat*inv_Shat*inv_k2_d2
              dg = dr*(1.0+cw2*(6.0*(r**5)-1.0))
              dfw= dg*glim*(1.0-g_6/(g_6+cw3_6))

              D = D+cw1*(dfw*Q0(6) + 2*fw)*Q0(6)/dist_i_2*cells(i,j,k)%volume
              ! --  end of source term -- !

              deltaU(1:6) = -matmul(PrecondInv,residue(i,j,k,1:6)) &
                - 0.5*((matmul(PrecondInv,DelIminusFlux) - LambdaTimesArea(1)*delQstar(i-1,j,k,1:6)) &
                     + (matmul(PrecondInv,DelJminusFlux) - LambdaTimesArea(2)*delQstar(i,j-1,k,1:6)) &
                     + (matmul(PrecondInv,DelKminusFlux) - LambdaTimesArea(3)*delQstar(i,j,k-1,1:6)) )

              delQstar(i,j,k,1:6) = deltaU(1:6)/D
            end do
          end do
        end do


        delQ=0.0
        !backward sweep
            do i=dims%imx-1,1,-1
          do j=dims%jmx-1,1,-1
        do k=dims%kmx-1,1,-1
              density = qp(i,j,k,1)
              C0  = (/Cells(i  ,j  ,k  )%Centerx,Cells(i  ,j  ,k  )%Centery,Cells(i  ,j  ,k  )%Centerz/)
              C1  = (/Cells(i-1,j  ,k  )%Centerx,Cells(i-1,j  ,k  )%Centery,Cells(i-1,j  ,k  )%Centerz/)
              C2  = (/Cells(i  ,j-1,k  )%Centerx,Cells(i  ,j-1,k  )%Centery,Cells(i  ,j-1,k  )%Centerz/)
              C3  = (/Cells(i  ,j  ,k-1)%Centerx,Cells(i  ,j  ,k-1)%Centery,Cells(i  ,j  ,k-1)%Centerz/)
              C4  = (/Cells(i+1,j  ,k  )%Centerx,Cells(i+1,j  ,k  )%Centery,Cells(i+1,j  ,k  )%Centerz/)
              C5  = (/Cells(i  ,j+1,k  )%Centerx,Cells(i  ,j+1,k  )%Centery,Cells(i  ,j+1,k  )%Centerz/)
              C6  = (/Cells(i  ,j  ,k+1)%Centerx,Cells(i  ,j  ,k+1)%Centery,Cells(i  ,j  ,k+1)%Centerz/)

              Q0  = qp(i  , j  , k  , 1:6)
              Q1  = qp(i-1, j  , k  , 1:6)
              Q2  = qp(i  , j-1, k  , 1:6)
              Q3  = qp(i  , j  , k-1, 1:6)
              Q4  = qp(i+1, j  , k  , 1:6)
              Q5  = qp(i  , j+1, k  , 1:6)
              Q6  = qp(i  , j  , k+1, 1:6)

              DQ0 = 0.0
              DQ4 = delQ(i+1, j  , k  , 1:6)
              DQ5 = delQ(i  , j+1, k  , 1:6)
              DQ6 = delQ(i  , j  , k+1, 1:6)

              Flist1(1) =  Ifaces(i,j,k)%A
              Flist1(2) = -Ifaces(i,j,k)%nx
              Flist1(3) = -Ifaces(i,j,k)%ny
              Flist1(4) = -Ifaces(i,j,k)%nz
              Flist1(5) = 0.5*(cells(i-1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))

              Flist2(1) =  Jfaces(i,j,k)%A
              Flist2(2) = -Jfaces(i,j,k)%nx
              Flist2(3) = -Jfaces(i,j,k)%ny
              Flist2(4) = -Jfaces(i,j,k)%nz
              Flist2(5) = 0.5*(cells(i  , j-1, k  )%volume + cells(i,j,k)%volume)
              Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))

              Flist3(1) =  Kfaces(i,j,k)%A
              Flist3(2) = -Kfaces(i,j,k)%nx
              Flist3(3) = -Kfaces(i,j,k)%ny
              Flist3(4) = -Kfaces(i,j,k)%nz
              Flist3(5) = 0.5*(cells(i  , j  , k-1)%volume + cells(i,j,k)%volume)
              Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))

              Flist4(1) =  Ifaces(i+1,j,k)%A
              Flist4(2) = +Ifaces(i+1,j,k)%nx
              Flist4(3) = +Ifaces(i+1,j,k)%ny
              Flist4(4) = +Ifaces(i+1,j,k)%nz
              Flist4(5) = 0.5*(cells(i+1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))

              Flist5(1) =  Jfaces(i,j+1,k)%A
              Flist5(2) = +Jfaces(i,j+1,k)%nx
              Flist5(3) = +Jfaces(i,j+1,k)%ny
              Flist5(4) = +Jfaces(i,j+1,k)%nz
              Flist5(5) = 0.5*(cells(i  , j+1, k  )%volume + cells(i,j,k)%volume)
              Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))

              Flist6(1) =  Kfaces(i,j,k+1)%A
              Flist6(2) = +Kfaces(i,j,k+1)%nx
              Flist6(3) = +Kfaces(i,j,k+1)%ny
              Flist6(4) = +Kfaces(i,j,k+1)%nz
              Flist6(5) = 0.5*(cells(i  , j  , k+1)%volume + cells(i,j,k)%volume)
              Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))

              NewIminusFlux     = SAFlux(Q4, Q0, DQ4, Flist4)
              NewJminusFlux     = SAFlux(Q5, Q0, DQ5, Flist5)
              NewKminusFlux     = SAFlux(Q6, Q0, DQ6, Flist6)
              OldIminusFlux     = SAFlux(Q4, Q0, DQ0, Flist4)
              OldJminusFlux     = SAFlux(Q5, Q0, DQ0, Flist5)
              OldKminusFlux     = SAFlux(Q6, Q0, DQ0, Flist6)

              !---preconditioning---
              r  = Q0(1)
              u  = Q0(2)
              v  = Q0(3)
              w  = Q0(4)
              p  = Q0(5)
              tv = Q0(6)
              VMag     = sqrt(u*u + v*v + w*w)
              SoundMag = sqrt(gm*p/r)
              M        = VMag/SoundMag 
              H  = (gm*p/(r*(gm-1.0))) + 0.5*(VMag)
              eps = min(1.0, max(M*M, Minf*Minf))
              factor = (1.0-eps)*(gm-1.0)/(SoundMag*SoundMag)

              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0, eps)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0, eps)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0, eps)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0, eps)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0, eps)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0, eps)


              PrecondInv(1,1) = 1.0 - factor*1*VMag*VMag/2.0
              PrecondInv(2,1) = 0.0 - factor*u*VMag*VMag/2.0
              PrecondInv(3,1) = 0.0 - factor*v*VMag*VMag/2.0
              PrecondInv(4,1) = 0.0 - factor*w*VMag*VMag/2.0
              PrecondInv(5,1) = 0.0 - factor*H*VMag*VMag/2.0
              PrecondInv(6,1) = 0.0 - factor*tv*VMag*VMag/2.0
              PrecondInv(1,2) = 0.0 - factor*1*(-u)
              PrecondInv(2,2) = 1.0 - factor*u*(-u)
              PrecondInv(3,2) = 0.0 - factor*v*(-u)
              PrecondInv(4,2) = 0.0 - factor*w*(-u)
              PrecondInv(5,2) = 0.0 - factor*H*(-u)
              PrecondInv(6,2) = 0.0 - factor*tv*(-u)
              PrecondInv(1,3) = 0.0 - factor*1*(-v)
              PrecondInv(2,3) = 0.0 - factor*u*(-v)
              PrecondInv(3,3) = 1.0 - factor*v*(-v)
              PrecondInv(4,3) = 0.0 - factor*w*(-v)
              PrecondInv(5,3) = 0.0 - factor*H*(-v)
              PrecondInv(6,3) = 0.0 - factor*tv*(-v)
              PrecondInv(1,4) = 0.0 - factor*1*(-w)
              PrecondInv(2,4) = 0.0 - factor*u*(-w)
              PrecondInv(3,4) = 0.0 - factor*v*(-w)
              PrecondInv(4,4) = 1.0 - factor*w*(-w)
              PrecondInv(5,4) = 0.0 - factor*H*(-w)
              PrecondInv(6,4) = 0.0 - factor*tv*(-w)
              PrecondInv(1,5) = 0.0 - factor*1*(1.)
              PrecondInv(2,5) = 0.0 - factor*u*(1.)
              PrecondInv(3,5) = 0.0 - factor*v*(1.)
              PrecondInv(4,5) = 0.0 - factor*w*(1.)
              PrecondInv(5,5) = 1.0 - factor*H*(1.)
              PrecondInv(6,5) = 0.0 - factor*tv*(1.)
              PrecondInv(1,6) = 0.0 - factor*1*(-1.)
              PrecondInv(2,6) = 0.0 - factor*u*(-1.)
              PrecondInv(3,6) = 0.0 - factor*v*(-1.)
              PrecondInv(4,6) = 0.0 - factor*w*(-1.)
              PrecondInv(5,6) = 0.0 - factor*H*(-1.)
              PrecondInv(6,6) = 1.0 - factor*tv*(-1.)
              !---end preconditioning


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux

              D = (cells(i,j,k)%volume/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)

              ! -- source term derivatives -- !
              Omega = sqrt( ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
                           + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
                           + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
                            )&
                       )
              dist_i = dist(i,j,k)
              dist_i_2 = dist_i*dist_i
              k2 = kappa_sa*kappa_sa
              nu   = mu(i,j,k)/density
              Ji   = Q0(6)/nu
              Ji_2 = Ji*Ji
              Ji_3 = Ji_2*ji


              ! ___ functions ___
              fv1  = (Ji_3)/((Ji_3) + (cv1_3))
              fv2  = 1.0 - Ji/(1.0 + (Ji*fv1))
              ! ___ Shear stress for production ___
              S = Omega
              inv_k2_d2 = 1.0/(k2*dist_i_2)
              Shat      = S + Q0(6)*fv2*inv_k2_d2
              Shat      = max(Shat, 1.0e-10)
              inv_Shat  = 1.0/Shat
              dfv1 = 3.0*Ji_2*cv1_3/(nu*(Ji_3+cv1_3)**2)
              dfv2 = -((1.0/nu) - Ji_2*dfv1)/((1.0+Ji*fv1)**2)
              dShat = (fv2+Q0(6)*dfv2)*inv_k2_d2

              D = D - cb1*(Q0(6)*dShat+Shat)*cells(i,j,k)%volume

              ! ___ Destruction term___ !
              r    = min(Q0(6)*inv_Shat*inv_k2_d2, 10.0)
              g    = r + cw2*((r**6) - r)
              g_6  = g**6
              glim = ((1.0+cw3_6)/(g_6+cw3_6))**(1.0/6.0)
              fw   = g*glim
              dr = (Shat-Q0(6)*dShat)*inv_Shat*inv_Shat*inv_k2_d2
              dg = dr*(1.0+cw2*(6.0*(r**5)-1.0))
              dfw= dg*glim*(1.0-g_6/(g_6+cw3_6))

              D = D+cw1*(dfw*Q0(6) + 2*fw)*Q0(6)/dist_i_2*cells(i,j,k)%volume
              ! --  end of source term -- !

              delQ(i,j,k,1:6) = delQstar(i,j,k,1:6) &
                - 0.5*((matmul(PrecondInv,DelIminusFlux) - LambdaTimesArea(4)*delQ(i+1,j,k,1:6)) &
                     + (matmul(PrecondInv,DelJminusFlux) - LambdaTimesArea(5)*delQ(i,j+1,k,1:6)) &
                     + (matmul(PrecondInv,DelKminusFlux) - LambdaTimesArea(6)*delQ(i,j,k+1,1:6)) )/D

            end do
          end do
        end do
        
        do k=1,dims%kmx-1
          do j = 1,dims%jmx-1
            do i = 1,dims%imx-1
              conservativeQ(1) = qp(i,j,k,1)
              conservativeQ(2) = qp(i,j,k,1) * qp(i,j,k,2)
              conservativeQ(3) = qp(i,j,k,1) * qp(i,j,k,3)
              conservativeQ(4) = qp(i,j,k,1) * qp(i,j,k,4)
              conservativeQ(5) = (qp(i,j,k,5) / (gm-1.0)) + ( 0.5 * qp(i,j,k,1) * sum( qp(i,j,k,2:4)**2) )
              conservativeQ(6) = qp(i,j,k,1) * qp(i,j,k,6)
              
              ! add new change into conservative solution
              conservativeQ(1:6) = conservativeQ(1:6) + delQ(i,j,k,1:6)

              ! convert back conservative to primitive
              qp(i,j,k,1) = conservativeQ(1)
              qp(i,j,k,2) = conservativeQ(2) / conservativeQ(1)
              qp(i,j,k,3) = conservativeQ(3) / conservativeQ(1)
              qp(i,j,k,4) = conservativeQ(4) / conservativeQ(1)
              qp(i,j,k,5) = (gm-1.0) * ( conservativeQ(5) - (0.5 * sum(conservativeQ(2:4)**2) / conservativeQ(1)) )
              qp(i,j,k,6) = conservativeQ(6) / conservativeQ(1)
              qp(i,j,k,6) = max(qp(i,j,k,6), 1.e-8)
            end do
          end do
        end do

    end subroutine update_SA_variables


    function SAFlux(ql, qr, du, inputs)
      !< Calculate the total flux through face for turbulent flow (SA)
      !--------------------------------------
      ! calculate the total flux through face
      !---------------------------------------
      implicit none
      real, dimension(1:n_var), intent(in) :: ql !left state
      real, dimension(1:n_var), intent(in) :: qr !right state
      !conservative form of updated neighbour
      real, dimension(1:n_var), intent(in) :: du
      real, dimension(1:7)    , intent(in) :: inputs
      real, dimension(1:n_var)             :: Flux
      real, dimension(1:n_var)             :: SAFlux
      real, dimension(1:n_var)             :: U ! conservative variables
      real, dimension(1:n_var)             :: W ! new primitive variables
      real, dimension(1:n_var)             :: P ! primitive variables of right cell

      !for extraction of the inputs
      real :: area
      real :: nx
      real :: ny
      real :: nz
      real :: volume
      real :: mmu
      real :: tmu


      real    :: dudx
      real    :: dudy
      real    :: dudz
      real    :: dvdx
      real    :: dvdy
      real    :: dvdz
      real    :: dwdx
      real    :: dwdy
      real    :: dwdz
      real    :: dTdx
      real    :: dTdy
      real    :: dTdz
      real    :: dtvdx
      real    :: dtvdy
      real    :: dtvdz
      real    :: T1, T2
      real    :: uface
      real    :: vface
      real    :: wface
      real    :: trace
      real    :: Tauxx
      real    :: Tauyy
      real    :: Tauzz
      real    :: Tauxy
      real    :: Tauxz
      real    :: Tauyz
      real    :: Qx
      real    :: Qy
      real    :: Qz
      real    :: HalfRhoUsquare
      real    :: RhoHt
      real    :: K_heat
      real    :: FaceNormalVelocity
      real    :: mu
      real    :: muCap

      area   = inputs(1)
      nx     = inputs(2)
      ny     = inputs(3)
      nz     = inputs(4)
      volume = inputs(5)
      mmu    = inputs(6)
      tmu    = inputs(7)


      !save the old stat in P
      P = qr

      ! find conservative variable
      U(1)   =   ql(1)
      U(2)   =   ql(1) * ql(2)
      U(3)   =   ql(1) * ql(3)
      U(4)   =   ql(1) * ql(4)
      U(5)   = ( ql(5) / (gm-1.0) ) + ( 0.5 * ql(1) * sum(ql(2:4)**2) )
      U(6)   =   ql(1) * ql(6)

      U(1:n_var) = U(1:n_var) + du(1:n_var)
      

      W(1)   =   U(1)
      W(2)   =   U(2) / U(1)
      W(3)   =   U(3) / U(1)
      W(4)   =   U(4) / U(1)
      W(5)   = (gm-1.0) * ( U(5) - ( 0.5 * SUM(U(2:4)**2) / U(1) ) )
      W(6)   =   U(6) / U(1)
      W(6) = max(W(6), 1e-8)

      FaceNormalVelocity = (W(2) * nx) + (W(3) * ny) + (W(4) * nz)
      uface = 0.5 * ( W(2) + P(2) )
      vface = 0.5 * ( W(3) + P(3) )
      wface = 0.5 * ( W(4) + P(4) )


      Flux(1) =   W(1) * FaceNormalVelocity
      Flux(2) = ( W(2) * Flux(1) ) + ( W(5) * nx )
      Flux(3) = ( W(3) * Flux(1) ) + ( W(5) * ny )
      Flux(4) = ( W(4) * Flux(1) ) + ( W(5) * nz )

      HalfRhoUsquare = 0.5 * W(1) * ( W(2)*W(2) + W(3)*W(3) + W(4)*W(4) )
      RhoHt          = ( (gm/(gm-1.0)) * W(5) ) + HalfRhoUsquare
      Flux(5)        = RhoHt * FaceNormalVelocity
      Flux(6) = ( W(6) * Flux(1) )   


      ! viscous terms
      muCap = 0.25*(P(1)+W(1))*(P(6) + W(6))
      mu = mmu + tmu
      T1     =    W(5) / ( W(1) * R_gas )
      T2     =    P(5) / ( P(1) * R_gas )
      dTdx   =  ( T2   - T1   ) * nx * Area / Volume
      dTdy   =  ( T2   - T1   ) * ny * Area / Volume
      dTdz   =  ( T2   - T1   ) * nz * Area / Volume
      dudx   =  ( P(2) - W(2) ) * nx * Area / Volume
      dudy   =  ( P(2) - W(2) ) * ny * Area / Volume
      dudz   =  ( P(2) - W(2) ) * nz * Area / Volume
      dvdx   =  ( P(3) - W(3) ) * nx * Area / Volume
      dvdy   =  ( P(3) - W(3) ) * ny * Area / Volume
      dvdz   =  ( P(3) - W(3) ) * nz * Area / Volume
      dwdx   =  ( P(4) - W(4) ) * nx * Area / Volume
      dwdy   =  ( P(4) - W(4) ) * ny * Area / Volume
      dwdz   =  ( P(4) - W(4) ) * nz * Area / Volume
      dtvdx  =  ( P(6) - W(6) ) * nx * Area / Volume
      dtvdy  =  ( P(6) - W(6) ) * ny * Area / Volume
      dtvdz  =  ( P(6) - W(6) ) * nz * Area / Volume

      trace = dudx + dvdy + dwdz
      Tauxx =  2. * mu * (dudx - trace/3.0)
      Tauyy =  2. * mu * (dvdy - trace/3.0)
      Tauzz =  2. * mu * (dwdz - trace/3.0)
      Tauxy = mu * (dvdx + dudy)
      Tauxz = mu * (dwdx + dudz)
      Tauyz = mu * (dwdy + dvdz)

      K_heat = ( mmu / Pr  + tmu/tpr) * gm * R_gas / ( gm - 1.0 )
      Qx = K_heat*dTdx
      Qy = K_heat*dTdy
      Qz = K_heat*dTdz
      tmu = 0.5*(W(6) + P(6))

      Flux(2) = Flux(2) - ( Tauxx * nx + Tauxy * ny + Tauxz * nz )
      Flux(3) = Flux(3) - ( Tauxy * nx + Tauyy * ny + Tauyz * nz )
      Flux(4) = Flux(4) - ( Tauxz * nx + Tauyz * ny + Tauzz * nz )
      Flux(5) = Flux(5) - ( Tauxx * uface + Tauxy * vface + Tauxz * wface + Qx ) * nx
      Flux(5) = Flux(5) - ( Tauxy * uface + Tauyy * vface + Tauyz * wface + Qy ) * ny
      Flux(5) = Flux(5) - ( Tauxz * uface + Tauyz * vface + Tauzz * wface + Qz ) * nz
      Flux(6) = Flux(6) + (mmu + muCap)*(dtvdx*nx + dtvdy*ny + dtvdz*nz)/sigma_sa

      Flux    = Flux * Area
      SAFlux = Flux

    end function SAFlux 


    subroutine update_lctm2015(qp, residue, cells, Ifaces, Jfaces, Kfaces, dims)
      !< Update the RANS/transition (LCTM2015) equation with LU-SGS
      implicit none
      type(extent), intent(in) :: dims
      real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(inout) :: qp
      real, dimension(:, :, :, :), intent(in)  :: residue
      !< Store residue at each cell-center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      integer :: i,j,k
        real, dimension(1:8)     :: deltaU
        real, dimension(1:8)     :: D
        real, dimension(1:8)     :: conservativeQ
        real, dimension(1:8)     :: OldIminusFlux
        real, dimension(1:8)     :: OldJminusFlux
        real, dimension(1:8)     :: OldKminusFlux
        real, dimension(1:8)     :: NewIminusFlux
        real, dimension(1:8)     :: NewJminusFlux
        real, dimension(1:8)     :: NewKminusFlux
        real, dimension(1:8)     :: DelIminusFlux
        real, dimension(1:8)     :: DelJminusFlux
        real, dimension(1:8)     :: DelKminusFlux
        real, dimension(1:6)     :: LambdaTimesArea
        real, dimension(1:8)     :: Q0 ! state at cell
        real, dimension(1:8)     :: Q1 ! state at neighbours 
        real, dimension(1:8)     :: Q2
        real, dimension(1:8)     :: Q3
        real, dimension(1:8)     :: Q4
        real, dimension(1:8)     :: Q5
        real, dimension(1:8)     :: Q6
        real, dimension(1:8)     :: DQ0! change in state
        real, dimension(1:8)     :: DQ1
        real, dimension(1:8)     :: DQ2
        real, dimension(1:8)     :: DQ3
        real, dimension(1:8)     :: DQ4
        real, dimension(1:8)     :: DQ5
        real, dimension(1:8)     :: DQ6

        real, dimension(1:8)     :: Flist1
        real, dimension(1:8)     :: Flist2
        real, dimension(1:8)     :: Flist3
        real, dimension(1:8)     :: Flist4
        real, dimension(1:8)     :: Flist5
        real, dimension(1:8)     :: Flist6
        real, dimension(1:3)     :: C0
        real, dimension(1:3)     :: C1
        real, dimension(1:3)     :: C2
        real, dimension(1:3)     :: C3
        real, dimension(1:3)     :: C4
        real, dimension(1:3)     :: C5
        real, dimension(1:3)     :: C6
        real                     :: beta
        real                     :: eps
        real                     :: M
        real                     :: VMag
        real                     :: SoundMag
        real                     :: u,v,w,r,p,kk,ww,H,im
        real                     :: factor
        real, dimension(1:8,1:8) :: PrecondInv
        real, dimension(1:8,1:8) :: Identity

        ! intermittency
        real :: Fonset1
        real :: Fonset2
        real :: Fonset3
        real :: Fonset
        real :: Rev
        Real :: RT
        real :: Fturb
        real :: Re_theta
        real :: TuL
        real :: strain
        real :: vort
        real :: Dp, De
        real :: Fpg
        real :: dvdy
        real :: lamd
        real :: density
        Dp = 0.0
        De = 0.0

        !Identity matrix
        Identity = 0.0
        do i = 1,8
            Identity(i,i) = 1.0
        end do


        !intialize delQ
        delQstar = 0.0

        !forward sweep
        do k=1,dims%kmx-1
          do j=1,dims%jmx-1
            do i=1,dims%imx-1
              density = qp(i,j,k,1)
              C0  = (/Cells(i  ,j  ,k  )%Centerx,Cells(i  ,j  ,k  )%Centery,Cells(i  ,j  ,k  )%Centerz/)
              C1  = (/Cells(i-1,j  ,k  )%Centerx,Cells(i-1,j  ,k  )%Centery,Cells(i-1,j  ,k  )%Centerz/)
              C2  = (/Cells(i  ,j-1,k  )%Centerx,Cells(i  ,j-1,k  )%Centery,Cells(i  ,j-1,k  )%Centerz/)
              C3  = (/Cells(i  ,j  ,k-1)%Centerx,Cells(i  ,j  ,k-1)%Centery,Cells(i  ,j  ,k-1)%Centerz/)
              C4  = (/Cells(i+1,j  ,k  )%Centerx,Cells(i+1,j  ,k  )%Centery,Cells(i+1,j  ,k  )%Centerz/)
              C5  = (/Cells(i  ,j+1,k  )%Centerx,Cells(i  ,j+1,k  )%Centery,Cells(i  ,j+1,k  )%Centerz/)
              C6  = (/Cells(i  ,j  ,k+1)%Centerx,Cells(i  ,j  ,k+1)%Centery,Cells(i  ,j  ,k+1)%Centerz/)

              Q0  = qp(i  , j  , k  , 1:8)
              Q1  = qp(i-1, j  , k  , 1:8)
              Q2  = qp(i  , j-1, k  , 1:8)
              Q3  = qp(i  , j  , k-1, 1:8)
              Q4  = qp(i+1, j  , k  , 1:8)
              Q5  = qp(i  , j+1, k  , 1:8)
              Q6  = qp(i  , j  , k+1, 1:8)

              DQ0 = 0.0
              DQ1 = delQstar(i-1, j  , k  , 1:8)
              DQ2 = delQstar(i  , j-1, k  , 1:8)
              DQ3 = delQstar(i  , j  , k-1, 1:8)

              Flist1(1) =  Ifaces(i,j,k)%A
              Flist1(2) = -Ifaces(i,j,k)%nx
              Flist1(3) = -Ifaces(i,j,k)%ny
              Flist1(4) = -Ifaces(i,j,k)%nz
              Flist1(5) = 0.5*(cells(i-1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))
              Flist1(8) = 0.5*(sst_F1(i-1, j  , k  ) + sst_F1(i,j,k))

              Flist2(1) =  Jfaces(i,j,k)%A
              Flist2(2) = -Jfaces(i,j,k)%nx
              Flist2(3) = -Jfaces(i,j,k)%ny
              Flist2(4) = -Jfaces(i,j,k)%nz
              Flist2(5) = 0.5*(cells(i  , j-1, k  )%volume + cells(i,j,k)%volume)
              Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))
              Flist2(8) = 0.5*(sst_F1(i  , j-1, k  ) + sst_F1(i,j,k))

              Flist3(1) =  Kfaces(i,j,k)%A
              Flist3(2) = -Kfaces(i,j,k)%nx
              Flist3(3) = -Kfaces(i,j,k)%ny
              Flist3(4) = -Kfaces(i,j,k)%nz
              Flist3(5) = 0.5*(cells(i  , j  , k-1)%volume + cells(i,j,k)%volume)
              Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))
              Flist3(8) = 0.5*(sst_F1(i  , j  , k-1) + sst_F1(i,j,k))

              Flist4(1) =  Ifaces(i+1,j,k)%A
              Flist4(2) = +Ifaces(i+1,j,k)%nx
              Flist4(3) = +Ifaces(i+1,j,k)%ny
              Flist4(4) = +Ifaces(i+1,j,k)%nz
              Flist4(5) = 0.5*(cells(i+1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))
              Flist4(8) = 0.5*(sst_F1(i+1, j  , k  ) + sst_F1(i,j,k))

              Flist5(1) =  Jfaces(i,j+1,k)%A
              Flist5(2) = +Jfaces(i,j+1,k)%nx
              Flist5(3) = +Jfaces(i,j+1,k)%ny
              Flist5(4) = +Jfaces(i,j+1,k)%nz
              Flist5(5) = 0.5*(cells(i  , j+1, k  )%volume + cells(i,j,k)%volume)
              Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))
              Flist5(8) = 0.5*(sst_F1(i  , j+1, k  ) + sst_F1(i,j,k))

              Flist6(1) =  Kfaces(i,j,k+1)%A
              Flist6(2) = +Kfaces(i,j,k+1)%nx
              Flist6(3) = +Kfaces(i,j,k+1)%ny
              Flist6(4) = +Kfaces(i,j,k+1)%nz
              Flist6(5) = 0.5*(cells(i  , j  , k+1)%volume + cells(i,j,k)%volume)
              Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))
              Flist6(8) = 0.5*(sst_F1(i  , j  , k+1) + sst_F1(i,j,k))

              NewIminusFlux     = lctm2015Flux(Q1, Q0, DQ1, Flist1)
              NewJminusFlux     = lctm2015Flux(Q2, Q0, DQ2, Flist2)
              NewKminusFlux     = lctm2015Flux(Q3, Q0, DQ3, Flist3)
              OldIminusFlux     = lctm2015Flux(Q1, Q0, DQ0, Flist1)
              OldJminusFlux     = lctm2015Flux(Q2, Q0, DQ0, Flist2)
              OldKminusFlux     = lctm2015Flux(Q3, Q0, DQ0, Flist3)

              !---preconditioning---
              r  = Q0(1)
              u  = Q0(2)
              v  = Q0(3)
              w  = Q0(4)
              p  = Q0(5)
              kk = Q0(6)
              ww = Q0(7)
              im = Q0(8)
              VMag     = sqrt(u*u + v*v + w*w)
              SoundMag = sqrt(gm*p/r)
              M        = VMag/SoundMag 
              H  = (gm*p/(r*(gm-1.0))) + 0.5*(VMag)
              eps = min(1.0, max(M*M, Minf*Minf))
              factor = (1.0-eps)*(gm-1.0)/(SoundMag*SoundMag)

              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0, eps)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0, eps)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0, eps)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0, eps)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0, eps)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0, eps)


              !preconditing start
              PrecondInv(:,1) = VMag*VMag/2.0
              PrecondInv(:,2) = -u
              PrecondInv(:,3) = -v
              PrecondInv(:,4) = -w
              PrecondInv(:,5) =  1.0
              PrecondInv(:,6) = -1.0
              PrecondInv(:,7) =  0.0
              PrecondInv(:,8) =  0.0
              PrecondInv(2,:) = u*PrecondInv(2,:)
              PrecondInv(3,:) = v*PrecondInv(3,:)
              PrecondInv(4,:) = w*PrecondInv(4,:)
              PrecondInv(5,:) = H*PrecondInv(5,:)
              PrecondInv(6,:) =kk*PrecondInv(6,:)
              PrecondInv(7,:) =ww*PrecondInv(7,:)
              PrecondInv(8,:) =im*PrecondInv(8,:)
              PrecondInv = Identity - (factor*PrecondInv)
              !---end preconditioning


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux


              D = (cells(i,j,k)%volume/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
              beta = sst_F1(i,j,k)*beta1 + (1.0-sst_F1(i,j,k))*beta2
              !D(6) = (D(6) + bstar*qp(i,j,k,7)*cells(i,j,k)%volume)
              D(6) = (D(6) + (bstar*qp(i,j,k,7))*cells(i,j,k)%volume)
              D(7) = (D(7) + 2.0*beta*qp(i,j,k,7)*cells(i,j,k)%volume)
              !gamma
              vort = sqrt(     ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
                              + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
                              + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
                               )&
                         )

              strain = sqrt(     ((gradw_y(i,j,k) + gradv_z(i,j,k))**2 &
                                + (gradu_z(i,j,k) + gradw_x(i,j,k))**2 &
                                + (gradv_x(i,j,k) + gradu_y(i,j,k))**2 &
                                + 2*(gradu_x(i,j,k))**2 &
                                + 2*(gradv_y(i,j,k))**2 &
                                + 2*(gradw_z(i,j,k))**2 &
                                 )&
                           )
              dvdy = DCCVnX(i,j,k)*CCnormalX(i,j,k) &
                   + DCCVnY(i,j,k)*CCnormalY(i,j,k) &
                   + DCCVnZ(i,j,k)*CCnormalZ(i,j,k)
              lamd =(-7.57e-3)*(dvdy*dist(i,j,k)*dist(i,j,k)*density/mu(i,j,k)) + 0.0128
              lamd = min(max(lamd, -1.0), 1.0)
              if(lamd>=0.0)then
                  Fpg = min(1.0 + 14.68*lamd, 1.5)
              else
                  Fpg = min(1.0 - 7.34*lamd, 3.0)
              end if
              Fpg = max(Fpg, 0.0)
              TuL = min(100.0*sqrt(2.0*qp(i,j,k,6)/3.0)/(qp(i,j,k,7)*dist(i,j,k)),100.0)
              Re_theta = 100.0 + 1000.0*exp(-TuL*Fpg)
              Rev = density*dist(i,j,k)*dist(i,j,k)*strain/mu(i,j,k)
              RT = density*qp(i,j,k,6)/(mu(i,j,k)*qp(i,j,k,7))
              Fturb = exp(-(0.5*Rt)**4)
              Fonset1 = Rev/(2.2*Re_theta)
              Fonset2 = min(Fonset1, 2.0)
              Fonset3 = max(1.0 - (RT/3.5)**3, 0.0)
              Fonset  = max(Fonset2 - Fonset3, 0.0)
              Dp = 100*density*strain*Fonset*(1.0-2.0*Q0(8))
              De = 0.06*vort*Fturb*density*(2.0*50.0*Q0(8) - 1.0)
              D(8) = (D(8) + (-Dp + DE )*cells(i,j,k)%volume)
              !storing D in Iflux array for backward sweep
              !F_p(i,j,k,1) = D

              deltaU(1:8) = -matmul(PrecondInv,residue(i,j,k,1:8)) &
                - 0.5*((matmul(PrecondInv,DelIminusFlux) - LambdaTimesArea(1)*delQstar(i-1,j,k,1:8)) &
                     + (matmul(PrecondInv,DelJminusFlux) - LambdaTimesArea(2)*delQstar(i,j-1,k,1:8)) &
                     + (matmul(PrecondInv,DelKminusFlux) - LambdaTimesArea(3)*delQstar(i,j,k-1,1:8)) )

              delQstar(i,j,k,1:8) = deltaU(1:8)/D
            end do
          end do
        end do

        delQ=0.0
        !backward sweep
            do i=dims%imx-1,1,-1
          do j=dims%jmx-1,1,-1
        do k=dims%kmx-1,1,-1
              density = qp(i,j,k,1)
              C0  = (/Cells(i  ,j  ,k  )%Centerx,Cells(i  ,j  ,k  )%Centery,Cells(i  ,j  ,k  )%Centerz/)
              C1  = (/Cells(i-1,j  ,k  )%Centerx,Cells(i-1,j  ,k  )%Centery,Cells(i-1,j  ,k  )%Centerz/)
              C2  = (/Cells(i  ,j-1,k  )%Centerx,Cells(i  ,j-1,k  )%Centery,Cells(i  ,j-1,k  )%Centerz/)
              C3  = (/Cells(i  ,j  ,k-1)%Centerx,Cells(i  ,j  ,k-1)%Centery,Cells(i  ,j  ,k-1)%Centerz/)
              C4  = (/Cells(i+1,j  ,k  )%Centerx,Cells(i+1,j  ,k  )%Centery,Cells(i+1,j  ,k  )%Centerz/)
              C5  = (/Cells(i  ,j+1,k  )%Centerx,Cells(i  ,j+1,k  )%Centery,Cells(i  ,j+1,k  )%Centerz/)
              C6  = (/Cells(i  ,j  ,k+1)%Centerx,Cells(i  ,j  ,k+1)%Centery,Cells(i  ,j  ,k+1)%Centerz/)

              Q0  = qp(i  , j  , k  , 1:8)
              Q1  = qp(i-1, j  , k  , 1:8)
              Q2  = qp(i  , j-1, k  , 1:8)
              Q3  = qp(i  , j  , k-1, 1:8)
              Q4  = qp(i+1, j  , k  , 1:8)
              Q5  = qp(i  , j+1, k  , 1:8)
              Q6  = qp(i  , j  , k+1, 1:8)

              DQ0 = 0.0
              DQ4 = delQ(i+1, j  , k  , 1:8)
              DQ5 = delQ(i  , j+1, k  , 1:8)
              DQ6 = delQ(i  , j  , k+1, 1:8)

              Flist1(1) =  Ifaces(i,j,k)%A
              Flist1(2) = -Ifaces(i,j,k)%nx
              Flist1(3) = -Ifaces(i,j,k)%ny
              Flist1(4) = -Ifaces(i,j,k)%nz
              Flist1(5) = 0.5*(cells(i-1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))
              Flist1(8) = 0.5*(sst_F1(i-1, j  , k  ) + sst_F1(i,j,k))

              Flist2(1) =  Jfaces(i,j,k)%A
              Flist2(2) = -Jfaces(i,j,k)%nx
              Flist2(3) = -Jfaces(i,j,k)%ny
              Flist2(4) = -Jfaces(i,j,k)%nz
              Flist2(5) = 0.5*(cells(i  , j-1, k  )%volume + cells(i,j,k)%volume)
              Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))
              Flist2(8) = 0.5*(sst_F1(i  , j-1, k  ) + sst_F1(i,j,k))

              Flist3(1) =  Kfaces(i,j,k)%A
              Flist3(2) = -Kfaces(i,j,k)%nx
              Flist3(3) = -Kfaces(i,j,k)%ny
              Flist3(4) = -Kfaces(i,j,k)%nz
              Flist3(5) = 0.5*(cells(i  , j  , k-1)%volume + cells(i,j,k)%volume)
              Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))
              Flist3(8) = 0.5*(sst_F1(i  , j  , k-1) + sst_F1(i,j,k))

              Flist4(1) =  Ifaces(i+1,j,k)%A
              Flist4(2) = +Ifaces(i+1,j,k)%nx
              Flist4(3) = +Ifaces(i+1,j,k)%ny
              Flist4(4) = +Ifaces(i+1,j,k)%nz
              Flist4(5) = 0.5*(cells(i+1, j  , k  )%volume + cells(i,j,k)%volume)
              Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))
              Flist4(8) = 0.5*(sst_F1(i+1, j  , k  ) + sst_F1(i,j,k))

              Flist5(1) =  Jfaces(i,j+1,k)%A
              Flist5(2) = +Jfaces(i,j+1,k)%nx
              Flist5(3) = +Jfaces(i,j+1,k)%ny
              Flist5(4) = +Jfaces(i,j+1,k)%nz
              Flist5(5) = 0.5*(cells(i  , j+1, k  )%volume + cells(i,j,k)%volume)
              Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))
              Flist5(8) = 0.5*(sst_F1(i  , j+1, k  ) + sst_F1(i,j,k))

              Flist6(1) =  Kfaces(i,j,k+1)%A
              Flist6(2) = +Kfaces(i,j,k+1)%nx
              Flist6(3) = +Kfaces(i,j,k+1)%ny
              Flist6(4) = +Kfaces(i,j,k+1)%nz
              Flist6(5) = 0.5*(cells(i  , j  , k+1)%volume + cells(i,j,k)%volume)
              Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))
              Flist6(8) = 0.5*(sst_F1(i  , j  , k+1) + sst_F1(i,j,k))

              NewIminusFlux     = lctm2015Flux(Q4, Q0, DQ4, Flist4)
              NewJminusFlux     = lctm2015Flux(Q5, Q0, DQ5, Flist5)
              NewKminusFlux     = lctm2015Flux(Q6, Q0, DQ6, Flist6)
              OldIminusFlux     = lctm2015Flux(Q4, Q0, DQ0, Flist4)
              OldJminusFlux     = lctm2015Flux(Q5, Q0, DQ0, Flist5)
              OldKminusFlux     = lctm2015Flux(Q6, Q0, DQ0, Flist6)

              !---preconditioning---
              r  = Q0(1)
              u  = Q0(2)
              v  = Q0(3)
              w  = Q0(4)
              p  = Q0(5)
              VMag     = sqrt(u*u + v*v + w*w)
              SoundMag = sqrt(gm*p/r)
              M        = VMag/SoundMag 
              H  = (gm*p/(r*(gm-1.0))) + 0.5*(VMag)
              eps = min(1.0, max(M*M, Minf*Minf))
              factor = (1.0-eps)*(gm-1.0)/(SoundMag*SoundMag)

              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0, eps)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0, eps)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0, eps)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0, eps)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0, eps)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0, eps)


              !preconditing start
              PrecondInv(:,1) = VMag*VMag/2.0
              PrecondInv(:,2) = -u
              PrecondInv(:,3) = -v
              PrecondInv(:,4) = -w
              PrecondInv(:,5) =  1.0
              PrecondInv(:,6) = -1.0
              PrecondInv(:,7) =  0.0
              PrecondInv(:,8) =  0.0
              PrecondInv(2,:) = u*PrecondInv(2,:)
              PrecondInv(3,:) = v*PrecondInv(3,:)
              PrecondInv(4,:) = w*PrecondInv(4,:)
              PrecondInv(5,:) = H*PrecondInv(5,:)
              PrecondInv(6,:) =kk*PrecondInv(6,:)
              PrecondInv(7,:) =ww*PrecondInv(7,:)
              PrecondInv(8,:) =im*PrecondInv(8,:)
              PrecondInv = Identity - (factor*PrecondInv)
              !---end preconditioning


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux

              D = (cells(i,j,k)%volume/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
              beta = sst_F1(i,j,k)*beta1 + (1.0-sst_F1(i,j,k))*beta2
              !D(6) = (D(6) + bstar*qp(i,j,k,7)*cells(i,j,k)%volume)
              D(6) = (D(6) + (bstar*qp(i,j,k,7))*cells(i,j,k)%volume)
              D(7) = (D(7) + 2.0*beta*qp(i,j,k,7)*cells(i,j,k)%volume)
              !gamma
              vort = sqrt(     ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
                              + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
                              + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
                               )&
                         )

              strain = sqrt(     ((gradw_y(i,j,k) + gradv_z(i,j,k))**2 &
                                + (gradu_z(i,j,k) + gradw_x(i,j,k))**2 &
                                + (gradv_x(i,j,k) + gradu_y(i,j,k))**2 &
                                + 2*(gradu_x(i,j,k))**2 &
                                + 2*(gradv_y(i,j,k))**2 &
                                + 2*(gradw_z(i,j,k))**2 &
                                 )&
                           )
              dvdy = DCCVnX(i,j,k)*CCnormalX(i,j,k) &
                   + DCCVnY(i,j,k)*CCnormalY(i,j,k) &
                   + DCCVnZ(i,j,k)*CCnormalZ(i,j,k)
              lamd =(-7.57e-3)*(dvdy*dist(i,j,k)*dist(i,j,k)*density/mu(i,j,k)) + 0.0128
              lamd = min(max(lamd, -1.0), 1.0)
              if(lamd>=0.0)then
                  Fpg = min(1.0 + 14.68*lamd, 1.5)
              else
                  Fpg = min(1.0 - 7.34*lamd, 3.0)
              end if
              Fpg = max(Fpg, 0.0)
              TuL = min(100.0*sqrt(2.0*qp(i,j,k,6)/3.0)/(qp(i,j,k,7)*dist(i,j,k)),100.0)
              Re_theta = 100.0 + 1000.0*exp(-TuL*Fpg)
              Rev = density*dist(i,j,k)*dist(i,j,k)*strain/mu(i,j,k)
              RT = density*qp(i,j,k,6)/(mu(i,j,k)*qp(i,j,k,7))
              Fturb = exp(-(0.5*Rt)**4)
              Fonset1 = Rev/(2.2*Re_theta)
              Fonset2 = min(Fonset1, 2.0)
              Fonset3 = max(1.0 - (RT/3.5)**3, 0.0)
              Fonset  = max(Fonset2 - Fonset3, 0.0)
              Dp = 100*density*strain*Fonset*(1.0-2.0*Q0(8))
              De = 0.06*vort*Fturb*density*(2.0*50.0*Q0(8) - 1.0)
              D(8) = (D(8) + (-Dp + DE )*cells(i,j,k)%volume)

              delQ(i,j,k,1:8) = delQstar(i,j,k,1:8) &
                - 0.5*((matmul(PrecondInv,DelIminusFlux) - LambdaTimesArea(4)*delQ(i+1,j,k,1:8)) &
                     + (matmul(PrecondInv,DelJminusFlux) - LambdaTimesArea(5)*delQ(i,j+1,k,1:8)) &
                     + (matmul(PrecondInv,DelKminusFlux) - LambdaTimesArea(6)*delQ(i,j,k+1,1:8)) )/D

            end do
          end do
        end do
        
        do k=1,dims%kmx-1
          do j = 1,dims%jmx-1
            do i = 1,dims%imx-1
              conservativeQ(1) = qp(i,j,k,1)
              conservativeQ(2) = qp(i,j,k,1) * qp(i,j,k,2)
              conservativeQ(3) = qp(i,j,k,1) * qp(i,j,k,3)
              conservativeQ(4) = qp(i,j,k,1) * qp(i,j,k,4)
              conservativeQ(5) = (qp(i,j,k,5) / (gm-1.0)) + ( 0.5 * qp(i,j,k,1) * sum( qp(i,j,k,2:4)**2) )
              conservativeQ(6) = qp(i,j,k,1) * qp(i,j,k,6)
              conservativeQ(7) = qp(i,j,k,1) * qp(i,j,k,7)
              conservativeQ(8) = qp(i,j,k,1) * qp(i,j,k,8)
              
              ! add new change into conservative solution
              conservativeQ(1:n_var) = conservativeQ(1:n_var) + delQ(i,j,k,1:n_var)

              ! convert back conservative to primitive
              qp(i,j,k,1) = conservativeQ(1)
              qp(i,j,k,2) = conservativeQ(2) / conservativeQ(1)
              qp(i,j,k,3) = conservativeQ(3) / conservativeQ(1)
              qp(i,j,k,4) = conservativeQ(4) / conservativeQ(1)
              qp(i,j,k,5) = (gm-1.0) * ( conservativeQ(5) - (0.5 * sum(conservativeQ(2:4)**2) / conservativeQ(1)) )
              if(conservativeQ(6)>0.0)then
              qp(i,j,k,6) = conservativeQ(6) / conservativeQ(1)
              end if
              if(conservativeQ(7)>0.0)then
              qp(i,j,k,7) = conservativeQ(7) / conservativeQ(1)
              end if
              qp(i,j,k,8) = conservativeQ(8) / conservativeQ(1)
              qp(i,j,k,8) = max(qp(i,j,k,8), 0.0)
              !qp(i,j,k,8) = min(qp(i,j,k,8), 1.0)
            end do
          end do
        end do
    end subroutine update_lctm2015

    function lctm2015flux(ql, qr, du, inputs)
      !< Calculate the total flux through face for turbulent/transition flow (LCTM2015)
      !--------------------------------------
      ! calculate the total flux through face
      !---------------------------------------
      implicit none
      real, dimension(1:n_var), intent(in) :: ql !left state
      real, dimension(1:n_var), intent(in) :: qr !right state
      !conservative form of updated neighbour
      real, dimension(1:n_var), intent(in) :: du
      real, dimension(1:8)    , intent(in) :: inputs
      real, dimension(1:n_var)             :: Flux
      real, dimension(1:n_var)             :: lctm2015flux
      real, dimension(1:n_var)             :: U ! conservative variables
      real, dimension(1:n_var)             :: W ! new primitive variables
      real, dimension(1:n_var)             :: P ! primitive variables of right cell

      !for extraction of the inputs
      real :: area
      real :: nx
      real :: ny
      real :: nz
      real :: volume
      real :: mmu
      real :: tmu


      real    :: dudx
      real    :: dudy
      real    :: dudz
      real    :: dvdx
      real    :: dvdy
      real    :: dvdz
      real    :: dwdx
      real    :: dwdy
      real    :: dwdz
      real    :: dTdx
      real    :: dTdy
      real    :: dTdz
      real    :: dtkdx
      real    :: dtkdy
      real    :: dtkdz
      real    :: dtwdx
      real    :: dtwdy
      real    :: dtwdz
      real    :: dtgmdx
      real    :: dtgmdy
      real    :: dtgmdz
      real    :: T1, T2
      real    :: uface
      real    :: vface
      real    :: wface
      real    :: trace
      real    :: Tauxx
      real    :: Tauyy
      real    :: Tauzz
      real    :: Tauxy
      real    :: Tauxz
      real    :: Tauyz
      real    :: Qx
      real    :: Qy
      real    :: Qz
      real    :: HalfRhoUsquare
      real    :: RhoHt
      real    :: K_heat
      real    :: FaceNormalVelocity
      real    :: mu
      real    :: sigma_k
      real    :: sigma_w
      real    :: F1

      area   = inputs(1)
      nx     = inputs(2)
      ny     = inputs(3)
      nz     = inputs(4)
      volume = inputs(5)
      mmu    = inputs(6)
      tmu    = inputs(7)
      F1     = inputs(8)


      !save the old stat in P
      P = qr

      ! find conservative variable
      U(1)   =   ql(1)
      U(2)   =   ql(1) * ql(2)
      U(3)   =   ql(1) * ql(3)
      U(4)   =   ql(1) * ql(4)
      U(5)   = ( ql(5) / (gm-1.0) ) + ( 0.5 * ql(1) * sum(ql(2:4)**2) )
      U(6)   =   ql(1) * ql(6)
      U(7)   =   ql(1) * ql(7)
      U(8)   =   ql(1) * ql(8)

      U(1:n_var) = U(1:n_var) + du(1:n_var)
      

      W(1)   =   U(1)
      W(2)   =   U(2) / U(1)
      W(3)   =   U(3) / U(1)
      W(4)   =   U(4) / U(1)
      W(5)   = (gm-1.0) * ( U(5) - ( 0.5 * SUM(U(2:4)**2) / U(1) ) )
      W(6)   =   U(6) / U(1)
      W(7)   =   U(7) / U(1)
      W(8)   =   U(8) / U(1)
      W(6)   = W(6) + 0.5*(1.-sign(1.,W(6)))*(ql(6)-W(6))
      W(7)   = W(7) + 0.5*(1.-sign(1.,W(7)))*(ql(7)-W(7))
      W(8) = max(W(8), 0.0)
      !W(8) = min(W(8), 1.0)

      FaceNormalVelocity = (W(2) * nx) + (W(3) * ny) + (W(4) * nz)
      uface = 0.5 * ( W(2) + P(2) )
      vface = 0.5 * ( W(3) + P(3) )
      wface = 0.5 * ( W(4) + P(4) )


      Flux(1) =   W(1) * FaceNormalVelocity
      Flux(2) = ( W(2) * Flux(1) ) + ( W(5) * nx )
      Flux(3) = ( W(3) * Flux(1) ) + ( W(5) * ny )
      Flux(4) = ( W(4) * Flux(1) ) + ( W(5) * nz )

      HalfRhoUsquare = 0.5 * W(1) * ( W(2)*W(2) + W(3)*W(3) + W(4)*W(4) )
      RhoHt          = ( (gm/(gm-1.0)) * W(5) ) + HalfRhoUsquare
      Flux(5)        = RhoHt * FaceNormalVelocity
      Flux(6) = ( W(6) * Flux(1) )   
      Flux(7) = ( W(7) * Flux(1) )   
      Flux(8) = ( W(8) * Flux(1) )   


      ! viscous terms
      mu = mmu + tmu
      T1     =    W(5) / ( W(1) * R_gas )
      T2     =    P(5) / ( P(1) * R_gas )
      dTdx   =  ( T2   - T1   ) * nx * Area / Volume
      dTdy   =  ( T2   - T1   ) * ny * Area / Volume
      dTdz   =  ( T2   - T1   ) * nz * Area / Volume
      dudx   =  ( P(2) - W(2) ) * nx * Area / Volume
      dudy   =  ( P(2) - W(2) ) * ny * Area / Volume
      dudz   =  ( P(2) - W(2) ) * nz * Area / Volume
      dvdx   =  ( P(3) - W(3) ) * nx * Area / Volume
      dvdy   =  ( P(3) - W(3) ) * ny * Area / Volume
      dvdz   =  ( P(3) - W(3) ) * nz * Area / Volume
      dwdx   =  ( P(4) - W(4) ) * nx * Area / Volume
      dwdy   =  ( P(4) - W(4) ) * ny * Area / Volume
      dwdz   =  ( P(4) - W(4) ) * nz * Area / Volume
      dtkdx  =  ( P(6) - W(6) ) * nx * Area / Volume
      dtkdy  =  ( P(6) - W(6) ) * ny * Area / Volume
      dtkdz  =  ( P(6) - W(6) ) * nz * Area / Volume
      dtwdx  =  ( P(7) - W(7) ) * nx * Area / Volume
      dtwdy  =  ( P(7) - W(7) ) * ny * Area / Volume
      dtwdz  =  ( P(7) - W(7) ) * nz * Area / Volume
      dtgmdx  =  ( P(8) - W(8) ) * nx * Area / Volume
      dtgmdy  =  ( P(8) - W(8) ) * ny * Area / Volume
      dtgmdz  =  ( P(8) - W(8) ) * nz * Area / Volume

      trace = dudx + dvdy + dwdz
      Tauxx =  2. * mu * (dudx - trace/3.0)
      Tauyy =  2. * mu * (dvdy - trace/3.0)
      Tauzz =  2. * mu * (dwdz - trace/3.0)
      Tauxy = mu * (dvdx + dudy)
      Tauxz = mu * (dwdx + dudz)
      Tauyz = mu * (dwdy + dvdz)

      K_heat = ( mmu / Pr  + tmu/tpr) * gm * R_gas / ( gm - 1.0 )
      Qx = K_heat*dTdx
      Qy = K_heat*dTdy
      Qz = K_heat*dTdz

      sigma_k = sigma_k1*F1 + sigma_k2*(1.0 - F1)
      sigma_w = sigma_w1*F1 + sigma_w2*(1.0 - F1)

      Flux(2) = Flux(2) - ( Tauxx * nx + Tauxy * ny + Tauxz * nz )
      Flux(3) = Flux(3) - ( Tauxy * nx + Tauyy * ny + Tauyz * nz )
      Flux(4) = Flux(4) - ( Tauxz * nx + Tauyz * ny + Tauzz * nz )
      Flux(5) = Flux(5) - ( Tauxx * uface + Tauxy * vface + Tauxz * wface + Qx ) * nx
      Flux(5) = Flux(5) - ( Tauxy * uface + Tauyy * vface + Tauyz * wface + Qy ) * ny
      Flux(5) = Flux(5) - ( Tauxz * uface + Tauyz * vface + Tauzz * wface + Qz ) * nz
      Flux(6) = Flux(6) + (mmu + sigma_k*tmu)*(dtkdx*nx + dtkdy*ny + dtkdz*nz)
      Flux(7) = Flux(7) + (mmu + sigma_w*tmu)*(dtwdx*nx + dtwdy*ny + dtwdz*nz)
      Flux(8) = Flux(8) + (mmu + tmu)*(dtgmdx*nx + dtgmdy*ny + dtgmdz*nz)

      Flux    = Flux * Area
      lctm2015flux = Flux
    end function lctm2015flux


end module plusgs
