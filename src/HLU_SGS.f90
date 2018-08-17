module HLU_SGS
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
  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx
  use global_vars, only : R_gas
  use global_vars, only : Pr
  use global_vars, only : tPr
  use global_vars, only : time_step_accuracy

  use global_vars, only : volume
  use global_vars, only : xnx, xny, xnz !face unit normal x
  use global_vars, only : ynx, yny, ynz !face unit normal y
  use global_vars, only : znx, zny, znz !face unit normal z
  use global_vars, only : xA, yA, zA    !face area
    
  use global_vars, only : n_var
  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx
  use global_vars, only : gm
  use global_vars, only : sst_n_var
  use global_vars, only : qp
  use global_vars, only : qp_inf
  use global_vars, only : density
  use global_vars, only : x_speed
  use global_vars, only : y_speed
  use global_vars, only : z_speed
  use global_vars, only : pressure
  use global_vars, only : dist
  use global_vars, only : mu
  use global_vars, only : mu_t
  use global_vars, only : tk_inf
  use global_vars, only : tkl_inf

  use global_vars, only : delta_t
  use global_vars, only : turbulence
  use global_vars, only : process_id

  use global_vars, only: F_p
  use global_vars, only: G_p
  use global_vars, only: H_p
  use global_vars, only: mass_residue
  use global_vars, only: x_mom_residue
  use global_vars, only: y_mom_residue
  use global_vars, only: z_mom_residue
  use global_vars, only: energy_residue
  use global_vars, only: TKE_residue
  use global_vars, only: omega_residue
  use global_vars, only: kl_residue
  use global_vars, only: residue
  use global_vars, only: mu_ref

  use global_vars, only: gradu_x
  use global_vars, only: gradu_y
  use global_vars, only: gradu_z
  use global_vars, only: gradv_x
  use global_vars, only: gradv_y
  use global_vars, only: gradv_z
  use global_vars, only: gradw_x
  use global_vars, only: gradw_y
  use global_vars, only: gradw_z
  use geometry   , only: CellCenter

  use utils, only: alloc
  use utils, only:  dealloc 
  use utils, only:  dmsg
  use utils, only:  DEBUG_LEVEL

  use string

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


#include "error.inc"
#include "debug.h"
#include "mpi.inc"


  real, dimension(:,:,:,:), allocatable :: Diagonal ! diagonal terms
  real, dimension(:,:,:,:), allocatable :: delQ ! initial delta
  real, dimension(:,:,:,:), allocatable :: delQk12 ! initial delta
  real, dimension(:,:,:,:), allocatable :: delQstar ! initial delta
  real, dimension(:,:,:,:), allocatable :: delQP ! Plus, New iteration 
  real, dimension(:,:,:,:), allocatable :: delQM ! Minus, Last iteration 
  real, dimension(:,:,:), allocatable, target :: dummy
  real, dimension(:,:,:), pointer :: tmu
  real, dimension(:,:,:), pointer :: mmu
  real, dimension(:,:,:), pointer :: blendF
  integer :: SGSIter = 5
        
  !parallel communication
  integer :: ibuf_size
  integer :: jbuf_size
  integer :: kbuf_size
  real, dimension(:), allocatable :: imin_send_buf
  real, dimension(:), allocatable :: jmin_send_buf
  real, dimension(:), allocatable :: kmin_send_buf
  real, dimension(:), allocatable :: imin_recv_buf
  real, dimension(:), allocatable :: jmin_recv_buf
  real, dimension(:), allocatable :: kmin_recv_buf
  real, dimension(:), allocatable :: imax_send_buf
  real, dimension(:), allocatable :: jmax_send_buf
  real, dimension(:), allocatable :: kmax_send_buf
  real, dimension(:), allocatable :: imax_recv_buf
  real, dimension(:), allocatable :: jmax_recv_buf
  real, dimension(:), allocatable :: kmax_recv_buf

  public :: update_with_HLU_SGS
  public :: setup_HLU_SGS
  public :: destroy_HLU_SGS

  contains

    subroutine setup_HLU_SGS()
      implicit none
      character(len=*), parameter :: &
        errmsg="module: DP_LUSGS, subrouinte setup"

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
      call alloc(delQk12, 0, imx, 0, jmx, 0, kmx, 1, n_var)
      call alloc(Diagonal, 0, imx, 0, jmx, 0, kmx, 1, n_var)

      if(mu_ref==0.0 .or. turbulence=='none') then
        call alloc(dummy, -2, imx+2, -2, jmx+2, -2, kmx+2)
        dummy = 0.0
      end if
      if(mu_ref==0.0)then
        mmu(-2:imx+2,-2:jmx+2,-2:kmx+2) => dummy(:,:,:)
      else
        mmu(-2:imx+2,-2:jmx+2,-2:kmx+2) => mu(:,:,:)
      end if
      if(trim(turbulence)=='none')then
        tmu(-2:imx+2,-2:jmx+2,-2:kmx+2) => dummy(:,:,:)
      else
        tmu(-2:imx+2,-2:jmx+2,-2:kmx+2) => mu_t(:,:,:)
      end if
      if((turbulence/='sst' .or. turbulence/='sst2003') .and. turbulence/='none') then
        call alloc(dummy, -2, imx+2, -2, jmx+2, -2, kmx+2)
        dummy = 0.0
        blendF(-2:imx+2, -2:jmx+2, -2:kmx+2) => dummy(:,:,:)
      elseif(turbulence=='none')then
        blendF(-2:imx+2, -2:jmx+2, -2:kmx+2) => dummy(:,:,:)
      else
        blendF(-2:imx+2,-2:jmx+2,-2:kmx+2) => sst_F1(:,:,:)
      end if
    end subroutine setup_HLU_SGS

    subroutine destroy_HLU_SGS()
      implicit none
      call dealloc(imin_send_buf)
      call dealloc(jmin_send_buf)
      call dealloc(kmin_send_buf)
      call dealloc(imin_recv_buf)
      call dealloc(jmin_recv_buf)
      call dealloc(kmin_recv_buf)
      call dealloc(imax_send_buf)
      call dealloc(jmax_send_buf)
      call dealloc(kmax_send_buf)
      call dealloc(imax_recv_buf)
      call dealloc(jmax_recv_buf)
      call dealloc(kmax_recv_buf)
      call dealloc(delQ)
      call dealloc(delQstar)
      call dealloc(delQk12)
      call dealloc(Diagonal)
      call dealloc(dummy)
    end subroutine destroy_HLU_SGS

    subroutine update_with_HLU_SGS()
      implicit none

          call update_laminar_variables()
          !call update_SST_variables()
!      select case(trim(turbulence))
!        case('none')
!          call update_laminar_variables()
!
!        case('sst', 'sst2003')
!          call update_SST_variables()
!
!        case('kkl')
!          call update_KKL_variables()
!
!        case('sa', 'saBC')
!          call update_SA_variables()
!
!        case Default
!          Fatal_error
!
!      end select

    end subroutine update_with_HLU_SGS


    subroutine UpdateState(qp,delQ,n_var)
      implicit none
      integer :: i,j,k
      integer :: n_var
      real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(inout)  :: qp
      real, dimension(0:imx,0:jmx,0:kmx,1:n_var), intent(in) :: delQ
      real, dimension(1:n_var) :: ConservativeQ

        do k=1,kmx-1
          do j = 1,jmx-1
            do i = 1,imx-1

              ! convert primitive to conservative
              conservativeQ(1) = qp(i,j,k,1)
              conservativeQ(2:n_var) = qp(i,j,k,1) * qp(i,j,k,2:n_var)
              conservativeQ(5) = (qp(i,j,k,5) / (gm-1.0)) + ( 0.5 * qp(i,j,k,1) * sum( qp(i,j,k,2:4)**2) )
              
              ! add new change into conservative solution
              conservativeQ(1:n_var) = conservativeQ(1:n_var) + delQ(i,j,k,1:n_var)

              ! convert back conservative to primitive
              qp(i,j,k,1) = conservativeQ(1)
              qp(i,j,k,2:n_var) = conservativeQ(2:n_var) / conservativeQ(1)
              qp(i,j,k,5) = (gm-1.0) * ( conservativeQ(5) - (0.5 * sum(conservativeQ(2:4)**2) / conservativeQ(1)) )
              ! Check for negative values
              qp(i,j,k,6:n_var) = max(qp(i,j,k,6:n_var), 1.e-10)
            end do
          end do
        end do

    end subroutine UpdateState


    subroutine update_laminar_variables()
      implicit none
      integer :: i,j,k, SubIter
        real, dimension(1:n_var)     :: deltaU
        real, dimension(1:n_var)     :: D
        real, dimension(1:n_var)     :: conservativeQ
        real, dimension(1:n_var)     :: OldIminusFlux
        real, dimension(1:n_var)     :: OldJminusFlux
        real, dimension(1:n_var)     :: OldKminusFlux
        real, dimension(1:n_var)     :: NewIminusFlux
        real, dimension(1:n_var)     :: NewJminusFlux
        real, dimension(1:n_var)     :: NewKminusFlux
        real, dimension(1:n_var)     :: DelIminusFlux
        real, dimension(1:n_var)     :: DelJminusFlux
        real, dimension(1:n_var)     :: DelKminusFlux
        real, dimension(1:n_var)     :: LOldIminusFlux
        real, dimension(1:n_var)     :: LOldJminusFlux
        real, dimension(1:n_var)     :: LOldKminusFlux
        real, dimension(1:n_var)     :: LNewIminusFlux
        real, dimension(1:n_var)     :: LNewJminusFlux
        real, dimension(1:n_var)     :: LNewKminusFlux
        real, dimension(1:n_var)     :: LDelIminusFlux
        real, dimension(1:n_var)     :: LDelJminusFlux
        real, dimension(1:n_var)     :: LDelKminusFlux
        real, dimension(1:n_var)     :: UOldIminusFlux
        real, dimension(1:n_var)     :: UOldJminusFlux
        real, dimension(1:n_var)     :: UOldKminusFlux
        real, dimension(1:n_var)     :: UNewIminusFlux
        real, dimension(1:n_var)     :: UNewJminusFlux
        real, dimension(1:n_var)     :: UNewKminusFlux
        real, dimension(1:n_var)     :: UDelIminusFlux
        real, dimension(1:n_var)     :: UDelJminusFlux
        real, dimension(1:n_var)     :: UDelKminusFlux
        real, dimension(1:6)     :: LambdaTimesArea
        real, dimension(1:n_var)     :: Q0 ! state at cell
        real, dimension(1:n_var)     :: Q1 ! state at neighbours 
        real, dimension(1:n_var)     :: Q2
        real, dimension(1:n_var)     :: Q3
        real, dimension(1:n_var)     :: Q4
        real, dimension(1:n_var)     :: Q5
        real, dimension(1:n_var)     :: Q6
        real, dimension(1:n_var)     :: DQ0! change in state
        real, dimension(1:n_var)     :: DQ1
        real, dimension(1:n_var)     :: DQ2
        real, dimension(1:n_var)     :: DQ3
        real, dimension(1:n_var)     :: DQ4
        real, dimension(1:n_var)     :: DQ5
        real, dimension(1:n_var)     :: DQ6
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
        real :: beta
        real :: fv1
        real :: fv2
        real :: fw
        real :: g
        real :: Scap
        real :: r
        real :: S_v
        real :: D_v
        real :: P_v
        real :: lamda
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
        real :: nu_t
        real :: glim
        real :: g_6
        real :: gamma_BC
        real :: dfv1
        real :: dfv2
        real :: dfw
        real :: dShat
        real :: dr
        real :: dg



        !intialize delQ
        delQ = 0.0
        delQstar = 0.0
        delQk12 = 0.0

        !Diagonal forward sweep
        if (time_step_accuracy == "dplusgs" .or. time_step_accuracy == "hlusgs") then
          do k=1,kmx-1
            do j=1,jmx-1
              do i=1,imx-1
                C0  = CellCenter(i  ,j  ,k  ,:)
                C1  = CellCenter(i-1,j  ,k  ,:)
                C2  = CellCenter(i  ,j-1,k  ,:)
                C3  = CellCenter(i  ,j  ,k-1,:)
                C4  = CellCenter(i+1,j  ,k  ,:)
                C5  = CellCenter(i  ,j+1,k  ,:)
                C6  = CellCenter(i  ,j  ,k+1,:)

                Q0  = qp(i  , j  , k  , 1:n_var)
                Q1  = qp(i-1, j  , k  , 1:n_var)
                Q2  = qp(i  , j-1, k  , 1:n_var)
                Q3  = qp(i  , j  , k-1, 1:n_var)
                Q4  = qp(i+1, j  , k  , 1:n_var)
                Q5  = qp(i  , j+1, k  , 1:n_var)
                Q6  = qp(i  , j  , k+1, 1:n_var)


                Flist1(1) =   xA(i,j,k)
                Flist1(2) = -xnx(i,j,k)
                Flist1(3) = -xny(i,j,k)
                Flist1(4) = -xnz(i,j,k)
                Flist1(5) = 0.5*(volume(i-1, j  , k  ) + volume(i,j,k))
                Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
                Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))
                Flist1(8) = 0.5*(blendF(i-1, j  , k  ) + blendF(i,j,k))

                Flist2(1) =   yA(i,j,k)
                Flist2(2) = -ynx(i,j,k)
                Flist2(3) = -yny(i,j,k)
                Flist2(4) = -ynz(i,j,k)
                Flist2(5) = 0.5*(volume(i  , j-1, k  ) + volume(i,j,k))
                Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
                Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))
                Flist2(8) = 0.5*(blendF(i  , j-1, k  ) + blendF(i,j,k))

                Flist3(1) =   zA(i,j,k)
                Flist3(2) = -znx(i,j,k)
                Flist3(3) = -zny(i,j,k)
                Flist3(4) = -znz(i,j,k)
                Flist3(5) = 0.5*(volume(i  , j  , k-1) + volume(i,j,k))
                Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
                Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))
                Flist3(8) = 0.5*(blendF(i  , j  , k-1) + blendF(i,j,k))

                Flist4(1) =   xA(i+1,j,k)
                Flist4(2) = +xnx(i+1,j,k)
                Flist4(3) = +xny(i+1,j,k)
                Flist4(4) = +xnz(i+1,j,k)
                Flist4(5) = 0.5*(volume(i+1, j  , k  ) + volume(i,j,k))
                Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
                Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))
                Flist4(8) = 0.5*(blendF(i+1, j  , k  ) + blendF(i,j,k))

                Flist5(1) =   yA(i,j+1,k)
                Flist5(2) = +ynx(i,j+1,k)
                Flist5(3) = +yny(i,j+1,k)
                Flist5(4) = +ynz(i,j+1,k)
                Flist5(5) = 0.5*(volume(i  , j+1, k  ) + volume(i,j,k))
                Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
                Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))
                Flist5(8) = 0.5*(blendF(i  , j+1, k  ) + blendF(i,j,k))

                Flist6(1) =   zA(i,j,k+1)
                Flist6(2) = +znx(i,j,k+1)
                Flist6(3) = +zny(i,j,k+1)
                Flist6(4) = +znz(i,j,k+1)
                Flist6(5) = 0.5*(volume(i  , j  , k+1) + volume(i,j,k))
                Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
                Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))
                Flist6(8) = 0.5*(blendF(i  , j  , k+1) + blendF(i,j,k))


                LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0)
                LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0)
                LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0)
                LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0)
                LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0)
                LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0)

                D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)

                select case(trim(turbulence))
                  case("none")
                    !do nothing
                    continue
                  case("sa")
                    ! -- source term derivatives -- !
                    Omega = sqrt( ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
                                 + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
                                 + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
                                  )&
                             )
                    dist_i = dist(i,j,k)
                    dist_i_2 = dist_i*dist_i
                    k2 = kappa_sa*kappa_sa
                    nu   = mu(i,j,k)/density(i,j,k)
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
                    D(6) = D(6) - cb1*(Q0(6)*dShat+Shat)*Volume(i,j,k)
                    ! ___ Destruction term___ !
                    r    = min(Q0(6)*inv_Shat*inv_k2_d2, 10.0)
                    g    = r + cw2*((r**6) - r)
                    g_6  = g**6
                    glim = ((1.0+cw3_6)/(g_6+cw3_6))**(1.0/6.0)
                    fw   = g*glim
                    dr = (Shat-Q0(6)*dShat)*inv_Shat*inv_Shat*inv_k2_d2
                    dg = dr*(1.0+cw2*(6.0*(r**5)-1.0))
                    dfw= dg*glim*(1.0-g_6/(g_6+cw3_6))
                    D(6) = D(6)+cw1*(dfw*Q0(6) + 2*fw)*Q0(6)/dist_i_2*volume(i,j,k)
                    ! --  end of source term -- !
                  case("sst", 'sst2003')
                    beta = sst_F1(i,j,k)*beta1 + (1.0-sst_F1(i,j,k))*beta2
                    D(6) = (D(6) + ( bstar*qp(i,j,k,7))*volume(i,j,k))
                    D(7) = (D(7) + 2.0*beta*qp(i,j,k,7)*volume(i,j,k))
                  case("kkl")
                    D(6) = D(6) + (2.5*(cmu**(0.75))*Q0(1)*(Q0(6)**(1.5))*volume(i,j,k)/Q0(7))
                    D(6) = D(6) + (2*mmu(i,j,k)*volume(i,j,k)/(dist(i,j,k)**2))
                    D(7) = D(7) + (6*mmu(i,j,k)*volume(i,j,k)/(dist(i,j,k)**2))
                  case Default
                    ! do nothing
                    continue
                end select 
                Diagonal(i,j,k,:) = D

                delQ(i,j,k,1:n_var) = -residue(i,j,k,1:n_var)/D

              end do
            end do
          end do
        end if


        do SubIter = 1,SGSIter
          ! Parallel data communication
          ! call apply_interface( Change_state, layers_of_ghost_cells)
          ! Layer should be 1 as delQ is not allocated further than 1 ghost layer
          call apply_interface(delQ, 1)

          !forward sweep
          do k=1,kmx-1
            do j=1,jmx-1
              do i=1,imx-1
                C0  = CellCenter(i  ,j  ,k  ,:)
                C1  = CellCenter(i-1,j  ,k  ,:)
                C2  = CellCenter(i  ,j-1,k  ,:)
                C3  = CellCenter(i  ,j  ,k-1,:)
                C4  = CellCenter(i+1,j  ,k  ,:)
                C5  = CellCenter(i  ,j+1,k  ,:)
                C6  = CellCenter(i  ,j  ,k+1,:)

                Q0  = qp(i  , j  , k  , 1:n_var)
                Q1  = qp(i-1, j  , k  , 1:n_var)
                Q2  = qp(i  , j-1, k  , 1:n_var)
                Q3  = qp(i  , j  , k-1, 1:n_var)
                Q4  = qp(i+1, j  , k  , 1:n_var)
                Q5  = qp(i  , j+1, k  , 1:n_var)
                Q6  = qp(i  , j  , k+1, 1:n_var)

                DQ0 = 0.0
                DQ1 = delQ(i-1, j  , k  , 1:n_var)
                DQ2 = delQ(i  , j-1, k  , 1:n_var)
                DQ3 = delQ(i  , j  , k-1, 1:n_var)

                Flist1(1) =   xA(i,j,k)
                Flist1(2) = -xnx(i,j,k)
                Flist1(3) = -xny(i,j,k)
                Flist1(4) = -xnz(i,j,k)
                Flist1(5) = 0.5*(volume(i-1, j  , k  ) + volume(i,j,k))
                Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
                Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))
                Flist1(8) = 0.5*(blendF(i-1, j  , k  ) + blendF(i,j,k))

                Flist2(1) =   yA(i,j,k)
                Flist2(2) = -ynx(i,j,k)
                Flist2(3) = -yny(i,j,k)
                Flist2(4) = -ynz(i,j,k)
                Flist2(5) = 0.5*(volume(i  , j-1, k  ) + volume(i,j,k))
                Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
                Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))
                Flist2(8) = 0.5*(blendF(i  , j-1, k  ) + blendF(i,j,k))

                Flist3(1) =   zA(i,j,k)
                Flist3(2) = -znx(i,j,k)
                Flist3(3) = -zny(i,j,k)
                Flist3(4) = -znz(i,j,k)
                Flist3(5) = 0.5*(volume(i  , j  , k-1) + volume(i,j,k))
                Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
                Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))
                Flist3(8) = 0.5*(blendF(i  , j  , k-1) + blendF(i,j,k))

                Flist4(1) =   xA(i+1,j,k)
                Flist4(2) = +xnx(i+1,j,k)
                Flist4(3) = +xny(i+1,j,k)
                Flist4(4) = +xnz(i+1,j,k)
                Flist4(5) = 0.5*(volume(i+1, j  , k  ) + volume(i,j,k))
                Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
                Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))
                Flist4(8) = 0.5*(blendF(i+1, j  , k  ) + blendF(i,j,k))

                Flist5(1) =   yA(i,j+1,k)
                Flist5(2) = +ynx(i,j+1,k)
                Flist5(3) = +yny(i,j+1,k)
                Flist5(4) = +ynz(i,j+1,k)
                Flist5(5) = 0.5*(volume(i  , j+1, k  ) + volume(i,j,k))
                Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
                Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))
                Flist5(8) = 0.5*(blendF(i  , j+1, k  ) + blendF(i,j,k))

                Flist6(1) =   zA(i,j,k+1)
                Flist6(2) = +znx(i,j,k+1)
                Flist6(3) = +zny(i,j,k+1)
                Flist6(4) = +znz(i,j,k+1)
                Flist6(5) = 0.5*(volume(i  , j  , k+1) + volume(i,j,k))
                Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
                Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))
                Flist6(8) = 0.5*(blendF(i  , j  , k+1) + blendF(i,j,k))

                NewIminusFlux     = Flux(Q1, Q0, DQ1, Flist1)
                NewJminusFlux     = Flux(Q2, Q0, DQ2, Flist2)
                NewKminusFlux     = Flux(Q3, Q0, DQ3, Flist3)
                OldIminusFlux     = Flux(Q1, Q0, DQ0, Flist1)
                OldJminusFlux     = Flux(Q2, Q0, DQ0, Flist2)
                OldKminusFlux     = Flux(Q3, Q0, DQ0, Flist3)

                LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0)
                LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0)
                LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0)
                LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0)
                LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0)
                LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0)


                ! multiply above flux with area to get correct values
                DelIminusFlux =  NewIminusFlux - OldIminusFlux
                DelJminusFlux =  NewJminusFlux - OldJminusFlux
                DelKminusFlux =  NewKminusFlux - OldKminusFlux


                !D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
                D = Diagonal(i,j,k,:) 
                !storing D in Iflux array for backward sweep
                !F_p(i,j,k,1) = D

                deltaU(1:n_var) = -residue(i,j,k,1:n_var) &
                  - 0.5*((DelIminusFlux - LambdaTimesArea(1)*delQ(i-1,j,k,1:n_var)) &
                       + (DelJminusFlux - LambdaTimesArea(2)*delQ(i,j-1,k,1:n_var)) &
                       + (DelKminusFlux - LambdaTimesArea(3)*delQ(i,j,k-1,1:n_var)) &
                       +  delQk12(i,j,k,1:n_var))

                delQ(i,j,k,1:n_var) = deltaU(1:n_var)/D

                delQk12(i,j,k,1:n_var) = &
                        ((DelIminusFlux - LambdaTimesArea(1)*delQ(i-1,j,k,1:n_var)) &
                       + (DelJminusFlux - LambdaTimesArea(2)*delQ(i,j-1,k,1:n_var)) &
                       + (DelKminusFlux - LambdaTimesArea(3)*delQ(i,j,k-1,1:n_var)) )
              end do
            end do
          end do

          !backward sweep
              do i=imx-1,1,-1
            do j=jmx-1,1,-1
          do k=kmx-1,1,-1
                C0  = CellCenter(i  ,j  ,k  ,:)
                C1  = CellCenter(i-1,j  ,k  ,:)
                C2  = CellCenter(i  ,j-1,k  ,:)
                C3  = CellCenter(i  ,j  ,k-1,:)
                C4  = CellCenter(i+1,j  ,k  ,:)
                C5  = CellCenter(i  ,j+1,k  ,:)
                C6  = CellCenter(i  ,j  ,k+1,:)

                Q0  = qp(i  , j  , k  , 1:n_var)
                Q1  = qp(i-1, j  , k  , 1:n_var)
                Q2  = qp(i  , j-1, k  , 1:n_var)
                Q3  = qp(i  , j  , k-1, 1:n_var)
                Q4  = qp(i+1, j  , k  , 1:n_var)
                Q5  = qp(i  , j+1, k  , 1:n_var)
                Q6  = qp(i  , j  , k+1, 1:n_var)

                DQ0 = 0.0
                DQ4 = delQ(i+1, j  , k  , 1:n_var)
                DQ5 = delQ(i  , j+1, k  , 1:n_var)
                DQ6 = delQ(i  , j  , k+1, 1:n_var)

                Flist1(1) =   xA(i,j,k)
                Flist1(2) = -xnx(i,j,k)
                Flist1(3) = -xny(i,j,k)
                Flist1(4) = -xnz(i,j,k)
                Flist1(5) = 0.5*(volume(i-1, j  , k  ) + volume(i,j,k))
                Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
                Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))
                Flist1(8) = 0.5*(blendF(i-1, j  , k  ) + blendF(i,j,k))

                Flist2(1) =   yA(i,j,k)
                Flist2(2) = -ynx(i,j,k)
                Flist2(3) = -yny(i,j,k)
                Flist2(4) = -ynz(i,j,k)
                Flist2(5) = 0.5*(volume(i  , j-1, k  ) + volume(i,j,k))
                Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
                Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))
                Flist2(8) = 0.5*(blendF(i  , j-1, k  ) + blendF(i,j,k))

                Flist3(1) =   zA(i,j,k)
                Flist3(2) = -znx(i,j,k)
                Flist3(3) = -zny(i,j,k)
                Flist3(4) = -znz(i,j,k)
                Flist3(5) = 0.5*(volume(i  , j  , k-1) + volume(i,j,k))
                Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
                Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))
                Flist3(8) = 0.5*(blendF(i  , j  , k-1) + blendF(i,j,k))

                Flist4(1) =   xA(i+1,j,k)
                Flist4(2) = +xnx(i+1,j,k)
                Flist4(3) = +xny(i+1,j,k)
                Flist4(4) = +xnz(i+1,j,k)
                Flist4(5) = 0.5*(volume(i+1, j  , k  ) + volume(i,j,k))
                Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
                Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))
                Flist4(8) = 0.5*(blendF(i+1, j  , k  ) + blendF(i,j,k))

                Flist5(1) =   yA(i,j+1,k)
                Flist5(2) = +ynx(i,j+1,k)
                Flist5(3) = +yny(i,j+1,k)
                Flist5(4) = +ynz(i,j+1,k)
                Flist5(5) = 0.5*(volume(i  , j+1, k  ) + volume(i,j,k))
                Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
                Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))
                Flist5(8) = 0.5*(blendF(i  , j+1, k  ) + blendF(i,j,k))

                Flist6(1) =   zA(i,j,k+1)
                Flist6(2) = +znx(i,j,k+1)
                Flist6(3) = +zny(i,j,k+1)
                Flist6(4) = +znz(i,j,k+1)
                Flist6(5) = 0.5*(volume(i  , j  , k+1) + volume(i,j,k))
                Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
                Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))
                Flist6(8) = 0.5*(blendF(i  , j  , k+1) + blendF(i,j,k))

                NewIminusFlux     = Flux(Q4, Q0, DQ4, Flist4)
                NewJminusFlux     = Flux(Q5, Q0, DQ5, Flist5)
                NewKminusFlux     = Flux(Q6, Q0, DQ6, Flist6)
                OldIminusFlux     = Flux(Q4, Q0, DQ0, Flist4)
                OldJminusFlux     = Flux(Q5, Q0, DQ0, Flist5)
                OldKminusFlux     = Flux(Q6, Q0, DQ0, Flist6)

                LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0)
                LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0)
                LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0)
                LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0)
                LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0)
                LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0)


                ! multiply above flux with area to get correct values
                DelIminusFlux =  NewIminusFlux - OldIminusFlux
                DelJminusFlux =  NewJminusFlux - OldJminusFlux
                DelKminusFlux =  NewKminusFlux - OldKminusFlux

                !D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
                D = Diagonal(i,j,k,:) 

                delQ(i,j,k,1:n_var) = (-residue(i,j,k,1:n_var) &
                  - 0.5*((DelIminusFlux - LambdaTimesArea(4)*delQ(i+1,j,k,1:n_var)) &
                       + (DelJminusFlux - LambdaTimesArea(5)*delQ(i,j+1,k,1:n_var)) &
                       + (DelKminusFlux - LambdaTimesArea(6)*delQ(i,j,k+1,1:n_var)) &
                       + delQk12(i,j,k,1:n_var)))/D

                delQk12(i,j,k,1:n_var) = &
                        ((DelIminusFlux - LambdaTimesArea(4)*delQ(i+1,j,k,1:n_var)) &
                       + (DelJminusFlux - LambdaTimesArea(5)*delQ(i,j+1,k,1:n_var)) &
                       + (DelKminusFlux - LambdaTimesArea(6)*delQ(i,j,k+1,1:n_var)) )
              end do
            end do
          end do
        end do


        call UpdateState(qp,delQ,n_var)

    end subroutine update_laminar_variables


    function Flux(ql, qr, du, inputs)
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
      real, dimension(1:n_var)             :: U ! conservative variables
      real, dimension(1:n_var)             :: W ! new primitive variables
      real, dimension(1:n_var)             :: P ! primitive variables of right cell
      real, dimension(1:n_var)             :: dQdx ! derivative of primitive variables wrt x
      real, dimension(1:n_var)             :: dQdy ! derivative of primitive variables wrt y
      real, dimension(1:n_var)             :: dQdz ! derivative of primitive variables wrt z

      !for extraction of the inputs
      real :: area
      real :: nx
      real :: ny
      real :: nz
      real :: volume
      real :: mmu
      real :: tmu
      real :: F1


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
      real    :: muCap
      real    :: sigma_k_sst
      real    :: sigma_w_sst

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
      U(1)       =   ql(1)
      U(2:n_var) =   ql(1) * ql(2:n_var)
      U(5)       = ( ql(5) / (gm-1.0) ) + ( 0.5 * ql(1) * sum(ql(2:4)**2) )

      U(1:n_var) = U(1:n_var) + du(1:n_var)
      

      W(1)       =   U(1)
      W(2) =   U(2) / U(1)
      W(3) =   U(3) / U(1)
      W(4) =   U(4) / U(1)
      W(2:n_var) =   U(2:n_var) / U(1)
      W(5)       = (gm-1.0) * ( U(5) - ( 0.5 * SUM(U(2:4)**2) / U(1) ) )
      w(5:n_var) = max(w(5:n_var),1.0e-10)

      FaceNormalVelocity = (W(2) * nx) + (W(3) * ny) + (W(4) * nz)
      uface = 0.5 * ( W(2) + P(2) )
      vface = 0.5 * ( W(3) + P(3) )
      wface = 0.5 * ( W(4) + P(4) )


      Flux(1)       =   W(1) * FaceNormalVelocity
      Flux(2:n_var) = W(2:n_var) * Flux(1)
      Flux(2)       = Flux(2) + ( W(5) * nx )
      Flux(3)       = Flux(3) + ( W(5) * ny )
      Flux(4)       = Flux(4) + ( W(5) * nz )

      HalfRhoUsquare = 0.5 * W(1) * ( W(2)*W(2) + W(3)*W(3) + W(4)*W(4) )
      RhoHt          = ( (gm/(gm-1.0)) * W(5) ) + HalfRhoUsquare
      Flux(5)        = RhoHt * FaceNormalVelocity


      ! viscous terms
      mu = mmu + tmu
      ! Storing Temperature in place of pressure
      !T1     =    W(5) / ( W(1) * R_gas )
      !T2     =    P(5) / ( P(1) * R_gas )
!      W(5)   =   W(5) / ( W(1) * R_gas )
!      P(5)   =   P(5) / ( P(1) * R_gas )
!      dTdx   =  ( T2   - T1   ) * nx * Area / Volume
!      dTdy   =  ( T2   - T1   ) * ny * Area / Volume
!      dTdz   =  ( T2   - T1   ) * nz * Area / Volume
!      dudx   =  ( P(2) - W(2) ) * nx * Area / Volume
!      dudy   =  ( P(2) - W(2) ) * ny * Area / Volume
!      dudz   =  ( P(2) - W(2) ) * nz * Area / Volume
!      dvdx   =  ( P(3) - W(3) ) * nx * Area / Volume
!      dvdy   =  ( P(3) - W(3) ) * ny * Area / Volume
!      dvdz   =  ( P(3) - W(3) ) * nz * Area / Volume
!      dwdx   =  ( P(4) - W(4) ) * nx * Area / Volume
!      dwdy   =  ( P(4) - W(4) ) * ny * Area / Volume
!      dwdz   =  ( P(4) - W(4) ) * nz * Area / Volume
      dQdx   =  ( P - W ) * nx * Area / Volume
      dQdy   =  ( P - W ) * ny * Area / Volume
      dQdz   =  ( P - W ) * nz * Area / Volume
      !finding Dt/Dx_i
      dQdx(5) = (dQdx(5) - (W(5)/W(1))*dQdx(1))/(R_gas*W(1))
      dQdy(5) = (dQdy(5) - (W(5)/W(1))*dQdy(1))/(R_gas*W(1))
      dQdz(5) = (dQdz(5) - (W(5)/W(1))*dQdz(1))/(R_gas*W(1))

      trace = dQdx(2) + dQdy(3) + dQdz(4)
      Tauxx =  2. * mu * (dQdx(2) - trace/3.0)
      Tauyy =  2. * mu * (dQdy(3) - trace/3.0)
      Tauzz =  2. * mu * (dQdz(4) - trace/3.0)
      Tauxy = mu * (dQdx(3) + dQdy(2))
      Tauxz = mu * (dQdx(4) + dQdz(2))
      Tauyz = mu * (dQdy(4) + dQdz(3))

      K_heat = ( mmu / Pr  + tmu/tpr) * gm * R_gas / ( gm - 1.0 )
      Qx = K_heat*dQdx(5)
      Qy = K_heat*dQdy(5)
      Qz = K_heat*dQdz(5)

      Flux(2) = Flux(2) - ( Tauxx * nx + Tauxy * ny + Tauxz * nz )
      Flux(3) = Flux(3) - ( Tauxy * nx + Tauyy * ny + Tauyz * nz )
      Flux(4) = Flux(4) - ( Tauxz * nx + Tauyz * ny + Tauzz * nz )
      Flux(5) = Flux(5) - ( Tauxx * uface + Tauxy * vface + Tauxz * wface + Qx ) * nx
      Flux(5) = Flux(5) - ( Tauxy * uface + Tauyy * vface + Tauyz * wface + Qy ) * ny
      Flux(5) = Flux(5) - ( Tauxz * uface + Tauyz * vface + Tauzz * wface + Qz ) * nz



      select case(trim(Turbulence))

        case("none")
          ! do nothing 
          continue
        case("sa")
          muCap = 0.25*(P(1)+W(1))*(P(6) + W(6))
          tmu = 0.5*(W(6) + P(6))
          Flux(6) = Flux(6) - (mmu + muCap)*(dQdx(6)*nx + dQdy(6)*ny + dQdz(6)*nz)/sigma_sa
        case("sst", 'sst2003')
          sigma_k_sst = sigma_k1*F1 + sigma_k2*(1.0 - F1) 
          sigma_w_sst = sigma_w1*F1 + sigma_w2*(1.0 - F1) 
          Flux(6) = Flux(6) - (mmu + sigma_k_sst*tmu)*(dQdx(6)*nx + dQdy(6)*ny + dQdz(6)*nz)
          Flux(7) = Flux(7) - (mmu + sigma_w_sst*tmu)*(dQdx(7)*nx + dQdy(7)*ny + dQdz(7)*nz)
        case("kkl")
          Flux(6) = Flux(6) - (mmu + sigma_k * tmu)*(dQdx(6)*nx + dQdy(6)*ny + dQdz(6)*nz)
          Flux(7) = Flux(7) - (mmu + sigma_phi*tmu)*(dQdx(7)*nx + dQdy(7)*ny + dQdz(7)*nz)
        case DEFAULT
          ! do nothing
          continue
      end select



      Flux    = Flux * Area

    end function Flux 


    function SpectralRadius(ql, qr, inputs, c1, c2)
      implicit none
      real, dimension(1:n_var), intent(in) :: ql
      real, dimension(1:n_var), intent(in) :: qr
      real, dimension(1:7)    , intent(in) :: inputs
      real, dimension(1:3)    , intent(in) :: c1
      real, dimension(1:3)    , intent(in) :: c2

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
      vis = gm * (mm/pr + tm/tpr) / ( rho * distance )
      SpectralRadius = ( NormalSpeed + SpeedOfSound + vis) * Area

    end function SpectralRadius


    subroutine update_SST_variables()
      implicit none
      integer :: i,j,k,SubIter
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



        !intialize delQ
        delQ = 0.0
        delQK12 = 0.0

        do SubIter = 1,SGSIter
          !call apply_interface(delQ, 1)
        !forward sweep
        do k=1,kmx-1
          do j=1,jmx-1
            do i=1,imx-1
              C0  = CellCenter(i  ,j  ,k  ,:)
              C1  = CellCenter(i-1,j  ,k  ,:)
              C2  = CellCenter(i  ,j-1,k  ,:)
              C3  = CellCenter(i  ,j  ,k-1,:)
              C4  = CellCenter(i+1,j  ,k  ,:)
              C5  = CellCenter(i  ,j+1,k  ,:)
              C6  = CellCenter(i  ,j  ,k+1,:)

              Q0  = qp(i  , j  , k  , 1:7)
              Q1  = qp(i-1, j  , k  , 1:7)
              Q2  = qp(i  , j-1, k  , 1:7)
              Q3  = qp(i  , j  , k-1, 1:7)
              Q4  = qp(i+1, j  , k  , 1:7)
              Q5  = qp(i  , j+1, k  , 1:7)
              Q6  = qp(i  , j  , k+1, 1:7)

              DQ0 = 0.0
              DQ1 = delQ(i-1, j  , k  , 1:7)
              DQ2 = delQ(i  , j-1, k  , 1:7)
              DQ3 = delQ(i  , j  , k-1, 1:7)

              Flist1(1) =   xA(i,j,k)
              Flist1(2) = -xnx(i,j,k)
              Flist1(3) = -xny(i,j,k)
              Flist1(4) = -xnz(i,j,k)
              Flist1(5) = 0.5*(volume(i-1, j  , k  ) + volume(i,j,k))
              Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))
              Flist1(8) = 0.5*(sst_F1(i-1, j  , k  ) + sst_F1(i,j,k))

              Flist2(1) =   yA(i,j,k)
              Flist2(2) = -ynx(i,j,k)
              Flist2(3) = -yny(i,j,k)
              Flist2(4) = -ynz(i,j,k)
              Flist2(5) = 0.5*(volume(i  , j-1, k  ) + volume(i,j,k))
              Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))
              Flist2(8) = 0.5*(sst_F1(i  , j-1, k  ) + sst_F1(i,j,k))

              Flist3(1) =   zA(i,j,k)
              Flist3(2) = -znx(i,j,k)
              Flist3(3) = -zny(i,j,k)
              Flist3(4) = -znz(i,j,k)
              Flist3(5) = 0.5*(volume(i  , j  , k-1) + volume(i,j,k))
              Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))
              Flist3(8) = 0.5*(sst_F1(i  , j  , k-1) + sst_F1(i,j,k))

              Flist4(1) =   xA(i+1,j,k)
              Flist4(2) = +xnx(i+1,j,k)
              Flist4(3) = +xny(i+1,j,k)
              Flist4(4) = +xnz(i+1,j,k)
              Flist4(5) = 0.5*(volume(i+1, j  , k  ) + volume(i,j,k))
              Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))
              Flist4(8) = 0.5*(sst_F1(i+1, j  , k  ) + sst_F1(i,j,k))

              Flist5(1) =   yA(i,j+1,k)
              Flist5(2) = +ynx(i,j+1,k)
              Flist5(3) = +yny(i,j+1,k)
              Flist5(4) = +ynz(i,j+1,k)
              Flist5(5) = 0.5*(volume(i  , j+1, k  ) + volume(i,j,k))
              Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))
              Flist5(8) = 0.5*(sst_F1(i  , j+1, k  ) + sst_F1(i,j,k))

              Flist6(1) =   zA(i,j,k+1)
              Flist6(2) = +znx(i,j,k+1)
              Flist6(3) = +zny(i,j,k+1)
              Flist6(4) = +znz(i,j,k+1)
              Flist6(5) = 0.5*(volume(i  , j  , k+1) + volume(i,j,k))
              Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))
              Flist6(8) = 0.5*(sst_F1(i  , j  , k+1) + sst_F1(i,j,k))

              NewIminusFlux     = SSTFlux(Q1, Q0, DQ1, Flist1)
              NewJminusFlux     = SSTFlux(Q2, Q0, DQ2, Flist2)
              NewKminusFlux     = SSTFlux(Q3, Q0, DQ3, Flist3)
              OldIminusFlux     = SSTFlux(Q1, Q0, DQ0, Flist1)
              OldJminusFlux     = SSTFlux(Q2, Q0, DQ0, Flist2)
              OldKminusFlux     = SSTFlux(Q3, Q0, DQ0, Flist3)

              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0)


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux


              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
              beta = sst_F1(i,j,k)*beta1 + (1.0-sst_F1(i,j,k))*beta2
              D(6) = (D(6) + bstar*qp(i,j,k,7)*volume(i,j,k))
              D(7) = (D(7) + 2.0*beta*qp(i,j,k,7)*volume(i,j,k))
              !storing D in Iflux array for backward sweep
              !F_p(i,j,k,1) = D

                deltaU(1:n_var) = -residue(i,j,k,1:n_var) &
                  - 0.5*((DelIminusFlux - LambdaTimesArea(1)*delQ(i-1,j,k,1:n_var)) &
                       + (DelJminusFlux - LambdaTimesArea(2)*delQ(i,j-1,k,1:n_var)) &
                       + (DelKminusFlux - LambdaTimesArea(3)*delQ(i,j,k-1,1:n_var)) &
                       +  delQk12(i,j,k,1:n_var))

                delQ(i,j,k,1:n_var) = deltaU(1:n_var)/D

                delQk12(i,j,k,1:n_var) = &
                        ((DelIminusFlux - LambdaTimesArea(1)*delQ(i-1,j,k,1:n_var)) &
                       + (DelJminusFlux - LambdaTimesArea(2)*delQ(i,j-1,k,1:n_var)) &
                       + (DelKminusFlux - LambdaTimesArea(3)*delQ(i,j,k-1,1:n_var)) )
            end do
          end do
        end do

        !delQ=0.0
        !backward sweep
            do i=imx-1,1,-1
          do j=jmx-1,1,-1
        do k=kmx-1,1,-1
              C0  = CellCenter(i  ,j  ,k  ,:)
              C1  = CellCenter(i-1,j  ,k  ,:)
              C2  = CellCenter(i  ,j-1,k  ,:)
              C3  = CellCenter(i  ,j  ,k-1,:)
              C4  = CellCenter(i+1,j  ,k  ,:)
              C5  = CellCenter(i  ,j+1,k  ,:)
              C6  = CellCenter(i  ,j  ,k+1,:)

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

              Flist1(1) =   xA(i,j,k)
              Flist1(2) = -xnx(i,j,k)
              Flist1(3) = -xny(i,j,k)
              Flist1(4) = -xnz(i,j,k)
              Flist1(5) = 0.5*(volume(i-1, j  , k  ) + volume(i,j,k))
              Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))
              Flist1(8) = 0.5*(sst_F1(i-1, j  , k  ) + sst_F1(i,j,k))

              Flist2(1) =   yA(i,j,k)
              Flist2(2) = -ynx(i,j,k)
              Flist2(3) = -yny(i,j,k)
              Flist2(4) = -ynz(i,j,k)
              Flist2(5) = 0.5*(volume(i  , j-1, k  ) + volume(i,j,k))
              Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))
              Flist2(8) = 0.5*(sst_F1(i  , j-1, k  ) + sst_F1(i,j,k))

              Flist3(1) =   zA(i,j,k)
              Flist3(2) = -znx(i,j,k)
              Flist3(3) = -zny(i,j,k)
              Flist3(4) = -znz(i,j,k)
              Flist3(5) = 0.5*(volume(i  , j  , k-1) + volume(i,j,k))
              Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))
              Flist3(8) = 0.5*(sst_F1(i  , j  , k-1) + sst_F1(i,j,k))

              Flist4(1) =   xA(i+1,j,k)
              Flist4(2) = +xnx(i+1,j,k)
              Flist4(3) = +xny(i+1,j,k)
              Flist4(4) = +xnz(i+1,j,k)
              Flist4(5) = 0.5*(volume(i+1, j  , k  ) + volume(i,j,k))
              Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))
              Flist4(8) = 0.5*(sst_F1(i+1, j  , k  ) + sst_F1(i,j,k))

              Flist5(1) =   yA(i,j+1,k)
              Flist5(2) = +ynx(i,j+1,k)
              Flist5(3) = +yny(i,j+1,k)
              Flist5(4) = +ynz(i,j+1,k)
              Flist5(5) = 0.5*(volume(i  , j+1, k  ) + volume(i,j,k))
              Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))
              Flist5(8) = 0.5*(sst_F1(i  , j+1, k  ) + sst_F1(i,j,k))

              Flist6(1) =   zA(i,j,k+1)
              Flist6(2) = +znx(i,j,k+1)
              Flist6(3) = +zny(i,j,k+1)
              Flist6(4) = +znz(i,j,k+1)
              Flist6(5) = 0.5*(volume(i  , j  , k+1) + volume(i,j,k))
              Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))
              Flist6(8) = 0.5*(sst_F1(i  , j  , k+1) + sst_F1(i,j,k))

              NewIminusFlux     = SSTFlux(Q4, Q0, DQ4, Flist4)
              NewJminusFlux     = SSTFlux(Q5, Q0, DQ5, Flist5)
              NewKminusFlux     = SSTFlux(Q6, Q0, DQ6, Flist6)
              OldIminusFlux     = SSTFlux(Q4, Q0, DQ0, Flist4)
              OldJminusFlux     = SSTFlux(Q5, Q0, DQ0, Flist5)
              OldKminusFlux     = SSTFlux(Q6, Q0, DQ0, Flist6)

              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0)


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux

              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
              beta = sst_F1(i,j,k)*beta1 + (1.0-sst_F1(i,j,k))*beta2
              D(6) = (D(6) + bstar*qp(i,j,k,7)*volume(i,j,k))
              D(7) = (D(7) + 2.0*beta*qp(i,j,k,7)*volume(i,j,k))

                delQ(i,j,k,1:n_var) = (-residue(i,j,k,1:n_var) &
                  - 0.5*((DelIminusFlux - LambdaTimesArea(4)*delQ(i+1,j,k,1:n_var)) &
                       + (DelJminusFlux - LambdaTimesArea(5)*delQ(i,j+1,k,1:n_var)) &
                       + (DelKminusFlux - LambdaTimesArea(6)*delQ(i,j,k+1,1:n_var)) &
                       + delQk12(i,j,k,1:n_var)))/D

                delQk12(i,j,k,1:n_var) = &
                        ((DelIminusFlux - LambdaTimesArea(4)*delQ(i+1,j,k,1:n_var)) &
                       + (DelJminusFlux - LambdaTimesArea(5)*delQ(i,j+1,k,1:n_var)) &
                       + (DelKminusFlux - LambdaTimesArea(6)*delQ(i,j,k+1,1:n_var)) )

            end do
          end do
        end do
        end do
        
        do k=1,kmx-1
          do j = 1,jmx-1
            do i = 1,imx-1
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
              qp(i,j,k,6) = conservativeQ(6) / conservativeQ(1)
              qp(i,j,k,7) = conservativeQ(7) / conservativeQ(1)
              qp(i,j,k,6) = max(qp(i,j,k,6), 1.e-8)
              qp(i,j,k,7) = max(qp(i,j,k,7), 1.e-8)
            end do
          end do
        end do

    end subroutine update_SST_variables


    function SSTFlux(ql, qr, du, inputs)
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
      W(6) = max(W(6), 1e-8)
      W(7) = max(W(7), 1e-8)

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

    subroutine update_KKL_variables()
      implicit none
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



        !intialize delQ
        delQstar = 0.0

        !forward sweep
        do k=1,kmx-1
          do j=1,jmx-1
            do i=1,imx-1
              C0  = CellCenter(i  ,j  ,k  ,:)
              C1  = CellCenter(i-1,j  ,k  ,:)
              C2  = CellCenter(i  ,j-1,k  ,:)
              C3  = CellCenter(i  ,j  ,k-1,:)
              C4  = CellCenter(i+1,j  ,k  ,:)
              C5  = CellCenter(i  ,j+1,k  ,:)
              C6  = CellCenter(i  ,j  ,k+1,:)

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

              Flist1(1) =   xA(i,j,k)
              Flist1(2) = -xnx(i,j,k)
              Flist1(3) = -xny(i,j,k)
              Flist1(4) = -xnz(i,j,k)
              Flist1(5) = 0.5*(volume(i-1, j  , k  ) + volume(i,j,k))
              Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))

              Flist2(1) =   yA(i,j,k)
              Flist2(2) = -ynx(i,j,k)
              Flist2(3) = -yny(i,j,k)
              Flist2(4) = -ynz(i,j,k)
              Flist2(5) = 0.5*(volume(i  , j-1, k  ) + volume(i,j,k))
              Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))

              Flist3(1) =   zA(i,j,k)
              Flist3(2) = -znx(i,j,k)
              Flist3(3) = -zny(i,j,k)
              Flist3(4) = -znz(i,j,k)
              Flist3(5) = 0.5*(volume(i  , j  , k-1) + volume(i,j,k))
              Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))

              Flist4(1) =   xA(i+1,j,k)
              Flist4(2) = +xnx(i+1,j,k)
              Flist4(3) = +xny(i+1,j,k)
              Flist4(4) = +xnz(i+1,j,k)
              Flist4(5) = 0.5*(volume(i+1, j  , k  ) + volume(i,j,k))
              Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))

              Flist5(1) =   yA(i,j+1,k)
              Flist5(2) = +ynx(i,j+1,k)
              Flist5(3) = +yny(i,j+1,k)
              Flist5(4) = +ynz(i,j+1,k)
              Flist5(5) = 0.5*(volume(i  , j+1, k  ) + volume(i,j,k))
              Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))

              Flist6(1) =   zA(i,j,k+1)
              Flist6(2) = +znx(i,j,k+1)
              Flist6(3) = +zny(i,j,k+1)
              Flist6(4) = +znz(i,j,k+1)
              Flist6(5) = 0.5*(volume(i  , j  , k+1) + volume(i,j,k))
              Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))

              NewIminusFlux     = KKLFlux(Q1, Q0, DQ1, Flist1)
              NewJminusFlux     = KKLFlux(Q2, Q0, DQ2, Flist2)
              NewKminusFlux     = KKLFlux(Q3, Q0, DQ3, Flist3)
              OldIminusFlux     = KKLFlux(Q1, Q0, DQ0, Flist1)
              OldJminusFlux     = KKLFlux(Q2, Q0, DQ0, Flist2)
              OldKminusFlux     = KKLFlux(Q3, Q0, DQ0, Flist3)

              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0)


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux


              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
              D(6) = D(6) + (2.5*(cmu**(0.75))*Q0(1)*(Q0(6)**(1.5))*volume(i,j,k)/Q0(7))
              D(6) = D(6) + (2*mmu(i,j,k)*volume(i,j,k)/(dist(i,j,k)**2))
              D(7) = D(7) + (6*mmu(i,j,k)*volume(i,j,k)/(dist(i,j,k)**2))
              !storing D in Iflux array for backward sweep
              !F_p(i,j,k,1) = D

              deltaU(1:7) = -residue(i,j,k,1:7) &
                - 0.5*((DelIminusFlux - LambdaTimesArea(1)*delQstar(i-1,j,k,1:7)) &
                     + (DelJminusFlux - LambdaTimesArea(2)*delQstar(i,j-1,k,1:7)) &
                     + (DelKminusFlux - LambdaTimesArea(3)*delQstar(i,j,k-1,1:7)) )

              delQstar(i,j,k,1:7) = deltaU(1:7)/D
            end do
          end do
        end do

        delQ=0.0
        !backward sweep
            do i=imx-1,1,-1
          do j=jmx-1,1,-1
        do k=kmx-1,1,-1
              C0  = CellCenter(i  ,j  ,k  ,:)
              C1  = CellCenter(i-1,j  ,k  ,:)
              C2  = CellCenter(i  ,j-1,k  ,:)
              C3  = CellCenter(i  ,j  ,k-1,:)
              C4  = CellCenter(i+1,j  ,k  ,:)
              C5  = CellCenter(i  ,j+1,k  ,:)
              C6  = CellCenter(i  ,j  ,k+1,:)

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

              Flist1(1) =   xA(i,j,k)
              Flist1(2) = -xnx(i,j,k)
              Flist1(3) = -xny(i,j,k)
              Flist1(4) = -xnz(i,j,k)
              Flist1(5) = 0.5*(volume(i-1, j  , k  ) + volume(i,j,k))
              Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))

              Flist2(1) =   yA(i,j,k)
              Flist2(2) = -ynx(i,j,k)
              Flist2(3) = -yny(i,j,k)
              Flist2(4) = -ynz(i,j,k)
              Flist2(5) = 0.5*(volume(i  , j-1, k  ) + volume(i,j,k))
              Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))

              Flist3(1) =   zA(i,j,k)
              Flist3(2) = -znx(i,j,k)
              Flist3(3) = -zny(i,j,k)
              Flist3(4) = -znz(i,j,k)
              Flist3(5) = 0.5*(volume(i  , j  , k-1) + volume(i,j,k))
              Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))

              Flist4(1) =   xA(i+1,j,k)
              Flist4(2) = +xnx(i+1,j,k)
              Flist4(3) = +xny(i+1,j,k)
              Flist4(4) = +xnz(i+1,j,k)
              Flist4(5) = 0.5*(volume(i+1, j  , k  ) + volume(i,j,k))
              Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))

              Flist5(1) =   yA(i,j+1,k)
              Flist5(2) = +ynx(i,j+1,k)
              Flist5(3) = +yny(i,j+1,k)
              Flist5(4) = +ynz(i,j+1,k)
              Flist5(5) = 0.5*(volume(i  , j+1, k  ) + volume(i,j,k))
              Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))

              Flist6(1) =   zA(i,j,k+1)
              Flist6(2) = +znx(i,j,k+1)
              Flist6(3) = +zny(i,j,k+1)
              Flist6(4) = +znz(i,j,k+1)
              Flist6(5) = 0.5*(volume(i  , j  , k+1) + volume(i,j,k))
              Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))

              NewIminusFlux     = KKLFlux(Q4, Q0, DQ4, Flist4)
              NewJminusFlux     = KKLFlux(Q5, Q0, DQ5, Flist5)
              NewKminusFlux     = KKLFlux(Q6, Q0, DQ6, Flist6)
              OldIminusFlux     = KKLFlux(Q4, Q0, DQ0, Flist4)
              OldJminusFlux     = KKLFlux(Q5, Q0, DQ0, Flist5)
              OldKminusFlux     = KKLFlux(Q6, Q0, DQ0, Flist6)

              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0)


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux

              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
              D(6) = D(6) + (2.5*(cmu**(0.75))*Q0(1)*(Q0(6)**(1.5))*volume(i,j,k)/Q0(7))
              D(6) = D(6) + (2*mmu(i,j,k)*volume(i,j,k)/(dist(i,j,k)**2))
              D(7) = D(7) + (6*mmu(i,j,k)*volume(i,j,k)/(dist(i,j,k)**2))


              delQ(i,j,k,1:7) = delQstar(i,j,k,1:7) &
                - 0.5*((DelIminusFlux - LambdaTimesArea(4)*delQ(i+1,j,k,1:7)) &
                     + (DelJminusFlux - LambdaTimesArea(5)*delQ(i,j+1,k,1:7)) &
                     + (DelKminusFlux - LambdaTimesArea(6)*delQ(i,j,k+1,1:7)) )/D

            end do
          end do
        end do
        
        do k=1,kmx-1
          do j = 1,jmx-1
            do i = 1,imx-1
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
              qp(i,j,k,6) = conservativeQ(6) / conservativeQ(1)
              qp(i,j,k,7) = conservativeQ(7) / conservativeQ(1)
              qp(i,j,k,6) = max(qp(i,j,k,6), 1.e-8)
              qp(i,j,k,7) = max(qp(i,j,k,7), 1.e-8)
            end do
          end do
        end do

    end subroutine update_KKL_variables


    function KKLFlux(ql, qr, du, inputs)
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
      real, dimension(1:n_var)             :: KKLFlux
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
      real    :: dtkldx
      real    :: dtkldy
      real    :: dtkldz
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
      U(7)   =   ql(1) * ql(7)

      U(1:n_var) = U(1:n_var) + du(1:n_var)
      

      W(1)   =   U(1)
      W(2)   =   U(2) / U(1)
      W(3)   =   U(3) / U(1)
      W(4)   =   U(4) / U(1)
      W(5)   = (gm-1.0) * ( U(5) - ( 0.5 * SUM(U(2:4)**2) / U(1) ) )
      W(6)   =   U(6) / U(1)
      W(7)   =   U(7) / U(1)
      W(6) = max(W(6), 1e-8)
      W(7) = max(W(7), 1e-8)

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
      dtkldx  =  ( P(7) - W(7) ) * nx * Area / Volume
      dtkldy  =  ( P(7) - W(7) ) * ny * Area / Volume
      dtkldz  =  ( P(7) - W(7) ) * nz * Area / Volume

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
      Flux(6) = Flux(6) + (mmu + sigma_k*tmu)*(dtkdx*nx + dtkdy*ny + dtkdz*nz)
      Flux(7) = Flux(7) + (mmu + sigma_phi*tmu)*(dtkldx*nx + dtkldy*ny + dtkldz*nz)

      Flux    = Flux * Area
      KKLFlux = Flux

    end function KKLFlux 


!    subroutine update_SA_variables()
!      implicit none
!      integer :: i,j,k
!      real, dimension(1:6)     :: deltaU
!      real, dimension(1:6)     :: D
!      real, dimension(1:6)     :: conservativeQ
!      real, dimension(1:6)     :: OldIminusFlux
!      real, dimension(1:6)     :: OldJminusFlux
!      real, dimension(1:6)     :: OldKminusFlux
!      real, dimension(1:6)     :: NewIminusFlux
!      real, dimension(1:6)     :: NewJminusFlux
!      real, dimension(1:6)     :: NewKminusFlux
!      real, dimension(1:6)     :: DelIminusFlux
!      real, dimension(1:6)     :: DelJminusFlux
!      real, dimension(1:6)     :: DelKminusFlux
!      real, dimension(1:6)     :: LambdaTimesArea
!      real, dimension(1:6)     :: Q0 ! state at cell
!      real, dimension(1:6)     :: Q1 ! state at neighbours 
!      real, dimension(1:6)     :: Q2
!      real, dimension(1:6)     :: Q3
!      real, dimension(1:6)     :: Q4
!      real, dimension(1:6)     :: Q5
!      real, dimension(1:6)     :: Q6
!      real, dimension(1:6)     :: DQ0! change in state
!      real, dimension(1:6)     :: DQ1
!      real, dimension(1:6)     :: DQ2
!      real, dimension(1:6)     :: DQ3
!      real, dimension(1:6)     :: DQ4
!      real, dimension(1:6)     :: DQ5
!      real, dimension(1:6)     :: DQ6
!      real, dimension(1:7)     :: Flist1
!      real, dimension(1:7)     :: Flist2
!      real, dimension(1:7)     :: Flist3
!      real, dimension(1:7)     :: Flist4
!      real, dimension(1:7)     :: Flist5
!      real, dimension(1:7)     :: Flist6
!      real, dimension(1:3)     :: C0
!      real, dimension(1:3)     :: C1
!      real, dimension(1:3)     :: C2
!      real, dimension(1:3)     :: C3
!      real, dimension(1:3)     :: C4
!      real, dimension(1:3)     :: C5
!      real, dimension(1:3)     :: C6
!      real :: fv1
!      real :: fv2
!      real :: fw
!      real :: g
!      real :: Scap
!      real :: r
!      real :: S_v
!      real :: D_v
!      real :: P_v
!      real :: lamda
!      real :: dist_i
!      real :: dist_i_2
!      real :: Ji
!      real :: Ji_2
!      real :: Ji_3
!      real :: S
!      real :: Omega
!      real :: k2
!      real :: inv_k2_d2
!      real :: Shat
!      real :: inv_Shat
!      real :: nu
!      real :: nu_t
!      real :: glim
!      real :: g_6
!      real :: gamma_BC
!      real :: dfv1
!      real :: dfv2
!      real :: dfw
!      real :: dShat
!      real :: dr
!      real :: dg
!
!
!
!        !intialize delQ
!        delQstar = 0.0
!
!        !forward sweep
!        do k=1,kmx-1
!          do j=1,jmx-1
!            do i=1,imx-1
!              C0  = CellCenter(i  ,j  ,k  ,:)
!              C1  = CellCenter(i-1,j  ,k  ,:)
!              C2  = CellCenter(i  ,j-1,k  ,:)
!              C3  = CellCenter(i  ,j  ,k-1,:)
!              C4  = CellCenter(i+1,j  ,k  ,:)
!              C5  = CellCenter(i  ,j+1,k  ,:)
!              C6  = CellCenter(i  ,j  ,k+1,:)
!
!              Q0  = qp(i  , j  , k  , 1:6)
!              Q1  = qp(i-1, j  , k  , 1:6)
!              Q2  = qp(i  , j-1, k  , 1:6)
!              Q3  = qp(i  , j  , k-1, 1:6)
!              Q4  = qp(i+1, j  , k  , 1:6)
!              Q5  = qp(i  , j+1, k  , 1:6)
!              Q6  = qp(i  , j  , k+1, 1:6)
!
!              DQ0 = 0.0
!              DQ1 = delQstar(i-1, j  , k  , 1:6)
!              DQ2 = delQstar(i  , j-1, k  , 1:6)
!              DQ3 = delQstar(i  , j  , k-1, 1:6)
!
!              Flist1(1) =   xA(i,j,k)
!              Flist1(2) = -xnx(i,j,k)
!              Flist1(3) = -xny(i,j,k)
!              Flist1(4) = -xnz(i,j,k)
!              Flist1(5) = 0.5*(volume(i-1, j  , k  ) + volume(i,j,k))
!              Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
!              Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))
!
!              Flist2(1) =   yA(i,j,k)
!              Flist2(2) = -ynx(i,j,k)
!              Flist2(3) = -yny(i,j,k)
!              Flist2(4) = -ynz(i,j,k)
!              Flist2(5) = 0.5*(volume(i  , j-1, k  ) + volume(i,j,k))
!              Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
!              Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))
!
!              Flist3(1) =   zA(i,j,k)
!              Flist3(2) = -znx(i,j,k)
!              Flist3(3) = -zny(i,j,k)
!              Flist3(4) = -znz(i,j,k)
!              Flist3(5) = 0.5*(volume(i  , j  , k-1) + volume(i,j,k))
!              Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
!              Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))
!
!              Flist4(1) =   xA(i+1,j,k)
!              Flist4(2) = +xnx(i+1,j,k)
!              Flist4(3) = +xny(i+1,j,k)
!              Flist4(4) = +xnz(i+1,j,k)
!              Flist4(5) = 0.5*(volume(i+1, j  , k  ) + volume(i,j,k))
!              Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
!              Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))
!
!              Flist5(1) =   yA(i,j+1,k)
!              Flist5(2) = +ynx(i,j+1,k)
!              Flist5(3) = +yny(i,j+1,k)
!              Flist5(4) = +ynz(i,j+1,k)
!              Flist5(5) = 0.5*(volume(i  , j+1, k  ) + volume(i,j,k))
!              Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
!              Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))
!
!              Flist6(1) =   zA(i,j,k+1)
!              Flist6(2) = +znx(i,j,k+1)
!              Flist6(3) = +zny(i,j,k+1)
!              Flist6(4) = +znz(i,j,k+1)
!              Flist6(5) = 0.5*(volume(i  , j  , k+1) + volume(i,j,k))
!              Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
!              Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))
!
!              NewIminusFlux     = SAFlux(Q1, Q0, DQ1, Flist1)
!              NewJminusFlux     = SAFlux(Q2, Q0, DQ2, Flist2)
!              NewKminusFlux     = SAFlux(Q3, Q0, DQ3, Flist3)
!              OldIminusFlux     = SAFlux(Q1, Q0, DQ0, Flist1)
!              OldJminusFlux     = SAFlux(Q2, Q0, DQ0, Flist2)
!              OldKminusFlux     = SAFlux(Q3, Q0, DQ0, Flist3)
!
!              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0)
!              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0)
!              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0)
!              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0)
!              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0)
!              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0)
!
!
!              ! multiply above flux with area to get correct values
!              DelIminusFlux =  NewIminusFlux - OldIminusFlux
!              DelJminusFlux =  NewJminusFlux - OldJminusFlux
!              DelKminusFlux =  NewKminusFlux - OldKminusFlux
!
!
!              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
!              !storing D in Iflux array for backward sweep
!              !F_p(i,j,k,1) = D
!              ! -- source term derivatives -- !
!              Omega = sqrt( ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
!                           + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
!                           + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
!                            )&
!                       )
!              dist_i = dist(i,j,k)
!              dist_i_2 = dist_i*dist_i
!              k2 = kappa_sa*kappa_sa
!              nu   = mu(i,j,k)/density(i,j,k)
!              Ji   = Q0(6)/nu
!              Ji_2 = Ji*Ji
!              Ji_3 = Ji_2*ji
!
!
!              ! ___ functions ___
!              fv1  = (Ji_3)/((Ji_3) + (cv1_3))
!              fv2  = 1.0 - Ji/(1.0 + (Ji*fv1))
!
!              ! ___ Shear stress for production ___
!              S = Omega
!              inv_k2_d2 = 1.0/(k2*dist_i_2)
!              Shat      = S + Q0(6)*fv2*inv_k2_d2
!              Shat      = max(Shat, 1.0e-10)
!              inv_Shat  = 1.0/Shat
!              dfv1 = 3.0*Ji_2*cv1_3/(nu*(Ji_3+cv1_3)**2)
!              dfv2 = -((1.0/nu) - Ji_2*dfv1)/((1.0+Ji*fv1)**2)
!              dShat = (fv2+Q0(6)*dfv2)*inv_k2_d2
!
!              D = D - cb1*(Q0(6)*dShat+Shat)*Volume(i,j,k)
!
!              ! ___ Destruction term___ !
!              r    = min(Q0(6)*inv_Shat*inv_k2_d2, 10.0)
!              g    = r + cw2*((r**6) - r)
!              g_6  = g**6
!              glim = ((1.0+cw3_6)/(g_6+cw3_6))**(1.0/6.0)
!              fw   = g*glim
!              dr = (Shat-Q0(6)*dShat)*inv_Shat*inv_Shat*inv_k2_d2
!              dg = dr*(1.0+cw2*(6.0*(r**5)-1.0))
!              dfw= dg*glim*(1.0-g_6/(g_6+cw3_6))
!
!              D = D+cw1*(dfw*Q0(6) + 2*fw)*Q0(6)/dist_i_2*volume(i,j,k)
!              ! --  end of source term -- !
!
!              deltaU(1:6) = -residue(i,j,k,1:6) &
!                - 0.5*((DelIminusFlux - LambdaTimesArea(1)*delQstar(i-1,j,k,1:6)) &
!                     + (DelJminusFlux - LambdaTimesArea(2)*delQstar(i,j-1,k,1:6)) &
!                     + (DelKminusFlux - LambdaTimesArea(3)*delQstar(i,j,k-1,1:6)) )
!
!              delQstar(i,j,k,1:6) = deltaU(1:6)/D
!            end do
!          end do
!        end do
!
!        !call apply_interface(delQstar, 1)
!
!        delQ=0.0
!        !backward sweep
!            do i=imx-1,1,-1
!          do j=jmx-1,1,-1
!        do k=kmx-1,1,-1
!              C0  = CellCenter(i  ,j  ,k  ,:)
!              C1  = CellCenter(i-1,j  ,k  ,:)
!              C2  = CellCenter(i  ,j-1,k  ,:)
!              C3  = CellCenter(i  ,j  ,k-1,:)
!              C4  = CellCenter(i+1,j  ,k  ,:)
!              C5  = CellCenter(i  ,j+1,k  ,:)
!              C6  = CellCenter(i  ,j  ,k+1,:)
!
!              Q0  = qp(i  , j  , k  , 1:6)
!              Q1  = qp(i-1, j  , k  , 1:6)
!              Q2  = qp(i  , j-1, k  , 1:6)
!              Q3  = qp(i  , j  , k-1, 1:6)
!              Q4  = qp(i+1, j  , k  , 1:6)
!              Q5  = qp(i  , j+1, k  , 1:6)
!              Q6  = qp(i  , j  , k+1, 1:6)
!
!              DQ0 = 0.0
!              DQ4 = delQ(i+1, j  , k  , 1:6)
!              DQ5 = delQ(i  , j+1, k  , 1:6)
!              DQ6 = delQ(i  , j  , k+1, 1:6)
!
!              Flist1(1) =   xA(i,j,k)
!              Flist1(2) = -xnx(i,j,k)
!              Flist1(3) = -xny(i,j,k)
!              Flist1(4) = -xnz(i,j,k)
!              Flist1(5) = 0.5*(volume(i-1, j  , k  ) + volume(i,j,k))
!              Flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
!              Flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))
!
!              Flist2(1) =   yA(i,j,k)
!              Flist2(2) = -ynx(i,j,k)
!              Flist2(3) = -yny(i,j,k)
!              Flist2(4) = -ynz(i,j,k)
!              Flist2(5) = 0.5*(volume(i  , j-1, k  ) + volume(i,j,k))
!              Flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
!              Flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))
!
!              Flist3(1) =   zA(i,j,k)
!              Flist3(2) = -znx(i,j,k)
!              Flist3(3) = -zny(i,j,k)
!              Flist3(4) = -znz(i,j,k)
!              Flist3(5) = 0.5*(volume(i  , j  , k-1) + volume(i,j,k))
!              Flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
!              Flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))
!
!              Flist4(1) =   xA(i+1,j,k)
!              Flist4(2) = +xnx(i+1,j,k)
!              Flist4(3) = +xny(i+1,j,k)
!              Flist4(4) = +xnz(i+1,j,k)
!              Flist4(5) = 0.5*(volume(i+1, j  , k  ) + volume(i,j,k))
!              Flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
!              Flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))
!
!              Flist5(1) =   yA(i,j+1,k)
!              Flist5(2) = +ynx(i,j+1,k)
!              Flist5(3) = +yny(i,j+1,k)
!              Flist5(4) = +ynz(i,j+1,k)
!              Flist5(5) = 0.5*(volume(i  , j+1, k  ) + volume(i,j,k))
!              Flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
!              Flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))
!
!              Flist6(1) =   zA(i,j,k+1)
!              Flist6(2) = +znx(i,j,k+1)
!              Flist6(3) = +zny(i,j,k+1)
!              Flist6(4) = +znz(i,j,k+1)
!              Flist6(5) = 0.5*(volume(i  , j  , k+1) + volume(i,j,k))
!              Flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
!              Flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))
!
!              NewIminusFlux     = SAFlux(Q4, Q0, DQ4, Flist4)
!              NewJminusFlux     = SAFlux(Q5, Q0, DQ5, Flist5)
!              NewKminusFlux     = SAFlux(Q6, Q0, DQ6, Flist6)
!              OldIminusFlux     = SAFlux(Q4, Q0, DQ0, Flist4)
!              OldJminusFlux     = SAFlux(Q5, Q0, DQ0, Flist5)
!              OldKminusFlux     = SAFlux(Q6, Q0, DQ0, Flist6)
!
!              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0)
!              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0)
!              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0)
!              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0)
!              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0)
!              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0)
!
!
!              ! multiply above flux with area to get correct values
!              DelIminusFlux =  NewIminusFlux - OldIminusFlux
!              DelJminusFlux =  NewJminusFlux - OldJminusFlux
!              DelKminusFlux =  NewKminusFlux - OldKminusFlux
!
!              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
!
!              ! -- source term derivatives -- !
!              Omega = sqrt( ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
!                           + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
!                           + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
!                            )&
!                       )
!              dist_i = dist(i,j,k)
!              dist_i_2 = dist_i*dist_i
!              k2 = kappa_sa*kappa_sa
!              nu   = mu(i,j,k)/density(i,j,k)
!              Ji   = Q0(6)/nu
!              Ji_2 = Ji*Ji
!              Ji_3 = Ji_2*ji
!
!
!              ! ___ functions ___
!              fv1  = (Ji_3)/((Ji_3) + (cv1_3))
!              fv2  = 1.0 - Ji/(1.0 + (Ji*fv1))
!              ! ___ Shear stress for production ___
!              S = Omega
!              inv_k2_d2 = 1.0/(k2*dist_i_2)
!              Shat      = S + Q0(6)*fv2*inv_k2_d2
!              Shat      = max(Shat, 1.0e-10)
!              inv_Shat  = 1.0/Shat
!              dfv1 = 3.0*Ji_2*cv1_3/(nu*(Ji_3+cv1_3)**2)
!              dfv2 = -((1.0/nu) - Ji_2*dfv1)/((1.0+Ji*fv1)**2)
!              dShat = (fv2+Q0(6)*dfv2)*inv_k2_d2
!
!              D = D - cb1*(Q0(6)*dShat+Shat)*Volume(i,j,k)
!
!              ! ___ Destruction term___ !
!              r    = min(Q0(6)*inv_Shat*inv_k2_d2, 10.0)
!              g    = r + cw2*((r**6) - r)
!              g_6  = g**6
!              glim = ((1.0+cw3_6)/(g_6+cw3_6))**(1.0/6.0)
!              fw   = g*glim
!              dr = (Shat-Q0(6)*dShat)*inv_Shat*inv_Shat*inv_k2_d2
!              dg = dr*(1.0+cw2*(6.0*(r**5)-1.0))
!              dfw= dg*glim*(1.0-g_6/(g_6+cw3_6))
!
!              D = D+cw1*(dfw*Q0(6) + 2*fw)*Q0(6)/dist_i_2*volume(i,j,k)
!              ! --  end of source term -- !
!
!              delQ(i,j,k,1:6) = delQstar(i,j,k,1:6) &
!                - 0.5*((DelIminusFlux - LambdaTimesArea(4)*delQ(i+1,j,k,1:6)) &
!                     + (DelJminusFlux - LambdaTimesArea(5)*delQ(i,j+1,k,1:6)) &
!                     + (DelKminusFlux - LambdaTimesArea(6)*delQ(i,j,k+1,1:6)) )/D
!
!            end do
!          end do
!        end do
!        
!        do k=1,kmx-1
!          do j = 1,jmx-1
!            do i = 1,imx-1
!              conservativeQ(1) = qp(i,j,k,1)
!              conservativeQ(2) = qp(i,j,k,1) * qp(i,j,k,2)
!              conservativeQ(3) = qp(i,j,k,1) * qp(i,j,k,3)
!              conservativeQ(4) = qp(i,j,k,1) * qp(i,j,k,4)
!              conservativeQ(5) = (qp(i,j,k,5) / (gm-1.0)) + ( 0.5 * qp(i,j,k,1) * sum( qp(i,j,k,2:4)**2) )
!              conservativeQ(6) = qp(i,j,k,1) * qp(i,j,k,6)
!              
!              ! add new change into conservative solution
!              conservativeQ(1:6) = conservativeQ(1:6) + delQ(i,j,k,1:6)
!
!              ! convert back conservative to primitive
!              qp(i,j,k,1) = conservativeQ(1)
!              qp(i,j,k,2) = conservativeQ(2) / conservativeQ(1)
!              qp(i,j,k,3) = conservativeQ(3) / conservativeQ(1)
!              qp(i,j,k,4) = conservativeQ(4) / conservativeQ(1)
!              qp(i,j,k,5) = (gm-1.0) * ( conservativeQ(5) - (0.5 * sum(conservativeQ(2:4)**2) / conservativeQ(1)) )
!              qp(i,j,k,6) = conservativeQ(6) / conservativeQ(1)
!              qp(i,j,k,6) = max(qp(i,j,k,6), 1.e-8)
!            end do
!          end do
!        end do
!
!    end subroutine update_SA_variables


!    function Flux(ql, qr, du, inputs)
!      !--------------------------------------
!      ! calculate the total flux through face
!      !---------------------------------------
!      implicit none
!      real, dimension(1:n_var), intent(in) :: ql !left state
!      real, dimension(1:n_var), intent(in) :: qr !right state
!      !conservative form of updated neighbour
!      real, dimension(1:n_var), intent(in) :: du
!      real, dimension(1:7)    , intent(in) :: inputs
!      real, dimension(1:n_var)             :: Flux
!      real, dimension(1:n_var)             :: SAFlux
!      real, dimension(1:n_var)             :: U ! conservative variables
!      real, dimension(1:n_var)             :: W ! new primitive variables
!      real, dimension(1:n_var)             :: P ! primitive variables of right cell
!
!      !for extraction of the inputs
!      real :: area
!      real :: nx
!      real :: ny
!      real :: nz
!      real :: volume
!      real :: mmu
!      real :: tmu
!
!
!      real    :: dudx
!      real    :: dudy
!      real    :: dudz
!      real    :: dvdx
!      real    :: dvdy
!      real    :: dvdz
!      real    :: dwdx
!      real    :: dwdy
!      real    :: dwdz
!      real    :: dTdx
!      real    :: dTdy
!      real    :: dTdz
!      real    :: dtvdx
!      real    :: dtvdy
!      real    :: dtvdz
!      real    :: T1, T2
!      real    :: uface
!      real    :: vface
!      real    :: wface
!      real    :: trace
!      real    :: Tauxx
!      real    :: Tauyy
!      real    :: Tauzz
!      real    :: Tauxy
!      real    :: Tauxz
!      real    :: Tauyz
!      real    :: Qx
!      real    :: Qy
!      real    :: Qz
!      real    :: HalfRhoUsquare
!      real    :: RhoHt
!      real    :: K_heat
!      real    :: FaceNormalVelocity
!      real    :: mu
!      real    :: muCap
!
!      area   = inputs(1)
!      nx     = inputs(2)
!      ny     = inputs(3)
!      nz     = inputs(4)
!      volume = inputs(5)
!      mmu    = inputs(6)
!      tmu    = inputs(7)
!
!
!      !save the old stat in P
!      P = qr
!
!      ! find conservative variable
!      U(1)   =   ql(1)
!      U(2:n_var)   =   ql(1) * ql(2:n_var)
!     ! U(3)   =   ql(1) * ql(3)
!     ! U(4)   =   ql(1) * ql(4)
!      U(5)   = ( ql(5) / (gm-1.0) ) + ( 0.5 * ql(1) * sum(ql(2:4)**2) )
!     ! U(6)   =   ql(1) * ql(6)
!
!      U(1:n_var) = U(1:n_var) + du(1:n_var)
!      
!
!      W(1)   =   U(1)
!      W(2)   =   U(2) / U(1)
!      W(3)   =   U(3) / U(1)
!      W(4)   =   U(4) / U(1)
!      W(2:n_var) = U(2:n_var) / U(1)
!      W(5)   = (gm-1.0) * ( U(5) - ( 0.5 * SUM(U(2:4)**2) / U(1) ) )
!      !W(6)   =   U(6) / U(1)
!      W(5:) = max(W(5:), 1e-8)
!
!      FaceNormalVelocity = (W(2) * nx) + (W(3) * ny) + (W(4) * nz)
!      uface = 0.5 * ( W(2) + P(2) )
!      vface = 0.5 * ( W(3) + P(3) )
!      wface = 0.5 * ( W(4) + P(4) )
!
!
!      Flux(1) =   W(1) * FaceNormalVelocity
!      Flux(2) = ( W(2) * Flux(1) ) + ( W(5) * nx )
!      Flux(3) = ( W(3) * Flux(1) ) + ( W(5) * ny )
!      Flux(4) = ( W(4) * Flux(1) ) + ( W(5) * nz )
!
!      HalfRhoUsquare = 0.5 * W(1) * ( W(2)*W(2) + W(3)*W(3) + W(4)*W(4) )
!      RhoHt          = ( (gm/(gm-1.0)) * W(5) ) + HalfRhoUsquare
!      Flux(5)        = RhoHt * FaceNormalVelocity
!      Flux(6) = ( W(6) * Flux(1) )   
!
!
!      ! viscous terms
!!      muCap = 0.25*(P(1)+W(1))*(P(6) + W(6))
!      mu = mmu + tmu
!      T1     =    W(5) / ( W(1) * R_gas )
!      T2     =    P(5) / ( P(1) * R_gas )
!      dTdx   =  ( T2   - T1   ) * nx * Area / Volume
!      dTdy   =  ( T2   - T1   ) * ny * Area / Volume
!      dTdz   =  ( T2   - T1   ) * nz * Area / Volume
!      dudx   =  ( P(2) - W(2) ) * nx * Area / Volume
!      dudy   =  ( P(2) - W(2) ) * ny * Area / Volume
!      dudz   =  ( P(2) - W(2) ) * nz * Area / Volume
!      dvdx   =  ( P(3) - W(3) ) * nx * Area / Volume
!      dvdy   =  ( P(3) - W(3) ) * ny * Area / Volume
!      dvdz   =  ( P(3) - W(3) ) * nz * Area / Volume
!      dwdx   =  ( P(4) - W(4) ) * nx * Area / Volume
!      dwdy   =  ( P(4) - W(4) ) * ny * Area / Volume
!      dwdz   =  ( P(4) - W(4) ) * nz * Area / Volume
!      dtvdx  =  ( P(6) - W(6) ) * nx * Area / Volume
!      dtvdy  =  ( P(6) - W(6) ) * ny * Area / Volume
!      dtvdz  =  ( P(6) - W(6) ) * nz * Area / Volume
!
!      trace = dudx + dvdy + dwdz
!      Tauxx =  2. * mu * (dudx - trace/3.0)
!      Tauyy =  2. * mu * (dvdy - trace/3.0)
!      Tauzz =  2. * mu * (dwdz - trace/3.0)
!      Tauxy = mu * (dvdx + dudy)
!      Tauxz = mu * (dwdx + dudz)
!      Tauyz = mu * (dwdy + dvdz)
!
!      K_heat = ( mmu / Pr  + tmu/tpr) * gm * R_gas / ( gm - 1.0 )
!      Qx = K_heat*dTdx
!      Qy = K_heat*dTdy
!      Qz = K_heat*dTdz
!      tmu = 0.5*(W(6) + P(6))
!
!      Flux(2) = Flux(2) - ( Tauxx * nx + Tauxy * ny + Tauxz * nz )
!      Flux(3) = Flux(3) - ( Tauxy * nx + Tauyy * ny + Tauyz * nz )
!      Flux(4) = Flux(4) - ( Tauxz * nx + Tauyz * ny + Tauzz * nz )
!      Flux(5) = Flux(5) - ( Tauxx * uface + Tauxy * vface + Tauxz * wface + Qx ) * nx
!      Flux(5) = Flux(5) - ( Tauxy * uface + Tauyy * vface + Tauyz * wface + Qy ) * ny
!      Flux(5) = Flux(5) - ( Tauxz * uface + Tauyz * vface + Tauzz * wface + Qz ) * nz
!!      Flux(6) = Flux(6) - (mmu + muCap)*(dtvdx*nx + dtvdy*ny + dtvdz*nz)/sigma_sa
!
!      Flux    = Flux * Area
!      SAFlux = Flux
!
!    end function Flux 


    subroutine apply_interface(state, layers)
      implicit none
      integer, intent(in):: layers
      real, dimension(1-layers:imx-1+layers,1-layers:jmx-1+layers,1-layers:kmx-1+layers,1:n_var), intent(inout):: state
      integer:: i,j,k,n,l
      integer:: status(MPI_STATUS_SIZE)
      integer:: ierr
      integer:: tag=1
      integer:: count=0


      !--- IMIN ---!
      DebugCall("apply_interface")
      if(imin_id>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do j=1,jmx-1
                count=count+1
                imin_send_buf(count) = state(l,j,k,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(imin_send_buf,ibuf_size, MPI_DOUBLE_PRECISION, imin_id,tag,&
                          imin_recv_buf,ibuf_size, MPI_DOUBLE_PRECISION, imin_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        if(dir_switch(1)==0)then
          count=0
          do n=1,n_var
            do l=1,layers
              do k=Pklo(1),Pkhi(1),PkDir(1)
                do j=Pjlo(1),Pjhi(1),PjDir(1)
                  count=count+1
                 state(1-l,j,k,n) = imin_recv_buf(count)
                end do
              end do
            end do
          end do
        else
          count=0
          do n=1,n_var
            do l=1,layers
              do j=Pjlo(1),Pjhi(1),PjDir(1)
                do k=Pklo(1),Pkhi(1),PkDir(1)
                  count=count+1
                  state(1-l,j,k,n) = imin_recv_buf(count)
                end do
              end do
            end do
          end do
        end if
      end if

      !--- IMAX ---!
      if(imax_id>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do j=1,jmx-1
                count=count+1
                imax_send_buf(count) = state(imx-l,j,k,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(imax_send_buf,ibuf_size, MPI_DOUBLE_PRECISION, imax_id,tag,&
                          imax_recv_buf,ibuf_size, MPI_DOUBLE_PRECISION, imax_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        if(dir_switch(2)==0)then
          count=0
          do n=1,n_var
            do l=1,layers
              do k=Pklo(2),Pkhi(2),PkDir(2)
                do j=Pjlo(2),Pjhi(2),PjDir(2)
                  count=count+1
                   state(imx+l-1,j,k,n) = imax_recv_buf(count)
                end do
              end do
            end do
          end do
        else
          count=0
          do n=1,n_var
            do l=1,layers
              do j=Pjlo(2),Pjhi(2),Pjdir(2)
                do k=Pklo(2),Pkhi(2),PkDir(2)
                  count=count+1
                   state(imx+l-1,j,k,n) = imax_recv_buf(count)
                end do
              end do
            end do
          end do
        end if
      end if



      !--- JMIN ---!
      if(jmin_id>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do i=1,imx-1
                count=count+1
                jmin_send_buf(count) = state(i,l,k,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(jmin_send_buf,jbuf_size, MPI_DOUBLE_PRECISION, jmin_id,tag,&
                          jmin_recv_buf,jbuf_size, MPI_DOUBLE_PRECISION, jmin_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        if(dir_switch(3)==0)then
          count=0
          do n=1,n_var
            do l=1,layers
              do k=Pklo(3),Pkhi(3),PkDir(3)
                do i=Pilo(3),Pihi(3),PiDir(3)
                  count=count+1
                  state(i,1-l,k,n) = jmin_recv_buf(count)
                end do
              end do
            end do
          end do
        else
          count=0
          do n=1,n_var
            do l=1,layers
              do i=Pilo(3),Pihi(3),PiDir(3)
                do k=Pklo(3),Pkhi(3),PkDir(3)
                  count=count+1
                  state(i,1-l,k,n) = jmin_recv_buf(count)
                end do
              end do
            end do
          end do
        end if
      end if

      !--- JMAX ---!
      if(jmax_id>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do i=1,imx-1
                count=count+1
                jmax_send_buf(count) = state(i,jmx-l,k,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(jmax_send_buf,jbuf_size, MPI_DOUBLE_PRECISION, jmax_id,tag,&
                          jmax_recv_buf,jbuf_size, MPI_DOUBLE_PRECISION, jmax_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        if(dir_switch(4)==0)then
          count=0
          do n=1,n_var
            do l=1,layers
              do k=Pklo(4),Pkhi(4),PkDir(4)
                do i=Pilo(4),Pihi(4),PiDir(4)
                  count=count+1
                  state(i,jmx+l-1,k,n) = jmax_recv_buf(count)
                end do
              end do
            end do
          end do
        else
          count=0
          do n=1,n_var
            do l=1,layers
              do i=Pilo(4),Pihi(4),PiDir(4)
                do k=Pklo(4),Pkhi(4),PkDir(4)
                  count=count+1
                  state(i,jmx+l-1,k,n) = jmax_recv_buf(count)
                end do
              end do
            end do
          end do
        end if
      end if

      !--- KMIN ---!
      if(kmin_id>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do j=1,jmx-1
              do i=1,imx-1
                count=count+1
                kmin_send_buf(count) = state(i,j,l,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(kmin_send_buf,kbuf_size, MPI_DOUBLE_PRECISION, kmin_id,tag,&
                          kmin_recv_buf,kbuf_size, MPI_DOUBLE_PRECISION, kmin_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        if(dir_switch(5)==0)then
          count=0
          do n=1,n_var
            do l=1,layers
              do j=Pjlo(5),Pjhi(5),PjDir(5)
                do i=Pilo(5),Pihi(5),PiDir(5)
                  count=count+1
                  state(i,j,1-l,n) = kmin_recv_buf(count)
                end do
              end do
            end do
          end do
        else
          count=0
          do n=1,n_var
            do l=1,layers
              do i=Pilo(5),Pihi(5),PiDir(5)
                do j=Pjlo(5),Pjhi(5),PjDir(5)
                  count=count+1
                  state(i,j,1-l,n) = kmin_recv_buf(count)
                end do
              end do
            end do
          end do
        end if
      end if

      !--- KMAX ---!
      if(kmax_id>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do j=1,jmx-1
              do i=1,imx-1
                count=count+1
                kmax_send_buf(count) = state(i,j,kmx-l,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(kmax_send_buf,kbuf_size, MPI_DOUBLE_PRECISION, kmax_id,tag,&
                          kmax_recv_buf,kbuf_size, MPI_DOUBLE_PRECISION, kmax_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        if(dir_switch(6)==0)then
          count=0
          do n=1,n_var
            do l=1,layers
              do j=Pjlo(6),Pjhi(6),PjDir(6)
                do i=Pilo(6),Pihi(6),PiDir(6)
                  count=count+1
                  state(i,j,kmx+l-1,n) = kmax_recv_buf(count)
                end do
              end do
            end do
          end do
        else
          count=0
          do n=1,n_var
            do l=1,layers
              do i=Pilo(6),Pihi(6),PiDir(6)
                do j=Pjlo(6),Pjhi(6),PjDir(6)
                  count=count+1
                  state(i,j,kmx+l-1,n) = kmax_recv_buf(count)
                end do
              end do
            end do
          end do
        end if
      end if
      call apply_periodic_bc(state, 1)
    end subroutine apply_interface

    subroutine apply_periodic_bc(state, layers)
      implicit none
      integer, intent(in) :: layers
      real, dimension(1-layers:imx-1+layers,1-layers:jmx-1+layers,1-layers:kmx-1+layers,1:n_var), intent(inout) :: state
      integer:: i,j,k,n,l
      integer:: status(MPI_STATUS_SIZE)
      integer:: ierr
      integer:: tag=1
      integer:: count=0

      call dmsg(1, 'interface', 'apply_periodic_boundary_condition')
      if(PbcId(1)>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do j=1,jmx-1
                count=count+1
                imin_send_buf(count) = state(l,j,k,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(imin_send_buf,ibuf_size, MPI_DOUBLE_PRECISION, PbcId(1),tag,&
                          imin_recv_buf,ibuf_size, MPI_DOUBLE_PRECISION, PbcId(1),tag,&
                          MPI_COMM_WORLD,status,ierr)
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do j=1,jmx-1
                count=count+1
                state(1-l,j,k,n) = imin_recv_buf(count)
              end do
            end do
          end do
        end do
      end if

      if(PbcId(2)>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do j=1,jmx-1
                count=count+1
                imax_send_buf(count) = state(imx-l,j,k,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(imax_send_buf,ibuf_size, MPI_DOUBLE_PRECISION, PbcId(2),tag,&
                          imax_recv_buf,ibuf_size, MPI_DOUBLE_PRECISION, PbcId(2),tag,&
                          MPI_COMM_WORLD,status,ierr)
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do j=1,jmx-1
                count=count+1
                 state(imx+l-1,j,k,n) = imax_recv_buf(count)
              end do
            end do
          end do
        end do
      end if
      !--- JMIN ---!
      if(PbcId(3)>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do i=1,imx-1
                count=count+1
                jmin_send_buf(count) = state(i,l,k,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(jmin_send_buf,jbuf_size, MPI_DOUBLE_PRECISION, PbcId(3),tag,&
                          jmin_recv_buf,jbuf_size, MPI_DOUBLE_PRECISION, PbcId(3),tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do i=1,imx-1
                count=count+1
                state(i,1-l,k,n) = jmin_recv_buf(count)
              end do
            end do
          end do
        end do
      end if

      !--- JMAX ---!
      if(PbcId(4)>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do i=1,imx-1
                count=count+1
                jmax_send_buf(count) = state(i,jmx-l,k,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(jmax_send_buf,jbuf_size, MPI_DOUBLE_PRECISION, PbcId(4),tag,&
                          jmax_recv_buf,jbuf_size, MPI_DOUBLE_PRECISION, PbcId(4),tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do i=1,imx-1
                count=count+1
                state(i,jmx+l-1,k,n) = jmax_recv_buf(count)
              end do
            end do
          end do
        end do
      end if

      !--- KMIN ---!
      if(PbcId(5)>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do j=1,jmx-1
              do i=1,imx-1
                count=count+1
                kmin_send_buf(count) = state(i,j,l,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(kmin_send_buf,kbuf_size, MPI_DOUBLE_PRECISION, PbcId(5),tag,&
                          kmin_recv_buf,kbuf_size, MPI_DOUBLE_PRECISION, PbcId(5),tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        count=0
        do n=1,n_var
          do l=1,layers
            do j=1,jmx-1
              do i=1,imx-1
                count=count+1
                state(i,j,1-l,n) = kmin_recv_buf(count)
              end do
            end do
          end do
        end do
      end if

      !--- KMAX ---!
      if(PbcId(6)>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do j=1,jmx-1
              do i=1,imx-1
                count=count+1
                kmax_send_buf(count) = state(i,j,kmx-l,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(kmax_send_buf,kbuf_size, MPI_DOUBLE_PRECISION, PbcId(6),tag,&
                          kmax_recv_buf,kbuf_size, MPI_DOUBLE_PRECISION, PbcId(6),tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        count=0
        do n=1,n_var
          do l=1,layers
            do j=1,jmx-1
              do i=1,imx-1
                count=count+1
                state(i,j,kmx+l-1,n) = kmax_recv_buf(count)
              end do
            end do
          end do
        end do
      end if


    end subroutine apply_periodic_bc
end module HLU_SGS
