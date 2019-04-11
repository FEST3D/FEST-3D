module solvergmres

  use global_vars, only : R_gas
  use global_vars, only : Pr
  use global_vars, only : tPr
  use global_vars, only : gm
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

  use global_vars, only : time_stepping_method
  use global_vars, only : time_step_accuracy
  use global_vars, only : global_time_step
  use global_vars, only : delta_t
  use global_vars, only : turbulence
  use global_vars, only : process_id
  use global_vars, only : total_process

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

  use geometry   , only: CellCenter

  use utils, only: alloc
  use utils, only:  dealloc 
  use utils, only:  dmsg
  use utils, only:  DEBUG_LEVEL

  use string

  !subroutine for residual calculation
  use face_interpolant,               only: interpolant
  use global_vars,                    only : mu_ref
  use interface1,                      only: apply_interface
  use bc_primitive,                   only: populate_ghost_primitive
  use face_interpolant,               only: compute_face_interpolant
  use boundary_state_reconstruction,  only: reconstruct_boundary_state
  use scheme,                         only: compute_fluxes
  use summon_grad_evaluation,         only: evaluate_all_gradients
  use viscosity                      ,only: calculate_viscosity
  use viscous,                        only: compute_viscous_fluxes
  use scheme,                         only: compute_residue
  use source,                         only: add_source_term_residue

  use time,                           only : compute_time_step
#include "error.inc"
#include "mpi.inc"

  real, dimension(:,:,:,:), allocatable :: Avj
  real, dimension(:,:,:,:), allocatable :: MinvAvj
  real, dimension(:), allocatable :: wj, vi, GlobalVj
  real, dimension(:,:,:,:), allocatable :: vj

  real, dimension(:,:,:,:), allocatable :: r0
  real, dimension(:,:,:,:), allocatable :: Minvr0 
  real, dimension(:,:,:,:), allocatable :: Ax0, Ax
  real, dimension(:,:,:,:), allocatable :: r, v1

  real, dimension(:,:), allocatable :: kspvectors


  real, dimension(:,:,:,:), allocatable :: old_state
  real, dimension(:,:,:,:), allocatable :: conser
  real, dimension(:,:,:,:), allocatable :: old_residue


  real, dimension(:,:,:,:), allocatable :: delQ
  real, dimension(:,:,:,:), allocatable :: delQstar
  real, dimension(:,:,:), allocatable, target :: dummy
  real, dimension(:,:,:), pointer :: tmu
  real, dimension(:,:,:), pointer :: mmu


  integer :: IDimMax
  integer :: JDimMax
  integer :: KDimMax
  integer, dimension(:), allocatable :: IDim
  integer, dimension(:), allocatable :: JDim
  integer, dimension(:), allocatable :: KDim
  integer :: SendSize
  integer :: GlobalSize
  integer, dimension(:), allocatable :: RecvSize
  integer, dimension(:), allocatable :: Displacement
  integer, dimension(:), allocatable :: StartIndex

  real    :: NormSum
  real, dimension(:), allocatable :: SendBuffer
  real, dimension(:), allocatable :: RecvBuffer
  real, dimension(:), allocatable :: Globalx
  real, dimension(:), allocatable :: Globalv1

contains

  subroutine setup_solvergmres(m)
    implicit none
    integer, intent(in) :: m
    integer :: ierr
    integer :: i

    call alloc(IDim, 1, total_process)
    call alloc(JDim, 1, total_process)
    call alloc(KDim, 1, total_process)

    call MPI_ALLGATHER(imx,1,MPI_INTEGER, IDim, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER(jmx,1,MPI_INTEGER, JDim, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER(kmx,1,MPI_INTEGER, KDim, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    IDimMax = maxval(IDim)
    JDimMax = maxval(JDim)
    KDimMax = maxval(KDim)

    SendSize = (imx-1)*(jmx-1)*(kmx-1)*n_var
    call alloc(RecvSize, 1, total_process)
    call alloc(Displacement, 1, total_process)
    call alloc(StartIndex, 1, total_process)
    RecvSize = 0
    do i = 1,total_process
      RecvSize(i) = (IDim(i)-1)*(JDim(i)-1)*(KDim(i)-1)*n_var
    end do
    Displacement(1) = 0
    do i = 2,total_process
      Displacement(i) = Displacement(i-1) + RecvSize(i-1)
    end do
    StartIndex = Displacement + 1
    GlobalSize = Sum(RecvSize)

    call alloc(SendBuffer, 1, SendSize)
    call alloc(RecvBuffer, 1, GlobalSize)
    call alloc(Globalx,    1, GlobalSize)
    call alloc(Globalv1,   1,GlobalSize)
    call alloc(wj,         1,GlobalSize)
    call alloc(vi,         1,GlobalSize)
    call alloc(Globalvj,   1,GlobalSize)
    Globalx = 0.0
    GlobalV1 =0.0
    wj = 0.0
    vi = 0.0

      call alloc(Avj, 1,imx-1,1,jmx-1,1,kmx-1,1,n_var)
      call alloc(MinvAvj, 1,imx-1,1,jmx-1,1,kmx-1,1,n_var)
      !call alloc(wj, 1,imx-1,1,jmx-1,1,kmx-1,1,n_var)
      !call alloc(vi, 1,imx-1,1,jmx-1,1,kmx-1,1,n_var)
      call alloc(vj, 1,imx-1,1,jmx-1,1,kmx-1,1,n_var)

      call alloc(r0, 1,imx-1,1,jmx-1,1,kmx-1,1,n_var)
      call alloc(Minvr0, 1,imx-1,1,jmx-1,1,kmx-1,1,n_var)
      call alloc(Ax0, 1,imx-1,1,jmx-1,1,kmx-1,1,n_var)
      call alloc(Ax, 1,imx-1,1,jmx-1,1,kmx-1,1,n_var)
      call alloc(r, 1,imx-1,1,jmx-1,1,kmx-1,1,n_var)
      call alloc(v1, 1,imx-1,1,jmx-1,1,kmx-1,1,n_var)

      !call alloc(kspvectors, 1,IDimMax-1,1,JDimMax-1,1,KDimMax-1,1,n_var,1,total_process,1,m+1)
      call alloc(kspvectors, 1,GlobalSize,1,m+1)
      Kspvectors = 0.0

      call alloc(old_state, -2,imx+2,-2,jmx+2,-2,kmx+2,1,n_var)
      call alloc(conser, -2,imx+2,-2,jmx+2,-2,kmx+2,1,n_var)
      call alloc(old_residue, 1,imx-1,1,jmx-1,1,kmx-1,1,n_var)

      call alloc(delQ, 0, imx, 0, jmx, 0, kmx, 1, n_var)
      call alloc(delQstar, 0, imx, 0, jmx, 0, kmx, 1, n_var)

      if(mu_ref==0.0 .or. turbulence=='none') then
        call alloc(dummy, 0, imx, 0, jmx, 0, kmx)
        dummy = 0.0
      end if
      if(mu_ref==0.0)then
        mmu => dummy
      else
        mmu => mu
      end if
      if(trim(turbulence)=='none')then
        tmu => dummy
      else
        tmu => mu_t
      end if

  end subroutine setup_solvergmres

  subroutine destroy_solvergmres()
    implicit none

      call dealloc(SendBuffer)
      call dealloc(RecvBuffer)
      call dealloc(Globalx)
      call dealloc(Globalv1)

      call dealloc(Avj)
      call dealloc(MinvAvj)
      call dealloc(wj)
      call dealloc(vi)
      call dealloc(vj)

      call dealloc(r0)
      call dealloc(Minvr0)
      call dealloc(Ax0)
      call dealloc(Ax)
      call dealloc(r)
      call dealloc(v1)

      call dealloc(kspvectors)

      call dealloc(old_state)
      call dealloc(conser)
      call dealloc(old_residue)

      call dealloc(delQ)
      call dealloc(delQstar)
      call dealloc(dummy)
  end subroutine destroy_solvergmres

  !===================================================
!  subroutine findnorm(norm,v1)
!
!    real, dimension(:,:,:,:), intent(in) :: v1
!    real, intent(out) :: norm
!
!    norm = SUM(v1**2)
!
!    norm = sqrt(norm)
!
!  end subroutine findnorm
!  !===================================================
!  subroutine innerproduct(v1dotv2,v1,v2)
!
!    real, dimension(:,:,:,:), intent(in) :: v1
!    real, dimension(:,:,:,:), intent(in) :: v2
!    real, intent(out) :: v1dotv2
!
!    v1dotv2 = SUM(v1*v2)
!
!  end subroutine innerproduct
  !===================================================
  subroutine arnoldialgorithm(v1,m,Hmat,kspvecs,lucky)

      integer,intent(in) :: m
      real, dimension(:), intent(in) :: v1
      real, intent(inout) :: Hmat(m+1,m)
      real, dimension(:,:), intent(inout) :: kspvecs
      logical,intent(inout) :: lucky


      !real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var) :: Avj
      !real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var) :: MinvAvj
      !real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var) :: wj
      !real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var) :: vi,vj
      integer :: i,j,k

      kspvecs(:,1)=v1
      lucky = .false.

      do j=1,m
      
          GlobalVj = kspvecs(:,j)
          call GetLocalData(GlobalVj, vj, StartIndex(process_id+1))

      call findAX(Avj,vj)
      call precond(MinvAvj,Avj)
      !MinvAvj = Avj
      !paralle call for MinvAvj
      call ParallelCommunication(wj,MinvAvj)

      do i=1,j
        ! vi should to global
        vi = kspvecs(:,i) 
        !call innerproduct(Hmat(i,j),MinvAvj,vi)
        !call innerproduct(Hmat(i,j),wj,vi)
        Hmat(i,j) = SUM(wj*vi)
      enddo

      !wj=MinvAvj
      ! check if can remove this assignment and work only with one variable
      ! make Wj as global variable

      do i=1,j
        vi = kspvecs(:,i)
        wj=wj-Hmat(i,j)*vi
      enddo

      !call findnorm(Hmat(j+1,j),wj)
      Hmat(j+1,j) = Sqrt(SUM(wj**2))

      if(Hmat(j+1,j) > 0.d0) then
        kspvecs(:,j+1)=wj(:)/Hmat(j+1,j)
      else
        lucky=.true.
      endif
    enddo


  end subroutine arnoldialgorithm
  !===================================================
  subroutine leastsquaresminimize(y,Hmat,m,beta)

      integer,intent(in)   :: m
      real,intent(inout) :: Hmat(m+1,m)
      real,intent(out)   :: y(m)
      real,intent(in)    :: beta

      real :: c,s,h_up,h_down,dtr
      real :: val1,val2
      real :: beta_e1(m+1)

      integer :: i,j

      beta_e1(:) = 0.d0
      beta_e1(1) = beta;

      do i=1,m

        h_up   = Hmat(i,i)
      h_down = Hmat(i+1,i)

      dtr = sqrt(h_up*h_up+h_down*h_down)

      c = h_up/dtr
      s = h_down/dtr

      do j=1,m
        
        h_up   = Hmat(i,j)
        h_down = Hmat(i+1,j)

        Hmat(i,j)   =  c*h_up  + s*h_down
        Hmat(i+1,j) = -s*h_up  + c*h_down 

      enddo

      val1 =  c*beta_e1(i)  + s*beta_e1(i+1); 
      val2 = -s*beta_e1(i)  + c*beta_e1(i+1);

      beta_e1(i)   = val1
      beta_e1(i+1) = val2 
    enddo


    y(m) = beta_e1(m)/Hmat(m,m)

    do i=m-1,1,-1

      y(i)=beta_e1(i)

      do j=i+1,m
        y(i)=y(i)-Hmat(i,j)*y(j)
      enddo

      y(i) = y(i)/Hmat(i,i)

    enddo

  end subroutine leastsquaresminimize
  !===================================================
  subroutine performgmres(b,x0,x,m,nrestarts)

      integer,intent(in) :: m,nrestarts
      real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(in) :: b, x0
      real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(out) :: x

      integer :: i,j,k

      !real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var) :: r0
      !real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var) :: Minvr0 
      !real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var) :: Ax0, Ax
      !real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var) :: r, v1
      real :: beta
      real :: y(m)

      !real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var,1:m+1) :: kspvectors
      real :: Hmat(m+1,m)

      logical :: lucky
      integer :: ierr

      call findAX(Ax0,x0)
      r0 = b-Ax0
      call precond(Minvr0,r0)
      r=Minvr0
      !r = r0
      x=x0
      call ParallelCommunication(Globalx, x)

      Hmat       = 0.d0
      kspvectors = 0.d0
      lucky      = .false.

      do i=1,nrestarts
        
        !call findnorm(beta,r)
        beta = Sum(r**2)
        ! finding global norm of residue
        call MPI_ALLREDUCE(beta, NormSum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        beta = sqrt(NormSum)
      if(process_id == 0 ) then 
        print *,"restart iteration:",i,"residual norm:",beta
      end if

      v1=r/beta
      
      ! send local V1 data to globalx
      call ParallelCommunication(Globalv1, v1)

        call arnoldialgorithm(Globalv1,m,Hmat,kspvectors,lucky)

      !if(lucky .eqv. .true.) then
      !  print *,"lucky condition"
      !  exit
      !endif

        call leastsquaresminimize(y,Hmat,m,beta)

      do j=1,m
        Globalx=Globalx+y(j)*kspvectors(:,j)
      enddo

      call GetLocalData(Globalx, x, StartIndex(process_id+1))

      call findAX(Ax,x)
      r=b-Ax
      call precond(Minvr0,r)
      r=Minvr0
    enddo

    call updateQp(x)
        
        
        
  end subroutine performgmres
  !===================================================


  subroutine GetLocalData(Globalvector, LocalMatrix, StartIndex)
    implicit none
    integer, intent(in) :: StartIndex
    real, dimension(:), intent(in) :: GlobalVector
    real, dimension(:,:,:,:), intent(out) :: LocalMatrix

    integer :: i,j,k,l, count

    count = StartIndex
    do l = 1,n_var
      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1
            LocalMatrix(i,j,k,l) = GlobalVector(count)
            count = count + 1
          end do
        end do
      end do
    end do

  end subroutine GetLocalData


  subroutine ParallelCommunication(Global, local)
    implicit none
    real, dimension(:), intent(out):: Global
    real, dimension(:,:,:,:)  , intent(in) :: local

    integer :: i,j,k,l,m
    integer :: count
    integer :: ierr

    ! collect data
    count = 0
    do l = 1,n_var
      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1
            count = count+1
            SendBuffer(count) = local(i,j,k,l)
          end do
        end do
      end do
    end do
    !print*, "Send:", SendSize
    !print*, "Recv:", RecvSize
    !print*, "Disp", Displacement

    !send data
    call MPI_ALLGATHERv(SendBuffer, SendSize, MPI_DOUBLE_PRECISION, Global, RecvSize, Displacement, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)

!    ! relocate data
!    Global = 0.0
!    count = 0
!    do m = 1,total_process
!      do l = 1,n_var
!        do k = 1,KDim(m)-1
!          do j = 1,JDim(m)-1
!            do i = 1,IDim(m)-1
!              count = count+1
!              Global(count) = RecvBuffer(count)
!            end do
!          end do
!        end do
!      end do
!    end do
  end subroutine ParallelCommunication

  subroutine findAX(Ax,x)

    implicit none
    real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(in) :: x
    real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(out) :: Ax
    !real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var) :: old_state
    !real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var) :: conser
    !real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var) :: old_residue
    real :: sigma
    real :: norm
    real :: Emach=1e-14
    integer :: i,j,k

    old_state = qp
    old_residue = residue
    !call innerproduct(sigma,qp,x)
    sigma = SUM(qp*x)
    !call findnorm(norm, x)
    norm = sqrt(SUM(x**2))
    if(norm>0 .and. abs(sigma)>0)then

      !sigma = sqrt(Emach)*sigma/norm
      sigma = sqrt(Emach)/norm
      call find_conservative_state(conser, qp)
      conser(1:imx-1,1:jmx-1,1:kmx-1,1:n_var) =  conser(1:imx-1,1:jmx-1,1:kmx-1,1:n_var) + sigma*x
      call find_primitive_state(conser, qp)
      !if(process_id==0)then
      !  print *, "variables"
      !  do k = -2,kmx+2
      !    do j = -2,jmx+2
      !      do i = -2,imx+2
      !        print *, i,j,k,sigma, qp(i,j,k,1)
      !      end do
      !    end do
      !  end do
      !endif

        call apply_interface()
        call populate_ghost_primitive()
        call compute_face_interpolant()
        call reconstruct_boundary_state(interpolant)
        call compute_fluxes()
        if (mu_ref /= 0.0) then
          call evaluate_all_gradients()
          call calculate_viscosity()
          call compute_viscous_fluxes(F_p, G_p, H_p)
        end if
        call compute_residue()
        call add_source_term_residue()

        residue = -residue

        Ax = -(residue - old_residue)/sigma
      else
        Ax=0.0
      end if
      qp = old_state
      residue = old_residue

  end subroutine findAX


  subroutine findAX1(Ax,x)

    implicit none
    real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(in) :: x
    real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(out) :: Ax
    real :: sigma(n_var)
    real :: norm(n_var)
    real :: Emach=1e-6
    integer :: i

    old_state = qp
    old_residue = residue

    call find_conservative_state(conser, qp)

    do i = 1,n_var
      sigma(i) = SUM(qp(1:imx-1,1:jmx-1,1:kmx-1,i)*x(1:imx-1,1:jmx-1,1:kmx-1,i))
      norm(i) = Sqrt(SUM(x(1:imx-1,1:jmx-1,1:kmx-1,:)**2))
      if(norm(i)>0.0) then
        !sigma(i) = Sqrt(Emach)*sigma(i)/norm(i)
        sigma(i) = Sqrt(Emach)/norm(i)
      else 
        sigma(i) = 0.0
      end if
      conser(1:imx-1,1:jmx-1,1:kmx-1,i) =  conser(1:imx-1,1:jmx-1,1:kmx-1,i)  &
                                        + sigma(i)*x(1:imx-1,1:jmx-1,1:kmx-1,i)
    end do

    call find_primitive_state(conser, qp)

    !compute residue
        call apply_interface()
        call populate_ghost_primitive()
        call compute_face_interpolant()
        !call reconstruct_boundary_state(interpolant)
        call reconstruct_boundary_state('none')
        call compute_fluxes()
        if (mu_ref /= 0.0) then
          !check if not recaculating gradient increase speed
          call evaluate_all_gradients()
          call calculate_viscosity()
          call compute_viscous_fluxes(F_p, G_p, H_p)
        end if
        call compute_residue()
        call add_source_term_residue()

        residue = -residue

      do i = 1,n_var
        if(sigma(i)>0.0)then
          Ax(1:imx-1,1:jmx-1,1:kmx-1,i) = -(residue(1:imx-1,1:jmx-1,1:kmx-1,i) - old_residue(1:imx-1,1:jmx-1,1:kmx-1,i))/sigma(i)
        else
          Ax(1:imx-1,1:jmx-1,1:kmx-1,i) = 0.0
        end if
      end do

      qp = old_state
      residue = old_residue

  end subroutine findAX1



  subroutine find_conservative_state(Q,P)
    implicit none
    real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(in) :: P
    real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(out) :: Q

    integer :: i,j,k
    real, dimension(1:n_var)             :: U ! conservative variables
    real, dimension(1:n_var)             :: W ! new primitive variables

    do k = 1,kmx-1
      do j = 1,jmx-1
        do i = 1,imx-1
          W(1:n_var) = P(i,j,k,1:n_var)

          U(1)   =   W(1)
          U(2:n_var)   =   W(1) * W(2:n_var)
          U(5)   = ( W(5) / (gm-1.0) ) + ( 0.5 * W(1) * sum(W(2:4)**2) )

          Q(i,j,k,1:n_var) = U(1:n_var)

        end do
      end do
    end do

    
  end subroutine find_conservative_state


  subroutine find_primitive_state(Q,P)
    implicit none
    real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(out) :: P
    real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(in) :: Q

    integer :: i,j,k
    real, dimension(1:n_var)             :: U ! conservative variables
    real, dimension(1:n_var)             :: W ! new primitive variables

    do k = 1,kmx-1
      do j = 1,jmx-1
        do i = 1,imx-1
          U(1:n_var) = Q(i,j,k,1:n_var)

          W(1)   =   U(1)
          W(2:n_var)   =   U(2:n_var) / U(1)
          W(5)   = (gm-1.0) * ( U(5) - ( 0.5 * SUM(U(2:4)**2) / U(1) ) )

          P(i,j,k,1:n_var) = W(1:n_var)

        end do
      end do
    end do

    
  end subroutine find_primitive_state


    subroutine precond(output, residue)
      implicit none
      real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(in) :: residue
      real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(out) :: output
      !real, dimension(0:imx,0:jmx,0:kmx,1:n_var) :: delQ
      !real, dimension(0:imx,0:jmx,0:kmx,1:n_var) :: delQstar
      !real, dimension(0:imx,0:jmx,0:kmx), target :: dummy
      !real, dimension(:,:,:), pointer :: tmu
      !real, dimension(:,:,:), pointer :: mmu
      
      integer :: i,j,k, SubIter
        real, dimension(1:5)     :: deltaU
        real                     :: D
        real, dimension(1:5)     :: LOldIminusFlux
        real, dimension(1:5)     :: LOldJminusFlux
        real, dimension(1:5)     :: LOldKminusFlux
        real, dimension(1:5)     :: LNewIminusFlux
        real, dimension(1:5)     :: LNewJminusFlux
        real, dimension(1:5)     :: LNewKminusFlux
        real, dimension(1:5)     :: LDelIminusFlux
        real, dimension(1:5)     :: LDelJminusFlux
        real, dimension(1:5)     :: LDelKminusFlux
        real, dimension(1:5)     :: UOldIminusFlux
        real, dimension(1:5)     :: UOldJminusFlux
        real, dimension(1:5)     :: UOldKminusFlux
        real, dimension(1:5)     :: UNewIminusFlux
        real, dimension(1:5)     :: UNewJminusFlux
        real, dimension(1:5)     :: UNewKminusFlux
        real, dimension(1:5)     :: UDelIminusFlux
        real, dimension(1:5)     :: UDelJminusFlux
        real, dimension(1:5)     :: UDelKminusFlux
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

              Q0  = qp(i  , j  , k  , 1:5)
              Q1  = qp(i-1, j  , k  , 1:5)
              Q2  = qp(i  , j-1, k  , 1:5)
              Q3  = qp(i  , j  , k-1, 1:5)
              Q4  = qp(i+1, j  , k  , 1:5)
              Q5  = qp(i  , j+1, k  , 1:5)
              Q6  = qp(i  , j  , k+1, 1:5)


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


              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0)

              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
              !storing D in Iflux array for backward sweep
              !F_p(i,j,k,1) = D

              deltaU(1:5) = -residue(i,j,k,1:5)

              delQstar(i,j,k,1:5) = deltaU(1:5)/D

            end do
          end do
        end do

        do SubIter = 1,1
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
              DQ4 = delQstar(i+1, j  , k  , 1:5)
              DQ5 = delQstar(i  , j+1, k  , 1:5)
              DQ6 = delQstar(i  , j  , k+1, 1:5)

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

              LNewIminusFlux     = Flux(Q1, Q0, DQ1, Flist1)
              LNewJminusFlux     = Flux(Q2, Q0, DQ2, Flist2)
              LNewKminusFlux     = Flux(Q3, Q0, DQ3, Flist3)
              LOldIminusFlux     = Flux(Q1, Q0, DQ0, Flist1)
              LOldJminusFlux     = Flux(Q2, Q0, DQ0, Flist2)
              LOldKminusFlux     = Flux(Q3, Q0, DQ0, Flist3)
              UNewIminusFlux     = Flux(Q4, Q0, DQ4, Flist4)
              UNewJminusFlux     = Flux(Q5, Q0, DQ5, Flist5)
              UNewKminusFlux     = Flux(Q6, Q0, DQ6, Flist6)
              UOldIminusFlux     = Flux(Q4, Q0, DQ0, Flist4)
              UOldJminusFlux     = Flux(Q5, Q0, DQ0, Flist5)
              UOldKminusFlux     = Flux(Q6, Q0, DQ0, Flist6)

              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0)


              ! multiply above flux with area to get correct values
              LDelIminusFlux =  LNewIminusFlux - LOldIminusFlux
              LDelJminusFlux =  LNewJminusFlux - LOldJminusFlux
              LDelKminusFlux =  LNewKminusFlux - LOldKminusFlux

              ! multiply above flux with area to get correct values
              UDelIminusFlux =  UNewIminusFlux - UOldIminusFlux
              UDelJminusFlux =  UNewJminusFlux - UOldJminusFlux
              UDelKminusFlux =  UNewKminusFlux - UOldKminusFlux

              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)

              delQ(i,j,k,1:5) = (-residue(i,j,k,1:5) &
                - 0.5*((LDelIminusFlux + UDelIminusFlux) &
                     + (LDelJminusFlux + UDelJminusFlux) &
                     + (LDelKminusFlux + UDelKminusFlux)))/D

              output(i,j,k,1:5) = delQ(i,j,k,1:5)
            end do
          end do
        end do
        end do
        

    end subroutine precond


    subroutine precond2(output, residue)
      implicit none
      real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(in) :: residue
      real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(out) :: output
      !real, dimension(0:imx,0:jmx,0:kmx,1:n_var) :: delq
      !real, dimension(0:imx,0:jmx,0:kmx,1:n_var) :: delqstar
      !real, dimension(0:imx,0:jmx,0:kmx), target :: dummy
      !real, dimension(:,:,:), pointer :: tmu
      !real, dimension(:,:,:), pointer :: mmu
      
      integer :: i,j,k
        real, dimension(1:5)     :: deltau
        real                     :: d
        real, dimension(1:5)     :: LOldIminusFlux
        real, dimension(1:5)     :: LOldJminusFlux
        real, dimension(1:5)     :: LOldKminusFlux
        real, dimension(1:5)     :: LNewIminusFlux
        real, dimension(1:5)     :: LNewJminusFlux
        real, dimension(1:5)     :: LNewKminusFlux
        real, dimension(1:5)     :: LDelIminusFlux
        real, dimension(1:5)     :: LDelJminusFlux
        real, dimension(1:5)     :: LDelKminusFlux
        real, dimension(1:5)     :: UOldIminusFlux
        real, dimension(1:5)     :: UOldJminusFlux
        real, dimension(1:5)     :: UOldKminusFlux
        real, dimension(1:5)     :: UNewIminusFlux
        real, dimension(1:5)     :: UNewJminusFlux
        real, dimension(1:5)     :: UNewKminusFlux
        real, dimension(1:5)     :: UDelIminusFlux
        real, dimension(1:5)     :: UDelJminusFlux
        real, dimension(1:5)     :: UDelKminusFlux
        real, dimension(1:5)     :: oldiminusflux
        real, dimension(1:5)     :: oldjminusflux
        real, dimension(1:5)     :: oldkminusflux
        real, dimension(1:5)     :: newiminusflux
        real, dimension(1:5)     :: newjminusflux
        real, dimension(1:5)     :: newkminusflux
        real, dimension(1:5)     :: deliminusflux
        real, dimension(1:5)     :: deljminusflux
        real, dimension(1:5)     :: delkminusflux
        real, dimension(1:6)     :: lambdatimesarea
        real, dimension(1:5)     :: q0 ! state at cell
        real, dimension(1:5)     :: q1 ! state at neighbours 
        real, dimension(1:5)     :: q2
        real, dimension(1:5)     :: q3
        real, dimension(1:5)     :: q4
        real, dimension(1:5)     :: q5
        real, dimension(1:5)     :: q6
        real, dimension(1:5)     :: dq0! change in state
        real, dimension(1:5)     :: dq1
        real, dimension(1:5)     :: dq2
        real, dimension(1:5)     :: dq3
        real, dimension(1:5)     :: dq4
        real, dimension(1:5)     :: dq5
        real, dimension(1:5)     :: dq6
        real, dimension(1:7)     :: flist1
        real, dimension(1:7)     :: flist2
        real, dimension(1:7)     :: flist3
        real, dimension(1:7)     :: flist4
        real, dimension(1:7)     :: flist5
        real, dimension(1:7)     :: flist6
        real, dimension(1:3)     :: c0
        real, dimension(1:3)     :: c1
        real, dimension(1:3)     :: c2
        real, dimension(1:3)     :: c3
        real, dimension(1:3)     :: c4
        real, dimension(1:3)     :: c5
        real, dimension(1:3)     :: c6


      !if(mu_ref==0.0 .or. turbulence=='none') then
      !  dummy = 0.0
      !end if
      !if(mu_ref==0.0)then
      !  mmu => dummy
      !else
      !  mmu => mu
      !end if
      !if(trim(turbulence)=='none')then
      !  tmu => dummy
      !else
      !  tmu => mu_t
      !end if


        !intialize delq
        delqstar = 0.0

        !forward sweep
        do k=0,kmx
          do j=0,jmx
            do i=0,imx
              C0  = CellCenter(i  ,j  ,k  ,:)
              C1  = CellCenter(i-1,j  ,k  ,:)
              C2  = CellCenter(i  ,j-1,k  ,:)
              C3  = CellCenter(i  ,j  ,k-1,:)
              C4  = CellCenter(i+1,j  ,k  ,:)
              C5  = CellCenter(i  ,j+1,k  ,:)
              C6  = CellCenter(i  ,j  ,k+1,:)

              Q0  = qp(i  , j  , k  , 1:5)
              Q1  = qp(i-1, j  , k  , 1:5)
              Q2  = qp(i  , j-1, k  , 1:5)
              Q3  = qp(i  , j  , k-1, 1:5)
              Q4  = qp(i+1, j  , k  , 1:5)
              Q5  = qp(i  , j+1, k  , 1:5)
              Q6  = qp(i  , j  , k+1, 1:5)


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


              LambdaTimesArea(1)= SpectralRadius(Q1, Q0, Flist1, C1, C0)
              LambdaTimesArea(2)= SpectralRadius(Q2, Q0, Flist2, C2, C0)
              LambdaTimesArea(3)= SpectralRadius(Q3, Q0, Flist3, C3, C0)
              LambdaTimesArea(4)= SpectralRadius(Q4, Q0, Flist4, C4, C0)
              LambdaTimesArea(5)= SpectralRadius(Q5, Q0, Flist5, C5, C0)
              LambdaTimesArea(6)= SpectralRadius(Q6, Q0, Flist6, C6, C0)

              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
              !storing D in Iflux array for backward sweep
              !F_p(i,j,k,1) = D

              deltaU(1:5) = -residue(i,j,k,1:5)

              delQstar(i,j,k,1:5) = deltaU(1:5)/D
              delQ(i,j,k,1:5) = deltaU(1:5)/D

            end do
          end do
        end do

        !forward sweep
        do k=1,kmx-1
          do j=1,jmx-1
            do i=1,imx-1
              c0  = cellcenter(i  ,j  ,k  ,:)
              c1  = cellcenter(i-1,j  ,k  ,:)
              c2  = cellcenter(i  ,j-1,k  ,:)
              c3  = cellcenter(i  ,j  ,k-1,:)
              c4  = cellcenter(i+1,j  ,k  ,:)
              c5  = cellcenter(i  ,j+1,k  ,:)
              c6  = cellcenter(i  ,j  ,k+1,:)

              q0  = qp(i  , j  , k  , 1:5)
              q1  = qp(i-1, j  , k  , 1:5)
              q2  = qp(i  , j-1, k  , 1:5)
              q3  = qp(i  , j  , k-1, 1:5)
              q4  = qp(i+1, j  , k  , 1:5)
              q5  = qp(i  , j+1, k  , 1:5)
              q6  = qp(i  , j  , k+1, 1:5)

              dq0 = 0.0
              dq1 = delqstar(i-1, j  , k  , 1:5)
              dq2 = delqstar(i  , j-1, k  , 1:5)
              dq3 = delqstar(i  , j  , k-1, 1:5)

              flist1(1) =   xa(i,j,k)
              flist1(2) = -xnx(i,j,k)
              flist1(3) = -xny(i,j,k)
              flist1(4) = -xnz(i,j,k)
              flist1(5) = 0.5*(volume(i-1, j  , k  ) + volume(i,j,k))
              flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))

              flist2(1) =   ya(i,j,k)
              flist2(2) = -ynx(i,j,k)
              flist2(3) = -yny(i,j,k)
              flist2(4) = -ynz(i,j,k)
              flist2(5) = 0.5*(volume(i  , j-1, k  ) + volume(i,j,k))
              flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))

              flist3(1) =   za(i,j,k)
              flist3(2) = -znx(i,j,k)
              flist3(3) = -zny(i,j,k)
              flist3(4) = -znz(i,j,k)
              flist3(5) = 0.5*(volume(i  , j  , k-1) + volume(i,j,k))
              flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))

              flist4(1) =   xa(i+1,j,k)
              flist4(2) = +xnx(i+1,j,k)
              flist4(3) = +xny(i+1,j,k)
              flist4(4) = +xnz(i+1,j,k)
              flist4(5) = 0.5*(volume(i+1, j  , k  ) + volume(i,j,k))
              flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))

              flist5(1) =   ya(i,j+1,k)
              flist5(2) = +ynx(i,j+1,k)
              flist5(3) = +yny(i,j+1,k)
              flist5(4) = +ynz(i,j+1,k)
              flist5(5) = 0.5*(volume(i  , j+1, k  ) + volume(i,j,k))
              flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))

              flist6(1) =   za(i,j,k+1)
              flist6(2) = +znx(i,j,k+1)
              flist6(3) = +zny(i,j,k+1)
              flist6(4) = +znz(i,j,k+1)
              flist6(5) = 0.5*(volume(i  , j  , k+1) + volume(i,j,k))
              flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))

              newiminusflux     = flux(q1, q0, dq1, flist1)
              newjminusflux     = flux(q2, q0, dq2, flist2)
              newkminusflux     = flux(q3, q0, dq3, flist3)
              oldiminusflux     = flux(q1, q0, dq0, flist1)
              oldjminusflux     = flux(q2, q0, dq0, flist2)
              oldkminusflux     = flux(q3, q0, dq0, flist3)

              lambdatimesarea(1)= spectralradius(q1, q0, flist1, c1, c0)
              lambdatimesarea(2)= spectralradius(q2, q0, flist2, c2, c0)
              lambdatimesarea(3)= spectralradius(q3, q0, flist3, c3, c0)
              lambdatimesarea(4)= spectralradius(q4, q0, flist4, c4, c0)
              lambdatimesarea(5)= spectralradius(q5, q0, flist5, c5, c0)
              lambdatimesarea(6)= spectralradius(q6, q0, flist6, c6, c0)


              ! multiply above flux with area to get correct values
              deliminusflux =  newiminusflux - oldiminusflux
              deljminusflux =  newjminusflux - oldjminusflux
              delkminusflux =  newkminusflux - oldkminusflux


              d = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*sum(lambdatimesarea)
              !storing d in iflux array for backward sweep
              !f_p(i,j,k,1) = d

              deltau(1:5) = +residue(i,j,k,1:5) &
                - 0.5*((deliminusflux - lambdatimesarea(1)*delqstar(i-1,j,k,1:5)) &
                     + (deljminusflux - lambdatimesarea(2)*delqstar(i,j-1,k,1:5)) &
                     + (delkminusflux - lambdatimesarea(3)*delqstar(i,j,k-1,1:5)) )

              delqstar(i,j,k,1:5) = deltau(1:5)/d
            end do
          end do
        end do

        delq=0.0
        !backward sweep
            do i=imx-1,1,-1
          do j=jmx-1,1,-1
        do k=kmx-1,1,-1
              c0  = cellcenter(i  ,j  ,k  ,:)
              c1  = cellcenter(i-1,j  ,k  ,:)
              c2  = cellcenter(i  ,j-1,k  ,:)
              c3  = cellcenter(i  ,j  ,k-1,:)
              c4  = cellcenter(i+1,j  ,k  ,:)
              c5  = cellcenter(i  ,j+1,k  ,:)
              c6  = cellcenter(i  ,j  ,k+1,:)

              q0  = qp(i  , j  , k  , 1:5)
              q1  = qp(i-1, j  , k  , 1:5)
              q2  = qp(i  , j-1, k  , 1:5)
              q3  = qp(i  , j  , k-1, 1:5)
              q4  = qp(i+1, j  , k  , 1:5)
              q5  = qp(i  , j+1, k  , 1:5)
              q6  = qp(i  , j  , k+1, 1:5)

              dq0 = 0.0
              dq4 = delq(i+1, j  , k  , 1:5)
              dq5 = delq(i  , j+1, k  , 1:5)
              dq6 = delq(i  , j  , k+1, 1:5)

              flist1(1) =   xa(i,j,k)
              flist1(2) = -xnx(i,j,k)
              flist1(3) = -xny(i,j,k)
              flist1(4) = -xnz(i,j,k)
              flist1(5) = 0.5*(volume(i-1, j  , k  ) + volume(i,j,k))
              flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))

              flist2(1) =   ya(i,j,k)
              flist2(2) = -ynx(i,j,k)
              flist2(3) = -yny(i,j,k)
              flist2(4) = -ynz(i,j,k)
              flist2(5) = 0.5*(volume(i  , j-1, k  ) + volume(i,j,k))
              flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))

              flist3(1) =   za(i,j,k)
              flist3(2) = -znx(i,j,k)
              flist3(3) = -zny(i,j,k)
              flist3(4) = -znz(i,j,k)
              flist3(5) = 0.5*(volume(i  , j  , k-1) + volume(i,j,k))
              flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))

              flist4(1) =   xa(i+1,j,k)
              flist4(2) = +xnx(i+1,j,k)
              flist4(3) = +xny(i+1,j,k)
              flist4(4) = +xnz(i+1,j,k)
              flist4(5) = 0.5*(volume(i+1, j  , k  ) + volume(i,j,k))
              flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))

              flist5(1) =   ya(i,j+1,k)
              flist5(2) = +ynx(i,j+1,k)
              flist5(3) = +yny(i,j+1,k)
              flist5(4) = +ynz(i,j+1,k)
              flist5(5) = 0.5*(volume(i  , j+1, k  ) + volume(i,j,k))
              flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))

              flist6(1) =   za(i,j,k+1)
              flist6(2) = +znx(i,j,k+1)
              flist6(3) = +zny(i,j,k+1)
              flist6(4) = +znz(i,j,k+1)
              flist6(5) = 0.5*(volume(i  , j  , k+1) + volume(i,j,k))
              flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))

              newiminusflux     = flux(q4, q0, dq4, flist4)
              newjminusflux     = flux(q5, q0, dq5, flist5)
              newkminusflux     = flux(q6, q0, dq6, flist6)
              oldiminusflux     = flux(q4, q0, dq0, flist4)
              oldjminusflux     = flux(q5, q0, dq0, flist5)
              oldkminusflux     = flux(q6, q0, dq0, flist6)

              lambdatimesarea(1)= spectralradius(q1, q0, flist1, c1, c0)
              lambdatimesarea(2)= spectralradius(q2, q0, flist2, c2, c0)
              lambdatimesarea(3)= spectralradius(q3, q0, flist3, c3, c0)
              lambdatimesarea(4)= spectralradius(q4, q0, flist4, c4, c0)
              lambdatimesarea(5)= spectralradius(q5, q0, flist5, c5, c0)
              lambdatimesarea(6)= spectralradius(q6, q0, flist6, c6, c0)


              ! multiply above flux with area to get correct values
              deliminusflux =  newiminusflux - oldiminusflux
              deljminusflux =  newjminusflux - oldjminusflux
              delkminusflux =  newkminusflux - oldkminusflux

              d = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*sum(lambdatimesarea)

              delq(i,j,k,1:5) = delqstar(i,j,k,1:5) &
                - 0.5*((deliminusflux - lambdatimesarea(4)*delq(i+1,j,k,1:5)) &
                     + (deljminusflux - lambdatimesarea(5)*delq(i,j+1,k,1:5)) &
                     + (delkminusflux - lambdatimesarea(6)*delq(i,j,k+1,1:5)) )/d

              output(i,j,k,1:5) = delq(i,j,k,1:5)
            end do
          end do
        end do
        

    end subroutine precond2


    subroutine precond1(output, residue)
      implicit none
      real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(in) :: residue
      real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(out) :: output
      !real, dimension(0:imx,0:jmx,0:kmx,1:n_var) :: delq
      !real, dimension(0:imx,0:jmx,0:kmx,1:n_var) :: delqstar
      !real, dimension(0:imx,0:jmx,0:kmx), target :: dummy
      !real, dimension(:,:,:), pointer :: tmu
      !real, dimension(:,:,:), pointer :: mmu
      
      integer :: i,j,k
        real, dimension(1:5)     :: deltau
        real                     :: d
        real, dimension(1:5)     :: oldiminusflux
        real, dimension(1:5)     :: oldjminusflux
        real, dimension(1:5)     :: oldkminusflux
        real, dimension(1:5)     :: newiminusflux
        real, dimension(1:5)     :: newjminusflux
        real, dimension(1:5)     :: newkminusflux
        real, dimension(1:5)     :: deliminusflux
        real, dimension(1:5)     :: deljminusflux
        real, dimension(1:5)     :: delkminusflux
        real, dimension(1:6)     :: lambdatimesarea
        real, dimension(1:5)     :: q0 ! state at cell
        real, dimension(1:5)     :: q1 ! state at neighbours 
        real, dimension(1:5)     :: q2
        real, dimension(1:5)     :: q3
        real, dimension(1:5)     :: q4
        real, dimension(1:5)     :: q5
        real, dimension(1:5)     :: q6
        real, dimension(1:5)     :: dq0! change in state
        real, dimension(1:5)     :: dq1
        real, dimension(1:5)     :: dq2
        real, dimension(1:5)     :: dq3
        real, dimension(1:5)     :: dq4
        real, dimension(1:5)     :: dq5
        real, dimension(1:5)     :: dq6
        real, dimension(1:7)     :: flist1
        real, dimension(1:7)     :: flist2
        real, dimension(1:7)     :: flist3
        real, dimension(1:7)     :: flist4
        real, dimension(1:7)     :: flist5
        real, dimension(1:7)     :: flist6
        real, dimension(1:3)     :: c0
        real, dimension(1:3)     :: c1
        real, dimension(1:3)     :: c2
        real, dimension(1:3)     :: c3
        real, dimension(1:3)     :: c4
        real, dimension(1:3)     :: c5
        real, dimension(1:3)     :: c6


      !if(mu_ref==0.0 .or. turbulence=='none') then
      !  dummy = 0.0
      !end if
      !if(mu_ref==0.0)then
      !  mmu => dummy
      !else
      !  mmu => mu
      !end if
      !if(trim(turbulence)=='none')then
      !  tmu => dummy
      !else
      !  tmu => mu_t
      !end if


        !intialize delq
        delqstar = 0.0

        !forward sweep
        do k=1,kmx-1
          do j=1,jmx-1
            do i=1,imx-1
              c0  = cellcenter(i  ,j  ,k  ,:)
              c1  = cellcenter(i-1,j  ,k  ,:)
              c2  = cellcenter(i  ,j-1,k  ,:)
              c3  = cellcenter(i  ,j  ,k-1,:)
              c4  = cellcenter(i+1,j  ,k  ,:)
              c5  = cellcenter(i  ,j+1,k  ,:)
              c6  = cellcenter(i  ,j  ,k+1,:)

              q0  = qp(i  , j  , k  , 1:5)
              q1  = qp(i-1, j  , k  , 1:5)
              q2  = qp(i  , j-1, k  , 1:5)
              q3  = qp(i  , j  , k-1, 1:5)
              q4  = qp(i+1, j  , k  , 1:5)
              q5  = qp(i  , j+1, k  , 1:5)
              q6  = qp(i  , j  , k+1, 1:5)

              dq0 = 0.0
              dq1 = delqstar(i-1, j  , k  , 1:5)
              dq2 = delqstar(i  , j-1, k  , 1:5)
              dq3 = delqstar(i  , j  , k-1, 1:5)

              flist1(1) =   xa(i,j,k)
              flist1(2) = -xnx(i,j,k)
              flist1(3) = -xny(i,j,k)
              flist1(4) = -xnz(i,j,k)
              flist1(5) = 0.5*(volume(i-1, j  , k  ) + volume(i,j,k))
              flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))

              flist2(1) =   ya(i,j,k)
              flist2(2) = -ynx(i,j,k)
              flist2(3) = -yny(i,j,k)
              flist2(4) = -ynz(i,j,k)
              flist2(5) = 0.5*(volume(i  , j-1, k  ) + volume(i,j,k))
              flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))

              flist3(1) =   za(i,j,k)
              flist3(2) = -znx(i,j,k)
              flist3(3) = -zny(i,j,k)
              flist3(4) = -znz(i,j,k)
              flist3(5) = 0.5*(volume(i  , j  , k-1) + volume(i,j,k))
              flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))

              flist4(1) =   xa(i+1,j,k)
              flist4(2) = +xnx(i+1,j,k)
              flist4(3) = +xny(i+1,j,k)
              flist4(4) = +xnz(i+1,j,k)
              flist4(5) = 0.5*(volume(i+1, j  , k  ) + volume(i,j,k))
              flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))

              flist5(1) =   ya(i,j+1,k)
              flist5(2) = +ynx(i,j+1,k)
              flist5(3) = +yny(i,j+1,k)
              flist5(4) = +ynz(i,j+1,k)
              flist5(5) = 0.5*(volume(i  , j+1, k  ) + volume(i,j,k))
              flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))

              flist6(1) =   za(i,j,k+1)
              flist6(2) = +znx(i,j,k+1)
              flist6(3) = +zny(i,j,k+1)
              flist6(4) = +znz(i,j,k+1)
              flist6(5) = 0.5*(volume(i  , j  , k+1) + volume(i,j,k))
              flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))

              newiminusflux     = flux(q1, q0, dq1, flist1)
              newjminusflux     = flux(q2, q0, dq2, flist2)
              newkminusflux     = flux(q3, q0, dq3, flist3)
              oldiminusflux     = flux(q1, q0, dq0, flist1)
              oldjminusflux     = flux(q2, q0, dq0, flist2)
              oldkminusflux     = flux(q3, q0, dq0, flist3)

              lambdatimesarea(1)= spectralradius(q1, q0, flist1, c1, c0)
              lambdatimesarea(2)= spectralradius(q2, q0, flist2, c2, c0)
              lambdatimesarea(3)= spectralradius(q3, q0, flist3, c3, c0)
              lambdatimesarea(4)= spectralradius(q4, q0, flist4, c4, c0)
              lambdatimesarea(5)= spectralradius(q5, q0, flist5, c5, c0)
              lambdatimesarea(6)= spectralradius(q6, q0, flist6, c6, c0)


              ! multiply above flux with area to get correct values
              deliminusflux =  newiminusflux - oldiminusflux
              deljminusflux =  newjminusflux - oldjminusflux
              delkminusflux =  newkminusflux - oldkminusflux


              d = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*sum(lambdatimesarea)
              !storing d in iflux array for backward sweep
              !f_p(i,j,k,1) = d

              deltau(1:5) = +residue(i,j,k,1:5) &
                - 0.5*((deliminusflux - lambdatimesarea(1)*delqstar(i-1,j,k,1:5)) &
                     + (deljminusflux - lambdatimesarea(2)*delqstar(i,j-1,k,1:5)) &
                     + (delkminusflux - lambdatimesarea(3)*delqstar(i,j,k-1,1:5)) )

              delqstar(i,j,k,1:5) = deltau(1:5)/d
            end do
          end do
        end do

        delq=0.0
        !backward sweep
            do i=imx-1,1,-1
          do j=jmx-1,1,-1
        do k=kmx-1,1,-1
              c0  = cellcenter(i  ,j  ,k  ,:)
              c1  = cellcenter(i-1,j  ,k  ,:)
              c2  = cellcenter(i  ,j-1,k  ,:)
              c3  = cellcenter(i  ,j  ,k-1,:)
              c4  = cellcenter(i+1,j  ,k  ,:)
              c5  = cellcenter(i  ,j+1,k  ,:)
              c6  = cellcenter(i  ,j  ,k+1,:)

              q0  = qp(i  , j  , k  , 1:5)
              q1  = qp(i-1, j  , k  , 1:5)
              q2  = qp(i  , j-1, k  , 1:5)
              q3  = qp(i  , j  , k-1, 1:5)
              q4  = qp(i+1, j  , k  , 1:5)
              q5  = qp(i  , j+1, k  , 1:5)
              q6  = qp(i  , j  , k+1, 1:5)

              dq0 = 0.0
              dq4 = delq(i+1, j  , k  , 1:5)
              dq5 = delq(i  , j+1, k  , 1:5)
              dq6 = delq(i  , j  , k+1, 1:5)

              flist1(1) =   xa(i,j,k)
              flist1(2) = -xnx(i,j,k)
              flist1(3) = -xny(i,j,k)
              flist1(4) = -xnz(i,j,k)
              flist1(5) = 0.5*(volume(i-1, j  , k  ) + volume(i,j,k))
              flist1(6) = 0.5*(   mmu(i-1, j  , k  ) +    mmu(i,j,k))
              flist1(7) = 0.5*(   tmu(i-1, j  , k  ) +    tmu(i,j,k))

              flist2(1) =   ya(i,j,k)
              flist2(2) = -ynx(i,j,k)
              flist2(3) = -yny(i,j,k)
              flist2(4) = -ynz(i,j,k)
              flist2(5) = 0.5*(volume(i  , j-1, k  ) + volume(i,j,k))
              flist2(6) = 0.5*(   mmu(i  , j-1, k  ) +    mmu(i,j,k))
              flist2(7) = 0.5*(   tmu(i  , j-1, k  ) +    tmu(i,j,k))

              flist3(1) =   za(i,j,k)
              flist3(2) = -znx(i,j,k)
              flist3(3) = -zny(i,j,k)
              flist3(4) = -znz(i,j,k)
              flist3(5) = 0.5*(volume(i  , j  , k-1) + volume(i,j,k))
              flist3(6) = 0.5*(   mmu(i  , j  , k-1) +    mmu(i,j,k))
              flist3(7) = 0.5*(   tmu(i  , j  , k-1) +    tmu(i,j,k))

              flist4(1) =   xa(i+1,j,k)
              flist4(2) = +xnx(i+1,j,k)
              flist4(3) = +xny(i+1,j,k)
              flist4(4) = +xnz(i+1,j,k)
              flist4(5) = 0.5*(volume(i+1, j  , k  ) + volume(i,j,k))
              flist4(6) = 0.5*(   mmu(i+1, j  , k  ) +    mmu(i,j,k))
              flist4(7) = 0.5*(   tmu(i+1, j  , k  ) +    tmu(i,j,k))

              flist5(1) =   ya(i,j+1,k)
              flist5(2) = +ynx(i,j+1,k)
              flist5(3) = +yny(i,j+1,k)
              flist5(4) = +ynz(i,j+1,k)
              flist5(5) = 0.5*(volume(i  , j+1, k  ) + volume(i,j,k))
              flist5(6) = 0.5*(   mmu(i  , j+1, k  ) +    mmu(i,j,k))
              flist5(7) = 0.5*(   tmu(i  , j+1, k  ) +    tmu(i,j,k))

              flist6(1) =   za(i,j,k+1)
              flist6(2) = +znx(i,j,k+1)
              flist6(3) = +zny(i,j,k+1)
              flist6(4) = +znz(i,j,k+1)
              flist6(5) = 0.5*(volume(i  , j  , k+1) + volume(i,j,k))
              flist6(6) = 0.5*(   mmu(i  , j  , k+1) +    mmu(i,j,k))
              flist6(7) = 0.5*(   tmu(i  , j  , k+1) +    tmu(i,j,k))

              newiminusflux     = flux(q4, q0, dq4, flist4)
              newjminusflux     = flux(q5, q0, dq5, flist5)
              newkminusflux     = flux(q6, q0, dq6, flist6)
              oldiminusflux     = flux(q4, q0, dq0, flist4)
              oldjminusflux     = flux(q5, q0, dq0, flist5)
              oldkminusflux     = flux(q6, q0, dq0, flist6)

              lambdatimesarea(1)= spectralradius(q1, q0, flist1, c1, c0)
              lambdatimesarea(2)= spectralradius(q2, q0, flist2, c2, c0)
              lambdatimesarea(3)= spectralradius(q3, q0, flist3, c3, c0)
              lambdatimesarea(4)= spectralradius(q4, q0, flist4, c4, c0)
              lambdatimesarea(5)= spectralradius(q5, q0, flist5, c5, c0)
              lambdatimesarea(6)= spectralradius(q6, q0, flist6, c6, c0)


              ! multiply above flux with area to get correct values
              deliminusflux =  newiminusflux - oldiminusflux
              deljminusflux =  newjminusflux - oldjminusflux
              delkminusflux =  newkminusflux - oldkminusflux

              d = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*sum(lambdatimesarea)

              delq(i,j,k,1:5) = delqstar(i,j,k,1:5) &
                - 0.5*((deliminusflux - lambdatimesarea(4)*delq(i+1,j,k,1:5)) &
                     + (deljminusflux - lambdatimesarea(5)*delq(i,j+1,k,1:5)) &
                     + (delkminusflux - lambdatimesarea(6)*delq(i,j,k+1,1:5)) )/d

              output(i,j,k,1:5) = delq(i,j,k,1:5)
            end do
          end do
        end do
        

    end subroutine precond1


    function Flux(ql, qr, du, inputs)
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



    subroutine updateQp(delQ)
      implicit none
      real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(inout) :: delQ
      integer :: i,j,k
      real, dimension(1:n_var) :: conservativeQ
        do k=1,kmx-1
          do j = 1,jmx-1
            do i = 1,imx-1
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

              !write(*,'(5ES27.16e3)') qp(i,j,k,:)
            end do
          end do
        end do
    end subroutine updateQp
end module solvergmres

