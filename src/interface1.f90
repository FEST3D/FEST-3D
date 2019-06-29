!< This module handles the MPI Communication calls for interface boundary conditions
module interface1
  !< This module handles the MPI Communication calls for interface boundary conditions
  use global_vars, only: imx
  use global_vars, only: jmx
  use global_vars, only: kmx
  use global_vars, only: n_var
  use global_vars, only: qp
  use global_vars, only: imin_id
  use global_vars, only: jmin_id
  use global_vars, only: kmin_id
  use global_vars, only: imax_id
  use global_vars, only: jmax_id
  use global_vars, only: kmax_id
  use global_vars, only: layers
  use global_vars, only: process_id
  use global_vars, only : dir_switch
  use global_vars, only: PbcId

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
  use mapping, only : read_interface_map

  use utils      , only:   alloc
  use utils      , only: dealloc
  use utils      , only: dmsg

#include "error.inc"
#include "mpi.inc"
  private
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

  public :: setup_interface
  public :: destroy_interface
  public :: apply_interface

  contains

    subroutine setup_interface()
      !< Allocate memory for the data communication between processors
      implicit none
      character(len=*), parameter :: &
        errmsg="module: interface, subrouinte setup"
      ibuf_size = (jmx-1)*(kmx-1)*n_var*layers
      jbuf_size = (imx-1)*(kmx-1)*n_var*layers
      kbuf_size = (imx-1)*(jmx-1)*n_var*layers
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
    end subroutine setup_interface


    subroutine destroy_interface()
      !< Deallocate all the memory being used  for data communication between processors
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
    end subroutine destroy_interface


    subroutine apply_interface()
      !< MPISEND_RECV call to exchange interface infromation between
      !< connected blocks.
      implicit none
      integer:: i,j,k,n,l
      integer:: status(MPI_STATUS_SIZE)
      integer:: ierr
      integer:: tag=1
      integer:: count=0


      !----------------------------------------------------------
      ! call pattern is change for first block = 0
      ! to avoid O-Grid infinite loop for mpi communication call
      !-----------------------------------------------------------
      if(mod(process_id,2)==0)then
      !--- IMAX ---!
      if(imax_id>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do j=1,jmx-1
                count=count+1
                imax_send_buf(count) = qp(imx-l,j,k,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(imax_send_buf,ibuf_size, MPI_DOUBLE_PRECISION, imax_id,tag,&
                          imax_recv_buf,ibuf_size, MPI_DOUBLE_PRECISION, imax_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
!        if(mpi_class(2)==0)then
!          call MPI_SEND(imax_send_buf,ibuf_size,MPI_DOUBLE_PRECISION,imax_id,tag,MPI_COMM_WORLD,ierr)
!          call MPI_RECV(imax_recv_buf,ibuf_size,MPI_DOUBLE_PRECISION,imax_id,tag,MPI_COMM_WORLD,status, ierr)
!        elseif(mpi_class(2)==1)then
!          call MPI_RECV(imax_recv_buf,ibuf_size,MPI_DOUBLE_PRECISION,imax_id,tag,MPI_COMM_WORLD,status, ierr)
!          call MPI_SEND(imax_send_buf,ibuf_size,MPI_DOUBLE_PRECISION,imax_id,tag,MPI_COMM_WORLD,ierr)
!        else
!          Fatal_error
!        end if
        ! redistribute data
        if(dir_switch(2)==0)then
          count=0
          do n=1,n_var
            do l=1,layers
              do k=Pklo(2),Pkhi(2),PkDir(2)
                do j=Pjlo(2),Pjhi(2),PjDir(2)
                  count=count+1
                   qp(imx+l-1,j,k,n) = imax_recv_buf(count)
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
                   qp(imx+l-1,j,k,n) = imax_recv_buf(count)
                end do
              end do
            end do
          end do
        end if
      end if

      !--- IMIN ---!
      call dmsg(1, 'interface', 'apply_interface')
      if(imin_id>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do j=1,jmx-1
                count=count+1
                imin_send_buf(count) = qp(l,j,k,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(imin_send_buf,ibuf_size, MPI_DOUBLE_PRECISION, imin_id,tag,&
                          imin_recv_buf,ibuf_size, MPI_DOUBLE_PRECISION, imin_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
!        if(mpi_class(1)==0)then
!          call MPI_SEND(imin_send_buf,ibuf_size,MPI_DOUBLE_PRECISION,imin_id,tag,MPI_COMM_WORLD,ierr)
!          call MPI_RECV(imin_recv_buf,ibuf_size,MPI_DOUBLE_PRECISION,imin_id,tag,MPI_COMM_WORLD,status, ierr)
!        elseif(mpi_class(1)==1)then
!          call MPI_RECV(imin_recv_buf,ibuf_size,MPI_DOUBLE_PRECISION,imin_id,tag,MPI_COMM_WORLD,status, ierr)
!          call MPI_SEND(imin_send_buf,ibuf_size,MPI_DOUBLE_PRECISION,imin_id,tag,MPI_COMM_WORLD,ierr)
!        else
!          Fatal_error
!        end if
        ! redistribute data
        if(dir_switch(1)==0)then
          count=0
          do n=1,n_var
            do l=1,layers
              do k=Pklo(1),Pkhi(1),PkDir(1)
                do j=Pjlo(1),Pjhi(1),PjDir(1)
                  count=count+1
                  qp(1-l,j,k,n) = imin_recv_buf(count)
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
                  qp(1-l,j,k,n) = imin_recv_buf(count)
                end do
              end do
            end do
          end do
        end if
      end if

      else
      !--- IMIN ---!
      call dmsg(1, 'interface', 'apply_interface')
      if(imin_id>=0)then
        !collect data
        count=0
        do n=1,n_var
          do l=1,layers
            do k=1,kmx-1
              do j=1,jmx-1
                count=count+1
                imin_send_buf(count) = qp(l,j,k,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(imin_send_buf,ibuf_size, MPI_DOUBLE_PRECISION, imin_id,tag,&
                          imin_recv_buf,ibuf_size, MPI_DOUBLE_PRECISION, imin_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
!        if(mpi_class(1)==0)then
!          call MPI_SEND(imin_send_buf,ibuf_size,MPI_DOUBLE_PRECISION,imin_id,tag,MPI_COMM_WORLD,ierr)
!          call MPI_RECV(imin_recv_buf,ibuf_size,MPI_DOUBLE_PRECISION,imin_id,tag,MPI_COMM_WORLD,status, ierr)
!        elseif(mpi_class(1)==1)then
!          call MPI_RECV(imin_recv_buf,ibuf_size,MPI_DOUBLE_PRECISION,imin_id,tag,MPI_COMM_WORLD,status, ierr)
!          call MPI_SEND(imin_send_buf,ibuf_size,MPI_DOUBLE_PRECISION,imin_id,tag,MPI_COMM_WORLD,ierr)
!        else
!          Fatal_error
!        end if
        ! redistribute data
        if(dir_switch(1)==0)then
          count=0
          do n=1,n_var
            do l=1,layers
              do k=Pklo(1),Pkhi(1),PkDir(1)
                do j=Pjlo(1),Pjhi(1),PjDir(1)
                  count=count+1
                  qp(1-l,j,k,n) = imin_recv_buf(count)
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
                  qp(1-l,j,k,n) = imin_recv_buf(count)
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
                imax_send_buf(count) = qp(imx-l,j,k,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(imax_send_buf,ibuf_size, MPI_DOUBLE_PRECISION, imax_id,tag,&
                          imax_recv_buf,ibuf_size, MPI_DOUBLE_PRECISION, imax_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
!        if(mpi_class(2)==0)then
!          call MPI_SEND(imax_send_buf,ibuf_size,MPI_DOUBLE_PRECISION,imax_id,tag,MPI_COMM_WORLD,ierr)
!          call MPI_RECV(imax_recv_buf,ibuf_size,MPI_DOUBLE_PRECISION,imax_id,tag,MPI_COMM_WORLD,status, ierr)
!        elseif(mpi_class(2)==1)then
!          call MPI_RECV(imax_recv_buf,ibuf_size,MPI_DOUBLE_PRECISION,imax_id,tag,MPI_COMM_WORLD,status, ierr)
!          call MPI_SEND(imax_send_buf,ibuf_size,MPI_DOUBLE_PRECISION,imax_id,tag,MPI_COMM_WORLD,ierr)
!        else
!          Fatal_error
!        end if
        ! redistribute data
        if(dir_switch(2)==0)then
          count=0
          do n=1,n_var
            do l=1,layers
              do k=Pklo(2),Pkhi(2),PkDir(2)
                do j=Pjlo(2),Pjhi(2),PjDir(2)
                  count=count+1
                   qp(imx+l-1,j,k,n) = imax_recv_buf(count)
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
                   qp(imx+l-1,j,k,n) = imax_recv_buf(count)
                end do
              end do
            end do
          end do
        end if
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
                jmin_send_buf(count) = qp(i,l,k,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(jmin_send_buf,jbuf_size, MPI_DOUBLE_PRECISION, jmin_id,tag,&
                          jmin_recv_buf,jbuf_size, MPI_DOUBLE_PRECISION, jmin_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
!        if(mpi_class(3)==0)then
!          call MPI_SEND(jmin_send_buf,jbuf_size,MPI_DOUBLE_PRECISION,jmin_id,tag,MPI_COMM_WORLD,ierr)
!          call MPI_RECV(jmin_recv_buf,jbuf_size,MPI_DOUBLE_PRECISION,jmin_id,tag,MPI_COMM_WORLD,status, ierr)
!        elseif(mpi_class(3)==1)then
!          call MPI_RECV(jmin_recv_buf,jbuf_size,MPI_DOUBLE_PRECISION,jmin_id,tag,MPI_COMM_WORLD,status, ierr)
!          call MPI_SEND(jmin_send_buf,jbuf_size,MPI_DOUBLE_PRECISION,jmin_id,tag,MPI_COMM_WORLD,ierr)
!        else
!          Fatal_error
!        end if
        ! redistribute data
        if(dir_switch(3)==0)then
          count=0
          do n=1,n_var
            do l=1,layers
              do k=Pklo(3),Pkhi(3),PkDir(3)
                do i=Pilo(3),Pihi(3),PiDir(3)
                  count=count+1
                  qp(i,1-l,k,n) = jmin_recv_buf(count)
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
                  qp(i,1-l,k,n) = jmin_recv_buf(count)
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
                jmax_send_buf(count) = qp(i,jmx-l,k,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(jmax_send_buf,jbuf_size, MPI_DOUBLE_PRECISION, jmax_id,tag,&
                          jmax_recv_buf,jbuf_size, MPI_DOUBLE_PRECISION, jmax_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
!        if(mpi_class(4)==0)then
!          call MPI_SEND(jmax_send_buf,jbuf_size,MPI_DOUBLE_PRECISION,jmax_id,tag,MPI_COMM_WORLD,ierr)
!          call MPI_RECV(jmax_recv_buf,jbuf_size,MPI_DOUBLE_PRECISION,jmax_id,tag,MPI_COMM_WORLD,status, ierr)
!        elseif(mpi_class(4)==1)then
!          call MPI_RECV(jmax_recv_buf,jbuf_size,MPI_DOUBLE_PRECISION,jmax_id,tag,MPI_COMM_WORLD,status, ierr)
!          call MPI_SEND(jmax_send_buf,jbuf_size,MPI_DOUBLE_PRECISION,jmax_id,tag,MPI_COMM_WORLD,ierr)
!        else
!          Fatal_error
!        end if

        ! redistribute data
        if(dir_switch(4)==0)then
          count=0
          do n=1,n_var
            do l=1,layers
              do k=Pklo(4),Pkhi(4),PkDir(4)
                do i=Pilo(4),Pihi(4),PiDir(4)
                  count=count+1
                  qp(i,jmx+l-1,k,n) = jmax_recv_buf(count)
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
                  qp(i,jmx+l-1,k,n) = jmax_recv_buf(count)
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
                kmin_send_buf(count) = qp(i,j,l,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(kmin_send_buf,kbuf_size, MPI_DOUBLE_PRECISION, kmin_id,tag,&
                          kmin_recv_buf,kbuf_size, MPI_DOUBLE_PRECISION, kmin_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
!        if(mpi_class(5)==0)then
!          call MPI_SEND(kmin_send_buf,kbuf_size,MPI_DOUBLE_PRECISION,kmin_id,tag,MPI_COMM_WORLD,ierr)
!          call MPI_RECV(kmin_recv_buf,kbuf_size,MPI_DOUBLE_PRECISION,kmin_id,tag,MPI_COMM_WORLD,status, ierr)
!        elseif(mpi_class(5)==1)then
!          call MPI_RECV(kmin_recv_buf,kbuf_size,MPI_DOUBLE_PRECISION,kmin_id,tag,MPI_COMM_WORLD,status, ierr)
!          call MPI_SEND(kmin_send_buf,kbuf_size,MPI_DOUBLE_PRECISION,kmin_id,tag,MPI_COMM_WORLD,ierr)
!        else
!          Fatal_error
!        end if
        ! redistribute data
        if(dir_switch(5)==0)then
          count=0
          do n=1,n_var
            do l=1,layers
              do j=Pjlo(5),Pjhi(5),PjDir(5)
                do i=Pilo(5),Pihi(5),PiDir(5)
                  count=count+1
                  qp(i,j,1-l,n) = kmin_recv_buf(count)
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
                  qp(i,j,1-l,n) = kmin_recv_buf(count)
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
                kmax_send_buf(count) = qp(i,j,kmx-l,n)
              end do
            end do
          end do
        end do
        call MPI_SENDRECV(kmax_send_buf,kbuf_size, MPI_DOUBLE_PRECISION, kmax_id,tag,&
                          kmax_recv_buf,kbuf_size, MPI_DOUBLE_PRECISION, kmax_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
!        if(mpi_class(6)==0)then
!          call MPI_SEND(kmax_send_buf,kbuf_size,MPI_DOUBLE_PRECISION,kmax_id,tag,MPI_COMM_WORLD,ierr)
!          call MPI_RECV(kmax_recv_buf,kbuf_size,MPI_DOUBLE_PRECISION,kmax_id,tag,MPI_COMM_WORLD,status, ierr)
!        elseif(mpi_class(6)==1)then
!          call MPI_RECV(kmax_recv_buf,kbuf_size,MPI_DOUBLE_PRECISION,kmax_id,tag,MPI_COMM_WORLD,status, ierr)
!          call MPI_SEND(kmax_send_buf,kbuf_size,MPI_DOUBLE_PRECISION,kmax_id,tag,MPI_COMM_WORLD,ierr)
!        else
!          Fatal_error
!        end if
        ! redistribute data
        if(dir_switch(6)==0)then
          count=0
          do n=1,n_var
            do l=1,layers
              do j=Pjlo(6),Pjhi(6),PjDir(6)
                do i=Pilo(6),Pihi(6),PiDir(6)
                  count=count+1
                  qp(i,j,kmx+l-1,n) = kmax_recv_buf(count)
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
                  qp(i,j,kmx+l-1,n) = kmax_recv_buf(count)
                end do
              end do
            end do
          end do
        end if
      end if
      call apply_periodic_bc()
    end subroutine apply_interface


    subroutine apply_periodic_bc()
      !<If a block is connected to another block in perodic
      !<fashion, this subroutine will take care of that boundary condition.
      implicit none
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
                imin_send_buf(count) = qp(l,j,k,n)
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
                qp(1-l,j,k,n) = imin_recv_buf(count)
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
                imax_send_buf(count) = qp(imx-l,j,k,n)
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
                 qp(imx+l-1,j,k,n) = imax_recv_buf(count)
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
                jmin_send_buf(count) = qp(i,l,k,n)
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
                qp(i,1-l,k,n) = jmin_recv_buf(count)
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
                jmax_send_buf(count) = qp(i,jmx-l,k,n)
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
                qp(i,jmx+l-1,k,n) = jmax_recv_buf(count)
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
                kmin_send_buf(count) = qp(i,j,l,n)
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
                qp(i,j,1-l,n) = kmin_recv_buf(count)
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
                kmax_send_buf(count) = qp(i,j,kmx-l,n)
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
                qp(i,j,kmx+l-1,n) = kmax_recv_buf(count)
              end do
            end do
          end do
        end do
      end if


    end subroutine apply_periodic_bc

end module interface1
