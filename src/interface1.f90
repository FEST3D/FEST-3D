!< This module handles the MPI Communication calls for interface boundary conditions
module interface1
  !< This module handles the MPI Communication calls for interface boundary conditions
  use vartypes
  use mpi
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

#include "debug.h"
#include "error.h"
  private
  integer :: layers = 3
  integer :: ibuf_size
  !< Size of the buffer for I face interface
  integer :: jbuf_size
  !< Size of the buffer for J face interface
  integer :: kbuf_size
  !< Size of the buffer for K face interface
  real(wp), dimension(:), allocatable :: imin_send_buf
  !< Array to store data to send data for Imin face
  real(wp), dimension(:), allocatable :: jmin_send_buf
  !< Array to store data to send data for Jmin face
  real(wp), dimension(:), allocatable :: kmin_send_buf
  !< Array to store data to send data for Kmin face
  real(wp), dimension(:), allocatable :: imin_recv_buf
  !< Array to store data to receive data for Imin face
  real(wp), dimension(:), allocatable :: jmin_recv_buf
  !< Array to store data to receive data for Jmin face
  real(wp), dimension(:), allocatable :: kmin_recv_buf
  !< Array to store data to receive data for Kmin face
  real(wp), dimension(:), allocatable :: imax_send_buf
  !< Array to store data to send data for Imax face
  real(wp), dimension(:), allocatable :: jmax_send_buf
  !< Array to store data to send data for Jmax face
  real(wp), dimension(:), allocatable :: kmax_send_buf
  !< Array to store data to send data for Kmax face
  real(wp), dimension(:), allocatable :: imax_recv_buf
  !< Array to store data to receive data for Imax face
  real(wp), dimension(:), allocatable :: jmax_recv_buf
  !< Array to store data to receive data for Jmax face
  real(wp), dimension(:), allocatable :: kmax_recv_buf
  !< Array to store data to receive data for Kmax face

  public :: setup_interface
  public :: apply_interface

  contains

    subroutine setup_interface(control, dims)
      !< Allocate memory for the data communication between processors
      implicit none
      type(controltype), intent(in) :: control
      !< Control parameters: n_var
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      integer :: imx, jmx, kmx, n_var
      character(len=*), parameter ::errmsg="module: interface, subrouinte setup"
      !< Error message

      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx

      n_var = control%n_var

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



    subroutine apply_interface(qp, control, bc, dims)
      !< MPISEND_RECV call to exchange interface infromation between
      !< connected blocks.
      implicit none
      type(controltype), intent(in) :: control
      !< Control parameters
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2,1:dims%n_var) :: qp
      !< Store primitive variable at cell center
      type(boundarytype), intent(in) :: bc
      !< boundary conditions and fixed values
      integer:: i,j,k,n,l
      integer:: status(MPI_STATUS_SIZE)
      integer:: ierr
      integer:: tag=1
      integer:: count=0
      integer :: imx, jmx, kmx, n_var

      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx

      n_var = control%n_var

      !----------------------------------------------------------
      ! call pattern is change for first block = 0
      ! to avoid O-Grid infinite loop for mpi communication call
      !-----------------------------------------------------------
      if(mod(control%process_id,2)==0)then
      !--- IMAX ---!
      if(bc%imax_id>=0)then
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
        call MPI_SENDRECV(imax_send_buf,ibuf_size, MPI_DOUBLE_PRECISION, bc%imax_id,tag,&
                          imax_recv_buf,ibuf_size, MPI_DOUBLE_PRECISION, bc%imax_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        if(bc%dir_switch(2)==0)then
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
      DebugCall('apply_interface')
      if(bc%imin_id>=0)then
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
        call MPI_SENDRECV(imin_send_buf,ibuf_size, MPI_DOUBLE_PRECISION, bc%imin_id,tag,&
                          imin_recv_buf,ibuf_size, MPI_DOUBLE_PRECISION, bc%imin_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        if(bc%dir_switch(1)==0)then
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
      DebugCall('apply_interface')
      if(bc%imin_id>=0)then
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
        call MPI_SENDRECV(imin_send_buf,ibuf_size, MPI_DOUBLE_PRECISION, bc%imin_id,tag,&
                          imin_recv_buf,ibuf_size, MPI_DOUBLE_PRECISION, bc%imin_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        if(bc%dir_switch(1)==0)then
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
      if(bc%imax_id>=0)then
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
        call MPI_SENDRECV(imax_send_buf,ibuf_size, MPI_DOUBLE_PRECISION, bc%imax_id,tag,&
                          imax_recv_buf,ibuf_size, MPI_DOUBLE_PRECISION, bc%imax_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        if(bc%dir_switch(2)==0)then
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
      if(bc%jmin_id>=0)then
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
        call MPI_SENDRECV(jmin_send_buf,jbuf_size, MPI_DOUBLE_PRECISION, bc%jmin_id,tag,&
                          jmin_recv_buf,jbuf_size, MPI_DOUBLE_PRECISION, bc%jmin_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        if(bc%dir_switch(3)==0)then
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
      if(bc%jmax_id>=0)then
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
        call MPI_SENDRECV(jmax_send_buf,jbuf_size, MPI_DOUBLE_PRECISION, bc%jmax_id,tag,&
                          jmax_recv_buf,jbuf_size, MPI_DOUBLE_PRECISION, bc%jmax_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        if(bc%dir_switch(4)==0)then
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
      if(bc%kmin_id>=0)then
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
        call MPI_SENDRECV(kmin_send_buf,kbuf_size, MPI_DOUBLE_PRECISION, bc%kmin_id,tag,&
                          kmin_recv_buf,kbuf_size, MPI_DOUBLE_PRECISION, bc%kmin_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        if(bc%dir_switch(5)==0)then
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
      if(bc%kmax_id>=0)then
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
        call MPI_SENDRECV(kmax_send_buf,kbuf_size, MPI_DOUBLE_PRECISION, bc%kmax_id,tag,&
                          kmax_recv_buf,kbuf_size, MPI_DOUBLE_PRECISION, bc%kmax_id,tag,&
                          MPI_COMM_WORLD,status,ierr)
        ! redistribute data
        if(bc%dir_switch(6)==0)then
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
      call apply_periodic_bc(qp, control, bc, dims)
    end subroutine apply_interface


    subroutine apply_periodic_bc(qp, control, bc, dims)
      !<If a block is connected to another block in perodic
      !<fashion, this subroutine will take care of that boundary condition.
      implicit none
      type(controltype), intent(in) :: control
      !< Control parameters
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      type(boundarytype), intent(in) :: bc
      !< boundary conditions and fixed values
      real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2,1:dims%n_var) :: qp
      !< Store primitive variable at cell center
      integer:: i,j,k,n,l
      integer:: status(MPI_STATUS_SIZE)
      integer:: ierr
      integer:: tag=1
      integer:: count=0
      integer :: imx, jmx, kmx, n_var

      DebugCall('apply_periodic_boundary_condition')

      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx

      n_var = control%n_var

      if(bc%PbcId(1)>=0)then
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
        call MPI_SENDRECV(imin_send_buf,ibuf_size, MPI_DOUBLE_PRECISION, bc%PbcId(1),tag,&
                          imin_recv_buf,ibuf_size, MPI_DOUBLE_PRECISION, bc%PbcId(1),tag,&
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

      if(bc%PbcId(2)>=0)then
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
        call MPI_SENDRECV(imax_send_buf,ibuf_size, MPI_DOUBLE_PRECISION, bc%PbcId(2),tag,&
                          imax_recv_buf,ibuf_size, MPI_DOUBLE_PRECISION, bc%PbcId(2),tag,&
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
      if(bc%PbcId(3)>=0)then
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
        call MPI_SENDRECV(jmin_send_buf,jbuf_size, MPI_DOUBLE_PRECISION, bc%PbcId(3),tag,&
                          jmin_recv_buf,jbuf_size, MPI_DOUBLE_PRECISION, bc%PbcId(3),tag,&
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
      if(bc%PbcId(4)>=0)then
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
        call MPI_SENDRECV(jmax_send_buf,jbuf_size, MPI_DOUBLE_PRECISION, bc%PbcId(4),tag,&
                          jmax_recv_buf,jbuf_size, MPI_DOUBLE_PRECISION, bc%PbcId(4),tag,&
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
      if(bc%PbcId(5)>=0)then
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
        call MPI_SENDRECV(kmin_send_buf,kbuf_size, MPI_DOUBLE_PRECISION, bc%PbcId(5),tag,&
                          kmin_recv_buf,kbuf_size, MPI_DOUBLE_PRECISION, bc%PbcId(5),tag,&
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
      if(bc%PbcId(6)>=0)then
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
        call MPI_SENDRECV(kmax_send_buf,kbuf_size, MPI_DOUBLE_PRECISION, bc%PbcId(6),tag,&
                          kmax_recv_buf,kbuf_size, MPI_DOUBLE_PRECISION, bc%PbcId(6),tag,&
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
