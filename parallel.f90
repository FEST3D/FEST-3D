module parallel
  ! contains routines to execute parallelly
  use mpi
  use utils, only: alloc, dealloc, dmsg, DEBUG_LEVEL
  use grid, only: imx,jmx,kmx
  use layout, only: process_id,imin_id,imax_id,jmin_id,jmax_id,&
       kmin_id,kmax_id
  use state 
  implicit none

  real, public, dimension(:),allocatable :: imin_send_buf,imin_recv_buf,&
       imax_send_buf,imax_recv_buf,jmin_send_buf,jmin_recv_buf, &
       jmax_send_buf,jmax_recv_buf,kmin_send_buf,kmin_recv_buf,&
       kmax_send_buf,kmax_recv_buf
  public :: allocate_buffer_cells
  public :: send_recv_i
  public :: send_recv_j
  public :: send_recv_k
  public :: send_recv

contains

  subroutine send_recv(layers)
    integer, intent(in) :: layers
    call dmsg(1, 'parallel', 'send_recv')
    call send_recv_i(layers)
    call send_recv_j(layers)
    call send_recv_k(layers)

  end subroutine send_recv

  subroutine send_recv_i(layers)
    !-----------------------------------------------------------
    !send and receive data left and right 
    ! i.e along i direction 
    !-----------------------------------------------------------
    implicit none
    integer, intent(in) :: layers
    integer :: j,k,l,count =1
    integer :: ierr,buf
    integer :: status(MPI_STATUS_SIZE)
    call dmsg(1, 'parallel', 'send_recv_i')
    if(imin_id >= 0 ) then
       !print *, "left id is", process_id,left_id
       ! first send message
       ! left_send_buf(1) = 10;
       count = 1
       do k = 0, kmx
          do j = 0, jmx

             do l = 1, layers
                imin_send_buf(count) = density(l,j,k)
                count = count+1
             end do
             !             imin_send_buf(count) = density(1,j,k)
             !             count = count+1
             !             imin_send_buf(count) = density(2,j,k)
             !             count = count+1
             !             imin_send_buf(count) = density(3,j,k)
             !             count = count+1             
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                imin_send_buf(count) = pressure(l,j,k)
                count = count+1
             end do
             !             imin_send_buf(count) = pressure(1,j,k)
             !             count = count+1
             !             imin_send_buf(count) = pressure(2,j,k)
             !             count = count+1
             !             imin_send_buf(count) = pressure(3,j,k)
             !             count = count+1 
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                imin_send_buf(count) = x_speed(l,j,k)
                count = count+1
             end do
             !             imin_send_buf(count) = x_speed(1,j,k)
             !             count = count+1
             !             imin_send_buf(count) = x_speed(2,j,k)
             !             count = count+1
             !             imin_send_buf(count) = x_speed(3,j,k)
             !             count = count+1 
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                imin_send_buf(count) = y_speed(l,j,k)
                count = count+1
             end do

             !             imin_send_buf(count) = y_speed(1,j,k)
             !             count = count+1
             !             imin_send_buf(count) = y_speed(2,j,k)
             !             count = count+1
             !             imin_send_buf(count) = y_speed(3,j,k)
             !             count = count+1 
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                imin_send_buf(count) = z_speed(l,j,k)
                count = count+1
             end do
             !             imin_send_buf(count) = z_speed(1,j,k)
             !             count = count+1
             !             imin_send_buf(count) = z_speed(2,j,k)
             !             count = count+1
             !             imin_send_buf(count) = z_speed(3,j,k)
             !             count = count+1 
          end do
       end do
       !print *,"count is ", count 
       ! send message to left process
       buf = (jmx+1)*(kmx+1)*5*layers ! three cells       
       !do k = 1, buf
       !print *,'left send - ', process_id, k ,left_send_buf(k)
       !end do       
       call MPI_SEND(imin_send_buf,buf,MPI_DOUBLE_PRECISION,imin_id,1,MPI_COMM_WORLD, ierr)
       call MPI_RECV(imin_recv_buf,buf,MPI_DOUBLE_PRECISION,imin_id,1,MPI_COMM_WORLD,status,ierr)
       ! updating solution
       !do k = 1, buf
       !print *,'left recv - ', process_id, k ,left_recv_buf(k)
       !end do       
       count = 1
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                density(1-l,j,k) = imin_recv_buf(count)
                count = count+1
             end do
             !              density(0,j,k) = imin_recv_buf(count)               
             !              count = count+1
             !              density(-1,j,k) = imin_recv_buf(count)               
             !              count = count+1
             !              density(-2,j,k) = imin_recv_buf(count)               
             !              count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                pressure(1-l,j,k) = imin_recv_buf(count)
                count = count+1
             end do

             !             pressure(0,j,k) = imin_recv_buf(count)
             !             count = count+1
             !             pressure(-1,j,k) = imin_recv_buf(count)
             !             count = count+1
             !             pressure(-2,j,k) = imin_recv_buf(count)
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                x_speed(1-l,j,k) = imin_recv_buf(count)
                count = count+1
             end do
             !             x_speed(0,j,k) = imin_recv_buf(count)             
             !             count = count+1
             !             x_speed(-1,j,k) = imin_recv_buf(count)             
             !             count = count+1
             !             x_speed(-2,j,k) = imin_recv_buf(count)             
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                y_speed(1-l,j,k) = imin_recv_buf(count)
                count = count+1
             end do
             !             y_speed(0,j,k) = imin_recv_buf(count)
             !             count = count+1
             !             y_speed(-1,j,k) = imin_recv_buf(count)
             !             count = count+1
             !             y_speed(-2,j,k) = imin_recv_buf(count)
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                z_speed(1-l,j,k) = imin_recv_buf(count)
                count = count+1
             end do
             !             z_speed(0,j,k) = imin_recv_buf(count)
             !             count = count+1
             !             z_speed(-1,j,k) = imin_recv_buf(count)
             !             count = count+1
             !             z_speed(-2,j,k) = imin_recv_buf(count)
             !             count = count+1
          end do
       end do
    end if


    if(imax_id >= 0) then
       !print *,"right id is ", process_id,right_id
       buf = (jmx+1)*(kmx+1)*5*layers
       call MPI_RECV(imax_recv_buf,buf,MPI_DOUBLE_PRECISION,imax_id,1,MPI_COMM_WORLD,status,ierr)
       ! updating solution       
       count = 1
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                density(imx+l-1,j,k) = imax_recv_buf(count)
                count = count+1
             end do

             !             density(imx,j,k) = imax_recv_buf(count)
             !             count = count+1
             !             density(imx+1,j,k) = imax_recv_buf(count)
             !             count = count+1
             !             density(imx+2,j,k) = imax_recv_buf(count)
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                pressure(imx+l-1,j,k) = imax_recv_buf(count)
                count = count+1
             end do
             !             pressure(imx,j,k) = imax_recv_buf(count)
             !             count = count+1
             !             pressure(imx+1,j,k) = imax_recv_buf(count)
             !             count = count+1
             !             pressure(imx+2,j,k) = imax_recv_buf(count)
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                x_speed(imx+l-1,j,k) = imax_recv_buf(count)
                count = count+1
             end do
             !             x_speed(imx,j,k) = imax_recv_buf(count)
             !             count = count+1
             !             x_speed(imx+1,j,k) = imax_recv_buf(count)
             !             count = count+1
             !             x_speed(imx+2,j,k) = imax_recv_buf(count)
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                y_speed(imx+l-1,j,k) = imax_recv_buf(count)
                count = count+1
             end do
             !             y_speed(imx,j,k) = imax_recv_buf(count)
             !             count = count+1
             !             y_speed(imx+1,j,k) = imax_recv_buf(count)
             !             count = count+1
             !             y_speed(imx+2,j,k) = imax_recv_buf(count)
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                z_speed(imx+l-1,j,k) = imax_recv_buf(count)
                count = count+1
             end do
             !             z_speed(imx,j,k) = imax_recv_buf(count)
             !             count = count+1
             !             z_speed(imx+1,j,k) = imax_recv_buf(count)
             !             count = count+1
             !             z_speed(imx+2,j,k) = imax_recv_buf(count)
             !             count = count+1
          end do
       end do
       ! creating right send  buffer
       count = 1
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                imax_send_buf(count) = density(imx-l,j,k)
                count = count+1
             end do
             !             imax_send_buf(count) = density(imx-1,j,k)             
             !             count = count+1
             !             imax_send_buf(count) = density(imx-2,j,k)             
             !             count = count+1
             !             imax_send_buf(count) = density(imx-3,j,k)             
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                imax_send_buf(count) = pressure(imx-l,j,k)
                count = count+1
             end do
             !             imax_send_buf(count) = pressure(imx-1,j,k)
             !             count = count+1
             !             imax_send_buf(count) = pressure(imx-2,j,k)
             !             count = count+1
             !             imax_send_buf(count) = pressure(imx-3,j,k)
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                imax_send_buf(count) = x_speed(imx-l,j,k)
                count = count+1
             end do
             !             imax_send_buf(count) = x_speed(imx-1,j,k)             
             !             count = count+1
             !             imax_send_buf(count) = x_speed(imx-2,j,k)             
             !             count = count+1
             !             imax_send_buf(count) = x_speed(imx-3,j,k)             
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                imax_send_buf(count) = y_speed(imx-l,j,k)
                count = count+1
             end do
             !             imax_send_buf(count) = y_speed(imx-1,j,k)
             !             count = count+1
             !             imax_send_buf(count) = y_speed(imx-2,j,k)
             !             count = count+1
             !             imax_send_buf(count) = y_speed(imx-3,j,k)
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             do l = 1, layers
                imax_send_buf(count) = z_speed(imx-l,j,k)
                count = count+1
             end do
             !             imax_send_buf(count) = z_speed(imx-1,j,k)
             !             count = count+1
             !             imax_send_buf(count) = z_speed(imx-2,j,k)
             !             count = count+1
             !             imax_send_buf(count) = z_speed(imx-3,j,k)
             !             count = count+1
          end do
       end do
       !do k = 1, buf
       !print *,'right send - ', process_id, k ,right_send_buf(k)
       !end do
       !print *, "total size is", buf 
       call MPI_SEND(imax_send_buf,buf,MPI_DOUBLE_PRECISION,imax_id,1,MPI_COMM_WORLD, ierr)    
    end if
  end subroutine send_recv_i


  subroutine send_recv_j(layers)
    !-----------------------------------------------------------
    !send and receive data left and right 
    ! i.e along i direction 
    !-----------------------------------------------------------
    implicit none
    integer, intent(in) :: layers
    integer :: i,k,l,count =1
    integer :: ierr,buf
    integer :: status(MPI_STATUS_SIZE)
    call dmsg(1, 'parallel', 'send_recv_i')
    if(jmin_id >= 0 ) then
       !print *, "left id is", process_id,left_id
       ! first send message
       ! left_send_buf(1) = 10;
       count = 1
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                jmin_send_buf(count) = density(i,l,k)
                count = count+1
             end do
             !             jmin_send_buf(count) = density(i,1,k)
             !             count = count+1
             !             jmin_send_buf(count) = density(i,2,k)
             !             count = count+1
             !             jmin_send_buf(count) = density(i,3,k)
             !             count = count+1             
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                jmin_send_buf(count) = pressure(i,l,k)
                count = count+1
             end do
             !             jmin_send_buf(count) = pressure(i,1,k)
             !             count = count+1
             !             jmin_send_buf(count) = pressure(i,2,k)
             !             count = count+1
             !             jmin_send_buf(count) = pressure(i,3,k)
             !             count = count+1 
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                jmin_send_buf(count) = x_speed(i,l,k)
                count = count+1
             end do
             !             jmin_send_buf(count) = x_speed(i,1,k)
             !             count = count+1
             !             jmin_send_buf(count) = x_speed(i,2,k)
             !             count = count+1
             !             jmin_send_buf(count) = x_speed(i,3,k)
             !             count = count+1 
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                jmin_send_buf(count) = y_speed(i,l,k)
                count = count+1
             end do
             !             jmin_send_buf(count) = y_speed(i,1,k)
             !             count = count+1
             !             jmin_send_buf(count) = y_speed(i,2,k)
             !             count = count+1
             !             jmin_send_buf(count) = y_speed(i,3,k)
             !             count = count+1 
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                jmin_send_buf(count) = z_speed(i,l,k)
                count = count+1
             end do
             !             jmin_send_buf(count) = z_speed(i,1,k)
             !             count = count+1
             !             jmin_send_buf(count) = z_speed(i,2,k)
             !             count = count+1
             !             jmin_send_buf(count) = z_speed(i,3,k)
             !             count = count+1 
          end do
       end do
       !print *,"count is ", count 
       ! send message to left process
       buf = (imx+1)*(kmx+1)*5*layers ! three layers       
       !do k = 1, buf
       !print *,'left send - ', process_id, k ,left_send_buf(k)
       !end do       
       call MPI_SEND(jmin_send_buf,buf,MPI_DOUBLE_PRECISION,jmin_id,1,MPI_COMM_WORLD, ierr)
       call MPI_RECV(jmin_recv_buf,buf,MPI_DOUBLE_PRECISION,jmin_id,1,MPI_COMM_WORLD,status,ierr)
       ! updating solution
       !do k = 1, buf
       !print *,'left recv - ', process_id, k ,left_recv_buf(k)
       !end do       
       count = 1
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                density(i,1-l,k) = jmin_recv_buf(count) 
                count = count+1
             end do
             !             density(i,0,k) = jmin_recv_buf(count)               
             !             count = count+1
             !             density(i,-1,k) = jmin_recv_buf(count)               
             !             count = count+1
             !             density(i,-2,k) = jmin_recv_buf(count)               
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                pressure(i,1-l,k) = jmin_recv_buf(count) 
                count = count+1
             end do
             !             pressure(i,0,k) = jmin_recv_buf(count)
             !             count = count+1
             !             pressure(i,-1,k) = jmin_recv_buf(count)
             !             count = count+1
             !             pressure(i,-2,k) = jmin_recv_buf(count)
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                x_speed(i,1-l,k) = jmin_recv_buf(count) 
                count = count+1
             end do
             !             x_speed(i,0,k) = jmin_recv_buf(count)             
             !             count = count+1
             !             x_speed(i,-1,k) = jmin_recv_buf(count)             
             !             count = count+1
             !             x_speed(i,-2,k) = jmin_recv_buf(count)             
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                y_speed(i,1-l,k) = jmin_recv_buf(count) 
                count = count+1
             end do
             !             y_speed(i,0,k) = jmin_recv_buf(count)             
             !             count = count+1
             !             y_speed(i,-1,k) = jmin_recv_buf(count)             
             !             count = count+1
             !             y_speed(i,-2,k) = jmin_recv_buf(count)             
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                z_speed(i,1-l,k) = jmin_recv_buf(count) 
                count = count+1
             end do
             !             z_speed(i,0,k) = jmin_recv_buf(count)             
             !             count = count+1
             !             z_speed(i,-1,k) = jmin_recv_buf(count)             
             !             count = count+1
             !             z_speed(i,-2,k) = jmin_recv_buf(count)             
             !             count = count+1
          end do
       end do
    end if


    if(jmax_id >= 0) then
       !print *,"right id is ", process_id,right_id
       buf = (imx+1)*(kmx+1)*5*layers
       call MPI_RECV(jmax_recv_buf,buf,MPI_DOUBLE_PRECISION,jmax_id,1,MPI_COMM_WORLD,status,ierr)
       ! updating solution       
       count = 1
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                density(i,jmx+l-1,k) = jmax_recv_buf(count) 
                count = count+1
             end do
             !             density(i,jmx,k) = jmax_recv_buf(count)
             !             count = count+1
             !             density(i,jmx+1,k) = jmax_recv_buf(count)
             !             count = count+1
             !             density(i,jmx+2,k) = jmax_recv_buf(count)
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                pressure(i,jmx+l-1,k) = jmax_recv_buf(count) 
                count = count+1
             end do
             !             pressure(i,jmx,k) = jmax_recv_buf(count)
             !             count = count+1
             !             pressure(i,jmx+1,k) = jmax_recv_buf(count)
             !             count = count+1
             !             pressure(i,jmx+2,k) = jmax_recv_buf(count)
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                x_speed(i,jmx+l-1,k) = jmax_recv_buf(count) 
                count = count+1
             end do
             !             x_speed(i,jmx,k) = jmax_recv_buf(count)
             !             count = count+1
             !             x_speed(i,jmx+1,k) = jmax_recv_buf(count)
             !             count = count+1
             !             x_speed(i,jmx+2,k) = jmax_recv_buf(count)
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                y_speed(i,jmx+l-1,k) = jmax_recv_buf(count) 
                count = count+1
             end do
             !             y_speed(i,jmx,k) = jmax_recv_buf(count)
             !             count = count+1
             !             y_speed(i,jmx+1,k) = jmax_recv_buf(count)
             !             count = count+1
             !             y_speed(i,jmx+2,k) = jmax_recv_buf(count)
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                z_speed(i,jmx+l-1,k) = jmax_recv_buf(count) 
                count = count+1
             end do
             !             z_speed(i,jmx,k) = jmax_recv_buf(count)
             !             count = count+1
             !             z_speed(i,jmx+1,k) = jmax_recv_buf(count)
             !             count = count+1
             !             z_speed(i,jmx+2,k) = jmax_recv_buf(count)
             !             count = count+1
          end do
       end do
       ! creating right send  buffer
       count = 1
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                jmax_send_buf(count) = density(i,jmx-l,k) 
                count = count+1
             end do
             !             jmax_send_buf(count) = density(i,jmx-1,k)             
             !             count = count+1
             !             jmax_send_buf(count) = density(i,jmx-2,k)             
             !             count = count+1
             !             jmax_send_buf(count) = density(i,jmx-3,k)             
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                jmax_send_buf(count) = pressure(i,jmx-l,k) 
                count = count+1
             end do
             !             jmax_send_buf(count) = pressure(i,jmx-1,k)             
             !             count = count+1
             !             jmax_send_buf(count) = pressure(i,jmx-2,k)             
             !             count = count+1
             !             jmax_send_buf(count) = pressure(i,jmx-3,k)             
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                jmax_send_buf(count) = x_speed(i,jmx-l,k) 
                count = count+1
             end do
             !             jmax_send_buf(count) = x_speed(i,jmx-1,k)             
             !             count = count+1
             !             jmax_send_buf(count) = x_speed(i,jmx-2,k)             
             !             count = count+1
             !             jmax_send_buf(count) = x_speed(i,jmx-3,k)             
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                jmax_send_buf(count) = y_speed(i,jmx-l,k) 
                count = count+1
             end do
             !             jmax_send_buf(count) = y_speed(i,jmx-1,k)             
             !             count = count+1
             !             jmax_send_buf(count) = y_speed(i,jmx-2,k)             
             !             count = count+1
             !             jmax_send_buf(count) = y_speed(i,jmx-3,k)             
             !             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             do l = 1, layers
                jmax_send_buf(count) = z_speed(i,jmx-l,k) 
                count = count+1
             end do
             !             jmax_send_buf(count) = z_speed(i,jmx-1,k)             
             !             count = count+1
             !             jmax_send_buf(count) = z_speed(i,jmx-2,k)             
             !             count = count+1
             !             jmax_send_buf(count) = z_speed(i,jmx-3,k)             
             !             count = count+1
          end do
       end do
       !do k = 1, buf
       !print *,'right send - ', process_id, k ,right_send_buf(k)
       !end do
       !print *, "total size is", buf 
       call MPI_SEND(jmax_send_buf,buf,MPI_DOUBLE_PRECISION,jmax_id,1,MPI_COMM_WORLD, ierr)    
    end if
  end subroutine send_recv_j

  subroutine send_recv_k(layers)
    !-----------------------------------------------------------
    !send and receive data left and right 
    ! i.e along i direction 
    !-----------------------------------------------------------
    implicit none
    integer, intent(in) :: layers
    integer :: j,i,l,count =1
    integer :: ierr,buf
    integer :: status(MPI_STATUS_SIZE)
    call dmsg(1, 'parallel', 'send_recv_k')
    if(kmin_id >= 0 ) then
       !print *, "left id is", process_id,left_id
       ! first send message
       ! left_send_buf(1) = 10;
       count = 1
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                kmin_send_buf(count) = density(i,j,l) 
                count = count+1
             end do
             !             kmin_send_buf(count) = density(i,j,1)
             !             count = count+1
             !             kmin_send_buf(count) = density(i,j,2)
             !             count = count+1
             !             kmin_send_buf(count) = density(i,j,3)
             !             count = count+1             
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                kmin_send_buf(count) = pressure(i,j,l) 
                count = count+1
             end do
             !             kmin_send_buf(count) = pressure(i,j,1)
             !             count = count+1
             !             kmin_send_buf(count) = pressure(i,j,2)
             !             count = count+1
             !             kmin_send_buf(count) = pressure(i,j,3)
             !             count = count+1  
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                kmin_send_buf(count) = x_speed(i,j,l) 
                count = count+1
             end do
             !             kmin_send_buf(count) = x_speed(i,j,1)
             !             count = count+1
             !             kmin_send_buf(count) = x_speed(i,j,2)
             !             count = count+1
             !             kmin_send_buf(count) = x_speed(i,j,3)
             !             count = count+1  
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                kmin_send_buf(count) = y_speed(i,j,l) 
                count = count+1
             end do
             !             kmin_send_buf(count) = y_speed(i,j,1)
             !             count = count+1
             !             kmin_send_buf(count) = y_speed(i,j,2)
             !             count = count+1
             !             kmin_send_buf(count) = y_speed(i,j,3)
             !             count = count+1  
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                kmin_send_buf(count) = z_speed(i,j,l) 
                count = count+1
             end do
             !             kmin_send_buf(count) = z_speed(i,j,1)
             !             count = count+1
             !             kmin_send_buf(count) = z_speed(i,j,2)
             !             count = count+1
             !             kmin_send_buf(count) = z_speed(i,j,3)
             !             count = count+1  
          end do
       end do
       !print *,"count is ", count 
       ! send message to left process
       buf = (jmx+1)*(imx+1)*5*layers ! three cells       
       !do k = 1, buf
       !print *,'left send - ', process_id, k ,left_send_buf(k)
       !end do       
       call MPI_SEND(kmin_send_buf,buf,MPI_DOUBLE_PRECISION,kmin_id,1,MPI_COMM_WORLD, ierr)
       call MPI_RECV(kmin_recv_buf,buf,MPI_DOUBLE_PRECISION,kmin_id,1,MPI_COMM_WORLD,status,ierr)
       ! updating solution
       !do k = 1, buf
       !print *,'left recv - ', process_id, k ,left_recv_buf(k)
       !end do       
       count = 1
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                density(i,j,1-l) = kmin_recv_buf(count) 
                count = count+1
             end do
             !             density(i,j,0) = kmin_recv_buf(count)               
             !             count = count+1
             !             density(i,j,-1) = kmin_recv_buf(count)               
             !             count = count+1
             !             density(i,j,-2) = kmin_recv_buf(count)               
             !             count = count+1
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                pressure(i,j,1-l) = kmin_recv_buf(count) 
                count = count+1
             end do
             !             pressure(i,j,0) = kmin_recv_buf(count)               
             !             count = count+1
             !             pressure(i,j,-1) = kmin_recv_buf(count)               
             !             count = count+1
             !             pressure(i,j,-2) = kmin_recv_buf(count)               
             !             count = count+1
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                x_speed(i,j,1-l) = kmin_recv_buf(count) 
                count = count+1
             end do
             !             x_speed(i,j,0) = kmin_recv_buf(count)               
             !             count = count+1
             !             x_speed(i,j,-1) = kmin_recv_buf(count)               
             !             count = count+1
             !             x_speed(i,j,-2) = kmin_recv_buf(count)               
             !             count = count+1
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                y_speed(i,j,1-l) = kmin_recv_buf(count) 
                count = count+1
             end do
             !             y_speed(i,j,0) = kmin_recv_buf(count)               
             !             count = count+1
             !             y_speed(i,j,-1) = kmin_recv_buf(count)               
             !             count = count+1
             !             y_speed(i,j,-2) = kmin_recv_buf(count)               
             !             count = count+1
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                z_speed(i,j,1-l) = kmin_recv_buf(count) 
                count = count+1
             end do
             !             z_speed(i,j,0) = kmin_recv_buf(count)               
             !             count = count+1
             !             z_speed(i,j,-1) = kmin_recv_buf(count)               
             !             count = count+1
             !             z_speed(i,j,-2) = kmin_recv_buf(count)               
             !             count = count+1
          end do
       end do
    end if


    if(kmax_id >= 0) then
       !print *,"right id is ", process_id,right_id
       buf = (jmx+1)*(imx+1)*5*layers
       call MPI_RECV(kmax_recv_buf,buf,MPI_DOUBLE_PRECISION,kmax_id,1,MPI_COMM_WORLD,status,ierr)
       ! updating solution       
       count = 1
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                density(i,j,kmx+l-1) = kmax_recv_buf(count)
                count = count+1
             end do
             !             density(i,j,kmx) = kmax_recv_buf(count)
             !             count = count+1
             !             density(i,j,kmx+1) = kmax_recv_buf(count)
             !             count = count+1
             !             density(i,j,kmx+2) = kmax_recv_buf(count)
             !             count = count+1
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                pressure(i,j,kmx+l-1) = kmax_recv_buf(count)
                count = count+1
             end do
             !             pressure(i,j,kmx) = kmax_recv_buf(count)
             !             count = count+1
             !             pressure(i,j,kmx+1) = kmax_recv_buf(count)
             !             count = count+1
             !             pressure(i,j,kmx+2) = kmax_recv_buf(count)
             !             count = count+1
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                x_speed(i,j,kmx+l-1) = kmax_recv_buf(count)
                count = count+1
             end do
             !             x_speed(i,j,kmx) = kmax_recv_buf(count)
             !             count = count+1
             !             x_speed(i,j,kmx+1) = kmax_recv_buf(count)
             !             count = count+1
             !             x_speed(i,j,kmx+2) = kmax_recv_buf(count)
             !             count = count+1
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                y_speed(i,j,kmx+l-1) = kmax_recv_buf(count)
                count = count+1
             end do
             !             y_speed(i,j,kmx) = kmax_recv_buf(count)
             !             count = count+1
             !             y_speed(i,j,kmx+1) = kmax_recv_buf(count)
             !             count = count+1
             !             y_speed(i,j,kmx+2) = kmax_recv_buf(count)
             !             count = count+1
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                z_speed(i,j,kmx+l-1) = kmax_recv_buf(count)
                count = count+1
             end do
             !             z_speed(i,j,kmx) = kmax_recv_buf(count)
             !             count = count+1
             !             z_speed(i,j,kmx+1) = kmax_recv_buf(count)
             !             count = count+1
             !             z_speed(i,j,kmx+2) = kmax_recv_buf(count)
             !             count = count+1
          end do
       end do
       ! creating right send  buffer
       count = 1
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                kmax_send_buf(count) = density(i,j,kmx-l)
                count = count+1
             end do
             !             kmax_send_buf(count) = density(i,j,kmx-1)             
             !             count = count+1
             !             kmax_send_buf(count) = density(i,j,kmx-2)             
             !             count = count+1
             !             kmax_send_buf(count) = density(i,j,kmx-3)             
             !             count = count+1
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                kmax_send_buf(count) = pressure(i,j,kmx-l)
                count = count+1
             end do
             !             kmax_send_buf(count) = pressure(i,j,kmx-1)             
             !             count = count+1
             !             kmax_send_buf(count) = pressure(i,j,kmx-2)             
             !             count = count+1
             !             kmax_send_buf(count) = pressure(i,j,kmx-3)             
             !             count = count+1
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                kmax_send_buf(count) = x_speed(i,j,kmx-l)
                count = count+1
             end do
             !             kmax_send_buf(count) = x_speed(i,j,kmx-1)             
             !             count = count+1
             !             kmax_send_buf(count) = x_speed(i,j,kmx-2)             
             !             count = count+1
             !             kmax_send_buf(count) = x_speed(i,j,kmx-3)             
             !             count = count+1
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                kmax_send_buf(count) = y_speed(i,j,kmx-l)
                count = count+1
             end do
             !             kmax_send_buf(count) = y_speed(i,j,kmx-1)             
             !             count = count+1
             !             kmax_send_buf(count) = y_speed(i,j,kmx-2)             
             !             count = count+1
             !             kmax_send_buf(count) = y_speed(i,j,kmx-3)             
             !             count = count+1
          end do
       end do
       do i = 0, imx
          do j = 0, jmx
             do l = 1, layers
                kmax_send_buf(count) = z_speed(i,j,kmx-l)
                count = count+1
             end do
             !             kmax_send_buf(count) = z_speed(i,j,kmx-1)             
             !             count = count+1
             !             kmax_send_buf(count) = z_speed(i,j,kmx-2)             
             !             count = count+1
             !             kmax_send_buf(count) = z_speed(i,j,kmx-3)             
             !             count = count+1
          end do
       end do
       !do k = 1, buf
       !print *,'right send - ', process_id, k ,right_send_buf(k)
       !end do
       !print *, "total size is", buf 
       call MPI_SEND(kmax_send_buf,buf,MPI_DOUBLE_PRECISION,kmax_id,1,MPI_COMM_WORLD, ierr)    
    end if
  end subroutine send_recv_k

  subroutine allocate_buffer_cells(layers)
    implicit none
    integer :: buf
    integer, intent(in) :: layers
    call dmsg(1, 'parallel', 'allocating memory for buffer cells for MP')  
    !left buffer jmx+1 * kmx + 1 * (dens, x , y, z speeds, pressure)

    buf = (jmx+1)*(kmx+1)*5*layers ! size of buffer cells left - right
    call alloc(imin_send_buf, 1,buf, &
         errmsg='Error: Unable to allocate memory for buffer ' // &
         'variable left_buf.')
    call alloc(imax_send_buf, 1,buf, &
         errmsg='Error: Unable to allocate memory for buffer ' // &
         'variable right_buf.')
    call alloc(imin_recv_buf, 1,buf, &
         errmsg='Error: Unable to allocate memory for buffer ' // &
         'variable left_buf.')
    call alloc(imax_recv_buf, 1,buf, &
         errmsg='Error: Unable to allocate memory for buffer ' // &
         'variable right_buf.')

    buf = (imx+1)*(kmx+1)*5*layers ! size of buffer top - bottom
    call alloc(jmin_send_buf, 1,buf, &
         errmsg='Error: Unable to allocate memory for buffer ' // &
         'variable top_buf.')
    call alloc(jmax_send_buf, 1,buf, &
         errmsg='Error: Unable to allocate memory for buffer ' // &
         'variable bottom_buf.')
    call alloc(jmin_recv_buf, 1,buf, &
         errmsg='Error: Unable to allocate memory for buffer ' // &
         'variable top_buf.')
    call alloc(jmax_recv_buf, 1,buf, &
         errmsg='Error: Unable to allocate memory for buffer ' // &
         'variable bottom_buf.')

    buf = (imx+1)*(jmx+1)*5*layers ! size of buffer front back

    call alloc(kmin_send_buf, 1,buf, &
         errmsg='Error: Unable to allocate memory for buffer ' // &
         'variable front_buf.')
    call alloc(kmax_send_buf, 1,buf, &
         errmsg='Error: Unable to allocate memory for buffer ' // &
         'variable back_buf.')
    call alloc(kmin_recv_buf, 1,buf, &
         errmsg='Error: Unable to allocate memory for buffer ' // &
         'variable front_buf.')
    call alloc(kmax_recv_buf, 1,buf, &
         errmsg='Error: Unable to allocate memory for buffer ' // &
         'variable back_buf.')                     

  end subroutine allocate_buffer_cells


end module parallel
