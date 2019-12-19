!--------------------------------------
! aim: 1) calculate mass balance
!      2) set flux to zero if wall
!---------------------------------------
  id = 9999  ! this number does not have meaning but it must be diffent than any boundary condition
  if(i_f==1 .and. i==1  ) then
    id=imin_id; normal_tag=-1
  end if
  if(i_f==1 .and. i==imx) then
    id=imax_id; normal_tag= 1
  end if
  if(j_f==1 .and. j==1  ) then
    id=jmin_id; normal_tag=-1
  end if
  if(j_f==1 .and. j==jmx) then
    id=jmax_id; normal_tag= 1
  end if
  if(k_f==1 .and. k==1  ) then
    id=kmin_id; normal_tag=-1
  end if
  if(k_f==1 .and. k==kmx) then
    id=kmax_id; normal_tag= 1
  end if

    select case(id)
      case(-5,-6,-7)
        F_plus(1)=0.0
        F_minus(1)=0.0
!        if(i==1)then
!          x_qp_left(1,:,:,2:4)=0.0
!          x_qp_right(1,:,:,2:4)=0.0
!        elseif(j==1)then
!          y_qp_left(:,1,:,2:4)=0.0
!          y_qp_right(:,1,:,2:4)=0.0
!        elseif(k==1)then
!          z_qp_left(:,:,1,2:4)=0.0
!          z_qp_right(:,:,1,2:4)=0.0
!        elseif(i==imx)then
!          x_qp_left(imx,:,:,2:4)=0.0
!          x_qp_right(imx,:,:,2:4)=0.0
!        elseif(j==jmx)then
!          y_qp_left(:,jmx,:,2:4)=0.0
!          y_qp_right(:,jmx,:,2:4)=0.0
!        elseif(k==kmx)then
!          z_qp_left(:,:,kmx,2:4)=0.0
!          z_qp_right(:,:,kmx,2:4)=0.0
!        end if
      case(-4,-3,-2,-1,-8,-9)
        merror = merror - normal_tag*(F_plus(1)+F_minus(1))*fA(i,j,k)
       ! print*, merror, normal_tag, (F_plus(1)+F_minus(1))*fA(i,j,k)
      case(0:)
        continue
      case DEFAULT
        continue
    end select


