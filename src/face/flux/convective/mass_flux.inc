!--------------------------------------
! aim: 1) calculate mass balance
!      2) set flux to zero if wall
!---------------------------------------
  id = 9999  ! this number does not have meaning but it must be diffent than any boundary condition
  if(i_f==1 .and. i==1  ) id=imin_id
  if(i_f==1 .and. i==imx) id=imax_id
  if(j_f==1 .and. j==1  ) id=jmin_id
  if(j_f==1 .and. j==jmx) id=jmax_id
  if(k_f==1 .and. k==1  ) id=kmin_id
  if(k_f==1 .and. k==kmx) id=kmax_id

    select case(id)
      case(-5,-6,-7)
        F_plus(1)=0.
        F_minus(1)=0
      case(-1,-3)
        merror = merror + (F_plus(1)+F_minus(1))*fA(i,j,k)
      case(-2,-4,-8)
        merror = merror - (F_plus(1)+F_minus(1))*fA(i,j,k)
      case(0:)
        continue
      case DEFAULT
        continue
    end select


