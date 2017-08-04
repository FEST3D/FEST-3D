
  case('Mass_abs')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(0)

  case('Resnorm_abs')
    write(RESNORM_FILE_UNIT, frm, advance='no') sqrt(sum(Res_abs(1:)**2))

  case('Viscous_abs')
    write(RESNORM_FILE_UNIT, frm, advance='no') sqrt(sum(Res_abs(1:5)**2))

  case('Turbulent_abs')
    write(RESNORM_FILE_UNIT, frm, advance='no') sqrt(sum(Res_abs(6:)**2))

  case('Continuity_abs')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(1)

  case('X-mom_abs')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(2)

  case('Z-mom_abs')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(3)

  case('Y-mom_abs')
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
    write(RESNORM_FILE_UNIT, frm, advance='no') sqrt(sum(Res_rel(6:)**2))

  case('Continuity_rel')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(1)

  case('X-mom_rel')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(2)

  case('Z-mom_rel')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(3)

  case('Y-mom_rel')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(4)

  case('Energy_rel')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(5)

  case('TKE_abs')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(6)

  case('tv_abs')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(6)

  case('Dissipation_abs')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(7)

  case('Omega_abs')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(7)

  case('Kl_abs')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_abs(7)

  case('TKE_rel')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(6)

  case('tv_rel')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(6)

  case('Dissipation_rel')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(7)

  case('Omega_rel')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(7)

  case('Kl_rel')
    write(RESNORM_FILE_UNIT, frm, advance='no') Res_rel(7)

  case DEFAULT
    ! making absolute resnorm default
    write(RESNORM_FILE_UNIT, frm, advance='no') sqrt(sum(Res_abs(1:)**2))
    Issure_warning
