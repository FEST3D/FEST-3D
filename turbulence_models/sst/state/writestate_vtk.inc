!should be included in state file, in subroutine name same as this file name 

            !-----------------------------------------------------------
            ! Write the state of the system to a file
            !-----------------------------------------------------------

            
            call dmsg(1, 'turbulence_models/sst/state', 'twritestate_vtk.inc')

            !!!!!!!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! k and omega will be append in existing file
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            ! Writing k (Turbulent Kinectic enrgy)
            write(OUT_FILE_UNIT, '(a)') 'SCALARS k FLOAT'
            write(OUT_FILE_UNIT, '(a)') 'LOOKUP_TABLE default'
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                write(OUT_FILE_UNIT, fmt='(f0.16)') tk(i, j, k)
              end do
             end do
            end do
            write(OUT_FILE_UNIT, *) 

            ! Writing Omega
            write(OUT_FILE_UNIT, '(a)') 'SCALARS Omega FLOAT'
            write(OUT_FILE_UNIT, '(a)') 'LOOKUP_TABLE default'
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                write(OUT_FILE_UNIT, fmt='(f0.16)') tw(i, j, k)
              end do
             end do
            end do
            write(OUT_FILE_UNIT, *) 
