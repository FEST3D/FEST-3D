!should be included in state file, in subroutine name same as this file name 

            !-----------------------------------------------------------
            ! Read the k and omega state of the system from a file
            !-----------------------------------------------------------

            call dmsg(1, 'turbulence_models/sst/state', 'readstate_vtk.inc')

            !only read k and omega from file
            !as rest is read in laminar solver
            read(OUT_FILE_UNIT, *) ! Skip Blank line
            read(OUT_FILE_UNIT, *) ! Skip SCALARS k
            read(OUT_FILE_UNIT, *) ! Skip LOOKUP_TABLE
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                read(OUT_FILE_UNIT, *) tk(i, j, k)
              end do
             end do
            end do

            read(OUT_FILE_UNIT, *) ! Skip Blank line
            read(OUT_FILE_UNIT, *) ! Skip SCALARS omega
            read(OUT_FILE_UNIT, *) ! Skip LOOKUP_TABLE
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                read(OUT_FILE_UNIT, *) tw(i, j, k)
              end do
             end do
            end do
