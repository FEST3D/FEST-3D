Set(Source 
vartypes.f90
lusgs.f90 
plusgs.f90 
layout.f90 
resnorm.f90 
interface1.f90 
time.f90 
viscous.f90 
grid.f90 
viscosity.f90 
solver.f90 
mapping.f90 
update.f90 
geometry.f90 
source.f90 
CC.f90
gradients.f90 
state.f90 
global/global_sa.f90 
global/global_kkl.f90 
global/global_sst.f90 
boundary/copy_bc.f90 
boundary/FT_bc.f90 
boundary/read_bc.f90 
boundary/boundary_state_reconstruction.f90 
boundary/bc_primitive.f90 
boundary/bc.f90 
read_write/read/read_output_tec.f90 
read_write/read/read_output_vtk.f90 
read_write/read/read_output.f90 
read_write/read/read.f90 
read_write/write/write_output_vtk.f90 
read_write/write/write_output.f90 
read_write/write/dump_solution.f90 
read_write/write/write_output_tec.f90 
read_write/write/write_output_tec_node.f90 
wall/wall.f90 
wall/wall_dist.f90 
convergence.f90 
read_write/check_output_control.f90 
utils.f90 
face/state/weno.f90 
face/state/weno_NM.f90 
face/state/muscl.f90 
face/state/face_interpolant.f90 
face/state/ppm.f90 
face/flux/convective/van_leer.f90 
face/flux/convective/ldfss0.f90 
face/flux/convective/ausm.f90 
face/flux/convective/ausmP.f90
face/flux/convective/ausmUP.f90
face/flux/convective/slau.f90
face/flux/convective/scheme.f90 )


#FILE(GLOB_RECURSE Source *.f90)

Add_library(lib OBJECT ${Source})
Add_executable(${EXE} ${SRC}/main.f90 $<TARGET_OBJECTS:lib>)
target_link_libraries(${EXE} ${lib} ${MPI_Fortran_LIBRARIES})
