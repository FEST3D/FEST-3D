#compiler options
program := FEST3D
subprog := wall
FC 			 = $(IFC)
FFLAGS	 = $(IFFLAGS)
VPATH 	:= src
BINDIR	:= bin
OBJDIR	:= obj
MODDIR	:= mod
GFC 		:= mpif90
IFC 		:= mpif90 -f90=ifort
GFFLAGS	:= -O3 -Wall -Wextra -Wconversion \
					-Wno-compare-reals \
					-fdefault-real-8 \
					-Waliasing \
					-Wsurprising \
					-Wintrinsic-shadow \
					-pedantic-errors \
					-fbounds-check \
					-cpp
#					-Werror \

IFFLAGS := -O3 -free -r8 -traceback \
           -cpp \
					 -zero \
					 -u \
					 -W1 \
					 -warn \
					 -warn error \
					 -warn declaration \
					 -warn uncalled \
					 -warn unused
subpath:= resnorm \
				 resnorm/viscous \
				 resnorm/turbulent

#files = global.f90 \
#				global_vars.f90 \
#				utils.f90 \
#				layout.f90 \
#				bitwise.f90 \
#				string.f90 \
#				read.f90 \
#				write_output_vtk.f90 \
#				grid.f90 \
#				geometry.f90 \
#				state.f90 \
#				ppm.f90 \
#				muscl.f90 \
#				face_interpolant.f90 \
#				van_leer.f90 \
#				ausm.f90 \
#				ldfss0.f90 \
#				scheme.f90 \
#				wall_find.f90 \
#				surfnode.f90 \
#				wall_dist.f90 \
#				source.f90 \
#				viscous.f90 \
#				parallel.f90 \
#				boundary_conditions.f90 \
#				boundary_state_reconstruction.f90 \
#				res_viscous.f90 \
#				res_turbulent.f90 \
#				resnorm/resnorm.f90 \
#				dump_solution.f90 \
#				solver.f90 \
#				main.f90
include files

subfiles = global.f90 \
					 global_vars.f90 \
					 utils.f90 \
					 grid.f90 \
					 wall_find.f90 \
					 surfnode.f90 \
					 wall_main.f90 


objects = $(addprefix $(OBJDIR)/, $(files:.f90=.o))
modules = $(addprefix $(MODDIR)/, $(files:.f90=.mod))
src     = $(addprefix $(VPATH)/, $(files))
subobjects = $(addprefix $(OBJDIR)/, $(notdir $(subfiles:.f90=.o)))
submodules = $(addprefix $(MODDIR)/, $(notdir $(subfiles:.f90=.mod)))
exe = $(BINDIR)/$(program)
pre = $(BINDIR)/$(subprog)

#path to python program that make dependency
MAKEDEPEND = python fort_depend.py
# $(DEP_FILE) is a .dep file generated by fort_depend.py
DEP_FILE = makefile.dep

all: $(exe) $(DEP_FILE) $(pre)

# Make dependencies
.PHONY: depend
depend: $(DEP_FILE)

# The .dep file depends on the source files, so it automatically gets updated
# when you change your source
$(DEP_FILE): $(OBJECTS)
	@echo "Making dependencies!"
	$(MAKEDEPEND) -w -o $(DEP_FILE) -f $(src)

include $(DEP_FILE)
         
$(exe) : $(objects)
ifeq ($(FC), mpif90)
	$(FC)  $(objects) -o $(exe) -J$(MODDIR)/
else
	$(FC)  $(objects) -o $(exe) -module $(MODDIR)
endif

$(OBJDIR)/%.o : %.f90 | setup 
ifeq ($(FC), mpif90)
	$(FC) $(FFLAGS) -c  $< -o $@ -J$(MODDIR)/
else
	$(FC) $(FFLAGS) -c  $< -o $@ -module $(MODDIR)/
endif

$(pre) : $(subobjects)
ifeq ($(FC), mpif90)
	$(FC)  $(subobjects) -o $(pre) -J$(MODDIR)/
else
	$(FC)   $(FFLAGS) -c src/wall_main.f90  -module $(MODDIR)
	$(FC)  $(subobjects) -o $(pre) -module $(MODDIR)
endif

setup:
	@if test ! -d $(OBJDIR); then mkdir $(OBJDIR); else : ; fi
	@if test ! -d $(MODDIR); then mkdir $(MODDIR); else : ; fi
	@if test ! -d $(BINDIR); then mkdir $(BINDIR); else : ; fi
	@if test ! -d $(OBJDIR)/resnorm/viscous; then mkdir $(OBJDIR)/resnorm/viscous; else : ; fi
	@if test ! -d $(OBJDIR)/resnorm/turbulent; then mkdir $(OBJDIR)/resnorm/turbulent; else : ; fi

clean:
	rm -f $(objects) $(modules) $(exe)

tar:
	tar -zcvf fest3d.tar.gz $(OBJDIR)/ $(MODDIR)/ $(BINDIR)/ $(VPATH)/ setup.sh makefile \
		makefile.dep fort_depend.py archive/ 

