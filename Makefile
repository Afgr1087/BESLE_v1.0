################################################################################
# MAKEFILE - Just edit if you know what you are doing!! 
################################################################################

# Directory of MUMPS library
topdir = ./MUMPS
libdir = $(topdir)/lib

# Compiler
FC = mpif90 

include $(topdir)/Makefile.inc

LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)

LIBDMUMPS = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

# Flags for debugging or for maximum performance (comment as necessary)
OPTF = -O3 -ffast-math -Dintel_ -DALLOW_NON_INIT 

FLAGS1 = -g -O0 -fbounds-check -fopenmp -fimplicit-none -fopenmp #-fdefault-integer-8
FLAGS2 =  -Waliasing -Wsurprising -fPIC -fno-range-check#-Wall -pedantic -Wextra 
FCFLAGS = $(FLAGS1) $(FLAGS2) 

# List of executables to be built 
PROGRAMS = BESLE

# "make" builds all
all: $(PROGRAMS) move
recompile: move_back $(PROGRAMS) move2
fix: $(PROGRAMS) move2

BESLE: Global_variables.o Set_parameters.o Global_functions.o Input.o Failure_Analysis.o Discretization.o BC_Module.o Int_points.o Fourier_Coefficients.o DRM.o Comp_H_G.o Vectors_b.o General_Matrices.o General_Sparce_Vectors.o Module_Solver.o Response_D_T.o Response_S_S.o Homogenization.o Output.o

# Building prog from prog.o; $^ (GNU extension) is
%: %.o
	@echo ' '
	@echo 'Building target: $@'
	@echo 'Invoking: GNU Fortran Linker'
	$(FC) $(FCFLAGS) -o $@ $(OPTF) $^ $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS) 
	@echo 'Finished building target: $@'
	@echo ' '
	
# Building prog.o from prog.f90 
%.o: src/%.f90
	@echo ' '
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	$(FC) $(FCFLAGS) $(OPTF) -I. -I$(topdir)/include -I$(topdir)/src -c $< 
	@echo 'Finished building: $<'
	@echo ' '


	
# Utility targets
.PHONY: clean veryclean

move:
	mkdir ./OOC
	mkdir ./obj
	mkdir ./Results
	mv *.mod obj/
	mv *.o obj/

move2:
	mv *.mod obj/
	mv *.o obj/
	rm -f ./OOC/mumps*

move_back:
	mv obj/*.mod ./
	mv obj/*.o ./ 

clean:
	rm -rf obj Results OOC 
	rm -f *.o *.mod *.MOD $(PROGRAMS) 
