################################################################################
# MAKEFILE - Just edit if you know what you are doing!! 
################################################################################

# Compiler
FC = gfortran 

# Flags for debugging or for maximum performance (comment as necessary)
FLAGS1 = -g -O0 -fbounds-check -fopenmp -fimplicit-none -fopenmp #-fdefault-integer-8
FLAGS2 =  -Waliasing -Wsurprising -fPIC -fno-range-check#-Wall -pedantic -Wextra 

FCFLAGS = $(FLAGS1) $(FLAGS2) 

# Libraries needed for linking
LDFLAGS = -llapack

# List of executables to be built 
PROGRAMS = F_Coeff

# "make" builds all
all: $(PROGRAMS) move	
recompile: move_back $(PROGRAMS) move
fix: $(PROGRAMS) move

F_Coeff: Global_variables.o Set_parameters.o Fourier_Coefficient.o Output_data.o

# Building prog from prog.o; $^ (GNU extension) is
%: %.o
	@echo 'Building target: $@'
	@echo 'Invoking: GNU Fortran Linker'
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
	@echo 'Finished building target: $@'
	@echo ' '

# Building prog.o from prog.f90 
%.o: src/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	$(FC) $(FCFLAGS) -c $<
	@echo 'Finished building: $<'
	@echo ' '

# Utility targets
.PHONY: clean veryclean

move:
	mkdir ./obj
	mv *.mod obj/
	mv *.o obj/

move_back:
	mv obj/*.mod ./
	mv obj/*.o ./

clean:
	rm -f *.o *.mod *.MOD $(PROGRAMS)
	rm -rf obj 

