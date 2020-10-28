
# This a common configuration file that includes definitions used by all
# the Makefiles.

# C++ compiler
CXX=g++
CC = gcc
# Flags for the C++ compiler
CFLAGS=-Wall -ansi -pedantic -O3 -Wformat

# Relative include and library paths for compilation of the examples
E_INC=-I voro++/src
INC=-I include
INC_t=-I triangle
E_LIB=-L voro++/src

# SRC is the directory in which the C source files are, and BIN is the
#   directory where you want to put the executable programs.  By default,
#   both are the current directory.

SRC = ./src/
SRC_t= ./triangle/
BIN = ./
# TRILIBDEFS is a list of definitions used to compile an object code version
#   of Triangle (triangle.o) to be called by another program.  The file
#   "triangle.h" contains detailed information on how to call triangle.o.
TRILIBDEFS = -DTRILIBRARY

# Installation directory
PREFIX=/usr/local

# Install command
INSTALL=install

# Flags for install command for executable
IFLAGS_EXEC=-m 0755

# Flags for install command for non-executable files
IFLAGS=-m 0644

# List of executables
PROGRAM = Structure Mesh

all: $(PROGRAM) move
recompile: move_back $(PROGRAM) move2
fix: $(PROGRAM) move2
   
# Makefile rules in C++

Structure: Mesh_Generator.o Set_parameters.o Polycrystalline_Material.o Input_data.o Ouput_data.o 
	$(CXX) $(CFLAGS) $(INC) $(E_INC) $(E_LIB) -o Structure Mesh_Generator.o Set_parameters.o Polycrystalline_Material.o Input_data.o Ouput_data.o ./voro++/src/libvoro++.a

Mesh_Generator.o: $(SRC)Mesh_Generator.cc $(SRC)Set_parameters.cc $(SRC)Polycrystalline_Material.cc $(SRC)Input_data.cc $(SRC)Ouput_data.cc 
	$(CXX) $(CFLAGS) $(INC) $(E_INC) $(E_LIB) -c $(SRC)Mesh_Generator.cc 

Set_parameters.o: $(SRC)Set_parameters.cc
	$(CXX) $(CFLAGS) $(INC) $(E_INC) $(E_LIB) -c $(SRC)Set_parameters.cc

Polycrystalline_Material.o: $(SRC)Polycrystalline_Material.cc
	$(CXX) $(CFLAGS) $(INC) $(E_INC) $(E_LIB) -c $(SRC)Polycrystalline_Material.cc 

Input_data.o: $(SRC)Input_data.cc
	$(CXX) $(CFLAGS) $(INC) $(E_INC) $(E_LIB) -c $(SRC)Input_data.cc 

Ouput_data.o: $(SRC)Ouput_data.cc 
	$(CXX) $(CFLAGS) $(INC) $(E_INC) $(E_LIB) -c $(SRC)Ouput_data.cc 

# Make rule for Triangle Mesh generator library in C

$(BIN)Mesh: $(SRC)Main2.c $(BIN)triangle.o $(BIN)Input_data_mesh.o $(BIN)2D_Data.o $(BIN)transf_coordinates.o $(BIN)Output_data.o
	$(CC) $(INC_t) -o $(BIN)Mesh $(SRC)Main2.c \
		$(BIN)triangle.o $(BIN)Input_data_mesh.o $(BIN)2D_Data.o $(BIN)transf_coordinates.o $(BIN)Output_data.o -lm

$(BIN)triangle.o: $(SRC_t)triangle.c $(SRC_t)triangle.h 
	$(CC) $(TRILIBDEFS) -c -o $(BIN)triangle.o \
		$(SRC_t)triangle.c

$(BIN)Input_data_mesh.o: $(SRC)Input_data_mesh.c 
	$(CC) -c -o $(BIN)Input_data_mesh.o \
		$(SRC)Input_data_mesh.c

$(BIN)2D_Data.o: $(SRC)2D_Data.c $(SRC)transf_coordinates.c
	$(CC) -c -o $(BIN)2D_Data.o \
		$(SRC)2D_Data.c

$(BIN)transf_coordinates.o: $(SRC)transf_coordinates.c
	$(CC) -c -o $(BIN)transf_coordinates.o \
		$(SRC)transf_coordinates.c

$(BIN)Output_data.o: $(SRC)Output_data.c
	$(CC) -c -o $(BIN)Output_data.o \
		$(SRC)Output_data.c

# Clean executables and .o files
move:
	mkdir ./IO_files
	mkdir ./Export
	mkdir ./obj
	mv *.o obj/

move2:
	mv *.o obj/

move_back:
	mv obj/*.o ./

clean:
	rm -f obj/*.o $(PROGRAM) 
	rm -rf obj IO_files Export
	rm -f *.o $(PROGRAMS)
.PHONY: all clean
