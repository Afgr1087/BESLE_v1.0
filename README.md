# BESLE_v1.0

A MPI-parallelized Fortran 3D

Boundary Element Software for 3D Linear Elasticity

By Andres F. Galvis 

Department of Computational Mechanics

School of Mechanical Engineering

University of Campinas

09/10/2020

====================================================================================================

BESLE is a software for carrying out 3D simulations of solids under 
quasi-static, inertial and high-rate boundary conditions for transient 
analysis. It uses the 3D elastostatic and dynamic boundary element formulation
with the fundamental solution based on double Fourier series. The software 
supports simulations of isotropic and anisotropic bodies considering, if 
necessary, several domains for modeling heterogeneous materials composed of 
different constituents. A distinguishing feature of the BESLE is that it 
carries out customized solid simulations in a straightforward configuration 
using only the surface mesh information. The software accounts with external 
sub-package for creating a material database, several options for mesh generation, 
diverse ways for the configuration of the boundary conditions. 

BESLE is composed by various Fortran 90 modules that can be compiled and 
modified from a Setup subroutine utility that is configured following the 
instruction provided in the User's Guide. This software employs external free
codes and software such as MUMPS, ScaLAPACK, SCOTCH, BLAS, LAPACK, Voro++ and 
triangle. The user must download and install them according to the instruction 
given in the INSTALL file. In order to run the different examples shown in the 
use User's Guide, this application requires the previous installation of MPICH2.
Moreover, this software accounts only with a parallel version, therefore, at least
2 threads must be indicated by the execution. For simulations with a large number
of degrees of freedom, MUMPS is installed in its 64-bit version.

Please read this README file and the documentation for a complete list of 
functionalities. Documentation and publications related to BESLE can also 
be found at https://github.com/Afgr1087/BESLE_v1.0.git. Please refer to INSTALL 
for installation instructions.

  This version of BESLE is provided to you free of charge. It is
  released under the GNU General Public License v3.0.


Contents:
-----------------------------------------------------------------------------------

INSTALL     README     VERSION		Makefile    
Make.inc/   doc/       src/			Mesh/    
Material/   Cracks/    Body_forces/ Examples/    

Make.inc/   contains all the Makefile.inc necessary for the installation of MUMPS, 
	        SCOTCH and ScaLAPACK.

doc/        contains the users' guide in pdf format.

src/        source *.f90 files of the main BESLE code.

Mesh/	    options for the mesh generation are provided in this folder
	  	    i) General using 3dxMax and ii) Polycrystal and iii) Matlab box.

Material/   is the folder that contains a serial Frotran code to generate the 
	  		material database.

Cracks/		contains the *.dat file with pre-cracks.

Body_forces/ contains the *.dat file with a body force vector.

Examples/	 contains the mesh, material, and simulation examples.


Additional contents after installation:
-----------------------------------------------------------------------------------

MUMPS/		The MUMPS solver must be installed in its 64-bit version	

Polycrystal/voro++/ 	3D Voronoi structure generator.

Polycrystal/triangle/	two-dimensional quality mesh generator and Delaunay triangulator.


The libraries voro++ and triangle have to be installed in case the user wants to
reproduce artificial polycrystalline structures. If it is not the case, their 
installation could be skipped. After BESLE compilation, temporary folders are generated 
as follows.

obj/ 	contains the *.o and *.mod files produced by the general compilation of BESLE.

OOC/	is the folder that MUMPS uses to write data while the system of equations is 
		being solved. After a success solution, MUMPS deletes these files itself.

Results/	contains the .vtk files for visualization using Paraview.

-----------------------------------------------------------------------------------

Acknowledgments
===============
The authors would like to thank to the University of Campinas (Brazil), Brunel University London (UK), 
and the University of Portsmouth (UK) for the facilities and structure provided to develop this version 
of BESLE. The project was funded by the National Council for Scientific and Technological Development 
CNPq (Grant Numbers: 312493/2013-4, 154283/2014-2 and 312536/2017-8), and the Brazilian Coordination for 
the Improvement of Higher Education Personnel CAPES (Grant Number: 435214/2019-01). 

Luiz C. Wrobel also thanks the CNPq for his personal financial support (Grant Number: 303770/2019-8).

This material is based upon work supported by the Air Force Office of Scientific 
Research-AFOSR under Award Numbers FA9550-18-1-0113 and FA9550-20-1-0133.

Computational sources were provided by the Center for Computational Engineering and Science-CCES 
at the University of Campinas funded by the Sao Paulo Research Foundation FAPESP (Grant Number: 2013/08293-7).

====================================================================================================

====================================================================================================
