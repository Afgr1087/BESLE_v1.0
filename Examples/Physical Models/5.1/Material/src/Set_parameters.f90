!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup

		n_materials=100

		fileplace="Data/Anisotropic/Orientations/Zn/"
		file_name="Angles.dat"

		Material="Anisotropic" 
		Lattice="hcp" 

		C11=165000.d0; C12=31100.d0; 
        C13=50000.d0; C33=61800.d0; 
        C44=33600.d0

	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
