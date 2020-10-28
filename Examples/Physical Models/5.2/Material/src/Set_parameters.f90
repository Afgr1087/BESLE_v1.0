!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup

		n_materials=100

		fileplace="Data/Anisotropic/Orientations/Fe/"
		file_name="Angles.dat"

		Material="Anisotropic" 
		Lattice="cubic" 

		C11=230000.d0; 
		C12=135000.d0; 
		C44=117000.d0;

	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
