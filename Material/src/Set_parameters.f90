!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup
	
		fileplace="Data/Isotropic/"
		file_name="Material_9.dat"

		Material="Isotropic" 

		E=200000.d0; nu=0.d0;

	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
