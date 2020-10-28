!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup

		fileplace="Data/Isotropic/"
		file_name="Material_7.dat"

		Material="Isotropic" 

		E=210000.d0; nu=0.3d0;

	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
