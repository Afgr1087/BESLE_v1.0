!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup

		fileplace = "Data/Others/Composite/"
		file_name = "Mat_1.dat"

		Material = 'Isotropic' 

		E = 241000.d0;
		nu = 0.26d0;

	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
