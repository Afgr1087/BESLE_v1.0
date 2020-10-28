!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup

		fileplace = "Data/Others/Composite/"
		file_name = "Mat_2.dat"

		Material = 'Isotropic' 
	

		E = 5100.d0;
		nu = 0.40d0;	

	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
