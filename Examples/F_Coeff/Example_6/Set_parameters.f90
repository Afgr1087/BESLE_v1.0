!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup

	n_materials = 8

	fileplace = "Data/Anisotropic/Multiple/Simulation_1/"
	file_name = "C_tensors.dat"

	Material = "Multiple_aniso" 

	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
