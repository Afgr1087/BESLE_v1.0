!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup

	n_materials = 8

	fileplace = "Data/Isotropic/Multiple/Simulation_1/"
	file_name = "E_v_constants.dat"

	Material = "Multiple_iso" 
	
	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
