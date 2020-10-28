!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup
		
		mesh_file = 'Input_data'
		fileplace_mesh = 'Mesh/General/Meshes/DAT/'
	    		
		n_Materials = 15
		material_coefficients_database = 'Vertebra'		
		fileplace_material = "Material/Data/Others/"
		
		quasi_static = 1
		static_steps = 4
		
		bc_groups = 4	
		fileplace_BCs = 'Mesh/General/BCs/BCs_DAT/'	
			    
		scale_results = 10
		results_file = 'Results'
		fileplace_results = 'Results/'			
		
	END SUBROUTINE Setup 
	
END MODULE Set_parameters
!====================================================================================!!
