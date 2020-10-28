!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup
		
		mesh_file = 'Input_data'
		fileplace_mesh = 'Mesh/General/Meshes/DAT/'
	    		
		n_Materials = 57
		material_coefficients_database = 'Composite'
		scale_size_1 = 0.001d0			

		quasi_static = 1
		static_steps = 50
		
		bc_groups = 2	
		fileplace_BCs = 'Mesh/General/BCs/BCs_DAT/'	
		
		scale_results = 200
		results_file = 'Results'
		fileplace_results = 'Results/'			
		
	END SUBROUTINE Setup 
	
END MODULE Set_parameters
!====================================================================================!!
