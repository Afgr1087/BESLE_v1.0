!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup
		
		Mesh_file = 'Mesh'
		fileplace_mesh = 'Mesh/Polycrystal/Export/'
			
		n_Materials = 96
		material_coefficients_database = 'Zn'
		fileplace_material = "Material/Data/Anisotropic/Orientations/"
		scale_size_1 = 0.001d0			

		quasi_static = 1
		static_steps = 20
		
	    load_profile='ramp'
		box_face_5 = (/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
		box_face_6 = (/0.d0,1.d0,0.d0,1.d0,0.00005d0,0.d0/)

		crack=1

		scale_results = 10
		results_file = 'Results'
		fileplace_results = 'Results/'	
		
		
	END SUBROUTINE Setup 
	
END MODULE Set_parameters
!====================================================================================!!
