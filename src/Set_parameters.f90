!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup
		
		mesh_file = 'Transient'
		fileplace_mesh = 'Mesh/Box/'
		scale_size_1 = 1e-3		
		
		material_coefficients_file = 'Material_9'
		fileplace_material = "Material/Data/Isotropic/"
		scale_prop_mat = 1e6
				
		transient = 1
		time_steps = 200
		dt = 1e-6
			
		density = 7850.d0
			
		load_profile='harmonic'
		omega = 0.99*79286.645975178
		phase = 0
		box_face_4 = (/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
		box_face_2 = (/100.d0*(1e6),1.d0,0.d0,1.d0,0.d0,1.d0/)

		results_file = 'Results'
		fileplace_results = 'Results/'
				
	END SUBROUTINE Setup 
	
END MODULE Set_parameters
!====================================================================================!!
