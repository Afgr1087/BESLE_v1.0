!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup
		
		mesh_file = 'Mesh'
		fileplace_mesh = 'Mesh/Polycrystal/Export/'
		
		n_Materials = 96
		material_coefficients_database = 'Fe'
		fileplace_material = "Material/Data/Anisotropic/Orientations/"
		scale_prop_mat = 1.d0*(1e3)
		scale_size_1 = 0.001d0			
		
		transient = 1
		time_steps = 50
		dt = 5e-9
		density = 7874.d0/(1000**3)

		box_face_5 = (/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
		box_face_6 = (/0.d0,1.d0,0.d0,1.d0,200.d0*(1e3),1.d0/)

		scale_results = 200
		results_file = 'Results'
		fileplace_results = 'Results/'	
		
		
	END SUBROUTINE Setup 
	
END MODULE Set_parameters
!====================================================================================!!
