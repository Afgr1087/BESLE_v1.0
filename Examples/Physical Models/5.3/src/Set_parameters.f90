!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup
		
		mesh_file = 'Input_data'
		fileplace_mesh = 'Mesh/Box/'
		
		material_coefficients_file = 'Material_7'
		fileplace_material = "Material/Data/Isotropic/"
						
		density = 7874.d0/(1000**3)

		bodyforces = 1
		bodyforce = (/0.d0,0.d0,-7.72439d0*(1e-5)/)

		quasi_static = 1
		static_steps = 1
		
		box_face_5 = (/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)

		scale_results = 4e7
		results_file = 'Results'
		fileplace_results = 'Results/'		
		
		
	END SUBROUTINE Setup 
	
END MODULE Set_parameters
!====================================================================================!!
