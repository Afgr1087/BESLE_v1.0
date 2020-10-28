!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup

		filename_in = 'Meshes/Obj/mesh.obj' 
    	filename_out = 'Input_data.dat'
    	BCSfilesPlace_in = 'BCs/BCs_Obj/'

		AnalysisType="quasi-elastostatic" 
    	Nsteps=50  

		NumBCS = 2   

  		BC(1)%DirType="xyz"              
		BC(1)%BCxType="displacement"     
		BC(1)%BCyType="displacement"      
		BC(1)%BCzType="displacement"
		BC(1)%FuncxType="linear"
		BC(1)%FuncyType="linear"
		BC(1)%FunczType="linear"
		BC(1)%LParam%LinearBCxMaxVal=0
		BC(1)%LParam%LinearBCyMaxVal=0
		BC(1)%LParam%LinearBCzMaxVal=0

		BC(2)%DirType="xyz"              
		BC(2)%BCxType="traction"     
		BC(2)%BCyType="free"      
		BC(2)%BCzType="traction"
		BC(2)%FuncxType="linear"
		BC(2)%FunczType="linear"
		BC(2)%LParam%LinearBCxMaxVal=100
		BC(2)%LParam%LinearBCzMaxVal=100

    
!-----------------------------------------------------------------------

	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
