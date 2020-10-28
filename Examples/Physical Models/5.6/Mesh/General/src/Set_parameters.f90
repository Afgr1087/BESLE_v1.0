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
    	Nsteps=4  

		NumBCS = 4   
    	MeshNumPress = 4 

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
		BC(2)%BCxType="displacement"     
		BC(2)%BCyType="displacement"      
		BC(2)%BCzType="traction"
		BC(2)%FuncxType="linear"
		BC(2)%FuncyType="linear"
		BC(2)%FunczType="sine"
		BC(2)%LParam%LinearBCxMaxVal=0
		BC(2)%LParam%LinearBCyMaxVal=0
		BC(2)%SParam%SineAz=0.268d0
		BC(2)%SParam%SineBz=1.d0
		BC(2)%SParam%SineCz=0.d0
		BC(2)%SParam%SineSOz=0.d0
		BC(2)%SParam%SineSfz=pi/2.d0

		BC(3)%DirType="normal"              
		BC(3)%BCnType="traction"
		BC(3)%FuncnType="sine"
		BC(3)%SParam%SineAn=0.201d0
		BC(3)%SParam%SineBn=1.d0
		BC(3)%SParam%SineCn=0.d0
		BC(3)%SParam%SineSOn=0.d0
		BC(3)%SParam%SineSfn=pi/2.d0
		
		BC(4)%DirType="normal"              
		BC(4)%BCnType="traction"
		BC(4)%FuncnType="sine"
		BC(4)%SParam%SineAn=0.201d0
		BC(4)%SParam%SineBn=1.d0
		BC(4)%SParam%SineCn=0.d0
		BC(4)%SParam%SineSOn=0.d0
		BC(4)%SParam%SineSfn=pi/2.d0

    
!-----------------------------------------------------------------------

	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
