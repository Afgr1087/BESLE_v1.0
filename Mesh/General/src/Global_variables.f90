!===============================================================================!
!-------------------------------------------------------------------------------!
!            MODULE FOR THE GLOBAL VARIABLES OF THE OVERALL PROBLEM             !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Global_variables
!=======================================================================
!=======================================================================
!reading mesh file section
!=======================================================================
!=======================================================================
!-----------------------------------------------------------------------
!calsses and variables declaration
!this program is made in single precision, it means each value uses 
!32 bits (8 decimals)
!-----------------------------------------------------------------------
!class "point" contains three cordinates, each one is a real of 32bits 
	type point
    	real(8), dimension(3)  :: Coord
    end type point
!-----------------------------------------------------------------------
!class "element", it is a triangle, which contains 3 vertices 
!(the incidence), 3 nodes (the x, y, z coordinates of each node) and 
!1 normal (x, y, z coordinates of the normal vector)  
    type element
    	integer(4), dimension(3)        :: Vertices
        type(point), dimension(3)      :: Nodes
        type(point)                    :: Normal
    end type element
!-----------------------------------------------------------------------
!class "reg", it is a region, it contains a list of vertices(points), 
!a list of elements (element) and two counters for total number of 
!vertices and elements inside of the region
    type Reg
    	integer(4) :: n_Vertices, n_Elements
        type(point), allocatable   :: VerticesCoor(:)
        type(element), allocatable :: Elements(:)
    end type Reg
!-----------------------------------------------------------------------
    type ElaPara
	real(8) staticBCxVal, staticBCyVal, staticBCzVal, staticBCnVal, staticBClVal
    end type
!-----------------------------------------------------------------------
    type LinPara
	real(8) LinearBCxMaxVal, LinearBCyMaxVal, LinearBCzMaxVal, LinearBCnMaxVal, LinearBClMaxVal
    end type
!-----------------------------------------------------------------------
    type QuadPara
	real(8) QuadSOx, QuadSOy, QuadSOz, QuadSOn, QuadSOl, QuadSfx, QuadSfy, QuadSfz, QuadSfn, QuadSfl
	real(8) QuadAx, QuadAy, QuadAz, QuadAn, QuadAl, QuadBx, QuadBy, QuadBz, QuadBn, QuadBl
	real(8) QuadCx, QuadCy, QuadCz, QuadCn, QuadCl, QuadSAux
    end type
!-----------------------------------------------------------------------
    type SinePara
	real(8) SineSOx, SineSOy, SineSOz, SineSOn, SineSOl, SineSfx, SineSfy, SineSfz, SineSfn, SineSfl
	real(8) SineAx, SineAy, SineAz, SineAn, SineAl, SineBx, SineBy, SineBz, SineBn, SineBl
	real(8) SineCx, SineCy, SineCz, SineCn, SineCl, SineSAux
    end type
!-----------------------------------------------------------------------
!class "BoundaryCondition" is a group of elements that are contained in 
!a boundary condition. The class contains the number of elements, the 
!elements(vertices, nodes and normals) and the values of the 
!BC (x, y and z)
    type BoundaryCondition
    	integer(4) :: n_BCElements,n_BCVertices,BCMode
	integer(4), allocatable :: BCElements(:)
	integer(4) :: NSteps
	real(8) :: Area
	real(8), allocatable :: BCxVal(:), BCyVal(:), BCzVal(:)
	integer(4), dimension(3)  :: IDtype
	type(ElaPara) :: EParam
	type(LinPara) :: LParam
	type(QuadPara) :: QParam
	type(SinePara) :: SParam
	Character(20) BCxType, BCyType, BCzType, BCnType, DirType
        Character(20) FuncxType, FuncyType, FunczType, FuncnType, FunclType
	real(8), dimension(3) :: LoadDirection
	type(point), allocatable   :: VerticesCoor(:)
	type(element), allocatable :: Elements(:)
    end type
!-----------------------------------------------------------------------
!object region 
    type(Reg), allocatable :: Region(:)
    type(BoundaryCondition) :: BC(1000000)
!-----------------------------------------------------------------------
    integer(4) stat,i,j,k
    integer(4) :: nvaux, neaux, nraux
    integer(4) :: m, n, num_points_mesh
    integer(4) :: num_elem_mesh, size_eq
    integer(4) :: n_Regions, MeshNumPress = 5 
    character(30), dimension(3) :: block
    character(100) line
    character(150) lineaux
    character(12) cp
    character(150) filename_in,filename_out
    character(150) file_input_mesh
    character(150) filename_mesh
    real(8), parameter :: PI= 3.141592653589793d0
!normal variables
    integer(4) :: nv1, nv2, nv3
    real(8), dimension(3) :: V1, V2, V3
    real(8), dimension(3) :: NA, NB, NC
    real(8), dimension(3) :: Norm
!formatting file variable
    real(8), allocatable :: POINTS(:,:), NORMAL_VECTORS(:,:)
    real(8), allocatable :: POINTS_MESH(:,:),BOND_COND(:,:)
    integer(4), allocatable :: ELEM(:,:), FACES(:,:), El_reg(:,:)
    integer(4), allocatable :: EQ_NODE(:,:)
    integer(4), allocatable :: ELEM_MESH(:,:)
    real(8), allocatable  :: POINTS_int(:,:)
    real(8), dimension(3) :: BCsVertReal
    integer(4), dimension(3) :: BCsVertInt
    integer(4) :: NT_Elements, NT_Vertices, neaux2
!rearrange process
    integer(4) :: E1NV1, E1NV2, E1NV3
    integer(4) :: E2NV1, E2NV2, E2NV3

    real(8) :: E1V1X, E1V1Y, E1V1Z, Er1, Er2, Er3, Tol
    real(8) :: E1V2X, E1V2Y, E1V2Z, Er4, Er5, Er6, Tol2
    real(8) :: E1V3X, E1V3Y, E1V3Z, Er7, Er8, Er9

    real(8) :: E2V1X, E2V1Y, E2V1Z
    real(8) :: E2V2X, E2V2Y, E2V2Z
    real(8) :: E2V3X, E2V3Y, E2V3Z 
       
    real(8) :: x1, x2, x3
    real(8) :: y1, y2, y3
    real(8) :: z1, z2, z3

    real(8) :: x11, x22, x33
    real(8) :: y11, y22, y33
    real(8) :: z11, z22, z33

    real(8) :: scale_mult, scale_div
!Bondary conditions
    integer(4) :: steps, NumBCS
    character(150) BCSfilesPlace_in
    character(150) BCSfilename_in,BCSfilename_out, BCSfilename_aux
!bcs global id
    integer(4) :: V1BC, V2BC, V3BC
    real(8) :: C1V1BC, C2V1BC, C3V1BC
    real(8) :: C1V2BC, C2V2BC, C3V2BC
    real(8) :: C1V3BC, C2V3BC, C3V3BC
    integer(4) :: V1M, V2M, V3M
    real(8) :: C1V1M, C2V1M, C3V1M
    real(8) :: C1V2M, C2V2M, C3V2M
    real(8) :: C1V3M, C2V3M, C3V3M
    real(8) :: EArea, AuxArea1, AuxArea2, AuxArea3
    integer(4) :: NT_BCElements
!analysis setting
    Character(20) AnalysisType
    integer(4) Nsteps
    real(8) TimeInterval
    real(8), allocatable :: Time(:)
    



!===============================================================================!
END MODULE Global_variables
