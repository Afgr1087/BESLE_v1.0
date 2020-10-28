!===============================================================================!
!-------------------------------------------------------------------------------!
!            MODULE FOR THE GLOBAL VARIABLES OF THE OVERALL PROBLEM             !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Global_variables

!===============================================================================!
!===============================================================================!
!===============================================================================!
!=========================Configuration variables===============================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!

	!==================================================================!
    ! Mesh 
    ! (0) 3D Matlab Mesh
    ! (1) 3D Polycrystalline 
    ! (2) 3D General mesh .obj
	CHARACTER(LEN=300) :: Mesh_file 
	CHARACTER(LEN=300) :: fileplace_mesh
	REAL(8) :: scale_size_1 = 1.d0 !POINTS*scale_size_1    

	!==================================================================!
    ! Material: (0) Only one material
	! Material: (1) Different material for each domain
    INTEGER :: n_Materials = 0           
	CHARACTER(LEN=300) :: Material_coefficients_file !-> (0)
	CHARACTER(LEN=300) :: Material_coefficients_database !-> (1)
	CHARACTER(LEN=300) :: fileplace_material
	REAL(8) :: scale_prop_mat = 1.d0 ! C*scale_prop_mat
	REAL(8) :: density 
	!==================================================================!
	! BCs: (0) automatic BCs for a box
    ! BCs: (1) Input file of customize BCs
	REAL(8) :: box_face_1(6)=-1.d0,box_face_2(6)=-1.d0
    REAL(8) :: box_face_3(6)=-1.d0,box_face_4(6)=-1.d0
	REAL(8) :: box_face_5(6)=-1.d0,box_face_6(6)=-1.d0
	REAL(8), ALLOCATABLE :: T_D_BC(:)
	! If (1) 
	INTEGER :: BC_groups=0
	CHARACTER(LEN=300) :: fileplace_BCs, load_profile='Heaviside'
	REAL(8) :: omega=0.d0, phase=0.d0;


	INTEGER :: IDcomp=30
	!==================================================================!
	! Number of realizations
    INTEGER :: N_r = 1 

    ! Failure analysis
    INTEGER :: Crack = 0
    INTEGER :: Failure = 0

    ! Homogenization
    INTEGER :: H_fields = 0

	!==================================================================!
	! Output solutions
    CHARACTER(LEN=300):: fileplace_results
	INTEGER :: Working_memory = 13, Out_of_core=1
	REAL(8) :: scale_results=1;
	!==================================================================!
    ! Static analysis
    INTEGER :: Quasi_static=0
    INTEGER :: Static_steps=0 ! Number of load increments

	!==================================================================!
    ! Body forces [x, y, z] (N/mm^3) MPa (0) No, (1) Yes
	INTEGER :: Bodyforces = 0, type_bf = 0
	REAL(8) :: BodyForce(3)
	CHARACTER(LEN=300):: fileplace_Bodyforce, Bodyforce_file 
	REAL(8), ALLOCATABLE :: Body_Forces(:,:)
	!==================================================================!
    ! Dynamic trasient analysis
	INTEGER :: Transient = 0
	REAL(8) :: Dt ! Time step
	INTEGER :: Time_steps = 0 ! Number of time steps
	REAL(8) :: T_t ! Total time

	!==================================================================!
    ! Position of physical nodes
    REAL(8) :: lambda = 0.155d0
    !==================================================================!
    ! Non-singular integration points 
    INTEGER :: npoints_rule_NSI = 13
    !Singular integration points 
    INTEGER :: npoints_rule_SI = 13
    ! Type of integration: 1 -> Entire triangle
    !                      4 -> 4 triangles
    !                      8 -> 8 triangles
    !	                  16 -> 16 triangles
    INTEGER :: Int_type = 4
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
!===============================================================================!
    
!-------------------------------------------------------------------------------!
! Variables for INPUT DATA
!-------------------------------------------------------------------------------! 
    INTEGER, ALLOCATABLE :: FACES(:,:), ELEM(:,:)!, SEGMENTS(:,:)
    INTEGER, ALLOCATABLE :: El_reg(:,:), Mat_reg(:,:)
    REAL(8), ALLOCATABLE :: POINTS(:,:), NORMAL_VECTORS(:,:), Volumes(:)
    
    INTEGER :: nreg
    REAL(8) :: BC_Face(6,6)
	REAL(8) :: dmin ! Minimum distance of overall structure
	! Type of boundary condition: constrained specific points
    INTEGER :: Config_BC(1) = (/0/)
!-------------------------------------------------------------------------------!
! Variables for Fourier Coefficients
!-------------------------------------------------------------------------------!
    
    REAL(KIND=8) :: Pi = 3.141592653589793d0
    REAL(KIND=8) :: error = 1e-6
    
    REAL(KIND=8), ALLOCATABLE :: Angles(:,:)

    TYPE Cell2
        REAL(KIND=8) :: re(3,3), im(3,3)
        INTEGER :: pos_mat(2)
    END TYPE Cell2

    TYPE Cell8
        REAL(KIND=8) :: C_reg(6,6)
    END TYPE Cell8
    TYPE(Cell8), DIMENSION(:), ALLOCATABLE :: C_tot
   
!-------------------------------------------------------------------------------!
! Variables for Discretization
!-------------------------------------------------------------------------------!
    INTEGER, PARAMETER :: nnos_el=3
    REAL(KIND=8) :: Mat_G(3,3) 

    TYPE Cell4
        REAL(8), ALLOCATABLE :: PHY_NODES(:,:), GEO_NODES(:,:)
    END TYPE Cell4
    TYPE(Cell4), DIMENSION(:), ALLOCATABLE :: GEO_PHY_NODES

    TYPE Cell5
        INTEGER, ALLOCATABLE :: ELEM_PHY_NODES(:,:), ELEM_GEO_NODES(:,:)
    END TYPE Cell5
    TYPE(Cell5), DIMENSION(:), ALLOCATABLE :: GEO_PHY_ELEM

    INTEGER ::  contrain_nodes(4,2)
    INTEGER, ALLOCATABLE :: Subregions(:,:), Interfaces(:,:), Non_interfaces(:,:)
    INTEGER, ALLOCATABLE :: Non_int_plot(:,:)
    INTEGER, ALLOCATABLE :: ELEM_GEO_NODES_global(:,:), bc_new(:,:), Int_reg(:,:)
    REAL(8), ALLOCATABLE :: BC_elem(:,:), GEO_NODES_global(:,:), PHY_NODES_global(:,:)
    REAL(8), ALLOCATABLE :: bct_0(:), bct_1(:), bct_2(:), bcs_3(:)
    REAL(8), ALLOCATABLE :: bcd_0(:), bcd_1(:), bcd_2(:), bcsd_3(:)
	INTEGER, ALLOCATABLE :: nodes_conectivity(:,:), Nodes_BCs(:,:)

	TYPE Cell14
		INTEGER :: nel, type_bc(3)
		INTEGER, ALLOCATABLE :: el(:)
        REAL(8), ALLOCATABLE :: bc(:,:)
    END TYPE Cell14
    TYPE(Cell14), DIMENSION(:), ALLOCATABLE :: BCs(:)

 
!-------------------------------------------------------------------------------!
! Variables for Int_points
!-------------------------------------------------------------------------------!
    ! Non-singular integration type, 4 sub-triangles (Four_tri) or
    ! one triangle (One_tri)
    REAL(8) :: N_Geo_nodes2(9,1), dNda_i(2,3)    
    INTEGER :: N_int_points(3), npoints_Nsi
    REAL(8), ALLOCATABLE :: weights_ns(:,:)
    TYPE Cell3
        REAL(KIND=8), ALLOCATABLE :: weights_s(:,:)
    END TYPE Cell3
    TYPE(Cell3) :: Int_data(nnos_el)
    REAL(8) :: N_Phy_nodes(9,1)
!-------------------------------------------------------------------------------!
! Variables for Comp_H_G
!-------------------------------------------------------------------------------!    
    TYPE Cell6
        REAL(8), ALLOCATABLE :: H(:,:), G(:,:)
    END TYPE Cell6
    TYPE(Cell6), DIMENSION(:), ALLOCATABLE :: H_G
!-------------------------------------------------------------------------------!
! Variables for DRM
!-------------------------------------------------------------------------------! 
    
    TYPE Cell7
        REAL(8), ALLOCATABLE :: M(:,:)
    END TYPE Cell7
    TYPE(Cell7), DIMENSION(:), ALLOCATABLE :: U_T_E_M
!-------------------------------------------------------------------------------!
! Variables for General_Matrices
!-------------------------------------------------------------------------------! 
    REAL(8), ALLOCATABLE :: b(:,:), bhat(:,:), V(:,:), b_cte(:,:), b_bf(:,:)
    REAL(8), ALLOCATABLE :: b_step(:,:), G_geral(:,:), B_geral(:,:)
    !REAL(8), ALLOCATABLE :: acsr_aux(:)
    !INTEGER, ALLOCATABLE :: ia_aux(:), ja_aux(:)
!-------------------------------------------------------------------------------!
! Variables for General_Matrices SPARCE and General_Matrices
!-------------------------------------------------------------------------------! 
	INTEGER ::  N_total
	INTEGER(8) :: nzA, nzG, nzV

	TYPE Cell13	
	INTEGER, ALLOCATABLE :: col_A(:), col_G(:), col_V(:)
	INTEGER, ALLOCATABLE :: row_A(:)
	REAL(8), ALLOCATABLE :: At(:), G_t(:), V_t(:)
    END TYPE Cell13
    TYPE(Cell13), DIMENSION(:), ALLOCATABLE :: S_E  

	INTEGER, ALLOCATABLE :: rows_G(:), rows_V(:)
	INTEGER(8), ALLOCATABLE :: scounts_G(:), scounts_V(:)
	INTEGER(8), ALLOCATABLE :: scounts_A(:)
	INTEGER(8), ALLOCATABLE :: displs_G(:), displs_V(:)
 

!-------------------------------------------------------------------------------!
! Variables for Module_Solver
!-------------------------------------------------------------------------------! 
    REAL(8), ALLOCATABLE :: x_BEM(:,:)
!-------------------------------------------------------------------------------!
! Variables for Response_D_T
!-------------------------------------------------------------------------------! 
    REAL(8), ALLOCATABLE :: Displacement(:,:), Traction(:,:), u_t(:,:)
!-------------------------------------------------------------------------------!
! Variables for Response_S_S
!-------------------------------------------------------------------------------! 
    REAL(8), ALLOCATABLE :: Strain(:,:), Stress(:,:)
    REAL(8), ALLOCATABLE :: NORMAL_VEC_Def(:,:)
!-------------------------------------------------------------------------------!
! Variables for validation
!-------------------------------------------------------------------------------! 
    REAL(8), ALLOCATABLE :: Disp_p(:), Stress_p(:,:), Strain_p(:,:), Stress_h(:)
    REAL(8), ALLOCATABLE :: T_n_max(:), T_n_min(:), S_n_max(:), Strain_h(:)

    
!================================================================================
! FAILURE ANALYSIS

    REAL(8) :: xmin_f, xmax_f, ymin_f, ymax_f, zmin_f, zmax_f  
        
    INTEGER :: n_TSL, Iteration
    TYPE Cell9
        INTEGER, ALLOCATABLE :: GBs(:)
	REAL(8) :: E0c(5), theta_z_x(2), R(3,3)
	REAL(8) :: R_inv(3,3)
    END TYPE Cell9
    TYPE(Cell9), DIMENSION(:), ALLOCATABLE :: TSL

    TYPE Cell10
	INTEGER :: alloc
	REAL(8) :: Es1(3), Es2(3), En(3), delta(9), Criterion(3)
	REAL(8), ALLOCATABLE :: H_cohe(:,:), G_cohe(:,:), M_cohe(:,:), M_cohe_G(:,:)
	REAL(8), ALLOCATABLE :: H_local(:,:), G_local(:,:), M_local(:,:)
    END TYPE Cell10
    TYPE(Cell10), DIMENSION(:), ALLOCATABLE :: Cohesive 

	TYPE Cell11
	REAL(8) :: S_l(3,9), e_l(3,9)
    END TYPE Cell11
    TYPE(Cell11), DIMENSION(:), ALLOCATABLE :: Stress_strain_local

	REAL(8) :: eta_load, load, step_load
	INTEGER :: contact, Reduce_A

	INTEGER, ALLOCATABLE :: Damage_zones(:)
	REAL(8), ALLOCATABLE :: Damage_coeff(:,:) 
	INTEGER :: Damage_elements

	INTEGER :: adhesive
	REAL(8) :: current_load

	! To read coordinates of pre-cracks
    CHARACTER(*),PARAMETER :: fileplace7 = "Cracks/"
    CHARACTER(*),PARAMETER :: fileplace8 = "Results/Paraview_plot/failure_def/"
    CHARACTER(*),PARAMETER :: fileplace9 = "Results/Paraview_plot/failure/"
    CHARACTER(LEN=200)::Test_name


	CHARACTER(*),PARAMETER :: fileplace_results_matlab = "Results/Matlab_plot/Historic/"
	CHARACTER(LEN=100) :: Results_file
!===============================================================================!
END MODULE Global_variables
