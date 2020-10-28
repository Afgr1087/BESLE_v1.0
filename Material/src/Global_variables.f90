!===============================================================================!
!-------------------------------------------------------------------------------!
!            MODULE FOR THE GLOBAL VARIABLES OF THE OVERALL PROBLEM             !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Global_variables
!-------------------------------------------------------------------------------!

    CHARACTER(LEN=300):: fileplace
    CHARACTER(LEN=300)::file_name
    INTEGER :: n_materials=1, step
    INTEGER :: z_x_z = 0
    INTEGER :: x_y_z = 1
    CHARACTER(LEN=300) :: Lattice, Material
    REAL(KIND=8) :: C(6,6), C_original(6,6) ! Compliance tensor
    REAL(KIND=8) :: C11=0.d0, C12=0.d0, C44=0.d0
    REAL(KIND=8) :: C15=0.d0, C16=0.d0, C14=0.d0, C33=0.d0, C13=0.d0
    REAL(KIND=8) :: C21=0.d0, C22=0.d0, C23=0.d0, C24=0.d0, C25=0.d0, C26=0.d0
    REAL(KIND=8) :: C31=0.d0, C32=0.d0, C34=0.d0, C35=0.d0, C36=0.d0
    REAL(KIND=8) :: C41=0.d0, C42=0.d0, C43=0.d0, C45=0.d0, C46=0.d0
    REAL(KIND=8) :: C51=0.d0, C52=0.d0, C53=0.d0, C54=0.d0, C55=0.d0, C56=0.d0
    REAL(KIND=8) :: C61=0.d0, C62=0.d0, C63=0.d0, C64=0.d0, C65=0.d0, C66=0.d0
    REAL(KIND=8) :: E, nu ! Case of isotropic material
    REAL(KIND=8) :: Pi = 3.1415926535897932d0
    REAL(KIND=8) :: phi_1, phi, phi_2
	REAL(KIND=8) :: theta_x=0.d0, theta_y=0.d0, theta_z=0.d0
    REAL(KIND=8) :: error = 1e-6, max_val
    TYPE Cell
        REAL(KIND=8) :: comp(3,3)
    END TYPE Cell
    TYPE(Cell), DIMENSION(:,:), ALLOCATABLE :: Rt_mat, Rc_mat
    TYPE(Cell), DIMENSION(:,:), ALLOCATABLE :: It_mat, Ic_mat
    TYPE(Cell), DIMENSION(:), ALLOCATABLE :: R0m_vet, Rm0_vet
    TYPE(Cell), DIMENSION(:), ALLOCATABLE :: I0m_vet, Im0_vet
    REAL(8), ALLOCATABLE :: Rt_matrix(:,:), Rc_matrix(:,:)
    REAL(8), ALLOCATABLE :: It_matrix(:,:), Ic_matrix(:,:)
    REAL(8), ALLOCATABLE :: R0m_vector(:,:), Rm0_vector(:,:)
    REAL(8), ALLOCATABLE :: I0m_vector(:,:), Im0_vector(:,:)
    REAL(KIND=8) :: R00_matrix(3,3)

	REAL(8), ALLOCATABLE :: C_tensors(:,:), E_v_constants(:,:)

!===============================================================================!
END MODULE Global_variables
