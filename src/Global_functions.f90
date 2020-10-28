!===============================================================================!
!-------------------------------------------------------------------------------!
!                         MODULE OF GLOBAL SUBROUTINES                          !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Global_functions
!-------------------------------------------------------------------------------!
USE Global_variables
!USE mpi
USE omp_lib
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
    SUBROUTINE Gauss_Legendre(x1,x2,n,x,w)

        INTEGER :: m, n, i, j, cont
        REAL(8) :: eps, xm, xl, x1, x2, x(:), w(:), z, p1, p2, p3, pp, z1

        eps = 3e-15
        m = NINT((n+1.d0)/2.d0)
        xm = 0.5d0*(x2+x1)
        xl = 0.5d0*(x2-x1)
        
        DO i=1,m
            z = DCOS(Pi*(i-0.25d0)/(n+0.5d0))
            cont = 0
            DO WHILE (1 == 1)
                cont = cont + 1
                p1 = 1.d0
                p2 = 0.d0
                DO j=1,n
                    p3 = p2
                    p2 = p1
                    p1 = ((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
                END DO
                pp=n*(z*p1-p2)/(z*z-1.d0)
                z1=z;
                z=z1-p1/pp
                IF (ABS(z-z1).LT.eps) THEN
                    EXIT
                END IF
                IF (cont == 20) EXIT
            END DO
            x(i)=xm-xl*z;
            x(n+1-i)=xm+xl*z;
            w(i)=2.d0*xl/((1-z*z)*pp*pp);
            w(n+1-i)=w(i);
        END DO

    END SUBROUTINE Gauss_Legendre
!===============================================================================!
    SUBROUTINE Shape_functions_cont(qsi,eta,tp,N)

        !--------------------------------------------------------------------------
        ! Input data:
        ! qsi, eta - Gauss points where the functions are calculated.
        ! Ouput data:
        ! [N] - Shape functions
        !--------------------------------------------------------------------------

        INTEGER :: tp
        REAL(8) :: qsi, eta, N(:)

        IF (tp.EQ.1) THEN
            N(1) = qsi
            N(2) = eta
            N(3) = 1.d0-qsi-eta
        ELSEIF(tp.EQ.2) THEN
            N(1) = qsi*(2.d0*qsi-1.d0)
            N(2) = eta*(2.d0*eta-1.d0)
            N(3) = (1.d0-qsi-eta)*(1.d0-2.d0*qsi-2.d0*eta)
            N(4) = 4.d0*qsi*eta
            N(5) = 4.d0*eta*(1.d0-qsi-eta)
            N(6) = 4.d0*qsi*(1.d0-qsi-eta)
        END IF

    END SUBROUTINE Shape_functions_cont
!===============================================================================!
SUBROUTINE Rule_points_weights(npoints_rule,rp,rw)  
    
    INTEGER :: npoints_rule
    REAL(8):: rp(:,:), rw(:)
    
    IF (npoints_rule.EQ.7) THEN
       
        rp(1,:) = (/0.3333333333d0,  0.3333333333d0,  0.3333333333d0/)
        rp(2,:) = (/0.7974269853d0,  0.1012865073d0,  0.1012865073d0/)
        rp(3,:) = (/0.1012865073d0,  0.7974269853d0,  0.1012865073d0/)
        rp(4,:) = (/0.1012865073d0,  0.1012865073d0,  0.7974269853d0/)
        rp(5,:) = (/0.0597158717d0,  0.4701420641d0,  0.4701420641d0/)
        rp(6,:) = (/0.4701420641d0,  0.0597158717d0,  0.4701420641d0/)
        rp(7,:) = (/0.4701420641d0,  0.4701420641d0,  0.0597158717d0/)

        rw(1) = 0.255d0; rw(2) = 0.1259391805d0; rw(3) = 0.1259391805d0
        rw(4) = 0.1259391805d0; rw(5) = 0.1329341527d0; rw(6) = 0.1329341527d0
        rw(7) = 0.1329341527d0
    ELSEIF (npoints_rule.EQ.13) THEN
        
        rp(1,:) = (/0.0651301029d0,  0.0651301029d0,  0.8697397942d0/)
        rp(2,:) = (/0.8697297941d0,  0.0651301029d0,  0.0651401030d0/)
        rp(3,:) = (/0.0651301029d0,  0.8697297941d0,  0.0651401030d0/)
        rp(4,:) = (/0.3128654960d0,  0.0486903154d0,  0.6384441886d0/)
        rp(5,:) = (/0.6384441885d0,  0.3128654960d0,  0.0486903155d0/)
        rp(6,:) = (/0.0486903154d0,  0.6384441885d0,  0.3128654961d0/)
        rp(7,:) = (/0.6384441885d0,  0.0486903154d0,  0.3128654961d0/)
        rp(8,:) = (/0.3128654960d0,  0.6384441885d0,  0.0486903155d0/)
        rp(9,:) = (/0.0486903154d0,  0.3128654960d0,  0.6384441886d0/)
        rp(10,:) = (/0.2603459660d0,  0.2603459660d0,  0.4793080680d0/)
        rp(11,:) = (/0.4793080678d0,  0.2603459660d0,  0.2603459662d0/)
        rp(12,:) = (/0.2603459660d0,  0.4793080678d0,  0.2603459662d0/)
        rp(13,:) = (/0.3333333333d0,  0.3333333333d0,  0.3333333333d0/)

        rw(1) = 0.0533472356d0; rw(2) = 0.0533472356d0; rw(3) = 0.0533472356d0
        rw(4) = 0.0771137608d0; rw(5) = 0.0771137608d0; rw(6) = 0.0771137608d0
        rw(7) = 0.0771137608d0; rw(8) = 0.0771137608d0; rw(9) = 0.0771137608d0
        rw(10) = 0.1756152574d0; rw(11) = 0.1756152574d0; rw(12) = 0.1756152574d0
        rw(13) = -0.1495700444d0
    END IF

END SUBROUTINE Rule_points_weights
!===============================================================================!
SUBROUTINE Shape_functions(qsi,eta,N)  
    
    REAL(8) :: qsi, eta, N(:), h_1, h_2, h_3
    REAL(8) :: N_1, N_2, N_3 
    
    ! shape functions of original continuous triangle element
    
    h_1 = qsi
    h_2 = eta
    h_3 = 1.d0-qsi-eta

    ! shape functions of discontinuous triangle element
    N_1 = Mat_G(1,1)*h_1+Mat_G(2,1)*h_2+Mat_G(3,1)*h_3
    N_2 = Mat_G(1,2)*h_1+Mat_G(2,2)*h_2+Mat_G(3,2)*h_3
    N_3 = Mat_G(1,3)*h_1+Mat_G(2,3)*h_2+Mat_G(3,3)*h_3
    
    ! Shape functions
    N = (/N_1,N_2,N_3/)
        
END SUBROUTINE Shape_functions
!===============================================================================!
SUBROUTINE Dev_shape_functions(dNdqsi,dNdeta)

    REAL(8) :: dNdqsi(3), dNdeta(3)
    REAL(8) :: dh1dqsi, dh2dqsi, dh3dqsi
    REAL(8) :: dh1deta, dh2deta, dh3deta
    REAL(8) :: dN1dqsi, dN2dqsi, dN3dqsi
    REAL(8) :: dN1deta, dN2deta, dN3deta

    ! derivative shape functions of original continuous triangle element
    
    dh1dqsi = 1.d0; dh2dqsi = 0.d0; dh3dqsi = -1.d0

    dh1deta = 0.d0; dh2deta = 1.d0; dh3deta = -1.d0

    ! derivatives shape functions of discontinuous triangle element
    dN1dqsi = Mat_G(1,1)*dh1dqsi+Mat_G(2,1)*dh2dqsi+Mat_G(3,1)*dh3dqsi
    dN2dqsi = Mat_G(1,2)*dh1dqsi+Mat_G(2,2)*dh2dqsi+Mat_G(3,2)*dh3dqsi
    dN3dqsi = Mat_G(1,3)*dh1dqsi+Mat_G(2,3)*dh2dqsi+Mat_G(3,3)*dh3dqsi
        
    dN1deta = Mat_G(1,1)*dh1deta+Mat_G(2,1)*dh2deta+Mat_G(3,1)*dh3deta
    dN2deta = Mat_G(1,2)*dh1deta+Mat_G(2,2)*dh2deta+Mat_G(3,2)*dh3deta
    dN3deta = Mat_G(1,3)*dh1deta+Mat_G(2,3)*dh2deta+Mat_G(3,3)*dh3deta
        
    ! derivatives of shape functions
    dNdqsi = (/dN1dqsi,dN2dqsi,dN3dqsi/)
    dNdeta = (/dN1deta,dN2deta,dN3deta/)

END SUBROUTINE Dev_shape_functions
!===============================================================================!
SUBROUTINE Dev_shape_funct_cont(qsi,eta,tp,dNdqsi,dNdeta)
    INTEGER :: tp
    REAL(8) :: qsi, eta, dNdqsi(:), dNdeta(:)
    REAL(8) :: dh1dqsi, dh2dqsi, dh3dqsi, dh4dqsi, dh5dqsi, dh6dqsi
    REAL(8) :: dh1deta, dh2deta, dh3deta, dh4deta, dh5deta, dh6deta

    IF (tp.EQ.1) THEN ! Three-node triangle continuous element
        dh1dqsi = 1.d0; dh2dqsi = 0.d0; dh3dqsi = -1.d0
        dh1deta = 0.d0; dh2deta = 1.d0; dh3deta = -1.d0

        dNdqsi = (/dh1dqsi,dh2dqsi,dh3dqsi/) 
        dNdeta = (/dh1deta,dh2deta,dh3deta/)
    ELSEIF (tp.EQ.2) THEN ! Six-node triangle continuous elment
        dh1dqsi = 4.d0*qsi-1.d0;          dh2dqsi = 0.d0
        dh3dqsi = 4.d0*qsi+4.d0*eta-3.d0; dh4dqsi = 4.d0*eta
        dh5dqsi = -4.d0*eta;              dh6dqsi = 4.d0*(1.d0-eta)

        dh1deta = 0.d0;                   dh2deta  = 4.d0*eta-1.d0
        dh3deta = 4.d0*qsi+4.d0*eta-3.d0; dh4deta  = 4.d0*qsi
        dh5deta = 4.d0*(1.d0-qsi);        dh6deta  = -4.d0*qsi

        dNdqsi = (/dh1dqsi,dh2dqsi,dh3dqsi,dh4dqsi,dh5dqsi,dh6dqsi/)  
        dNdeta = (/dh1deta,dh2deta,dh3deta,dh4deta,dh5deta,dh6deta/)     
    END IF

END SUBROUTINE Dev_shape_funct_cont
!===============================================================================!
SUBROUTINE Calc_jac_tri3(X1,X2,X3,dNdqsi,dNdeta,J)

    REAL(8) :: X1(:), X2(:), X3(:), dNdqsi(:), dNdeta(:), J
    REAL(8) :: dxdqsi, dydqsi, dxdeta, dydeta

    dxdqsi = X1(1)*dNdqsi(1) + X2(1)*dNdqsi(2) + X3(1)*dNdqsi(3)
    dydqsi = X1(2)*dNdqsi(1) + X2(2)*dNdqsi(2) + X3(2)*dNdqsi(3)

    dxdeta = X1(1)*dNdeta(1) + X2(1)*dNdeta(2) + X3(1)*dNdeta(3)
    dydeta = X1(2)*dNdeta(1) + X2(2)*dNdeta(2) + X3(2)*dNdeta(3)

    J = 0.5d0*DSQRT((dxdqsi*dydeta - dydqsi*dxdeta)**2);

END SUBROUTINE Calc_jac_tri3
!===============================================================================!
SUBROUTINE Shape_functions_quad4(qsi,eta,N)

    REAL(8) :: qsi,eta, N(:)
    
        N(1) = (1.d0/4.d0)*(1.d0 - qsi)*(1.d0 - eta)
        N(2) = (1.d0/4.d0)*(1.d0 + qsi)*(1.d0 - eta)
        N(3) = (1.d0/4.d0)*(1.d0 + qsi)*(1.d0 + eta)
        N(4) = (1.d0/4.d0)*(1.d0 - qsi)*(1.d0 + eta)

END SUBROUTINE Shape_functions_quad4
!===============================================================================!
SUBROUTINE Dev_shape_functions_quad4(qsi,eta,dNdqsi,dNdeta)

    REAL(8) :: qsi,eta, dNdqsi(:), dNdeta(:)
    
    dNdqsi(1) = (1.d0/4.d0)*(-(1.d0-eta))
    dNdqsi(2) = (1.d0/4.d0)*(1.d0 - eta)
    dNdqsi(3) = (1.d0/4.d0)*(1.d0 + eta)
    dNdqsi(4) = (1.d0/4.d0)*(-(1.d0+eta))

    dNdeta(1) = (1.d0/4.d0)*(-(1.d0 - qsi))
    dNdeta(2) = (1.d0/4.d0)*(-(1.d0 + qsi))
    dNdeta(3) = (1.d0/4.d0)*(1.d0 + qsi)
    dNdeta(4) = (1.d0/4.d0)*(1.d0 - qsi)

END SUBROUTINE Dev_shape_functions_quad4
!===============================================================================!
SUBROUTINE Calc_jac_quad4(X1,X2,X3,X4,dNdqsi,dNdeta,J)

    REAL(8) :: X1(:),X2(:),X3(:),X4(:),dNdqsi(:),dNdeta(:), J
    REAL(8) :: dxdqsi, dydqsi, dxdeta, dydeta

    dxdqsi = X1(1)*dNdqsi(1) + X2(1)*dNdqsi(2) + X3(1)*dNdqsi(3) + X4(1)*dNdqsi(4); 
    dydqsi = X1(2)*dNdqsi(1) + X2(2)*dNdqsi(2) + X3(2)*dNdqsi(3) + X4(2)*dNdqsi(4);

    dxdeta = X1(1)*dNdeta(1) + X2(1)*dNdeta(2) + X3(1)*dNdeta(3) + X4(1)*dNdeta(4); 
    dydeta = X1(2)*dNdeta(1) + X2(2)*dNdeta(2) + X3(2)*dNdeta(3) + X4(2)*dNdeta(4); 

    J = DSQRT((dxdqsi*dydeta - dydqsi*dxdeta)**2)

END SUBROUTINE Calc_jac_quad4
!===============================================================================!
SUBROUTINE find_vector(Vector,val,pos)
    
    
    INTEGER :: Vector(nnos_el), i, pos, val
    pos = 0
    DO i = 1,nnos_el
        IF (val.EQ.Vector(i)) THEN
            pos = i
        END IF
    END DO

END SUBROUTINE find_vector
!===============================================================================!
SUBROUTINE find_vector_g(Vector,size_v,val,pos)
    
    
    INTEGER :: size_v, i, pos, val
    INTEGER, ALLOCATABLE :: Vector(:)
    pos = 0
    DO i = 1,size_v
        IF (val.EQ.Vector(i)) THEN
            pos = i
        END IF
    END DO

END SUBROUTINE find_vector_g
!===============================================================================!
SUBROUTINE spherical_coordinates(xsp,xfp,theta,phi,r)

    REAL(8) :: xsp(:), xfp(:), x(3), theta, phi, r, x1, x2, x3

    x = xfp - xsp
    x1 = x(1); x2 = x(2); x3 = x(3)
    CALL Atan_2(x1,x2,theta)
    CALL Atan_2(x3,DSQRT(x1**2+x2**2),phi)
    r = DSQRT(x1**2 + x2**2 + x3**2)
    
END SUBROUTINE spherical_coordinates
!===============================================================================!
SUBROUTINE Atan_2(x,y,ang)

    REAL(8) :: x, y, x1, x2, ang, Tol

    Tol = 1e-10

    IF ((x.GT.0).AND.(y.GT.0)) THEN ! First quadrant
        x1 = ABS(x); x2 = ABS(y)
        ang = DATAN(x2/x1)
    END IF
    IF ((x.LT.0).AND.(y.GT.0)) THEN ! Second quadrant
        x1 = ABS(x); x2 = ABS(y)
        ang = Pi - ATAN(x2/x1)
    END IF
    IF ((x.LT.0).AND.(y.LT.0)) THEN ! Third quadrant
        x1 = ABS(x); x2 = ABS(y)
        ang = -Pi + DATAN(x2/x1)
    END IF
    IF ((x.GT.0).AND.(y.LT.0)) THEN ! Fourth quadrant
        x1 = ABS(x); x2 = ABS(y)
        ang = -DATAN(x2/x1)
    END IF
    IF ((x.GT.0).AND.(ABS(y).LT.Tol)) THEN ! in positive x-axis
        ang = 0.d0
    END IF
    IF ((x.LT.0).AND.(ABS(y).LT.Tol)) THEN ! in negative x-axis
        ang = Pi
    END IF
    IF ((ABS(x).LT.Tol).AND.(y.GT.0)) THEN ! in positive y-axis
        ang = Pi/2
    END IF
    IF ((ABS(x).LT.Tol).AND.(y.LT.0)) THEN ! in negative y-axis
        ang = 3.d0*Pi/2
    END IF
    
END SUBROUTINE Atan_2
!===============================================================================!
SUBROUTINE Inv_matrix(Matrix,Matrix2)

    INTEGER :: M, N, LDA, INFO, LWORK
    INTEGER, ALLOCATABLE :: IPIV(:)
    REAL(8), ALLOCATABLE :: Matrix(:,:), Matrix2(:,:), WORK(:)

    ! compute an LU factorization of a general M-by-N matrix F_reg
    M = SIZE(Matrix,1)
    N = SIZE(Matrix,2)
    LDA = M
    ALLOCATE(IPIV(M))
    IPIV = 0
    INFO = 0
    CALL DGETRF(M, N, Matrix, LDA, IPIV, INFO)
    
    ! compute the inverse of a matrix using the LU	factorization
    LWORK = N
    ALLOCATE(WORK(LWORK))
    WORK = 0.d0
    CALL DGETRI(N, Matrix, LDA, IPIV, WORK, LWORK, INFO)

    Matrix2 = Matrix

    DEALLOCATE(IPIV, WORK)

END SUBROUTINE Inv_matrix
!===============================================================================!
SUBROUTINE Matmul_parallel(nt1,me1,Matrix_A,Matrix_B,Matrix_C,NEWCOMM1)
    
	include 'mpif.h'
    INTEGER :: nt1, me1
    REAL(8), ALLOCATABLE :: Matrix_A(:,:) ,Matrix_B(:,:) ,Matrix_C(:,:)  
    
    INTEGER, ALLOCATABLE :: scounts(:), displs(:), scounts_B(:), displs_B(:)
    INTEGER :: factor2, factor4, m, n, sub_nt, val, NEWCOMM1, size_mat_A
    REAL(8) :: factor1, factor3, nelem_real

    INTEGER :: root=0, mpierr, World_group, New_group, NEWCOMM, sub_me, i
    INTEGER, ALLOCATABLE :: ranks(:)    

    INTEGER :: NFB, NCB, row_i_A, row_f_A, nelem

    REAL(8), ALLOCATABLE :: mat_A(:,:), sub_mat_B(:,:), sub_mat_C(:,:), mat_C(:,:)

    CALL MPI_BARRIER(NEWCOMM1,mpierr);

    ! Division of Matrix A and B
    IF (me1.EQ.0) THEN
        
		! Matrix B
		nelem = SIZE(Matrix_B,2) 
		ALLOCATE(scounts_B(nt1), displs_B(nt1))
		CALL Par_Div_32(nt1, nelem, scounts_B, displs_B)
		
        nelem = SIZE(Matrix_A,1) 

		ALLOCATE(scounts(nt1), displs(nt1))
		CALL Par_Div_32(nt1, nelem, scounts, displs)
		
        ! ---------------------------------------------------------------------------
		sub_nt = nt1
    END IF
	
    !---------------------------------------------------------------------------------
    ! New group according to the size of Matrix_A

    CALL MPI_Bcast(sub_nt,1,MPI_INTEGER,root,NEWCOMM1,mpierr)
    
    CALL MPI_COMM_GROUP (NEWCOMM1, World_group, mpierr)
     
    ALLOCATE(ranks(sub_nt))
    ranks = (/(i,i=0,sub_nt-1)/)
      
    CALL MPI_GROUP_INCL (World_group, sub_nt, ranks, New_group, mpierr)
    CALL MPI_COMM_CREATE (NEWCOMM1, New_group, NEWCOMM, mpierr)
        
    sub_me = -1; sub_nt = -1
    
    IF (MPI_COMM_NULL .NE. NEWCOMM) THEN
        CALL MPI_Comm_rank(NEWCOMM, sub_me, mpierr);
        CALL MPI_Comm_size(NEWCOMM, sub_nt, mpierr);
        
        IF (sub_me.NE.0) THEN
            ALLOCATE(scounts(sub_nt), displs(sub_nt))
			ALLOCATE(scounts_B(sub_nt), displs_B(sub_nt))
        END IF
        CALL MPI_Bcast(scounts,sub_nt,MPI_INTEGER,root,NEWCOMM,mpierr)
        CALL MPI_Bcast(displs,sub_nt,MPI_INTEGER,root,NEWCOMM,mpierr)
		CALL MPI_Bcast(scounts_B,sub_nt,MPI_INTEGER,root,NEWCOMM,mpierr)
        CALL MPI_Bcast(displs_B,sub_nt,MPI_INTEGER,root,NEWCOMM,mpierr)
		

        ! Distribution of mat_A to all blocks
        IF (sub_me.EQ.0) THEN
            NFB = SIZE(Matrix_B,1); NCB = SIZE(Matrix_B,2)
			
        END IF

        CALL MPI_Bcast(NFB,1,MPI_INTEGER,root,NEWCOMM,mpierr)
        CALL MPI_Bcast(NCB,1,MPI_INTEGER,root,NEWCOMM,mpierr)
		
		
        CALL Distr_mat2(sub_nt,sub_me,NEWCOMM,scounts_B,displs_B,NFB,NCB, &
                                                              Matrix_B,sub_mat_B)

        n = 0
1       n = n + 1 ! Loop for divisions of Matrix_A
		
            ! Select mat_A (row matrix) from Matrix_A
            row_i_A = displs(n) + 1
            row_f_A = displs(n) + scounts(n) 
            
            ALLOCATE(mat_A(scounts(n),NFB))
            IF (sub_me.EQ.0) THEN
                mat_A = Matrix_A(row_i_A:row_f_A,:)
				size_mat_A = SIZE(mat_A)
            END IF
			CALL MPI_Bcast(size_mat_A,1,MPI_INTEGER,root,NEWCOMM,mpierr)
            CALL MPI_Bcast(mat_A,size_mat_A,MPI_REAL8,root,NEWCOMM,mpierr)

            ! Result = mat_A * sub_mat_B
            ALLOCATE(sub_mat_C(scounts(n),scounts_B(sub_me+1)))

            sub_mat_C = MATMUL(mat_A,sub_mat_B)

            ! Join results in the host: mat_C
            IF (sub_me.EQ.0) THEN
                ALLOCATE(mat_C(scounts(n),displs_B(sub_nt)+scounts_B(sub_nt)))
				
            END IF
			
            CALL Join_mat2(sub_nt,sub_me,NEWCOMM,scounts_B,displs_B,sub_mat_C,mat_C)
            
			

            IF (sub_me.EQ.0) THEN
				
                Matrix_C(row_i_A:row_f_A,:) = mat_C
        		
                DEALLOCATE(mat_C)

            END IF
            DEALLOCATE(mat_A,sub_mat_C)

        IF (n.LE.sub_nt-1) THEN
            GOTO 1 ! Go to the next element
        ELSE
            DEALLOCATE(scounts,displs,scounts_B,displs_B)
            GOTO 2 ! All processors finish the task
        END IF

2       DEALLOCATE(sub_mat_B)
        
    END IF
    !---------------------------------------------------------------------------------

    IF (MPI_COMM_NULL .NE. NEWCOMM) THEN
        CALL MPI_COMM_FREE(NEWCOMM, mpierr)
    END IF
        
    CALL MPI_Group_free(World_group, mpierr);
    CALL MPI_Group_free(New_group, mpierr);

    DEALLOCATE(ranks)

END SUBROUTINE Matmul_parallel
!===============================================================================!
SUBROUTINE Phy_to_Geo(field_g)

	INTEGER :: i, nelem, j, var
	REAL(8) :: N(3), phi_x(3), phi_y(3), phi_z(3), geo_x(3), geo_y(3), geo_z(3)
    REAL(8), ALLOCATABLE :: field_g(:,:)	

	nelem = SIZE(field_g,1)

	DO i=1,nelem/3

        phi_x = field_g(3*i-2:3*i,1)
        phi_y = field_g(3*i-2:3*i,2)
        phi_z = field_g(3*i-2:3*i,3)
        
        DO j=1,nnos_el

            N = N_Geo_nodes2(3*j-2:3*j,1)
            
            geo_x(j) = phi_x(1)*N(1) + phi_x(2)*N(2) + phi_x(3)*N(3)
            geo_y(j) = phi_y(1)*N(1) + phi_y(2)*N(2) + phi_y(3)*N(3) 
            geo_z(j) = phi_z(1)*N(1) + phi_z(2)*N(2) + phi_z(3)*N(3)

        END DO

        field_g(3*i-2:3*i,1) = geo_x
        field_g(3*i-2:3*i,2) = geo_y
        field_g(3*i-2:3*i,3) = geo_z
            
    END DO


END SUBROUTINE Phy_to_Geo
!===============================================================================!
SUBROUTINE Par_Div(nt_sub, n_A, scounts, displs)

	INTEGER :: i
	INTEGER(8) :: m, n, n_A, val, nt_sub
	INTEGER(8), ALLOCATABLE :: scounts(:)
	INTEGER(8), ALLOCATABLE :: displs(:)
	
	m = n_A/nt_sub;
    
    IF (MOD(n_A,nt_sub).EQ.0) THEN
    	scounts = n_A/(nt_sub);
    ELSE    	
        n = n_A - m*(nt_sub-1);
		scounts(1:nt_sub-1) = m
		
		! The last thread have more than one element compated to the other threads
		IF (n-m.GT.1) THEN 
			! the additional number of elements are greater than thenumber of 
			! remaining threads
							
			DO i=1,n-m
				scounts(i) = scounts(i) + 1
			END DO
			n = m

		END IF

		scounts(nt_sub) = n
         
    END IF

	val = 0
    displs(1) = 0

    DO i=1,nt_sub-1
    val = val + scounts(i)
    displs(i+1) = val
    END DO

END SUBROUTINE Par_Div
!===============================================================================!
SUBROUTINE Par_Div_32(nt_sub, n_A, scounts, displs)

	INTEGER :: i, nt_sub
	INTEGER :: m, n, n_A, val
	INTEGER, ALLOCATABLE :: scounts(:)
	INTEGER, ALLOCATABLE :: displs(:)

	m = INT(n_A/nt_sub);
	   
    IF (MOD(n_A,nt_sub).EQ.0) THEN
    	scounts = n_A/(nt_sub);
    ELSE    	
        n = n_A - m*(nt_sub-1);
		scounts(1:nt_sub-1) = m
		
		! The last thread have more than one element compated to the other threads
		IF (n-m.GT.1) THEN 
			! the additional number of elements are greater than thenumber of 
			! remaining threads
							
			DO i=1,n-m
				scounts(i) = scounts(i) + 1
			END DO
			n = m

		END IF

		scounts(nt_sub) = n
         
    END IF

	val = 0
    displs(1) = 0
    DO i=1,nt_sub-1
    val = val + scounts(i)
    displs(i+1) = val
    END DO

END SUBROUTINE Par_Div_32
!===============================================================================!
SUBROUTINE Matmul_parallel2(nt,me,Matrix_A,Matrix_B,Matrix_C)
    include 'mpif.h'
    INTEGER :: nt, me
    REAL(8), ALLOCATABLE :: Matrix_A(:,:) ,Matrix_B(:,:) ,Matrix_C(:,:)  
    
    INTEGER, ALLOCATABLE :: scounts(:), displs(:)
    INTEGER :: factor2, factor4, m, n, sub_nt, val
    REAL(8) :: factor1, factor3, nelem_real

    INTEGER :: root=0, mpierr, World_group, New_group, NEWCOMM, sub_me, i
    INTEGER, ALLOCATABLE :: ranks(:)    

    INTEGER :: NFB, NCB, row_i_A, row_f_A, nelem

    REAL(8), ALLOCATABLE :: mat_A(:,:), sub_mat_B(:,:), sub_mat_C(:,:), mat_C(:,:)

    CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr);

    ! Division of Matrix A and B
    IF (me.EQ.0) THEN
        
        nelem = SIZE(Matrix_A,1) 
        

		ALLOCATE(scounts(nt), displs(nt))
		CALL Par_Div_32(nt, nelem, scounts, displs)
		
        ! ---------------------------------------------------------------------------
		sub_nt = nt
        ! ---------------------------------------------------------------------------
        
    END IF

    !---------------------------------------------------------------------------------
    ! New group according to the size of Matrix_A

    CALL MPI_Bcast(sub_nt,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
    
    CALL MPI_COMM_GROUP (MPI_COMM_WORLD, World_group, mpierr)
     
    ALLOCATE(ranks(sub_nt))
    ranks = (/(i,i=0,sub_nt-1)/)
      
    CALL MPI_GROUP_INCL (World_group, sub_nt, ranks, New_group, mpierr)
    CALL MPI_COMM_CREATE (MPI_COMM_WORLD, New_group, NEWCOMM, mpierr)
        
    sub_me = -1; sub_nt = -1
    
    IF (MPI_COMM_NULL .NE. NEWCOMM) THEN
        CALL MPI_Comm_rank(NEWCOMM, sub_me, mpierr);
        CALL MPI_Comm_size(NEWCOMM, sub_nt, mpierr);
        
        IF (sub_me.NE.0) THEN
            ALLOCATE(scounts(sub_nt), displs(sub_nt))
        END IF
        CALL MPI_Bcast(scounts,sub_nt,MPI_INTEGER,root,NEWCOMM,mpierr)
        CALL MPI_Bcast(displs,sub_nt,MPI_INTEGER,root,NEWCOMM,mpierr)

        ! Distribution of mat_A to all blocks
        IF (sub_me.EQ.0) THEN
            NFB = SIZE(Matrix_B,1); NCB = SIZE(Matrix_B,2)
        END IF

        CALL MPI_Bcast(NFB,1,MPI_INTEGER,root,NEWCOMM,mpierr)
        CALL MPI_Bcast(NCB,1,MPI_INTEGER,root,NEWCOMM,mpierr)
		
		
        CALL Distr_mat2(sub_nt,sub_me,NEWCOMM,scounts,displs,NFB,NCB, &
                                                              Matrix_B,sub_mat_B)

        n = 0
1       n = n + 1 ! Loop for divisions of Matrix_A
		
            ! Select mat_A (row matrix) from Matrix_A
            row_i_A = displs(n) + 1
            row_f_A = displs(n) + scounts(n) 
            
            ALLOCATE(mat_A(scounts(n),displs(sub_nt)+scounts(sub_nt)))
            IF (sub_me.EQ.0) THEN
                mat_A = Matrix_A(row_i_A:row_f_A,:)
            END IF

            CALL MPI_Bcast(mat_A,SIZE(mat_A),MPI_REAL8,root,NEWCOMM,mpierr)

            ! Result = mat_A * sub_mat_B
            ALLOCATE(sub_mat_C(scounts(n),scounts(sub_me+1)))
            sub_mat_C = MATMUL(mat_A,sub_mat_B)

            ! Join results in the host: mat_C
            IF (sub_me.EQ.0) THEN
                ALLOCATE(mat_C(scounts(n),displs(sub_nt)+scounts(sub_nt)))
            END IF
			
            CALL Join_mat2(sub_nt,sub_me,NEWCOMM,scounts,displs,sub_mat_C,mat_C)
            
			

            IF (sub_me.EQ.0) THEN

                Matrix_C(row_i_A:row_f_A,:) = mat_C
        
                DEALLOCATE(mat_C)

            END IF
            DEALLOCATE(mat_A,sub_mat_C)

        IF (n.LE.sub_nt-1) THEN
            GOTO 1 ! Go to the next element
        ELSE
            DEALLOCATE(scounts,displs)
            GOTO 2 ! All processors finish the task
        END IF

2       DEALLOCATE(sub_mat_B)
        
    END IF
    !---------------------------------------------------------------------------------

    IF (MPI_COMM_NULL .NE. NEWCOMM) THEN
        CALL MPI_COMM_FREE(NEWCOMM, mpierr)
    END IF
        
    CALL MPI_Group_free(World_group, mpierr);
    CALL MPI_Group_free(New_group, mpierr);

    DEALLOCATE(ranks)

END SUBROUTINE Matmul_parallel2
!===============================================================================!
SUBROUTINE Distr_mat2(sub_nt,sub_me,NEWCOMM,scounts,displs,NFA,NCA, &
                                                                    mat_A,sub_mat_A)
    include 'mpif.h'
    INTEGER :: sub_nt, root = 0, mpierr, sub_me, i, NEWCOMM, status
    INTEGER :: NCA, NFA
    INTEGER :: sub_nt2, World_group2, New_group2, NEWCOMM2, sub_me2
    INTEGER :: ndims, oldsize(2), newsize(2), starts(2), newtype
    REAL(8), ALLOCATABLE :: mat_A(:,:), sub_mat_A(:,:), sub_mat_A_aux(:,:)
    INTEGER, ALLOCATABLE :: scounts(:), displs(:)
    INTEGER, ALLOCATABLE :: ranks(:), scounts2(:), displs2(:)
    !----------------------------------------------------------------------------------------------

    CALL MPI_BARRIER(NEWCOMM,mpierr);

    ! --------------------------------------------------------------------------------

    ! Submatrices A
    
    ALLOCATE(sub_mat_A(NFA,scounts(sub_me+1)))

    ! Send blocks from 1 to sub_nt-1 to all processors -------------------------------
    sub_nt2 = sub_nt - 1
    CALL MPI_COMM_GROUP (NEWCOMM, World_group2, mpierr)
        
    ALLOCATE(ranks(sub_nt2))
    ranks = (/(i,i=0,sub_nt2-1)/)
      
    CALL MPI_GROUP_INCL (World_group2, sub_nt2, ranks, New_group2, mpierr)
    CALL MPI_COMM_CREATE (NEWCOMM, New_group2, NEWCOMM2, mpierr)
        
    sub_me2 = -1; sub_nt2 = -1

    IF (MPI_COMM_NULL .NE. NEWCOMM2) THEN
        CALL MPI_Comm_rank(NEWCOMM2, sub_me2, mpierr);
        CALL MPI_Comm_size(NEWCOMM2, sub_nt2, mpierr);

        ALLOCATE(scounts2(sub_nt2), displs2(sub_nt2))
        scounts2 = scounts(1:sub_nt2)
        displs2 = displs(1:sub_nt2)
        !---------------------------------------------------------------------------------------------
        ndims = 2
        starts = (/0,0/)

        ! Make a newtype array from the global array of matriz A =============
        oldsize = (/NFA,displs(sub_nt)/)
        newsize = (/NFA,scounts2(sub_me2+1)/)
        
        CALL MPI_TYPE_CREATE_SUBARRAY(ndims, oldsize, newsize, &
                                           starts, MPI_ORDER_FORTRAN, MPI_REAL8, newtype, mpierr)
        CALL MPI_TYPE_COMMIT(newtype,mpierr)
        
        ! Sending all blocks of same size to all processors inside the group
        CALL MPI_SCATTERV(mat_A,scounts2*NFA,displs2*NFA,MPI_REAL8,sub_mat_A, & 
                                         scounts2(sub_me2+1)*NFA,newtype,root,NEWCOMM2,mpierr) 
        CALL MPI_TYPE_FREE(newtype,mpierr)
        
        DEALLOCATE(ranks,scounts2, displs2)

    END IF
    
    IF (MPI_COMM_NULL .NE. NEWCOMM2) THEN
        CALL MPI_COMM_FREE(NEWCOMM2, mpierr)
    END IF
        
    CALL MPI_Group_free(World_group2, mpierr);
    CALL MPI_Group_free(New_group2, mpierr);

    ! Send the last block to processor sub_nt ------------------------------------------
    IF (sub_me.EQ.0) THEN

        ! Fot sub_mat_A
        ALLOCATE(sub_mat_A_aux(NFA,scounts(sub_nt)))   
        sub_mat_A_aux = mat_A(:,displs(sub_nt)+1:displs(sub_nt)+scounts(sub_nt))
        CALL MPI_SEND(sub_mat_A_aux, scounts(sub_nt)*NFA, MPI_REAL8, sub_nt-1, 100, NEWCOMM, mpierr ) 
        DEALLOCATE(sub_mat_A_aux)

    ELSEIF (sub_me .EQ. sub_nt-1) THEN

        ! Fot sub_mat_A
        CALL MPI_RECV(sub_mat_A, scounts(sub_nt)*NFA, MPI_REAL8, 0, 100, NEWCOMM, status, mpierr ) 
       
    END IF


END SUBROUTINE Distr_mat2
!===============================================================================!
SUBROUTINE Join_mat2(sub_nt,sub_me,NEWCOMM,scounts,displs, sub_mat_C,mat_C)
    include 'mpif.h'
    INTEGER :: sub_nt, root = 0, mpierr, sub_me, i, NEWCOMM
    
    INTEGER :: NFC
    
    REAL(8), ALLOCATABLE :: sub_mat_C(:,:), mat_C(:,:), sub_mat_C_aux(:,:)

    INTEGER, ALLOCATABLE :: scounts(:), displs(:)
    INTEGER :: sub_nt2, World_group2, New_group2, NEWCOMM2, sub_me2
    INTEGER, ALLOCATABLE :: ranks(:), scounts2(:), displs2(:)
    INTEGER :: status(MPI_STATUS_SIZE)
    !----------------------------------------------------------------------------------------------

    CALL MPI_BARRIER(NEWCOMM,mpierr);

    ! --------------------------------------------------------------------------------
    
    NFC = SIZE(sub_mat_C,1)      

    ! --------------------------------------------------------------------------------

    ! Send blocks from 1 to sub_nt-1 to all processors -------------------------------
    sub_nt2 = sub_nt - 1
    CALL MPI_COMM_GROUP (NEWCOMM, World_group2, mpierr)
        
    ALLOCATE(ranks(sub_nt2))
    ranks = (/(i,i=0,sub_nt2-1)/)
      
    CALL MPI_GROUP_INCL (World_group2, sub_nt2, ranks, New_group2, mpierr)
    CALL MPI_COMM_CREATE (NEWCOMM, New_group2, NEWCOMM2, mpierr)
        
    sub_me2 = -1; sub_nt2 = -1
    
    IF (MPI_COMM_NULL .NE. NEWCOMM2) THEN
        CALL MPI_Comm_rank(NEWCOMM2, sub_me2, mpierr);
        CALL MPI_Comm_size(NEWCOMM2, sub_nt2, mpierr);
        
        ALLOCATE(scounts2(sub_nt2), displs2(sub_nt2))
        scounts2 = scounts(1:sub_nt2)
        displs2 = displs(1:sub_nt2)
        !---------------------------------------------------------------------------------------------
                		
		CALL MPI_GATHERV(sub_mat_C,scounts2(sub_me+1)*NFC,MPI_REAL8,mat_C, &
                                            scounts2*NFC, displs2*NFC,MPI_REAL8,root,NEWCOMM2,mpierr)
		  
        !---------------------------------------------------------------------------------------------
        DEALLOCATE(ranks,scounts2, displs2)

    END IF
    
		
    IF (MPI_COMM_NULL .NE. NEWCOMM2) THEN
        CALL MPI_COMM_FREE(NEWCOMM2, mpierr)
    END IF
  
    CALL MPI_Group_free(World_group2, mpierr);
    CALL MPI_Group_free(New_group2, mpierr);
    	

    ! Send the last block to processor sub_nt ------------------------------------------
    IF (sub_me .EQ. sub_nt-1) THEN
        
        CALL MPI_SEND(sub_mat_C, scounts(sub_nt)*NFC, MPI_REAL8, 0, 101, NEWCOMM, mpierr ) 
                
    ELSEIF (sub_me .EQ. 0) THEN

        ALLOCATE(sub_mat_C_aux(NFC,scounts(sub_nt)))
        CALL MPI_RECV(sub_mat_C_aux, scounts(sub_nt)*NFC, MPI_REAL8, sub_nt-1, 101, NEWCOMM, status, mpierr ) 
        mat_C(:,displs(sub_nt)+1:displs(sub_nt)+scounts(sub_nt)) = sub_mat_C_aux
        DEALLOCATE(sub_mat_C_aux)

    END IF    



    !-----------------------------------------------------------------------------------------------
    
END SUBROUTINE Join_mat2
!===============================================================================!
SUBROUTINE Inv_matrix5(Matrix,Matrix2,sub_me,sub_nt,NEWCOMM)

	include 'mpif.h'

    INTEGER :: n, m1, m2,sub_me,sub_nt, mpierr,NEWCOMM
    REAL(8), ALLOCATABLE :: Matrix(:,:), Matrix2(:,:)
    REAL(8), ALLOCATABLE :: A11(:,:), A12(:,:), A21(:,:), A22(:,:)
    REAL(8), ALLOCATABLE :: P1(:,:), P2(:,:), P3(:,:), P4(:,:), P5(:,:), P6(:,:)
    REAL(8), ALLOCATABLE :: C11(:,:), C12(:,:), C21(:,:), C22(:,:)

	CALL MPI_BARRIER(NEWCOMM,mpierr);
	
	IF (sub_me.EQ.0) THEN

		! Division of Matrix into quarters
		n = SIZE(Matrix,1)
	  
		IF (MOD(n,2).EQ.0) THEN
		    m1 = n/2;
		    m2 = m1 + 1;
		ELSE
		    m1 = (n+1)/2;
		    m2 = m1 + 1;
		END IF
		
		ALLOCATE(A11(m1,m1),A12(m1,n-m1),A21(n-m1,m1),A22(n-m1,n-m1))
		ALLOCATE(P1(m1,m1),P2(n-m1,m1),P3(m1,n-m1))
		ALLOCATE(P4(n-m1,n-m1),P5(n-m1,n-m2),P6(n-m1,n-m1))

		A11 = Matrix(1:m1,1:m1); A12 = Matrix(1:m1,m2:n);
		A21 = Matrix(m2:n,1:m1); A22 = Matrix(m2:n,m2:n);
		
		!CALL Inv_matrix(A11,P1)

	END IF
		
		CALL Inv_matrix6(A11,P1,sub_me,sub_nt,NEWCOMM)
		
!!$OMP PARALLEL SHARED(A12,A21,A22,P1,P2,P3,P4,P5)
!!$OMP SECTIONS 
!!$OMP SECTION
!    P2 = MATMUL(A21,P1);

	CALL Matmul_parallel(sub_nt,sub_me,A21,P1,P2,NEWCOMM)
	
!!$OMP SECTION
!    P3 = MATMUL(P1,A12);

	CALL Matmul_parallel(sub_nt,sub_me,P1,A12,P3,NEWCOMM) 
	
!!$OMP SECTION
!    P4 = MATMUL(A21,MATMUL(P1,A12)); 

	CALL Matmul_parallel(sub_nt,sub_me,A21,P3,P4,NEWCOMM)

!!$OMP SECTION
!    P5 = MATMUL(A21,MATMUL(P1,A12))-A22; 
!!$OMP END SECTIONS NOWAIT
!!$OMP END PARALLEL
    
	IF (sub_me.EQ.0) THEN

		P5 = P4 - A22

   		!CALL Inv_matrix(P5,P6)
		
    
    	DEALLOCATE(A11,A12,A21,A22)
    	ALLOCATE(C11(m1,m1),C12(m1,n-m1),C21(n-m1,m1),C22(n-m1,n-m1))
		
	END IF
		
		CALL Inv_matrix6(P5,P6,sub_me,sub_nt,NEWCOMM)

!!$OMP PARALLEL SHARED(P1,P2,P3,P6,C11,C12,C21,C22)
!!$OMP SECTIONS 
!!$OMP SECTION
    !C12 = MATMUL(P3,P6)
	CALL Matmul_parallel(sub_nt,sub_me,P3,P6,C12,NEWCOMM)

!!$OMP SECTION 
!    C21 = MATMUL(P6,P2);
	CALL Matmul_parallel(sub_nt,sub_me,P6,P2,C21,NEWCOMM)
!!$OMP SECTION
!    C11 = P1-MATMUL(P3,MATMUL(P6,P2))
	CALL Matmul_parallel(sub_nt,sub_me,P3,C21,C11,NEWCOMM)
!!$OMP SECTION
!    C22 = -P6; 
!!$OMP END SECTIONS NOWAIT
!!$OMP END PARALLEL
    
	IF (sub_me.EQ.0) THEN

		C11 = P1 - C11
		C22 = -P6;

		DEALLOCATE(P1,P2,P3,P4,P5,P6)
		   
		Matrix2(1:m1,1:m1) = C11; Matrix2(1:m1,m2:n) = C12
		Matrix2(m2:n,1:m1) = C21; Matrix2(m2:n,m2:n) = C22
		    
		DEALLOCATE(C11,C12,C21,C22)

	END IF
	
END SUBROUTINE Inv_matrix5
!===============================================================================!
SUBROUTINE Inv_matrix6(Matrix,Matrix2,sub_me,sub_nt,NEWCOMM)

	include 'mpif.h'

    INTEGER :: n, m1, m2,sub_me,sub_nt, mpierr,NEWCOMM
    REAL(8), ALLOCATABLE :: Matrix(:,:), Matrix2(:,:)
    REAL(8), ALLOCATABLE :: A11(:,:), A12(:,:), A21(:,:), A22(:,:)
    REAL(8), ALLOCATABLE :: P1(:,:), P2(:,:), P3(:,:), P4(:,:), P5(:,:), P6(:,:)
    REAL(8), ALLOCATABLE :: C11(:,:), C12(:,:), C21(:,:), C22(:,:)

	CALL MPI_BARRIER(NEWCOMM,mpierr);
	
	IF (sub_me.EQ.0) THEN

		! Division of Matrix into quarters
		n = SIZE(Matrix,1)
	  
		IF (MOD(n,2).EQ.0) THEN
		    m1 = n/2;
		    m2 = m1 + 1;
		ELSE
		    m1 = (n+1)/2;
		    m2 = m1 + 1;
		END IF
		
		ALLOCATE(A11(m1,m1),A12(m1,n-m1),A21(n-m1,m1),A22(n-m1,n-m1))
		ALLOCATE(P1(m1,m1),P2(n-m1,m1),P3(m1,n-m1))
		ALLOCATE(P4(n-m1,n-m1),P5(n-m1,n-m2),P6(n-m1,n-m1))

		A11 = Matrix(1:m1,1:m1); A12 = Matrix(1:m1,m2:n);
		A21 = Matrix(m2:n,1:m1); A22 = Matrix(m2:n,m2:n);
		
		!CALL Inv_matrix(A11,P1)

	END IF
		
		CALL Inv_matrix7(A11,P1,sub_me,sub_nt,NEWCOMM)
		
!!$OMP PARALLEL SHARED(A12,A21,A22,P1,P2,P3,P4,P5)
!!$OMP SECTIONS 
!!$OMP SECTION
!    P2 = MATMUL(A21,P1);

	CALL Matmul_parallel(sub_nt,sub_me,A21,P1,P2,NEWCOMM)
	
!!$OMP SECTION
!    P3 = MATMUL(P1,A12);

	CALL Matmul_parallel(sub_nt,sub_me,P1,A12,P3,NEWCOMM) 
	
!!$OMP SECTION
!    P4 = MATMUL(A21,MATMUL(P1,A12)); 

	CALL Matmul_parallel(sub_nt,sub_me,A21,P3,P4,NEWCOMM)

!!$OMP SECTION
!    P5 = MATMUL(A21,MATMUL(P1,A12))-A22; 
!!$OMP END SECTIONS NOWAIT
!!$OMP END PARALLEL
    
	IF (sub_me.EQ.0) THEN

		P5 = P4 - A22

   		!CALL Inv_matrix(P5,P6)
		
    
    	DEALLOCATE(A11,A12,A21,A22)
    	ALLOCATE(C11(m1,m1),C12(m1,n-m1),C21(n-m1,m1),C22(n-m1,n-m1))
		
	END IF
		
		CALL Inv_matrix7(P5,P6,sub_me,sub_nt,NEWCOMM)

!!$OMP PARALLEL SHARED(P1,P2,P3,P6,C11,C12,C21,C22)
!!$OMP SECTIONS 
!!$OMP SECTION
    !C12 = MATMUL(P3,P6)
	CALL Matmul_parallel(sub_nt,sub_me,P3,P6,C12,NEWCOMM)

!!$OMP SECTION 
!    C21 = MATMUL(P6,P2);
	CALL Matmul_parallel(sub_nt,sub_me,P6,P2,C21,NEWCOMM)
!!$OMP SECTION
!    C11 = P1-MATMUL(P3,MATMUL(P6,P2))
	CALL Matmul_parallel(sub_nt,sub_me,P3,C21,C11,NEWCOMM)
!!$OMP SECTION
!    C22 = -P6; 
!!$OMP END SECTIONS NOWAIT
!!$OMP END PARALLEL
    
	IF (sub_me.EQ.0) THEN

		C11 = P1 - C11
		C22 = -P6;

		DEALLOCATE(P1,P2,P3,P4,P5,P6)
		   
		Matrix2(1:m1,1:m1) = C11; Matrix2(1:m1,m2:n) = C12
		Matrix2(m2:n,1:m1) = C21; Matrix2(m2:n,m2:n) = C22
		    
		DEALLOCATE(C11,C12,C21,C22)

	END IF
	
END SUBROUTINE Inv_matrix6
!===============================================================================!
SUBROUTINE Inv_matrix7(Matrix,Matrix2,sub_me,sub_nt,NEWCOMM)

	include 'mpif.h'

    INTEGER :: n, m1, m2,sub_me,sub_nt, mpierr,NEWCOMM
    REAL(8), ALLOCATABLE :: Matrix(:,:), Matrix2(:,:)
    REAL(8), ALLOCATABLE :: A11(:,:), A12(:,:), A21(:,:), A22(:,:)
    REAL(8), ALLOCATABLE :: P1(:,:), P2(:,:), P3(:,:), P4(:,:), P5(:,:), P6(:,:)
    REAL(8), ALLOCATABLE :: C11(:,:), C12(:,:), C21(:,:), C22(:,:)

	CALL MPI_BARRIER(NEWCOMM,mpierr);
	
	IF (sub_me.EQ.0) THEN

		! Division of Matrix into quarters
		n = SIZE(Matrix,1)
	  
		IF (MOD(n,2).EQ.0) THEN
		    m1 = n/2;
		    m2 = m1 + 1;
		ELSE
		    m1 = (n+1)/2;
		    m2 = m1 + 1;
		END IF
		
		ALLOCATE(A11(m1,m1),A12(m1,n-m1),A21(n-m1,m1),A22(n-m1,n-m1))
		ALLOCATE(P1(m1,m1),P2(n-m1,m1),P3(m1,n-m1))
		ALLOCATE(P4(n-m1,n-m1),P5(n-m1,n-m2),P6(n-m1,n-m1))

		A11 = Matrix(1:m1,1:m1); A12 = Matrix(1:m1,m2:n);
		A21 = Matrix(m2:n,1:m1); A22 = Matrix(m2:n,m2:n);
		
		!CALL Inv_matrix(A11,P1)

	END IF
		
	CALL Inv_matrix8(A11,P1,sub_me,sub_nt,NEWCOMM)
		
!!$OMP PARALLEL SHARED(A12,A21,A22,P1,P2,P3,P4,P5)
!!$OMP SECTIONS 
!!$OMP SECTION
!    P2 = MATMUL(A21,P1);

	CALL Matmul_parallel(sub_nt,sub_me,A21,P1,P2,NEWCOMM)
	
!!$OMP SECTION
!    P3 = MATMUL(P1,A12);

	CALL Matmul_parallel(sub_nt,sub_me,P1,A12,P3,NEWCOMM) 
	
!!$OMP SECTION
!    P4 = MATMUL(A21,MATMUL(P1,A12)); 

	CALL Matmul_parallel(sub_nt,sub_me,A21,P3,P4,NEWCOMM)

!!$OMP SECTION
!    P5 = MATMUL(A21,MATMUL(P1,A12))-A22; 
!!$OMP END SECTIONS NOWAIT
!!$OMP END PARALLEL
    
	IF (sub_me.EQ.0) THEN

		P5 = P4 - A22

   		!CALL Inv_matrix(P5,P6)
		
    
    	DEALLOCATE(A11,A12,A21,A22)
    	ALLOCATE(C11(m1,m1),C12(m1,n-m1),C21(n-m1,m1),C22(n-m1,n-m1))
		
	END IF
		
		CALL Inv_matrix8(P5,P6,sub_me,sub_nt,NEWCOMM)

!!$OMP PARALLEL SHARED(P1,P2,P3,P6,C11,C12,C21,C22)
!!$OMP SECTIONS 
!!$OMP SECTION
    !C12 = MATMUL(P3,P6)
	CALL Matmul_parallel(sub_nt,sub_me,P3,P6,C12,NEWCOMM)

!!$OMP SECTION 
!    C21 = MATMUL(P6,P2);
	CALL Matmul_parallel(sub_nt,sub_me,P6,P2,C21,NEWCOMM)
!!$OMP SECTION
!    C11 = P1-MATMUL(P3,MATMUL(P6,P2))
	CALL Matmul_parallel(sub_nt,sub_me,P3,C21,C11,NEWCOMM)
!!$OMP SECTION
!    C22 = -P6; 
!!$OMP END SECTIONS NOWAIT
!!$OMP END PARALLEL
    
	IF (sub_me.EQ.0) THEN

		C11 = P1 - C11
		C22 = -P6;

		DEALLOCATE(P1,P2,P3,P4,P5,P6)
		   
		Matrix2(1:m1,1:m1) = C11; Matrix2(1:m1,m2:n) = C12
		Matrix2(m2:n,1:m1) = C21; Matrix2(m2:n,m2:n) = C22
		    
		DEALLOCATE(C11,C12,C21,C22)

	END IF
	
END SUBROUTINE Inv_matrix7
!===============================================================================!
SUBROUTINE Inv_matrix8(Matrix,Matrix2,sub_me,sub_nt,NEWCOMM)

	include 'mpif.h'

    INTEGER :: n, m1, m2,sub_me,sub_nt, mpierr,NEWCOMM
    REAL(8), ALLOCATABLE :: Matrix(:,:), Matrix2(:,:)
    REAL(8), ALLOCATABLE :: A11(:,:), A12(:,:), A21(:,:), A22(:,:)
    REAL(8), ALLOCATABLE :: P1(:,:), P2(:,:), P3(:,:), P4(:,:), P5(:,:), P6(:,:)
    REAL(8), ALLOCATABLE :: C11(:,:), C12(:,:), C21(:,:), C22(:,:)

	CALL MPI_BARRIER(NEWCOMM,mpierr);
	
	IF (sub_me.EQ.0) THEN

		! Division of Matrix into quarters
		n = SIZE(Matrix,1)
	  
		IF (MOD(n,2).EQ.0) THEN
		    m1 = n/2;
		    m2 = m1 + 1;
		ELSE
		    m1 = (n+1)/2;
		    m2 = m1 + 1;
		END IF
		
		ALLOCATE(A11(m1,m1),A12(m1,n-m1),A21(n-m1,m1),A22(n-m1,n-m1))
		ALLOCATE(P1(m1,m1),P2(n-m1,m1),P3(m1,n-m1))
		ALLOCATE(P4(n-m1,n-m1),P5(n-m1,n-m2),P6(n-m1,n-m1))

		A11 = Matrix(1:m1,1:m1); A12 = Matrix(1:m1,m2:n);
		A21 = Matrix(m2:n,1:m1); A22 = Matrix(m2:n,m2:n);
		
		CALL Inv_matrix(A11,P1)

	END IF
		
	!CALL Inv_matrix7(A11,P1,sub_me,sub_nt,NEWCOMM)
		
!!$OMP PARALLEL SHARED(A12,A21,A22,P1,P2,P3,P4,P5)
!!$OMP SECTIONS 
!!$OMP SECTION
!    P2 = MATMUL(A21,P1);

	CALL Matmul_parallel(sub_nt,sub_me,A21,P1,P2,NEWCOMM)
	
!!$OMP SECTION
!    P3 = MATMUL(P1,A12);

	CALL Matmul_parallel(sub_nt,sub_me,P1,A12,P3,NEWCOMM) 
	
!!$OMP SECTION
!    P4 = MATMUL(A21,MATMUL(P1,A12)); 

	CALL Matmul_parallel(sub_nt,sub_me,A21,P3,P4,NEWCOMM)

!!$OMP SECTION
!    P5 = MATMUL(A21,MATMUL(P1,A12))-A22; 
!!$OMP END SECTIONS NOWAIT
!!$OMP END PARALLEL
    
	IF (sub_me.EQ.0) THEN

		P5 = P4 - A22

   		CALL Inv_matrix(P5,P6)
		
    
    	DEALLOCATE(A11,A12,A21,A22)
    	ALLOCATE(C11(m1,m1),C12(m1,n-m1),C21(n-m1,m1),C22(n-m1,n-m1))
		
	END IF
		
		!CALL Inv_matrix7(P5,P6,sub_me,sub_nt,NEWCOMM)

!!$OMP PARALLEL SHARED(P1,P2,P3,P6,C11,C12,C21,C22)
!!$OMP SECTIONS 
!!$OMP SECTION
    !C12 = MATMUL(P3,P6)
	CALL Matmul_parallel(sub_nt,sub_me,P3,P6,C12,NEWCOMM)

!!$OMP SECTION 
!    C21 = MATMUL(P6,P2);
	CALL Matmul_parallel(sub_nt,sub_me,P6,P2,C21,NEWCOMM)
!!$OMP SECTION
!    C11 = P1-MATMUL(P3,MATMUL(P6,P2))
	CALL Matmul_parallel(sub_nt,sub_me,P3,C21,C11,NEWCOMM)
!!$OMP SECTION
!    C22 = -P6; 
!!$OMP END SECTIONS NOWAIT
!!$OMP END PARALLEL
    
	IF (sub_me.EQ.0) THEN

		C11 = P1 - C11
		C22 = -P6;

		DEALLOCATE(P1,P2,P3,P4,P5,P6)
		   
		Matrix2(1:m1,1:m1) = C11; Matrix2(1:m1,m2:n) = C12
		Matrix2(m2:n,1:m1) = C21; Matrix2(m2:n,m2:n) = C22
		    
		DEALLOCATE(C11,C12,C21,C22)

	END IF
	
END SUBROUTINE Inv_matrix8
!===============================================================================!
SUBROUTINE connectivity_plot(Field)

	INTEGER :: i, no, nob, cont, j
	REAL(8), ALLOCATABLE :: Field(:,:), field_no(:)

	ALLOCATE(field_no(SIZE(Field,2)))

	DO i=1,SIZE(nodes_conectivity,1)

		no = nodes_conectivity(i,1)

		field_no = Field(no,:)

		cont = nodes_conectivity(i,2)

		DO j = 1,cont

			nob = nodes_conectivity(i,j+3)		

			field_no = field_no + Field(nob,:)

		END DO

		field_no = field_no/(cont+1)

		Field(no,:) = field_no

		DO j = 1,cont

			nob = nodes_conectivity(i,j+3)

			Field(nob,:) = field_no

		END DO

	END DO

	DEALLOCATE(field_no)


END SUBROUTINE connectivity_plot
!===============================================================================!
SUBROUTINE connectivity_plot_disp(Field)

	INTEGER :: i, no, nob, cont, j, el, dof, k
	REAL(8), ALLOCATABLE :: Field(:,:), field_no(:)
	REAL(8) :: Type_cond, val_cond, field_dof

	dof = 3;

	ALLOCATE(field_no(SIZE(Field,2)))

	DO i=1,SIZE(nodes_conectivity,1)

		el = nodes_conectivity(i,3) 

		no = nodes_conectivity(i,1)

		field_no = Field(no,:)

		cont = nodes_conectivity(i,2)
		
		IF (el .NE. 0) THEN ! displacement is imposed at any dof of "el"

			DO j=1,dof

				IF (nreg.EQ.1) THEN

					Type_cond = INT(BC_elem(el,2*j))
					val_cond = BC_elem(el,2*j+1)

				ELSE

					Type_cond = INT(BC_elem(el,2*j+2))
					val_cond = BC_elem(el,2*j+3)

				END IF


				IF (Type_cond.EQ.0) THEN ! average is not necessary

					Field(no,j) = val_cond;
					
					DO k=1,cont

						nob = nodes_conectivity(i,k+3)

						Field(nob,j) = val_cond; 
						
					END DO

				ELSE ! average is necessary

					field_dof = Field(no,j)
					DO k = 1,cont

						nob = nodes_conectivity(i,k+3)		

						field_dof = field_dof + Field(nob,j)

					END DO

					field_dof = field_dof/(cont+1)

					Field(no,j) = field_dof

					DO k = 1,cont

						nob = nodes_conectivity(i,k+3)

						Field(nob,j) = field_dof

					END DO

				END IF


			END DO


		ELSEIF(el .EQ. 0) THEN 

			DO j = 1,cont

				nob = nodes_conectivity(i,j+3)		

				field_no = field_no + Field(nob,:)

			END DO

			field_no = field_no/(cont+1)

			Field(no,:) = field_no

			DO j = 1,cont

				nob = nodes_conectivity(i,j+3)

				Field(nob,:) = field_no

			END DO

		END IF

	END DO

	DEALLOCATE(field_no)


END SUBROUTINE connectivity_plot_disp
!===============================================================================!
SUBROUTINE Check_Size(nt)

	INTEGER:: nt
    INTEGER ::i, nel_reg,dof=3
    INTEGER :: nint_reg, nelem

	INTEGER(8) :: ndofs_reg1, ndofs_reg2
	INTEGER(8) :: n_A, n_G, n_V, nt_sub, max_size, val, nzA
	
	INTEGER(8), ALLOCATABLE :: s_A_aux(:) 
	INTEGER(8), ALLOCATABLE :: displs_A(:) 

	IF (nreg .GT. 1) THEN

		nelem = SIZE(ELEM,1) ! Number of elements
		nt_sub = nt - 1 ! Threads without the master "0"

		n_V = 0
		n_G = 0
		DO i=1,nreg

		    nel_reg = El_reg(i,1)
		    ndofs_reg1 = nel_reg + El_reg(i,7) + El_reg(i,8) 
		    n_V = n_V + ndofs_reg1**2 

		    nint_reg = Int_reg(i,2)
		    ndofs_reg2 = nint_reg + El_reg(i,7) + El_reg(i,8)       
		    n_G = n_G + (ndofs_reg2*ndofs_reg1)

		END DO
		
		n_A = n_V + n_G 
		n_G = n_V - n_G

		n_A = n_A*(nnos_el*dof)**2 - Reduce_A ! Total size of matrix A
		n_G = n_G*(nnos_el*dof)**2 ! Total size of matrix G

		ALLOCATE(s_A_aux(nt-1))
		s_A_aux = 0; 

		ALLOCATE(displs_A(nt-1))
		displs_A = 0; 
		
		N_total = nelem*nnos_el*dof ! rows or columns  of general system

		! Function for parallel distribution
		! Par_Div(threads, vector_size, scounts)
		CALL Par_Div(nt_sub, n_A, s_A_aux, displs_A) ! A_values
		
		max_size = MAXVAL(s_A_aux)

		IF (max_size .GT. (2**31-1)*0.7) THEN

			WRITE(*,*) "   THE PROBLEM EXCEEDS THE MAXIMUM NUMBER OF POSSIBLE DEGREES OF FREEDOM"
			WRITE(*,*) ""
			
			DO i=1,nt-1
				val = s_A_aux(i)-INT((2**31-1)*0.7)
				WRITE(*,'(A15,I4,A20,I14)') "Rank: ",i-1, "Additional DOF: ", val


			END DO
			WRITE(*,*) ""

			call exit(123)

		END IF

		DEALLOCATE (s_A_aux,displs_A) 

	ELSE
	
		nt_sub = nt - 1 ! Threads without the master "0"
		nelem = SIZE(ELEM,1) ! Number of elements
		N_total = nelem*nnos_el*dof ! rows or columns  of general system
		nzA = N_total

		ALLOCATE(s_A_aux(nt-1))
		s_A_aux = 0; 

		ALLOCATE(displs_A(nt-1))
		displs_A = 0; 


		! Function for parallel distribution
		! Par_Div(threads, vector_size, scounts)
		CALL Par_Div(nt_sub, nzA, s_A_aux, displs_A) ! A_values	


		max_size = MAXVAL(s_A_aux)

		IF (max_size .GT. (2**31-1)*0.7) THEN

			WRITE(*,*) "   THE PROBLEM EXCEEDS THE MAXIMUM NUMBER OF POSIBLE DEGREES OF FREEDOM"
			WRITE(*,*) ""
			
			DO i=1,nt-1
				val = s_A_aux(i)-INT((2**31-1)*0.7)
				WRITE(*,'(A15,I4,A20,I14)') "Rank: ",i-1, "Additional DOF: ", val


			END DO
			WRITE(*,*) ""

			call exit(123)

		END IF

		DEALLOCATE (s_A_aux,displs_A)


	END IF

END SUBROUTINE Check_Size
!===============================================================================!
SUBROUTINE Check_Size_nt(nt)

	INTEGER:: nt 

	
		IF (nt .EQ. 1) THEN
			WRITE(*,*) ""
			WRITE(*,*) "                ***** MORE THAN ONE THREAD IS REQUIRED ***** "
			WRITE(*,*) ""
			
			call exit(123)

		END IF


END SUBROUTINE Check_Size_nt
!===============================================================================!
END MODULE Global_functions
