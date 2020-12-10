!===============================================================================!
!-------------------------------------------------------------------------------!
!          MODULE TO COMPUTE THE INTEGRATION POINTS FOR A                       !
!       DISCONTINUOUS QUADRATIC SIX-NODE TRIANGLE BEM ELEMENT                   !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Int_points
!-------------------------------------------------------------------------------!
USE Global_variables
USE Global_functions
!USE mpi
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE Integration_Points(me,time)

    include 'mpif.h'
    INTEGER :: i, cont, j,ii,jj
    INTEGER, PARAMETER :: ngp=4
    INTEGER :: V_divisions(3)
    CHARACTER(4), ALLOCATABLE :: Tri_Quad(:) 
    REAL(8), ALLOCATABLE :: rp(:,:), rw(:), N_nsc(:,:), N_sc(:,:)
    REAL(8) :: qsi, eta, w, N(3), dNdqsi(3), dNdeta(3), x, y
    REAL(8) :: X1(2), X2(2), X3(2), N2(3), dN2dqsi(3), dN2deta(3), Je
    REAL(8) :: gp(ngp), gw(ngp), X4(2), w1, w2, N3(4), dN3dqsi(4), dN3deta(4)
    REAL(8) :: positions_Geo_Nodes(3,2), S_F(3)

    REAL(8) :: t1, t2, duration, time
    INTEGER :: me, mpierr, root=0
    
    t1 = MPI_Wtime();

    CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr);

    IF (me.EQ.0) THEN

        WRITE (*,*) '----------------------------------------------------------------------------'
        WRITE (*,*) ' '
        WRITE (*,'(A)',advance='no') ' 0. Matrices of integration points ...'

    END IF

    ! Matrix G to compute the shape functions of discontinuous elements
    CALL Matrix_G(N_Phy_nodes)

    ! Derivative of the shape functions
    CALL Dev_shape_functions(dNdqsi,dNdeta)
    dNda_i(1,:) = dNdqsi
    dNda_i(2,:) = dNdeta

    !====================================================================
    ! Coordinates of the geometrical nodes in the parametric element
    positions_Geo_Nodes(1,:) = (/1.d0  , 0.d0/)
    positions_Geo_Nodes(2,:) = (/0.d0  , 1.d0/)
    positions_Geo_Nodes(3,:) = (/0.d0  , 0.d0/)
    S_F = 0.d0
    DO i=1,nnos_el
        ! Parametric coordinates
        qsi = positions_Geo_Nodes(i,1); eta = positions_Geo_Nodes(i,2)
        ! Shape functions of three-node triangle discontinuous element
        CALL Shape_functions(qsi,eta,S_F) 
        N_Geo_nodes2(3*i-2:3*i,1) = S_F

    END DO
    !====================================================================

    ! NON-SINGULAR INTEGRATION - [weights_ns]
    
    ! Number of integration points in the entire triangle or in each 
    ! subtriangles depending on the case
    ! npoints_rule_NSI (7) or (13)  
    
    ALLOCATE(rp(npoints_rule_NSI,3), rw(npoints_rule_NSI))
    ! rp: rule points, rw: weights
    CALL Rule_points_weights(npoints_rule_NSI,rp,rw)   
    
    IF (Int_type.EQ.1) THEN ! Entire triangle
        
        ALLOCATE(weights_ns(npoints_rule_NSI,12)) ! Weights matrix

        DO i=1,npoints_rule_NSI
            qsi = rp(i,1); eta = rp(i,2); w = rw(i)
            ! Shape functions of three-node triangle discontinous element
            CALL Shape_functions(qsi,eta,N)
            ! Derivative of the shape functions
            CALL Dev_shape_functions(dNdqsi,dNdeta)
            ! Non-singular weights matrix
            weights_ns(i,1:2) = (/qsi,eta/)
            weights_ns(i,3:5) = N; weights_ns(i,6:8) = dNdqsi; 
            weights_ns(i,9:11) = dNdeta; weights_ns(i,12) = 0.5d0*w  
        END DO
    
    ! Triangle divided into 4 sub-triangles
    ELSEIF(Int_type.EQ.4) THEN
        
        ALLOCATE(weights_ns(4*npoints_rule_NSI,12)) ! Weights matrix
        ! Points of the new sub-triangles
		ALLOCATE(N_sc(4,6))

        N_sc(1,:) = (/0.5d0,0.5d0,   0.d0, 1.d0,   0.d0 ,0.5d0/)
        N_sc(2,:) = (/0.5d0, 0.d0,   0.d0,0.5d0,   0.5d0,0.5d0/)
        N_sc(3,:) = (/0.5d0, 0.d0,   0.d0,0.5d0,   0.d0 , 0.d0/)
        N_sc(4,:) = (/0.5d0,0.5d0,   1.d0, 0.d0,   0.5d0, 0.d0/)
        cont = 0;

        DO i=1,4
            ! Coordinates
            X1 = N_sc(i,1:2); X2 = N_sc(i,3:4); X3 = N_sc(i,5:6)
            DO j=1,npoints_rule_NSI
                qsi = rp(j,1); eta = rp(j,2); w = rw(j)
                ! shape functions of three-node triangle continous element
                CALL Shape_functions_cont(qsi,eta,1,N2)
                ! Derivative of the shape functions of three-node triangle continous element
                CALL Dev_shape_funct_cont(qsi,eta,1,dN2dqsi,dN2deta)
                ! Jacobian
                CALL Calc_jac_tri3(X1,X2,X3,dN2dqsi,dN2deta,Je)
                ! x and y coordinates
                x = N2(1)*X1(1)+N2(2)*X2(1)+N2(3)*X3(1)
                y = N2(1)*X1(2)+N2(2)*X2(2)+N2(3)*X3(2)  
                ! shape functions of three-node triangle discontinous element
                CALL Shape_functions(x,y,N)
                ! Derivative of the shape functions
                CALL Dev_shape_functions(dNdqsi,dNdeta)
                ! Non-singular weights matrix
                cont = cont + 1
                weights_ns(cont,1:2) = (/x,y/)
                weights_ns(cont,3:5) = N;        weights_ns(cont,6:8) = dNdqsi; 
                weights_ns(cont,9:11) = dNdeta; weights_ns(cont,12) = Je*w 
            END DO

        END DO

	! Triangle divided into 8 sub-triangles
    ELSEIF(Int_type.EQ.8) THEN

		ALLOCATE(weights_ns(8*npoints_rule_NSI,12)) ! Weights matrix
        ! Points of the new sub-triangles
		ALLOCATE(N_sc(8,6))

        N_sc(1,:) = (/  0.d0,  0.5d0,   0.25d0, 0.75d0,     0.d0,   1.d0/)
        N_sc(2,:) = (/  0.d0,  0.5d0,    0.5d0,  0.5d0,   0.25d0, 0.75d0/)
        N_sc(3,:) = (/0.25d0, 0.25d0,    0.5d0,  0.5d0,     0.d0,  0.5d0/)
        N_sc(4,:) = (/  0.d0,   0.d0,   0.25d0, 0.25d0,     0.d0,  0.5d0/)
		N_sc(5,:) = (/  0.d0,   0.d0,    0.5d0,   0.d0,   0.25d0, 0.25d0/)
        N_sc(6,:) = (/ 0.5d0,   0.d0,    0.5d0,  0.5d0,   0.25d0, 0.25d0/)
        N_sc(7,:) = (/ 0.5d0,   0.d0,   0.75d0, 0.25d0,    0.5d0,  0.5d0/)
        N_sc(8,:) = (/ 0.5d0,   0.d0,   0.75d0, 0.25d0,     1.d0,   0.d0/)

        cont = 0;

		DO i=1,8
            ! Coordinates
            X1 = N_sc(i,1:2); X2 = N_sc(i,3:4); X3 = N_sc(i,5:6)
            DO j=1,npoints_rule_NSI
                qsi = rp(j,1); eta = rp(j,2); w = rw(j)
                ! shape functions of three-node triangle continous element
                CALL Shape_functions_cont(qsi,eta,1,N2)
                ! Derivative of the shape functions of three-node triangle continous element
                CALL Dev_shape_funct_cont(qsi,eta,1,dN2dqsi,dN2deta)
                ! Jacobian
                CALL Calc_jac_tri3(X1,X2,X3,dN2dqsi,dN2deta,Je)
                ! x and y coordinates
                x = N2(1)*X1(1)+N2(2)*X2(1)+N2(3)*X3(1)
                y = N2(1)*X1(2)+N2(2)*X2(2)+N2(3)*X3(2)  
                ! shape functions of three-node triangle discontinous element
                CALL Shape_functions(x,y,N)
                ! Derivative of the shape functions
                CALL Dev_shape_functions(dNdqsi,dNdeta)
                ! Non-singular weights matrix
                cont = cont + 1
                weights_ns(cont,1:2) = (/x,y/)
                weights_ns(cont,3:5) = N;        weights_ns(cont,6:8) = dNdqsi; 
                weights_ns(cont,9:11) = dNdeta; weights_ns(cont,12) = Je*w 
            END DO

        END DO

	! Triangle divided into 16 sub-triangles
    ELSEIF(Int_type.EQ.16) THEN

		ALLOCATE(weights_ns(16*npoints_rule_NSI,12)) ! Weights matrix
        ! Points of the new sub-triangles
		ALLOCATE(N_sc(16,6))

        N_sc(1,:) = (/  0.d0, 0.75d0,   0.25d0, 0.75d0,     0.d0,   1.d0/)
        N_sc(2,:) = (/  0.d0,  0.5d0,   0.25d0, 0.75d0,     0.d0, 0.75d0/)
        N_sc(3,:) = (/  0.d0,  0.5d0,   0.25d0,  0.5d0,   0.25d0, 0.75d0/)
        N_sc(4,:) = (/0.25d0,  0.5d0,    0.5d0,  0.5d0,   0.25d0, 0.75d0/)
		N_sc(5,:) = (/  0.d0, 0.25d0,   0.25d0, 0.25d0,     0.d0,  0.5d0/)
        N_sc(6,:) = (/0.25d0, 0.25d0,   0.25d0,  0.5d0,     0.d0,  0.5d0/)
        N_sc(7,:) = (/0.25d0, 0.25d0,    0.5d0,  0.5d0,   0.25d0,  0.5d0/)
        N_sc(8,:) = (/0.25d0, 0.25d0,    0.5d0, 0.25d0,    0.5d0,  0.5d0/)
        N_sc(9,:) = (/ 0.5d0, 0.25d0,   0.75d0, 0.25d0,    0.5d0,  0.5d0/)
		N_sc(10,:) = (/  0.d0,   0.d0,   0.25d0, 0.25d0,     0.d0, 0.25d0/)
		N_sc(11,:) = (/  0.d0,   0.d0,   0.25d0, 0.25d0,   0.25d0,   0.d0/)
		N_sc(12,:) = (/0.25d0,   0.d0,    0.5d0,   0.d0,   0.25d0, 0.25d0/)
		N_sc(13,:) = (/0.25d0, 0.25d0,    0.5d0,   0.d0,    0.5d0, 0.25d0/)
		N_sc(14,:) = (/ 0.5d0,   0.d0,   0.75d0, 0.25d0,    0.5d0, 0.25d0/)
		N_sc(15,:) = (/ 0.5d0,   0.d0,   0.75d0,   0.d0,   0.75d0, 0.25d0/)
        N_sc(16,:) = (/0.75d0,   0.d0,     1.d0,   0.d0,   0.75d0, 0.25d0/)

        cont = 0;

		DO i=1,16
            ! Coordinates
            X1 = N_sc(i,1:2); X2 = N_sc(i,3:4); X3 = N_sc(i,5:6)
            DO j=1,npoints_rule_NSI
                qsi = rp(j,1); eta = rp(j,2); w = rw(j)
                ! shape functions of three-node triangle continous element
                CALL Shape_functions_cont(qsi,eta,1,N2)
                ! Derivative of the shape functions of three-node triangle continous element
                CALL Dev_shape_funct_cont(qsi,eta,1,dN2dqsi,dN2deta)
                ! Jacobian
                CALL Calc_jac_tri3(X1,X2,X3,dN2dqsi,dN2deta,Je)
                ! x and y coordinates
                x = N2(1)*X1(1)+N2(2)*X2(1)+N2(3)*X3(1)
                y = N2(1)*X1(2)+N2(2)*X2(2)+N2(3)*X3(2)  
                ! shape functions of three-node triangle discontinous element
                CALL Shape_functions(x,y,N)
                ! Derivative of the shape functions
                CALL Dev_shape_functions(dNdqsi,dNdeta)
                ! Non-singular weights matrix
                cont = cont + 1
                weights_ns(cont,1:2) = (/x,y/)
                weights_ns(cont,3:5) = N;        weights_ns(cont,6:8) = dNdqsi; 
                weights_ns(cont,9:11) = dNdeta; weights_ns(cont,12) = Je*w 
            END DO

        END DO

    END IF
    
    ! Number of points for non-singular integration
    npoints_Nsi = SIZE(weights_ns,1)

    DEALLOCATE(rp,rw)
    
    ! SINGULAR INTEGRATION - int_data%

    ! Number of integration points for each sub-triangles
    ! npoints_rule_SI (7) or (13) 

    ALLOCATE(rp(npoints_rule_SI,3), rw(npoints_rule_SI))
    ! rp: rule points, rw: weights
    CALL Rule_points_weights(npoints_rule_SI,rp,rw)
    ! gp: Gauss points, gw: weights
    CALL Gauss_Legendre(-1.d0,1.d0,ngp,gp,gw)
    
    ! Number of divisions in the element
    V_divisions = (/7,7,8/)
    ! Total number of integration points per element
     N_int_points(1) = npoints_rule_SI + 6*(ngp**2)
     N_int_points(2) = npoints_rule_SI + 6*(ngp**2)
     N_int_points(3) = 8*(ngp**2)
    
    DO i=1,nnos_el
        
        ! Configuration points when the singularity is at node (i)
        ALLOCATE(N_nsc(V_divisions(i),8))
        ALLOCATE(Tri_Quad(V_divisions(i)))
        CALL calc_pts_node(i,N_nsc,Tri_Quad)

        ALLOCATE(Int_data(i)%weights_s(N_int_points(i),12))
        cont=0
            
        DO j=1,V_divisions(i) ! loop for all triangles and quads
            ! Coordinates
            X1 = N_nsc(j,1:2); X2 = N_nsc(j,3:4); 
            X3 = N_nsc(j,5:6); X4 = N_nsc(j,7:8);    
            
            IF (Tri_Quad(j).EQ.'Tria') THEN ! Triangle shape function
                DO ii=1, npoints_rule_SI
                    qsi = rp(ii,1); eta = rp(ii,2); w = rw(ii);
                    ! Shape functions of three-node triangle continuous element
                    CALL  Shape_functions_cont(qsi,eta,1,N2)
                    ! Derivative of the shape functions of three-node triangle continous element
                    CALL Dev_shape_funct_cont(qsi,eta,1,dN2dqsi,dN2deta)
                    ! Jacobian
                    CALL Calc_jac_tri3(X1,X2,X3,dN2dqsi,dN2deta,Je)
                    ! x and y coordinates
                    x = N2(1)*X1(1)+N2(2)*X2(1)+N2(3)*X3(1)
                    y = N2(1)*X1(2)+N2(2)*X2(2)+N2(3)*X3(2)
                    ! shape functions of three-node triangle discontinous element
                    CALL Shape_functions(x,y,N)      
                    ! Derivative of the shape functions
                    CALL Dev_shape_functions(dNdqsi,dNdeta)
                    
                    cont = cont + 1;
                    Int_data(i)%weights_s(cont,1:2) = (/x,y/)
                    Int_data(i)%weights_s(cont,3:5) = N
                    Int_data(i)%weights_s(cont,6:8) = dNdqsi
                    Int_data(i)%weights_s(cont,9:11) = dNdeta
                    Int_data(i)%weights_s(cont,12) = Je*w
                END DO
            ELSEIF (Tri_Quad(j).EQ.'Quad') THEN ! Quad shape function
                DO ii=1,ngp ! Coordinate a1
                    qsi = gp(ii); w1 = gw(ii)
                    DO jj=1,ngp ! Coordinate a2
                        eta = gp(jj); w2 = gw(jj)       
                        !Shape functions of constant quad4 element
                        CALL Shape_functions_quad4(qsi,eta,N3)
                        !Derivative shape functions of constant quad4 element
                        CALL Dev_shape_functions_quad4(qsi,eta,dN3dqsi,dN3deta)
                        ! Jacobian
                        CALL Calc_jac_quad4(X1,X2,X3,X4,dN3dqsi,dN3deta,Je)
                        ! x and y coordinates
                        x = N3(1)*X1(1)+N3(2)*X2(1)+N3(3)*X3(1)+N3(4)*X4(1);
                        y = N3(1)*X1(2)+N3(2)*X2(2)+N3(3)*X3(2)+N3(4)*X4(2);
                        ! shape functions of six-node triangle discontinous element
                        CALL Shape_functions(x,y,N)      
                        ! Derivative of the shape functions
                        CALL Dev_shape_functions(dNdqsi,dNdeta)
                        
                        cont = cont + 1;
                        Int_data(i)%weights_s(cont,1:2) = (/x,y/)
                        Int_data(i)%weights_s(cont,3:5) = N
                        Int_data(i)%weights_s(cont,6:8) = dNdqsi
                        Int_data(i)%weights_s(cont,9:11) = dNdeta
                        Int_data(i)%weights_s(cont,12) = w1*w2*Je
                    END DO
                END DO
            END IF

        END DO

        DEALLOCATE(N_nsc,Tri_Quad)  
    END DO    

    DEALLOCATE(rp,rw)
    
    !---------------------------------------------------------------------------------
    t2 = MPI_Wtime() 
    duration = t2-t1
    CALL MPI_REDUCE(duration,time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,mpierr);
    IF (me.EQ.0) THEN
        WRITE(*,'(A,F15.4,A)') 'COMPLETED! (Time :',time,'s)'
        WRITE (*,*) ''
    END IF
    CALL MPI_Bcast(time,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
    !---------------------------------------------------------------------------------

END SUBROUTINE Integration_Points
!===============================================================================!
SUBROUTINE Matrix_G(N_Phy_nodes)

    REAL(8) :: positions_Phy_Nodes(3,3), N_Phy_nodes(:,:), L(3,3)    
    REAL(8) :: qsi, eta, S_F1(3)
    INTEGER :: M=3, N=3, LDA=3, IPIV(3), LWORK=3, INFO, i, j
    REAL(8) :: WORK(3)

    ! Coordinates of the physical nodes in the parametric element
    positions_Phy_Nodes(1,:) = (/1.0d0, 1.d0-2*lambda, lambda/) ! Node 1
    positions_Phy_Nodes(2,:) = (/2.0d0, lambda, 1.d0-2*lambda/) ! Node 2
    positions_Phy_Nodes(3,:) = (/3.0d0, lambda, lambda/) ! Node 3

    DO i=1,nnos_el !Loop over the six physical nodes
        ! Parametric coordinates
        qsi = positions_Phy_Nodes(i,2); eta = positions_Phy_Nodes(i,3)
        ! Shape functions of three-node triangle continuous element
        CALL Shape_functions_cont(qsi,eta,1,S_F1)
        N_Phy_nodes(3*i-2:3*i,1) = S_F1
    END DO

    ! Matrix G to evaluate the discontinuous shape functions
    DO i=1,nnos_el !loop for shape functions
        DO j=1,nnos_el ! loop for nodes
            qsi = positions_Phy_Nodes(j,2)
            eta = positions_Phy_Nodes(j,3)
            ! shape functions of three-node triangle continous element
            CALL Shape_functions_cont(qsi,eta,1,S_F1)
            L(j,i) = S_F1(i); ! Matrix [L]
        END DO
    END DO

    CALL DGETRF(M,N,L,LDA,IPIV,INFO)
    CALL DGETRI(N,L,LDA,IPIV,WORK,LWORK,INFO)

    Mat_G = L

END SUBROUTINE Matrix_G
!===============================================================================!
SUBROUTINE calc_pts_node(node,N_nsc,Tri_Quad)
    
    INTEGER :: node
    CHARACTER*(*) :: Tri_Quad(:) 
    REAL(8) :: N_nsc(:,:)
    REAL(8) :: X0(2),X1(2),X2(2),X3(2),X4(2),X5(2),X6(2),X7(2),X8(2)
    REAL(8) :: X9(2),X10(2),X11(2),X12(2),X13(2),X14(2),X15(2),X16(2)
    REAL(8) :: X17(2),X18(2),X19(2),X20(2),X21(2),X22(2),X23(2),X24(2)
    REAL(8) :: X25(2)

     X1 =  (/1.d0           ,            0.d0/);   X2 =  (/1.05d0-2*lambda , 2*lambda-0.05d0/); 
     X3 =  (/1.05d0-2*lambda,            0.d0/);   X4 =  (/1.d0-2*lambda   ,          lambda/);
     X5 =  (/0.95d0-2*lambda,            0.d0/);   X6 =  (/0.95d0-2*lambda , 0.05d0+2*lambda/);
     X7 =  (/0.75d0-2*lambda,            0.d0/);   X8 =  (/0.75d0-2*lambda , 2*lambda+0.25d0/);   
     X9 =  (/0.d0           ,            0.d0/);   X10 = (/0.d0            ,            1.d0/);
     X11 = (/2*lambda-0.05d0, 1.05d0-2*lambda/);   X12 = (/0.d0            , 1.05d0-2*lambda/);
     X13 = (/lambda         ,   1.d0-2*lambda/);   X14 = (/0.d0            , 0.95d0-2*lambda/);     
     X15 = (/0.05d0+2*lambda, 0.95d0-2*lambda/);   X16 = (/0.d0            , 0.75d0-2*lambda/);
     X17 = (/2*lambda+0.25d0, 0.75d0-2*lambda/);   X18 = (/lambda          ,          lambda/);
     X19 = (/lambda+0.05d0  ,            0.d0/);   X20 = (/lambda+0.05d0   ,    lambda+0.05d0/);   
     X21 = (/0.d0           ,   lambda+0.05d0/);   X22 = (/lambda+0.25d0   ,            0.d0/);
     X23 = (/lambda+0.25d0  ,   lambda+0.25d0/);   X24 = (/0.d0            ,   lambda+0.25d0/);
     X25 = (/0.5d0          ,           0.5d0/);   X0  = (/0.d0            ,            0.d0/); 

     SELECT CASE(node)

        CASE(1)     
            N_nsc(1,1:8) = (/X1,X2,X3,X0/)
            N_nsc(2,1:8) = (/X4,X4,X3,X2/)
            N_nsc(3,1:8) = (/X4,X4,X5,X3/)
            N_nsc(4,1:8) = (/X4,X4,X6,X5/)
            N_nsc(5,1:8) = (/X4,X4,X6,X2/)
            N_nsc(6,1:8) = (/X7,X5,X6,X8/)
            N_nsc(7,1:8) = (/X9,X7,X8,X10/)
            Tri_Quad(:) = (/'Tria','Quad','Quad','Quad','Quad','Quad','Quad'/) 
        CASE(2) 
            N_nsc(1,1:8) = (/X10,X11,X12,X0/)
            N_nsc(2,1:8) = (/X13,X13,X11,X12/)
            N_nsc(3,1:8) = (/X13,X13,X12,X14/)
            N_nsc(4,1:8) = (/X13,X13,X14,X15/)
            N_nsc(5,1:8) = (/X13,X13,X15,X11/)
            N_nsc(6,1:8) = (/X16,X17,X15,X14/)
            N_nsc(7,1:8) = (/X9,X1,X17,X16/)
            Tri_Quad(:) = (/'Tria','Quad','Quad','Quad','Quad','Quad','Quad'/) 

        CASE(3) 
            N_nsc(1,1:8) = (/X18,X18,X9,X19/)
            N_nsc(2,1:8) = (/X18,X18,X19,X20/)
            N_nsc(3,1:8) = (/X18,X18,X20,X21/)
            N_nsc(4,1:8) = (/X18,X18,X21,X9/)
            N_nsc(5,1:8) = (/X19,X22,X23,X20/)
            N_nsc(6,1:8) = (/X21,X20,X23,X24/)
            N_nsc(7,1:8) = (/X22,X1,X25,X23/)
            N_nsc(8,1:8) = (/X24,X23,X25,X10/)
            Tri_Quad(:) = (/'Quad','Quad','Quad','Quad','Quad','Quad','Quad','Quad'/)

     END SELECT

END SUBROUTINE calc_pts_node
!===============================================================================!
END MODULE Int_points
