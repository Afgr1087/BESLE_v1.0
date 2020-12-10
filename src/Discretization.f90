!===============================================================================!
!-------------------------------------------------------------------------------!
!                     MODULE TO DISCRETIZE THE GEOMETRY
!
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Discretization
!-------------------------------------------------------------------------------!
USE Global_variables
USE Global_functions
USE Set_parameters
USE Failure_Analysis
!USE mpi
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE Discretization_Matrices(ts1)
    
    INTEGER::t0,t1,rate
    REAL(8),INTENT(OUT)::ts1
    
    INTEGER :: nel_reg, i, nnos_total, nelem_global, nnos_global
    INTEGER, ALLOCATABLE :: ELEM_PHY_NODES_reg(:,:), ELEM_GEO_NODES_reg(:,:)
    REAL(8), ALLOCATABLE :: GEO_NODES_reg(:,:), PHY_NODES_reg(:,:)
    
    CALL SYSTEM_CLOCK(t0, rate) ! Initializates the time counter
    WRITE (*,*) '----------------------------------------------------------------------------'
    WRITE (*,*) ' '
    WRITE (*,'(A)',advance='no') ' 3. Discretization matrices        ...'
           
    !POINTS = POINTS*scale_size
  
    nelem_global = SIZE(ELEM,1)
    nnos_global = nelem_global*nnos_el

    ALLOCATE(GEO_PHY_NODES(nreg), GEO_PHY_ELEM(nreg))
    ALLOCATE(GEO_NODES_global(nnos_global,4))
    ALLOCATE(PHY_NODES_global(nnos_global,4))
    ALLOCATE(ELEM_GEO_NODES_global(nelem_global,nnos_el+1))
       
    DO i=1,nreg ! Loops over regions

        nel_reg = El_reg(i,1) ! Number of elemens in regions i
        nnos_total = nel_reg*nnos_el ! Number nodes in the region i

        ALLOCATE(GEO_PHY_NODES(i)%PHY_NODES(nnos_total,4))
        ALLOCATE(GEO_PHY_NODES(i)%GEO_NODES(nnos_total,4))
        ALLOCATE(GEO_PHY_ELEM(i)%ELEM_PHY_NODES(nel_reg,nnos_el+1))
        ALLOCATE(GEO_PHY_ELEM(i)%ELEM_GEO_NODES(nel_reg,nnos_el+1))

        
        CALL Elements_Nodes(i,nel_reg,nnos_total,PHY_NODES_reg, &
                               GEO_NODES_reg,ELEM_PHY_NODES_reg,ELEM_GEO_NODES_reg)
        
        ! Saving a cell with discretizarion data for each region
        GEO_PHY_NODES(i)%PHY_NODES = PHY_NODES_reg       
        GEO_PHY_NODES(i)%GEO_NODES = GEO_NODES_reg 
        GEO_PHY_ELEM(i)%ELEM_PHY_NODES = ELEM_PHY_NODES_reg
        GEO_PHY_ELEM(i)%ELEM_GEO_NODES = ELEM_GEO_NODES_reg

        DEALLOCATE(ELEM_PHY_NODES_reg, ELEM_GEO_NODES_reg)
        DEALLOCATE(GEO_NODES_reg, PHY_NODES_reg)

    END DO     

    CALL SYSTEM_CLOCK(t1)       ! Time counter stopped
    WRITE (*,'(A,F15.3,A)') 'COMPLETED! (Time :', REAL(t1-t0)/rate,'s)'
    WRITE (*,*) ''
    ts1 = REAL(t1-t0)/rate

END SUBROUTINE Discretization_Matrices
!===============================================================================!
SUBROUTINE Elements_Nodes(reg,nel_reg,nnos_total,PHY_NODES_reg, &
                               GEO_NODES_reg,ELEM_PHY_NODES_reg,ELEM_GEO_NODES_reg)

    INTEGER :: i, j, index_geo, index_phy, no1, no2, no3, index_geo_global
    INTEGER :: nel_reg, reg, el, nnos_total, index_elem_global
    INTEGER :: index_elem_geo_global, index_phy_global
    INTEGER, ALLOCATABLE :: ELEM_PHY_NODES_reg(:,:), ELEM_GEO_NODES_reg(:,:)
    REAL(8), ALLOCATABLE :: GEO_NODES_reg(:,:), PHY_NODES_reg(:,:)
    REAL(8) :: X1(3), X2(3), X3(3), x, y, z, N_Phy(3)
    REAL(8) :: X_c(3,3)

	!--------------------------------------------------------------------------
	! Matrix of geometrical and physical nodes
	
    ALLOCATE(GEO_NODES_reg(nnos_total,4), PHY_NODES_reg(nnos_total,4))
    GEO_NODES_reg = 0.d0; PHY_NODES_reg = 0.d0

    ALLOCATE(ELEM_PHY_NODES_reg(nel_reg,nnos_el+1), ELEM_GEO_NODES_reg(nel_reg,nnos_el+1))
    ELEM_PHY_NODES_reg = 0; ELEM_GEO_NODES_reg = 0

    index_geo = 0 ! counter of geometrical nodes
    index_phy = 0 ! counter of physical nodes

    index_geo_global = 0 ! counter of geometrical nodes global
    index_phy_global = 0 ! counter of physical nodes global
    index_elem_global = 0 ! counter of global element
    index_elem_geo_global = 0 ! counter of geometrical nodes global  

    IF (reg.GT.1) THEN
        DO i=1,reg-1
            index_geo_global = index_geo_global + El_reg(i,1)*nnos_el ! counter of geometrical nodes global
            index_phy_global = index_phy_global + El_reg(i,1)*nnos_el
            index_elem_global = index_elem_global + El_reg(i,1) ! counter of global element
            index_elem_geo_global = index_elem_geo_global + El_reg(i,1)*nnos_el ! counter of geometrical nodes global 
        END DO
    END IF
    

    DO i=1,nel_reg

        index_elem_global = index_elem_global + 1        

        el = Subregions(reg,i) ! element

        ! Matrix ELEM_PHY_NODES
        ELEM_PHY_NODES_reg(i,1) = i

        ! Matrix of geometrical nodes
        no1 = ELEM(el,1); no2 = ELEM(el,2); no3 = ELEM(el,3) ! Nodes 1,2 and 3
        ! XYZ coordinates of nodes 1, 2 and 3
        X1 = POINTS(no1,1:3); X2 = POINTS(no2,1:3); X3 = POINTS(no3,1:3);

        ! First three geometrical nodes
        index_geo = index_geo + 1
        GEO_NODES_reg(index_geo,1) = index_geo
        GEO_NODES_reg(index_geo,2:4) = X1

        index_geo = index_geo + 1
        GEO_NODES_reg(index_geo,1) = index_geo
        GEO_NODES_reg(index_geo,2:4) = X2

        index_geo = index_geo + 1
        GEO_NODES_reg(index_geo,1) = index_geo
        GEO_NODES_reg(index_geo,2:4) = X3

        ! First three geometrical nodes global-------------------------
        index_geo_global = index_geo_global + 1
        GEO_NODES_global(index_geo_global,1) = index_geo_global
        GEO_NODES_global(index_geo_global,2:4) = X1

        index_geo_global = index_geo_global + 1
        GEO_NODES_global(index_geo_global,1) = index_geo_global
        GEO_NODES_global(index_geo_global,2:4) = X2

        index_geo_global = index_geo_global + 1
        GEO_NODES_global(index_geo_global,1) = index_geo_global
        GEO_NODES_global(index_geo_global,2:4) = X3
        ! -------------------------------------------------------------

        ! Coordinates
        X_c = GEO_NODES_reg(3*i-2:3*i,2:4)
        X1 = X_c(1,:); X2 = X_c(2,:); X3 = X_c(3,:)

        ! Matrix of physical nodes
        DO j=1,nnos_el !3 new physical nodes
            index_phy = index_phy + 1
            ! Shape functions
            N_Phy = N_Phy_nodes(3*j-2:3*j,1)

            ! x, y, z coordinates of the new physical nodes
            x = N_Phy(1)*X1(1)+N_Phy(2)*X2(1)+N_Phy(3)*X3(1)
            
            y = N_Phy(1)*X1(2)+N_Phy(2)*X2(2)+N_Phy(3)*X3(2)

            z = N_Phy(1)*X1(3)+N_Phy(2)*X2(3)+N_Phy(3)*X3(3)

            ! New physical nodes
            PHY_NODES_reg(index_phy,1) = index_phy
            PHY_NODES_reg(index_phy,2:4) = (/x,y,z/)

            index_phy_global = index_phy_global +1
            PHY_NODES_global(index_phy_global,1) = index_phy_global
            PHY_NODES_global(index_phy_global,2:4) = (/x,y,z/)

            ! Matrix of ELEM_PHY_NODES
            ELEM_PHY_NODES_reg(i,j+1) = index_phy

            ! Matrix of ELEM_GEO_NODES global -------------------------
            index_elem_geo_global = index_elem_geo_global + 1
            ELEM_GEO_NODES_global(index_elem_global,1) = index_elem_global;
            ELEM_GEO_NODES_global(index_elem_global,j+1) = index_elem_geo_global    
            ! -------------------------------------------------------------
        END DO

    END DO

    ELEM_GEO_NODES_reg = ELEM_PHY_NODES_reg(:,1:4)

END SUBROUTINE Elements_Nodes
!===============================================================================!
SUBROUTINE Subregions_Interfaces(nt,me,time)
    include 'mpif.h'
    INTEGER :: max_val, ind1, ind2, i, n_el, j, cont, face_A, face_B, dof
    INTEGER :: nelem, el_A, no11, no21, no31, el_B, no12, no22, no32
    REAL(8) :: x11, y11, z11, x21, y21, z21, x31, y31, z31
    REAL(8) :: x12, y12, z12, x22, y22, z22, x32, y32, z32
    REAL(8) :: Er1x, Er1y, Er1z, Er2x, Er2y, Er2z, Er3x, Er3y, Er3z, Tol
    INTEGER, ALLOCATABLE :: Interfaces_aux(:,:), ELEM_reg(:,:)
    INTEGER :: nelem_global, ind_el_1, ind_el_2, ind_nodes_1, ind_nodes_2, nnos
    INTEGER :: reg_A, reg_B, el, val
    
    REAL(8), ALLOCATABLE :: pointsr2(:), pointsr3(:), pointsr4(:)   
    INTEGER :: me,nt,mpierr, root=0, npoints, scounts(nt), m, n
    INTEGER :: displs(nt), size_int, div
    INTEGER, ALLOCATABLE :: no11v(:), no21v(:), no31v(:), reg_Av(:), face_Av(:)
    INTEGER, ALLOCATABLE :: el_Av(:)
    INTEGER, ALLOCATABLE :: no12v(:), no22v(:), no32v(:), regv(:), facev(:)
    INTEGER, ALLOCATABLE :: Interfaces_aux2(:), Interfaces_aux3(:), Interfaces2(:,:)

	INTEGER :: reg_A_aux, reg_B_aux, cont2, cont1 ,iinterface, Mat_Reg_A, Mat_Reg_B
	INTEGER, ALLOCATABLE :: vec_aux(:), Mat_aux(:,:)

	INTEGER :: no1,no2,no3
	REAL(8) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, d_el(3)
	REAL(8), ALLOCATABLE :: vec_dmin(:)

    REAL(8) :: t1, t2, duration, time
            
    t1 = MPI_Wtime();

    

    CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr);

    IF (me.EQ.0) THEN
        
        IF (nreg.GT.1) THEN

            WRITE (*,*) '----------------------------------------------------------------------------'
            WRITE (*,*) ' '
            WRITE (*,'(A)',advance='no') ' 2. Subregions and Interfaces      '
       
        ELSE

            WRITE (*,*) '----------------------------------------------------------------------------'
            WRITE (*,*) ' '
            WRITE (*,'(A)',advance='no') ' 2. Subregion matrix               '

        END IF

        nelem = SIZE(ELEM,1) ! Total number of elements

		ELEM(:,4) = (/ (i,i=1,nelem) /) ! Each element is a face

        !------------ Subregions ---------------------
        max_val = MAXVAL(El_reg(:,1))
        ALLOCATE(Subregions(nreg,max_val))
        Subregions = 0
        ALLOCATE(ELEM_reg(nelem,2))
        ELEM_reg = 0
        
        ind1 = 0; ind2 = 0
        DO i=1,nreg
            n_el = El_reg(i,1)
            ind1 = ind2 + 1; 
            ind2 = ind2 + n_el
            Subregions(i,1:n_el) = (/ (j,j=ind1,ind2) /) 
            DO j=1,n_el
                el = Subregions(i,j)
                ELEM_reg(el,:) = (/el,i/)
            END DO 
        END DO

		! ---------------------------------------------------------------------------
		! find the minimum segment of an element
		! for each element computes the distance of its three segments and
		! select the minimum. Then, compare and select the minimum over all elements
		

		ALLOCATE(vec_dmin(nelem))
		
		DO i=1,nelem

			no1 = ELEM(i,1); no2 = ELEM(i,2); no3 = ELEM(i,3);
			
			x1 = POINTS(no1,1); y1 = POINTS(no1,2); z1 = POINTS(no1,3);	
			x2 = POINTS(no2,1); y2 = POINTS(no2,2); z2 = POINTS(no2,3);	
			x3 = POINTS(no3,1); y3 = POINTS(no3,2); z3 = POINTS(no3,3);	

			d_el(1) = DSQRT((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
			d_el(2) = DSQRT((x2-x3)**2 + (y2-y3)**2 + (z2-z3)**2)
			d_el(3) = DSQRT((x1-x3)**2 + (y1-y3)**2 + (z1-z3)**2)					
				
			vec_dmin(i) = MINVAL(d_el)
			
		END DO

		dmin = MINVAL(vec_dmin)

    END IF

	CALL MPI_Bcast(dmin,1,MPI_REAL8,root,MPI_COMM_WORLD,mpierr)
    
    !---------------------------------------------
    !-------------Interfaces----------------------
    
    CALL MPI_Bcast(nreg,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)

    IF (nreg.GT.1) THEN ! Just for more than one region

        IF (me.EQ.0) THEN    
        
            ! Global number of elements
            nelem_global = 0
            ind_nodes_2 = 0
            dof = 3

            DO i=1,nreg
                nelem = El_reg(i,1) ! elements per region
                nelem_global = nelem_global + nelem ! global elements
                ind_el_2 = 0
                DO j=(nelem_global-nelem+1),nelem_global
                    ind_el_1 = ind_el_2 + 1
                    ind_el_2 = ind_el_2 + nnos_el*dof
                    ELEM(j,6:7) = (/ind_el_1,ind_el_2/)
                END DO

                nnos = nelem*nnos_el*dof
                ind_nodes_1 = ind_nodes_2 + 1 
                ind_nodes_2 = ind_nodes_2 + nnos      
                
                El_reg(i,3:4) = (/ind_nodes_1,ind_nodes_2/)

            END DO
           
            ! Allocate points to send to all processors     
            npoints = SIZE(POINTS,1)   
            ALLOCATE(pointsr2(npoints))
            pointsr2 = POINTS(:,1)
                    
            ALLOCATE(pointsr3(npoints))
            pointsr3 = POINTS(:,2)
            
            ALLOCATE(pointsr4(npoints))
            pointsr4 = POINTS(:,3)

            !Division of Matrix ELEM  
			IF (MOD(nelem_global,nt).EQ.0) THEN
                scounts = nelem_global/nt;
            ELSE
                m = INT(nelem_global/nt);
                n = nelem_global - m*(nt-1);
                scounts(1:nt-1) = m
                scounts(nt) = n                
            END IF
            val = 0
            displs(1) = 0
            DO i=1,nt-1
                val = val + scounts(i)
                displs(i+1) = val
            END DO
            
            ind2 = 0 ! counter of elements
        END IF

        Tol =  2e-5 * dmin
		
		IF (me.EQ.0) THEN 

			IF (dmin .LT. 1e-15) THEN
				WRITE(*,*) " "
				WRITE(*,*) " "
				WRITE(*,*) "       A MINIMUM SEGMENT DISTANCE GREATER THAN ZERO IS REQUIRED"
				WRITE(*,*) " "				
				call exit(123)

			END IF

		END IF

        ! ELEMENTS B =============================================================
        ! Sending regions and faces to all processors ----------------------------
        CALL MPI_Bcast(scounts,nt,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
        CALL MPI_Bcast(displs,nt,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
        ALLOCATE(regv(scounts(me+1)),facev(scounts(me+1)))    
        CALL MPI_SCATTERV(ELEM_reg(:,2),scounts,displs,MPI_INTEGER,regv, & 
                                 scounts(me+1),MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        
        CALL MPI_SCATTERV(ELEM(:,4),scounts,displs,MPI_INTEGER,facev, & 
                                 scounts(me+1),MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        ! Sending points to all processors ---------------------------------------
        
        CALL MPI_Bcast(nelem_global,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
        CALL MPI_Bcast(npoints,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
        IF (me.NE.0) THEN
            ALLOCATE(pointsr2(npoints),pointsr3(npoints),pointsr4(npoints))  
        END IF     
        CALL MPI_Bcast(pointsr2,npoints,MPI_REAL8,root,MPI_COMM_WORLD,mpierr)
        CALL MPI_Bcast(pointsr3,npoints,MPI_REAL8,root,MPI_COMM_WORLD,mpierr)
        CALL MPI_Bcast(pointsr4,npoints,MPI_REAL8,root,MPI_COMM_WORLD,mpierr)
        ! Sending group of elements to each processor -----------------------------
            
        ALLOCATE(no12v(scounts(me+1)),no22v(scounts(me+1)),no32v(scounts(me+1)))    
        CALL MPI_SCATTERV(ELEM(:,1),scounts,displs,MPI_INTEGER,no12v,scounts(me+1),& 
                                           MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        CALL MPI_SCATTERV(ELEM(:,2),scounts,displs,MPI_INTEGER,no22v,scounts(me+1),& 
                                            MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        CALL MPI_SCATTERV(ELEM(:,3),scounts,displs,MPI_INTEGER,no32v,scounts(me+1),& 
                                            MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
       
        ALLOCATE(Interfaces_aux(nelem_global,7))
        Interfaces_aux = 0
        !========================================================================
        ! ELEMENT A =============================================================
        ! Sending the same element to all processors ----------------------------
        !========================================================================
        div = 1 ! division counter for element A
1       n = scounts(div)
        IF (me.EQ.0) THEN                
            ind1 = ind2 + 1
            ind2 = ind2 + n  
            ALLOCATE(el_Av(n),reg_Av(n),face_Av(n))
            ALLOCATE(no11v(n),no21v(n),no31v(n)) 
            el_Av = (/ (i,i=ind1,ind2) /) 
            reg_Av = ELEM_reg(ind1:ind2,2) ! region A
            face_Av = ELEM(ind1:ind2,4) ! Face A  
            no11v = ELEM(ind1:ind2,1)
            no21v = ELEM(ind1:ind2,2)
            no31v = ELEM(ind1:ind2,3)
        END IF

        IF (me.NE.0) THEN
            ALLOCATE(el_Av(n),reg_Av(n),face_Av(n))
            ALLOCATE(no11v(n),no21v(n),no31v(n))
        END IF

        CALL MPI_Bcast(el_Av,scounts(div),MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
        CALL MPI_Bcast(reg_Av,scounts(div),MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
        CALL MPI_Bcast(face_Av,scounts(div),MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
        CALL MPI_Bcast(no11v,scounts(div),MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
        CALL MPI_Bcast(no21v,scounts(div),MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
        CALL MPI_Bcast(no31v,scounts(div),MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)

        !========================================================================
        ! MATRIX OF INTERFACES ==================================================
        !========================================================================

        DO i = 1,scounts(div)
        
            el_A = el_Av(i)
            reg_A = reg_Av(i);
            face_A = face_Av(i)

            ! Nodes of Element A
            no11 = no11v(i); no21 = no21v(i); no31 = no31v(i);
            ! Coordinates node 1 elemen A
            x11=pointsr2(no11); y11=pointsr3(no11); z11=pointsr4(no11);
            ! Coordinates node 2 elemen A
            x21=pointsr2(no21); y21=pointsr3(no21); z21=pointsr4(no21);
            ! Coordinates node 3 elemen A
            x31=pointsr2(no31); y31=pointsr3(no31); z31=pointsr4(no31);

            ! Comparison between element A and elements B
        
            DO j=1,scounts(me+1)

                IF (nt.GT.1) THEN
                    IF (me.EQ.nt-1) THEN
                        el_B = scounts(me)*(me)+j
                    END IF
                    IF (me.LT.nt-1) THEN
                        el_B = scounts(me+1)*(me)+j  
                    END IF    
                END IF
                IF (nt.EQ.1) THEN
                    el_B = scounts(me+1)*(me)+j        
                END IF  
                reg_B = regv(j) ! region A
                face_B = facev(j) ! Face A

                IF (reg_A.NE.reg_B) THEN

                    ! Nodes of Element B
                    no12=no12v(j); no22=no22v(j); no32=no32v(j);
                    ! Coordinates node 1 elemen B
                    x12=pointsr2(no12); y12=pointsr3(no12); z12=pointsr4(no12);
                    ! Coordinates node 2 elemen B
                    x22=pointsr2(no22); y22=pointsr3(no22); z22=pointsr4(no22);
                    ! Coordinates node 3 elemen B
                    x32=pointsr2(no32); y32=pointsr3(no32); z32=pointsr4(no32);

                    ! Error node 1 in x, y and z
                    Er1x=ABS(x11-x12); Er1y=ABS(y11-y12); Er1z=ABS(z11-z12);
                    ! Error node 2 in x, y and z
                    Er2x=ABS(x21-x22); Er2y=ABS(y21-y22); Er2z=ABS(z21-z22);
                    ! Error node 3 in x, y and z
                    Er3x=ABS(x31-x32); Er3y=ABS(y31-y32); Er3z=ABS(z31-z32);

                    IF(((Er1x.LT.Tol).AND.(Er1y.LT.Tol).AND.(Er1z.LT.Tol)).AND. & 
                       ((Er2x.LT.Tol).AND.(Er2y.LT.Tol).AND.(Er2z.LT.Tol)).AND. & 
                       ((Er3x.LT.Tol).AND.(Er3y.LT.Tol).AND.(Er3z.LT.Tol))) THEN 
                    
                        
                        Interfaces_aux(el_A,1) = el_A ! Number of interface
                        Interfaces_aux(el_A,2) = reg_A ! Region A
                        Interfaces_aux(el_A,3) = reg_B ! Region B
                        Interfaces_aux(el_A,4) = el_A ! Element A
                        Interfaces_aux(el_A,5) = el_B ! Element B
                        Interfaces_aux(el_A,6) = face_A ! face A
                        Interfaces_aux(el_A,7) = face_B ! face B
						
                    END IF                    
                END IF 
            END DO
        END DO
    

        IF (div.LE.nt-1) THEN
            div = div + 1 
            DEALLOCATE(no11v,no21v,no31v)
            DEALLOCATE(reg_Av,face_Av,el_Av)
            GOTO 1 ! Go to the next element
        ELSE
            DEALLOCATE(no11v,no21v,no31v)
            DEALLOCATE(reg_Av,face_Av,el_Av)
            DEALLOCATE(pointsr2,pointsr3,pointsr4)
            DEALLOCATE(no12v,no22v,no32v)
            GOTO 2 ! All processors finish the task
        END IF

        ! Number of interfaces computed by each processor
    2   cont = 0
        DO i=1,nelem_global
            val = Interfaces_aux(i,1)
            IF (val.NE.0) THEN
                cont = cont + 1
            END IF            
        END DO
        
        ! Join the number of interfaces to all processors
        CALL MPI_ALLGATHER(7*cont,1,MPI_INTEGER,scounts,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)

        ! Saving interfaces in a overall vector to be send to me = 0
        ALLOCATE(Interfaces_aux2(cont*7))
        Interfaces_aux2 = 0
        cont = 0
        DO i=1,nelem_global
            val = Interfaces_aux(i,1)
            IF (val.NE.0) THEN
                cont = cont + 1
                Interfaces_aux2(7*cont-6:7*cont) = Interfaces_aux(i,:)  
            END IF            
        END DO 

        DEALLOCATE(Interfaces_aux) 

        ! Plus all interfaces computed by each processor to obtain the total number
        ! interfaces and send to me = 0
        CALL MPI_REDUCE(7*cont,size_int,1,MPI_INTEGER,MPI_SUM,root,MPI_COMM_WORLD,mpierr)
        
        IF (me.EQ.0) THEN
			
            ALLOCATE(Interfaces_aux3(size_int))
            Interfaces_aux3 = 0
        END IF

        ! Send interfaces to me=0
        val = 0
        displs(1) = 0
        DO i=1,nt-1
            val = val + scounts(i)
            displs(i+1) = val
        END DO

        
        CALL MPI_GATHERV(Interfaces_aux2,scounts(me+1),MPI_INTEGER,Interfaces_aux3,scounts,&
                                    displs,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
        
        DEALLOCATE(Interfaces_aux2)

        ! Final matrix of interfaces in the me = 0
        IF (me.EQ.0) THEN
            
            ALLOCATE(Interfaces2(nelem_global,7))
            Interfaces2 = 0
            
            ALLOCATE(Int_reg(nreg,3))

            ind2 = 0
            Int_reg = 0
            DO i=1,size_int/7
                ind1 = ind2 + 1
                ind2 = ind2 + 7
                el_A = Interfaces_aux3(ind1+3)
                
                Interfaces2(el_A,:) = Interfaces_aux3(ind1:ind2)
                
                ! Matrix of interfaces per region
                reg_A = Interfaces2(el_A,2)
                Int_reg(reg_A,1) = reg_A
                cont = Int_reg(reg_A,2)
                cont = cont + 1
                Int_reg(reg_A,2) = cont
				
            END DO
            

            ALLOCATE(Interfaces(size_int/7,20))
			Interfaces = 0
			
            cont = 0
            
            DO i=1,nelem_global
                val = Interfaces2(i,1)
                IF (val.NE.0) THEN
                    cont = cont + 1
                    Interfaces(cont,1) = cont
                    Interfaces(cont,2:7) = Interfaces2(i,2:7)
                    !WRITE(*,*) Interfaces(cont,1:7)

                    el_A = Interfaces(cont,4)
                    el_B = Interfaces(cont,5)
                    ! Indentification of interfaces in the matrix of elements
                    ELEM(el_A,12:13) = (/el_B,cont/)
                END IF
            END DO

            DEALLOCATE(Interfaces2)
			
        END IF

    END IF

	IF (me.EQ.0) THEN
		IF (Crack.EQ.1) THEN
			! Include pre-cracks in the physical model
			CALL Pre_cracks
		    cont = SIZE(Interfaces,1)
		END IF 
	END IF 

	!! Number of TSL in the overall model
	IF (me.EQ.0) THEN
		IF (Failure.EQ.1) THEN

			ALLOCATE(Interfaces2(size_int/7,7))
			Interfaces2 = Interfaces
			ALLOCATE(vec_aux(size_int/7), Mat_aux(size_int/7,size_int/7))
			vec_aux = 0; Mat_aux = 0; cont1=0
			
			DO i=1,cont
				reg_A = Interfaces2(i,2)
				reg_B = Interfaces(i,3)
				vec_aux(1) = i
				
 				cont2 = 0
				DO j=1,cont

					reg_A_aux = Interfaces2(j,2)
					reg_B_aux = Interfaces2(j,3)

					IF (i.NE.j) THEN

						IF (((reg_A.EQ.reg_A_aux).AND.(reg_B.EQ.reg_B_aux)).OR. &
  						    ((reg_A.EQ.reg_B_aux).AND.(reg_B.EQ.reg_A_aux))) THEN

							Interfaces2(j,:) = 0

							cont2 = cont2 + 1

							vec_aux(1+cont2) = j											

						END IF

					END IF
					
				END DO	

				IF (cont2.GT.0) THEN

					cont1 = cont1 + 1
					Mat_aux(cont1,1) = cont2+1
					Mat_aux(cont1,2:cont2+2) = vec_aux(1:cont2+1)

					vec_aux = 0

					DO j=2,cont2+2
						iinterface = Mat_aux(cont1,j)
						Interfaces(iinterface,8) = cont1
					END DO

					iinterface = Mat_aux(cont1,2)
					reg_A = Interfaces2(iinterface,2); 
					reg_B = Interfaces2(iinterface,3);
				
					Mat_Reg_A = Mat_reg(reg_A,2); 
					Mat_Reg_B = Mat_reg(reg_B,2);	
					IF ((Mat_Reg_A.EQ.0).AND.(Mat_Reg_B.EQ.0)) THEN
				    
						adhesive = cont1

				 	END IF

				END IF

			END DO

			ALLOCATE(TSL(cont1))
			n_TSL = cont1

			WRITE(*,*) ''
			WRITE (*,'(A,I6,A)') '     -- Number of Energy criteria: ',n_TSL 
			WRITE(*,*) ''

			DO i=1,cont1
				cont2 = Mat_aux(i,1)
				ALLOCATE(TSL(i)%GBs(cont2))
				TSL(i)%GBs = Mat_aux(i,2:cont2+1)
				
			END DO

			DEALLOCATE(vec_aux,Interfaces2,Mat_aux)

		END IF
	END IF

    !---------------------------------------------------------------------------------
    t2 = MPI_Wtime() 
    duration = t2-t1
    CALL MPI_REDUCE(duration,time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,mpierr);
    IF (me.EQ.0) THEN
  	
        IF ((Failure.EQ.1).OR.(Crack.EQ.1)) THEN

            WRITE(*,'(A,F15.4,A)') '		                   ...COMPLETED! (Time :',time,'s)'
            WRITE (*,*) ''

        ELSE

            WRITE(*,'(A,F15.4,A)') '...COMPLETED! (Time :',time,'s)'
            WRITE (*,*) ''

        END IF

    END IF
    CALL MPI_Bcast(time,1,MPI_DOUBLE,root,MPI_COMM_WORLD,mpierr)
    !---------------------------------------------------------------------------------

END SUBROUTINE Subregions_Interfaces
!===============================================================================!
END MODULE Discretization
