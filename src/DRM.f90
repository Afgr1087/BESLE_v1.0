!===============================================================================!
!-------------------------------------------------------------------------------!
!          MODULE TO COMPUTE THE MATRICES [Û],[T^],[F]                          !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE DRM
!-------------------------------------------------------------------------------!
USE Global_variables
USE Global_functions
!USE mpi
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE compute_Uhat_E(nt,me,Uhat_reg,E_reg,That_reg,reg,C, &
                                                ELEMv,NORMAL_VECTORSv)

    INTEGER :: me, reg, nelem, factor2, factor4, m, nt, sub_nt, n, val, i 
    INTEGER, ALLOCATABLE :: scounts(:), displs(:)
    REAL(8) :: factor1, factor3, nelem_real

    REAL(8), ALLOCATABLE :: Uhat_reg(:,:), E_reg(:,:), That_reg(:,:)
    REAL(8) :: C(6,6)

    INTEGER, ALLOCATABLE :: ELEMv(:)
    REAL(8),ALLOCATABLE :: NORMAL_VECTORSv(:,:) 
	
    IF (me.EQ.0) THEN

		!nelem = El_reg(reg,1)
		!ALLOCATE(scounts(nt), displs(nt))
		!CALL Par_Div_32(nt, nelem, scounts, displs)
		
		!WRITE(*,*) scounts
		!DEALLOCATE(scounts, displs)
		!WRITE(*,*) ''
        
        nelem = El_reg(reg,1) !Number of elements
        nelem_real = El_reg(reg,1)
        !Division of Matrix ELEM ----------------------------------------------------
        factor1 = nelem_real/nt
        factor2 = NINT(factor1)
        IF (factor2.LT.factor1) THEN
            factor2 = factor2 + 1
        END IF
        factor3 = nelem_real/factor2
        factor4 = INT(factor3)*factor2
        m = factor2
        
        IF (MOD(nelem,factor2).LT.(1e-20)) THEN
            sub_nt = INT(factor3)
            ALLOCATE(scounts(sub_nt), displs(sub_nt))
            scounts = m            
        ELSE
            n = nelem - factor4
            sub_nt = INT(factor3) + 1
            ALLOCATE(scounts(sub_nt), displs(sub_nt))
            scounts(1:sub_nt-1) = m
            scounts(sub_nt) = n
            
        END IF

        val = 0
        displs(1) = 0
        DO i=1,sub_nt-1
            val = val + scounts(i)
            displs(i+1) = val
        END DO
		!WRITE(*,*) sub_nt
		!WRITE(*,*) scounts
		
        ! ---------------------------------------------------------------------------
		
        
        ! ---------------------------------------------------------------------------
        
    END IF

    !--------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------
     
    CALL Parallel_Uhat_E(reg,sub_nt,me,scounts,displs,Uhat_reg,E_reg,That_reg, &
                                                                C,ELEMv,NORMAL_VECTORSv)
    
    !--------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------

END SUBROUTINE compute_Uhat_E
!===============================================================================!
SUBROUTINE Parallel_Uhat_E(reg,sub_nt,me,scounts,displs,Uhat_reg,E_reg,That_reg, &
                                                                C,ELEMv,NORMAL_VECTORSv)
    include 'mpif.h'
    INTEGER :: reg, sub_nt, me
    INTEGER, ALLOCATABLE :: scounts(:), displs(:), scounts2(:), displs2(:) 
    REAL(8) :: C(6,6) 
    REAL(8), ALLOCATABLE :: Uhat_reg(:,:), E_reg(:,:), That_reg(:,:)
    REAL(8), ALLOCATABLE :: F_reg(:,:)

    INTEGER :: dof, nelem, nnos
    INTEGER, ALLOCATABLE :: ELEM_PHY_NODES_reg(:,:), Subregionsv(:)
    REAL(8), ALLOCATABLE :: PHY_NODES_reg2(:), PHY_NODES_reg3(:), PHY_NODES_reg4(:)

    INTEGER :: root=0, mpierr, World_group, New_group, NEWCOMM, sub_me, i
    INTEGER, ALLOCATABLE :: ranks(:)

    INTEGER, ALLOCATABLE :: ELEM_PHY_NODES_reg2(:), ELEM_PHY_NODES_reg3(:)
    INTEGER, ALLOCATABLE :: ELEM_PHY_NODES_reg4(:) 
    INTEGER, ALLOCATABLE :: ELEM_PHY_NODES_regv(:,:)
    REAL(8), ALLOCATABLE :: PHY_NODES_regv(:,:)
    REAL(8), ALLOCATABLE :: Uhat_regv(:),F_regv(:),That_regv(:)
    INTEGER, ALLOCATABLE :: spv(:), vec_row(:), vec_col(:)

    INTEGER :: ind1, ind2, div, n, nnosv, sp, elem_fp, j, no_i, no_f, el, face, fp, k
    REAL(8) :: normal(3)
    REAL(8) :: Xk(3), Xj(3), uc(3,3), F(3,3), Tc(3,3), ucv(9), Fv(9), Tcv(9)

    INTEGER, ALLOCATABLE :: ELEMv(:)
    REAL(8),ALLOCATABLE :: NORMAL_VECTORSv(:,:)

    INTEGER :: ind3, ind4, ind5, ind6, ind_ini, ind_fin

    INTEGER :: val, nnosvall, row, col
    INTEGER, ALLOCATABLE :: vec_rowall(:), vec_colall(:)
    REAL(8), ALLOCATABLE :: Uhat_regvall(:),F_regvall(:),That_regvall(:)

	CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr);

    dof = 3 ! Degrees of freedom 
    IF (me.EQ.0) THEN
		
        nelem = El_reg(reg,1) !Number of elements
        nnos = nelem*nnos_el ! Total number of nodes
        
        ! Matrix [Û], [T^] and [F]
        ALLOCATE(Uhat_reg(dof*nnos,dof*nnos))
        Uhat_reg(:,:) = 0.d0
        ALLOCATE(F_reg(dof*nnos,dof*nnos))    
        F_reg(:,:) = 0.d0
        ALLOCATE(That_reg(dof*nnos,dof*nnos))
        That_reg(:,:) = 0.d0

        ! Matrices ELEM_PHY_NODES per region
        ALLOCATE(ELEM_PHY_NODES_reg(nelem,nnos_el))
        ELEM_PHY_NODES_reg = GEO_PHY_ELEM(reg)%ELEM_PHY_NODES(:,2:4)

        ! Matrices PHY_NODES per region
        ALLOCATE(PHY_NODES_reg2(nnos))
        PHY_NODES_reg2 = GEO_PHY_NODES(reg)%PHY_NODES(:,2)
        ALLOCATE(PHY_NODES_reg3(nnos))
        PHY_NODES_reg3 = GEO_PHY_NODES(reg)%PHY_NODES(:,3)
        ALLOCATE(PHY_NODES_reg4(nnos)) 
        PHY_NODES_reg4 = GEO_PHY_NODES(reg)%PHY_NODES(:,4)

        ALLOCATE(Subregionsv(nelem))
        Subregionsv = Subregions(reg,1:nelem)

    END IF

    ! Sending number of elements and nodes to all processors
    CALL MPI_Bcast(nelem,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
    CALL MPI_Bcast(nnos,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)    

    !---------------------------------------------------------------------------------
    ! New group to compute the matrices [Û], [T^] and [F] for each region 

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
        !----------------------------------------------------------------------------
         ! Sending elements in the regions
        IF (sub_me.NE.0) THEN
            ALLOCATE(Subregionsv(nelem))
        END IF
        
        CALL MPI_Bcast(Subregionsv,nelem,MPI_INTEGER,root,NEWCOMM,mpierr) 
            
        !---------------------------------------------------------------------------------
        ! If the number of elements is greater than the number of processors, make the
        ! parallelization for two loops 
        IF (sub_me.NE.0) THEN
            ALLOCATE(scounts(sub_nt), displs(sub_nt))
        END IF
        CALL MPI_Bcast(scounts,sub_nt,MPI_INTEGER,root,NEWCOMM,mpierr)
        CALL MPI_Bcast(displs,sub_nt,MPI_INTEGER,root,NEWCOMM,mpierr)
        
        ! Sending ELEM_PHY_NODES_reg to all processors
        ALLOCATE(ELEM_PHY_NODES_reg2(scounts(sub_me+1)))
        ALLOCATE(ELEM_PHY_NODES_reg3(scounts(sub_me+1)))
        ALLOCATE(ELEM_PHY_NODES_reg4(scounts(sub_me+1)))
           
        CALL MPI_SCATTERV(ELEM_PHY_NODES_reg(:,1),scounts,displs,MPI_INTEGER,ELEM_PHY_NODES_reg2, & 
                                     scounts(sub_me+1),MPI_INTEGER,0,NEWCOMM,mpierr)
        CALL MPI_SCATTERV(ELEM_PHY_NODES_reg(:,2),scounts,displs,MPI_INTEGER,ELEM_PHY_NODES_reg3, & 
                                     scounts(sub_me+1),MPI_INTEGER,0,NEWCOMM,mpierr)
        CALL MPI_SCATTERV(ELEM_PHY_NODES_reg(:,3),scounts,displs,MPI_INTEGER,ELEM_PHY_NODES_reg4, & 
                                     scounts(sub_me+1),MPI_INTEGER,0,NEWCOMM,mpierr)
        
        ! Sending PHY_NODES_reg to all processors
        IF (sub_me.NE.0) THEN
            ALLOCATE(PHY_NODES_reg2(nnos),PHY_NODES_reg3(nnos),PHY_NODES_reg4(nnos))  
        END IF  
        CALL MPI_Bcast(PHY_NODES_reg2,nnos,MPI_REAL8,root,NEWCOMM,mpierr)
        CALL MPI_Bcast(PHY_NODES_reg3,nnos,MPI_REAL8,root,NEWCOMM,mpierr)
        CALL MPI_Bcast(PHY_NODES_reg4,nnos,MPI_REAL8,root,NEWCOMM,mpierr)

        ALLOCATE(ELEM_PHY_NODES_regv(scounts(sub_me+1),3))
        ELEM_PHY_NODES_regv(:,1) = ELEM_PHY_NODES_reg2
        ELEM_PHY_NODES_regv(:,2) = ELEM_PHY_NODES_reg3
        ELEM_PHY_NODES_regv(:,3) = ELEM_PHY_NODES_reg4

        ALLOCATE(PHY_NODES_regv(nnos,3))
        PHY_NODES_regv(:,1) = PHY_NODES_reg2
        PHY_NODES_regv(:,2) = PHY_NODES_reg3
        PHY_NODES_regv(:,3) = PHY_NODES_reg4

        IF (sub_me.EQ.0) THEN
            DEALLOCATE(ELEM_PHY_NODES_reg)
        END IF
        DEALLOCATE (PHY_NODES_reg2,PHY_NODES_reg3, PHY_NODES_reg4)
        DEALLOCATE(ELEM_PHY_NODES_reg2,ELEM_PHY_NODES_reg3,ELEM_PHY_NODES_reg4)
        !---------------------------------------------------------------------------------

        ind2 = 0; ! Counters for sp 
        div = 1 ! division counter for nnos
1       n = nnos_el*scounts(div)
        nnosv = (n*dof)*(nnos_el*scounts(sub_me+1)*dof)
        ALLOCATE(Uhat_regv(nnosv),F_regv(nnosv),That_regv(nnosv))
        ALLOCATE(vec_row(n*dof*nnos_el*scounts(sub_me+1)))
        ALLOCATE(vec_col(n*dof*scounts(sub_me+1)))
            
        ind1 = ind2 + 1
        ind2 = ind2 + n
        ALLOCATE(spv(n))
        spv = (/ (i,i=ind1,ind2) /) ! Source points
        
        ind4 = 0
        ind5 = 0
        ind_fin = 0

        DO i = 1,n ! source points
        
            sp = spv(i)
            Xk = PHY_NODES_regv(sp,1:3)
            
            DO j = 1,scounts(sub_me+1) ! Elements (fielt points)

                IF (sub_nt.GT.1) THEN
                    IF (sub_me.EQ.sub_nt-1) THEN
                        elem_fp = scounts(sub_me)*(sub_me)+j
                    ELSE
                        elem_fp = scounts(sub_me+1)*(sub_me)+j            
                    END IF
                END IF
                IF (sub_nt.EQ.1) THEN
                    elem_fp = scounts(sub_me+1)*(sub_me)+j
                END IF

                no_i = ELEM_PHY_NODES_regv(j,1)  
                no_f = ELEM_PHY_NODES_regv(j,3)          
                
                el = Subregionsv(elem_fp) ! element
                face = ELEMv(el) ! Face of the element
                normal = NORMAL_VECTORSv(face,:) ! Normal vector of the face
                
                ! Vector of rows in the matrix 
                ind3 = ind4 + 1
                ind4 = ind4 + 9
                ind6 = 0
                DO fp = no_i,no_f ! point

                    Xj = PHY_NODES_regv(fp,1:3)
                    
                    ! Matrix [U]
                    CALL compute_Uhat(Xk,Xj,uc)  

                    ! Matrix [F]
                    CALL compute_F(Xk,Xj,F,sp,fp,C)

                    ! Matrix [T]
                    CALL compute_That(Xk,Xj,Tc,normal,sp,fp,C)

                    ! vector format
                    DO k=1,3
                        ucv(3*k-2:3*k) =  uc(k,:)
                        Fv(3*k-2:3*k) = F(k,:)
                        Tcv(3*k-2:3*k) =  Tc(k,:)   
                    END DO    

                    vec_row(ind3+ind6:ind4) = (/3*fp-2,3*fp-1,3*fp/)
                    ind6 = ind6 + 3
                    
                    ! Vector of rows in the matrix 
                    ind5 = ind5 + 1
                    vec_col(ind5) = 3*sp-2

                    ! Vector to be send
                    ind_ini = ind_fin + 1
                    ind_fin = ind_fin + dof*nnos_el 

                    Uhat_regv(ind_ini:ind_fin) = ucv
                    F_regv(ind_ini:ind_fin) = Fv
                    That_regv(ind_ini:ind_fin) = Tcv

                END DO


            END DO
        END DO

        ALLOCATE(scounts2(sub_nt),displs2(sub_nt))
        
        ! Sending vector_row to the overall in the me=0 -----------------------------------------
        CALL MPI_ALLGATHER(ind4,1,MPI_INTEGER,scounts2,1,MPI_INTEGER,NEWCOMM,mpierr)
        
        val = 0
        displs2(1) = 0
        DO i=1,sub_nt-1
            val = val + scounts2(i)
            displs2(i+1) = val
        END DO
        
        CALL MPI_REDUCE(SIZE(vec_row),nnosvall,1,MPI_INTEGER,MPI_SUM,root,NEWCOMM,mpierr)

        IF (sub_me.EQ.0) THEN 
            ALLOCATE(vec_rowall(nnosvall))
        END IF
        
        CALL MPI_GATHERV(vec_row,scounts2(sub_me+1),MPI_INTEGER,vec_rowall,scounts2,&
                                       displs2,MPI_INTEGER,root,NEWCOMM,mpierr)
        
        ! Sending vector_col to the overall in the me=0 -----------------------------------------
        CALL MPI_ALLGATHER(ind5,1,MPI_INTEGER,scounts2,1,MPI_INTEGER,NEWCOMM,mpierr)
            
        val = 0
        displs2(1) = 0
        DO i=1,sub_nt-1
            val = val + scounts2(i)
            displs2(i+1) = val
        END DO

        CALL MPI_REDUCE(SIZE(vec_col),nnosvall,1,MPI_INTEGER,MPI_SUM,root,NEWCOMM,mpierr)
              
        IF (sub_me.EQ.0) THEN 
            ALLOCATE(vec_colall(nnosvall))
        END IF

        CALL MPI_GATHERV(vec_col,scounts2(sub_me+1),MPI_INTEGER,vec_colall,scounts2,&
                                       displs2,MPI_INTEGER,root,NEWCOMM,mpierr)

        ! Sending vector Uhat_regv, F_regv and That_regv to the overall  --------------------
        CALL MPI_ALLGATHER(ind_fin,1,MPI_INTEGER,scounts2,1,MPI_INTEGER,NEWCOMM,mpierr)
            
        val = 0
        displs2(1) = 0
        DO i=1,sub_nt-1
            val = val + scounts2(i)
            displs2(i+1) = val
        END DO

        CALL MPI_REDUCE(nnosv,nnosvall,1,MPI_INTEGER,MPI_SUM,root,NEWCOMM,mpierr)

        IF (sub_me.EQ.0) THEN 
            ALLOCATE(Uhat_regvall(nnosvall))
            ALLOCATE(F_regvall(nnosvall))
            ALLOCATE(That_regvall(nnosvall))
        END IF

        CALL MPI_GATHERV(Uhat_regv,scounts2(sub_me+1),MPI_REAL8,Uhat_regvall,scounts2,&
                                       displs2,MPI_REAL8,root,NEWCOMM,mpierr)
        CALL MPI_GATHERV(F_regv,scounts2(sub_me+1),MPI_REAL8,F_regvall,scounts2,&
                                       displs2,MPI_REAL8,root,NEWCOMM,mpierr)
        CALL MPI_GATHERV(That_regv,scounts2(sub_me+1),MPI_REAL8,That_regvall,scounts2,&
                                       displs2,MPI_REAL8,root,NEWCOMM,mpierr)
        
        !------------------------------------------------------------------------------------------
        ! Fill matrices Uhat_reg, F_reg and That_reg in the me=0
        IF (sub_me.EQ.0) THEN
           ind4 = 0
           DO i=1,SIZE(vec_colall)
                
                col = vec_colall(i)

                DO j=3*i-2,3*i
                    
                    row = vec_rowall(j)
                    
                    ind3 = ind4 + 1
                    ind4 = ind4 + dof
                              
                    Uhat_reg(row,col:col+dof-1) = Uhat_regvall(ind3:ind4)
                    F_reg(row,col:col+dof-1) = F_regvall(ind3:ind4)
                    That_reg(row,col:col+dof-1) = That_regvall(ind3:ind4)

                 END DO
            END DO 
            DEALLOCATE(Uhat_regvall,F_regvall,That_regvall,vec_colall,vec_rowall)   
                 
        END IF

        DEALLOCATE(scounts2, displs2)

        IF (div.LE.sub_nt-1) THEN
            div = div + 1 
            DEALLOCATE(spv,Uhat_regv,F_regv,That_regv,vec_row,vec_col)
            
            GOTO 1 ! Go to the next element
        ELSE
            DEALLOCATE(spv,Uhat_regv,F_regv,That_regv,vec_row,vec_col)
            GOTO 2 ! All processors finish the task
        END IF

        !------------------------------------------------------------------------------------------
        
2       DEALLOCATE(scounts, displs,Subregionsv,ELEM_PHY_NODES_regv,PHY_NODES_regv)
        
        IF (sub_me.EQ.0) THEN
      		
            ! Matrix [E] = INV[F] -------------------------------------------------------!
            ALLOCATE(E_reg(dof*nnos,dof*nnos))    
            E_reg = 0.d0
			!CALL Inv_matrix(F_reg,E_reg)
		END IF
		
        CALL Inv_matrix_i(F_reg,E_reg,sub_me,sub_nt,NEWCOMM)
		
		IF (sub_me.EQ.0) THEN
      		
            DEALLOCATE(F_reg)
                
        END IF

    END IF

    IF (MPI_COMM_NULL .NE. NEWCOMM) THEN
        CALL MPI_COMM_FREE(NEWCOMM, mpierr)
    END IF
        
    CALL MPI_Group_free(World_group, mpierr);
    CALL MPI_Group_free(New_group, mpierr);

    DEALLOCATE(ranks)

END SUBROUTINE Parallel_Uhat_E
!===============================================================================!
SUBROUTINE Parallel_Uhat_E2(Uhat_reg,E_reg,That_reg,reg,C)
    
    INTEGER :: sp, fp, nelem,dof, nnos, reg, face, elem_f, el
    INTEGER :: no_i, no_f
    REAL(8) ::Xk(3), Xj(3), uc(3,3), Tc(3,3)
    REAL(8), ALLOCATABLE :: Uhat_reg(:,:), PHY_NODES_reg(:,:), E_reg(:,:)
    REAL(8), ALLOCATABLE :: That_reg(:,:)
    REAL(8), ALLOCATABLE :: F_reg(:,:)
    REAL(8) :: F(3,3), normal(3)

    REAL(KIND=8) :: C(6,6) ! Compliance tensor

    !INTEGER :: TID, OMP_GET_THREAD_NUM, NTHREADS, OMP_GET_NUM_THREADS

    nelem = El_reg(reg,1) !Number of elements
    dof = 3 ! Degrees of freedom
    nnos = nelem*nnos_el ! Total number of nodes

    ALLOCATE(PHY_NODES_reg(nelem,nnos_el+1))
    PHY_NODES_reg = GEO_PHY_NODES(reg)%PHY_NODES
    
    ! Matrix [Û] and [F]------------------------------------------------------!
    ALLOCATE(Uhat_reg(dof*nnos,dof*nnos))
    Uhat_reg(:,:) = 0.d0
    ALLOCATE(F_reg(dof*nnos,dof*nnos))    
    F_reg(:,:) = 0.d0
    ALLOCATE(That_reg(dof*nnos,dof*nnos))
    That_reg(:,:) = 0.d0
    
    DO sp = 1,nnos ! source points
            
        Xk = PHY_NODES_reg(sp,2:4)
        
        fp = 0
        DO elem_f = 1,nelem ! elements field

            no_i = GEO_PHY_ELEM(reg)%ELEM_PHY_NODES(elem_f,2)  
            no_f = GEO_PHY_ELEM(reg)%ELEM_PHY_NODES(elem_f,4)          
            
            el = Subregions(reg,elem_f) ! element
            face = ELEM(el,4) ! Face of the element
            normal = NORMAL_VECTORS(face,:) ! Normal vector of the face
            
            DO fp = no_i,no_f ! point
                
                Xj = PHY_NODES_reg(fp,2:4)
                
                ! Matrix [U]
                CALL compute_Uhat(Xk,Xj,uc)

                ! Matrix [F]
                CALL compute_F(Xk,Xj,F,sp,fp,C)
        
                ! Matrix [T]
                CALL compute_That(Xk,Xj,Tc,normal,sp,fp,C)

                Uhat_reg(3*fp-2:3*fp,3*sp-2:3*sp) = uc
                
                F_reg(3*fp-2:3*fp,3*sp-2:3*sp) = F
                
                That_reg(3*fp-2:3*fp,3*sp-2:3*sp) = Tc
                
            END DO
            
        END DO
        
    END DO
    
    ! Matrix [E] = INV[F] -------------------------------------------------------!
    ALLOCATE(E_reg(dof*nnos,dof*nnos))    
    E_reg(:,:) = 0.d0
    CALL Inv_matrix(F_reg,E_reg)

    ! Matrix [E] = INV[F] Using G. Sharma et. al.-------------------------------!
    !CALL Inv_matrix2(F_reg,E_reg)  

    DEALLOCATE(PHY_NODES_reg, F_reg)
    
END SUBROUTINE Parallel_Uhat_E2
!===============================================================================!
SUBROUTINE compute_Uhat(Xk,Xj,uc)
  
    REAL(8) :: delta(3,3), Xk(3), Xj(3), R(3), nr, uc(3,3)

    delta(:,:) = 0.d0
    delta(1,1) = 1.d0; delta(2,2) = 1.d0; delta(3,3) = 1.d0

    R = Xj - Xk
    nr = SQRT(R(1)**2 + R(2)**2 + R(3)**2)
    uc = delta*(nr**2 + nr**3) 

END SUBROUTINE compute_Uhat
!===============================================================================!
SUBROUTINE compute_F(Xk,Xj,F,sp,fp,C)

    INTEGER i, n, j, k ,l,sp,fp
    REAL(8) :: delta(3,3), Xk(3), Xj(3), R(3), nr, dr(3), F(3,3)

    REAL(KIND=8) :: C(6,6) ! Compliance tensor

    delta(:,:) = 0.d0
    delta(1,1) = 1.d0; delta(2,2) = 1.d0; delta(3,3) = 1.d0

    R = Xj - Xk
    nr = SQRT(R(1)**2 + R(2)**2 + R(3)**2)    

    IF (sp.NE.fp) THEN 
        dr = R/nr
    ELSE
        dr = (/0.d0, 0.d0, 0.d0/)
    END IF
    F = 0.d0
    DO i=1,3
        DO n=1,3
            DO j=1,3
                DO k=1,3
                    DO l=1,3
                        F(i,n) = F(i,n) + Cmnrs(i,j,k,l,C)*delta(k,n)* &
                        (2.d0*delta(j,l) + 3.d0*nr*delta(j,l) + 3.d0*nr*dr(j)*dr(l))
                    END DO
                END DO
            END DO
        END DO
     END DO

END SUBROUTINE compute_F
!===============================================================================!
SUBROUTINE compute_That(Xk,Xj,Tc,normal,sp,fp,C)
                        
    INTEGER i, n, j, k ,l,sp,fp
    REAL(8) :: delta(3,3), Xk(3), Xj(3), R(3), nr, dr(3), Tc(3,3), normal(3)

    REAL(KIND=8) :: C(6,6) ! Compliance tensor

    R = Xj - Xk
    nr = SQRT(R(1)**2 + R(2)**2 + R(3)**2)    

    delta(:,:) = 0.d0
    delta(1,1) = 1.d0; delta(2,2) = 1.d0; delta(3,3) = 1.d0

    IF (sp.NE.fp) THEN 
        dr = R/nr
    ELSE
        dr = (/0.d0, 0.d0, 0.d0/)
    END IF
    Tc = 0.d0

    DO i=1,3
        DO n=1,3
            DO j=1,3
                DO k=1,3
                    DO l=1,3
                        Tc(i,n) = Tc(i,n) + Cmnrs(i,j,k,l,C)*delta(k,n)* &
                        (2.d0*nr + 3.d0*nr**2)*dr(l)*normal(j)
                    END DO
                END DO
            END DO
        END DO
     END DO

END SUBROUTINE compute_That
!========================================================================================!
FUNCTION Cmnrs(m,n,r,s,C)  ! [1]
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    REAL(8) :: Cmnrs
    INTEGER, INTENT(IN) :: m, n, r, s

    INTEGER::Maux(3,3), val_i, val_j

    REAL(KIND=8) :: C(6,6) ! Compliance tensor

    Maux(1,:) = (/1, 6, 5/)
    Maux(2,:) = (/6, 2, 4/)
    Maux(3,:) = (/5, 4, 3/)

    val_i = Maux(m,n)
    val_j = Maux(r,s)

    Cmnrs = C(val_i,val_j) 

END FUNCTION Cmnrs
!===============================================================================!
SUBROUTINE Inv_matrix_i(Matrix,Matrix2,sub_me,sub_nt,NEWCOMM)

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
		
		!CALL Inv_matrix8(A11,P1)

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

   		!CALL Inv_matrix8(P5,P6)
		
    
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
	
END SUBROUTINE Inv_matrix_i
!===============================================================================!
END MODULE DRM
