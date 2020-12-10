!===============================================================================!
!-------------------------------------------------------------------------------!
!          MODULE TO COMPUTE THE MATRICES [H] AND [G]                           !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Comp_H_G
!-------------------------------------------------------------------------------!
USE Global_variables
USE Global_functions
USE Fourier_Coefficients
USE DRM

!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE Matrices_H_G(nt,me,time)
    
    !INTEGER::t0,t1,rate
    !REAL(8),INTENT(OUT)::ts1
    include 'mpif.h'
    INTEGER :: i, dof, nnos_reg, num, seed, nelem_global, size_reg,j
    REAL(8), ALLOCATABLE :: G_reg(:,:), H_reg(:,:), Uhat_reg(:,:), E_reg(:,:)
    REAL(8), ALLOCATABLE :: That_reg(:,:), M_reg(:,:), Ht_reg(:,:)
    REAL(8), ALLOCATABLE :: GT_reg(:,:), HU_reg(:,:) 
    CHARACTER(LEN=20)::num_str, Test_name

    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIt_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIc_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RI0m_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIm0_cell
    INTEGER, ALLOCATABLE :: pos_vet(:)
    REAL(KIND=8):: R00(3,3)

    INTEGER :: size_coef_cells(4)

    REAL(KIND=8) :: C(6,6) ! Compliance tensor
    REAL(KIND=8) :: S(6,6) ! Stiffness tensor

    INTEGER :: nt,me, mpierr, root=0
    REAL(8) :: t1, t2, duration, time, Anglesv(nreg,4)

    INTEGER :: nfaces, nnormals
    INTEGER, ALLOCATABLE :: ELEMv(:)
    REAL(8), ALLOCATABLE :: Normal_aux(:), NORMAL_VECTORSv(:,:)
    

	

    t1 = MPI_Wtime();

    CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr);

    IF (me.EQ.0) THEN 
        
        WRITE (*,*) '----------------------------------------------------------------------------'
        WRITE (*,*) ' '
        WRITE (*,'(A)',advance='no') ' 5. Matrices [H] and [G]'
        WRITE (*,*) ''    
    
        dof = 3 ! Degrees of freedom

        ALLOCATE(H_G(nreg), Angles(nreg,4))
        IF ((Transient.EQ.1).OR.(Bodyforces.EQ.1)) THEN
            ALLOCATE(U_T_E_M(nreg))
        END IF
        Angles = 0.d0
        CALL SYSTEM_CLOCK(seed)
        CALL SRAND(seed)

        ! Faces 
        nelem_global = SIZE(ELEM,1)
        ALLOCATE(ELEMv(nelem_global))
        ELEMv = ELEM(:,4)

        ! Normal vectors
        nfaces = SIZE(ELEM,1)
        ALLOCATE(Normal_aux(nfaces*3))
        DO i=1,nfaces
            Normal_aux(3*i-2:3*i) = NORMAL_VECTORS(i,:)
        END DO
        nnormals = nfaces*dof
    END IF

    !--------------------------------------------------------------------------
    ! Sending faces
    CALL MPI_Bcast(nelem_global,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
    IF (me.NE.0) THEN
        ALLOCATE(ELEMv(nelem_global))
    END IF
    CALL MPI_Bcast(ELEMv,nelem_global,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
    !--------------------------------------------------------------------------
    ! Sending normal vectors
    CALL MPI_Bcast(nnormals,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
    IF (me.NE.0) THEN
        ALLOCATE(Normal_aux(nnormals))
    END IF
    CALL MPI_Bcast(Normal_aux,nnormals,MPI_REAL8,root,MPI_COMM_WORLD,mpierr) 
    ALLOCATE(NORMAL_VECTORSv(nnormals/3,3)) 
    DO i=1,nnormals/3
        NORMAL_VECTORSv(i,:) = Normal_aux(3*i-2:3*i)
    END DO
    DEALLOCATE(Normal_aux)

    !--------------------------------------------------------------------------
    !===================================================================================
    !===================================================================================
    !===================================================================================
    ! ================= OVER REGIONS ===================================================

    ALLOCATE(C_tot(nreg))

    i=1 ! Loop over regions 

    ! Reading Fourier coefficient from directory----------------------------------------
   
3   IF (me.EQ.0) THEN

        size_reg = El_reg(i,1)*nnos_el*dof

        WRITE (*,*) ''
        WRITE (*,*) ' Region: ', i, '    Elements: ', El_reg(i,1), '    Size: ', size_reg 
        WRITE (*,*) '    Computing matrices H and G ...'
        

        IF (n_Materials.GE.1) THEN 

            num = Mat_reg(i,2)
            
            WRITE(num_str,"(I10)") num
            Test_name = "Mat_"//TRIM(ADJUSTL(num_str))//".dat" 
            
        END IF

    END IF

    CALL MPI_Bcast(num,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
    WRITE(num_str,"(I10)") num
    Test_name = "Mat_"//TRIM(ADJUSTL(num_str))//".dat" 
    Anglesv(i,1) = i

    CALL Fourier_Coeff(me,Test_name,i,RIt_cell, RIc_cell,RI0m_cell, &
                       RIm0_cell,R00,pos_vet,size_coef_cells,C,S,Anglesv)
    
    !-----------------------------------------------------------------------------------!

    ! Compute matrices [H] and [G]------------------------------------------------------!
    
    CALL Compute_H_G(nt,me,ELEMv,NORMAL_VECTORSv,H_reg,G_reg,i,RIt_cell, &
                        RIc_cell, RI0m_cell,RIm0_cell,R00,pos_vet,size_coef_cells,C)

    IF (me.EQ.0) THEN 

        Angles(i,:) = Anglesv(i,:)  
        nnos_reg = El_reg(i,1)*nnos_el            

        ALLOCATE(H_G(i)%H(nnos_reg*dof,nnos_reg*dof))
        ALLOCATE(H_G(i)%G(nnos_reg*dof,nnos_reg*dof))
        H_G(i)%H = 0.d0; H_G(i)%G = 0.d0;

        H_G(i)%G = G_reg 
        IF (Transient.EQ.0) THEN
            H_G(i)%H = H_reg;  
        END IF
        
    END IF 
    
    DEALLOCATE(RIt_cell, RIc_cell, RI0m_cell, RIm0_cell, pos_vet)

    !===================================================================================
    !--------------------------------------------------------------------!  
    ! analysis of body forces / Koegl and Gaul - 2000 (CMES)
    !--------------------------------------------------------------------! 

        IF ((Transient.EQ.1).OR.(Bodyforces.EQ.1)) THEN

            ! Compute matrix [Û], [T^] and [E]-------------------------------!
            IF (me.EQ.0) THEN         
                WRITE (*,*) '    Computing matrices of DRM  ...'
				
            END IF
			
            CALL compute_Uhat_E(nt,me,Uhat_reg,E_reg,That_reg,i,C, &
                                                    ELEMv,NORMAL_VECTORSv)  
            
            ! Compute matrix [M]---------------------------------------------!
            ! M = (GT^ - HÛ)E

            ! First: GT^
            IF (me.EQ.0) THEN
				
                ALLOCATE(GT_reg(nnos_reg*dof,nnos_reg*dof))
				GT_reg = 0.d0 !MATMUL(G_reg,That_reg)
            END IF 
            
            CALL Matmul_parallel2(nt,me,G_reg,That_reg,GT_reg)
            
            ! Second: HÛ
            IF (me.EQ.0) THEN
				
                DEALLOCATE(That_reg)
                ALLOCATE(HU_reg(nnos_reg*dof,nnos_reg*dof))
				HU_reg = 0.d0 !MATMUL(H_reg,Uhat_reg)
            END IF	
			
            CALL Matmul_parallel2(nt,me,H_reg,Uhat_reg,HU_reg)
            
            ! Third: GT^ - HÛ
            IF (me.EQ.0) THEN
				
                DEALLOCATE(Uhat_reg)
                GT_reg = GT_reg-HU_reg
                DEALLOCATE(HU_reg)
                ALLOCATE(M_reg(nnos_reg*dof,nnos_reg*dof))
				M_reg = 0.d0 !MATMUL(GT_reg,E_reg)
            END IF
			
            ! Fourth: M = (GT^ - HÛ)E
            CALL Matmul_parallel2(nt,me,GT_reg,E_reg,M_reg)

            IF (me.EQ.0) THEN
				
                DEALLOCATE(E_reg)
                DEALLOCATE(GT_reg)
                !DEALLOCATE(M_reg)
                
            END IF        
            ! ----------------------------------------------------------------
        END IF

        

        IF (me.EQ.0) THEN

            IF ((Transient.EQ.1).OR.(Bodyforces.EQ.1)) THEN
                
                ! Compute matrix [M]---------------------------------------------!
                ALLOCATE(U_T_E_M(i)%M(nnos_reg*dof,nnos_reg*dof))
                !ALLOCATE(M_reg(nnos_reg*dof,nnos_reg*dof))

                !ALLOCATE(GT_reg(nnos_reg*dof,nnos_reg*dof))
                !ALLOCATE(HU_reg(nnos_reg*dof,nnos_reg*dof))
                
                !CALL Matmul_parallel(G_reg,That_reg,GT_reg)
                   
                !CALL Matmul_parallel(H_reg,Uhat_reg,HU_reg)
                
                !GT_reg = GT_reg-HU_reg
                !CALL Matmul_parallel(GT_reg,E_reg,M_reg)
                       
                !M_reg = MATMUL(MATMUL(G_reg,That_reg)-MATMUL(H_reg,Uhat_reg),E_reg)
                
                !DO j=1,10
                !    WRITE(*,*) M_reg(j,1:3)
                !END DO
                !DEALLOCATE(HU_reg, GT_reg)

                U_T_E_M(i)%M = M_reg
                
                !----------------------------------------------------------------!  

                ! Compute matrix [(2/Dt²)*density*V+H]---------------------------!
                IF (Transient.EQ.1) THEN
                    ALLOCATE(Ht_reg(nnos_reg*dof,nnos_reg*dof))
                    Ht_reg = (2.d0/Dt**2)*density*M_reg + H_reg
                    H_G(i)%H = Ht_reg
                    DEALLOCATE(Ht_reg)
                END IF    

                !DEALLOCATE(Uhat_reg, E_reg, That_reg, M_reg)
                DEALLOCATE(M_reg)
            END IF
            !--------------------------------------------------------------------!   
            DEALLOCATE(H_reg, G_reg)

        END IF
        !===================================================================================

    IF (i.LE.nreg-1) THEN
        i = i + 1 
        GOTO 3 ! Go to the next region
    END IF

    DEALLOCATE(ELEMv,NORMAL_VECTORSv)

    !---------------------------------------------------------------------------------
    t2 = MPI_Wtime() 
    duration = t2-t1
    CALL MPI_REDUCE(duration,time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,mpierr);
    IF (me.EQ.0) THEN
        WRITE (*,*) ''
        WRITE(*,'(A,F15.4,A)') '                                  ...COMPLETED! (Time :',time,'s)'
        WRITE (*,*) ''
    END IF
    CALL MPI_Bcast(time,1,MPI_DOUBLE,root,MPI_COMM_WORLD,mpierr)
    !---------------------------------------------------------------------------------

END SUBROUTINE Matrices_H_G
!===============================================================================!
SUBROUTINE Compute_H_G(nt,me,ELEMv,NORMAL_VECTORSv,H_reg,G_reg,reg,RIt_cell, &
                            RIc_cell, RI0m_cell, RIm0_cell,R00,pos_vet, size_coef_cells,C)
    include 'mpif.h'
    INTEGER :: nelem, dof, nnos, i, reg, me, nt, sub_nt!, section  
    INTEGER, ALLOCATABLE :: scounts(:), displs(:)
    REAL(8), ALLOCATABLE :: G_reg(:,:), H_reg(:,:), NORMAL_VECTORSv(:,:)
    
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIt_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIc_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RI0m_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIm0_cell
    INTEGER, ALLOCATABLE :: pos_vet(:), ELEMv(:)
    REAL(KIND=8):: R00(3,3)
    INTEGER :: size_coef_cells(4)
    REAL(KIND=8) :: C(6,6) ! Compliance tensor

    INTEGER :: factor2, factor4, m, n, val
    REAL(8) :: factor1, factor3, nelem_real

    IF (me.EQ.0) THEN
        
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
		
        ! ---------------------------------------------------------------------------
        
    END IF
    
    !--------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------

    CALL Parallel_compute_H_G(reg,sub_nt,me,scounts,displs,RIt_cell,RIc_cell,RI0m_cell,&
                              RIm0_cell,R00,pos_vet,size_coef_cells,C, &
                              ELEMv,NORMAL_VECTORSv,H_reg,G_reg)

    !--------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------

END SUBROUTINE Compute_H_G
!===============================================================================!
SUBROUTINE Parallel_compute_H_G(reg,sub_nt,me,scounts,displs,RIt_cell,RIc_cell,RI0m_cell,&
                              RIm0_cell,R00,pos_vet,size_coef_cells,C, &
                              ELEMv,NORMAL_VECTORSv,H_reg,G_reg)
    include 'mpif.h'
    INTEGER :: reg, me, dof, nelem, nnos, sub_nt, root=0, mpierr, World_group
    INTEGER :: ind2, n, div, nnosv, ind1, ind4, ind5, ind_fin, i, sp, j, elem_fp, pos
    INTEGER :: ind3, ind_ini, nnosvall, ind6, ind7, col, row
    INTEGER :: New_group, NEWCOMM, sub_me, val

    INTEGER, ALLOCATABLE :: ranks(:)
    INTEGER, ALLOCATABLE :: ELEM_PHY_NODES_reg(:,:), ELEM_PHY_NODES_reg2(:) 
    INTEGER, ALLOCATABLE :: ELEM_PHY_NODES_reg3(:), ELEM_PHY_NODES_reg4(:) 
    INTEGER, ALLOCATABLE :: ELEM_PHY_NODES_regv(:,:)
    INTEGER, ALLOCATABLE :: Subregionsv(:), spv(:), ELEMv(:)
    REAL(8), ALLOCATABLE :: G_reg(:,:), H_reg(:,:)
    REAL(8), ALLOCATABLE :: PHY_NODES_regv(:,:), PHY_NODES_reg2(:), PHY_NODES_reg3(:)
    REAL(8), ALLOCATABLE :: PHY_NODES_reg4(:), NORMAL_VECTORSv(:,:)
    REAL(8) :: g_elv(27), h_elv(27)
    INTEGER, ALLOCATABLE :: vec_row(:), vec_col(:), vec_rowall(:), vec_colall(:)
    REAL(8), ALLOCATABLE :: G_regv(:), H_regv(:), G_regvall(:), H_regvall(:)
    INTEGER, ALLOCATABLE :: scounts(:), displs(:),  scounts2(:), displs2(:) 

    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIt_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIc_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RI0m_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIm0_cell
    INTEGER, ALLOCATABLE :: pos_vet(:)
    REAL(KIND=8):: R00(3,3)
    INTEGER :: size_coef_cells(4)
    REAL(KIND=8) :: C(6,6) ! Compliance tensor

    dof = 3 ! Degrees of freedom 
    IF (me.EQ.0) THEN

        nelem = El_reg(reg,1) !Number of elements
        nnos = nelem*nnos_el ! Total number of nodes
        
        ! Matrices H and G per region
        ALLOCATE(G_reg(dof*nnos,dof*nnos),H_reg(dof*nnos,dof*nnos))
        G_reg(:,:) = 0.d0; H_reg(:,:) = 0.d0

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
    ! New group to compute the matrices [H] and [G] for each region 

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

        DEALLOCATE(ELEM_PHY_NODES_reg2,ELEM_PHY_NODES_reg3,ELEM_PHY_NODES_reg4)
        DEALLOCATE(PHY_NODES_reg2,PHY_NODES_reg3,PHY_NODES_reg4)
           
        ind2 = 0; ! Counters for sp 
        div = 1 ! division counter for nnos
1       n = nnos_el*scounts(div)

        nnosv = (n*dof)*(nnos_el*scounts(sub_me+1)*dof)
        ALLOCATE(G_regv(nnosv),H_regv(nnosv))
        ALLOCATE(vec_row(n*nnos_el*scounts(sub_me+1)))
        ALLOCATE(vec_col(n*scounts(sub_me+1)))
            
        ind1 = ind2 + 1
        ind2 = ind2 + n
        ALLOCATE(spv(n))
        spv = (/ (i,i=ind1,ind2) /) ! Source points
           
        ind4 = 0
        ind5 = 0
        ind_fin = 0
        DO i = 1,n ! source points

            sp = spv(i)

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

                ! If the source point is in the element
                CALL find_vector(ELEM_PHY_NODES_regv(j,1:3),sp,pos) 

                IF(pos.NE.0) THEN ! YES
                    ! Singular integration  
                    CALL H_G_singularv(elem_fp,j,sp,ELEM_PHY_NODES_regv,PHY_NODES_regv,pos,g_elv,h_elv, &
                                           RIt_cell, RIc_cell, RI0m_cell, RIm0_cell,R00, pos_vet, &
                                            size_coef_cells,C,ELEMv,NORMAL_VECTORSv,Subregionsv)     
                ELSE ! NO
                    ! Non-singular integration
                    CALL H_G_Non_singularv(elem_fp,j,sp,ELEM_PHY_NODES_regv,PHY_NODES_regv,g_elv,h_elv, &
                                              RIt_cell, RIc_cell, RI0m_cell, RIm0_cell,R00, pos_vet, &
                                              size_coef_cells,C,ELEMv,NORMAL_VECTORSv,Subregionsv)
                END IF
                    
                ! Vector of rows in the matrix 
                ind3 = ind4 + 1
                ind4 = ind4 + 3
                vec_row(ind3:ind4) = (/3*sp-2,3*sp-1,3*sp/)
        
                ! Vector of columns in the matrix 
                ind5 = ind5 + 1
                vec_col(ind5) = 9*elem_fp-8
                    
                ! Vector to be send
                ind_ini = ind_fin + 1
                ind_fin = ind_fin + 3*nnos_el*dof                 

                G_regv(ind_ini:ind_fin) = g_elv
                H_regv(ind_ini:ind_fin) = h_elv

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

        ! Sending vector G_regv and H_regv to the overall G and H in the me=0 --------------------
        CALL MPI_ALLGATHER(ind_fin,1,MPI_INTEGER,scounts2,1,MPI_INTEGER,NEWCOMM,mpierr)
            
        val = 0
        displs2(1) = 0
        DO i=1,sub_nt-1
            val = val + scounts2(i)
            displs2(i+1) = val
        END DO

        CALL MPI_REDUCE(nnosv,nnosvall,1,MPI_INTEGER,MPI_SUM,root,NEWCOMM,mpierr)

        IF (sub_me.EQ.0) THEN 
            ALLOCATE(G_regvall(nnosvall))
            ALLOCATE(H_regvall(nnosvall))
        END IF

        CALL MPI_GATHERV(G_regv,scounts2(sub_me+1),MPI_REAL8,G_regvall,scounts2,&
                                       displs2,MPI_REAL8,root,NEWCOMM,mpierr)
        CALL MPI_GATHERV(H_regv,scounts2(sub_me+1),MPI_REAL8,H_regvall,scounts2,&
                                       displs2,MPI_REAL8,root,NEWCOMM,mpierr)
        !------------------------------------------------------------------------------------------
        ! Fill matrices H_reg and G_reg in the me=0
        IF (sub_me.EQ.0) THEN
            ind7 = 0
            DO i=1,SIZE(vec_colall)

                col = vec_colall(i)
                    
                DO j=3*i-2,3*i
                    row = vec_rowall(j)
                        
                    ind6 = ind7 + 1
                    ind7 = ind7 + nnos_el*dof
                        
                    G_reg(row,col:col+nnos_el*dof-1) = G_regvall(ind6:ind7)
                    H_reg(row,col:col+nnos_el*dof-1) = H_regvall(ind6:ind7)

                END DO
            END DO 
            DEALLOCATE(G_regvall,H_regvall,vec_rowall,vec_colall)         
        END IF

        DEALLOCATE(scounts2,displs2)
        !------------------------------------------------------------------------------------------

        IF (div.LE.sub_nt-1) THEN
            div = div + 1 
            DEALLOCATE(spv)
            DEALLOCATE(G_regv,H_regv,vec_col,vec_row)
            GOTO 1 ! Go to the next element
        ELSE
            DEALLOCATE(spv)
            DEALLOCATE(G_regv,H_regv,vec_col,vec_row)
            GOTO 2 ! All processors finish the task
        END IF

2       DEALLOCATE(scounts,displs,ELEM_PHY_NODES_regv,PHY_NODES_regv,Subregionsv) 

        !Diagonal elements in matrix [H]
        IF (sub_me.EQ.0) THEN

            DO i=1,nnos
                H_reg(3*i-2:3*i,3*i-2:3*i) = 0.d0
                DO j=1,nnos
                    IF (j.NE.i) THEN
                    H_reg(3*i-2:3*i,3*i-2:3*i) = H_reg(3*i-2:3*i,3*i-2:3*i) &
                                                       & - H_reg(3*i-2:3*i,3*j-2:3*j)
                    END IF
                END DO
            END DO

            DEALLOCATE(ELEM_PHY_NODES_reg)

        END IF

        !---------------------------------------------------------------------------------

    END IF

    IF (MPI_COMM_NULL .NE. NEWCOMM) THEN
        CALL MPI_COMM_FREE(NEWCOMM, mpierr)
    END IF
        
    CALL MPI_Group_free(World_group, mpierr);
    CALL MPI_Group_free(New_group, mpierr);

    DEALLOCATE(ranks)

END SUBROUTINE Parallel_compute_H_G
!===============================================================================!
SUBROUTINE H_G_singularv(elem_fp_g,elem_fp,sp,ELEM_PHY_NODES_reg,PHY_NODES_reg,pos,g_el,h_el, &
                                  RIt_cell, RIc_cell, RI0m_cell, RIm0_cell, R00, pos_vet, &
                                  size_coef_cells,C,ELEMv,NORMAL_VECTORSv,Subregionsv) 

    INTEGER :: elem_fp, sp, pos, node1, node2, node3, face
    INTEGER :: npts, i, el, elem_fp_g
    INTEGER, ALLOCATABLE :: ELEM_PHY_NODES_reg(:,:)
    REAL(8) :: X1(3), X2(3), X3(3), g_el(27), XD(3), n(3)
    REAL(8) :: wgts(N_int_points(pos)), N1(N_int_points(pos)), h_el(27)
    REAL(8) :: N2(N_int_points(pos)), N3(N_int_points(pos))
    REAL(8) :: dN1dqsi(N_int_points(pos))
    REAL(8) :: dN2dqsi(N_int_points(pos)), dN3dqsi(N_int_points(pos))
    REAL(8) :: dN1deta(N_int_points(pos))
    REAL(8) :: dN2deta(N_int_points(pos)), dN3deta(N_int_points(pos))
    REAL(8) :: xc(N_int_points(pos))
    REAL(8) :: yc(N_int_points(pos)), zc(N_int_points(pos))
    REAL(8) :: dxdqsi(N_int_points(pos)), dydqsi(N_int_points(pos))
    REAL(8) :: dzdqsi(N_int_points(pos)), dxdeta(N_int_points(pos))
    REAL(8) :: dydeta(N_int_points(pos)), dzdeta(N_int_points(pos))
    REAL(8) :: g1(N_int_points(pos)), g2(N_int_points(pos))
    REAL(8) :: g3(N_int_points(pos)), J(N_int_points(pos))
    REAL(8) :: uast11(N_int_points(pos)), uast12(N_int_points(pos))
    REAL(8) :: uast22(N_int_points(pos)), uast33(N_int_points(pos))
    REAL(8) :: uast13(N_int_points(pos)), uast23(N_int_points(pos))
    REAL(8) :: tast11(N_int_points(pos)), tast12(N_int_points(pos))
    REAL(8) :: tast13(N_int_points(pos)), tast21(N_int_points(pos))
    REAL(8) :: tast22(N_int_points(pos)), tast23(N_int_points(pos))
    REAL(8) :: tast31(N_int_points(pos)), tast32(N_int_points(pos))
    REAL(8) :: tast33(N_int_points(pos)), xsp(3), xfp(3), theta, phi, r
    REAL(8) :: Uuv(3,3), Tij(3,3), g11_1, g12_1, g22_1, g33_1, g13_1, g23_1
    REAL(8) :: g11_2, g12_2, g22_2, g33_2, g13_2, g23_2, g11_3, g12_3, g22_3
    REAL(8) :: g33_3, g13_3, g23_3, h11_1, h12_1, h13_1, h21_1, h22_1, h23_1
    REAL(8) :: h31_1, h32_1, h33_1, h11_2, h12_2, h13_2, h21_2, h22_2, h23_2
    REAL(8) :: h31_2, h32_2, h33_2, h11_3, h12_3, h13_3, h21_3, h22_3, h23_3
    REAL(8) :: h31_3, h32_3, h33_3!, Error
    REAL(8), ALLOCATABLE :: PHY_NODES_reg(:,:)

    INTEGER :: size_coef_cells(4)

    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIt_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIc_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RI0m_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIm0_cell
    INTEGER, ALLOCATABLE :: pos_vet(:)
    REAL(KIND=8):: R00(3,3)

    REAL(KIND=8) :: C(6,6) ! Compliance tensor

    INTEGER, ALLOCATABLE :: Subregionsv(:), ELEMv(:)
    REAL(8),ALLOCATABLE :: NORMAL_VECTORSv(:,:)

    g_el = 0.d0; h_el = 0.d0 ! Initializating [g] and [h]
    el = Subregionsv(elem_fp_g) ! element
    face = ELEMv(el) ! Face of the element
    n = NORMAL_VECTORSv(face,:) ! Normal vector of the face
    
    ! Nodes of the element that contains the field point
    node1 = ELEM_PHY_NODES_reg(elem_fp,1); 
    node2 = ELEM_PHY_NODES_reg(elem_fp,2);
    node3 = ELEM_PHY_NODES_reg(elem_fp,3); 

    ! Coordinates 
    X1 = PHY_NODES_reg(node1,1:3); 
    X2 = PHY_NODES_reg(node2,1:3); 
    X3 = PHY_NODES_reg(node3,1:3); 

    XD = PHY_NODES_reg(sp,1:3) ! Source point

    ! Weights for the numerical integration
    wgts = Int_data(pos)%weights_s(:,12)

    ! Shape functions of discontinuous element
    N1 = Int_data(pos)%weights_s(:,3)
    N2 = Int_data(pos)%weights_s(:,4)
    N3 = Int_data(pos)%weights_s(:,5)
    
    ! Derivatives of the shape functions of discontinuous elements
    dN1dqsi = Int_data(pos)%weights_s(:,6);  
    dN2dqsi = Int_data(pos)%weights_s(:,7);
    dN3dqsi = Int_data(pos)%weights_s(:,8); 

    dN1deta = Int_data(pos)%weights_s(:,9); 
    dN2deta = Int_data(pos)%weights_s(:,10);
    dN3deta = Int_data(pos)%weights_s(:,11); 

    ! Coordinates of the field points 
    xc = N1*X1(1) + N2*X2(1) + N3*X3(1) 
    yc = N1*X1(2) + N2*X2(2) + N3*X3(2) 
    zc = N1*X1(3) + N2*X2(3) + N3*X3(3) 

    ! Jacobian
    dxdqsi = X1(1)*dN1dqsi + X2(1)*dN2dqsi + X3(1)*dN3dqsi 
    dydqsi = X1(2)*dN1dqsi + X2(2)*dN2dqsi + X3(2)*dN3dqsi 
    dzdqsi = X1(3)*dN1dqsi + X2(3)*dN2dqsi + X3(3)*dN3dqsi 

    dxdeta = X1(1)*dN1deta + X2(1)*dN2deta + X3(1)*dN3deta 
    dydeta = X1(2)*dN1deta + X2(2)*dN2deta + X3(2)*dN3deta
    dzdeta = X1(3)*dN1deta + X2(3)*dN2deta + X3(3)*dN3deta 

    g1 = dydqsi*dzdeta - dzdqsi*dydeta
    g2 = dzdqsi*dxdeta - dxdqsi*dzdeta
    g3 = dxdqsi*dydeta - dydqsi*dxdeta

    J = DSQRT(g1**2 + g2**2 + g3**2)

    ! Initializing u and t matrices of fundamental solutions

    npts = N_int_points(pos)
    uast11 = 0.d0; uast12 = 0.d0
    uast22 = 0.d0; uast33 = 0.d0
    uast13 = 0.d0; uast23 = 0.d0

    tast11 = 0.d0; tast12 = 0.d0; tast13 = 0.d0
    tast21 = 0.d0; tast22 = 0.d0; tast23 = 0.d0
    tast31 = 0.d0; tast32 = 0.d0; tast33 = 0.d0
    
    xsp = XD

    DO i=1,npts

        xfp = (/xc(i),yc(i),zc(i)/)
        CALL spherical_coordinates(xsp,xfp,theta,phi,r)
        
        Uuv = 0.d0; Tij = 0.d0
        CALL Fund_Solution_UT(theta,phi,r,n,Uuv,Tij, &
                              RIt_cell, RIc_cell, RI0m_cell, RIm0_cell,R00, pos_vet, &
                              size_coef_cells,C)

        uast11(i) = Uuv(1,1); uast12(i) = Uuv(1,2)
        uast22(i) = Uuv(2,2); uast33(i) = Uuv(3,3)
        uast13(i) = Uuv(1,3); uast23(i) = Uuv(2,3)
        
        tast11(i) = Tij(1,1); tast12(i) = Tij(1,2)
        tast13(i) = Tij(1,3); tast21(i) = Tij(2,1)
        tast22(i) = Tij(2,2); tast23(i) = Tij(2,3)
        tast31(i) = Tij(3,1); tast32(i) = Tij(3,2)
        tast33(i) = Tij(3,3);
            
    END DO
    
    ![G]j: = [ [G]j1, [G]j2, [G]j3, [G]j4, [G]j5, [G]j6 ]

    ! [G]j1 components
    g11_1 = SUM(uast11*wgts*J*N1); g12_1 = SUM(uast12*wgts*J*N1) 
    g22_1 = SUM(uast22*wgts*J*N1); g33_1 = SUM(uast33*wgts*J*N1) 
    g13_1 = SUM(uast13*wgts*J*N1); g23_1 = SUM(uast23*wgts*J*N1)
    ! [G]j2 components
    g11_2 = SUM(uast11*wgts*J*N2); g12_2 = SUM(uast12*wgts*J*N2) 
    g22_2 = SUM(uast22*wgts*J*N2); g33_2 = SUM(uast33*wgts*J*N2) 
    g13_2 = SUM(uast13*wgts*J*N2); g23_2 = SUM(uast23*wgts*J*N2)
    ! [G]j3 components
    g11_3 = SUM(uast11*wgts*J*N3); g12_3 = SUM(uast12*wgts*J*N3) 
    g22_3 = SUM(uast22*wgts*J*N3); g33_3 = SUM(uast33*wgts*J*N3) 
    g13_3 = SUM(uast13*wgts*J*N3); g23_3 = SUM(uast23*wgts*J*N3)

    ! Numeric integration of [g]
    g_el(1:9) = (/g11_1,g12_1,g13_1,g11_2,g12_2,g13_2,g11_3,g12_3,g13_3/)
    g_el(10:18) = (/g12_1,g22_1,g23_1,g12_2,g22_2,g23_2,g12_3,g22_3,g23_3/)
    g_el(19:27) = (/g13_1,g23_1,g33_1,g13_2,g23_2,g33_2,g13_3,g23_3,g33_3/)

    ![H]j: = [ [H]j1, [H]j2, [H]j3, [H]j4, [H]j5, [H]j6 ]
    
    ! [H]j1 components
    h11_1 = SUM(tast11*wgts*J*N1); h12_1 = SUM(tast12*wgts*J*N1); h13_1 = SUM(tast13*wgts*J*N1) 
    h21_1 = SUM(tast21*wgts*J*N1); h22_1 = SUM(tast22*wgts*J*N1); h23_1 = SUM(tast23*wgts*J*N1)
    h31_1 = SUM(tast31*wgts*J*N1); h32_1 = SUM(tast32*wgts*J*N1); h33_1 = SUM(tast33*wgts*J*N1)
    ! [H]j2 components
    h11_2 = SUM(tast11*wgts*J*N2); h12_2 = SUM(tast12*wgts*J*N2); h13_2 = SUM(tast13*wgts*J*N2) 
    h21_2 = SUM(tast21*wgts*J*N2); h22_2 = SUM(tast22*wgts*J*N2); h23_2 = SUM(tast23*wgts*J*N2)
    h31_2 = SUM(tast31*wgts*J*N2); h32_2 = SUM(tast32*wgts*J*N2); h33_2 = SUM(tast33*wgts*J*N2)
    ! [H]j3 components
    h11_3 = SUM(tast11*wgts*J*N3); h12_3 = SUM(tast12*wgts*J*N3); h13_3 = SUM(tast13*wgts*J*N3) 
    h21_3 = SUM(tast21*wgts*J*N3); h22_3 = SUM(tast22*wgts*J*N3); h23_3 = SUM(tast23*wgts*J*N3)
    h31_3 = SUM(tast31*wgts*J*N3); h32_3 = SUM(tast32*wgts*J*N3); h33_3 = SUM(tast33*wgts*J*N3)

    ! Numeric integration of [h]
    h_el(1:9) = (/h11_1,h12_1,h13_1,h11_2,h12_2,h13_2,h11_3,h12_3,h13_3/)
    h_el(10:18) = (/h21_1,h22_1,h23_1,h21_2,h22_2,h23_2,h21_3,h22_3,h23_3/)
    h_el(19:27) = (/h31_1,h32_1,h33_1,h31_2,h32_2,h33_2,h31_3,h32_3,h33_3/)

END SUBROUTINE H_G_singularv
!===============================================================================!
SUBROUTINE H_G_Non_singularv(elem_fp_g,elem_fp,sp,ELEM_PHY_NODES_reg,PHY_NODES_reg,g_el,h_el, &
                                  RIt_cell, RIc_cell, RI0m_cell, RIm0_cell, R00, pos_vet, &
                                  size_coef_cells,C,ELEMv,NORMAL_VECTORSv,Subregionsv)

    INTEGER :: elem_fp, sp, node1, node2, node3, face
    INTEGER :: npts, i, el, elem_fp_g
    INTEGER, ALLOCATABLE :: ELEM_PHY_NODES_reg(:,:)
    REAL(8) :: X1(3), X2(3), X3(3), g_el(27), XD(3), n(3)
    REAL(8) :: h_el(27), wgts(npoints_Nsi), N1(npoints_Nsi), N2(npoints_Nsi)
    REAL(8) :: N3(npoints_Nsi)
    REAL(8) :: dN1dqsi(npoints_Nsi), dN2dqsi(npoints_Nsi)
    REAL(8) :: dN3dqsi(npoints_Nsi)
    REAL(8) :: dN1deta(npoints_Nsi), dN2deta(npoints_Nsi)
    REAL(8) :: dN3deta(npoints_Nsi)
    REAL(8) :: xc(npoints_Nsi), yc(npoints_Nsi)
    REAL(8) :: zc(npoints_Nsi), dxdqsi(npoints_Nsi), dydqsi(npoints_Nsi)
    REAL(8) :: dzdqsi(npoints_Nsi), dxdeta(npoints_Nsi), dydeta(npoints_Nsi)
    REAL(8) :: dzdeta(npoints_Nsi), g1(npoints_Nsi), g2(npoints_Nsi)
    REAL(8) :: g3(npoints_Nsi), J(npoints_Nsi), uast11(npoints_Nsi)
    REAL(8) :: uast12(npoints_Nsi), uast22(npoints_Nsi), uast33(npoints_Nsi)
    REAL(8) :: uast13(npoints_Nsi), uast23(npoints_Nsi), tast11(npoints_Nsi)
    REAL(8) :: tast12(npoints_Nsi), tast13(npoints_Nsi), tast21(npoints_Nsi)
    REAL(8) :: tast22(npoints_Nsi), tast23(npoints_Nsi), tast31(npoints_Nsi)
    REAL(8) :: tast32(npoints_Nsi), tast33(npoints_Nsi), xsp(3), xfp(3) 
    REAL(8) :: theta, phi, r, Uuv(3,3), Tij(3,3)
    REAL(8) :: g11_1, g12_1, g22_1, g33_1, g13_1, g23_1
    REAL(8) :: g11_2, g12_2, g22_2, g33_2, g13_2, g23_2, g11_3, g12_3, g22_3
    REAL(8) :: g33_3, g13_3, g23_3, h11_1, h12_1, h13_1, h21_1, h22_1, h23_1
    REAL(8) :: h31_1, h32_1, h33_1, h11_2, h12_2, h13_2, h21_2, h22_2, h23_2
    REAL(8) :: h31_2, h32_2, h33_2, h11_3, h12_3, h13_3, h21_3, h22_3, h23_3
    REAL(8) :: h31_3, h32_3, h33_3!, Error
    REAL(8), ALLOCATABLE :: PHY_NODES_reg(:,:)

    INTEGER :: size_coef_cells(4)

    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIt_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIc_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RI0m_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIm0_cell
    INTEGER, ALLOCATABLE :: pos_vet(:)
    REAL(KIND=8):: R00(3,3)    

    REAL(KIND=8) :: C(6,6) ! Compliance tensor

    INTEGER, ALLOCATABLE :: Subregionsv(:), ELEMv(:)
    REAL(8),ALLOCATABLE :: NORMAL_VECTORSv(:,:)

    g_el = 0.d0; h_el=0.d0 ! Initializating [g] and [h]
    el = Subregionsv(elem_fp_g) ! element 
    face = ELEMv(el) ! Face of the element
    n = NORMAL_VECTORSv(face,:) ! Normal vector of the face

    ! Nodes of the element that contains the field point
    node1 = ELEM_PHY_NODES_reg(elem_fp,1); 
    node2 = ELEM_PHY_NODES_reg(elem_fp,2);
    node3 = ELEM_PHY_NODES_reg(elem_fp,3);
    ! Coordinates 
    X1 = PHY_NODES_reg(node1,1:3); 
    X2 = PHY_NODES_reg(node2,1:3); 
    X3 = PHY_NODES_reg(node3,1:3)

    XD = PHY_NODES_reg(sp,1:3) ! Source point

    ! Weight for numeric integration 
    wgts = weights_ns(:,12)
    
    ! Shape functions of discontinuous element
    N1 = weights_ns(:,3); N2 = weights_ns(:,4); N3 = weights_ns(:,5);
    
    ! Derivatives of the shape functions of discontinous element
    dN1dqsi = weights_ns(:,6);  
    dN2dqsi = weights_ns(:,7);
    dN3dqsi = weights_ns(:,8); 

    dN1deta = weights_ns(:,9); 
    dN2deta = weights_ns(:,10)
    dN3deta = weights_ns(:,11);

    ! Coordinates of the field points 
    xc = N1*X1(1) + N2*X2(1) + N3*X3(1)
    yc = N1*X1(2) + N2*X2(2) + N3*X3(2) 
    zc = N1*X1(3) + N2*X2(3) + N3*X3(3) 

    ! Jacobian
    dxdqsi = X1(1)*dN1dqsi + X2(1)*dN2dqsi + X3(1)*dN3dqsi 
    dydqsi = X1(2)*dN1dqsi + X2(2)*dN2dqsi + X3(2)*dN3dqsi 
    dzdqsi = X1(3)*dN1dqsi + X2(3)*dN2dqsi + X3(3)*dN3dqsi 

    dxdeta = X1(1)*dN1deta + X2(1)*dN2deta + X3(1)*dN3deta 
    dydeta = X1(2)*dN1deta + X2(2)*dN2deta + X3(2)*dN3deta 
    dzdeta = X1(3)*dN1deta + X2(3)*dN2deta + X3(3)*dN3deta 

    g1 = dydqsi*dzdeta - dzdqsi*dydeta
    g2 = dzdqsi*dxdeta - dxdqsi*dzdeta
    g3 = dxdqsi*dydeta - dydqsi*dxdeta

    J = DSQRT(g1**2 + g2**2 + g3**2)

    ! Initializing u and t matrices of fundamental solutions
    npts = npoints_Nsi
    uast11 = 0.d0; uast12 = 0.d0
    uast22 = 0.d0; uast33 = 0.d0
    uast13 = 0.d0; uast23 = 0.d0

    tast11 = 0.d0; tast12 = 0.d0; tast13 = 0.d0
    tast21 = 0.d0; tast22 = 0.d0; tast23 = 0.d0
    tast31 = 0.d0; tast32 = 0.d0; tast33 = 0.d0
    
    xsp = XD
    
    DO i=1,npts

        xfp = (/xc(i),yc(i),zc(i)/)
        CALL spherical_coordinates(xsp,xfp,theta,phi,r)

        Uuv = 0.d0; Tij = 0.d0
        CALL Fund_Solution_UT(theta,phi,r,n,Uuv,Tij, &
                              RIt_cell, RIc_cell, RI0m_cell, RIm0_cell, R00, pos_vet, &
                              size_coef_cells,C)

        uast11(i) = Uuv(1,1); uast12(i) = Uuv(1,2)
        uast22(i) = Uuv(2,2); uast33(i) = Uuv(3,3)
        uast13(i) = Uuv(1,3); uast23(i) = Uuv(2,3)
        
        tast11(i) = Tij(1,1); tast12(i) = Tij(1,2)
        tast13(i) = Tij(1,3); tast21(i) = Tij(2,1)
        tast22(i) = Tij(2,2); tast23(i) = Tij(2,3)
        tast31(i) = Tij(3,1); tast32(i) = Tij(3,2)
        tast33(i) = Tij(3,3);
            
    END DO
    
     ![G]j: = [ [G]j1, [G]j2, [G]j3, [G]j4, [G]j5, [G]j6 ]

    ! [G]j1 components
    g11_1 = SUM(uast11*wgts*J*N1); g12_1 = SUM(uast12*wgts*J*N1) 
    g22_1 = SUM(uast22*wgts*J*N1); g33_1 = SUM(uast33*wgts*J*N1) 
    g13_1 = SUM(uast13*wgts*J*N1); g23_1 = SUM(uast23*wgts*J*N1)
    ! [G]j2 components
    g11_2 = SUM(uast11*wgts*J*N2); g12_2 = SUM(uast12*wgts*J*N2) 
    g22_2 = SUM(uast22*wgts*J*N2); g33_2 = SUM(uast33*wgts*J*N2) 
    g13_2 = SUM(uast13*wgts*J*N2); g23_2 = SUM(uast23*wgts*J*N2)
    ! [G]j3 components
    g11_3 = SUM(uast11*wgts*J*N3); g12_3 = SUM(uast12*wgts*J*N3) 
    g22_3 = SUM(uast22*wgts*J*N3); g33_3 = SUM(uast33*wgts*J*N3) 
    g13_3 = SUM(uast13*wgts*J*N3); g23_3 = SUM(uast23*wgts*J*N3)

    ! Numeric integration of [g]
    g_el(1:9) = (/g11_1,g12_1,g13_1,g11_2,g12_2,g13_2,g11_3,g12_3,g13_3/)
    g_el(10:18) = (/g12_1,g22_1,g23_1,g12_2,g22_2,g23_2,g12_3,g22_3,g23_3/)
    g_el(19:27) = (/g13_1,g23_1,g33_1,g13_2,g23_2,g33_2,g13_3,g23_3,g33_3/)

    ![H]j: = [ [H]j1, [H]j2, [H]j3, [H]j4, [H]j5, [H]j6 ]
    
    ! [H]j1 components
    h11_1 = SUM(tast11*wgts*J*N1); h12_1 = SUM(tast12*wgts*J*N1); h13_1 = SUM(tast13*wgts*J*N1) 
    h21_1 = SUM(tast21*wgts*J*N1); h22_1 = SUM(tast22*wgts*J*N1); h23_1 = SUM(tast23*wgts*J*N1)
    h31_1 = SUM(tast31*wgts*J*N1); h32_1 = SUM(tast32*wgts*J*N1); h33_1 = SUM(tast33*wgts*J*N1)
    ! [H]j2 components
    h11_2 = SUM(tast11*wgts*J*N2); h12_2 = SUM(tast12*wgts*J*N2); h13_2 = SUM(tast13*wgts*J*N2) 
    h21_2 = SUM(tast21*wgts*J*N2); h22_2 = SUM(tast22*wgts*J*N2); h23_2 = SUM(tast23*wgts*J*N2)
    h31_2 = SUM(tast31*wgts*J*N2); h32_2 = SUM(tast32*wgts*J*N2); h33_2 = SUM(tast33*wgts*J*N2)
    ! [H]j3 components
    h11_3 = SUM(tast11*wgts*J*N3); h12_3 = SUM(tast12*wgts*J*N3); h13_3 = SUM(tast13*wgts*J*N3) 
    h21_3 = SUM(tast21*wgts*J*N3); h22_3 = SUM(tast22*wgts*J*N3); h23_3 = SUM(tast23*wgts*J*N3)
    h31_3 = SUM(tast31*wgts*J*N3); h32_3 = SUM(tast32*wgts*J*N3); h33_3 = SUM(tast33*wgts*J*N3)

    ! Numeric integration of [h]
    h_el(1:9) = (/h11_1,h12_1,h13_1,h11_2,h12_2,h13_2,h11_3,h12_3,h13_3/)
    h_el(10:18) = (/h21_1,h22_1,h23_1,h21_2,h22_2,h23_2,h21_3,h22_3,h23_3/)
    h_el(19:27) = (/h31_1,h32_1,h33_1,h31_2,h32_2,h33_2,h31_3,h32_3,h33_3/)

END SUBROUTINE H_G_Non_singularv
!===============================================================================!
SUBROUTINE Fund_Solution_UT(theta,phi,r,nfp,Uuv,Tij, &
                              RIt_cell, RIc_cell, RI0m_cell, RIm0_cell, R00, pos_vet, &
                              size_coef_cells,C)

    INTEGER :: nRIt, nRIc, nRI0m, nRIm0, i, m, n
    REAL(8) :: theta, phi, r, nfp(:), Uuv(:,:), Tij(:,:), w0_1, w1_1, w2_1, w0_2
    REAL(8) :: w1_2, w2_2, w0_3, w1_3, w2_3, sum1(3,3), sum2(3,3), sum3(3,3)
    REAL(8) :: Rt_mn(3,3), It_mn(3,3), Ga(3,3), Gt(3,3), Rc_mn(3,3), Ic_mn(3,3)
    REAL(8) :: Gb(3,3), Gc(3,3), R0m(3,3), I0m(3,3), ggb(3,3), ggc(3,3)
    REAL(8) :: Rm0(3,3), Im0(3,3), gga(3,3), ggt(3,3), Uuv_1(3,3), Uuv_2(3,3)
    REAL(8) :: Uuv_3(3,3), U11_1, U12_1, U13_1, U21_1, U22_1, U23_1, U31_1
    REAL(8) :: U32_1, U33_1, U11_2, U12_2, U13_2, U21_2, U22_2, U23_2, U31_2 
    REAL(8) :: U32_2, U33_2, U11_3, U12_3, U13_3, U21_3, U22_3, U23_3, U31_3
    REAL(8) :: U32_3, U33_3, S1_vet(6), S2_vet(6), S3_vet(6), S1_mat(3,3)
    REAL(8) :: S2_mat(3,3), S3_mat(3,3), Error_1, Error_2, Tol
    REAL(8) :: phi_aux, theta_aux

    INTEGER :: size_coef_cells(4)
    REAL(KIND=8) :: C(6,6) ! Compliance tensor

    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIt_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIc_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RI0m_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIm0_cell
    INTEGER, ALLOCATABLE :: pos_vet(:)
    REAL(KIND=8):: R00(3,3)
        
    nRIt = size_coef_cells(1); nRIc = size_coef_cells(2)
    nRI0m = size_coef_cells(3); nRIm0 = size_coef_cells(4)
    
    w0_1 = DSIN(phi)*DCOS(theta);
    w0_2 = DSIN(phi)*DSIN(theta);
    w0_3 = DCOS(phi);    

    w1_1 = -DSIN(theta)/DSIN(phi);
    w1_2 = DCOS(theta)/DSIN(phi);
    w1_3 = 0.d0;

    w2_1 = DCOS(phi)*cos(theta);
    w2_2 = DCOS(phi)*DSIN(theta);
    w2_3 = -DSIN(phi);

    ! Treatment of syngularity -------------------------
    phi_aux = phi; theta_aux = theta
    Tol = 1e-6
    Error_1 = ABS(phi_aux) - 0.d0;
    Error_2 = ABS(phi_aux) - (Pi);
    IF (ABS(Error_1).LT.Tol) THEN ! when phi = 0 
        theta_aux = 0.d0; 
        phi_aux = 1e-6
        w1_1 = -DSIN(theta_aux)/DSIN(phi_aux);
        theta_aux = Pi/2;
        w1_2 = DCOS(theta_aux)/DSIN(phi_aux);
    END IF
    IF (ABS(Error_2).LT.Tol) THEN ! When phi = Pi
        theta_aux = 0.d0; phi_aux = 1e-6
        w1_1 = -DSIN(theta_aux)/DSIN(phi_aux);
        theta_aux = Pi/2;
        w1_2 = DCOS(theta_aux)/DSIN(phi_aux);
    END IF
    ! ---------------------------------------------------


    sum1 = 0.d0; sum2 = 0.d0; sum3 = 0.d0
    
    DO i=1,nRIt
        
        Rt_mn = RIt_cell(i)%re; It_mn = RIt_cell(i)%im
        m = RIt_cell(i)%pos_mat(1); n = RIt_cell(i)%pos_mat(2)
        Ga = Rt_mn*DCOS(m*theta) - It_mn*DSIN(m*theta)
        Gt = Rt_mn*DSIN(m*theta) + It_mn*DCOS(m*theta)

        sum1 = sum1 + Ga*DCOS(n*phi)
        sum2 = sum2 + m*Gt*DCOS(n*phi)
        sum3 = sum3 + n*Ga*DSIN(n*phi)
        
    END DO
    
    DO i=1,nRIc

        Rc_mn = RIc_cell(i)%re; Ic_mn = RIc_cell(i)%im
        m = RIc_cell(i)%pos_mat(1); n = RIc_cell(i)%pos_mat(2)

        Gb = Rc_mn*DSIN(m*theta) + Ic_mn*DCOS(m*theta)
        Gc = Rc_mn*DCOS(m*theta) - Ic_mn*DSIN(m*theta)

        sum1 = sum1 - Gb*DSIN(n*phi)
        sum2 = sum2 + m*Gc*DSIN(n*phi)
        sum3 = sum3 + n*Gb*DCOS(n*phi)
        
    END DO
    
    DO i=1,nRI0m

        R0m = RI0m_cell(i)%re; I0m = RI0m_cell(i)%im
        m = pos_vet(i)

        ggb = R0m*DCOS(m*phi) - I0m*DSIN(m*phi)
        ggc = R0m*DSIN(m*phi) + I0m*DCOS(m*phi)

        sum1 = sum1 + ggb
        sum3 = sum3 + m*ggc
                
    END DO
    
    DO i=1,nRIm0

        Rm0 = RIm0_cell(i)%re; Im0 = RIm0_cell(i)%im
        m = pos_vet(i)

        gga = Rm0*DCOS(m*theta) - Im0*DSIN(m*theta)
        ggt = Rm0*DSIN(m*theta) + Im0*DCOS(m*theta)

        sum1 = sum1 + gga
        sum2 = sum2 + m*ggt
                
    END DO
    
    sum1 = sum1 + R00/2

    Uuv_1 = 1.d0/(2.d0*pi*(r**2))*(-w0_1*sum1 - w1_1*sum2 - w2_1*sum3);
    Uuv_2 = 1.d0/(2.d0*pi*(r**2))*(-w0_2*sum1 - w1_2*sum2 - w2_2*sum3);
    Uuv_3 = 1.d0/(2.d0*pi*(r**2))*(-w0_3*sum1 - w1_3*sum2 - w2_3*sum3);

    Uuv = (sum1)/(2.d0*pi*r);
    
    U11_1 = Uuv_1(1,1); U12_1 = Uuv_1(1,2); U13_1 = Uuv_1(1,3)
    U21_1 = Uuv_1(2,1); U22_1 = Uuv_1(2,2); U23_1 = Uuv_1(2,3)
    U31_1 = Uuv_1(3,1); U32_1 = Uuv_1(3,2); U33_1 = Uuv_1(3,3)

    U11_2 = Uuv_2(1,1); U12_2 = Uuv_2(1,2); U13_2 = Uuv_2(1,3)
    U21_2 = Uuv_2(2,1); U22_2 = Uuv_2(2,2); U23_2 = Uuv_2(2,3)
    U31_2 = Uuv_2(3,1); U32_2 = Uuv_2(3,2); U33_2 = Uuv_2(3,3)

    U11_3 = Uuv_3(1,1); U12_3 = Uuv_3(1,2); U13_3 = Uuv_3(1,3)
    U21_3 = Uuv_3(2,1); U22_3 = Uuv_3(2,2); U23_3 = Uuv_3(2,3)
    U31_3 = Uuv_3(3,1); U32_3 = Uuv_3(3,2); U33_3 = Uuv_3(3,3)

    S1_vet = MATMUL(C,(/U11_1,U21_2,U31_3,(U21_3 + U31_2),(U11_3 + U31_1),(U11_2 + U21_1)/))
    S2_vet = MATMUL(C,(/U12_1,U22_2,U32_3,(U22_3 + U32_2),(U12_3 + U32_1),(U12_2 + U22_1)/))
    S3_vet = MATMUL(C,(/U13_1,U23_2,U33_3,(U23_3 + U33_2),(U13_3 + U33_1),(U13_2 + U23_1)/))

    S1_mat(1,:) = (/S1_vet(1),  S1_vet(6),  S1_vet(5)/)
    S1_mat(2,:) = (/S1_vet(6),  S1_vet(2),  S1_vet(4)/)
    S1_mat(3,:) = (/S1_vet(5),  S1_vet(4),  S1_vet(3)/)

    S2_mat(1,:) = (/S2_vet(1),  S2_vet(6),  S2_vet(5)/)
    S2_mat(2,:) = (/S2_vet(6),  S2_vet(2),  S2_vet(4)/)
    S2_mat(3,:) = (/S2_vet(5),  S2_vet(4),  S2_vet(3)/)
      
    S3_mat(1,:) = (/S3_vet(1),  S3_vet(6),  S3_vet(5)/)
    S3_mat(2,:) = (/S3_vet(6),  S3_vet(2),  S3_vet(4)/)
    S3_mat(3,:) = (/S3_vet(5),  S3_vet(4),  S3_vet(3)/)
    
    Tij(1,:) = MATMUL(S1_mat,nfp)
    Tij(2,:) = MATMUL(S2_mat,nfp)
    Tij(3,:) = MATMUL(S3_mat,nfp)

END SUBROUTINE Fund_Solution_UT
!===============================================================================!
SUBROUTINE New_C_H_Matrices(nt,me,reg,ELEMv,NORMAL_VECTORSv)

	INTEGER :: nt, me, reg,i
	CHARACTER(LEN=20)::num_str, Test_name
	TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIt_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIc_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RI0m_cell
    TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIm0_cell
	INTEGER, ALLOCATABLE :: pos_vet(:)
	REAL(KIND=8):: R00(3,3)
	INTEGER :: size_coef_cells(4)
	REAL(KIND=8) :: C(6,6), S(6,6), Anglesv(nreg,4) ! Compliance tensor
	
	REAL(8), ALLOCATABLE :: G_reg(:,:), H_reg(:,:)
	INTEGER, ALLOCATABLE :: ELEMv(:)
    REAL(8), ALLOCATABLE :: NORMAL_VECTORSv(:,:)
	
    WRITE(num_str,"(I10)") reg
    Test_name = "Test_"//TRIM(ADJUSTL(num_str))//".dat" 
    Anglesv(reg,1) = reg
	
	! Read the Fourier coefficients-----------------------------------------------------! 
	CALL Fourier_Coeff(me,Test_name,reg,RIt_cell, RIc_cell,RI0m_cell, &
                       RIm0_cell,R00,pos_vet,size_coef_cells,C,S,Anglesv)
    
    ! Compute matrices [H] and [G]------------------------------------------------------!   
    CALL Compute_H_G(nt,me,ELEMv,NORMAL_VECTORSv,H_reg,G_reg,reg,RIt_cell, &
                        RIc_cell, RI0m_cell,RIm0_cell,R00,pos_vet,size_coef_cells,C)
	
	DEALLOCATE(RIt_cell, RIc_cell, RI0m_cell, RIm0_cell, pos_vet)

	IF (me.EQ.0) THEN 

        Angles(reg,:) = Anglesv(reg,:)  
        
        H_G(reg)%H = 0.d0; H_G(reg)%G = 0.d0;

        H_G(reg)%G = G_reg 
        IF (Transient.EQ.0) THEN
            H_G(reg)%H = H_reg;  
        END IF
	        
	DEALLOCATE(H_reg, G_reg)

    END IF 
	

	

	

END SUBROUTINE  New_C_H_Matrices
!===============================================================================!
END MODULE Comp_H_G
