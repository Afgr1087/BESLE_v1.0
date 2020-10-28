!===============================================================================!
!-------------------------------------------------------------------------------!
!          MODULE TO COMPUTE THE VECTOR {A} AND VECTOR {b} AFTER FAILURE        !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Sparse_Vectors_Failure
!-------------------------------------------------------------------------------!
USE Set_parameters
USE Global_variables
USE Vectors_b
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE Re_General_Assembly

	INTEGER, ALLOCATABLE :: ID_block_reg(:,:),ID_index_reg(:,:),columns(:,:)
	INTEGER :: n_V, n_G, i, nel_reg, ndofs_reg1, nint_reg, ndofs_reg2, n_rows
	INTEGER :: dof, nelem, j, k, n_A
	REAL(8), ALLOCATABLE :: A_values_aux(:), G_values_aux(:), V_values_aux(:)
    INTEGER, ALLOCATABLE :: columns_A_aux(:), rows_A_aux(:)
    INTEGER, ALLOCATABLE :: columns_G_aux(:)
    INTEGER, ALLOCATABLE :: columns_V_aux(:) 

	INTEGER :: Tol, cf, cf2, cf3, fi, fi2, fi3, nzA, nzG, nzV
	INTEGER :: col_blocks, ndf, cont_G, cont_V, reg, element, pos
	INTEGER :: ini_col, ind, ii
	REAL(8) :: val
	INTEGER :: ci, ci2, ci3

    ! To find the index of the sparse matrix
    ! this function is in Comp_H_G
    
    CALL Sparce_array_Failure(ID_block_reg,ID_index_reg,columns)

	dof = 3
	nelem = SIZE(ELEM,1) ! Number of elements

	n_V = 0
    n_G = 0
    DO i=1,nreg

        nel_reg = El_reg(i,2)
        ndofs_reg1 = nel_reg*nnos_el*dof    
        n_V = n_V + ndofs_reg1**2 

        nint_reg = Int_reg(i,2)
        ndofs_reg2 = nint_reg*nnos_el*dof        
        n_G = n_G + (ndofs_reg2*ndofs_reg1)

    END DO
    
    n_A = n_V + n_G
    n_G = n_V - n_G
    
    n_rows = nelem*nnos_el*dof+1

    ALLOCATE(A_values_aux(n_A), columns_A_aux(n_A), rows_A_aux(n_A))
    ALLOCATE(G_values_aux(n_G), columns_G_aux(n_G), rows_G(n_rows))
	
	! Matrix [V] Body forces analysis
    IF (Bodyforces.EQ.1) THEN
        ALLOCATE(V_values_aux(n_V), columns_V_aux(n_A), rows_V(n_rows))
        V_values_aux = 0.d0
        columns_V_aux = 0; 
        rows_V = 0
        CALL Bodyforces_vet
    END IF

    ! Matrix [V] Trasient analysis
    IF (Transient.EQ.1) THEN
        ALLOCATE(V_values_aux(n_A), columns_V_aux(n_A), rows_V(n_rows))
        V_values_aux = 0.d0
        columns_V_aux = 0; 
        rows_V = 0
    END IF

	A_values_aux = 0.d0; columns_A_aux = 0; rows_A_aux = 0
    G_values_aux = 0.d0; columns_G_aux = 0; rows_G = 0

	Tol = 1E-25
    nel_reg = 0;
    cf = 0; cf2 = 0; cf3 = 0;
    fi = 0; 
    fi2 = 0; fi3 = 0;
    nzA = 0; nzG = 0; nzV = 0;

	DO i=1,nreg
        col_blocks = ID_block_reg(i,2)
        nel_reg = El_reg(i,2)
        ndf = 0;
        DO j = 1,nel_reg*nnos_el*dof 
            ndf = ndf + 1 ! number of row in local region
            !cont_A = 0;     
            cont_G = 0; 
            cont_V = 0;   
            fi = fi + 1
            fi2 = fi2 + 1
            fi3 = fi3 + 1
            DO k=1,col_blocks 
                reg = ID_block_reg(i,k+2)
                element = ID_index_reg(i,k+2)
                pos = columns(i,k+1);
                IF (reg.EQ.(nreg+1)) THEN ! Boundary block

					ini_col = ELEM(element,6) !Initial index in matrix H of reg
                    
                    !-----------------------------------------------------------------
                    ci = cf + 1           ! initial column for A       
                    cf = cf + nnos_el*dof ! final column for A 
                    
                    ind = 0;
                    DO ii=ci,cf
                        val = H_G(i)%H(ndf,ini_col+ind)
                        IF (ABS(val) .GT. Tol) THEN
                            nzA = nzA + 1
                            A_values_aux(nzA) = H_G(i)%H(ndf,ini_col+ind); 
                            columns_A_aux(nzA) = pos + ind
                            rows_A_aux(nzA) = fi
                            !IF (cont_A .EQ. 0) THEN
                            !    fi = fi + 1
                            !    rows_A_aux(fi) = nzA
                            !    cont_A = 1
                            !END IF
                        END IF
                        ind = ind + 1
                    END DO
                    !-----------------------------------------------------------------
                    ci2 = cf2 + 1           ! initial column for G 
                    cf2 = cf2 + nnos_el*dof ! final column for G   

                    ind = 0;
                    
                    DO ii=ci2,cf2
                        val = H_G(i)%G(ndf,ini_col+ind)
                        IF (ABS(val) .GT. Tol) THEN 
                            nzG = nzG + 1
                            G_values_aux(nzG) = H_G(i)%G(ndf,ini_col+ind); 
                            columns_G_aux(nzG) = pos + ind                   
                            IF (cont_G .EQ. 0) THEN  
                                rows_G(fi2) = nzG
                                cont_G = 1
                            END IF
                        END IF
                        ind = ind + 1
                    END DO
                    !-----------------------------------------------------------------
                    ci3 = cf3 + 1           ! initial column for V  
                    cf3 = cf3 + nnos_el*dof ! final column for V

                    IF ((Bodyforces.EQ.1).OR.(Transient.EQ.1)) THEN
                        ind = 0;
                        DO ii=ci3,cf3  
                            val = U_T_E_M(i)%M(ndf,ini_col+ind)
                            IF (ABS(val) .GT. Tol) THEN
                                nzV = nzV + 1
                                V_values_aux(nzV) = U_T_E_M(i)%M(ndf,ini_col+ind); 
                                columns_V_aux(nzV) = pos + ind
                                IF (cont_V .EQ. 0) THEN 
                                    rows_V(fi3) = nzV
                                    cont_V = 1;
                                END IF
                            END IF
                            ind = ind + 1
                        END DO
                    END IF
                    !-----------------------------------------------------------------
                END IF
				IF ((i.LT.ABS(reg)).AND.(reg.NE.(nreg+1))) THEN ! reg_A < reg_B

                    IF (reg.GT.0) THEN ! [F]^a 

                        ini_col = ELEM(element,6) !Initial index in matrix H of reg
                       
                        !-----------------------------------------------------------------
                        ci = cf + 1           ! initial column for A       
                        cf = cf + nnos_el*dof ! final column for A 

                        ind = 0;
                        DO ii=ci,cf
                            val = H_G(i)%H(ndf,ini_col+ind)
                            IF (ABS(val) .GT. Tol) THEN
                                nzA = nzA + 1
                                A_values_aux(nzA) = H_G(i)%H(ndf,ini_col+ind); 
                                columns_A_aux(nzA) = pos + ind
                                rows_A_aux(nzA) = fi
                                !IF (cont_A .EQ. 0) THEN
                                !    fi = fi + 1
                                !    rows_A_aux(fi) = nzA
                                !    cont_A = 1
                                !END IF
                            END IF
                            ind = ind + 1
                        END DO
                        !-----------------------------------------------------------------
                        ci3 = cf3 + 1           ! initial column for V  
                        cf3 = cf3 + nnos_el*dof ! final column for V

                        IF ((Bodyforces.EQ.1).OR.(Transient.EQ.1)) THEN
                            ind = 0;
                            DO ii=ci3,cf3 
                                val = U_T_E_M(i)%M(ndf,ini_col+ind)
                                IF (ABS(val) .GT. Tol) THEN
                                    nzV = nzV + 1
                                    V_values_aux(nzV) = U_T_E_M(i)%M(ndf,ini_col+ind); 
                                    columns_V_aux(nzV) = pos + ind
                                    IF (cont_V .EQ. 0) THEN 
                                        rows_V(fi3) = nzV
                                        cont_V = 1;
                                    END IF
                                END IF
                                ind = ind + 1
                            END DO
                        END IF
                        !-----------------------------------------------------------------
                        
                    ELSE               ! [G]^a 

                        ini_col = ELEM(element,6) !Initial index in matrix H of reg
                                                
                        !-----------------------------------------------------------------
                        ci = cf + 1           ! initial column for A       
                        cf = cf + nnos_el*dof ! final column for A 

                        ind = 0;
                        DO ii=ci,cf
                            val = H_G(i)%H(ndf,ini_col+ind)
                            IF (ABS(val) .GT. Tol) THEN
                                nzA = nzA + 1
                                A_values_aux(nzA) = H_G(i)%G(ndf,ini_col+ind); 
                                columns_A_aux(nzA) = pos + ind
                                rows_A_aux(nzA) = fi
                                !IF (cont_A .EQ. 0) THEN
                                !    fi = fi + 1
                                !    rows_A_aux(fi) = nzA
                                !    cont_A = 1
                                !END IF
                            END IF
                            ind = ind + 1
                        END DO
                        !-----------------------------------------------------------------
                        
                    END IF

                END IF
                IF ((i.GT.ABS(reg)).AND.(reg.NE.(nreg+1))) THEN ! reg_A > reg_B

                    IF (reg.GT.0) THEN  ! [G]^b 

                        ini_col = ELEM(element,6) !Initial index in matrix H of reg
                                                              
                        !-----------------------------------------------------------------
                        ci = cf + 1           ! initial column for A       
                        cf = cf + nnos_el*dof ! final column for A 

                        ind = 0;
                        DO ii=ci,cf
                            val = H_G(i)%H(ndf,ini_col+ind)
                            IF (ABS(val) .GT. Tol) THEN
                                nzA = nzA + 1
                                A_values_aux(nzA) = -H_G(i)%G(ndf,ini_col+ind); 
                                columns_A_aux(nzA) = pos + ind
                                rows_A_aux(nzA) = fi
                                !IF (cont_A .EQ. 0) THEN
                                !    fi = fi + 1
                                !    rows_A_aux(fi) = nzA
                                !    cont_A = 1
                                !END IF
                            END IF
                            ind = ind + 1
                        END DO
                        !-----------------------------------------------------------------
                        
                    ELSE                ! [F]^b 

                        ini_col = ELEM(element,6) !Initial index in matrix H of reg
                       
                        !-----------------------------------------------------------------
                        ci = cf + 1           ! initial column for A       
                        cf = cf + nnos_el*dof ! final column for A 

                        ind = 0;
                        DO ii=ci,cf
                            val = H_G(i)%H(ndf,ini_col+ind)
                            IF (ABS(val) .GT. Tol) THEN
                                nzA = nzA + 1
                                A_values_aux(nzA) = H_G(i)%H(ndf,ini_col+ind); 
                                columns_A_aux(nzA) = pos + ind
                                rows_A_aux(nzA) = fi
                                !IF (cont_A .EQ. 0) THEN
                                !    fi = fi + 1
                                !    rows_A_aux(fi) = nzA
                                !    cont_A = 1
                                !END IF
                            END IF
                            ind = ind + 1
                        END DO
                        !-----------------------------------------------------------------
                        ci3 = cf3 + 1           ! initial column for V  
                        cf3 = cf3 + nnos_el*dof ! final column for V

                        IF ((Bodyforces.EQ.1).OR.(Transient.EQ.1)) THEN
                            ind = 0;
                            DO ii=ci3,cf3  
                                val = U_T_E_M(i)%M(ndf,ini_col+ind)
                                IF (ABS(val) .GT. Tol) THEN
                                    nzV = nzV + 1
                                    V_values_aux(nzV) = U_T_E_M(i)%M(ndf,ini_col+ind); 
                                    columns_V_aux(nzV) = pos + ind
                                    IF (cont_V .EQ. 0) THEN 
                                        rows_V(fi3) = nzV
                                        cont_V = 1;
                                    END IF
                                END IF
                                ind = ind + 1
                            END DO
                        END IF
                        !-----------------------------------------------------------------
                      
                    END IF
                END IF
                
            END DO
        END DO
    END DO

	DEALLOCATE(ID_block_reg,ID_index_reg,columns)

	! Matrix A --------------------------------------------------

    ALLOCATE(A_values(nzA))      
    A_values = A_values_aux(1:nzA)   
    DEALLOCATE(A_values_aux)
    
	ALLOCATE(columns_A(nzA))
    columns_A = columns_A_aux(1:nzA)    
    DEALLOCATE(columns_A_aux)
    
	ALLOCATE(rows_A(nzA))  
    rows_A = rows_A_aux(1:nzA)
    DEALLOCATE(rows_A_aux)

    ! Matrix G --------------------------------------------------

    rows_G(fi2+1) = nzG + 1
	ALLOCATE(G_values(nzG))
    G_values = G_values_aux(1:nzG)
    DEALLOCATE(G_values_aux)
    
	ALLOCATE(columns_G(nzG))
    columns_G = columns_G_aux(1:nzG)
    DEALLOCATE(columns_G_aux)

    ! Matrix V --------------------------------------------------
    IF ((Bodyforces.EQ.1).OR.(Transient.EQ.1)) THEN
        rows_V(fi3+1) = nzV + 1
        ALLOCATE(V_values(nzV),columns_V(nzV))
        V_values = V_values_aux(1:nzV)
        columns_V = columns_V_aux(1:nzV)
        DEALLOCATE(V_values_aux,columns_V_aux)
    END IF 
	
	! -------------------------------------------------------------------------

END SUBROUTINE Re_General_Assembly
!===============================================================================!
SUBROUTINE Sparce_array_Failure(ID_block_reg,ID_index_reg,columns)

	INTEGER :: dof, val1, val2, val3, i, nelem, size_nelem
	INTEGER, ALLOCATABLE :: vec1(:), vec2(:), vec3(:), ID_block_reg(:,:)
	INTEGER, ALLOCATABLE :: ID_index_reg(:,:),columns(:,:),Vec_indices(:) 
	INTEGER :: j, element, el, k, element_A, element_B, reg_A, reg_B
	INTEGER :: ndci, ndcf, ind, ind_B, ini_col, fin_col

	dof = 3 ! degrees of fredom

	! displs of regions
    val1 = 0
    El_reg(1,5) = 0
    val2 = 0
    El_reg(1,6) = 0
    val3 = 0
    Int_reg(1,3) = 0
    DO i=1,nreg-1
        val1 = val1 + El_reg(i,2)*nnos_el*dof
        El_reg(i+1,5) = val1
        val2 = val2 + El_reg(i,2)
        El_reg(i+1,6) = val2
        val3 = val3 + Int_reg(i,2)
        Int_reg(i+1,3) = val3
    END DO
	
	nelem = SIZE(ELEM,1) ! Number of elements

	
	
    ALLOCATE(vec1(nreg), vec2(nreg), vec3(nreg))
    vec1 = El_reg(:,2)
    vec2 = Int_reg(:,2)
    vec3 = vec1 + vec2
    size_nelem = MAXVAL(vec3)
    DEALLOCATE(vec1, vec2)
     
    ALLOCATE(ID_block_reg(nreg,size_nelem+2))
    ALLOCATE(ID_index_reg(nreg,size_nelem+2))
    ALLOCATE(columns(nreg,1+size_nelem))
    ALLOCATE(Vec_indices(nreg))  

	Vec_indices = 0
    columns = 0
    ID_block_reg = 0; ID_index_reg = 0; 

	DO i=1,nreg ! Loop over regions
        ID_block_reg(i,1:2) = (/i,vec3(i)/) ! [region, N_blocks]
        ID_index_reg(i,1:2) = (/i,vec3(i)/) ! [region, N_blocks]
             
        DO j=1,El_reg(i,2) ! loop over elements per regions    
            element = Subregions(i,j) ! element
            el = ELEM(element,12) ! see if the element is an interface
            IF (el.NE.0) THEN
                !-------------------------------------------------------------------------
                ! BLOCKS THAT BELONG TO THE INTERFACES, APPICATION OF DISPLACEMENT
                ! COMPATIBILITY AND TRACTION EQUILIBRIUM CONDITIONS
                !-------------------------------------------------------------------------
                
                k = ELEM(element,13) ! This is the interface
                
                element_A = Interfaces(k,4) ! Element in region A
                element_B = Interfaces(k,5) ! ELement in region B

                reg_A = Interfaces(k,2);
                reg_B = Interfaces(k,3);
                !-------------------------------------------------------------------------
                ! DISPLACEMENT COMPATIBILTY IN THE INTERFACES regA<regB OR regA>regB 
                ! CALCULATION OF BLOCK OF THE REGIONS A AND B
                !-------------------------------------------------------------------------
                ! index columns                      
                !initial index of the block in the general matrix 
                !final index of the block in the general matrix 
                ndci = El_reg(reg_A,5) + (j - 1)*nnos_el*dof + 1 
                ndcf = ndci + nnos_el*dof - 1
                            
                ! Index to use in the order results module
                ELEM(element_A,8:9) = (/ndci,ndcf/) 
                    
                ! DATA FOR SPARCE SYSTEM ==========================================
                
                ind = Vec_indices(reg_A)
                ind = ind  + 1
                Vec_indices(reg_A) = ind
                ind_B = Vec_indices(reg_B)
                ind_B = ind_B + 1
                Vec_indices(reg_B) = ind_B
 
                columns(reg_A,1+ind) = ndci
                columns(reg_B,1+ind_B) = ndci
 
                ID_block_reg(reg_A,2+ind) = reg_B
                ID_block_reg(reg_B,2+ind_B) = -reg_A
                                                        
                ID_index_reg(reg_A,2+ind) = element_A
                ID_index_reg(reg_B,2+ind_B) = element_B
				
            ELSE
                
                !-------------------------------------------------------------------------
                ! BLOCKS THAT BELONG TO BOUNDARY SEGMENTS OF THE GEOMETRY WITH KNOWN
                ! BOUNDARY CONDITIONS
                !-------------------------------------------------------------------------
                
                ! index columns   
                !initial index of the block in the general matrix 
                !final index of the block in the general matrix 
                ndci = El_reg(i,5) + (j - 1)*nnos_el*dof + 1 
                ndcf = ndci + nnos_el*dof - 1

                ! Index to use in the order results module
                ELEM(element,8:9) = (/ndci,ndcf/) 

                ini_col = ELEM(element,6) !Initial index in matrix H of reg
                fin_col = ELEM(element,7) !Final index in matrix H of reg

                ! DATA FOR SPARCE SYSTEM ========================================== 
                
                ind = Vec_indices(i)
                ind = ind  + 1
                Vec_indices(i) = ind
                        
                columns(i,1+ind) = ndci

                ID_block_reg(i,2+ind) = nreg+1

                ID_index_reg(i,2+ind) = element   
                ! =================================================================

            END IF
            
        END DO
        
    END DO

	DEALLOCATE(Vec_indices,vec3) 
	
END SUBROUTINE Sparce_array_Failure
!===============================================================================!
END MODULE Sparse_Vectors_Failure
