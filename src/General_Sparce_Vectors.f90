!===============================================================================!
!-------------------------------------------------------------------------------!
!          MODULE TO COMPUTE THE VECTOR {A} AND VECTOR {b}                      !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE General_Sparce_Vectors
!-------------------------------------------------------------------------------!
USE Global_variables
USE Global_functions
USE Vectors_b
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE Final_Sparce_Vectors(ts1,step,nt)

    INTEGER::t0,t1,rate, t2,t3, rate2, step, nt
    REAL(8),INTENT(OUT)::ts1
    
    INTEGER, ALLOCATABLE :: ID_block_reg(:,:),ID_index_reg(:,:),columns(:,:)

    INTEGER :: ndf, cf, i, col_blocks, nel_reg, j, dof, k, reg
    INTEGER :: element, ci, ini_col, nelem
    INTEGER :: cf2, ci2, ci3, cf3, nint_reg, n_rows
    INTEGER :: fi2, fi3, ind, ii, pos, fi !, cont_A
    INTEGER :: cont_V, cont_G, el_I

    REAL(8) :: val, Tol, Tol_V, Tol_G
	INTEGER :: cohe_val, condition_int, Int_el, ninterfaces, reg_A, cond_node 

	INTEGER :: proc_G, proc_V, proc_A
	INTEGER(8), ALLOCATABLE :: s_A_aux(:), s_G_aux(:), s_V_aux(:)  
	INTEGER(8), ALLOCATABLE :: displs_A(:) 
	INTEGER(8) :: nzA_loc, nzG_loc, nzV_loc, n_A, n_G, n_V
	INTEGER(8) :: ndofs_reg1, ndofs_reg2, val1, nt_sub

	TYPE Cell12	
	INTEGER, ALLOCATABLE :: col_A_aux(:), col_G_aux(:), col_V_aux(:)
	INTEGER, ALLOCATABLE :: row_A_aux(:)
	REAL(8), ALLOCATABLE :: A_aux(:), G_aux(:), V_aux(:)
    END TYPE Cell12
    TYPE(Cell12), DIMENSION(:), ALLOCATABLE :: S_E_aux


    CALL SYSTEM_CLOCK(t0, rate) ! Initializates the time counter
    WRITE (*,*) '----------------------------------------------------------------------------'
    WRITE (*,*) ' '
    WRITE (*,'(A)',advance='no') ' 6. Vectors {A} and {b}           ...'
    WRITE (*,*) ' '

    ! To find the index of the sparse matrix
    ! this function is in Comp_H_G
    
    CALL Sparce_array(step,ID_block_reg,ID_index_reg,columns)
    
    CALL SYSTEM_CLOCK(t2, rate2) ! Initializates the time counter
    WRITE (*,'(A)',advance='no') '       - 6.2 Vectors              ...'

    dof = 3;
    nelem = SIZE(ELEM,1) ! Number of elements
	ninterfaces = SIZE(Interfaces,1) ! Number of elements
	
	val1 = 0
	IF (step.EQ.1) THEN

    	ALLOCATE(b(nelem*nnos_el*dof,1),b_cte(nelem*nnos_el*dof,1))
    	ALLOCATE(b_step(nelem*nnos_el*dof,1), b_bf(nelem*nnos_el*dof,1)) 
		

		b = 0.d0; b_step = 0.d0; b_cte = 0.d0; b_bf = 0.d0

		Reduce_A = 0
		
	ELSEIF ((step.GT.1).AND.(Iteration.EQ.1)) THEN 
		
		DEALLOCATE(S_E)
		DEALLOCATE(rows_G)
		DEALLOCATE(scounts_A, scounts_G, displs_G)
		DEALLOCATE(b,b_step,b_cte)
		
		IF ((Bodyforces.EQ.1).OR.(Transient.EQ.1)) THEN
			DEALLOCATE(rows_V,scounts_V, displs_V)
		END IF 
		
		DO i=1,nreg
			val1 = val1 + El_reg(i,8)
		END DO
		
		ALLOCATE(b(nelem*nnos_el*dof+val1,1),b_cte(nelem*nnos_el*dof+val1,1))
		ALLOCATE(b_step(nelem*nnos_el*dof+val1,1))
  			
		b = 0.d0; b_step = 0.d0; b_cte = 0.d0

	END IF	

	! Distribution of vectors A, G, V throughout the processors  -------------------------

	nt_sub = nt - 1 ! Threads without the master "0"

	ALLOCATE(S_E_aux(nt_sub)) ! system of equations "aux"

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

	ALLOCATE(s_A_aux(nt-1), s_G_aux(nt-1))
	s_A_aux = 0; s_G_aux = 0

	ALLOCATE(displs_A(nt-1), displs_G(nt-1))
	displs_A = 0; displs_G = 0
	
	
	N_total = nelem*nnos_el*dof ! rows or columns  of general system

	! Function for parallel distribution
	! Par_Div(threads, vector_size, scounts)
	CALL Par_Div(nt_sub, n_A, s_A_aux, displs_A) ! A_values

	CALL Par_Div(nt_sub, n_G, s_G_aux, displs_G) ! G_values
	
	! Matrix [V] Body forces or transient analyses
    IF ((Bodyforces.EQ.1).OR.(Transient.EQ.1)) THEN
		ALLOCATE(s_V_aux(nt-1),displs_V(nt-1))
		s_V_aux = 0; displs_V = 0
		CALL Par_Div(nt_sub, n_V, s_V_aux, displs_V) ! V_values
		n_V = n_V*(nnos_el*dof)**2 ! Total size of matrix V
		s_V_aux = s_V_aux*(nnos_el*dof)**2 ! Size of each V_loc  
		displs_V = displs_V*(nnos_el*dof)**2 ! Positions of each V_loc   

		! Body forces analysis
		IF (Bodyforces.EQ.1) THEN
		    
		    ALLOCATE(bhat(nelem*nnos_el*dof,1))
		    CALL Bodyforces_vet

		END IF

		! Trasient analysis
	 	IF (Transient.EQ.1) THEN
			IF (step.EQ.1) THEN
		    	ALLOCATE(u_t(nelem*nnos_el*dof,3))
		    	u_t = 0.d0
			END IF		    
		END IF
   
    END IF
	! ----------------------------------------------------------------------------------
	! Allocating A, G and V auxiliar matrices

	n_rows = N_total+1+val1
	
	DO i=1,nt_sub

		ALLOCATE(S_E_aux(i)%A_aux(s_A_aux(i))) ! A_aux
		ALLOCATE(S_E_aux(i)%col_A_aux(s_A_aux(i))) ! columns_Aaux
		ALLOCATE(S_E_aux(i)%row_A_aux(s_A_aux(i))) ! rows_Aaux
		S_E_aux(i)%A_aux = 0.d0
		S_E_aux(i)%col_A_aux = 0
		S_E_aux(i)%row_A_aux = 0

		ALLOCATE(S_E_aux(i)%G_aux(s_G_aux(i))) ! G_aux
		ALLOCATE(S_E_aux(i)%col_G_aux(s_G_aux(i))) ! columns_Gaux
		
		S_E_aux(i)%G_aux = 0.d0
		S_E_aux(i)%col_G_aux = 0
		
		
		! Matrix [V] Body forces or transient analyses
    	IF ((Bodyforces.EQ.1).OR.(Transient.EQ.1)) THEN
			
			ALLOCATE(S_E_aux(i)%V_aux(s_A_aux(i))) ! V_aux			
			ALLOCATE(S_E_aux(i)%col_V_aux(s_A_aux(i))) ! columns_Vaux		
			
			S_E_aux(i)%V_aux = 0.d0
			S_E_aux(i)%col_V_aux = 0
			
		END IF

	END DO

	! Matrix [V] Body forces or transient analyses
    IF ((Bodyforces.EQ.1).OR.(Transient.EQ.1)) THEN
		ALLOCATE(rows_V(n_rows))
		rows_V = 0
	END IF
	
	ALLOCATE(rows_G(n_rows))
	rows_G = 0
	
	! ----------------------------------------------------------------------------------

	Tol = 1E-25
	Tol_V = 1E-25
	Tol_G = 1E-25
    nel_reg = 0;
    cf = 0; cf2 = 0; cf3 = 0;
    fi = 0; 
    fi2 = 0; fi3 = 0;
    nzA = 0; nzG = 0; nzV = 0;
    

	nzA_loc = 0
	proc_A = 1
	nzG_loc = 0
	proc_G = 1
	nzV_loc = 0
    proc_V = 1
	
    DO i=1,nreg
        col_blocks = ID_block_reg(i,2)
        nel_reg = El_reg(i,1)
		
        ndf = 0;
		
        DO j = 1,nel_reg*nnos_el*dof+El_reg(i,8)
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
		                cf = cf + ELEM(element,14)  ! final column for A 
		                
		                ind = 0;
		                DO ii=ci,cf

							IF (ndf.LE.SIZE(H_G(i)%H,1)) THEN
		                    	val = H_G(i)%H(ndf,ini_col+ind)
							ELSE
								val = 0.d0
							END IF
							
		                    IF (ABS(val) .GT. Tol) THEN

		                        nzA = nzA + 1
		                        							
								nzA_loc = nzA_loc + 1
								
								S_E_aux(proc_A)%A_aux(nzA_loc) = H_G(i)%H(ndf,ini_col+ind)
								S_E_aux(proc_A)%col_A_aux(nzA_loc) = pos + ind
								S_E_aux(proc_A)%row_A_aux(nzA_loc) = fi
								
								IF (nzA.EQ.(displs_A(proc_A)+s_A_aux(proc_A))) THEN
									proc_A = proc_A + 1		
									nzA_loc = 0						
								END IF

		                    END IF
		                    ind = ind + 1
		                END DO
		                !-----------------------------------------------------------------
		                ci2 = cf2 + 1           ! initial column for G 
		                cf2 = cf2 + ELEM(element,14)  ! final column for G   

		                ind = 0;
		                
		                DO ii=ci2,cf2

							IF (ndf.LE.SIZE(H_G(i)%G,1)) THEN
		                    	val = H_G(i)%G(ndf,ini_col+ind)
							ELSE
								val = 0.d0
							END IF

		                    IF (ABS(val) .GT. Tol_G) THEN 
		                        nzG = nzG + 1
		                        								                 
		                        IF (cont_G .EQ. 0) THEN  
		                            rows_G(fi2) = nzG
		                            cont_G = 1
		                        END IF

								nzG_loc = nzG_loc + 1
								
								S_E_aux(proc_G)%G_aux(nzG_loc) = H_G(i)%G(ndf,ini_col+ind);
								S_E_aux(proc_G)%col_G_aux(nzG_loc) = pos + ind
								
								IF (nzG.EQ.(displs_G(proc_G)+s_G_aux(proc_G))) THEN
									proc_G = proc_G + 1		
									nzG_loc = 0						
								END IF


		                    END IF
		                    ind = ind + 1
		                END DO
		                !-----------------------------------------------------------------
		                ci3 = cf3 + 1           ! initial column for V  
		                cf3 = cf3 + ELEM(element,14)  ! final column for V

		                IF ((Bodyforces.EQ.1).OR.(Transient.EQ.1)) THEN
		                    ind = 0;
		                    DO ii=ci3,cf3  
				
								IF (ndf.LE.SIZE(U_T_E_M(i)%M,1)) THEN
		                        	val = U_T_E_M(i)%M(ndf,ini_col+ind)
								ELSE
									val = 0.d0
								END IF

		                        IF (ABS(val) .GT. Tol_V ) THEN
		                            nzV = nzV + 1
		                            
		                            IF (cont_V .EQ. 0) THEN 
		                                rows_V(fi3) = nzV
		                                cont_V = 1;
		                            END IF

									nzV_loc = nzV_loc + 1
									
									S_E_aux(proc_V)%V_aux(nzV_loc) = U_T_E_M(i)%M(ndf,ini_col+ind); 
									S_E_aux(proc_V)%col_V_aux(nzV_loc) = pos + ind
								
									IF (nzV.EQ.(displs_V(proc_V)+s_V_aux(proc_V))) THEN
										proc_V = proc_V + 1		
										nzV_loc = 0						
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
                        cf = cf + ELEM(element,14)  ! final column for A 
						
                        ind = 0;

						cohe_val = 0
						Int_el = ELEM(element,13)
						condition_int = Interfaces(Int_el,18) 
                        DO ii=ci,cf

							IF (condition_int.EQ.1) THEN ! Cohesive zone element
								
								cohe_val = cohe_val + 1
								val = Cohesive(element)%H_cohe(ndf,cohe_val)
								
							ELSE							! Pristine element

								IF (ndf.LE.SIZE(H_G(i)%H,1)) THEN
									val = H_G(i)%H(ndf,ini_col+ind)
								ELSE
									val=0.d0
								END IF

							END IF

                            IF (ABS(val) .GT. Tol) THEN
                                nzA = nzA + 1
                                								
                                nzA_loc = nzA_loc + 1
								
								S_E_aux(proc_A)%A_aux(nzA_loc) = val
								S_E_aux(proc_A)%col_A_aux(nzA_loc) = pos + ind
								S_E_aux(proc_A)%row_A_aux(nzA_loc) = fi

								IF (nzA.EQ.(displs_A(proc_A)+s_A_aux(proc_A))) THEN
									proc_A = proc_A + 1		
									nzA_loc = 0						
								END IF

                            END IF
                            ind = ind + 1
                        END DO
			
                        !-----------------------------------------------------------------
                        ci3 = cf3 + 1           ! initial column for V  
                        cf3 = cf3 + ELEM(element,14)  ! final column for V

                        IF ((Bodyforces.EQ.1).OR.(Transient.EQ.1)) THEN
                            ind = 0;
							cohe_val = 0
                            DO ii=ci3,cf3 

								IF (condition_int.EQ.1) THEN ! Cohesive zone element
								
									cohe_val = cohe_val + 1
									val = Cohesive(element)%M_cohe(ndf,cohe_val)
								
								ELSE

									IF (ndf.LE.SIZE(U_T_E_M(i)%M,1)) THEN
				                    	val = U_T_E_M(i)%M(ndf,ini_col+ind)
									ELSE
										val = 0.d0
									END IF

								END IF

                                IF (ABS(val) .GT. Tol_V ) THEN
                                    nzV = nzV + 1
                                    
                                    IF (cont_V .EQ. 0) THEN 
                                        rows_V(fi3) = nzV
                                        cont_V = 1;
                                    END IF

									nzV_loc = nzV_loc + 1
									
									S_E_aux(proc_V)%V_aux(nzV_loc) = val;
									S_E_aux(proc_V)%col_V_aux(nzV_loc) = pos + ind
								
									IF (nzV.EQ.(displs_V(proc_V)+s_V_aux(proc_V))) THEN
										proc_V = proc_V + 1		
										nzV_loc = 0						
									END IF

                                END IF
                                ind = ind + 1
                            END DO
                        END IF
                        !-----------------------------------------------------------------
                        
                    ELSE  !              ! [G]^a 

                        ini_col = ELEM(element,6) !Initial index in matrix H of reg
                                                
                        !-----------------------------------------------------------------

						Int_el = ELEM(element,13)

						el_I = Interfaces(Int_el,5) 

                        ci = cf + 1           ! initial column for A       
                        cf = cf + ELEM(el_I,14)  ! final column for A 
						
						ind = 0;
						cohe_val = 0
						Int_el = ELEM(element,13)
						condition_int = Interfaces(Int_el,18) 
                        DO ii=ci,cf

							IF (condition_int.EQ.1) THEN ! Cohesive zone element
								
								cohe_val = cohe_val + 1
								val = Cohesive(element)%G_cohe(ndf,cohe_val)
								
							ELSE							! Pristine element
						
								IF (ndf.LE.SIZE(H_G(i)%G,1)) THEN
		                        	val = H_G(i)%G(ndf,ini_col+ind)
								ELSE
									val = 0.d0
								END IF
	
							END IF

                            IF (ABS(val) .GT. Tol) THEN
                                nzA = nzA + 1
                                								
                                nzA_loc = nzA_loc + 1
								
								S_E_aux(proc_A)%A_aux(nzA_loc) = val
								S_E_aux(proc_A)%col_A_aux(nzA_loc) = pos + ind
								S_E_aux(proc_A)%row_A_aux(nzA_loc) = fi

								IF (nzA.EQ.(displs_A(proc_A)+s_A_aux(proc_A))) THEN
									proc_A = proc_A + 1		
									nzA_loc = 0
								END IF

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
                        cf = cf + ELEM(element,14) ! final column for A 
						
						ind = 0;
						cohe_val = 0
						Int_el = ELEM(element,13)
						condition_int = Interfaces(Int_el,18) 
                        DO ii=ci,cf

							IF (condition_int.EQ.1) THEN ! Cohesive zone element
								
								cohe_val = cohe_val + 1
								val = Cohesive(element)%G_cohe(ndf,cohe_val)
								
							ELSE							! Pristine element
		
								IF (ndf.LE.SIZE(H_G(i)%G,1)) THEN
		                        	val = -H_G(i)%G(ndf,ini_col+ind)
								ELSE
									val = 0.d0
								END IF
	
							END IF
							
                            IF (ABS(val) .GT. Tol) THEN
                                nzA = nzA + 1
                                								
                                nzA_loc = nzA_loc + 1
								
								S_E_aux(proc_A)%A_aux(nzA_loc) = val
								S_E_aux(proc_A)%col_A_aux(nzA_loc) = pos + ind
								S_E_aux(proc_A)%row_A_aux(nzA_loc) = fi

								IF (nzA.EQ.(displs_A(proc_A)+s_A_aux(proc_A))) THEN
									proc_A = proc_A + 1		
									nzA_loc = 0							
								END IF

                            END IF
                            ind = ind + 1
                        END DO
                        !-----------------------------------------------------------------	


						!-----------------------------------------------------------------
						IF (Failure.EQ.1) THEN
							! Only for failure analysis
				            ci3 = cf3 + 1           ! initial column for V  
				            cf3 = cf3 + ELEM(element,14)  ! final column for V

				            IF ((Bodyforces.EQ.1).OR.(Transient.EQ.1)) THEN
				            	ind = 0;
								cohe_val = 0
				                DO ii=ci3,cf3 

									IF (condition_int.EQ.1) THEN ! Cohesive zone element
								
										cohe_val = cohe_val + 1
										val = Cohesive(element)%M_cohe_G(ndf,cohe_val)
																			
									ELSE

										IF (ndf.LE.SIZE(U_T_E_M(i)%M,1)) THEN
								            val = 0.d0
										ELSE
											val = 0.d0
										END IF

									END IF

				                    IF (ABS(val) .GT. Tol_V ) THEN
				                    	nzV = nzV + 1
				                        
				                        IF (cont_V .EQ. 0) THEN 
				                        	rows_V(fi3) = nzV
				                            cont_V = 1;
				                        END IF
										
										nzV_loc = nzV_loc + 1
									
										S_E_aux(proc_V)%V_aux(nzV_loc) = val; 
										S_E_aux(proc_V)%col_V_aux(nzV_loc) = pos + ind
									
										IF (nzV.EQ.(displs_V(proc_V)+s_V_aux(proc_V))) THEN
											proc_V = proc_V + 1		
											nzV_loc = 0						
										END IF
										
	
				                    END IF
				                    ind = ind + 1
				                END DO
				             END IF
						END IF
                        !-----------------------------------------------------------------
                        
                    ELSE                ! [F]^b 

                        ini_col = ELEM(element,6) !Initial index in matrix H of reg
                       
						Int_el = ELEM(element,13)

						el_I = Interfaces(Int_el,5) 

                        !-----------------------------------------------------------------
                        ci = cf + 1           ! initial column for A       
                        cf = cf + ELEM(el_I,14)  ! final column for A 
						
						ind = 0;
						cohe_val = 0
						
						condition_int = Interfaces(Int_el,18) 
                        DO ii=ci,cf

							IF (condition_int.EQ.1) THEN ! Cohesive zone element
								
								cohe_val = cohe_val + 1
								val = Cohesive(element)%H_cohe(ndf,cohe_val)
								
							ELSE							! Pristine element
								
								IF (ndf.LE.SIZE(H_G(i)%H,1)) THEN
		                        	val = H_G(i)%H(ndf,ini_col+ind)
								ELSE
									val = 0.d0
								END IF

							END IF

                            IF (ABS(val) .GT. Tol) THEN
                                nzA = nzA + 1
                               											
                                nzA_loc = nzA_loc + 1
								
								S_E_aux(proc_A)%A_aux(nzA_loc) = val
								S_E_aux(proc_A)%col_A_aux(nzA_loc) = pos + ind
								S_E_aux(proc_A)%row_A_aux(nzA_loc) = fi

								IF (nzA.EQ.(displs_A(proc_A)+s_A_aux(proc_A))) THEN
									proc_A = proc_A + 1		
									nzA_loc = 0							
								END IF

                            END IF
                            ind = ind + 1
                        END DO
	
                        !-----------------------------------------------------------------
                        ci3 = cf3 + 1           ! initial column for V  
                        cf3 = cf3 + ELEM(element,14)  ! final column for V

                        IF ((Bodyforces.EQ.1).OR.(Transient.EQ.1)) THEN
                            ind = 0;
							cohe_val = 0
                            DO ii=ci3,cf3  

								IF (condition_int.EQ.1) THEN ! Cohesive zone element
								
									cohe_val = cohe_val + 1
									val = Cohesive(element)%M_cohe(ndf,cohe_val)
									
								ELSE

									IF (ndf.LE.SIZE(U_T_E_M(i)%M,1)) THEN
				                    	val = U_T_E_M(i)%M(ndf,ini_col+ind)
									ELSE
										val = 0.d0
									END IF

								END IF

                                IF (ABS(val) .GT. Tol_V ) THEN
                                    nzV = nzV + 1
                                    
                                    IF (cont_V .EQ. 0) THEN 
                                        rows_V(fi3) = nzV
                                        cont_V = 1;
                                    END IF

									nzV_loc = nzV_loc + 1
									
									S_E_aux(proc_V)%V_aux(nzV_loc) = val; 
									S_E_aux(proc_V)%col_V_aux(nzV_loc) = pos + ind
									IF (nzV.EQ.(displs_V(proc_V)+s_V_aux(proc_V))) THEN
										proc_V = proc_V + 1		
										nzV_loc = 0						
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
	
	ALLOCATE(S_E(nt_sub)) 
    
    ! Matrix A --------------------------------------------------
    
	DO i=1,nt_sub-1

		ALLOCATE(S_E(i)%At(s_A_aux(i)))
		S_E(i)%At = S_E_aux(i)%A_aux
		DEALLOCATE(S_E_aux(i)%A_aux)

		ALLOCATE(S_E(i)%col_A(s_A_aux(i)))
		S_E(i)%col_A = S_E_aux(i)%col_A_aux
		DEALLOCATE(S_E_aux(i)%col_A_aux)

		ALLOCATE(S_E(i)%row_A(s_A_aux(i)))
		S_E(i)%row_A = S_E_aux(i)%row_A_aux
		DEALLOCATE(S_E_aux(i)%row_A_aux)		

	END DO

	val1 = n_A-nzA
	
	s_A_aux(nt_sub) = s_A_aux(nt_sub) - val1

	ALLOCATE(S_E(nt_sub)%At(s_A_aux(nt_sub)))
	S_E(nt_sub)%At = S_E_aux(nt_sub)%A_aux(1:s_A_aux(nt_sub))
	DEALLOCATE(S_E_aux(nt_sub)%A_aux)

	ALLOCATE(S_E(nt_sub)%col_A(s_A_aux(nt_sub)))
	S_E(nt_sub)%col_A = S_E_aux(nt_sub)%col_A_aux(1:s_A_aux(nt_sub))
	DEALLOCATE(S_E_aux(nt_sub)%col_A_aux)

	ALLOCATE(S_E(nt_sub)%row_A(s_A_aux(nt_sub)))
	S_E(nt_sub)%row_A = S_E_aux(nt_sub)%row_A_aux(1:s_A_aux(nt_sub))
	DEALLOCATE(S_E_aux(nt_sub)%row_A_aux)
    ! Matrix G --------------------------------------------------
    
	DO i=1,nt_sub-1

		ALLOCATE(S_E(i)%G_t(s_G_aux(i)))
		S_E(i)%G_t = S_E_aux(i)%G_aux
		DEALLOCATE(S_E_aux(i)%G_aux)

		ALLOCATE(S_E(i)%col_G(s_G_aux(i)))
		S_E(i)%col_G = S_E_aux(i)%col_G_aux
		DEALLOCATE(S_E_aux(i)%col_G_aux)

	END DO

	val1 = n_G-nzG
	s_G_aux(nt_sub) = s_G_aux(nt_sub) - val1

	ALLOCATE(S_E(nt_sub)%G_t(s_G_aux(nt_sub)))
	S_E(nt_sub)%G_t = S_E_aux(nt_sub)%G_aux(1:s_G_aux(nt_sub))
	DEALLOCATE(S_E_aux(nt_sub)%G_aux)

	ALLOCATE(S_E(nt_sub)%col_G(s_G_aux(nt_sub)))
	S_E(nt_sub)%col_G = S_E_aux(nt_sub)%col_G_aux(1:s_G_aux(nt_sub))
	DEALLOCATE(S_E_aux(nt_sub)%col_G_aux)

	rows_G(fi2+1) = nzG + 1

    ! Matrix V --------------------------------------------------
    IF ((Bodyforces.EQ.1).OR.(Transient.EQ.1)) THEN
        
		DO i=1,nt_sub-1

			ALLOCATE(S_E(i)%V_t(s_V_aux(i)))
			S_E(i)%V_t = S_E_aux(i)%V_aux
			DEALLOCATE(S_E_aux(i)%V_aux)

			ALLOCATE(S_E(i)%col_V(s_V_aux(i)))
			S_E(i)%col_V = S_E_aux(i)%col_V_aux
			DEALLOCATE(S_E_aux(i)%col_V_aux)		

		END DO

		val1 = n_V-nzV
		s_V_aux(nt_sub) = s_V_aux(nt_sub) - val1

		ALLOCATE(S_E(nt_sub)%V_t(s_V_aux(nt_sub)))
		S_E(nt_sub)%V_t = S_E_aux(nt_sub)%V_aux(1:s_V_aux(nt_sub))
		DEALLOCATE(S_E_aux(nt_sub)%V_aux)

		ALLOCATE(S_E(nt_sub)%col_V(s_V_aux(nt_sub)))
		S_E(nt_sub)%col_V = S_E_aux(nt_sub)%col_V_aux(1:s_V_aux(nt_sub))
		DEALLOCATE(S_E_aux(nt_sub)%col_V_aux)

		rows_V(fi3+1) = nzV + 1

    END IF 
    
    ! -------------------------------------------------------------------------
	! final size vectors

	ALLOCATE(scounts_A(nt_sub))
	scounts_A = s_A_aux
	DEALLOCATE(s_A_aux)

	ALLOCATE(scounts_G(nt_sub))
	scounts_G = s_G_aux
	DEALLOCATE(s_G_aux)
	
	IF ((Bodyforces.EQ.1).OR.(Transient.EQ.1)) THEN
	
		ALLOCATE(scounts_V(nt_sub))
		scounts_V = s_V_aux
		DEALLOCATE(s_V_aux)

	END IF 

	DEALLOCATE(S_E_aux,displs_A)

	! -------------------------------------------------------------------------

    CALL SYSTEM_CLOCK(t3)       ! Time counter stopped
    WRITE (*,'(A,F15.3,A)') 'COMPLETED! (Time :', REAL(t3-t2)/rate2,'s)'
    WRITE (*,*) ''
    
    ! -------------------------------------------------------------------------

    CALL SYSTEM_CLOCK(t1)       ! Time counter stopped
    WRITE(*,'(A,F15.4,A)') '                                  ...COMPLETED! (Time :',REAL(t1-t0)/rate,'s)'
    WRITE (*,*) ''
    ts1 = REAL(t1-t0)/rate

END SUBROUTINE Final_Sparce_Vectors
!===============================================================================!
SUBROUTINE Sparce_array(step,ID_block_reg,ID_index_reg,columns)

    INTEGER::t0,t1,rate, step
      
    INTEGER :: nelem, dof, ind, size_reg, val1, val2, el
    !REAL(8), ALLOCATABLE :: H_reg(:,:), G_reg(:,:), H_bc(:,:), G_bc(:,:)
    INTEGER :: i, j, k, element, element_A, element_B, reg_A
    INTEGER :: reg_B, ndci, ndcf, ini_col, fin_col, val3, el_aux
    REAL(8) :: Vec_aux_1(21), Vec_aux_2(9)

    INTEGER :: ind_B, size_nelem
    INTEGER, ALLOCATABLE :: ID_block_reg(:,:), ID_index_reg(:,:), columns(:,:)
    INTEGER, ALLOCATABLE :: Vec_indices(:), vec1(:), vec2(:), vec3(:) 

    CALL SYSTEM_CLOCK(t0, rate) ! Initializates the time counter
    WRITE (*,*) ' '
    WRITE (*,'(A)',advance='no') '       - 6.1 Sparce_array         ...'

    dof = 3 ! degrees of fredom
    
    ! displs of regions
    val1 = 0
    El_reg(1,5) = 0
    val2 = 0
    El_reg(1,6) = 0
    val3 = 0
    Int_reg(1,3) = 0
    DO i=1,nreg-1
        val1 = val1 + (El_reg(i,1)*nnos_el*dof)+El_reg(i,7) ! ???
        El_reg(i+1,5) = val1		
        val2 = val2 + El_reg(i,1)
        El_reg(i+1,6) = val2
        val3 = val3 + Int_reg(i,2)
        Int_reg(i+1,3) = val3
    END DO
    
    nelem = SIZE(ELEM,1) ! Number of elements
    ALLOCATE(vec1(nreg), vec2(nreg), vec3(nreg))
    vec1 = El_reg(:,1)
    vec2 = Int_reg(:,2)
	
    vec3 = vec1 + vec2
    size_nelem = MAXVAL(vec3)
    DEALLOCATE(vec1, vec2)
        
	IF (step.EQ.1) THEN

    	ALLOCATE(B_geral(nelem*nnos_el*dof,1))
		B_geral = 0.d0

	ELSEIF ((step.GT.1).AND.(Iteration.EQ.1)) THEN 

		DEALLOCATE(B_geral)

		val1 = 0
		DO i=1,nreg
			val1 = val1 + El_reg(i,8)
		END DO

		ALLOCATE(B_geral(nelem*nnos_el*dof+val1,1))
		B_geral = 0.d0;
  	
	END IF

    ALLOCATE(ID_block_reg(nreg,size_nelem+2))
    ALLOCATE(ID_index_reg(nreg,size_nelem+2))
    ALLOCATE(columns(nreg,1+size_nelem))
    ALLOCATE(Vec_indices(nreg))    

    Vec_indices = 0
    
    columns = 0
    ID_block_reg = 0; ID_index_reg = 0;
	ndcf = 0
    
    DO i=1,nreg ! Loop over regions
        ID_block_reg(i,1:2) = (/i,vec3(i)/) ! [region, N_blocks]
        ID_index_reg(i,1:2) = (/i,vec3(i)/) ! [region, N_blocks]
             
        DO j=1,El_reg(i,1) ! loop over elements per regions    
            element = Subregions(i,j) ! element
            el = ELEM(element,12) ! see if the element is an interface or CZ
			
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
				
                ndci = ndcf + 1 
                ndcf = ndcf + ELEM(element_A,14) 
	            
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
                ndci = ndcf + 1 
                ndcf = ndcf + ELEM(element,14)

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

				! Application of boundary conditions
                IF (step.EQ.1) THEN

		            CALL Apply_BC_Multidomain(H_G(i)%H,H_G(i)%G,element,ini_col)
		        

				    Vec_aux_1 = BC_elem(element,:)
				    Vec_aux_2(1:3) = (/Vec_aux_1(5),Vec_aux_1(7),Vec_aux_1(9)/)
				    Vec_aux_2(4:6) = (/Vec_aux_1(11),Vec_aux_1(13),Vec_aux_1(15)/)
				    Vec_aux_2(7:9) = (/Vec_aux_1(17),Vec_aux_1(19),Vec_aux_1(21)/)

				    B_geral(ndci:ndcf,1) = Vec_aux_2

				END IF
				
            END IF
            
        END DO
        
    END DO
    !ndci = 0
	!DO i=1,SIZE(B_geral,1)
	!	IF (B_geral(i,1).NE.0.d0) THEN
    !         
	!		WRITE(*,*) i,  ndci, B_geral(i,1)
	!	END IF
	!END DO


    DEALLOCATE(Vec_indices,vec3)
       
    ! ===========================================================

    CALL SYSTEM_CLOCK(t1)       ! Time counter stopped
    WRITE (*,'(A,F15.3,A)') 'COMPLETED! (Time :', REAL(t1-t0)/rate,'s)'
    WRITE (*,*) ''
        
END SUBROUTINE Sparce_array
!===============================================================================!
!SUBROUTINE Apply_BC_Multidomain(Hg,Gg,element,ini_col,H_bc,G_bc)
SUBROUTINE Apply_BC_Multidomain(Hg,Gg,element,ini_col)

    INTEGER :: i, dof, ind, ini_col, element, type_bc    
    REAL(8), ALLOCATABLE :: change(:,:), Hg(:,:), Gg(:,:)!, H_bc(:,:), G_bc(:,:)

    ALLOCATE(change(SIZE(Hg,1),1));
    change = 0.d0
    dof = 3

    DO i=1,nnos_el*dof ! select the type of BC 
        type_bc =  INT(BC_elem(element,2*i+2)) ! type of BC 
         
        IF ((type_bc .EQ. 0)) THEN ! Displacement is known
			
            ind = ini_col + i - 1
            change(:,1) = Gg(:,ind)
            Gg(:,ind) = -Hg(:,ind) 
            Hg(:,ind)  = -change(:,1)
        END IF
    END DO

    !H_bc = Hg; G_bc = Gg
    
END SUBROUTINE Apply_BC_Multidomain
!===============================================================================!
END MODULE General_Sparce_Vectors
