!===============================================================================!
!-------------------------------------------------------------------------------!
!                            MODULE TO COMPUTE VECTORS {b}                      !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Vectors_b
!-------------------------------------------------------------------------------!
USE Global_variables
USE BC_Module
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE b_update(ts1,step)

    INTEGER::t0,t1,rate
    REAL(8),INTENT(OUT)::ts1
    
    INTEGER :: n, step, i
    REAL(8), ALLOCATABLE :: bhat_aux(:,:), b_step_aux(:,:)

    CALL SYSTEM_CLOCK(t0, rate) ! Initializates the time counter
    WRITE (*,*) '----------------------------------------------------------------------------'
    WRITE (*,*) ' '
    WRITE (*,'(A)',advance='no') ' 7. Computing vector {b}          ...'

    ! Quasi_static and body forces analysis --------------------------------------
    IF (((Quasi_static.EQ.1).AND.(Bodyforces.EQ.0)).OR. &
       ((Quasi_static.EQ.1).AND.(Bodyforces.EQ.1))) THEN

        IF (step.GT.1) THEN
			
            IF (BC_groups.EQ.0) THEN

                CALL New_B_geral(step)
		        
            ELSE

                CALL New_B_geral_M(step)
            	
            END IF
			
        END IF
	
        ! Vector Gt
        IF (nreg.EQ. 1) THEN
            b_cte = MATMUL(G_geral,B_geral)
        ELSE
            b_cte = 0.d0; 
	    	CALL MATMUL_Mumps(1,displs_G,scounts_G,rows_G,B_geral,b_cte)
            !CALL MATMUL_sparse(G_values,rows_G,columns_G,B_geral,b_cte)
        END IF 

        ! Vector Gt + Vb
        IF (Bodyforces.EQ.1) THEN
			
            IF (step.GT.1) THEN
                CALL New_bhat(step)
            END IF
			
            IF (nreg.EQ. 1) THEN
                b_bf = MATMUL(V,bhat)
                b = b_cte + b_bf
            ELSE
                n = SIZE(bhat,1)
                ALLOCATE(bhat_aux(n,1))
                bhat_aux = bhat				
                b_bf = 0.d0				
				CALL MATMUL_Mumps(2,displs_V,scounts_V,rows_V,bhat_aux,b_bf)
                !CALL MATMUL_sparse(V_values,rows_V,columns_V,bhat_aux,bhat)
                b = b_cte + b_bf
				
                DEALLOCATE(bhat_aux)
            END IF        
        END IF  

        ! Conventional 
        IF ((Bodyforces.EQ.0).AND.(Quasi_static.EQ.1)) THEN
            b = b_cte
        END IF
        
    END IF

    ! For transient analysis ----------------------------------------------------
    ! Vector Gt + (1/Dt²)*density*V*(5u_t - 4u_t-dt + u_t-2t)
    IF (Transient.EQ.1) THEN
    
        IF (step.GT.1) THEN

            IF (BC_groups.EQ.0) THEN

                CALL New_B_geral(step)
		        
            ELSE

                CALL New_B_geral_M(step)
            	
            END IF
	    
        END IF



		IF (Bodyforces.EQ.1) THEN

			IF (step.GT.1) THEN
                CALL New_bhat(step)
            END IF

			IF (nreg.EQ. 1) THEN
                b_bf = MATMUL(V,bhat)
                
            ELSE
                n = SIZE(bhat,1)
                ALLOCATE(bhat_aux(n,1))
                bhat_aux = bhat				
                b_bf = 0.d0				
				CALL MATMUL_Mumps(2,displs_V,scounts_V,rows_V,bhat_aux,b_bf)
                !CALL MATMUL_sparse(V_values,rows_V,columns_V,bhat_aux,bhat)
                				
                DEALLOCATE(bhat_aux)
            END IF 
			

		END IF

			


        IF (nreg.EQ. 1) THEN
            b_cte = MATMUL(G_geral,B_geral)
            b_step(:,1) = (1.d0/Dt**2)*density*MATMUL(V,(5.d0*u_t(:,3)-4.d0*u_t(:,2)+u_t(:,1)))
            b = b_cte + b_step + b_bf

        ELSE
            b_cte = 0.d0; 
			CALL MATMUL_Mumps(1,displs_G,scounts_G,rows_G,B_geral,b_cte)
			
            !CALL MATMUL_sparse(G_values,rows_G,columns_G,B_geral,b_cte)
            n = SIZE(b_step,1)
            ALLOCATE(b_step_aux(n,1))
            b_step_aux = 0.d0
            b_step_aux(:,1) = (1.d0/Dt**2)*density*(5.d0*u_t(:,3)-4.d0*u_t(:,2)+u_t(:,1));            
            b_step = 0.d0
			
			CALL MATMUL_Mumps(2,displs_V,scounts_V,rows_V,b_step_aux,b_step)

            !CALL MATMUL_sparse(V_values,rows_V,columns_V,b_step_aux,b_step)
            b = b_cte + b_step + b_bf

			
            DEALLOCATE(b_step_aux)
        END IF 
        
    END IF
    !------------------------------------------------------------------------------
    CALL SYSTEM_CLOCK(t1)       ! Time counter stopped
    WRITE (*,'(A,F15.3,A)') 'COMPLETED! (Time :', REAL(t1-t0)/rate,'s)'
    WRITE (*,*) ''
    ts1 = REAL(t1-t0)/rate

END SUBROUTINE b_update
!===============================================================================!
SUBROUTINE Bodyforces_vet
    
    INTEGER :: i, ind1, ind2, Steps
    
	Steps = 1 ! Default
    IF (Transient.EQ.1) THEN
        Steps = Time_steps
    END IF

    IF (Quasi_static.EQ.1)THEN
    	Steps = Static_steps
    END IF		

	IF (type_bf.EQ.0) THEN 

		ind2 = 0
		DO i=1,SIZE(bhat,1)/3
		    
			ind1 = ind2 + 1
		    ind2 = ind2 + 3

		    bhat(ind1:ind2,1) = BodyForce
		        
		END DO

	ELSE

		OPEN(1,file=trim(fileplace_Bodyforce)//trim(Bodyforce_file)//'.dat',STATUS='OLD')
				
			ALLOCATE(Body_Forces(Steps,3))
			Body_Forces = 0.d0

			DO i = 1,Steps
				READ (1,*) Body_Forces(i,:) 
			END DO
			
		CLOSE(1)

			ind2 = 0
			DO i=1,SIZE(bhat,1)/3
				
				ind1 = ind2 + 1
				ind2 = ind2 + 3

				bhat(ind1:ind2,1) = Body_Forces(1,:)
				
			END DO
			
	END IF

		
       
END SUBROUTINE Bodyforces_vet
!===============================================================================!
SUBROUTINE New_bhat(step)

    INTEGER :: step, ind2, i, ind1
    
		
		IF (type_bf.EQ.1) THEN 

		    ind2 = 0
			DO i=1,SIZE(bhat,1)/3
				
				ind1 = ind2 + 1
				ind2 = ind2 + 3

				bhat(ind1:ind2,1) = Body_Forces(step,:)
				    
			END DO	
	
		END IF

END SUBROUTINE New_bhat
!===============================================================================!
SUBROUTINE New_B_geral(step)

    INTEGER :: i, nelem, j, m, ind, dof, Type_bc, step
    REAL(8) :: val_t

    INTEGER :: face, cont1, face_bc, nbc, el, cont2, cond
    REAL(8) :: bc(6), val_BC

    INTEGER :: element, k, ninterfaces, element_A, element_B, reg_A, reg_B
    INTEGER :: ndci, ndcf
    REAL(8) :: Vec_aux_1(21), Vec_aux_2(9)

	REAL(8) :: val, f, Tol

    ! Traction and displacement
    f = 1.d0
	Tol = 1e-7
    DO i=1,SIZE(T_D_BC)/7

        ind = (i-1)*7
        face = INT(T_D_BC(ind+1)) 
        cont1 = 1
        cont2 = 2
        
        DO j=1,3

            Type_bc = INT(T_D_BC(ind+2*j+1)) ! Displacement or traction
                
            ! Traction
            IF (Type_bc.EQ.1) THEN
        
                bc(cont1) = 1.d0
					
            ELSE ! Displacements

                bc(cont1) = 0.d0

            END IF

			bc(cont2) = BCs(i)%bc(j,step)
            

            cont1 = cont1 + 2
            cont2 = cont2 + 2

        END DO

        nbc = SIZE(bc_new,1)
                
        DO k=1,nbc
            face_bc = bc_new(k,1)
            el = bc_new(k,2)
            IF (face.EQ.face_bc) THEN
                IF (nreg.GT.1) THEN
                    BC_elem(el,4:21) = (/bc,bc,bc/)
                ELSE
                    BC_elem(el,2:19) = (/bc,bc,bc/)
                END IF
            END IF
        END DO 

    END DO

    ! If there are constrained nodes in the model ---------------------------
    cond = Config_BC(1) ! Boundary condition
    IF (cond.GT.0) THEN

        ! Apply the boundary conditions to those nodes
        CALL Constrain_node

    END IF
    !-----------------------------------------------------------------------


    !----------------------------------------------------------------
    IF (nreg.EQ.1) THEN

        nelem = El_reg(nreg,1) ! Number of elements
        dof = 3 ! degrees of fredom
        m = nnos_el*dof
      
        DO i=1,nelem
            DO j=1,m ! select the type of BC
                type_bc = INT(BC_elem(i,2*j)) ! type of BC
                IF (type_bc .LE. 1) THEN ! Traction is known
                   ind = i*m-m+j
                   val_t = BC_elem(i,2*j+1) ! Traction value 
                   B_geral(ind,1) = val_t
                END IF
            END DO
        END DO

    ELSE

        ninterfaces = SIZE(Interfaces,1)
        ndcf=0
        dof = 3 ! degrees of fredom

        DO i=1,nreg ! Loop over regions
            DO j=1,El_reg(i,1) ! loop over elements per regions    
                element = Subregions(i,j) ! element
                DO k=1,ninterfaces ! loop over interfaces
                    element_A = Interfaces(k,4) ! Element in region A
                    element_B = Interfaces(k,5) ! ELement in region B
                    IF (element.EQ.element_A) THEN 
                        reg_A = Interfaces(k,2);
                        reg_B = Interfaces(k,3);
                        IF (reg_A.LT.reg_B) THEN
                            ndci = ndcf + 1 !initial index of the block in the general matrix 
                            ndcf = ndcf + ELEM(element_A,14) !final index of the block in the general matrix 
                            EXIT
                        ELSE
                            ndci = ndcf + 1 !initial index of the block in the general matrix 
                            ndcf = ndcf + ELEM(element_A,14) !final index of the block in the general matrix 
                            EXIT
                        END IF
                    END IF
                    IF ((k.EQ.ninterfaces).AND.(element.NE.element_A).AND.(element.NE.element_B)) THEN

                        ndci = ndcf + 1 !initial index of the block in the general matrix 
                        ndcf = ndcf + ELEM(element,14) !final index of the block in the general matrix
                        
                        Vec_aux_1 = BC_elem(element,:)
                        
                        Vec_aux_2(1:3) = (/Vec_aux_1(5),Vec_aux_1(7),Vec_aux_1(9)/)
                        Vec_aux_2(4:6) = (/Vec_aux_1(11),Vec_aux_1(13),Vec_aux_1(15)/)
                        Vec_aux_2(7:9) = (/Vec_aux_1(17),Vec_aux_1(19),Vec_aux_1(21)/)

                        B_geral(ndci:ndcf,1) = Vec_aux_2

                    END IF
                END DO
            END DO
        END DO
   END IF   
	
END SUBROUTINE New_B_geral
!===============================================================================!
SUBROUTINE New_B_geral_M(step)

    INTEGER :: i, nelem, j, m, ind, dof, Type_bc, step, bc_xyz(3)
    REAL(8) :: val_t

    INTEGER :: face, face_bc, nbc, el, cond
    REAL(8) :: bc(6), val_BC

    INTEGER :: element, k, ninterfaces, element_A, element_B, reg_A, reg_B
    INTEGER :: ndci, ndcf
    REAL(8) :: Vec_aux_1(21), Vec_aux_2(9), bc_val_xyz(3), normal(3)

    ! Traction and displacement
  
    DO i=1,BC_groups

		DO j=1,BCs(i)%nel

			el = BCs(i)%el(j)

			IF (BCs(i)%type_bc(1).EQ.2) THEN ! normal displacement

				bc_xyz = (/0,0,0/)

				normal = NORMAL_VECTORS(el,:)

				bc_val_xyz = normal*BCs(i)%bc(step,1)

			ELSE IF (BCs(i)%type_bc(1).EQ.3) THEN ! normal traction

				bc_xyz = (/1,1,1/)

				normal = NORMAL_VECTORS(el,:)

				bc_val_xyz = normal*BCs(i)%bc(step,1)
				
			ELSE IF ( (BCs(i)%type_bc(1).EQ.0).OR. &
						(BCs(i)%type_bc(1).EQ.1) ) THEN

			   bc_xyz = BCs(i)%type_bc ! in Cartesian coordinates
			   bc_val_xyz = BCs(i)%bc(step,:)
			

			END IF

			DO k=1,nnos_el

        		Type_bc = bc_xyz(k)
                         
		        ! Traction
		        IF (Type_bc.EQ.1) THEN
		            
		            bc(2*k-1) = 1.d0
		                
		        ELSE ! Displacements
		            
		            bc(2*k-1) = 0.d0

		        END IF

				bc(2*k) = bc_val_xyz(k)
       		END DO

		    IF (nreg.GT.1) THEN
		        BC_elem(el,4:21) = (/bc,bc,bc/)
		    ELSE
		        BC_elem(el,2:19) = (/bc,bc,bc/)
		    END IF

		END DO

    END DO

    ! If there are constrained nodes in the model ---------------------------
    cond = Config_BC(1) ! Boundary condition
    IF (cond.GT.0) THEN

        ! Apply the boundary conditions to those nodes
        CALL Constrain_node

    END IF
    !-----------------------------------------------------------------------


    !----------------------------------------------------------------
    IF (nreg.EQ.1) THEN

        nelem = El_reg(nreg,2) ! Number of elements
        dof = 3 ! degrees of fredom
        m = nnos_el*dof
      
        DO i=1,nelem
            DO j=1,m ! select the type of BC
                type_bc = INT(BC_elem(i,2*j)) ! type of BC
                IF (type_bc .LE. 1) THEN ! Traction is known
                   ind = i*m-m+j
                   val_t = BC_elem(i,2*j+1) ! Traction value 
                   B_geral(ind,1) = val_t
                END IF
            END DO
        END DO

    ELSE

        ninterfaces = SIZE(Interfaces,1)
        ndcf=0
        dof = 3 ! degrees of fredom

        DO i=1,nreg ! Loop over regions
            DO j=1,El_reg(i,1) ! loop over elements per regions    
                element = Subregions(i,j) ! element
                DO k=1,ninterfaces ! loop over interfaces
                    element_A = Interfaces(k,4) ! Element in region A
                    element_B = Interfaces(k,5) ! ELement in region B
                    IF (element.EQ.element_A) THEN 
                        reg_A = Interfaces(k,2);
                        reg_B = Interfaces(k,3);
                        IF (reg_A.LT.reg_B) THEN
                            ndci = ndcf + 1 !initial index of the block in the general matrix 
                            ndcf = ndcf + nnos_el*dof !final index of the block in the general matrix 
                            EXIT
                        ELSE
                            ndci = ndcf + 1 !initial index of the block in the general matrix 
                            ndcf = ndcf + nnos_el*dof !final index of the block in the general matrix 
                            EXIT
                        END IF
                    END IF
                    IF ((k.EQ.ninterfaces).AND.(element.NE.element_A).AND.(element.NE.element_B)) THEN

                        ndci = ndcf + 1 !initial index of the block in the general matrix 
                        ndcf = ndcf + nnos_el*dof !final index of the block in the general matrix
                        
                        Vec_aux_1 = BC_elem(element,:)
                        
                        Vec_aux_2(1:3) = (/Vec_aux_1(5),Vec_aux_1(7),Vec_aux_1(9)/)
                        Vec_aux_2(4:6) = (/Vec_aux_1(11),Vec_aux_1(13),Vec_aux_1(15)/)
                        Vec_aux_2(7:9) = (/Vec_aux_1(17),Vec_aux_1(19),Vec_aux_1(21)/)

                        B_geral(ndci:ndcf,1) = Vec_aux_2

                    END IF
                END DO
            END DO
        END DO
   END IF 

	   

END SUBROUTINE New_B_geral_M
!===============================================================================!
SUBROUTINE MATMUL_sparse(Mat_A,rows,columns,vec_b,vec_c) ! PARDISO format

    INTEGER :: n,k,ncol_ini,ncol_fin,l,column
    REAL(8) :: Aij, bj, cj
    REAL(8), ALLOCATABLE :: Mat_A(:), vec_b(:,:), vec_c(:,:)
    INTEGER, ALLOCATABLE :: rows(:), columns(:)
    
    n = SIZE(rows)
    
    DO k=1,n-1 
        ncol_ini = rows(k)
        ncol_fin = rows(k+1)-1
        IF (ncol_ini.NE.0) THEN
            cj = 0.d0;
            DO l=ncol_ini,ncol_fin
                column = columns(l)
                Aij = Mat_A(l)
                bj = vec_b(column,1)
                cj = cj + Aij*bj
            END DO
            vec_c(k,1) = cj
        END IF
    END DO

END SUBROUTINE MATMUL_sparse
!===============================================================================!
! Matrix = 1 -> G and 2 -> V
SUBROUTINE MATMUL_Mumps(var,displs,scounts,rows,vec_b,vec_c)
 
    INTEGER :: var, n, proc, k, ncol_ini, ncol_fin, l, ind_loc, column
    REAL(8) :: Aij, bj, cj
    REAL(8), ALLOCATABLE :: vec_b(:,:), vec_c(:,:)
	INTEGER, ALLOCATABLE :: rows(:)
    INTEGER(8), ALLOCATABLE ::  scounts(:)
	INTEGER(8), ALLOCATABLE :: displs(:)
    
	
    n = SIZE(rows)
    
	proc = 1
    DO k=1,n-1 
        ncol_ini = rows(k)
        ncol_fin = rows(k+1)-1
        IF (ncol_ini.NE.0) THEN
            cj = 0.d0;
            DO l=ncol_ini,ncol_fin
	
				IF (l.GT.(displs(proc)+scounts(proc))) THEN

					proc = proc + 1
				
				END IF

				ind_loc = l - displs(proc) ! local position

				IF (var.EQ.1) THEN
					column = S_E(proc)%col_G(ind_loc)
					Aij = S_E(proc)%G_t(ind_loc)
				ELSEIF (var.EQ.2) THEN
					column = S_E(proc)%col_V(ind_loc)
					Aij = S_E(proc)%V_t(ind_loc)
				END IF

				bj = vec_b(column,1)
				cj = cj + Aij*bj

            END DO
            vec_c(k,1) = cj
        END IF
    END DO
	!WRITE(*,*) N_total
	!DO k=1,N_total
	!	IF (vec_b(k,1).NE. 0.d0) THEN
	!		WRITE(*,*) k,vec_b(k,1)
	!	END IF
	!END DO
	
END SUBROUTINE MATMUL_Mumps
!===============================================================================!
END MODULE Vectors_b
