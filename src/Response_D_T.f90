!===============================================================================!
!-------------------------------------------------------------------------------!
!                          MODULE ORGANIZE THE RESPONSE                         !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Response_D_T
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE Displacement_traction(ts1,step)

    INTEGER::t0,t1,rate
    REAL(8),INTENT(OUT)::ts1
    
    INTEGER :: nelem, dof, nnos_total, i, j, ind, ii, step
    REAL(8) :: val, type_bc
    REAL(8), ALLOCATABLE :: Disp_aux(:,:), Trac_aux(:,:), X_aux(:,:)
    INTEGER :: k, element, element_A, element_B, reg_A, reg_B, ndci, ndcf
    INTEGER :: ndfi, ndff, ndci2, ndcf2, ndfi2, ndff2, ninterfaces

	INTEGER :: Int_el, condition_int
    INTEGER :: ind_dof, no_A, no_B, iTSL, ind_no
	REAL(8) :: u_i, u_j, delta, T, Tol      

	INTEGER :: n_el, el, face_el, col_i, col_f, Int_B, k1
	REAL(8) :: Normal_face(3), n(3), f_local(3), epsilon_d
	
	INTEGER :: g_dof(9), condition_node, g_no(9), ind_g(9), ind_g2(9)
	REAL(8) :: dus1, dus2, dus, dun, ducrs1, ducrs2, ducrs, ducrn
	REAL(8) :: dus_ducrs, dun_ducrn, d_new, d, du, del

	INTEGER :: ndci_a, ndci_b, ind3, ind4, ind5, ind6, Iteration_yes

    CALL SYSTEM_CLOCK(t0, rate) ! Initializates the time counter
    WRITE (*,*) '----------------------------------------------------------------------------'
    WRITE (*,*) ' '
    WRITE (*,'(A)',advance='no') ' 9. [D] and [T] response          ...'    

    IF (nreg.EQ.1) THEN 

        nelem = El_reg(nreg,1)! Number of elements
        dof = 3 ! Degrees of freedom
        nnos_total = nelem*nnos_el ! number of nodes

        
        ALLOCATE(Disp_aux(nnos_total*dof,1),Trac_aux(nnos_total*dof,1))
        Disp_aux = 0.d0;  Trac_aux = 0.d0

        IF (step.EQ.1) THEN

            ALLOCATE(Displacement(nnos_total,4),Traction(nnos_total,4))
            Displacement = 0.d0; Traction = 0.d0

        END IF
       
        IF (Transient.EQ.1) THEN
            ALLOCATE(X_aux(nnos_total*dof,1))
            X_aux = x_BEM
        END IF

        DO i=1,nelem
            DO j=1,nnos_el*dof
                IF (BC_elem(i,2*j).LE.(1e-4)) THEN ! Displacement is known
                    val = BC_elem(i,2*j+1)
                    ind = nnos_el*dof*(i-1)+j
                    Disp_aux(ind,1) = val
                    Trac_aux(ind,1) = x_BEM(ind,1)
                    IF (Transient.EQ.1) THEN
                        X_aux(ind,1) = 0.d0
                    END IF
                ELSE ! Traction is known
                    val = BC_elem(i,2*j+1)
                    ind = nnos_el*dof*(i-1)+j
                    Trac_aux(ind,1) = val
                    Disp_aux(ind,1) = x_BEM(ind,1)
                END IF
            END DO
        END DO

        
        DO i=1,nnos_total

            Displacement(i,1) = i
            Displacement(i,2:4) = Disp_aux(3*i-2:3*i,1)

            Traction(i,1) = i
            Traction(i,2:4) = Trac_aux(3*i-2:3*i,1)

        END DO
      

        IF (Transient.EQ.1) THEN
            u_t(:,1) = u_t(:,2)
            u_t(:,2) = u_t(:,3)
            u_t(:,3) = X_aux(:,1)
        END IF

        !DO i=1,nnos_total
        !    WRITE(*,*) Displacement(i,:)
        !END DO

        !DO i=1,nnos_total
        !    WRITE(*,*) Traction(i,:)
        !END DO

		DO i = 1,nreg ! Loop over regions
		    n_el = El_reg(i,1) ! Elements per region
		    DO j = 1,n_el 
		        el = Subregions(i,j) ! element          
		        face_el = ELEM(el,4)
				
		        Normal_face = NORMAL_VECTORS(face_el,:)
		        
				!Normal vector of the deformed element
		        CALL Normal_vector(Normal_face,el,n)

		        NORMAL_VEC_Def(el,:) = n

		    END DO 
		END DO
        
        DEALLOCATE(Disp_aux,Trac_aux,x_BEM)

    ELSE

        nelem = SIZE(ELEM,1)
        dof = 3 ! Degrees of freedom

        ALLOCATE(Disp_aux(nelem*nnos_el*dof,1))
        ALLOCATE(Trac_aux(nelem*nnos_el*dof,1))
        Disp_aux = 0.d0;  Trac_aux = 0.d0

        IF (Transient.EQ.1) THEN
            ALLOCATE(X_aux(nelem*nnos_el*dof,1))
            X_aux = x_BEM
        END IF

        ninterfaces = SIZE(Interfaces,1)
        ndff=0
		
        DO i=1,nreg ! Loop over regions
            DO j=1,El_reg(i,1) ! loop over elements per regions    
                element = Subregions(i,j) ! element
                DO k=1,ninterfaces ! loop over interfaces
                    element_A = Interfaces(k,4) ! Element in region A
                    element_B = Interfaces(k,5) ! ELement in region B
                    IF (element.EQ.element_A) THEN 

                        reg_A = Interfaces(k,2);
                        reg_B = Interfaces(k,3);
						
						Int_el = ELEM(element_A,13)
						condition_int = Interfaces(Int_el,18)

						k1 = 0
                        !-------------------------------------------------------------------------
                        ! DISPLACEMENT COMPATIBILTY IN THE INTERFACES regA<regB
                        !-------------------------------------------------------------------------
                        IF (reg_A.LT.reg_B) THEN 
                            
                            ! Index
                            ndci = ELEM(element_A,8) ! index of block Fij 
                            ndcf = ELEM(element_A,9) ! index of block Fij
                                                  
                            ndfi = ndff + 1  
                            ndff = ndff + ELEM(element_A,14)
							
                            ! Index to use in the region B of the interface
                            ELEM(element_A,10:11) = (/ndfi,ndff/) 
                            
                            ! Displacements in the region A of the interface or CZ 
                            Disp_aux(ndfi:ndff,1) = x_BEM(ndci:ndcf,1)
							
                            EXIT
                        !-------------------------------------------------------------------------
                        ! TRACTION EQUILIBRIUM IN THE INTERFACES regA>reg_B
                        !-------------------------------------------------------------------------      
                        ELSE 
							
                            ! Index of block Gij
                            ndci = ELEM(element_A,8)   
                            ndcf = ELEM(element_A,9)                       
                            
							! Index of block Fji
		                    ndci2 = ELEM(element_B,8) 
		                    ndcf2 = ELEM(element_B,9)

							! index of block Gji
		                    ndfi2 = ELEM(element_B,10) 
		                    ndff2 = ELEM(element_B,11) 

							IF (condition_int.EQ.0) THEN ! Pristine interface
								 
								ndfi = ndff + 1  
                            	ndff = ndff + ELEM(element_A,14) 
								
								! Tractions in the region A of the interface
		                        Trac_aux(ndfi:ndff,1) = x_BEM(ndci:ndcf,1)

								IF (Transient.EQ.1) THEN
		                            X_aux(ndci:ndcf,1) = 0.d0
		                        END IF

								! Displacements in the region B of the interface 
		                        Disp_aux(ndfi:ndff,1) = x_BEM(ndci2:ndcf2,1)

								! Tractions in the region B of the interface
		                        Trac_aux(ndfi2:ndff2,1) = -x_BEM(ndci:ndcf,1)

								IF (Transient.EQ.1) THEN
		                            X_aux(ndci:ndcf,1) = 0.d0
		                        END IF

							ELSE ! Interface in cohesive, contact or separation condition
								
								ndfi = ndff   
                            	ndff = ndff + ELEM(element_B,14) 
								
								DO ii=1,ELEM(element_B,14) 

									ndfi = ndfi + 1
									
									condition_node = Interfaces(Int_el,8+ii)									
	
									IF ((condition_node.EQ.0).OR.(condition_node.EQ.2)) THEN ! pristine dof or contact du3=0
										
										! Tractions in the region B of the interface
		                        		Trac_aux(ndfi,1) = x_BEM(ndci+ii-1,1)
										
										IF (Transient.EQ.1) THEN
						                	X_aux(ndci+ii-1,1) = 0.d0
						                END IF

										! Displacements in the region B of the interface 
		                        		Disp_aux(ndfi,1) = x_BEM(ndci2+ii-1,1)

										! Tractions in the region A of the interface
								        Trac_aux(ndfi2+ii-1,1) = -x_BEM(ndci+ii-1,1)

										IF (Transient.EQ.1) THEN
								        	X_aux(ndci+ii-1,1) = 0.d0
								        END IF
										
									ELSEIF (condition_node.EQ.1) THEN! Separation
								

										! Tractions in the region A of the interface
		                        		Trac_aux(ndfi,1) = 0.d0

										IF (Transient.EQ.1) THEN
								        	X_aux(ndci+ii-1,1) = 0.d0
								        END IF

										! Displacements in the region B of the CZ				
										Disp_aux(ndfi,1) = x_BEM(ndci+ii-1,1)


										! Tractions in the region B of the interface
										Trac_aux(ndfi2+ii-1,1) = 0.d0

										IF (Transient.EQ.1) THEN
								        	X_aux(ndci+ii-1,1) = 0.d0
								        END IF

										
									ELSEIF (condition_node.EQ.3) THEN! Contact

										


									END IF


								END DO

							END IF

							EXIT
						END IF

					END IF
                    !-------------------------------------------------------------------------
                    ! BLOCKS THAT BELONG TO BOUNDARY SEGMENTS OF THE GEOMETRY WITH KNOWN
                    ! BOUNDARY CONDITIONS
                    !-------------------------------------------------------------------------
                    IF ((k.EQ.ninterfaces).AND.(element.NE.element_A).AND.(element.NE.element_B)) THEN    
                        
                        ! Index
                        ndci = ELEM(element,8)  
                        ndcf = ELEM(element,9)                       
                        ndfi = ndff + 1  
                        ndff = ndff + ELEM(element,14)
                        ! Boundary conditions          
                        DO ii=1,ELEM(element,14) 
                            type_bc = BC_elem(element,2*ii+2)
                            IF (type_bc .LE. 1e-4) THEN
                                val = BC_elem(element,2*ii+3)
                                ind = ndfi+ii-1
                                Disp_aux(ind,1) = val
                                Trac_aux(ind,1) = x_BEM(ndci+ii-1,1)
                                IF (Transient.EQ.1) THEN
                                    X_aux(ndci+ii-1,1) = 0.d0
                                END IF
                            ELSE 
                                val = BC_elem(element,2*ii+3)
                                ind = ndfi+ii-1
                                Trac_aux(ind,1) = val
                                Disp_aux(ind,1) = x_BEM(ndci+ii-1,1)         
                            END IF
                        END DO

                    END IF
                END DO
            END DO
        END DO

        nnos_total = nelem*nnos_el ! number of nodes
        
        IF (step.EQ.1) THEN
            ALLOCATE(Displacement(nnos_total,4),Traction(nnos_total,4)) 
	     
        END IF

        DO i=1,nnos_total

            Displacement(i,1) = i
            Displacement(i,2:4) = Disp_aux(3*i-2:3*i,1)
			
            Traction(i,1) = i
            Traction(i,2:4) = Trac_aux(3*i-2:3*i,1)

        END DO

		! Normal vectors in the deformed position --------------------------------------------
		
		DO i = 1,nreg ! Loop over regions
		    n_el = El_reg(i,1) ! Elements per region
		    DO j = 1,n_el 
		        el = Subregions(i,j) ! element          
		        face_el = ELEM(el,4)
				
		        Normal_face = NORMAL_VECTORS(face_el,:)
		        
				!Normal vector of the deformed element
		        CALL Normal_vector(Normal_face,el,n)

		        NORMAL_VEC_Def(el,:) = n

		    END DO 
		END DO
		

		! Compute traction at CZ using the CZ equation =======================================

		IF (Failure.EQ.1) THEN

			epsilon_d = (1e-4)*1.d0
			
			Iteration = 0
			Iteration_yes = 0
			contact = 0

			Tol = 1e-7

			ind_g = (/1,2,3,1,2,3,1,2,3/)
			ind_g2 = (/1,1,1,2,2,2,3,3,3/)
			
			DO i=1,ninterfaces

				condition_int = Interfaces(i,18)		
				iTSL = Interfaces(i,8)
				reg_A = Interfaces(i,2);
		        reg_B = Interfaces(i,3);

				IF ((condition_int.EQ.1).AND.(reg_A.LT.reg_B)) THEN

					element_A = Interfaces(i,4)
					element_B = Interfaces(i,5)

					ndci_a = ELEM(element_A,8)-1
					ndci_b = ELEM(element_B,8)-1
					
					DO j=1,nnos_el*dof

						ind_no = ind_g2(j)

						no_A = ELEM_GEO_NODES_global(element_A,ind_no+1)

						no_B = ELEM_GEO_NODES_global(element_B,ind_no+1)

						ind_dof = ind_g(j)

						u_i = Displacement(no_A,ind_dof+1)

						u_j = Displacement(no_B,ind_dof+1)

						Cohesive(element_A)%delta(j) = u_j-u_i
						Cohesive(element_B)%delta(j) = u_j-u_i

					END DO

					DO j=1,nnos_el

						condition_node = Interfaces(i,8+j*dof)

						IF (Cohesive(element_A)%delta(j*dof).LT.0.d0) THEN ! from separated to contact

							Iteration = 0
							contact = 1
					
							Int_el = ELEM(element_B,13)
							Interfaces(i,8+j*dof) = 2
							Interfaces(Int_el,8+j*dof) = 2					

						ELSEIF (condition_node.EQ.2) THEN ! from contacto to separation 

							no_A = ELEM_GEO_NODES_global(element_A,j+1)
							T = traction(no_A,dof+1)
							
							IF (T.GT.0.d0) THEN

								Iteration = 0
								contact = 1

								Int_el = ELEM(element_B,13)
								Interfaces(i,8+j*dof) = 1
								Interfaces(Int_el,8+j*dof) = 1

								IF (Transient.EQ.1) THEN

									ind3 = ndci_a + 1; 
									ind4 = ndci_a + nnos_el*dof
									ind5 = ndci_b + 1; 
									ind6 = ndci_b + nnos_el*dof
									u_t(ind5:ind6,:) = u_t(ind3:ind4,:)

								END IF

							END IF

						END IF	

					END DO

				END IF

			END DO

			IF (Iteration.EQ.1) THEN
				WRITE(*,*) ''
				WRITE(*,*) ''
				WRITE (*,'(A,I5,A)') '     -- Interface in separation or contact'	
				WRITE(*,*) ''
			END IF
							
			IF (Iteration.EQ.0) THEN

				g_dof = (/1,2,3,1,2,3,1,2,3/)
				g_no = (/1,1,1,2,2,2,3,3,3/)

				DO i=1,ninterfaces

				condition_int = Interfaces(i,18)		
				iTSL = Interfaces(i,8)
				reg_A = Interfaces(i,2);
		        reg_B = Interfaces(i,3);

					IF ((condition_int.EQ.1).AND.(reg_A.LT.reg_B)) THEN

						element_A = Interfaces(i,4)
						element_B = Interfaces(i,5)

						DO j=1,nnos_el
				
							no_A = ELEM_GEO_NODES_global(element_A,j+1)	

							f_local = Displacement(no_A,2:4)

							Displacement(no_A,2:4) = MATMUL(TSL(iTSL)%R_inv,f_local)
											
							f_local = Traction(no_A,2:4)

							Traction(no_A,2:4) = MATMUL(TSL(iTSL)%R_inv,f_local)

							! ---------------------------------------------------------------

							no_B = ELEM_GEO_NODES_global(element_B,j+1)	

							f_local = Displacement(no_B,2:4)

							Displacement(no_B,2:4) = MATMUL(TSL(iTSL)%R_inv,f_local)	

							f_local = Traction(no_B,2:4)

							Traction(no_B,2:4) = MATMUL(TSL(iTSL)%R_inv,f_local)
						
						END DO

					END IF

				END DO

			END IF

			IF (Transient.EQ.1) THEN
		        u_t(:,1) = u_t(:,2)
		        u_t(:,2) = u_t(:,3)
		        u_t(:,3) = X_aux(:,1)
		    END IF

		END IF
		
		! ====================================================================================
		
        
        IF ((Transient.EQ.1).AND.(Failure.EQ.0)) THEN
            u_t(:,1) = u_t(:,2)
            u_t(:,2) = u_t(:,3)
            u_t(:,3) = X_aux(:,1)
        END IF

        !DO i=1,nnos_total
        !    WRITE(*,*) Displacement(i,:)
        !END DO

        !DO i=1,nnos_total
        !    WRITE(*,*) Traction(i,:)
        !END DO

        DEALLOCATE(Disp_aux,Trac_aux,x_BEM)

    END IF
    
    CALL SYSTEM_CLOCK(t1)       ! Time counter stopped
	IF (Iteration.EQ.1) THEN
    	WRITE (*,'(A,F15.3,A)') '                                  ...COMPLETED! (Time :', REAL(t1-t0)/rate,'s)'
	ELSE
		WRITE (*,'(A,F15.3,A)') 'COMPLETED! (Time :', REAL(t1-t0)/rate,'s)'
	END IF
	
    WRITE (*,*) ''
    ts1 = REAL(t1-t0)/rate

END SUBROUTINE Displacement_traction
!===============================================================================!
SUBROUTINE Normal_vector(Nf,el,n)

    INTEGER :: el, no1, no2, no3
    REAL(8) :: x1, y1, z1, u1, v1, w1, x2, y2, z2, u2, v2, w2
    REAL(8) :: x3, y3, z3, u3, v3, w3, a(3), b(3), d(3), n(3)
    REAL(8) :: Nf(3), d_dot_Nf(3), vec

    ! Node 1 and coordinates ------------------
    no1 = ELEM_GEO_NODES_global(el,2)
    x1 = PHY_NODES_global(no1,2); 
    y1 = PHY_NODES_global(no1,3);         
    z1 = PHY_NODES_global(no1,4);
    ! Diplacement
    u1 = Displacement(no1,2); ! dx
    v1 = Displacement(no1,3); ! dy
    w1 = Displacement(no1,4); ! dz
    ! Current position
    x1 = x1 + u1; y1 = y1 + v1; z1 = z1 + w1;
                     
    ! Node 2 and coordinates ------------------
    no2 = ELEM_GEO_NODES_global(el,3)   
    x2 = PHY_NODES_global(no2,2); 
    y2 = PHY_NODES_global(no2,3);         
    z2 = PHY_NODES_global(no2,4);
    ! Diplacement
    u2 = Displacement(no2,2); ! dx
    v2 = Displacement(no2,3); ! dy
    w2 = Displacement(no2,4); ! dz
    ! Current position
    x2 = x2 + u2; y2 = y2 + v2; z2 = z2 + w2;

    ! Node 2 and coordinates --------------------
    no3 = ELEM_GEO_NODES_global(el,4)     
    x3 = PHY_NODES_global(no3,2); 
    y3 = PHY_NODES_global(no3,3);         
    z3 = PHY_NODES_global(no3,4); 
    ! Diplacement
    u3 = Displacement(no3,2); ! dx
    v3 = Displacement(no3,3); ! dy
    w3 = Displacement(no3,4); ! dz
    ! Current position
    x3 = x3 + u3; y3 = y3 + v3; z3 = z3 + w3;      

    ! Vector product d = a x b
    a = (/x1,y1,z1/) - (/x3,y3,z3/)
    b = (/x2,y2,z2/) - (/x3,y3,z3/)
    d = (/a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)/)
    n = d / DSQRT(d(1)**2 + d(2)**2 + d(3)**2)
    
    ! Projection of d over Normal_face
    vec = n(1)*Nf(1) + n(2)*Nf(2) + n(3)*Nf(3)
    d_dot_Nf = vec*Nf
    
    IF ((d_dot_Nf(1).LT.0).AND.(Nf(1).GT.0)) THEN
        n(1) = -n(1)
    END IF
    IF ((d_dot_Nf(1).GT.0).AND.(Nf(1).LT.0)) THEN
        n(1) = -n(1)
    END IF
    IF ((d_dot_Nf(2).LT.0).AND.(Nf(2).GT.0)) THEN
        n(2) = -n(2)
    END IF
    IF ((d_dot_Nf(2).GT.0).AND.(Nf(2).LT.0)) THEN
        n(2) = -n(2)
    END IF
    IF ((d_dot_Nf(3).LT.0).AND.(Nf(3).GT.0)) THEN
        n(3) = -n(3)
    END IF
    IF ((d_dot_Nf(3).GT.0).AND.(Nf(3).LT.0)) THEN
        n(3) = -n(3)
    END IF

END SUBROUTINE Normal_vector
!===============================================================================!
END MODULE Response_D_T
