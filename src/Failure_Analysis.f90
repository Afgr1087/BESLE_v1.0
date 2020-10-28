!===============================================================================!
!-------------------------------------------------------------------------------!
!                     MODULE FAILURE ANALYSIS
!
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Failure_Analysis
!-------------------------------------------------------------------------------!
USE Global_variables
USE Global_functions
!USE mpi
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!

!   1. Pre_cracks

!===============================================================================!
!===============================================================================!
!   1. Pre_cracks
!===============================================================================!
SUBROUTINE Pre_cracks
    
    INTEGER :: ninterfaces, i, j, el_A, no, cont, val, el_B, reg_A, ncracks, nelem
    INTEGER, ALLOCATABLE :: Interfaces_aux(:,:)
    REAL(8) :: xf1, xf2, xf3, yf1, yf2, yf3, zf1, zf2, zf3
    REAL(8), ALLOCATABLE :: i_Cracks(:,:)

        
    Test_name = "Input_cracks.dat"
    OPEN(1,file=fileplace7//Test_name,STATUS='OLD')

    REWIND(1)

    READ (1,*) ncracks

    ALLOCATE(i_Cracks(ncracks,6))
    DO i = 1,ncracks
        READ (1,*) i_Cracks(i,:)
    END DO

	nelem = SIZE(ELEM,1)

	ALLOCATE(Damage_zones(nelem),Damage_coeff(nelem,3)) 
	Damage_zones = 0
	Damage_coeff = 0.d0 

    CLOSE(1) 

    ninterfaces = SIZE(Interfaces,1)

    ALLOCATE(Interfaces_aux(ninterfaces,17))

    Interfaces_aux = Interfaces
    cont = 0
    DO i=1,ncracks
        
        xmin_f = i_Cracks(i,1);  xmax_f = i_Cracks(i,2); 
        ymin_f = i_Cracks(i,3);  ymax_f = i_Cracks(i,4); 
        zmin_f = i_Cracks(i,5);  zmax_f = i_Cracks(i,6);           
        
        DO j=1,ninterfaces

            el_A = Interfaces(j,4)
                
            no = ELEM(el_A,1); 
            xf1 = POINTS(no,1); yf1 = POINTS(no,2); zf1 = POINTS(no,3);

            no = ELEM(el_A,2); 
            xf2 = POINTS(no,1); yf2 = POINTS(no,2); zf2 = POINTS(no,3);

            no = ELEM(el_A,3); 
            xf3 = POINTS(no,1); yf3 = POINTS(no,2); zf3 = POINTS(no,3);
                    
                
            IF ((xf1.GT.xmin_f).AND.(xf1.LT.xmax_f).AND. &
                (yf1.GT.ymin_f).AND.(yf1.LT.ymax_f).AND. &
                (zf1.GT.zmin_f).AND.(zf1.LT.zmax_f).AND. &
                (xf2.GT.xmin_f).AND.(xf2.LT.xmax_f).AND. &
                (yf2.GT.ymin_f).AND.(yf2.LT.ymax_f).AND. &
                (zf2.GT.zmin_f).AND.(zf2.LT.zmax_f).AND. &
                (xf3.GT.xmin_f).AND.(xf3.LT.xmax_f).AND. &
                (yf3.GT.ymin_f).AND.(yf3.LT.ymax_f).AND. &
                (zf3.GT.zmin_f).AND.(zf3.LT.zmax_f)) THEN    
            
                Interfaces_aux(j,:) = 0                          
                cont = cont + 1
				Damage_zones(el_A) = el_A
                Damage_coeff(el_A,:) = 1.d0

            END IF
                
        END DO

    END DO
	Damage_elements = cont

	WRITE(*,*) ''
	WRITE(*,*) ''
	IF (ncracks.EQ.1) THEN
    	WRITE (*,'(A,I6,A)') '     -- Deleting ',cont/2 ,'   interfaces to create the crack'
	ELSEIF (ncracks.GT.1) THEN
		WRITE (*,'(A,I6,A)') '     -- Deleting ',cont/2 ,'   interfaces to create the cracks'
	END IF
    WRITE(*,*) ''

    Interfaces = Interfaces_aux
    Interfaces_aux = 0;
    cont = 0
    DO i=1,ninterfaces
        val = Interfaces(i,1)
        IF (val.NE.0) THEN
            cont = cont + 1
            Interfaces_aux(cont,1) = cont
            Interfaces_aux(cont,:) = Interfaces(i,:)
                    
        END IF
    END DO
    ninterfaces = cont

    DEALLOCATE(Interfaces)
    ALLOCATE(Interfaces(cont,20))
	
	Interfaces = 0
    Interfaces(1:cont,1:18) = Interfaces_aux(1:cont,1:18)
    DEALLOCATE(Interfaces_aux)

    ELEM(:,12:13) = 0;
    Int_reg = 0
            
    DO i=1,ninterfaces
                
        el_A = Interfaces(i,4)
        el_B = Interfaces(i,5)
                
        ! Indentification of interfaces in the matrix of elements
        ELEM(el_A,12:13) = (/el_B,i/)
                
        reg_A = Interfaces(i,2)  	      
        Int_reg(reg_A,1) = reg_A
        cont = Int_reg(reg_A,2)
        cont = cont + 1
        Int_reg(reg_A,2) = cont
    END DO        
    
    DEALLOCATE(i_Cracks)

END SUBROUTINE Pre_cracks
!===============================================================================!
SUBROUTINE Read_TSL

	INTEGER :: i,m
	CHARACTER(LEN=20)::Test_name

	INTEGER :: iinterface, face, el_a
	REAL(8) :: val_1(5)
	REAL(8) :: n_v(3), x1(3), ang_x, ang_z, theta, Tol
	
	REAL(8) :: R_d(3,3)
	INTEGER :: M1, N, LDA, INFO, LWORK
    INTEGER, ALLOCATABLE :: IPIV(:)
    REAL(8), ALLOCATABLE :: WORK(:)
		
    WRITE (*,'(A,I3,A)') '     -- Reading critical energies for failure analysis'
    WRITE(*,*) ''


    Test_name = "Ec.txt"
    OPEN(1,file=fileplace7//Test_name,STATUS='OLD')
    
	REWIND(1)
    
	Tol = (1e-7)*1.d0
	DO i=1,n_TSL

		! From local to global coordinates-------------------------------------------------
		
		iinterface = TSL(i)%GBs(1)
		face = Interfaces(iinterface,6)
		n_v = NORMAL_VECTORS(face,:)

		theta = DACOS(n_v(3))
		
		IF (theta.LT.Tol) THEN ! If theta = 0

			ang_z = 0.d0; ang_x = 0.d0

		ELSEIF (ABS(theta-Pi).LT.Tol) THEN ! If theta = pi

			ang_z = 0.d0; ang_x = Pi

		ELSEIF ((ABS(theta-Pi).NE.Tol).AND. &
									(theta.NE.0.d0)) THEN ! If No

			x1 = (/-n_v(2),n_v(1),0.d0/) ! z x N

			CALL Atan_2(x1(1),x1(2),ang_z)
		
			CALL Atan_2(n_v(3),DSQRT(n_v(1)**2+n_v(1)**2),ang_x)

		END IF
		
		TSL(i)%theta_z_x(1) = ang_z
		TSL(i)%theta_z_x(2) = ang_x
		
		R_d(1,1) = DCOS(ang_z)*DCOS(ang_z)-DCOS(ang_x)*DSIN(ang_z)*DSIN(ang_z)
		R_d(1,2) = DCOS(ang_z)*DSIN(ang_z)+DCOS(ang_x)*DCOS(ang_z)*DSIN(ang_z)
		R_d(1,3) = DSIN(ang_z)*DSIN(ang_x)
		R_d(2,1) = -DSIN(ang_z)*DCOS(ang_z)-DCOS(ang_x)*DSIN(ang_z)*DCOS(ang_z)
		R_d(2,2) = -DSIN(ang_z)*DSIN(ang_z)+DCOS(ang_x)*DCOS(ang_z)*DCOS(ang_z)
		R_d(2,3) = DCOS(ang_z)*DSIN(ang_x)
		R_d(3,1) = DSIN(ang_x)*DSIN(ang_z)
		R_d(3,2) = -DSIN(ang_x)*DCOS(ang_z)
		R_d(3,3) = DCOS(ang_x)

		TSL(i)%R = R_d

        M1 = SIZE(R_d,1)
        N = SIZE(R_d,2)
        LDA = M1
		
        ALLOCATE(IPIV(M1))
        IPIV = 0
        INFO = 0
        CALL DGETRF(M1, N, R_d, LDA, IPIV, INFO)
        
        LWORK = N
        ALLOCATE(WORK(LWORK))
        WORK = 0.d0
        CALL DGETRI(N, R_d, LDA, IPIV, WORK, LWORK, INFO)

		DEALLOCATE(IPIV,WORK)

		TSL(i)%R_inv = R_d

		! ---------------------------------------------------------------------------------
		! TSL values in the LCS -----------------------------------------------------------
	
		READ (1,'(5F30.15)') val_1
		
		val_1 = val_1/scale_prop_mat

		IF (i.NE.adhesive) THEN
		
			val_1 = val_1*(1e4)

		END IF
		
		! Critical energies 
		! 4 shear energies corresponding to each quadrant of the local plano x,y
		! 1 normal energy density
		! Es1, Es2, Es3, Es4, En
		TSL(i)%E0c = val_1 
			
	END DO
		

	CLOSE(1)

END SUBROUTINE Read_TSL
!===============================================================================!
SUBROUTINE Local_fields(nt,me,time,step)
	include 'mpif.h'
    
	INTEGER :: nt,me, mpierr, root=0, step
    REAL(8) :: t1, t2, duration, time

	INTEGER :: ninterfaces, nelem, i, dof, cont, reg_a, reg_b, el_a, el_b
	INTEGER :: iTSL, j, no, ind1, ind2, ind_i_a, ind_i_b, n_dof_a, n_dof_b
	INTEGER :: ndr_a, ndc_a, ndr_b, ndc_b, col_i, col_f, condition_node
	INTEGER :: condition, ndci_a, ndci_b
	REAL(8) :: eta_load, R_d(3,3), S_lk(3,3), S_lk_1(3,3), e_lk(3,3), e_lk_1(3,3)
	REAL(8) :: Eks1, Eks2, Ekn, Es, Es_crit, alpha, beta, En, En_crit
	REAL(8) :: Tol
	REAL(8), ALLOCATABLE :: H_aux1(:,:), G_aux1(:,:), M_aux1(:,:)
	REAL(8), ALLOCATABLE :: H_aux2(:,:), G_aux2(:,:), eta_load_aux(:), M_aux2(:,:)
	INTEGER :: ind3,ind4,ind5,ind6, reset

    t1 = MPI_Wtime();

    CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr);
	
	IF (me.EQ.0) THEN    

		WRITE (*,*) '----------------------------------------------------------------------------'
		WRITE (*,*) ''
        WRITE (*,'(A)',advance='no') ' 13. Checking interfaces'
		
		ninterfaces = SIZE(Interfaces,1)
		nelem = SIZE(ELEM,1)

		IF (step.EQ.1) THEN

			ALLOCATE(Cohesive(nelem),Stress_strain_local(nelem))
			
			DO i = 1,nelem
				Cohesive(i)%Es1 = 0.d0
				Cohesive(i)%Es2 = 0.d0
				Cohesive(i)%En = 0.d0
				Cohesive(i)%alloc = 0
				Cohesive(i)%Criterion = 0.d0
				Stress_strain_local(i)%s_l = 0.d0
				Stress_strain_local(i)%e_l = 0.d0

			END DO

			IF (Crack.EQ.0) THEN

				ALLOCATE(Damage_zones(nelem),Damage_coeff(nelem,3)) 
				Damage_zones = 0
				Damage_coeff = 0.d0 

			END IF

		END IF

		dof = 3

		eta_load = 0.d0
		El_reg(:,7:8) = 0
		Interfaces(:,19:20) = 0
		ELEM(:,14) = nnos_el*dof
		Iteration=0
        cont = 0
		Tol = (1e-7)*1.d0

		ALLOCATE(eta_load_aux(ninterfaces))
		eta_load_aux = 0.d0
		
		DO i=1,ninterfaces
			
		    reg_a = Interfaces(i,2)
		    reg_b = Interfaces(i,3)    
		
		    el_a = Interfaces(i,4)
		    el_b = Interfaces(i,5)
					
			ndci_a = ELEM(el_a,8)-1
			ndci_b = ELEM(el_b,8)-1

			iTSL = Interfaces(i,8)

			R_d = TSL(iTSL)%R

			DO j=1,nnos_el

				IF (reg_a.LT.reg_b) THEN

					no = ELEM_GEO_NODES_global(el_a,j+1)

					! sigma_global ---------------------------------------------------------
					S_lk(1,:) = (/Stress(no,1), Stress(no,6), Stress(no,5)/)
					S_lk(2,:) = (/Stress(no,6), Stress(no,2), Stress(no,4)/)
					S_lk(3,:) = (/Stress(no,5), Stress(no,4), Stress(no,3)/)	

					! sigma_local at k
					S_lk = MATMUL(S_lk,TRANSPOSE(R_d)) 
					S_lk = MATMUL(R_d,S_lk)

					! sigma_local at k-1
					S_lk_1 = Stress_strain_local(el_a)%S_l(:,3*j-2:3*j)

					! epsilon_global ---------------------------------------------------------
					e_lk(1,:) = (/Strain(no,1), Strain(no,6), Strain(no,5)/)
					e_lk(2,:) = (/Strain(no,6), Strain(no,2), Strain(no,4)/)
					e_lk(3,:) = (/Strain(no,5), Strain(no,4), Strain(no,3)/)	

					! sigma_local at k
					e_lk = MATMUL(e_lk,TRANSPOSE(R_d)) 
					e_lk = MATMUL(R_d,e_lk)

					! sigma_local at k-1
					e_lk_1 = Stress_strain_local(el_a)%e_l(:,3*j-2:3*j)
	
					! Strain energy density 
					Eks1 = ((S_lk(1,3)+S_lk_1(1,3))*(e_lk(1,3)-e_lk_1(1,3)))/2.d0
					Cohesive(el_a)%Es1(j) = Cohesive(el_a)%Es1(j) + Eks1
					Cohesive(el_b)%Es1(j) = Cohesive(el_b)%Es1(j) + Eks1
						
					Eks2 = ((S_lk(2,3)+S_lk_1(2,3))*(e_lk(2,3)-e_lk_1(2,3)))/2.d0
					Cohesive(el_a)%Es2(j) = Cohesive(el_a)%Es2(j) + Eks2
					Cohesive(el_b)%Es2(j) = Cohesive(el_b)%Es2(j) + Eks2
				
					Ekn = ((S_lk(3,3)+S_lk_1(3,3))*(e_lk(3,3)-e_lk_1(3,3)))/2.d0	
					Cohesive(el_a)%En(j) = Cohesive(el_a)%En(j) + Ekn
					Cohesive(el_b)%En(j) = Cohesive(el_b)%En(j) + Ekn

					! Update
					Stress_strain_local(el_a)%S_l(:,3*j-2:3*j) = S_lk
					Stress_strain_local(el_b)%S_l(:,3*j-2:3*j) = S_lk
					Stress_strain_local(el_a)%e_l(:,3*j-2:3*j) = e_lk
					Stress_strain_local(el_b)%e_l(:,3*j-2:3*j) = e_lk

				END IF
				
				! Generalized energy failure criterion ==================================================
				! select the more critical shear 

				S_lk = Stress_strain_local(el_a)%S_l(:,3*j-2:3*j)
				Eks1 = Cohesive(el_a)%Es1(j)
				Eks2 = Cohesive(el_a)%Es2(j)

				Es = 0.d0
				Es_crit = 1.d0 ! just for calculations

				IF ((S_lk(1,3).GT.0.d0).AND.(S_lk(2,3).GT.0.d0)) THEN

					alpha = Eks1/TSL(iTSL)%E0c(1)
					beta = Eks2/TSL(iTSL)%E0c(3)

					IF (alpha.GT.beta) THEN
						Es = Eks1
						Es_crit = TSL(iTSL)%E0c(1)
					ELSE
						Es = Eks2
						Es_crit = TSL(iTSL)%E0c(3)
					END IF

				ELSEIF ((S_lk(1,3).GT.0.d0).AND.(S_lk(2,3).LT.0.d0)) THEN

					alpha = Eks1/TSL(iTSL)%E0c(1)
					beta = Eks2/TSL(iTSL)%E0c(4)

					IF (alpha.GT.beta) THEN
						Es = Eks1
						Es_crit = TSL(iTSL)%E0c(1)
					ELSE
						Es = Eks2
						Es_crit = TSL(iTSL)%E0c(4)
					END IF

				ELSEIF ((S_lk(1,3).LT.0.d0).AND.(S_lk(2,3).LT.0.d0)) THEN

					alpha = Eks1/TSL(iTSL)%E0c(2)
					beta = Eks2/TSL(iTSL)%E0c(4)

					IF (alpha.GT.beta) THEN
						Es = Eks1
						Es_crit = TSL(iTSL)%E0c(2)
					ELSE
						Es = Eks2
						Es_crit = TSL(iTSL)%E0c(4)
					END IF

				ELSEIF ((S_lk(1,3).LT.0.d0).AND.(S_lk(2,3).GT.0.d0)) THEN
							
					alpha = Eks1/TSL(iTSL)%E0c(2)
					beta = Eks2/TSL(iTSL)%E0c(3)

					IF (alpha.GT.beta) THEN
						Es = Eks1
						Es_crit = TSL(iTSL)%E0c(2)
					ELSE
						Es = Eks2
						Es_crit = TSL(iTSL)%E0c(3)
					END IF

				END IF

				! critical traction
				En = 0.d0
				En_crit = 1.d0 ! just for calculations
				IF (S_lk(3,3).GT.0.d0) THEN

					En = Cohesive(el_a)%En(j)
					En_crit = TSL(iTSL)%E0c(5)

				ELSE

					En = 0.d0
					En_crit = TSL(iTSL)%E0c(5)

				END IF

				Cohesive(el_a)%Criterion(j) = (Es/Es_crit) + (En/En_crit)
				Damage_coeff(el_a,j) = Cohesive(el_a)%Criterion(j)

				ind1 = 8+(3*j-2); ind2 = 8+(3*j)
			    condition_node = SUM(Interfaces(i,ind1:ind2)) ! Only analyse pristine nodes

				IF (((Cohesive(el_a)%Criterion(j).GT.1.d0).OR. &
								(ABS(Cohesive(el_a)%Criterion(j)-1.d0).LT.Tol)).AND. &
																(condition_node.EQ.0)) THEN

					Interfaces(i,ind1:ind2) = 1
					Interfaces(i,18) = 1
					Iteration=1

					cont = cont + 1

					IF ((reg_a.LT.reg_b).AND.(Transient.EQ.1)) THEN

						ind3 = ndci_a + 1; 
						ind4 = ndci_a + nnos_el*dof
						ind5 = ndci_b + 1; 
						ind6 = ndci_b + nnos_el*dof
						u_t(ind5:ind6,:) = u_t(ind3:ind4,:)
						
						
					END IF

					Damage_zones(el_a) = el_a
					Damage_coeff(el_a,j) = 1.d0
					
					Cohesive(el_a)%Criterion(j) = 1.d0

				ELSEIF (condition_node.GT.0) THEN

					Cohesive(el_a)%Criterion(j) = 0.d0

					Damage_coeff(el_a,j) = 1.d0

				ENDIF 

				! ============================================================================================

			END DO

			eta_load_aux(i) = MAXVAL(Cohesive(el_a)%Criterion)

			

		END DO
		
		! New coefficient to modify the time step				
		eta_load = MAXVAL(eta_load_aux) 
		DEALLOCATE(eta_load_aux)

		IF (Iteration.EQ.1) THEN
			WRITE (*,*) ''
			WRITE (*,*) ''
			WRITE (*,'(A,I5,A)') '     -- Propagation of :', cont, ' interfaces'
			WRITE (*,*) ''
		ELSE

			WRITE (*,*) ''
			WRITE (*,*) ''
			WRITE (*,'(A,I5,A)') '     -- No propagation occurs'	
			WRITE (*,*) ''
		END IF	

		reset  = 0
		DO i=1,ninterfaces
			
		    reg_a = Interfaces(i,2)
		    reg_b = Interfaces(i,3)    
		
		    el_a = Interfaces(i,4)
		    el_b = Interfaces(i,5)

			iTSL = Interfaces(i,8)

			ndci_a = ELEM(el_a,8)-1
			ndci_b = ELEM(el_b,8)-1

			! If at least one node of the interfaces fails in one or more axes
			! It is required to evaluate the new sub-blocks 
			! H_cohe, G_cohe
			IF (Interfaces(i,18).EQ.1) THEN

				reset = reset + 1
				IF (reset.EQ.1) THEN
					Reduce_A = 0 
				END IF

				ind_i_a = ELEM(el_a,6)-1
				ind_i_b = ELEM(el_b,6)-1
				n_dof_a = El_reg(reg_a,1)*nnos_el*dof
				n_dof_b = El_reg(reg_b,1)*nnos_el*dof

				IF (reg_a.LT.reg_b) THEN 

					ndr_a = El_reg(reg_a,1)*nnos_el*dof ! rows block F of reg_A
					ndc_a = nnos_el*dof ! columns block F of reg_A
					ndr_b = El_reg(reg_b,1)*nnos_el*dof + El_reg(reg_b,8) ! rows block F of reg_B
					ndc_b = nnos_el*dof ! columns block F of reg_B
					
				ELSE

					ndr_a = El_reg(reg_a,1)*nnos_el*dof + El_reg(reg_a,8) ! rows block G of reg_B
					ndc_a = nnos_el*dof + Interfaces(i,20) ! columns block G of reg_A
					ndr_b = El_reg(reg_b,1)*nnos_el*dof  ! rows block G of reg_B
					ndc_b = nnos_el*dof + Interfaces(i,20) ! columns block G of reg_B
					
				END IF
				
				IF ((Cohesive(el_a)%alloc.EQ.1).AND.(Cohesive(el_b)%alloc.EQ.1)) THEN
					
					DEALLOCATE(Cohesive(el_a)%H_cohe)
					DEALLOCATE(Cohesive(el_b)%H_cohe)
					DEALLOCATE(Cohesive(el_a)%G_cohe)
					DEALLOCATE(Cohesive(el_b)%G_cohe)
			
					IF (Transient.EQ.1) THEN
						DEALLOCATE(Cohesive(el_a)%M_cohe)
						DEALLOCATE(Cohesive(el_b)%M_cohe)
						DEALLOCATE(Cohesive(el_b)%M_cohe_G)
					END IF

					DEALLOCATE(Cohesive(el_a)%H_local)
					DEALLOCATE(Cohesive(el_a)%G_local)
					DEALLOCATE(Cohesive(el_b)%H_local)
					DEALLOCATE(Cohesive(el_b)%G_local)

					IF (Transient.EQ.1) THEN
						DEALLOCATE(Cohesive(el_a)%M_local)
						DEALLOCATE(Cohesive(el_b)%M_local)
					END IF
			
					Cohesive(el_a)%alloc = 0
					Cohesive(el_b)%alloc = 0
		
				END IF

				IF ((Cohesive(el_a)%alloc.EQ.0).AND.(Cohesive(el_b)%alloc.EQ.0)) THEN

					IF (reg_a.LT.reg_b) THEN

					! === From global to local coordinates ===============
						
						ALLOCATE(H_aux1(n_dof_a,dof),G_aux1(n_dof_a,dof))
						ALLOCATE(H_aux2(n_dof_b,dof),G_aux2(n_dof_b,dof))
						ALLOCATE(Cohesive(el_a)%H_local(n_dof_a,nnos_el*dof))
						ALLOCATE(Cohesive(el_a)%G_local(n_dof_a,nnos_el*dof))
						ALLOCATE(Cohesive(el_b)%H_local(n_dof_b,nnos_el*dof))
						ALLOCATE(Cohesive(el_b)%G_local(n_dof_b,nnos_el*dof))
	
						IF (Transient.EQ.1) THEN
							ALLOCATE(M_aux1(n_dof_a,dof),M_aux2(n_dof_b,dof))
							ALLOCATE(Cohesive(el_a)%M_local(n_dof_a,nnos_el*dof))
							ALLOCATE(Cohesive(el_b)%M_local(n_dof_b,nnos_el*dof))
						END IF

						DO j=1,nnos_el

							col_i = ind_i_a + (3*j-2) 
							col_f = ind_i_a + 3*j
						
							H_aux1 = H_G(reg_a)%H(:,col_i:col_f)
							G_aux1 = H_G(reg_a)%G(:,col_i:col_f)
							IF (Transient.EQ.1) THEN
								M_aux1 = U_T_E_M(reg_a)%M(:,col_i:col_f)
							END IF
							 
							Cohesive(el_a)%H_local(:,3*j-2:3*j) = MATMUL(H_aux1,TSL(iTSL)%R)
							Cohesive(el_a)%G_local(:,3*j-2:3*j) = MATMUL(G_aux1,TSL(iTSL)%R)
							IF (Transient.EQ.1) THEN
								Cohesive(el_a)%M_local(:,3*j-2:3*j) = MATMUL(M_aux1,TSL(iTSL)%R)
							END IF


							col_i = ind_i_b + (3*j-2) 
							col_f = ind_i_b + 3*j

							H_aux2 = H_G(reg_b)%H(:,col_i:col_f)
							G_aux2 = H_G(reg_b)%G(:,col_i:col_f)
							IF (Transient.EQ.1) THEN
								M_aux2 = U_T_E_M(reg_b)%M(:,col_i:col_f)
							END IF
	
							Cohesive(el_b)%H_local(:,3*j-2:3*j) = MATMUL(H_aux2,TSL(iTSL)%R)
							Cohesive(el_b)%G_local(:,3*j-2:3*j) = MATMUL(G_aux2,TSL(iTSL)%R)
							IF (Transient.EQ.1) THEN
								Cohesive(el_b)%M_local(:,3*j-2:3*j) = MATMUL(M_aux2,TSL(iTSL)%R)
							END IF

						END DO

						DEALLOCATE(H_aux1,G_aux1,H_aux2,G_aux2)
						IF (Transient.EQ.1) THEN
							DEALLOCATE(M_aux1,M_aux2)
						END IF

						! ====================================================
						ALLOCATE(Cohesive(el_a)%H_cohe(ndr_a,ndc_a))
						ALLOCATE(Cohesive(el_b)%H_cohe(ndr_b,ndc_b))

						IF (Transient.EQ.1) THEN
							ALLOCATE(Cohesive(el_a)%M_cohe(ndr_a,ndc_a))
							ALLOCATE(Cohesive(el_b)%M_cohe(ndr_b,ndc_b))
							ALLOCATE(Cohesive(el_b)%M_cohe_G(ndr_b,ndc_b))
							Cohesive(el_a)%M_cohe = 0.d0
							Cohesive(el_b)%M_cohe = 0.d0
							Cohesive(el_b)%M_cohe_G = 0.d0
						END IF
						

						Cohesive(el_a)%H_cohe = 0.d0
						Cohesive(el_b)%H_cohe = 0.d0						

					ELSE
						ALLOCATE(Cohesive(el_a)%G_cohe(ndr_a,ndc_a))
						ALLOCATE(Cohesive(el_b)%G_cohe(ndr_b,ndc_b))
						Cohesive(el_a)%G_cohe = 0.d0
						Cohesive(el_b)%G_cohe = 0.d0
						Cohesive(el_a)%alloc = 1
						Cohesive(el_b)%alloc = 1
					END IF
					
				END IF

				DO j=1,nnos_el*dof
					
					! - Pristine interface
					! - loading Cohesive
					! - Separation					
				
					condition = Interfaces(i,8+j)
					
					IF (condition.EQ.0) THEN ! Pristine interface

						IF (reg_a.LT.reg_b) THEN

							Cohesive(el_a)%H_cohe(:,j) = Cohesive(el_a)%H_local(:,j) 

							Cohesive(el_b)%H_cohe(1:n_dof_b,j) = Cohesive(el_b)%H_local(:,j)

							IF (Transient.EQ.1) THEN

								Cohesive(el_a)%M_cohe(:,j) = Cohesive(el_a)%M_local(:,j)

								Cohesive(el_b)%M_cohe(1:n_dof_b,j) = Cohesive(el_b)%M_local(:,j)

							END IF
  						
						ELSE

							Cohesive(el_a)%G_cohe(1:n_dof_a,j) = -Cohesive(el_a)%G_local(:,j)

							Cohesive(el_b)%G_cohe(:,j) = Cohesive(el_b)%G_local(:,j)

							!IF (Transient.EQ.1) THEN

							!	Cohesive(el_a)%M_cohe(1:n_dof_a,j) = 0.d0

							!	Cohesive(el_b)%M_cohe(:,j) = 0.d0

							!END IF
							
						END IF

					ELSEIF (condition.EQ.1) THEN ! Separation
						
						IF (reg_a.LT.reg_b) THEN		

							Cohesive(el_a)%H_cohe(:,j) = Cohesive(el_a)%H_local(:,j) 

							Cohesive(el_b)%H_cohe(1:n_dof_b,j) = 0.d0
	
							Reduce_A = Reduce_A + n_dof_b

							IF (Transient.EQ.1) THEN

								Cohesive(el_a)%M_cohe(:,j) = Cohesive(el_a)%M_local(:,j)
								
								Cohesive(el_b)%M_cohe(1:n_dof_b,j) = 0.d0
	
								!Cohesive(el_b)%M_cohe_G(1:n_dof_b,j) = Cohesive(el_b)%M_local(:,j)

							END IF


						ELSE

							Cohesive(el_a)%G_cohe(1:n_dof_a,j) = Cohesive(el_a)%H_local(:,j) 

							Cohesive(el_b)%G_cohe(:,j) = 0.d0
	
							Reduce_A = Reduce_A + n_dof_b

							IF (Transient.EQ.1) THEN

								Cohesive(el_a)%M_cohe_G(1:n_dof_a,j) = Cohesive(el_a)%M_local(:,j)

								!Cohesive(el_b)%M_cohe(:,j) = 0.d0

							END IF

						END IF
													
					ELSEIF (condition.EQ.2) THEN ! Contact
						
						IF (reg_a.LT.reg_b) THEN

							Cohesive(el_a)%H_cohe(:,j) = Cohesive(el_a)%H_local(:,j) 

							Cohesive(el_b)%H_cohe(1:n_dof_b,j) = Cohesive(el_b)%H_local(:,j)

							IF (Transient.EQ.1) THEN

								Cohesive(el_a)%M_cohe(:,j) = Cohesive(el_a)%M_local(:,j)

								Cohesive(el_b)%M_cohe(1:n_dof_b,j) = Cohesive(el_b)%M_local(:,j)

								u_t(ndci_b+j,:) = 0.d0

							END IF
							
						ELSE

							Cohesive(el_a)%G_cohe(1:n_dof_a,j) = -Cohesive(el_a)%G_local(:,j)

							Cohesive(el_b)%G_cohe(:,j) = Cohesive(el_b)%G_local(:,j)

							IF (Transient.EQ.1) THEN

								Cohesive(el_a)%M_cohe(1:n_dof_a,j) = 0.d0

								Cohesive(el_b)%M_cohe(:,j) = 0.d0

							END IF
							
						END IF

					END IF

				END DO

			END IF

		END DO

	END IF

	IF (me.EQ.0) THEN
		
		Damage_elements = 0
		
		DO i=1,nelem
			IF (Damage_zones(i).NE.0) THEN

				Damage_elements = Damage_elements + 1
				
			END IF
		END DO 
		
	END IF

	

    t2 = MPI_Wtime() 
    duration = t2-t1
    CALL MPI_REDUCE(duration,time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,mpierr);
    IF (me.EQ.0) THEN
        WRITE (*,'(A,F15.3,A)') '                                  ...COMPLETED! (Time :',time,'s)'
        WRITE (*,*) ''
		WRITE (*,*) '----------------------------------------------------------------------------'
    END IF
    CALL MPI_Bcast(time,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
    !---------------------------------------------------------------------------------

END SUBROUTINE Local_fields
!===============================================================================!
END MODULE Failure_Analysis
