!===============================================================================!
!-------------------------------------------------------------------------------!
!          MODULE TO COMPUTE THE VECTOR {A} AND VECTOR {b}                      !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Cohesive_process
!-------------------------------------------------------------------------------!
USE Global_variables
USE Set_parameters
USE General_Sparce_Vectors
USE Response_D_T
USE Module_Solver
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE Iterative_process(me,nt,step,time)
	
	IMPLICIT NONE
    INCLUDE 'mpif.h'

	INTEGER :: me, step, mpierr, root=0, nt
	REAL(8) :: ts6, ts8, ts7, ts9

    REAL(8) :: t1, t2, duration, time
    
    t1 = MPI_Wtime();

	CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr);

1	IF (me.EQ.0) THEN

		CALL New_cohesive_matrices
	
		!Applying boundary conditions to obtain {A} and {b} SPARCE
		CALL Final_Sparce_Vectors(ts6,step,nt)
		
		!6. New vector b
        CALL b_update(ts7,step)

	END IF

	!SOLUTION [A]x = {b}
    CALL Solver_MUMPS(me,nt,ts8)

	IF (me.EQ.0) THEN

		! Displacement and traction response
        CALL Displacement_traction(ts9,step)

	END IF

	CALL MPI_Bcast(Iteration,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
	
	IF (Iteration.EQ.1) THEN
    	GOTO 1 ! Next iteration in the same step
    ELSE
		
        GOTO 2 ! Next step
		
    END IF

2	IF (me.EQ.0) THEN
		WRITE (*,*) '----------------------------------------------------------------------------'
		WRITE(*,*) ''
		WRITE (*,'(A,I3,A)') '     -- Next step load'
		WRITE(*,*) ''
		WRITE (*,*) '----------------------------------------------------------------------------'
		WRITE (*,*) '----------------------------------------------------------------------------'
	END IF

	!---------------------------------------------------------------------------------
    t2 = MPI_Wtime() 
    duration = t2-t1
    CALL MPI_REDUCE(duration,time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,mpierr);
    CALL MPI_Bcast(time,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
    !---------------------------------------------------------------------------------
	

END SUBROUTINE Iterative_process
!===============================================================================!
SUBROUTINE New_cohesive_matrices

	INTEGER :: i, j, ind_g(9), ninterfaces, reg_a, reg_b, el_a, el_b, dof
	INTEGER :: condition, iTSL, ind_i_a, ind_i_b, n_dof_a, n_dof_b,ind_g2(9)
	INTEGER :: ndr_a, ndc_a, ndr_b, ndc_b, k1, ind ,ndci_a, ndci_b
	REAL(8), ALLOCATABLE :: H_aux(:), G_aux(:)
	REAL(8) :: damage

	ind_g = (/1,2,3,1,2,3,1,2,3/)
	ind_g2 = (/1,1,1,2,2,2,3,3,3/)	
	ninterfaces = SIZE(Interfaces,1)

	dof = 3
	
	WRITE (*,*) '----------------------------------------------------------------------------'
	WRITE (*,*) '----------------------------------------------------------------------------'
	WRITE (*,*) '----------------------------------------------------------------------------'
	WRITE (*,*) ' '
	WRITE (*,'(A,I3,A)') '     -- Update block matrices'
	WRITE (*,*) ' '

	El_reg(:,7:8) = 0
	Interfaces(:,19:20) = 0
	ELEM(:,14) = nnos_el*dof 
		
	DO i=1,ninterfaces

		reg_a = Interfaces(i,2)
		reg_b = Interfaces(i,3)    
		
		el_a = Interfaces(i,4)
		el_b = Interfaces(i,5)

		iTSL = Interfaces(i,8)

		! If at least one node of the interfaces fails in one or more axes
		! It is required to evaluate the new sub-blocks 
		! H_cohe, G_cohe
		IF (Interfaces(i,18).EQ.1) THEN
	
			ind_i_a = ELEM(el_a,6)-1
			ind_i_b = ELEM(el_b,6)-1
			n_dof_a = El_reg(reg_a,2)*nnos_el*dof
			n_dof_b = El_reg(reg_b,2)*nnos_el*dof

			IF (reg_a.LT.reg_b) THEN 

				ndr_a = El_reg(reg_a,2)*nnos_el*dof ! rows block F of reg_A
				ndc_a = nnos_el*dof ! columns block F of reg_A
				ndr_b = El_reg(reg_b,2)*nnos_el*dof + El_reg(reg_b,8) ! rows block F of reg_B
				ndc_b = nnos_el*dof ! columns block F of reg_B
					
			ELSE

				ndr_a = El_reg(reg_a,2)*nnos_el*dof + El_reg(reg_a,8) ! rows block G of reg_B
				ndc_a = nnos_el*dof + Interfaces(i,20) ! columns block G of reg_A
				ndr_b = El_reg(reg_b,2)*nnos_el*dof  ! rows block G of reg_B
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
				END IF
			
				Cohesive(el_a)%alloc = 0
				Cohesive(el_b)%alloc = 0
		
			END IF
			
			IF ((Cohesive(el_a)%alloc.EQ.0).AND.(Cohesive(el_b)%alloc.EQ.0)) THEN

				IF (reg_a.LT.reg_b) THEN
				
					ALLOCATE(Cohesive(el_a)%H_cohe(ndr_a,ndc_a))
					ALLOCATE(Cohesive(el_b)%H_cohe(ndr_b,ndc_b))
				
					IF (Transient.EQ.1) THEN
						ALLOCATE(Cohesive(el_a)%M_cohe(ndr_a,ndc_a))
						ALLOCATE(Cohesive(el_b)%M_cohe(ndr_b,ndc_b))
						Cohesive(el_a)%M_cohe = 0.d0
						Cohesive(el_b)%M_cohe = 0.d0
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

			k1 = 0
			ind = Interfaces(i,19)-1	

			ndci_a = ELEM(el_a,8)-1
			ndci_b = ELEM(el_b,8)-1
			
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

						IF (Transient.EQ.1) THEN

							Cohesive(el_a)%M_cohe(1:n_dof_a,j) = 0.d0

							Cohesive(el_b)%M_cohe(:,j) = 0.d0

						END IF
							
					END IF

				ELSEIF (condition.EQ.1) THEN ! Separation
						
					IF (reg_a.LT.reg_b) THEN		

						Cohesive(el_a)%H_cohe(:,j) = Cohesive(el_a)%H_local(:,j) 

						Cohesive(el_b)%H_cohe(1:n_dof_b,j) = 0.d0

						IF (Transient.EQ.1) THEN

							Cohesive(el_a)%M_cohe(:,j) = Cohesive(el_a)%M_local(:,j)
								
							Cohesive(el_b)%M_cohe(1:n_dof_b,j) = 0.d0

						END IF


					ELSE

						Cohesive(el_a)%G_cohe(1:n_dof_a,j) = Cohesive(el_a)%H_local(:,j) 

						Cohesive(el_b)%G_cohe(:,j) = 0.d0

						IF (Transient.EQ.1) THEN

							Cohesive(el_a)%M_cohe(1:n_dof_a,j) = Cohesive(el_a)%M_local(:,j)

							Cohesive(el_b)%M_cohe(:,j) = 0.d0

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

END SUBROUTINE New_cohesive_matrices
!===============================================================================!
END MODULE Cohesive_process





















