!====================================================================================!
!====================================================================================!
!                                                                                    !
!               U N I C A M P - UNIVERSIDADE ESTADUAL DE CAMPINAS                    !
!                     F E M - FACULDADE DE ENGENHARIA MECÂNICA                       !
!                  D M C - DEPARTAMENTE DE MECÂNICA COMPUTACIONAL                    !
!                                                                                    !
! AUTHOR: Andrés Felipe Galvis Rodriguez                                             !
! SUPERVISOR: Prof. Dr. Paulo Sollero                                                !
!************************************************************************************!
! Program: 3D_BEM_Aniso                                                              !
! Description: 3D-BEM PROGRAM TO SOLVE LINEAR ELASTICITY PROBLEMS                    !
!               Using the fundamental solution based on double Fourier series        !
!               and 3-node triangular Discontinuous elements                         !
!====================================================================================!
!====================================================================================!
PROGRAM Main

!====================================================================================!
! MODULES
!------------------------------------------------------------------------------------!

    USE Global_variables
	USE Global_functions
    USE Set_parameters
    USE Input
    USE Discretization
    USE BC_Module
    USE Int_points
    USE Comp_H_G
    USE General_Matrices
    USE General_Sparce_Vectors
    USE Vectors_b
    USE Module_Solver
    USE Response_D_T
    USE Response_S_S
    USE Homogenization    
    USE Failure_Analysis
    USE Output
    
    !USE mpi
    !USE omp_lib

    IMPLICIT NONE
    include 'mpif.h'
    REAL(8) :: t1, t2, duration, global
    INTEGER :: me,nt,mpierr

    REAL(8) :: ts0, ts1, ts2, ts3, ts4, ts5, ts6, ts7, ts8, ts9, ts10, ts11
    REAL(8) :: ts12a, ts12b, ts13, ts14, tts
    REAL(8) :: t=0
    INTEGER :: r=1, Steps=1, step, root=0

    CALL MPI_INIT(mpierr) 
    !---------------------------------------------------------------------------------
    t1 = MPI_Wtime();
    !---------------------------------------------------------------------------------
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nt,mpierr) 
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,me,mpierr)   
    CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr);

    IF (me.EQ.0) THEN 
        CALL start_greet        
    END IF

	CALL Check_Size_nt(nt)

	CALL Setup

    !0. Matrices of integration points
    CALL Integration_Points(me,ts0)

    IF (me.EQ.0) THEN
        
        !1. Reading data from the file
        CALL Read_data(ts1,r)			

    END IF
	    
    ! 2. Matrices of subregions and interfaces
    CALL Subregions_Interfaces(nt,me,ts2)
    
    IF (me.EQ.0) THEN

        !3. Sub-regions, Interfacesm, Nodes and Boundary Conditions
        CALL Discretization_Matrices(ts3)

		IF (Failure.EQ.1) THEN
            
            ! Reading the traction-separation law
            CALL Read_TSL

        END IF 

		Steps = 1 ! Default
        IF (Transient.EQ.1) THEN
            t = 0.d0
			T_t = Dt*Time_steps
            Steps = Time_steps
        END IF

        IF (Quasi_static.EQ.1)THEN
            Steps = Static_steps
        END IF

        !4. Matrices of Boundary Conditions and Non_interfaces
        CALL BC_Ninterfaces_matrices(ts4,Steps)

		! check if the problems is adequate in size
		CALL Check_Size(nt)

    END IF

    !5. Compute matrices [H] and [G]
    CALL Matrices_H_G(nt,me,ts5)

	IF (me.EQ.0) THEN

		tts = tts + ts0 + ts1 + ts2 + ts3 + ts4 + ts5 + ts6

    END IF

    CALL MPI_Bcast(Steps,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)

    DO step=1,Steps

        IF (me.EQ.0) THEN

			! Only for the first step if there is not failure
			! If failure occurs, the system must be assembly
			
			IF ((Iteration.EQ.1).OR.(step.EQ.1)) THEN 
				
				IF (nreg.EQ.1) THEN   
				    !6. Applying boundary conditions to obtain [A] and {b}
				    CALL Final_Arrays_1D(ts6,nt)  
				ELSE
				    !6. Applying boundary conditions to obtain {A} and {b} SPARCE
				    CALL Final_Sparce_Vectors(ts6,step,nt)
				ENDIF				

			END IF

            !6. New vector b
            CALL b_update(ts7,step)
            
        END IF

        !7. SOLUTION [A]x = {b}
        CALL Solver_MUMPS(me,nt,ts8)
		
        IF (me.EQ.0) THEN
			
            !8. Displacement and traction response
            CALL Displacement_traction(ts9,step)

		END IF

		CALL MPI_Bcast(Iteration,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)		
		ts14 = 0.d0
		
		IF (me.EQ.0) THEN

            !9. Strain and Stress tensors
            CALL Strain_Stress(ts10,step)

            IF (H_fields.EQ.1) THEN

                !10. Compute the homogenized sigma and epsilon
                CALL Homogenized_sigma_epsilon(ts11,step)

            END IF

		END IF

		IF (Failure.EQ.1) THEN

            CALL Local_fields(nt,me,ts13,step)

        END IF
		
		IF (me.EQ.0) THEN
			
            !9. Output data 
	    	
			CALL Output_data_par(ts12b,step) ! paraview
			

        END IF

		IF (me.EQ.0) THEN

	    ! CURRENT processing time
            tts = tts + ts7 + ts8 + ts9 + ts10 + ts11 + ts12a + ts12b + ts13 + ts14

            IF (Transient.EQ.1) THEN
                t = t + Dt
                WRITE (*,*) ''
                WRITE (*,*) '****************************************************************************'
                WRITE (*,'(A,I10)') ' ******* STEP : ', step
                WRITE (*,'(A,F15.2,A)') ' ******* CURRENT PROCESSING TIME : ', tts, ' s'
                WRITE (*,*) '****************************************************************************'
                WRITE (*,*) ''
            
	    ELSEIF (Quasi_static.EQ.1) THEN

		WRITE (*,*) ''
                WRITE (*,*) '****************************************************************************'
                WRITE (*,'(A,I10)') ' ******* STEP : ', step  
                WRITE (*,'(A,F15.2,A)') ' ******* CURRENT PROCESSING TIME : ', tts, ' s'
                WRITE (*,*) '****************************************************************************'
                WRITE (*,*) ''

	    END IF  
		
	END IF

    END DO

    !---------------------------------------------------------------------------------
    t2 = MPI_Wtime() 
    duration = t2-t1
    CALL MPI_REDUCE(duration,global,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,mpierr);
    IF (me.EQ.0) THEN
        WRITE (*,*) '****************************************************************************'
        WRITE (*,*) ''
        WRITE(*,'(A,F15.4,A)') '              TOTAL PROCESSING TIME: ',global,' s'
        WRITE (*,*) ''
        WRITE (*,*) '*****************************************************************************' 
    END IF
    !---------------------------------------------------------------------------------

    CALL MPI_FINALIZE(mpierr)

END PROGRAM Main




