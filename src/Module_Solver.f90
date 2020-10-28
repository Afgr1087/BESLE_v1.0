!
!  This file is part of MUMPS 5.0.2, released
!  on Fri Jul 15 09:12:54 UTC 2016
!
      MODULE Module_Solver

      USE Global_variables 
      USE Global_functions 

      CONTAINS

      SUBROUTINE Solver_MUMPS(me,nt,time)

      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC) mumps_par
      INTEGER IERR, I
      !REAL(8), ALLOCATABLE :: b(:,:), A_values(:)
      !INTEGER, ALLOCATABLE :: rows_A(:), columns_A(:)

      INTEGER :: me, nt, sub_nt, j, root=0, status(MPI_STATUS_SIZE)
      REAL(8) :: t1, t2, duration, time
    
      t1 = MPI_Wtime();
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR);

      IF (me.EQ.0) THEN 
        
        WRITE (*,*) '----------------------------------------------------------------------------'
        WRITE (*,*) ' '
        WRITE (*,'(A)',advance='no') ' 8. Solver MUMPS'
        WRITE (*,*) ' '
      END IF

! Define a communicator for the package.
      mumps_par%COMM = MPI_COMM_WORLD

      IF (me.EQ.0) THEN 
           WRITE (*,*) ' '
      	   WRITE (*,*) '     -- Define a communicator for the package.'
      END IF
!  Initialize an instance of the package
!  for L U factorization (sym = 0, with working host)
      
      IF (me.EQ.0) THEN 
      	   WRITE (*,*) '     -- Initialize an instance of the package'
      	   WRITE (*,*) '     -- for L U factorization (sym = 0, with no working host).'
      END IF

      mumps_par%WRITE_PROBLEM = 'A_values'

      mumps_par%JOB = -1
      mumps_par%SYM = 0
      mumps_par%PAR = 0
        
      CALL DMUMPS(mumps_par)
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
                  "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),& 
                  "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
       GOTO 500
      END IF

    
      
!  Define problem on the host (processor 0)

      IF (me.EQ.0) THEN 
      	   WRITE (*,*) '     -- Define problem on all processors except the host.'
		   WRITE (*,*) '        Size: N = ', N_total, 'NZ = ', nzA
      END IF

      IF ( mumps_par%MYID .eq. 0 ) THEN 
        
        mumps_par%ICNTL(14) = 50
        mumps_par%ICNTL(22) = Out_of_core
		mumps_par%ICNTL(23) = 1024*Working_memory
        mumps_par%ICNTL(4) = 1
		mumps_par%ICNTL(5) = 0	
		mumps_par%ICNTL(18) = 3
        
        mumps_par%N = N_total
        mumps_par%NNZ = nzA
		
		ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
		mumps_par%RHS = b(:,1)

		
	  END IF
	 
	  sub_nt = nt-1
      
	  IF (mumps_par%MYID .EQ. 0 ) THEN
			
	  		DO j=1,sub_nt

				CALL MPI_Send(scounts_A(j),1,MPI_INTEGER,j,100,MPI_COMM_WORLD,IERR)

		    END DO
			
	  ELSE 

	    	CALL MPI_Recv(mumps_par%NNZ_loc,1,MPI_INTEGER,0,100,MPI_COMM_WORLD,status,IERR )

	  END IF

	  IF ( mumps_par%MYID .GT. 0 ) THEN
    	ALLOCATE(mumps_par%IRN_loc(mumps_par%NNZ_loc) )
		ALLOCATE(mumps_par%JCN_loc(mumps_par%NNZ_loc) )
		ALLOCATE(mumps_par%A_loc(mumps_par%NNZ_loc) )	
	  END IF

	   IF (mumps_par%MYID .EQ. 0 ) THEN
			
	  		DO j=1,sub_nt

				CALL MPI_Send(S_E(j)%row_A,scounts_A(j),MPI_INTEGER,j,101,MPI_COMM_WORLD,IERR)
				CALL MPI_Send(S_E(j)%col_A,scounts_A(j),MPI_INTEGER,j,102,MPI_COMM_WORLD,IERR)
				CALL MPI_Send(S_E(j)%At,scounts_A(j),MPI_REAL8,j,103,MPI_COMM_WORLD,IERR)

		    END DO
			
	  ELSE 

			CALL MPI_Recv(mumps_par%IRN_loc,mumps_par%NNZ_loc,MPI_INTEGER,0,101,MPI_COMM_WORLD,status,IERR )
			CALL MPI_Recv(mumps_par%JCN_loc,mumps_par%NNZ_loc,MPI_INTEGER,0,102,MPI_COMM_WORLD,status,IERR )
			CALL MPI_Recv(mumps_par%A_loc,mumps_par%NNZ_loc,MPI_REAL8,0,103,MPI_COMM_WORLD,status,IERR )

	  END IF
		
      mumps_par%OOC_TMPDIR = 'OOC'

!  Call package for solution
		
      IF (me.EQ.0) THEN 
      	   WRITE (*,*) '     -- Call package for solution.'
      END IF

      mumps_par%JOB = 6
	  
      CALL DMUMPS(mumps_par)
	
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
                  "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),& 
                  "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
       GOTO 500
      END IF
!  Solution has been assembled on the host

      IF (me.EQ.0) THEN 
      	   WRITE (*,*) '     -- Solution has been assembled on the host.'
      END IF

      IF ( mumps_par%MYID .eq. 0 ) THEN

        ALLOCATE( x_BEM ( mumps_par%N,1  ) )    
        x_BEM(:,1) = mumps_par%RHS
      END IF
!  Deallocate user data

      IF (me.EQ.0) THEN 
      	   WRITE (*,*) '     -- Deallocate user data.'
      END IF

	  IF ( mumps_par%MYID .GT. 0 ) THEN
	  	DEALLOCATE( mumps_par%IRN_loc )
      	DEALLOCATE( mumps_par%JCN_loc )
      	DEALLOCATE( mumps_par%A_loc   )
	  END IF
      IF ( mumps_par%MYID .eq. 0 )THEN
        DEALLOCATE( mumps_par%RHS )
      END IF
!  Destroy the instance (deallocate internal data structures)

      IF (me.EQ.0) THEN 
      	   WRITE (*,*) '     -- Destroy the instance (deallocate internal data structures).'
      END IF
      mumps_par%JOB = -2
      CALL DMUMPS(mumps_par)
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
                  "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),& 
                  "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
       GOTO 500
      END IF
    
      !---------------------------------------------------------------------------------
500   t2 = MPI_Wtime() 
      duration = t2-t1
      CALL MPI_REDUCE(duration,time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,IERR);
      IF (me.EQ.0) THEN
	    WRITE (*,*) ' '
        WRITE(*,'(A,F15.3,A)') '                                  ...COMPLETED! (Time :',time,'s)'
        WRITE (*,*) ''
      END IF
      CALL MPI_Bcast(time,1,MPI_INTEGER,root,MPI_COMM_WORLD,IERR)
      !---------------------------------------------------------------------------------

      END SUBROUTINE Solver_MUMPS
      END MODULE Module_Solver
