!===============================================================================!
!-------------------------------------------------------------------------------!
!          MODULE TO COMPUTE THE MATRIX [A] AND VECTOR {b}                      !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE General_Matrices
!-------------------------------------------------------------------------------!
USE Global_variables
USE Vectors_b
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE Final_Arrays_1D(ts1,nt)

    INTEGER::t0,t1,rate, nt, i 
	INTEGER(8) :: nt_sub
    REAL(8),INTENT(OUT)::ts1

    INTEGER :: reg, nelem, dof
    REAL(8), ALLOCATABLE :: H_reg(:,:), G_reg(:,:), A(:,:)

    INTEGER :: n
    INTEGER, ALLOCATABLE :: rowind(:), colind(:)!, ja_aux2(:), ia_aux2(:)
    REAL(8), ALLOCATABLE :: acoo(:)!, acsr_aux2(:)

	INTEGER(8), ALLOCATABLE :: displs_A(:)

    CALL SYSTEM_CLOCK(t0, rate) ! Initializates the time counter
    WRITE (*,*) '-----------------------------------------------------------'
    WRITE (*,*) ' '
    WRITE (*,'(A)',advance='no') ' 6. Vectors {A} and {b}           ...'

        reg = nreg
        nelem = El_reg(reg,1) ! Number of elements
        dof = 3 ! degrees of fredom
        ALLOCATE(A(nelem*nnos_el*dof,nelem*nnos_el*dof))
        ALLOCATE(H_reg(nelem*nnos_el*dof,nelem*nnos_el*dof))
        ALLOCATE(G_reg(nelem*nnos_el*dof,nelem*nnos_el*dof))
        ALLOCATE(G_geral(nelem*nnos_el*dof,nelem*nnos_el*dof))
        
        H_reg = H_G(reg)%H; G_reg = H_G(reg)%G

	
        ! Matrix [V] Body forces analysis
        IF (Bodyforces.EQ.1) THEN

			IF (Transient.EQ.0) THEN

            	ALLOCATE(V(nelem*nnos_el*dof,nelem*nnos_el*dof))
            	
				V = U_T_E_M(reg)%M

			END IF

			ALLOCATE(bhat(nelem*nnos_el*dof,1))
                        
            CALL Bodyforces_vet

            

        END IF

        ! Matrix [V] Trasient analysis
        IF (Transient.EQ.1) THEN
            ALLOCATE(u_t(nelem*nnos_el*dof,3))
            u_t = 0.d0
            ALLOCATE(V(nelem*nnos_el*dof,nelem*nnos_el*dof))
            V = U_T_E_M(reg)%M
        END IF

        CALL Apply_BC(H_reg,G_reg,reg,A)  

        CALL MakeSparse(A,rowind, colind, acoo,n) 

	!-----------------------------------------------------------------------

		nt_sub = nt - 1 ! Threads without the master "0"
	
		ALLOCATE(S_E(nt_sub),scounts_A(nt_sub),displs_A(nt_sub)) 

		N_total = nelem*nnos_el*dof ! rows or columns  of general system

		nzA = n

		! Function for parallel distribution
		! Par_Div(threads, vector_size, scounts)
		CALL Par_Div(nt_sub, nzA, scounts_A, displs_A) ! A_values		

		DO i=1,nt_sub-1

			ALLOCATE(S_E(i)%At(scounts_A(i)))
			S_E(i)%At = acoo(displs_A(i)+1:displs_A(i+1))
			
			ALLOCATE(S_E(i)%col_A(scounts_A(i)))
			S_E(i)%col_A = colind(displs_A(i)+1:displs_A(i+1))

			ALLOCATE(S_E(i)%row_A(scounts_A(i)))
			S_E(i)%row_A = rowind(displs_A(i)+1:displs_A(i+1))

		END DO

		ALLOCATE(S_E(nt_sub)%At(scounts_A(nt_sub)))
		S_E(i)%At = acoo(displs_A(nt_sub)+1:displs_A(nt_sub)+scounts_A(nt_sub))

		ALLOCATE(S_E(nt_sub)%col_A(scounts_A(nt_sub)))
		S_E(i)%col_A = colind(displs_A(nt_sub)+1:displs_A(nt_sub)+scounts_A(nt_sub))

		ALLOCATE(S_E(nt_sub)%row_A(scounts_A(nt_sub)))
		S_E(i)%row_A = rowind(displs_A(nt_sub)+1:displs_A(nt_sub)+scounts_A(nt_sub))

        
		DEALLOCATE(acoo,colind,rowind,displs_A)


        CALL SYSTEM_CLOCK(t1)       ! Time counter stopped
        WRITE (*,'(A,F15.3,A)') 'COMPLETED! (Time :', REAL(t1-t0)/rate,'s)'
        WRITE (*,*) ''
        ts1 = REAL(t1-t0)/rate

END SUBROUTINE Final_Arrays_1D
!===============================================================================!
SUBROUTINE MakeSparse(A,rowind,colind,acoo,cont)

REAL(8), ALLOCATABLE::acoo(:)
INTEGER, ALLOCATABLE::rowind(:),colind(:)

INTEGER :: ii, jj, cont, mm, nn
        
REAL(8), ALLOCATABLE :: val_aux(:), A(:,:)       
INTEGER, ALLOCATABLE :: i_row_aux(:),j_col_aux(:)

mm = SIZE(A,1); nn = SIZE(A,2)

cont = 0

ALLOCATE( val_aux(mm*nn), i_row_aux(mm*nn), j_col_aux(mm*nn) )
val_aux = 0.d0; i_row_aux = 0; j_col_aux = 0

DO ii = 1,mm
    DO jj = 1,nn
        IF (ABS(A(ii,jj)) .GE. 1E-20) THEN 
            cont = cont + 1
            val_aux(cont) = A(ii,jj)
            i_row_aux(cont) = ii
            j_col_aux(cont) = jj   
        END IF
    END DO
END DO

ALLOCATE( acoo(cont), rowind(cont), colind(cont) )
acoo = 0.d0; rowind = 0; colind = 0

acoo = val_aux(1:cont)
rowind = i_row_aux(1:cont)
colind = j_col_aux(1:cont)

DEALLOCATE(val_aux,i_row_aux,j_col_aux,A)
    
END SUBROUTINE MakeSparse
!===============================================================================!
SUBROUTINE MakeCSR(nrow, nnz, rowind, colind, acoo, rowcsr, colcsr, acsr)

INTEGER::nrow, nnz, rowind(:), colind(:)
REAL(8)::acoo(:)

INTEGER,ALLOCATABLE::rowcsr(:), colcsr(:)
REAL(8),ALLOCATABLE::acsr(:)

INTEGER::kk, jj, ii, iad, k0
REAL(8)::xx

ALLOCATE( rowcsr(nrow + 1) )
rowcsr = 0

! Row lengths
DO kk = 1, nnz
    rowcsr(rowind(kk)) = rowcsr(rowind(kk)) + 1
END DO

! Starting point of each row
kk = 1
DO jj = 1, nrow + 1
    k0 = rowcsr(jj)
    rowcsr(jj) = kk
    kk = kk + k0
END DO

!  Go through the structure once more.  Fill in output matrix.
ALLOCATE( acsr(nnz), colcsr(nnz) )
DO kk = 1, nnz
    ii = rowind(kk)
    jj = colind(kk)
    xx = acoo(kk)
    iad = rowcsr(ii)
    acsr(iad) = xx
    colcsr(iad) = jj
    rowcsr(ii) = iad + 1
END DO

!  Shift back rowcsr
DO jj = nrow, 1, -1
    rowcsr(jj + 1) = rowcsr(jj)
END DO
rowcsr(1) = 1

END SUBROUTINE MakeCSR
!===============================================================================!
SUBROUTINE Apply_BC(Hg,Gg,reg,A)
    
    INTEGER :: nelem, dof, m, i, j, ind, reg, type_bc
    REAL(8) :: val_d, val_t
    REAL(8),ALLOCATABLE :: change(:,:), Hg(:,:), Gg(:,:), A(:,:)    

    nelem = El_reg(reg,1) ! Number of elements
    dof = 3 ! degrees of fredom
    m = nnos_el*dof
    ALLOCATE(change(nelem*m,1),b_cte(nelem*m,1),B_geral(nelem*m,1)); 
    ALLOCATE(b(nelem*m,1), b_step(nelem*m,1), b_bf(nelem*m,1));
    change = 0.d0; b = 0.d0; b_cte = 0.d0; b_step = 0.d0; b_bf = 0.d0
    !WRITE(*,*) ''
    DO i=1,nelem
        DO j=1,m ! select the type of BC
            type_bc = INT(BC_elem(i,2*j)) ! type of BC
            !WRITE(*,*) i,j,type_bc,BC_elem(i,2*j+1)
            IF (type_bc .EQ. 0) THEN ! Displacement is known
                ind = i*m-m+j ! index column to be change in [H]
                change(:,1) = Gg(:,ind)
                Gg(:,ind) = -Hg(:,ind) 
                Hg(:,ind)  = -change(:,1)
                val_d = BC_elem(i,2*j+1) ! Displacement value
                B_geral(ind,1) = val_d  
            ELSE ! Traction is known
               ind = i*m-m+j
               val_t = BC_elem(i,2*j+1) ! Traction value 
               B_geral(ind,1) = val_t
            END IF
        END DO
    END DO

    A = Hg
    G_geral = Gg
    
    DEALLOCATE(Hg,Gg,change)

END SUBROUTINE Apply_BC
!===============================================================================!
END MODULE General_Matrices

