!===============================================================================!
!-------------------------------------------------------------------------------!
!                               MODULE FOR INPUT DATA                           !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Input
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!
CONTAINS
    SUBROUTINE start_greet
        WRITE(*,*) ''
        WRITE(*,*) '*****************************************************************************'
        WRITE(*,*) '*****************************************************************************'
        WRITE(*,*) '*****************************************************************************'
        WRITE(*,*) '*****************************************************************************'
        WRITE(*,*) '                                  BESLE:                               '
	    WRITE(*,*) '               Boundary Element Software for 3D Linear Elasticiy '
		WRITE(*,*) '        An MPI-parallelized Fortran package for heterogeneous materials                    '
        WRITE(*,*) ''
		WRITE(*,*) 'By Andrés F. Galvis'
		WRITE(*,*) 'Contributions: Daniel M. Prada'
		WRITE(*,*) '               Lucas S. Moura'
		WRITE(*,*) 'Coordinator: Paulo Sollero'
		WRITE(*,*) 'School of Mechanical Engineering'
		WRITE(*,*) 'University of Campinas'
		WRITE(*,*) 'Brazil'
		WRITE(*,*) '28/10/2020'
        WRITE(*,*) '*****************************************************************************'
        WRITE(*,*) '*****************************************************************************'
        WRITE(*,*) '*****************************************************************************'
        WRITE(*,*) '*****************************************************************************'
    END SUBROUTINE start_greet
!===============================================================================!
    SUBROUTINE Read_data(ts1,r)

        INTEGER :: i, m, n, num, r, dof=3
        INTEGER::t0,t1,rate
        REAL(8),INTENT(OUT)::ts1
        CHARACTER(LEN=20)::num_str
        REAL(8) :: max_x, max_y, max_z, min_x, min_y, min_z

        CALL SYSTEM_CLOCK(t0, rate) ! Initializates the time counter
        WRITE (*,*) '----------------------------------------------------------------------------'
        WRITE (*,*) ' '
        WRITE (*,'(A)',advance='no') ' 1. Reading data from *.dat file   ...'

		OPEN(1,file=trim(fileplace_mesh)//trim(Mesh_file)//'.dat',STATUS='OLD')

        REWIND(1)

		READ (1,*) nreg

!-----------------------------------------------------------------------
! POINTS : [Node, node1, node2, node3]
!-----------------------------------------------------------------------
! Read matrix of points

        READ (1,*) m
        ALLOCATE(POINTS(m,3))
        DO i = 1,m
            READ (1,*) POINTS(i,1:3)
        END DO
        
        POINTS = POINTS*scale_size_1

        max_x = MAXVAL(POINTS(:,1))
        max_y = MAXVAL(POINTS(:,2))
        max_z = MAXVAL(POINTS(:,3))

        min_x = MINVAL(POINTS(:,1))
        min_y = MINVAL(POINTS(:,2))
        min_z = MINVAL(POINTS(:,3))

        POINTS(:,1) = POINTS(:,1) - min_x
        POINTS(:,2) = POINTS(:,2) - min_y
        POINTS(:,3) = POINTS(:,3) - min_z

        max_x = MAXVAL(POINTS(:,1))
        max_y = MAXVAL(POINTS(:,2))
        max_z = MAXVAL(POINTS(:,3))

! -----------------------------------------------------------------------
! ELEM : [Element, Face, node1, node2, node3]
! -----------------------------------------------------------------------
! Read matrix f elements
        READ (1,*) m
        ALLOCATE(ELEM(m,14))
        ELEM = 0
		ELEM(:,14) = nnos_el*dof
        DO i = 1,m
            READ (1,*) ELEM(i,1:3)
        END DO
! -----------------------------------------------------------------------
! NORMAL_VECTORS : [Face, nx, ny, nz]
! -----------------------------------------------------------------------
! Read matrix of normal vectors in each face
        READ (1,*) m
        
        ALLOCATE(NORMAL_VECTORS(m,3))
        ALLOCATE(NORMAL_VEC_Def(m,3))
        NORMAL_VEC_Def = 0.d0
        DO i = 1,m
            READ (1,*) NORMAL_VECTORS(i,1:3)
        END DO
        
! -----------------------------------------------------------------------
! El_reg : [Region, Elements]
! -----------------------------------------------------------------------
! Read matrix of elements per region
        READ (1,*) m 
        
        ALLOCATE(El_reg(m,8))
        El_reg = 0
	
        DO i = 1,m
            READ (1,*) El_reg(i,1)
            
        END DO

! -----------------------------------------------------------------------
    
        CLOSE(1)

        IF (n_Materials.GE.1) THEN
            CALL materials_list(r)
        END IF

        CALL SYSTEM_CLOCK(t1)       ! Time counter stopped
        WRITE (*,'(A,F15.2,A)') 'COMPLETED! (Time :', REAL(t1-t0)/rate,'s)'
        WRITE (*,*) ''
        ts1 = REAL(t1-t0)/rate
        
    END SUBROUTINE Read_data
!===============================================================================!

    SUBROUTINE release_mem

        DEALLOCATE(FACES,POINTS,ELEM,NORMAL_VECTORS,El_reg,Volumes)
        DEALLOCATE(GEO_PHY_ELEM,Subregions)

        DEALLOCATE(GEO_PHY_NODES,GEO_NODES_global,PHY_NODES_global)
        DEALLOCATE(ELEM_GEO_NODES_global)

        DEALLOCATE(bcd_0,bcd_1,bcd_2)
        DEALLOCATE(bcsd_3)
        DEALLOCATE(bct_0,bct_1,bct_2)
        DEALLOCATE(bcs_3)
        DEALLOCATE(BC_elem, bc_new)
        
        DEALLOCATE(Angles)

        IF (nreg.GT.1) THEN
            
            DEALLOCATE(Non_interfaces,Non_int_plot)

            DEALLOCATE(Int_reg,Interfaces)

            IF (Transient.EQ.1) THEN
            !    DEALLOCATE(rows_V,V_values,columns_V)
                DEALLOCATE(b_step,u_t, Disp_p)
            END IF

            IF (Bodyforces.EQ.1) THEN
                !DEALLOCATE(rows_V,bhat,V_values,columns_V)
		DEALLOCATE(rows_V,bhat)
            END IF
            IF (Quasi_static.EQ.1) THEN
                DEALLOCATE(Disp_p)
            END IF
            !DEALLOCATE(A_values,columns_A, rows_A)
            !DEALLOCATE(G_values,columns_G, rows_G)

        ELSE

            !DEALLOCATE(acsr_aux,ja_aux,ia_aux)

            !DEALLOCATE(A,G_geral)
            !DEALLOCATE(G_geral)
            IF (Transient.EQ.1) THEN
                DEALLOCATE(V,b_step,u_t,Disp_p)
            END IF
            IF (Bodyforces.EQ.1) THEN
                DEALLOCATE(V,bhat)
            END IF
            IF (Quasi_static.EQ.1) THEN
                DEALLOCATE(Disp_p)
            END IF

        END IF

        DEALLOCATE(b,b_cte)
        DEALLOCATE(B_geral, Displacement, Traction)
           
    END SUBROUTINE release_mem
!===============================================================================!
SUBROUTINE materials_list(r)

    INTEGER :: i,num, r, val
    CHARACTER(LEN=30)::num_str, Test_name

    num = r
    WRITE(num_str,"(I10)") num
  
    !Test_name = "Model"//TRIM(ADJUSTL(num_str))//"/MaterialsM"//TRIM(ADJUSTL(num_str))//".txt"
    Test_name = "List_Mat.dat"

    OPEN(1,file=trim(fileplace_material)//trim(Material_coefficients_database)//'/'//trim(Test_name),STATUS='OLD')

        ALLOCATE(Mat_reg(nreg,2))
	    Mat_reg = 0	
        DO i = 1,nreg
			Mat_reg(i,1) = i
            READ (1,*) Mat_reg(i,2)
            
        END DO

    CLOSE(1)

! -----------------------------------------------------------------------

!===============================================================================!
END SUBROUTINE materials_list

END MODULE Input
