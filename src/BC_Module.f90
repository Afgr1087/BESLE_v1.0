!===============================================================================!
!-------------------------------------------------------------------------------!
!                      MODULE TO SPECIFY THE BOUNDARY CONDITIONS               !
!
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE BC_Module
!-------------------------------------------------------------------------------!
USE Global_variables
USE Global_functions
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE BC_Ninterfaces_matrices(ts1,Steps)

    INTEGER::t0,t1,rate,Steps
    REAL(8),INTENT(OUT)::ts1
    
    INTEGER :: npf, nfaces, i, j, pos, INFO, INFO2, ind1, n_el, nelem, el, face
    INTEGER ::  principal_face, cont1, ind, face_bc, ind2, ind3, ind4, k
    INTEGER :: cont, cond, vec_face(6), val, BC_info(6,2)
    REAL(8) :: val_max_x, val_max_y, val_max_z, bc(21), bc_aux(6), Tol
	REAL(8) :: vec_box_bc(6), vec_aux(6)
    !REAL(8), ALLOCATABLE :: BC_elem_aux(:,:)
    INTEGER, ALLOCATABLE :: bc_new_aux(:,:)

    INTEGER :: size_v
    INTEGER, ALLOCATABLE :: Vector(:)

    INTEGER :: n, m
	CHARACTER(LEN=20)::num_str
	CHARACTER(LEN=100)::f_name
    
    CALL SYSTEM_CLOCK(t0, rate) ! Initializates the time counter
    WRITE (*,*) '----------------------------------------------------------------------------'
    WRITE (*,*) ' '
    WRITE (*,'(A)',advance='no') ' 4. Boundary conditions            ...'

    IF (BC_groups.EQ.0) THEN

		vec_face = 0
		vec_box_bc = 0.d0
		vec_box_bc(1) = box_face_1(2)
		vec_box_bc(2) = box_face_2(2)
		vec_box_bc(3) = box_face_3(2)
		vec_box_bc(4) = box_face_4(2)
		vec_box_bc(5) = box_face_5(2)
		vec_box_bc(6) = box_face_6(2)
	
		cont1 = 0
		DO i=1,6
			IF (vec_box_bc(i).GT.-0.5d0) THEN
				cont1 = cont1 + 1
				vec_face(cont1) = i
			END IF
		END DO
		
		ALLOCATE(T_D_BC(cont1*7))
		ind2=0
		DO i=1,6

			val = vec_face(i)


			IF (val.GT.0) THEN
	
				ind1 = ind2 + 1
				ind2 = ind2 + 7

				SELECT CASE (val)		
					CASE (1)				
						T_D_BC(ind1:ind2) = (/1.d0,box_face_1/)
					CASE (2)				
						T_D_BC(ind1:ind2) = (/2.d0,box_face_2/)
					CASE (3)				
						T_D_BC(ind1:ind2) = (/3.d0,box_face_3/)
					CASE (4)				
						T_D_BC(ind1:ind2) = (/4.d0,box_face_4/)
					CASE (5)				
						T_D_BC(ind1:ind2) = (/5.d0,box_face_5/)
					CASE (6)			
						T_D_BC(ind1:ind2) = (/6.d0,box_face_6/)
				END SELECT
			END IF
				
		END DO
		
		
		CALL bct_vector(Steps)
		BC_info = 0
		
		DO i=1,cont1

			ind1 = 7*i-6
			ind2 = 7*i

			vec_aux = T_D_BC(ind1+1:ind2)
			val = T_D_BC(ind1)
		    
			cond = 0
			DO j=1,3

				cond = cond + INT(vec_aux(2*j))

			END DO

			IF (cond.EQ.0) THEN ! displacement imposed

				BC_info(val,:) = (/val, cond/) ! face and type of condition
				
			ENDIF
			

		END DO

	END IF
	
    !POINTS = POINTS*scale_size
    npf = 6 ! Number of principal faces
    val_max_x = MAXVAL(POINTS(:,1))
    val_max_y = MAXVAL(POINTS(:,2))
    val_max_z = MAXVAL(POINTS(:,3))
    Tol = (5e-5) !* dmin 
    
    IF (nreg.GT.1) THEN ! For more than one region  
        !-------------Non-Interfaces----------------------
        
        nfaces = SIZE(ELEM,1) ! number of faces
        ALLOCATE(Non_interfaces(nfaces,5))
        Non_interfaces = 0

        
        
        ind4 = 0
        
        size_v = SIZE(Interfaces,1)
        ALLOCATE(Vector(size_v))
        Vector = Interfaces(:,6)     
        
                
        ind2 = 0
        
        DO i=1,nfaces 
                
            n_el = 1 ! Elements per face (each element is a face)
            ! Searching if the ith face is an interface
            !CALL find_vector(Interfaces(:,6),j,pos)

            CALL find_vector_g(Vector,size_v,i,pos) 

            ind1 = ind2 + 1   
            ind2 = ind2 + n_el

            IF (pos.EQ.0) THEN ! The face is not an interface

                DO j=1,npf ! Principal face
        
                    ! Find the principal face to each face
                    CALL Principal_faces(j,ind2,val_max_x,val_max_y,val_max_z,INFO,INFO2)

                    IF (INFO.EQ.1) THEN ! The face belongs to principal face
                        ind4 = ind4 + 1
                        Non_interfaces(i,:) = (/i,j,n_el,ind1,ind2/) ! Face
                        EXIT
                    END IF 

                    IF((INFO.EQ.0).AND.(INFO2.EQ.0).AND.(j.EQ.npf)) THEN
                        
                        ind4 = ind4 + 1
                        Non_interfaces(i,:) = (/i,npf+1,n_el,ind1,ind2/) ! Face
                        
                    END IF

                END DO

            END IF
                
        END DO
        
        DEALLOCATE(Vector) 

        ALLOCATE(Non_int_plot(ind4,5))
        Non_int_plot = 0
        ind3 = 0
        DO i=1,nfaces
            IF (Non_interfaces(i,1).NE.0) THEN
                ind3 = ind3 + 1
                Non_int_plot(ind3,:) = Non_interfaces(i,:)
            END IF
        END DO

        !---------------------------------------------
        ! ------------ Bc_elem -----------------------

        IF (BC_groups.EQ.0) THEN

            nelem = SIZE(ELEM,1)
            ALLOCATE(BC_elem(nelem,21))
            BC_elem = 0.d0
            cont = 0

			ALLOCATE(Nodes_BCs(nelem*nnos_el,2))
			Nodes_BCs = 0
            

            DO el=1,nelem

                face = ELEM(el,4)
                principal_face = Non_interfaces(face,2)
                
                IF (principal_face.NE.0) THEN

                    IF (principal_face.GT.npf) THEN

                        bc_aux = (/1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0/)
                        bc(1:3) = (/el,face,principal_face/)
                        bc(4:21) = (/bc_aux,bc_aux,bc_aux/)

                        BC_elem(el,:) = bc

                        !if (principal_face.EQ.7) then
            
                        !    WRITE(*,*) principal_face,el
                            
                        !end if

                    ELSE

                        CALL bc_vector(principal_face,face,el,bc)    
           
                        BC_elem(el,:) = bc

						IF (BC_info(principal_face,1).NE.0) THEN
			
							Nodes_BCs(3*el-2,:) = (/3*el-2,el/)
							Nodes_BCs(3*el-1,:) = (/3*el-1,el/)
							Nodes_BCs(3*el,:) = (/3*el,el/)

						END IF
                       
                       !if (principal_face.EQ.5) then
            
                       !     WRITE(*,*) el
                            
                       !end if

                    END IF
                   
                END IF

            END DO

        ELSE 

            nelem = SIZE(ELEM,1)

			ALLOCATE(BCs(BC_groups))

			DO i = 1,BC_groups
			    
				WRITE(num_str,"(I10)") i
				f_name = "BCs_"//TRIM(ADJUSTL(num_str))//".dat"
		
				OPEN(1,file=trim(fileplace_BCs)//trim(f_name),STATUS='OLD')
				
					READ (1,'(1I10)') BCs(i)%nel  ! Number of elements
					
					m = BCs(i)%nel
					ALLOCATE(BCs(i)%el(m))
					DO j = 1,m
						READ (1,'(1I10)') BCs(i)%el(j)
					END DO	
					
					READ (1,'(3I10)') BCs(i)%type_bc ! Type of BCs 
					
					ALLOCATE(BCs(i)%bc(Steps,3))

					DO j = 1,Steps

						READ (1,'(3F27.15)') BCs(i)%bc(j,:)
						
					END DO

				CLOSE(1)

			END DO
			
			nelem = SIZE(ELEM,1)
            ALLOCATE(BC_elem(nelem,21))
            BC_elem=0.d0

			ALLOCATE(Nodes_BCs(nelem*nnos_el,2))
			Nodes_BCs = 0

			DO i=1,BC_groups

				DO j=1,BCs(i)%nel

					el = BCs(i)%el(j)
					face = ELEM(el,4)  

					CALL bct_vector_M(i,el,bc_aux)

					bc(1:3) = (/el,face,4/)
                    
                    bc(4:21) = (/bc_aux,bc_aux,bc_aux/)        
                    
                    BC_elem(el,:) = bc

					Nodes_BCs(3*el-2,:) = (/3*el-2,el/)
					Nodes_BCs(3*el-1,:) = (/3*el-1,el/)
					Nodes_BCs(3*el,:) = (/3*el,el/)

				END DO
			END DO

			DO i=1,nelem

            	el = INT(BC_elem(i,1))
					
               	IF (el .EQ. 0) THEN

                	face = ELEM(i,2)
            
                    bc_aux = (/1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0/)
                    bc(1:3) = (/i,face,4/)
                    bc(4:21) = (/bc_aux,bc_aux,bc_aux/)

                    BC_elem(i,:) = bc       
						
                END IF


        	END DO

        END IF

    ELSE      

		!---------------------------------------------
        ! ------------ Bc_elem -----------------------

        IF (BC_groups.EQ.0) THEN 
        
		    nelem = SIZE(ELEM,1)
		    ALLOCATE(BC_elem(nelem,19))

			ALLOCATE(Nodes_BCs(nelem*nnos_el,2))
			Nodes_BCs = 0
		    
		    INFO2 = 0
		    DO i=1,npf ! Principal face

		        DO j=1,nelem
		            el = j 
		            ! Find the principal face to each face
		            
		            CALL Principal_faces(i,el,val_max_x,val_max_y,val_max_z,INFO,INFO2)
		            
		            IF (INFO2.NE.0) THEN

		                principal_face = INFO2
		                CALL bc_vector(principal_face,face,j,bc)
						BC_elem(j,1) = principal_face
		                BC_elem(j,2:19) = bc(2:19)
						
						Nodes_BCs(3*j-2,:) = (/3*j-2,el/)
						Nodes_BCs(3*j-1,:) = (/3*j-1,el/)
						Nodes_BCs(3*j,:) = (/3*j,el/)

		            END IF 
		            
		        END DO

		    END DO

		ELSE

			nelem = SIZE(ELEM,1)

			ALLOCATE(BCs(BC_groups))

			DO i = 1,BC_groups
			    
				WRITE(num_str,"(I10)") i
				f_name = "BCs_"//TRIM(ADJUSTL(num_str))//".dat"
		
				OPEN(1,file=trim(fileplace_BCs)//trim(f_name),STATUS='OLD')
				
				READ (1,*) BCs(i)%nel  ! Number of elements
					
				m = BCs(i)%nel
				ALLOCATE(BCs(i)%el(m))
				DO j = 1,m
					READ (1,*) BCs(i)%el(j)
				END DO	
					
				READ (1,*) BCs(i)%type_bc ! Type of BCs   
					
				ALLOCATE(BCs(i)%bc(Steps,3))

				DO j = 1,Steps

					READ (1,*) BCs(i)%bc(j,:)
						
				END DO

				CLOSE(1)
	
			END DO

			nelem = SIZE(ELEM,1)
            ALLOCATE(BC_elem(nelem,21))
            BC_elem=0.d0

			ALLOCATE(Nodes_BCs(nelem*nnos_el,2))
			Nodes_BCs = 0

			DO i=1,BC_groups

				DO j=1,BCs(i)%nel

					el = BCs(i)%el(j)
					face = ELEM(el,4)  

					CALL bct_vector_M(i,el,bc_aux)

					bc(1) = el
                    
                    bc(2:19) = (/bc_aux,bc_aux,bc_aux/)        
                    
                    BC_elem(el,:) = bc

					Nodes_BCs(3*el-2,:) = (/3*el-2,el/)
					Nodes_BCs(3*el-1,:) = (/3*el-1,el/)
					Nodes_BCs(3*el,:) = (/3*el,el/)

				END DO
			END DO

			DO i=1,nelem

            	el = INT(BC_elem(i,1))
					
               	IF (el .EQ. 0) THEN
            
                    bc_aux = (/1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0/)
                    bc(1) = i
                    bc(2:19) = (/bc_aux,bc_aux,bc_aux/)

                    BC_elem(i,:) = bc       
						
                END IF


        	END DO


		END IF     

    END IF

	!===========================================================================
	CALL nodes_conectivity_geo
	!===========================================================================

	
    
	!===========================================================================

    ! Traction and displacement ---------------------------------------
    
    IF (BC_groups.EQ.0) THEN

        ALLOCATE(bc_new_aux(nelem,2) )
        
        bc_new_aux = 0

        cont1 = 0
        
        DO i=1,SIZE(T_D_BC)/7
            
            ind = (i-1)*7
            face_bc = INT(T_D_BC(ind+1)) 
            
            DO k=1,nelem
                IF (nreg.GT.1) THEN
                    face = INT(BC_elem(k,3))
                ELSE
                    face = INT(BC_elem(k,1))!ELEM(k,4)
                END IF  
                IF (face_bc.EQ.face) THEN  
                    cont1 = cont1 + 1         
                    bc_new_aux(cont1,:) = (/face_bc,k/)
                END IF
            END DO
      
        END DO                 
        
        ALLOCATE(bc_new(cont1,2))
        bc_new = bc_new_aux(1:cont1,:)
        
        DEALLOCATE(bc_new_aux) 

    END IF

    ! If there are constrained nodes in the model ---------------------------
    cond = Config_BC(1) ! Boundary condition
    IF (cond.GT.0) THEN

        ! Find nodes to be constrained
        CALL Find_constrained_node

        ! Apply the boundary conditions to those nodes
        CALL Constrain_node

    END IF
    !-----------------------------------------------------------------------

    !POINTS = POINTS/scale_size

    CALL SYSTEM_CLOCK(t1)       ! Time counter stopped
    WRITE (*,'(A,F15.3,A)') 'COMPLETED! (Time :', REAL(t1-t0)/rate,'s)'
    WRITE (*,*) ''
    ts1 = REAL(t1-t0)/rate

END SUBROUTINE BC_Ninterfaces_matrices
!===============================================================================!
SUBROUTINE Principal_faces(principal_face,el,val_max_x,val_max_y,val_max_z,INFO,INFO2)

    INTEGER :: el,INFO,no1,no2,no3, principal_face, INFO2
    REAL(8) :: val_max_x,val_max_y,val_max_z
    REAL(8) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, Tol
    REAL(8) :: Er11x, Er21x, Er31x, Er12x, Er22x, Er32x
    REAL(8) :: Er11y, Er21y, Er31y, Er12y, Er22y, Er32y
    REAL(8) :: Er11z, Er21z, Er31z, Er12z, Er22z, Er32z
   
    ! Nodes
    no1 = ELEM(el,1); no2 = ELEM(el,2); no3 = ELEM(el,3);
    
    ! Coordinates nodes
    x1=POINTS(no1,1); y1=POINTS(no1,2); z1=POINTS(no1,3);
    x2=POINTS(no2,1); y2=POINTS(no2,2); z2=POINTS(no2,3);
    x3=POINTS(no3,1); y3=POINTS(no3,2); z3=POINTS(no3,3);
   
    Tol = (5e-5) !* dmin 

    INFO = 0 
    INFO2 = 0 
    SELECT CASE(principal_face)
        CASE(4)
            ! Error, when x=0 
            Er11x=ABS(x1-0.d0);      Er21x=ABS(x2-0.d0);      Er31x=ABS(x3-0.d0);
            ! Coordinate x of nodes 1,2 and 3 are equal to x = 0
            IF ((Er11x.LT.Tol).AND.(Er21x.LT.Tol).AND.(Er31x.LT.Tol)) THEN
                INFO = 1
                INFO2 = 4
            END IF
        CASE(2)
            ! Error, when x=val_max_x
            Er12x=ABS(x1-val_max_x); Er22x=ABS(x2-val_max_x); Er32x=ABS(x3-val_max_x);
            ! Coordinate x of nodes 1,2 and 3 are equal to x = val_max_x
            IF ((Er12x.LT.Tol).AND.(Er22x.LT.Tol).AND.(Er32x.LT.Tol)) THEN
                INFO = 1
                INFO2 = 2
            END IF
        CASE(1)
            ! Error, when y=0 
            Er11y=ABS(y1-0.d0);      Er21y=ABS(y2-0.d0);      Er31y=ABS(y3-0.d0);
            ! Coordinate y of nodes 1,2 and 3 are equal to y = 0
            IF ((Er11y.LT.Tol).AND.(Er21y.LT.Tol).AND.(Er31y.LT.Tol)) THEN
                INFO = 1
                INFO2 = 1
            END IF
        CASE(3)
            ! Error, when y=val_max_y
            Er12y=ABS(y1-val_max_y); Er22y=ABS(y2-val_max_y); Er32y=ABS(y3-val_max_y);
            ! Coordinate y of nodes 1,2 and 3 are equal to y = val_max_x
            IF ((Er12y.LT.Tol).AND.(Er22y.LT.Tol).AND.(Er32y.LT.Tol)) THEN
                INFO = 1
                INFO2 = 3
            END IF
        CASE(5)
            ! Error, when z=0 
            Er11z=ABS(z1-0.d0);      Er21z=ABS(z2-0.d0);      Er31z=ABS(z3-0.d0);
            ! Coordinate y of nodes 1,2 and 3 are equal to z = 0
            IF ((Er11z.LT.Tol).AND.(Er21z.LT.Tol).AND.(Er31z.LT.Tol)) THEN
                INFO = 1
                INFO2 = 5
            END IF
        CASE(6)
            ! Error, when z=val_max_z
            Er12z=ABS(z1-val_max_z); Er22z=ABS(z2-val_max_z); Er32z=ABS(z3-val_max_z);
            ! Coordinate z of nodes 1,2 and 3 are equal to z = val_max_z
            IF ((Er12z.LT.Tol).AND.(Er22z.LT.Tol).AND.(Er32z.LT.Tol)) THEN
                INFO = 1
                INFO2 = 6
            END IF 
    END SELECT

END SUBROUTINE Principal_faces
!===============================================================================!
SUBROUTINE bc_vector(principal_face,face,el,bc) 

    INTEGER :: face, el,principal_face
    REAL(8) :: bc(21), bc_aux(6)

    IF (nreg.GT.1) THEN    
    
        bc_aux = BC_Face(principal_face,:)
        bc(1:3) = (/el,face,principal_face/)
        bc(4:21) = (/bc_aux,bc_aux,bc_aux/)

    ELSE

        bc_aux = BC_Face(principal_face,:)
        bc(1) = el
        bc(2:19) = (/bc_aux,bc_aux,bc_aux/)

    END IF

END SUBROUTINE bc_vector
!===============================================================================!
SUBROUTINE bct_vector(Steps)

    INTEGER :: i, face, ind, j, Type_bc, k, cont, Steps, cont1
    REAL(8) :: CTE, t, Applied_bc

    ! Type_BC_(x,y,z) --> 0 Displacement or 1 Traction
    ! BC_Face = [Type_BC_x, Val_BC_x, Type_BC_y, Val_BC_y, Type_BC_y, Val_BC_y] 
    BC_Face(1,:) = (/1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0/)
    BC_Face(2,:) = (/1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0/)
    BC_Face(3,:) = (/1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0/)
    BC_Face(4,:) = (/1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0/)
    BC_Face(5,:) = (/1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0/)
    BC_Face(6,:) = (/1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0/)

    ! displacement and traction ---------------------------------------------
	
	ALLOCATE(BCs(SIZE(T_D_BC)/7))
	

    DO i=1,SIZE(T_D_BC)/7 ! face
    
        ind = (i-1)*7
        face = INT(T_D_BC(ind+1)) 
        cont=2
        cont1=1

		ALLOCATE(BCs(i)%bc(3,Steps))
		BCs(i)%bc = 0.d0
 
        DO j=1,nnos_el! x,y,z

            Type_bc = INT(T_D_BC(ind+2*j+1)) ! Displacement or traction
		            
			Applied_bc = T_D_BC(ind+2*j)

            ! Traction
            IF (Type_bc.EQ.1) THEN               
                
                BC_Face(face,cont1)= 1.d0                
     
            ELSE ! Displacement

                BC_Face(face,cont1)= 0.d0                        

            END IF

			IF (load_profile.EQ.'ramp') THEN

				CTE = Applied_bc/Steps
		        t = 0.d0
		        DO k=1,Steps
		        	t = t + 1
					BCs(i)%bc(j,k) = CTE*k
		        END DO
				BC_Face(face,cont) = BCs(i)%bc(j,1)
		        cont=cont+2
		        cont1=cont1+2

			ELSEIF (load_profile.EQ.'Heaviside') THEN

				CTE = Applied_bc

		        DO k=1,Steps
					BCs(i)%bc(j,k) = CTE!*k
		        END DO
				BC_Face(face,cont) = BCs(i)%bc(j,1)
		        cont=cont+2
		        cont1=cont1+2

			ELSEIF (load_profile.EQ.'harmonic') THEN

				CTE = Applied_bc
				t = 0.d0;
		        DO k=1,Steps
					t = t + dt;
					BCs(i)%bc(j,k) = CTE*DSIN( omega*t + phase)
		        END DO
				BC_Face(face,cont) = BCs(i)%bc(j,1)
		        cont=cont+2
		        cont1=cont1+2

			END IF

        END DO
    END DO

END SUBROUTINE bct_vector

!===============================================================================!
SUBROUTINE bct_vector_M(group,el,bc)

    INTEGER :: i, nelem, j, m, Type_bc, group, el, bc_xyz(3)

    REAL(8) :: bc(6), bc_val_xyz(3), normal(3)

	bc = 0.d0

	IF (BCs(group)%type_bc(1).EQ.2) THEN ! normal displacement

		bc_xyz = (/0,0,0/)

		normal = NORMAL_VECTORS(el,:)

		bc_val_xyz = normal*BCs(group)%bc(1,1)

	ELSE IF (BCs(group)%type_bc(1).EQ.3) THEN ! normal traction

		bc_xyz = (/1,1,1/)

		normal = NORMAL_VECTORS(el,:)

		bc_val_xyz = normal*BCs(group)%bc(1,1)
		
	ELSE IF ( (BCs(group)%type_bc(1).EQ.0).OR. &
				(BCs(group)%type_bc(1).EQ.1) ) THEN

	   bc_xyz = BCs(group)%type_bc ! in Cartesian coordinates
	   bc_val_xyz = BCs(group)%bc(1,:)
	

	END IF


    ! Traction and displacement
          
        DO j=1,nnos_el

            Type_bc = bc_xyz(j)
                 	   
            ! Traction
            IF (Type_bc.EQ.1) THEN
                
                bc(2*j-1) = 1.d0
    
            ELSE ! Displacements
                
                bc(2*j-1) = 0.d0
				
            END IF

			bc(2*j) = bc_val_xyz(j)
			
        END DO
	 
END SUBROUTINE bct_vector_M



!===============================================================================!
SUBROUTINE Find_constrained_node

    INTEGER :: Condition_BC(6), FACES_point(6,4), i, nodes, face, cond, face_e
    INTEGER :: nelem, j, k, kk, no, cont, pos, el, face_el
    REAL(8) :: BC_nodes(6,3), Tol, x, y, z, x_el, y_el, z_el, d, val
    REAL(8) :: xmax, ymax, zmax, xmin, ymin, zmin


    xmax = MAXVAL(POINTS(:,1)); ymax = MAXVAL(POINTS(:,2)); zmax = MAXVAL(POINTS(:,3));
    xmin = MINVAL(POINTS(:,1)); ymin = MINVAL(POINTS(:,2)); zmin = MINVAL(POINTS(:,3));

    Tol = 0.8d0 * dmin  

    ! Configuration of BC [node_a,node_b,node_c,node_d]
    Condition_BC = (/1, 0, 0, 0, 0, 0/)

    ! Faces of each point [face_a,face_b,face_c,face_d]
    FACES_point(1,:) = (/5, 6, 0, 0/)
    FACES_point(2,:) = (/0, 0, 0, 0/)
    FACES_point(3,:) = (/0, 0, 0, 0/)
    FACES_point(4,:) = (/0, 0, 0, 0/)
    FACES_point(5,:) = (/0, 0, 0, 0/)
    FACES_point(6,:) = (/0, 0, 0, 0/)    

    ! Points to be found [x,y,z]   
    BC_nodes(1,:) = (/(xmax-xmin)/2,  ymin, (zmax-zmin)/2/)
    BC_nodes(2,:) = (/xmax,  (ymax-ymin)/2, (zmax-zmin)/2/)
    BC_nodes(3,:) = (/(xmax-xmin)/2 , ymax, (zmax-zmin)/2/)
    BC_nodes(4,:) = (/xmin,  (ymax-ymin)/2, (zmax-zmin)/2/)
    BC_nodes(5,:) = (/xmin, (ymax)/2,  0.25d0*zmax/)
    BC_nodes(6,:) = (/xmin, (ymax)/2,  0.75d0*zmax/)

    nelem = SIZE(ELEM_GEO_NODES_global,1)
    
    contrain_nodes = 0
    
    ! Find the nearest nodes
    WRITE(*,*) ''
    DO i=1,1!SIZE(Config_BC)
        
        cond = 1 !Config_BC(i)
        nodes = Condition_BC(i)
        
        
        DO j=1,nodes
            cont = 1
            val = 0
            ! Face of specific node
            face = FACES_point(cond,j)
            
            ! Coordinates of specific node
            x = BC_nodes(face,1); y = BC_nodes(face,2); z = BC_nodes(face,3); 
            
            ! Elements that belong to the face
            DO k=1,nelem     
                IF (nreg.GT.1) THEN
                    face_el = INT(BC_elem(k,3))
                ELSE
                    face_el = ELEM(k,4)
                END IF 
				face=4
				
                IF (face.EQ.face_el) THEN
                    
                    DO kk=1,3 ! Loop for three nodes of each element
                        ! Coordinate of the kth node in the jth element
                        no = ELEM_GEO_NODES_global(k,kk+1) 
                        x_el = PHY_NODES_global(no,2);
                        y_el = PHY_NODES_global(no,3);         
                        z_el = PHY_NODES_global(no,4);  
                        ! Distance between specific node and kth node in the jth element
                        d = DSQRT((x-x_el)**2+(y-y_el)**2+(z-z_el)**2) 
                     
                        ! Determining the minimun distance
                        IF (cont.EQ.1) THEN
                            val = d  
                        END IF
                        IF ((cont.GT.1).AND.(val.GT.d)) THEN
                            val = d
                            pos = kk ! position in the element 1,2 or 3
                            el = k ! element
                        END IF                          
                        cont = cont + 1
                    END DO
                END IF 
            END DO
            contrain_nodes(j,:) = (/el,pos/)
        END DO 
    END DO

	
    
END SUBROUTINE Find_constrained_node
!===============================================================================!
SUBROUTINE Constrain_node

    INTEGER :: cond, i, el, pos
    REAL(8) :: bc(6)

    cond = 1!Config_BC(1) ! Boundary condition
    
    ! apply the boundary condition in the BC_elem matrix

    SELECT CASE (cond)
       
        CASE(1) 
            bc = (/1.d0,0.d0,1.d0,0.d0,0.d0,0.d0/)       
                    
        !CASE(2) 
        !    bc = (/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0/) 
            
        !CASE(3) 
        !    bc = (/0.d0,0.d0,1.d0,0.d0,1.d0,0.d0/) 

        !CASE(4)
        !    bc = (/1.d0,0.d0,1.d0,0.d0,0.d0,0.d0/) 

        !CASE(5)
        !    bc = (/0.d0,0.d0,1.d0,0.d0,1.d0,0.d0/) 

        !CASE(6)
        !    bc = (/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0/) 
        
    END SELECT

    DO i=1,SIZE(contrain_nodes,1)
        el = contrain_nodes(i,1)
        pos = contrain_nodes(i,2)

        WRITE(*,*) el, pos

        IF (el.NE.0) THEN
            IF (nreg.GT.1) THEN
                
                SELECT CASE (pos)
       
                    CASE(1) ! Node 1
                        BC_elem(el,4:9) = bc       
                    
                    CASE(2) ! Node 2
                        BC_elem(el,10:15) = bc
            
                    CASE(3) ! Node 3
                        BC_elem(el,16:21) = bc
        
                END SELECT

            ELSE
                
                SELECT CASE (pos)
       
                    CASE(1) ! Node 1
                        BC_elem(el,2:7) = bc

                    CASE(2) ! Node 2
                        BC_elem(el,8:13) = bc    
            
                    CASE(3) ! Node 3
                        BC_elem(el,14:19) = bc       

                END SELECT

            END IF
        END IF
    END DO

END SUBROUTINE Constrain_node
!===============================================================================!
SUBROUTINE nodes_conectivity_geo

	INTEGER :: nelem, i, nel, j, el, k, no, cont, jj, elb, kk, nob, max_cont, cont2
	REAL(8) :: x, y, z, xb, yb, zb, d
	INTEGER, ALLOCATABLE :: ELEM_aux(:,:)
	INTEGER, ALLOCATABLE :: nodes_conectivity_aux(:,:), nodes_conectivity_aux2(:,:)

	nelem = SIZE(ELEM_GEO_NODES_global,1)
	
	ALLOCATE(ELEM_aux(nelem,3))
	ELEM_aux = ELEM_GEO_NODES_global(:,2:4)

	ALLOCATE(nodes_conectivity_aux(nelem*nnos_el,IDcomp))
	nodes_conectivity_aux = 0

	DO i=1,nreg

		nel = El_reg(i,1)
		
		DO j=1,nel
			
			el = Subregions(i,j)
			
			DO k=1,nnos_el

				no = ELEM_aux(el,k)

				IF (no.NE.0) THEN

					x = GEO_NODES_global(no,2); 
            		y = GEO_NODES_global(no,3); 
            		z = GEO_NODES_global(no,4); 

					cont = 0
					nodes_conectivity_aux(no,1) = no
					
					DO jj=1,nel

						elb = Subregions(i,jj)

						IF (el.NE.elb) THEN
					
							DO kk=1,nnos_el

								nob = ELEM_aux(elb,kk)

								IF (nob.NE.0) THEN
								
									xb = GEO_NODES_global(nob,2); 
		        					yb = GEO_NODES_global(nob,3); 
		        					zb = GEO_NODES_global(nob,4); 
									
									d = SQRT((x-xb)**2 + (y-yb)**2 + (z-zb)**2)						

									IF (d.LT.(1e-9)) THEN

										cont = cont + 1
										
										nodes_conectivity_aux(no,cont+3) = nob

										ELEM_aux(elb,kk) = 0

										IF (el .EQ. Nodes_BCs(no,2)) THEN 
											
											nodes_conectivity_aux(no,3) = el

										END IF

										IF (elb .EQ. Nodes_BCs(nob,2)) THEN 

											nodes_conectivity_aux(no,3) = elb

										END IF

									END IF

								END IF

							END DO

						END IF

					END DO

					nodes_conectivity_aux(no,2) = cont
					ELEM_aux(el,k) = 0

				END IF

			END DO

		END DO

	END DO

	max_cont = MAXVAL(nodes_conectivity_aux(:,2))

	ALLOCATE(nodes_conectivity_aux2(nelem*nnos_el,max_cont+3))
	nodes_conectivity_aux2 = 0
	cont2 = 0
	DO i = 1,nelem*nnos_el

		no = nodes_conectivity_aux(i,1)

		IF (no.NE.0) THEN	

			cont2 = cont2 + 1
			cont = nodes_conectivity_aux(i,2)
			nodes_conectivity_aux2(cont2,1:cont+3) = nodes_conectivity_aux(i,1:cont+3)

		END IF

	END DO

	DEALLOCATE(nodes_conectivity_aux)

	ALLOCATE(nodes_conectivity(cont2,max_cont+3))

	nodes_conectivity = nodes_conectivity_aux2(1:cont2,:)

	DEALLOCATE(nodes_conectivity_aux2, ELEM_aux)
	
END SUBROUTINE nodes_conectivity_geo
!===============================================================================!
END MODULE BC_Module

