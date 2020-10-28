!===============================================================================!
!-------------------------------------------------------------------------------!
!                   MODULE TO WRITE THE DATA TO PLOT THE RESPONSE               !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Output
!-------------------------------------------------------------------------------!
USE Global_variables
USE Global_functions
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE Output_data_par(ts1,step)

    INTEGER::t0,t1,rate
    REAL(8),INTENT(OUT)::ts1

    INTEGER :: step, npoints, nelem, no, no_par(4), i, j
    REAL(8) :: xi, yi, zi, dx, dy, dz, xf, yf, zf, dispt, N(3)
    REAL(8) :: phi_x(3), phi_y(3), phi_z(3), geo_x(3), geo_y(3), geo_z(3)
    REAL(8), ALLOCATABLE :: POINTS_def(:,:), Displacement_geo(:,:)
    REAL(8) phi_xx(3), phi_yy(3), phi_zz(3), phi_yz(3), phi_xz(3), phi_xy(3)
    REAL(8) geo_xx(3), geo_yy(3), geo_zz(3), geo_yz(3), geo_xz(3), geo_xy(3)
    REAL(8), ALLOCATABLE :: Stress_geo(:,:), Strain_geo(:,:), VMS(:,:)
        
    CHARACTER(LEN=20)::num_str


    CALL SYSTEM_CLOCK(t0, rate) ! Initializates the time counter
    
    WRITE (*,*) ' '
    WRITE (*,*) '     -- Paraview' 
    WRITE (*,*) ' ' 

	IF ((Static_steps.GT.1).OR.(Time_steps.GT.1)) THEN

    	WRITE(num_str,"(I10)") step
    	Test_name = TRIM(Results_file)//'_'//TRIM(ADJUSTL(num_str))//".vtk"

	ELSEIF ((Static_steps.EQ.1).OR.(Time_steps.EQ.1)) THEN

		Test_name = TRIM(Results_file)//".vtk"

	END IF

    OPEN(3,file=trim(fileplace_results)//trim(Test_name),STATUS='REPLACE')

    REWIND(3)

    WRITE(3,'(A)') '# vtk DataFile Version 2.3'
    WRITE(3,'(A)') Test_name
    WRITE(3,'(A)') 'ASCII'
    WRITE(3,'(A)') ''
    WRITE(3,'(A)') 'DATASET UNSTRUCTURED_GRID'
    !--------------------------------------------------------------
    npoints = SIZE(GEO_NODES_global,1)
	nelem = SIZE(ELEM_GEO_NODES_global,1)

	ALLOCATE(Displacement_geo(npoints,3))

    DO i=1,nelem

        phi_x = Displacement(3*i-2:3*i,2)
        phi_y = Displacement(3*i-2:3*i,3)
        phi_z = Displacement(3*i-2:3*i,4)
        
        DO j=1,nnos_el

            N = N_Geo_nodes2(3*j-2:3*j,1)
            
            geo_x(j) = phi_x(1)*N(1) + phi_x(2)*N(2) + phi_x(3)*N(3)
            geo_y(j) = phi_y(1)*N(1) + phi_y(2)*N(2) + phi_y(3)*N(3) 
            geo_z(j) = phi_z(1)*N(1) + phi_z(2)*N(2) + phi_z(3)*N(3)

        END DO

        Displacement_geo(3*i-2:3*i,1) = geo_x
        Displacement_geo(3*i-2:3*i,2) = geo_y
        Displacement_geo(3*i-2:3*i,3) = geo_z
            
    END DO

	CALL connectivity_plot_disp(Displacement_geo)
		    
    ALLOCATE(POINTS_def(npoints,3))
    
    WRITE(3,'(A,1I20,A)') 'POINTS', npoints, ' double'
    
    DO i=1,nelem

        DO j=1,nnos_el
                
            no = ELEM_GEO_NODES_global(i,j+1)
            xi = GEO_NODES_global(no,2); 
            yi = GEO_NODES_global(no,3); 
            zi = GEO_NODES_global(no,4); 

            dx = Displacement_geo(no,1);
            dy = Displacement_geo(no,2);
            dz = Displacement_geo(no,3);

            xf = xi + dx*scale_results; yf = yi + dy*scale_results; zf = zi + dz*scale_results;

            POINTS_def(no,:) = (/xf, yf, zf/)

        END DO
    END DO

    DO i=1,npoints
        WRITE(3,*) POINTS_def(i,:)
    END DO
    DEALLOCATE(POINTS_def)
    !----------------------------------------------------------------
    WRITE(3,'(A,1I20,1I20)') 'CELLS', nelem, nelem*4

    no_par(1) = 3
   
    DO i=1,nelem
        
        no_par(2:4) = ELEM_GEO_NODES_global(i,2:4) - 1
        
        WRITE(3,*) no_par

    END DO
    !----------------------------------------------------------------
    WRITE(3,'(A,1I20)') 'CELL_TYPES', nelem
    
    DO i=1,nelem
        WRITE(3,*) 5
    END DO
    !----------------------------------------------------------------
    !----------------------------------------------------------------
    !----------------------------------------------------------------
    WRITE(3,'(A,1I20)') 'POINT_DATA', npoints
    WRITE(3,'(A)') 'FIELD Results 4'
    WRITE(3,'(A,1I20,A)') 'Disp 3', npoints, ' double'


    DO i=1,npoints
        dx = Displacement_geo(i,1);
        dy = Displacement_geo(i,2);
        dz = Displacement_geo(i,3);        
		
        dispt = SQRT(dx**2 + dy**2 + dz**2)

        WRITE(3,*) Displacement_geo(i,1:3)!, dispt

    END DO

    DEALLOCATE(Displacement_geo)

    !----------------------------------------------------------------
    
    ALLOCATE(Stress_geo(npoints,6),VMS(npoints,1))

    DO i=1,nelem

        phi_xx = Stress(3*i-2:3*i,1)
        phi_yy = Stress(3*i-2:3*i,2)
        phi_zz = Stress(3*i-2:3*i,3)
        phi_yz = Stress(3*i-2:3*i,4)
        phi_xz = Stress(3*i-2:3*i,5)
        phi_xy = Stress(3*i-2:3*i,6)
        
        DO j=1,nnos_el

            N = N_Geo_nodes2(3*j-2:3*j,1)
            
            geo_xx(j) = phi_xx(1)*N(1) + phi_xx(2)*N(2) + phi_xx(3)*N(3)
            geo_yy(j) = phi_yy(1)*N(1) + phi_yy(2)*N(2) + phi_yy(3)*N(3) 
            geo_zz(j) = phi_zz(1)*N(1) + phi_zz(2)*N(2) + phi_zz(3)*N(3)
            geo_yz(j) = phi_yz(1)*N(1) + phi_yz(2)*N(2) + phi_yz(3)*N(3)
            geo_xz(j) = phi_xz(1)*N(1) + phi_xz(2)*N(2) + phi_xz(3)*N(3) 
            geo_xy(j) = phi_xy(1)*N(1) + phi_xy(2)*N(2) + phi_xy(3)*N(3)
            

        END DO

        Stress_geo(3*i-2:3*i,1) = geo_xx
        Stress_geo(3*i-2:3*i,2) = geo_yy
        Stress_geo(3*i-2:3*i,3) = geo_zz
        Stress_geo(3*i-2:3*i,4) = geo_yz
        Stress_geo(3*i-2:3*i,5) = geo_xz
        Stress_geo(3*i-2:3*i,6) = geo_xy

	!Stress_geo(3*i-2:3*i,1) = phi_xx
        !Stress_geo(3*i-2:3*i,2) = phi_yy
        !Stress_geo(3*i-2:3*i,3) = phi_zz
        !Stress_geo(3*i-2:3*i,4) = phi_yz
        !Stress_geo(3*i-2:3*i,5) = phi_xz
        !Stress_geo(3*i-2:3*i,6) = phi_xy

		

		VMS(3*i-2:3*i,1) = (0.5d0*((phi_xx-phi_yy)**2 + (phi_yy-phi_zz)**2 + &
				 (phi_zz-phi_xx)**2 + 6.d0*(phi_yz**2 + phi_xz**2 + phi_xy**2)))**(0.5d0)

    END DO

	CALL connectivity_plot(Stress_geo)

    WRITE(3,'(A,1I20,A)') 'Sigma_ij 6', npoints, ' double'

    DO i=1,npoints
        
        WRITE(3,*) Stress_geo(i,1:6)/scale_prop_mat

    END DO
    
    !WRITE(3,'(A,1I20,A)') 'Sigma_ij 3', npoints, ' double'
    
    !DO i=1,npoints
        
    !    WRITE(3,*) Stress_geo(i,4:6)*scale_prop_mat

    !END DO

	CALL connectivity_plot(VMS)

	WRITE(3,'(A,1I20,A)') 'sigma_VM 1', npoints, ' double'
    
    DO i=1,npoints
        
        WRITE(3,*) VMS(i,1)/scale_prop_mat

    END DO

    DEALLOCATE(Stress_geo)    

    !----------------------------------------------------------------
    ALLOCATE(Strain_geo(npoints,6))
    DO i=1,nelem

        phi_xx = Strain(3*i-2:3*i,1)
        phi_yy = Strain(3*i-2:3*i,2)
        phi_zz = Strain(3*i-2:3*i,3)
        phi_yz = Strain(3*i-2:3*i,4)
        phi_xz = Strain(3*i-2:3*i,5)
        phi_xy = Strain(3*i-2:3*i,6)
        
        DO j=1,nnos_el

            N = N_Geo_nodes2(3*j-2:3*j,1)
            
            geo_xx(j) = phi_xx(1)*N(1) + phi_xx(2)*N(2) + phi_xx(3)*N(3)
            geo_yy(j) = phi_yy(1)*N(1) + phi_yy(2)*N(2) + phi_yy(3)*N(3) 
            geo_zz(j) = phi_zz(1)*N(1) + phi_zz(2)*N(2) + phi_zz(3)*N(3)
            geo_yz(j) = phi_yz(1)*N(1) + phi_yz(2)*N(2) + phi_yz(3)*N(3)
            geo_xz(j) = phi_xz(1)*N(1) + phi_xz(2)*N(2) + phi_xz(3)*N(3) 
            geo_xy(j) = phi_xy(1)*N(1) + phi_xy(2)*N(2) + phi_xy(3)*N(3)
            

        END DO

        Strain_geo(3*i-2:3*i,1) = geo_xx
        Strain_geo(3*i-2:3*i,2) = geo_yy
        Strain_geo(3*i-2:3*i,3) = geo_zz
        Strain_geo(3*i-2:3*i,4) = geo_yz
        Strain_geo(3*i-2:3*i,5) = geo_xz
        Strain_geo(3*i-2:3*i,6) = geo_xy
            
    END DO

	CALL connectivity_plot(Strain_geo)

    WRITE(3,'(A,1I20,A)') 'Epsilon_ij 6', npoints, ' double'
    
    DO i=1,npoints
        
        WRITE(3,*) Strain_geo(i,1:6)

    END DO

    !----------------------------------------------------------------
    !WRITE(3,'(A,1I20,A)') 'Epsilon_ij 3', npoints, ' double'
    
    !DO i=1,npoints
        
    !    WRITE(3,*) Strain_geo(i,4:6)

    !END DO
    DEALLOCATE(Strain_geo) 
        

            

    CLOSE(3) 

    CALL SYSTEM_CLOCK(t1)       ! Time counter stopped
    WRITE (*,'(A,F15.3,A)') '                                  ...COMPLETED! (Time :', REAL(t1-t0)/rate,'s)'
    WRITE (*,*) ''
    ts1 = REAL(t1-t0)/rate

END SUBROUTINE Output_data_par
!====================================================================================!
END MODULE Output
