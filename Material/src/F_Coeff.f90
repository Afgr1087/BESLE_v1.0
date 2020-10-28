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
! Program: Fourier_Coefficients                                                      !
! Description: PROGRAM TO FIND THE FOURIER COEFFICIENTS                              !
!====================================================================================!
!====================================================================================!
PROGRAM Main
!====================================================================================!
! MODULES
!------------------------------------------------------------------------------------!
    USE Global_variables 
	USE Set_parameters
    USE Fourier_Coefficient
    USE Output_data
!====================================================================================!
! VARIABLES
!------------------------------------------------------------------------------------!
    IMPLICIT NONE
    
    INTEGER :: i
    REAL :: ts1, ts2
    CHARACTER(LEN=20)::num_str
    REAL(8), ALLOCATABLE :: Angles(:,:)
    
!====================================================================================!
! INPUT DATA
!------------------------------------------------------------------------------------!

    CALL start_greet
    
	CALL Setup

	IF (Material.EQ.'Multiple_iso') THEN

		OPEN(1,file=trim(fileplace)//file_name,STATUS='OLD')

		REWIND(1)

		ALLOCATE(E_v_constants(n_materials,2))

		DO i = 1,n_materials

			READ (1,*) E_v_constants(i,:)

		END DO
		
		CLOSE(1)

	END IF

	IF (Material.EQ.'Multiple_aniso') THEN

		OPEN(1,file=trim(fileplace)//file_name,STATUS='OLD')

		REWIND(1)

		ALLOCATE(C_tensors(n_materials,21))

		DO i = 1,n_materials

			READ (1,*) C_tensors(i,:)

		END DO
		
		CLOSE(1)

	END IF

	IF ((n_materials.GT.1).AND.(Material.EQ.'Anisotropic')) THEN
		
		OPEN(1,file=trim(fileplace)//file_name,STATUS='OLD')

		REWIND(1)
				
		ALLOCATE(Angles(n_materials,3))
		DO i = 1,n_materials

			IF (z_x_z.EQ.1) THEN

				x_y_z = 0

				READ (1,*) theta_x, theta_z
				Angles(i,:) = (/theta_z,theta_x,theta_z/)*(Pi/180.d0)					

			ELSEIF (x_y_z.EQ.1) THEN

				READ (1,*) Angles(i,:)
				Angles(i,:) = Angles(i,:)*(Pi/180.d0)

			END IF
			

		END DO
		
		CLOSE(1)

	ELSE IF ((n_materials.EQ.1).AND.(Material.EQ.'Anisotropic')) THEN

		ALLOCATE(Angles(n_materials,3))
	
		IF (z_x_z.EQ.1) THEN
			x_y_z = 0
		END IF


		IF (z_x_z.EQ.1) THEN
			Angles(1,:) = (/theta_z, theta_x, theta_z/)*(Pi/180.d0)
		ELSE IF (x_y_z.EQ.1) THEN
			Angles(1,:) = (/theta_x, theta_y, theta_z/)*(Pi/180.d0)
		END IF

	END IF


    DO i=1,n_materials

		step = i

		IF (n_materials.GT.1) THEN
	
		    WRITE(num_str,"(I10)") i-1
		    file_name = 'Mat_'//TRIM(ADJUSTL(num_str))//".dat"

		END IF

		IF ((n_materials.GE.1).AND.(Material.EQ.'Anisotropic')) THEN
        
		    phi_1 = Angles(i,1)
		    phi = Angles(i,2)
		    phi_2 = Angles(i,3)

		    IF (z_x_z.EQ.1) THEN

				WRITE (*,*) ' '
		        WRITE (*,*) '   fileplace   :', '   ', trim(fileplace)//trim(file_name)
				WRITE (*,*) ' '
		        WRITE (*,'(A,F7.2)') '    theta_z: ', phi_1*(180.d0/Pi)
		        WRITE (*,'(A,F7.2)') '    theta_x: ', phi*(180.d0/Pi)
		        WRITE (*,'(A,F7.2)') '    theta_z: ', phi_2*(180.d0/Pi)
				WRITE (*,*) ' '

		    ELSE IF (x_y_z.EQ.1) THEN

				WRITE (*,*) ' '
		        WRITE (*,*) '   filepleace   :', '   ', trim(fileplace)//trim(file_name)
				WRITE (*,*) ' '
		        WRITE (*,'(A,F7.2)') '    theta_x: ', phi_1*(180.d0/Pi)
		        WRITE (*,'(A,F7.2)') '    theta_y: ', phi*(180.d0/Pi)
		        WRITE (*,'(A,F7.2)') '    theta_z: ', phi_2*(180.d0/Pi)
				WRITE (*,*) ' '

			END IF

		END IF

		IF (Material.EQ.'Isotropic') THEN

			WRITE (*,*) ' '
		    WRITE (*,*) '   fileplace   :', '   ', trim(fileplace)//trim(file_name)
			WRITE (*,*) ' '

		END IF

		IF ((Material.EQ.'Multiple_aniso').OR.(Material.EQ.'Multiple_iso')) THEN

			WRITE (*,*) ' '
		    WRITE (*,*) '   fileplace   :', '   ', trim(fileplace)//trim(file_name)
			WRITE (*,*) ' '

		END IF

        !1. Computing the Fourier Coefficients
        CALL Fourier_Coeff(ts1)

        !2. Data for the main 3D-BEM program
        CALL Print_data(ts2)

    END DO

END PROGRAM Main
