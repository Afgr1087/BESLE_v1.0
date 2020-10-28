!===============================================================================!
!-------------------------------------------------------------------------------!
!            MODULE TO COMPUTE THE FOURIER COEFFICIENT FOR THE                  !
!                            FUNDAMENTAL SOLUTIONS                              !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Fourier_Coefficient
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE start_greet
    WRITE(*,*) ''
    WRITE(*,*) '*****************************************************************************'
    WRITE(*,*) '*****************************************************************************'
    WRITE(*,*) '*****************************************************************************'
    WRITE(*,*) '*****************************************************************************'
    WRITE(*,*) '              SOFTWARE FOR COMPUTING THE FOURIER COEFFICIENTS                '
    WRITE(*,*) '' 
	WRITE(*,*) 'By Andrés F. Galvis and Daniel M. Prada'
	WRITE(*,*) 'Department of Computational Mechanics'
	WRITE(*,*) 'School of Mechanical Engineering'
	WRITE(*,*) 'University of Campinas'
	WRITE(*,*) '18/02/2020'
    WRITE(*,*) '*****************************************************************************'
    WRITE(*,*) '*****************************************************************************'
    WRITE(*,*) '*****************************************************************************'
    WRITE(*,*) '*****************************************************************************'
    END SUBROUTINE start_greet
!===============================================================================!
    SUBROUTINE Fourier_Coeff(ts1)
        INTEGER::t0,t1,rate
        REAL,INTENT(OUT)::ts1
        INTEGER :: ngp=64, alpha = 16! Number of Gauss points and terms
        REAL(8) :: xgp(64), wgp(64)
        
        CALL SYSTEM_CLOCK(t0, rate) ! Initializates the time counter
        WRITE (*,*) '----------------------------------------------------------------------'
        WRITE (*,*) ' '
        WRITE (*,'(A)',advance='no') ' 1. Computing Fourier Coefficients       ...'

        CALL Compliance_tensor

        CALL Gauss_Legendre(-1.d0,1.d0,ngp,xgp,wgp)
      
        CALL Rot_C
        
        CALL RI_mat(alpha,xgp,wgp)

        CALL SYSTEM_CLOCK(t1)       ! Time counter stopped
        WRITE (*,'(A,F7.3,A)') 'COMPLETED! (Time :', REAL(t1-t0)/rate,'s)'
        WRITE (*,*) ''
        ts1 = REAL(t1-t0)/rate

    END SUBROUTINE Fourier_Coeff
!===============================================================================!
    SUBROUTINE Compliance_tensor

        REAL(8) :: cte
        
        IF (step.EQ.1) THEN
            IF (Material.EQ.'Anisotropic') THEN
                SELECT CASE(Lattice)

                    CASE('cubic')

                        C21 = C12; C13 = C12; C31 = C12; C23 = C12; C32 = C12;
                        C22 = C11; C33 = C11; C55 = C44; C66 = C44;
                        C14 = 0.d0; C15 = 0.d0; C16 = 0.d0; C24 = 0.d0; C25 = 0.d0; C26 = 0.d0;
                        C34 = 0.d0; C35 = 0.d0; C36 = 0.d0; C45 = 0.d0; C46 = 0.d0; C56 = 0.d0;
                        C41 = 0.d0; C42 = 0.d0; C43 = 0.d0; C51 = 0.d0; C52 = 0.d0; C53 = 0.d0;
                        C54 = 0.d0; C61 = 0.d0; C62 = 0.d0; C63 = 0.d0; C64 = 0.d0; C65 = 0.d0;
                    CASE('hcp')

                        C21 = C12; C31 = C13; C23 = C13; C32 = C13;
                        C22 = C11; C55 = C44; C66 = 0.5d0*(C11-C12);
                        C14 = 0.d0; C15 = 0.d0; C16 = 0.d0; C24 = 0.d0; C25 = 0.d0; C26 = 0.d0;
                        C34 = 0.d0; C35 = 0.d0; C36 = 0.d0; C45 = 0.d0; C46 = 0.d0; C56 = 0.d0;
                        C41 = 0.d0; C42 = 0.d0; C43 = 0.d0; C51 = 0.d0; C52 = 0.d0; C53 = 0.d0;
                        C54 = 0.d0; C61 = 0.d0; C62 = 0.d0; C63 = 0.d0; C64 = 0.d0; C65 = 0.d0;
                    CASE('trigonal')

                        C21 = C12; C31 = C13; C41 = C14; C23 = C13; C32 = C13; C22 = C11
                        C24 = -C14; C42 = -C14; C55 = C44; C56 = C14; C65 = C14;
                        C66 = 0.5d0*(C11-C12);
                        C15 = 0.d0; C16 = 0.d0; C25 = 0.d0; C26 = 0.d0;
                        C34 = 0.d0; C35 = 0.d0; C36 = 0.d0; C45 = 0.d0; C46 = 0.d0;
                        C43 = 0.d0; C51 = 0.d0; C52 = 0.d0; C53 = 0.d0;
                        C54 = 0.d0; C61 = 0.d0; C62 = 0.d0; C63 = 0.d0; C64 = 0.d0;

					CASE('full')
			
						C21=C12; C31=C13; C41=C14; C51=C15; C61=C16;
						         C32=C23; C42=C24; C52=C25; C62=C26;
										  C43=C34; C53=C35; C63=C36;
												   C54=C45; C64=C46;
												   			C65=C56;


                END SELECT
            ELSEIF (Material.EQ.'Isotropic') THEN
                cte = E/((1.d0+nu)*(1.d0-2.d0*nu));
                C11 = cte*(1.d0-nu); C12 = cte*(nu); C13 = cte*(nu); C44 = cte*((1.d0-2.d0*nu)/2);
                C21 = C12; C22 = C11; C23 = C13; C31 = C13; C32 = C13; C33 = C11;
                C55 = C44; C66 = C44;
                C14 = 0.d0; C15 = 0.d0; C16 = 0.d0; C24 = 0.d0; C25 = 0.d0; C26 = 0.d0;
                C34 = 0.d0; C35 = 0.d0; C36 = 0.d0; C45 = 0.d0; C46 = 0.d0; C56 = 0.d0;
                C41 = 0.d0; C42 = 0.d0; C43 = 0.d0; C51 = 0.d0; C52 = 0.d0; C53 = 0.d0;
                C54 = 0.d0; C61 = 0.d0; C62 = 0.d0; C63 = 0.d0; C64 = 0.d0; C65 = 0.d0;
            END IF
            
			C_original(1,:) = (/C11, C12, C13, C14, C15, C16/)
            C_original(2,:) = (/C21, C22, C23, C24, C25, C26/)
            C_original(3,:) = (/C31, C32, C33, C34, C35, C36/)
            C_original(4,:) = (/C41, C42, C43, C44, C45, C46/)
            C_original(5,:) = (/C51, C52, C53, C54, C55, C56/)
            C_original(6,:) = (/C61, C62, C63, C64, C65, C66/) 


        END IF

		IF (Material.EQ.'Multiple_aniso') THEN

		    C11 = C_tensors(step,1); C12 = C_tensors(step,2); C13 = C_tensors(step,3)
		    C14 = C_tensors(step,4); C15 = C_tensors(step,4); C16 = C_tensors(step,6)
			C22 = C_tensors(step,7); C23 = C_tensors(step,8); C24 = C_tensors(step,9)
			C25 = C_tensors(step,10); C26 = C_tensors(step,11); C33 = C_tensors(step,12)
			C34 = C_tensors(step,13); C35 = C_tensors(step,14); C36 = C_tensors(step,15)
			C44 = C_tensors(step,16); C45 = C_tensors(step,17); C46 = C_tensors(step,18)
			C55 = C_tensors(step,19); C56 = C_tensors(step,20); C66 = C_tensors(step,21)
		
			C21=C12; C31=C13; C41=C14; C51=C15; C61=C16;
				     C32=C23; C42=C24; C52=C25; C62=C26;
							  C43=C34; C53=C35; C63=C36;
									   C54=C45; C64=C46;
												C65=C56;

			C_original(1,:) = (/C11, C12, C13, C14, C15, C16/)
            C_original(2,:) = (/C21, C22, C23, C24, C25, C26/)
            C_original(3,:) = (/C31, C32, C33, C34, C35, C36/)
            C_original(4,:) = (/C41, C42, C43, C44, C45, C46/)
            C_original(5,:) = (/C51, C52, C53, C54, C55, C56/)
            C_original(6,:) = (/C61, C62, C63, C64, C65, C66/) 


		END IF

		IF (Material.EQ.'Multiple_iso') THEN
		
			E = E_v_constants(step,1); nu = E_v_constants(step,2); 

			cte = E/((1.d0+nu)*(1.d0-2.d0*nu));
            C11 = cte*(1.d0-nu); C12 = cte*(nu); C13 = cte*(nu); C44 = cte*((1.d0-2.d0*nu)/2);
            C21 = C12; C22 = C11; C23 = C13; C31 = C13; C32 = C13; C33 = C11;
            C55 = C44; C66 = C44;

			C_original(1,:) = (/C11, C12, C13, C14, C15, C16/)
            C_original(2,:) = (/C21, C22, C23, C24, C25, C26/)
            C_original(3,:) = (/C31, C32, C33, C34, C35, C36/)
            C_original(4,:) = (/C41, C42, C43, C44, C45, C46/)
            C_original(5,:) = (/C51, C52, C53, C54, C55, C56/)
            C_original(6,:) = (/C61, C62, C63, C64, C65, C66/) 

		END IF

        C = C_original
        
    END SUBROUTINE Compliance_tensor
!===============================================================================!
    SUBROUTINE Rot_C

        REAL(8) :: gamma(9,3), gamma_axis(3,3)
        REAL(8) :: g11, g12, g13, g21, g22, g23, g31, g32, g33
        REAL(8) :: K1(3,3), K2(3,3), K3(3,3), K4(3,3), K(6,6), K_t(6,6)
        REAL(8) :: C_aux1(6,6), C_aux2(6,6), m, n
        INTEGER :: i

                
        C_aux1 = C
				

        
        IF (z_x_z.EQ.1) THEN

            ! Rotation in z
            m = DCOS(phi_1); n = DSIN(phi_1)
            gamma(1,:) = (/   m,   n,  0.d0/)
            gamma(2,:) = (/  -n,   m,  0.d0/)
            gamma(3,:) = (/ 0.d0, 0.d0,1.d0 /)

            ! Rotation in x
            m = DCOS(phi); n = DSIN(phi)
            gamma(4,:) = (/1.d0,0.d0,0.d0/)
            gamma(5,:) = (/0.d0, m,   n  /)
            gamma(6,:) = (/0.d0,-n,   m  /)

            ! Rotation in z
            m = DCOS(phi_2); n = DSIN(phi_2)
            gamma(7,:) = (/   m,   n,  0.d0/)
            gamma(8,:) = (/  -n,   m,  0.d0/)
            gamma(9,:) = (/ 0.d0, 0.d0,1.d0 /)

        ELSE IF (x_y_z.EQ.1) THEN

            ! Rotation in x
            m = DCOS(phi_1); n = DSIN(phi_1)
            gamma(1,:) = (/1.d0,0.d0,0.d0/)
            gamma(2,:) = (/0.d0, m,   n  /)
            gamma(3,:) = (/0.d0,-n,   m  /)

            ! Rotation in y
            m = DCOS(phi); n = DSIN(phi)
            gamma(4,:) = (/ m,  0.d0, -n /)
            gamma(5,:) = (/0.d0,1.d0,0.d0/)
            gamma(6,:) = (/ n,  0.d0,  m /)

            ! Rotation in z
            m = DCOS(phi_2); n = DSIN(phi_2)
            gamma(7,:) = (/   m,   n,  0.d0/)
            gamma(8,:) = (/  -n,   m,  0.d0/)
            gamma(9,:) = (/ 0.d0, 0.d0,1.d0 /)

		END IF

        DO i=1,3 ! Three axes of rotation x,y and z

            gamma_axis = gamma(3*i-2:3*i,:)

            g11 = gamma_axis(1,1); g12 = gamma_axis(1,2); g13 = gamma_axis(1,3)
            g21 = gamma_axis(2,1); g22 = gamma_axis(2,2); g23 = gamma_axis(2,3)
            g31 = gamma_axis(3,1); g32 = gamma_axis(3,2); g33 = gamma_axis(3,3)

            K1(1,:) = (/g11**2,g12**2,g13**2/)
            K1(2,:) = (/g21**2,g22**2,g23**2/)
            K1(3,:) = (/g31**2,g32**2,g33**2/)

            K2(1,:) = (/g12*g13, g13*g11, g11*g12/)
            K2(2,:) = (/g22*g23, g23*g21, g21*g22/)
            K2(3,:) = (/g32*g33, g33*g31, g31*g32/)

            K3(1,:) = (/g21*g31, g22*g32, g23*g33/)
            K3(2,:) = (/g31*g11, g32*g12, g33*g13/)
            K3(3,:) = (/g11*g21, g12*g22, g13*g23/)

            K4(1,:) = (/g22*g33+g23*g32, g23*g31+g21*g33, g21*g32+g22*g31/)
            K4(2,:) = (/g32*g13+g33*g12, g33*g11+g31*g13, g31*g12+g32*g11/)
            K4(3,:) = (/g12*g23+g13*g22, g13*g21+g11*g23, g11*g22+g12*g21/)

            K(1:3,1:3) = K1; K(1:3,4:6) = 2.d0*K2
            K(4:6,1:3) = K3; K(4:6,4:6) = K4

            K_t = TRANSPOSE(K)

            C_aux2 = MATMUL(K,MATMUL(C_aux1,K_t))
            C_aux1 = C_aux2

        END DO

        C = C_aux1

    END SUBROUTINE Rot_C
!===============================================================================!
   SUBROUTINE Gauss_Legendre(x1,x2,n,x,w)

        INTEGER :: m, n, i, j, cont
        REAL(8) :: eps, xm, xl, x1, x2, x(:), w(:), z, p1, p2, p3, pp, z1

        eps = 3e-15
        m = NINT((n+1.d0)/2.d0)
        xm = 0.5d0*(x2+x1)
        xl = 0.5d0*(x2-x1)
        
        DO i=1,m
            z = DCOS(Pi*(i-0.25d0)/(n+0.5d0))
            cont = 0
            DO WHILE (1 == 1)
                cont = cont + 1
                p1 = 1.d0
                p2 = 0.d0
                DO j=1,n
                    p3 = p2
                    p2 = p1
                    p1 = ((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
                END DO
                pp=n*(z*p1-p2)/(z*z-1.d0)
                z1=z;
                z=z1-p1/pp
                IF (ABS(z-z1).LT.eps) THEN
                    EXIT
                END IF
                IF (cont == 20) EXIT
            END DO
            x(i)=xm-xl*z;
            x(n+1-i)=xm+xl*z;
            w(i)=2.d0*xl/((1-z*z)*pp*pp);
            w(n+1-i)=w(i);
        END DO

    END SUBROUTINE Gauss_Legendre
!===============================================================================!
    SUBROUTINE RI_mat(alpha,xgp,wgp)


        INTEGER :: alpha, m, n, i,j
        REAL(8) :: max_mat(alpha,alpha), max_vet(alpha), xgp(:), wgp(:)
        REAL(8) :: Ruv1(3,3), Ruv2(3,3), Iuv1(3,3), Iuv2(3,3)
        REAL(8) :: R0m(3,3), Rm0(3,3), I0m(3,3), Im0(3,3)
        
        ALLOCATE(Rt_mat(alpha,alpha), Rc_mat(alpha,alpha))
        ALLOCATE(It_mat(alpha,alpha), Ic_mat(alpha,alpha))
        ALLOCATE(R0m_vet(alpha), Rm0_vet(alpha))
        ALLOCATE(I0m_vet(alpha), Im0_vet(alpha))

        DO m=1,alpha
            DO n=1,alpha

                CALL lambda_uv(m,n,xgp,wgp,Ruv1,Iuv1)
                CALL lambda_uv(m,-n,xgp,wgp,Ruv2,Iuv2)

                Rt_mat(m,n)%comp = Ruv1 + Ruv2
                Rc_mat(m,n)%comp = Ruv1 - Ruv2

                It_mat(m,n)%comp = Iuv1 + Iuv2
                Ic_mat(m,n)%comp = Iuv1 - Iuv2

                max_mat(m,n) = MAX(NORM2(Rt_mat(m,n)%comp),NORM2(Rc_mat(m,n)%comp), &
                               & NORM2(It_mat(m,n)%comp),NORM2(Ic_mat(m,n)%comp))
            END DO

            CALL lambda_uv(0,m,xgp,wgp,R0m,I0m)
            CALL lambda_uv(m,0,xgp,wgp,Rm0,Im0)

            R0m_vet(m)%comp = R0m
            Rm0_vet(m)%comp = Rm0

            I0m_vet(m)%comp = I0m
            Im0_vet(m)%comp = Im0

            max_vet(m) = MAXVAL((/max_mat(m,:),NORM2(R0m_vet(m)%comp),NORM2(I0m_vet(m)%comp), &
                         & NORM2(Rm0_vet(m)%comp),NORM2(Im0_vet(m)%comp)/))

       END DO

       max_val = MAXVAL(max_vet)

       CALL lambda_uv(0,0,xgp,wgp,Ruv1,Iuv1)

       R00_matrix = Ruv1
     
       ALLOCATE(Rt_matrix(alpha*3,alpha*3))
       ALLOCATE(It_matrix(alpha*3,alpha*3))
       ALLOCATE(Rc_matrix(alpha*3,alpha*3))
       ALLOCATE(Ic_matrix(alpha*3,alpha*3))

       ALLOCATE(R0m_vector(alpha*3,3), Rm0_vector(alpha*3,3))
       ALLOCATE(I0m_vector(alpha*3,3), Im0_vector(alpha*3,3))
        
       DO i=1,alpha
           DO j=1,alpha
               Rt_matrix(3*i-2:3*i,3*j-2:3*j) = Rt_mat(i,j)%comp
               Rc_matrix(3*i-2:3*i,3*j-2:3*j) = Rc_mat(i,j)%comp
               It_matrix(3*i-2:3*i,3*j-2:3*j) = It_mat(i,j)%comp
               Ic_matrix(3*i-2:3*i,3*j-2:3*j) = Ic_mat(i,j)%comp
           END DO

           R0m_vector(3*i-2:3*i,:) = R0m_vet(i)%comp
           I0m_vector(3*i-2:3*i,:) = I0m_vet(i)%comp
           Rm0_vector(3*i-2:3*i,:) = Rm0_vet(i)%comp
           Im0_vector(3*i-2:3*i,:) = Im0_vet(i)%comp
       END DO

    END SUBROUTINE RI_mat
!===============================================================================!
    SUBROUTINE lambda_uv(m,n,xgp,wgp,Ruv1,Iuv1)

        INTEGER :: m, n, i, j
        REAL(8) :: xgp(:), wgp(:), Ruv1(:,:), Iuv1(:,:)
        COMPLEX(8) :: l_uv(3,3), f(3,3)
        l_uv = 0.d0; f = 0.d0

        DO i=1,SIZE(xgp)
            DO j=1,SIZE(xgp)
                CALL fuv(xgp(i)*Pi,xgp(j)*Pi,m,n,f)
                l_uv = l_uv + wgp(i)*wgp(j)*f
            END DO
        END DO

        l_uv = l_uv/4.d0

        Ruv1 = REAL(l_uv)
        Iuv1 = AIMAG(l_uv)

    END SUBROUTINE lambda_uv
!===============================================================================!
    SUBROUTINE fuv(theta,phi,m,n,f)

        INTEGER :: m, n
        REAL(8) :: theta, phi, H(3,3)
        COMPLEX(8) :: f(3,3), cte

        cte = DCMPLX(0.d0,-1.d0)
		
        CALL Huv(theta,phi,H)
        f = H*EXP(cte*(m*theta + n*phi))

    END SUBROUTINE fuv
!===============================================================================!
    SUBROUTINE Huv(theta,phi,H)

        REAL(8) :: theta, phi, m1, m2, n1, n2, n3
        REAL(8) :: K11, K12, K13, K22, K23, K33, detk
        REAL(8) :: R11, R12, R13, R21, R22, R23, R31, R32, R33
        REAL(8) :: V11, V12, V13, V22, V23, V33
        REAL(8) :: Q11, Q12, Q13, Q22, Q23, Q33
        REAL(8) :: fp0, fp1, fp2, fp3, fp4, fp5, fp6, a(7)
        REAL(8) :: beta1, beta2, beta3, cte, q0, q1, q2, q3, q4, H(3,3)
        REAL(8) :: H_bl11, H_bl12, H_bl13, H_bl22, H_bl23, H_bl33
        COMPLEX(8) :: x(6), p1, p2, p3, p1_conj, p2_conj, p3_conj
        COMPLEX(8) :: t1, t2, t3

        ! Compliance tensor
        C11 = C(1,1); C12 = C(1,2); C13 = C(1,3); C14 = C(1,4); C15 = C(1,5);
        C16 = C(1,6); C22 = C(2,2); C23 = C(2,3); C24 = C(2,4); C25 = C(2,5);
        C26 = C(2,6); C33 = C(3,3); C34 = C(3,4); C35 = C(3,5); C36 = C(3,6);
        C44 = C(4,4); C45 = C(4,5); C46 = C(4,6); C55 = C(5,5); C56 = C(5,6);
        C66 = C(6,6);
		
        m1 = -DSIN(theta); m2 = DCOS(theta); ! m3 = 0, values involving m3 were cut off
        n1 = DCOS(phi)*DCOS(theta); n2 = DCOS(phi)*DSIN(theta); n3 = -DSIN(phi);

        ! Matriz [K]
        K11 = C11*m1*m1 + 2.d0*C16*m1*m2 + C66*m2*m2;
        K12 = C16*m1*m1 + (C12 + C66)*m1*m2 + C26*m2*m2;
        K13 = C15*m1*m1 + (C14 + C56)*m1*m2 + C46*m2*m2;
        K22 = C66*m1*m1 + 2.d0*C26*m1*m2 + C22*m2*m2;
        K23 = C56*m1*m1 + (C46 + C25)*m1*m2 + C24*m2*m2;
        K33 = C55*m1*m1 + (C45 + C45)*m1*m2 + C44*m2*m2;

        detK = -K33*K12**2 + 2.d0*K12*K13*K23 - K22*K13**2 - K11*K23**2 + K11*K22*K33

        ! Matriz [Q]
        Q11 = C11*n1*n1 + 2.d0*C16*n1*n2 + 2.d0*C15*n1*n3 + C66*n2*n2 + 2.d0*C56*n2*n3 + &
               & C55*n3*n3;
        Q12 = C16*n1*n1 + (C12 + C66)*n1*n2 + (C14 + C56)*n1*n3 + C26*n2*n2 + &
               & (C46 + C25)*n2*n3 + C45*n3*n3;
        Q13 = C15*n1*n1 + (C14 + C56)*n1*n2 + (C13 + C55)*n1*n3 + C46*n2*n2 + &
               & (C36 + C45)*n2*n3 + C35*n3*n3;
        Q22 = C66*n1*n1 + 2.d0*C26*n1*n2 + 2.d0*C46*n1*n3 + C22*n2*n2 + 2*C24*n2*n3 + &
               & C44*n3*n3;
        Q23 = C56*n1*n1 + (C46 + C25)*n1*n2 + (C36 + C45)*n1*n3 + C24*n2*n2 + &
               & (C23 + C44)*n2*n3 + C34*n3*n3;
        Q33 = C55*n1*n1 + 2.d0*C45*n1*n2 + 2.d0*C35*n1*n3 + C44*n2*n2 + 2.d0*C34*n2*n3 + &
               & C33*n3*n3;

        ! [R] & [V]
        R11 = C11*n1*m1 + C16*n1*m2 + C16*n2*m1 + C66*n2*m2 + C15*n3*m1 + C56*n3*m2;
        R12 = C16*n1*m1 + C12*n1*m2 + C66*n2*m1 + C26*n2*m2 + C56*n3*m1 + C25*n3*m2;
        R13 = C15*n1*m1 + C14*n1*m2 + C56*n2*m1 + C46*n2*m2 + C55*n3*m1 + C45*n3*m2;
        R21 = C16*n1*m1 + C66*n1*m2 + C12*n2*m1 + C26*n2*m2 + C14*n3*m1 + C46*n3*m2;
        R22 = C66*n1*m1 + C26*n1*m2 + C26*n2*m1 + C22*n2*m2 + C46*n3*m1 + C24*n3*m2;
        R23 = C56*n1*m1 + C46*n1*m2 + C25*n2*m1 + C24*n2*m2 + C45*n3*m1 + C44*n3*m2;
        R31 = C15*n1*m1 + C56*n1*m2 + C14*n2*m1 + C46*n2*m2 + C13*n3*m1 + C36*n3*m2;
        R32 = C56*n1*m1 + C25*n1*m2 + C46*n2*m1 + C24*n2*m2 + C36*n3*m1 + C23*n3*m2;
        R33 = C55*n1*m1 + C45*n1*m2 + C45*n2*m1 + C44*n2*m2 + C35*n3*m1 + C34*n3*m2;

        V11 = 2.d0*R11; V12 = R12 + R21; V13 = R13 + R31; V22 = 2.d0*R22;
        V23 = R23 + R32; V33 = 2.d0*R33; ! V = R + R'

        ! STROH's EIGENVALUES
        fp6 =  -K33*K12**2 + 2.d0*K12*K13*K23 - K22*K13**2 - K11*K23**2 + K11*K22*K33;

        fp5 =  -V33*K12**2 + 2.d0*V23*K12*K13 + 2.d0*V13*K12*K23 - 2.d0*K33*V12*K12 - &
              & V22*K13**2 + 2.d0*V12*K13*K23 - 2.d0*K22*V13*K13 - V11*K23**2 - &
              & 2.d0*K11*V23*K23 + K11*K22*V33 + K11*K33*V22 + K22*K33*V11

        fp4 =  -Q33*K12**2 + 2.d0*Q23*K12*K13 + 2.d0*Q13*K12*K23 - 2.d0*V33*K12*V12 + &
              & 2.d0*K12*V13*V23 - 2.d0*K33*Q12*K12 - Q22*K13**2 + 2.d0*Q12*K13*K23 + &
              & 2.d0*K13*V12*V23 - 2.d0*V22*K13*V13 - 2.d0*K22*Q13*K13 - Q11*K23**2 + &
              & 2.d0*K23*V12*V13 - 2.d0*V11*K23*V23 - 2.d0*K11*Q23*K23 - K33*V12**2 - &
              & K22*V13**2 - K11*V23**2 + K11*K22*Q33 + K11*K33*Q22 + K22*K33*Q11 + &
              & K11*V22*V33 + K22*V11*V33 + K33*V11*V22

        fp3 =  2.d0*K12*Q13*V23 - V13**2*V22 - V12**2*V33 - V11*V23**2 + 2.d0*K12*Q23*V13 + &
             & 2.d0*K13*Q12*V23 - 2.d0*K13*Q13*V22 - 2.d0*K13*Q22*V13 + 2.d0*K13*Q23*V12 - &
             & 2.d0*K22*Q13*V13 + 2.d0*K23*Q12*V13 + 2.d0*K23*Q13*V12 - 2.d0*K11*Q23*V23 - &
             & 2.d0*K12*Q12*V33 - 2.d0*K12*Q33*V12 - 2.d0*K23*Q11*V23 - 2.d0*K23*Q23*V11 - &
             & 2.d0*K33*Q12*V12 + K11*Q22*V33 + K11*Q33*V22 + K22*Q11*V33 + &
             & K22*Q33*V11 + K33*Q11*V22 + K33*Q22*V11 + 2*V12*V13*V23 + V11*V22*V33

        fp2 = - K33*Q12**2 + 2.d0*K23*Q12*Q13 + 2.d0*K13*Q12*Q23 - 2.d0*V33*Q12*V12 + &
              & 2.d0*Q12*V13*V23 - 2.d0*K12*Q33*Q12 - K22*Q13**2 + 2.d0*K12*Q13*Q23 + &
              & 2.d0*Q13*V12*V23 - 2.d0*V22*Q13*V13 - 2.d0*K13*Q22*Q13 - K11*Q23**2 + &
              & 2.d0*Q23*V12*V13 - 2.d0*V11*Q23*V23 - 2.d0*K23*Q11*Q23 - Q33*V12**2 - &
              & Q22*V13**2 - Q11*V23**2 + K11*Q22*Q33 + K22*Q11*Q33 + &
              & K33*Q11*Q22 + Q11*V22*V33 + Q22*V11*V33 + Q33*V11*V22

        fp1 = - V33*Q12**2 + 2.d0*V23*Q12*Q13 + 2.d0*V13*Q12*Q23 - 2.d0*Q33*V12*Q12 - &
              & V22*Q13**2 + 2.d0*V12*Q13*Q23 - 2.d0*Q22*V13*Q13 - V11*Q23**2 - &
              & 2.d0*Q11*V23*Q23 + Q11*Q22*V33 + Q11*Q33*V22 + Q22*Q33*V11

        fp0 = - Q33*Q12**2 + 2.d0*Q12*Q13*Q23 - Q22*Q13**2 - Q11*Q23**2 + Q11*Q22*Q33

        ! Solving the sextic polynomial equation
        a = (/fp6,fp5,fp4,fp3,fp2,fp1,fp0/)
		
        CALL SexticRoots(a,x)
        
        ! Sorting Eigenvalues
        beta1 = AIMAG(x(1));
        beta2 = AIMAG(x(3));
        beta3 = AIMAG(x(5));

        p1 = x(1); p2 = x(3); p3 = x(5);
        p1_conj = x(2); p2_conj = x(4); p3_conj = x(6);

        cte = -1.d0/(2.d0*beta1*beta2*beta3)

        t1 = (p1 - p2_conj)*(p1 - p3_conj);
        t2 = (p2 - p3_conj)*(p2 - p1_conj);
        t3 = (p3 - p1_conj)*(p3 - p2_conj);

        q0 = REAL(1.d0/t1 + 1.d0/t2 + 1.d0/t3);
        q1 = REAL(p1/t1 + p2/t2 + p3/t3);
        q2 = REAL(p1**2/t1 + p2**2/t2 + p3**2/t3) - 1.d0;
        q3 = REAL(p1*p2_conj*p3_conj/t1 + p2*p3_conj*p1_conj/t2 + p3*p1_conj*p2_conj/t3);
        q4 = REAL(p1**2*p2_conj*p3_conj/t1 + p2**2*p3_conj*p1_conj/t2 + p3**2*p1_conj*p2_conj/t3);

        ! Matrix [H]

        H_bl11 = q0*(Q22*Q33 - Q23*Q23) + &
                 q1*(V22*Q33 + Q22*V33 - 2.d0*V23*Q23) + &
                 q2*(K22*Q33 + Q22*K33 + V22*V33 - 2.d0*K23*Q23 - V23*V23) + &
                 q3*(V22*K33 + K22*V33 - 2.d0*V23*K23) + &
                 q4*(K22*K33 - K23*K23);

        H_bl12 = q0*(Q23*Q13 - Q12*Q33) + &
                 q1*(V23*Q13 + Q23*V13 - V12*Q33 - Q12*V33) + &
                 q2*(K23*Q13 + Q23*K13 + V23*V13 - K12*Q33 - Q12*K33 - V12*V33) + &
                 q3*(V23*K13 + K23*V13 - V12*K33 - K12*V33) + &
                 q4*(K23*K13 - K12*K33);

        H_bl13 = q0*(Q12*Q23 - Q22*Q13) + &
                 q1*(V12*Q23 + Q12*V23 - V22*Q13 - Q22*V13) + &
                 q2*(K12*Q23 + Q12*K23 + V12*V23 - K22*Q13 - Q22*K13 - V22*V13) + &
                 q3*(V12*K23 + K12*V23 - V22*K13 - K22*V13) + &
                 q4*(K12*K23 - K22*K13);

        H_bl22 = q0*(Q33*Q11 - Q13*Q13) + &
                 q1*(V33*Q11 + Q33*V11 - 2.d0*V13*Q13)+ &
                 q2*(K33*Q11 + Q33*K11 + V33*V11 - 2.d0*K13*Q13 - V13*V13) + &
                 q3*(V33*K11 + K33*V11 - 2.d0*V13*K13) + &
                 q4*(K33*K11 - K13*K13);

        H_bl23 = q0*(Q13*Q12 - Q23*Q11) + &
                 q1*(V13*Q12 + Q13*V12 - V23*Q11 - Q23*V11) + &
                 q2*(K13*Q12 + Q13*K12 + V13*V12 - K23*Q11 - Q23*K11 - V23*V11) + &
                 q3*(V13*K12 + K13*V12 - V23*K11 - K23*V11) + &
                 q4*(K13*K12 - K23*K11);

        H_bl33 = q0*(Q11*Q22 - Q12*Q12) + &
                 q1*(V11*Q22 + Q11*V22 - 2.d0*V12*Q12) + &
                 q2*(K11*Q22 + Q11*K22 + V11*V22 - 2.d0*Q12*K12 - V12*V12) + &
                 q3*(V11*K22 + K11*V22 - 2.d0*V12*K12) + &
                 q4*(K11*K22 - K12*K12);

        H(1,:) = (/H_bl11, H_bl12, H_bl13/)
        H(2,:) = (/H_bl12, H_bl22, H_bl23/)
        H(3,:) = (/H_bl13, H_bl23, H_bl33/)

        H = H*cte/detK

    END SUBROUTINE Huv
!===============================================================================!
    SUBROUTINE SexticRoots(a,x)

        REAL(8) :: a(:), fp6, fp5, fp4, fp3, fp2, fp1, fp0
        REAL(8) :: Cm0, Cm1, Cm2, Cm3, Cm4, Cm5
        REAL(8) :: Comp_mat(6,6)
        COMPLEX(8) :: x(:)
        CHARACTER(1) :: JOBVL='N', JOBVR='N'
        INTEGER :: N=6, LDA=6, LDVL=1, LDVR=1, LWORK=3*6, INFO=0
        REAL(8) :: WR(6),WI(6),VL(1,6),VR(1,6),WORK(3*6)

        WR=0.d0; WI=0.d0; VL=0.d0; VR=0.d0; WORK=0.d0

        fp6 = a(1); fp5 = a(2);
        fp4 = a(3); fp3 = a(4); fp2 = a(5); fp1 = a(6); fp0 = a(7)

        Cm0 = -fp0/fp6
        Cm1 = -fp1/fp6
        Cm2 = -fp2/fp6
        Cm3 = -fp3/fp6
        Cm4 = -fp4/fp6
        Cm5 = -fp5/fp6


        Comp_mat(1,:) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, Cm0/)
        Comp_mat(2,:) = (/1.d0, 0.d0, 0.d0, 0.d0, 0.d0, Cm1/)
        Comp_mat(3,:) = (/0.d0, 1.d0, 0.d0, 0.d0, 0.d0, Cm2/)
        Comp_mat(4,:) = (/0.d0, 0.d0, 1.d0, 0.d0, 0.d0, Cm3/)
        Comp_mat(5,:) = (/0.d0, 0.d0, 0.d0, 1.d0, 0.d0, Cm4/)
        Comp_mat(6,:) = (/0.d0, 0.d0, 0.d0, 0.d0, 1.d0, Cm5/)
		
        CALL DGEEV(JOBVL,JOBVR,N,Comp_mat,LDA,WR,WI,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
		
        x(1) = DCMPLX(WR(1),WI(1))
        x(2) = DCMPLX(WR(2),WI(2))
        x(3) = DCMPLX(WR(3),WI(3))
        x(4) = DCMPLX(WR(4),WI(4))
        x(5) = DCMPLX(WR(5),WI(5))
        x(6) = DCMPLX(WR(6),WI(6))

    END SUBROUTINE SexticRoots
!===============================================================================!
END MODULE Fourier_Coefficient

