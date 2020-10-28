!===============================================================================!
!-------------------------------------------------------------------------------!
!                 MODULE TO EVALUATE STRAIN AND STRESS TENSORS                  !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Response_S_S
!-------------------------------------------------------------------------------!
USE Global_variables
USE Global_functions
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE Strain_Stress(ts1,step)

    INTEGER::t0,t1,rate
    REAL(8),INTENT(OUT)::ts1

    INTEGER :: step, i, n_el, el, j, nelem
    REAL(8) :: Normal_face(3), grad_u_el(3,9), Strain_el(3,6)
    REAL(8) :: sigma_el(3,6)
    
    CALL SYSTEM_CLOCK(t0, rate) ! Initializates the time counter
    WRITE (*,*) '----------------------------------------------------------------------------'
    WRITE (*,*) ' '
    WRITE (*,'(A)',advance='no') ' 10. Strain and Stress tensors    ...'

    IF (step.EQ.1) THEN
        nelem = SIZE(ELEM,1)   
        ALLOCATE(Strain(nelem*nnos_el,6),Stress(nelem*nnos_el,6))
    END IF
    
    DO i=1,nreg

        n_el = El_reg(i,1) ! Elements per region
        
        DO j = 1,n_el

            el = Subregions(i,j) ! element

            ! Normal vector of the deformed element
	    	Normal_face = NORMAL_VEC_Def(el,:)

            ! Gradient of displacement field 
            CALL gradient_u(i,el,Normal_face,grad_u_el,sigma_el)

            ! Strain and stress tensor 
            CALL strain_t(i,grad_u_el,Strain_el)    

            Strain(3*el-2:3*el,:) = Strain_el
            Stress(3*el-2:3*el,:) = sigma_el  

        END DO

    END DO

    !CALL constraint_elements

    !Stress = Stress / scale_prop_mat

    CALL SYSTEM_CLOCK(t1)       ! Time counter stopped
    WRITE (*,'(A,F15.3,A)') 'COMPLETED! (Time :', REAL(t1-t0)/rate,'s)'
    WRITE (*,*) ''
    ts1 = REAL(t1-t0)/rate

END SUBROUTINE Strain_Stress
!===============================================================================!
SUBROUTINE constraint_elements

    INTEGER :: nelem, i, j, type_bc, disp_val, element, ii, size_v, pos1, pos2
    INTEGER, ALLOCATABLE :: Vector1(:), Vector2(:)

    IF (nreg.EQ.1) THEN

        nelem =SIZE(ELEM,1)
        DO i=1,nelem
            DO j=1,nnos_el
                ! Displacement is known
                type_bc = INT(SQRT(BC_elem(i,2*j)**2 + BC_elem(i,2*j+2)**2 + BC_elem(i,2*j+4)**2))
                disp_val = SQRT(BC_elem(i,2*j+1)**2 + BC_elem(i,2*j+3)**2 + BC_elem(i,2*j+5)**2)
               
                IF ((type_bc.EQ.0).AND.(ABS(disp_val).LE.(1e-15))) THEN 
                    
                    Strain(3*i-2:3*i,:) = 0.d0    
                    Stress(3*i-2:3*i,:) = -Stress(3*i-2:3*i,:)

                END IF
            END DO
        END DO

    ELSE

        size_v = SIZE(Interfaces,1)
        ALLOCATE(Vector1(size_v),Vector2(size_v))
        Vector1 = Interfaces(:,4)
        Vector2 = Interfaces(:,5) 

        DO i=1,nreg ! Loop over regions
            DO j=1,El_reg(i,1) ! loop over elements per regions    
                element = Subregions(i,j) ! element
                      
                CALL find_vector_g(Vector1,size_v,element,pos1) 
                CALL find_vector_g(Vector2,size_v,element,pos2) 

                IF (pos1.EQ.0 .AND. pos2.EQ.0) THEN    
                        
                    ! Boundary conditions          
                    DO ii=1,nnos_el
                        ! Displacement is known
                        type_bc = INT(SQRT(BC_elem(element,2*ii+2)**2 + BC_elem(element,2*ii+4)**2 + &
                                                                                    BC_elem(element,2*ii+6)**2))
                        disp_val = SQRT(BC_elem(element,2*ii+3)**2 + BC_elem(element,2*ii+5)**2 + BC_elem(element,2*ii+7)**2)
                           
                        IF ((type_bc.EQ.0).AND.(ABS(disp_val).LE.(1e-15))) THEN 
                                
                            Strain(3*element-2:3*element,:) = 0.d0    
                            Stress(3*element-2:3*element,:) = -Stress(3*element-2:3*element,:)

                        END IF
                    END DO
                END IF
            END DO
        END DO

        DEALLOCATE(Vector1,Vector2)

    END IF


END SUBROUTINE constraint_elements
!===============================================================================!
SUBROUTINE strain_t(reg,grad_u_el,Strain_el) 

    INTEGER :: i, reg
    REAL(8) :: grad_u_el(3,9), Strain_el(3,6)

    DO i=1,nnos_el

        Strain_el(i,1) = grad_u_el(i,1) ! e11
        Strain_el(i,2) = grad_u_el(i,2) ! e22
        Strain_el(i,3) = grad_u_el(i,3) ! e33
        Strain_el(i,4) = 0.5d0*(grad_u_el(i,4) + grad_u_el(i,5)) ! e23  
        Strain_el(i,5) = 0.5d0*(grad_u_el(i,6) + grad_u_el(i,7)) ! e13  
        Strain_el(i,6) = 0.5d0*(grad_u_el(i,8) + grad_u_el(i,9)) ! e12

   END DO

END SUBROUTINE strain_t
!===============================================================================!
SUBROUTINE gradient_u(reg,el,n,grad_u_el,sigma_el)

    INTEGER :: el, no1, no2, no3, reg, i
    REAL(8) :: n(3), grad_u_el(3,9), x_el_aux(15,1)

    REAL(8) :: x1, y1, z1, u1, v1, w1, x2, y2, z2, u2, v2, w2
    REAL(8) :: x3, y3, z3, u3, v3, w3
    REAL(8) :: dN1dqsi, dN2dqsi, dN3dqsi, dN1deta, dN2deta, dN3deta
    REAL(8) :: dxdqsi, dydqsi, dzdqsi, dxdeta, dydeta, dzdeta 
    REAL(8) :: dudqsi, dvdqsi, dwdqsi, dudeta, dvdeta, dwdeta
    REAL(8) :: c11, c12, c13, c14, c15, c16, c21, c22, c23, c24, c25, c26
    REAL(8) :: c31, c32, c33, c34, c35, c36, c41, c42, c43, c44, c45, c46
    REAL(8) :: c51, c52, c53, c54, c55, c56, c61, c62, c63, c64, c65, c66
    REAL(8) :: I_C(6,15), N_0(3,15), O_D(6,15), Mat(15,15), t_el(3,3), vec(15,1)
    REAL(8) :: trac(3), sigma_el(3,6)

    INTEGER :: M, N1, LDA, INFO, LWORK
    INTEGER, ALLOCATABLE :: IPIV(:)
    REAL(8), ALLOCATABLE :: WORK(:)
    
    ! Node 1 and coordinates ------------------
    no1 = ELEM_GEO_NODES_global(el,2)
    x1 = PHY_NODES_global(no1,2); 
    y1 = PHY_NODES_global(no1,3);         
    z1 = PHY_NODES_global(no1,4);
    ! Diplacement
    u1 = Displacement(no1,2); ! dx
    v1 = Displacement(no1,3); ! dy
    w1 = Displacement(no1,4); ! dz
    ! Traction
    t_el(1,:) = Traction(no1,2:4); 
    ! Current position
    x1 = x1 + u1; y1 = y1 + v1; z1 = z1 + w1;
                     
    ! Node 2 and coordinates ------------------
    no2 = ELEM_GEO_NODES_global(el,3)   
    x2 = PHY_NODES_global(no2,2); 
    y2 = PHY_NODES_global(no2,3);         
    z2 = PHY_NODES_global(no2,4);
    ! Diplacement
    u2 = Displacement(no2,2); ! dx
    v2 = Displacement(no2,3); ! dy
    w2 = Displacement(no2,4); ! dz
    ! Traction
    t_el(2,:) = Traction(no2,2:4); 
    ! Current position
    x2 = x2 + u2; y2 = y2 + v2; z2 = z2 + w2;

    ! Node 2 and coordinates --------------------
    no3 = ELEM_GEO_NODES_global(el,4)     
    x3 = PHY_NODES_global(no3,2); 
    y3 = PHY_NODES_global(no3,3);         
    z3 = PHY_NODES_global(no3,4); 
    ! Diplacement
    u3 = Displacement(no3,2); ! dx
    v3 = Displacement(no3,3); ! dy
    w3 = Displacement(no3,4); ! dz
    ! Traction
    t_el(3,:) = Traction(no3,2:4); 
    ! Current position
    x3 = x3 + u3; y3 = y3 + v3; z3 = z3 + w3; 
      

    ! Derivatives of the shape functions of discontinous element
    dN1dqsi = dNda_i(1,1);  
    dN2dqsi = dNda_i(1,2); 
    dN3dqsi = dNda_i(1,3);  

    dN1deta = dNda_i(2,1);  
    dN2deta = dNda_i(2,2); 
    dN3deta = dNda_i(2,3); 

    ! First derivative of the current position
    dxdqsi = x1*dN1dqsi +  x2*dN2dqsi + x3*dN3dqsi
    dydqsi = y1*dN1dqsi +  y2*dN2dqsi + y3*dN3dqsi
    dzdqsi = z1*dN1dqsi +  z2*dN2dqsi + z3*dN3dqsi

    dxdeta = x1*dN1deta +  x2*dN2deta + x3*dN3deta
    dydeta = y1*dN1deta +  y2*dN2deta + y3*dN3deta
    dzdeta = z1*dN1deta +  z2*dN2deta + z3*dN3deta

    ! First derivative of the displacement respect to psi and eta
    dudqsi = u1*dN1dqsi +  u2*dN2dqsi + u3*dN3dqsi
    dvdqsi = v1*dN1dqsi +  v2*dN2dqsi + v3*dN3dqsi
    dwdqsi = w1*dN1dqsi +  w2*dN2dqsi + w3*dN3dqsi

    dudeta = u1*dN1deta +  u2*dN2deta + u3*dN3deta
    dvdeta = v1*dN1deta +  v2*dN2deta + v3*dN3deta
    dwdeta = w1*dN1deta +  w2*dN2deta + w3*dN3deta
    
    ! Stiffness coefficients
    c11 = C_tot(reg)%C_reg(1,1); c12 = C_tot(reg)%C_reg(1,2);
    c13 = C_tot(reg)%C_reg(1,3); c14 = C_tot(reg)%C_reg(1,4);
    c15 = C_tot(reg)%C_reg(1,5); c16 = C_tot(reg)%C_reg(1,6);
    c21 = C_tot(reg)%C_reg(2,1); c22 = C_tot(reg)%C_reg(2,2);
    c23 = C_tot(reg)%C_reg(2,3); c24 = C_tot(reg)%C_reg(2,4);
    c25 = C_tot(reg)%C_reg(2,5); c26 = C_tot(reg)%C_reg(2,6);
    c31 = C_tot(reg)%C_reg(3,1); c32 = C_tot(reg)%C_reg(3,2);
    c33 = C_tot(reg)%C_reg(3,3); c34 = C_tot(reg)%C_reg(3,4);
    c35 = C_tot(reg)%C_reg(3,5); c36 = C_tot(reg)%C_reg(3,6);
    c41 = C_tot(reg)%C_reg(4,1); c42 = C_tot(reg)%C_reg(4,2);
    c43 = C_tot(reg)%C_reg(4,3); c44 = C_tot(reg)%C_reg(4,4);
    c45 = C_tot(reg)%C_reg(4,5); c46 = C_tot(reg)%C_reg(4,6);
    c51 = C_tot(reg)%C_reg(5,1); c52 = C_tot(reg)%C_reg(5,2);
    c53 = C_tot(reg)%C_reg(5,3); c54 = C_tot(reg)%C_reg(5,4);
    c55 = C_tot(reg)%C_reg(5,5); c56 = C_tot(reg)%C_reg(5,6);
    c61 = C_tot(reg)%C_reg(6,1); c62 = C_tot(reg)%C_reg(6,2);
    c63 = C_tot(reg)%C_reg(6,3); c64 = C_tot(reg)%C_reg(6,4);
    c65 = C_tot(reg)%C_reg(6,5); c66 = C_tot(reg)%C_reg(6,6);
   
    ! Matriz [[I][C]]{x} = {0}
    I_C(1,:) = (/1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-c11,-c12,-c13,-0.5d0*c14,-0.5d0*c14,-0.5d0*c15,-0.5d0*c15,-0.5d0*c16,-0.5d0*c16/)
    I_C(2,:) = (/0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,-c21,-c22,-c23,-0.5d0*c24,-0.5d0*c24,-0.5d0*c25,-0.5d0*c25,-0.5d0*c26,-0.5d0*c26/)
    I_C(3,:) = (/0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,-c31,-c32,-c33,-0.5d0*c34,-0.5d0*c34,-0.5d0*c35,-0.5d0*c35,-0.5d0*c36,-0.5d0*c36/)
    I_C(4,:) = (/0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,-c41,-c42,-c43,-0.5d0*c44,-0.5d0*c44,-0.5d0*c45,-0.5d0*c45,-0.5d0*c46,-0.5d0*c46/)
    I_C(5,:) = (/0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,-c51,-c52,-c53,-0.5d0*c54,-0.5d0*c54,-0.5d0*c55,-0.5d0*c55,-0.5d0*c56,-0.5d0*c56/)
    I_C(6,:) = (/0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,-c61,-c62,-c63,-0.5d0*c64,-0.5d0*c64,-0.5d0*c65,-0.5d0*c65,-0.5d0*c66,-0.5d0*c66/)      

    ! Matrix [[N][0]]{x} = {t}
    N_0(1,:) = (/n(1),0.d0,0.d0,0.d0,n(3),n(2),0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
    N_0(2,:) = (/0.d0,n(2),0.d0,n(3),0.d0,n(1),0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
    N_0(3,:) = (/0.d0,0.d0,n(3),n(2),n(1),0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)

    ! Matrix [[0][D]]{x} = {grad_a u}
    O_D(1,:) = (/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,dxdqsi,0.d0,0.d0,0.d0,0.d0,dzdqsi,0.d0,dydqsi,0.d0/) 
    O_D(2,:) = (/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,dydqsi,0.d0,dzdqsi,0.d0,0.d0,0.d0,0.d0,dxdqsi/)
    O_D(3,:) = (/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,dzdqsi,0.d0,dydqsi,0.d0,dxdqsi,0.d0,0.d0/)
    O_D(4,:) = (/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,dxdeta,0.d0,0.d0,0.d0,0.d0,dzdeta,0.d0,dydeta,0.d0/) 
    O_D(5,:) = (/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,dydeta,0.d0,dzdeta,0.d0,0.d0,0.d0,0.d0,dxdeta/)
    O_D(6,:) = (/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,dzdeta,0.d0,dydeta,0.d0,dxdeta,0.d0,0.d0/)

    ! Full matrix
    Mat(1:6,1:15) = I_C
    Mat(7:9,1:15) = N_0
    Mat(10:15,1:15) = O_D
    
    ! Inverse
    M = 15
    N1 = 15
    LDA = M
    ALLOCATE(IPIV(M))
    IPIV = 0
    INFO = 0
    CALL DGETRF(M, N1, Mat, LDA, IPIV, INFO)
    LWORK = N1
    ALLOCATE(WORK(LWORK))
    WORK = 0.d0
    CALL DGETRI(N1, Mat, LDA, IPIV, WORK, LWORK, INFO)
    DEALLOCATE(IPIV,WORK)
    
    DO i=1,nnos_el

        trac = t_el(i,:)

        ! vector {b} = [{0}^T {t}^T {grad_a u}^T]
        vec(:,1) = (/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,trac(1),trac(2),trac(3), &
                                          dudqsi,dvdqsi,dwdqsi,dudeta,dvdeta,dwdeta/)    
        
        x_el_aux = MATMUL(Mat,vec)

        sigma_el(i,:) = x_el_aux(1:6,1)
        grad_u_el(i,:) = x_el_aux(7:15,1)
          
    END DO

END SUBROUTINE gradient_u
!===============================================================================!
END MODULE Response_S_S
