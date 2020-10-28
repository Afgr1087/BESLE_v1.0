!===============================================================================!
!-------------------------------------------------------------------------------!
!                      MODULE TO APPLY HOMOGENIZATION                           !
!
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Homogenization
!-------------------------------------------------------------------------------!
USE Global_variables
USE Global_functions
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
SUBROUTINE Homogenized_sigma_epsilon(ts1,step)

    INTEGER::t0,t1,rate, step
    REAL(8),INTENT(OUT)::ts1

    INTEGER :: i, n_el, j, el!, node, pos
    REAL(8) :: n(3), Volume, H_e_el (3,3), H_strain(3,3)
    REAL(8) :: H_strain_reg(3,3), H_s_el (3,3), H_stress_reg(3,3), H_stress(3,3)
    REAL(8) :: strain(6), sigma(6), max_x, max_y, max_z   

    CALL SYSTEM_CLOCK(t0, rate) ! Initializates the time counter
    WRITE (*,*) '----------------------------------------------------------------------------'
    WRITE (*,*) ' '
    WRITE (*,'(A)',advance='no') ' 11. Homogenization               '

    max_x = MAXVAL(POINTS(:,1))
    max_y = MAXVAL(POINTS(:,2))
    max_z = MAXVAL(POINTS(:,3))
    
    Volume = max_x*max_y*max_z
   
    H_strain = 0.d0
    H_stress = 0.d0

    DO i = 1,nreg ! Loop over regions
        n_el = El_reg(i,1) ! Elements per region
		
        !Volume = Volumes(i) ! Grain volume
        !Volume = scale_size**3
        H_strain_reg = 0.d0 ! Homogenized strain tensor / reg
        H_stress_reg = 0.d0 ! Homogenized stress tensor / reg
        DO j = 1,n_el 
            el = Subregions(i,j) 
            IF (el.NE.0) THEN 

                ! Normal vector of the deformed element
				n = NORMAL_VEC_Def(el,:)

                ! Homogenized strain tensor per element
                CALL Homogenized_strain_tensor(el,n,H_e_el)      
                !H_strain_reg = H_strain_reg + H_e_el

                ! Homogenized stress tensor per element
                CALL Homogenized_stress_tensor(el,H_s_el) 
                !H_stress_reg = H_stress_reg + H_s_el

				!IF (Damage_zones(el).EQ.0) THEN
				
					H_strain_reg = H_strain_reg + H_e_el
					H_stress_reg = H_stress_reg + H_s_el			
		
				!END IF

            END IF 
        END DO 

        ! Homogenized strain tensor per grain
        H_strain_reg = (1.d0/(2.d0*Volume))*H_strain_reg
        ! Homogenized strain tensor total
        H_strain = H_strain + H_strain_reg    

        ! Homogenized stress tensor per grain
        H_stress_reg = (1.d0/(2.d0*Volume))*H_stress_reg
        ! Homogenized stress tensor total
        H_stress = H_stress + H_stress_reg  

    END DO
    
    WRITE(*,*) ''
    WRITE(*,*) ''
    WRITE(*,*) '     -- Homogenized strain tensor'
    WRITE(*,*) ''
    DO i=1,3
        WRITE(*,'(A,3E15.4)') '        ', H_strain(i,1), H_strain(i,2), H_strain(i,3)
    END DO
    
    WRITE(*,*) ''
    WRITE(*,*) '     -- Homogenized stress tensor'
    WRITE(*,*) ''	
    DO i=1,3
        WRITE(*,'(A,3F15.4)') '        ', H_stress(i,1)/scale_prop_mat, H_stress(i,2)/scale_prop_mat, H_stress(i,3)/scale_prop_mat
    END DO
    WRITE(*,*) ''
    !Stress_h(step) = H_stress(3,3)
    
    !Strain_h(step) = H_strain(3,3)

    CALL SYSTEM_CLOCK(t1)       ! Time counter stopped
    WRITE (*,'(A,F15.3,A)') '                                  ...COMPLETED! (Time :', REAL(t1-t0)/rate,'s)'
    WRITE (*,*) ''
    ts1 = REAL(t1-t0)/rate

END SUBROUTINE Homogenized_sigma_epsilon
!===============================================================================!
SUBROUTINE Homogenized_stress_tensor(el,s_el)

    INTEGER :: el, no1, no2, no3, i, npts
    REAL(8) :: t1(3), t2(3), t3(3), X1(3), X2(3), X3(3), u1(3), u2(3), u3(3)
    
    REAL(8) :: wgts(npoints_Nsi), N1(npoints_Nsi), N2(npoints_Nsi), N3(npoints_Nsi)
    REAL(8) :: dN1dqsi(npoints_Nsi), dN2dqsi(npoints_Nsi), dN3dqsi(npoints_Nsi)
    REAL(8) :: dN1deta(npoints_Nsi), dN2deta(npoints_Nsi), dN3deta(npoints_Nsi)
    REAL(8) :: tx(npoints_Nsi), ty(npoints_Nsi), tz(npoints_Nsi) 
    REAL(8) :: Xx(npoints_Nsi), Xy(npoints_Nsi), Xz(npoints_Nsi)
    REAL(8) :: dxdqsi(npoints_Nsi), dydqsi(npoints_Nsi), dzdqsi(npoints_Nsi)
    REAL(8) :: dxdeta(npoints_Nsi), dydeta(npoints_Nsi), dzdeta(npoints_Nsi)
    REAL(8) :: g1(npoints_Nsi), g2(npoints_Nsi), g3(npoints_Nsi), Jac(npoints_Nsi)
    REAL(8) :: t(3), Xc(3)
    REAL(8) :: sigma11(npoints_Nsi), sigma12(npoints_Nsi), sigma13(npoints_Nsi)
    REAL(8) :: sigma21(npoints_Nsi), sigma22(npoints_Nsi), sigma23(npoints_Nsi)
    REAL(8) :: sigma31(npoints_Nsi), sigma32(npoints_Nsi), sigma33(npoints_Nsi)
    REAL(8) :: s11, s12, s13, s21, s22, s23, s31, s32, s33, s_el(3,3)

    ! Node 1 ---------------------------------
    no1 = ELEM_GEO_NODES_global(el,2)
    ! Traction and coordinates 
    t1 = Traction(no1,2:4); 
    X1 = PHY_NODES_global(no1,2:4);
    u1 = Displacement(no1,2:4);
    X1 = X1 + u1            
    ! Node 2 ---------------------------------
    no2 = ELEM_GEO_NODES_global(el,3) 
    ! Traction and coordinates 
    t2 = Traction(no2,2:4);  
    X2 = PHY_NODES_global(no2,2:4); 
    u2 = Displacement(no2,2:4);
    X2 = X2 + u2    
    ! Node 3 ---------------------------------
    no3 = ELEM_GEO_NODES_global(el,4) 
    ! Traction and coordinates 
    t3 = Traction(no3,2:4);   
    X3 = PHY_NODES_global(no3,2:4);     
    u3 = Displacement(no3,2:4);
    X3 = X3 + u3 
    ! Weight for numeric integration 
    wgts = weights_ns(:,12)
            
    ! Shape functions of discontinuous element
    N1 = weights_ns(:,3); N2 = weights_ns(:,4); N3 = weights_ns(:,5);
            
    ! Derivatives of the shape functions of discontinous element
    dN1dqsi = weights_ns(:,6);  
    dN2dqsi = weights_ns(:,7);
    dN3dqsi = weights_ns(:,8); 

    dN1deta = weights_ns(:,9); 
    dN2deta = weights_ns(:,10)
    dN3deta = weights_ns(:,11);

    ! Coordinates of the points in the triangular element 
    Xx = N1*X1(1) + N2*X2(1) + N3*X3(1)
    Xy = N1*X1(2) + N2*X2(2) + N3*X3(2) 
    Xz = N1*X1(3) + N2*X2(3) + N3*X3(3) 

    ! Value for traction in each coordinate
    tx = N1*t1(1) + N2*t2(1) + N3*t3(1)
    ty = N1*t1(2) + N2*t2(2) + N3*t3(2) 
    tz = N1*t1(3) + N2*t2(3) + N3*t3(3) 

    ! Jacobian
    dxdqsi = X1(1)*dN1dqsi + X2(1)*dN2dqsi + X3(1)*dN3dqsi 
    dydqsi = X1(2)*dN1dqsi + X2(2)*dN2dqsi + X3(2)*dN3dqsi 
    dzdqsi = X1(3)*dN1dqsi + X2(3)*dN2dqsi + X3(3)*dN3dqsi 

    dxdeta = X1(1)*dN1deta + X2(1)*dN2deta + X3(1)*dN3deta 
    dydeta = X1(2)*dN1deta + X2(2)*dN2deta + X3(2)*dN3deta 
    dzdeta = X1(3)*dN1deta + X2(3)*dN2deta + X3(3)*dN3deta 

    g1 = dydqsi*dzdeta - dzdqsi*dydeta
    g2 = dzdqsi*dxdeta - dxdqsi*dzdeta
    g3 = dxdqsi*dydeta - dydqsi*dxdeta

    Jac = DSQRT(g1**2 + g2**2 + g3**2)

    npts = npoints_Nsi

    sigma11 = 0.d0; sigma12 = 0.d0; sigma13 = 0.d0
    sigma21 = 0.d0; sigma22 = 0.d0; sigma23 = 0.d0
    sigma31 = 0.d0; sigma32 = 0.d0; sigma33 = 0.d0

    DO i = 1,npts ! Over all points in the element

        t = (/tx(i), ty(i), tz(i)/)
        Xc = (/Xx(i), Xy(i), Xz(i)/)
        
        sigma11(i) = t(1)*Xc(1) + t(1)*Xc(1)
        sigma12(i) = t(1)*Xc(2) + t(2)*Xc(1)
        sigma13(i) = t(1)*Xc(3) + t(3)*Xc(1)

        sigma21(i) = t(2)*Xc(1) + t(1)*Xc(2)
        sigma22(i) = t(2)*Xc(2) + t(2)*Xc(2)
        sigma23(i) = t(2)*Xc(3) + t(3)*Xc(2)

        sigma31(i) = t(3)*Xc(1) + t(1)*Xc(3)
        sigma32(i) = t(3)*Xc(2) + t(2)*Xc(3)
        sigma33(i) = t(3)*Xc(3) + t(3)*Xc(3)
   
    END DO
        
    ! Components of stress tensor
    s11 = SUM(sigma11*wgts*Jac);
    s12 = SUM(sigma12*wgts*Jac);
    s13 = SUM(sigma13*wgts*Jac); 
    s21 = SUM(sigma21*wgts*Jac); 
    s22 = SUM(sigma22*wgts*Jac);
    s23 = SUM(sigma23*wgts*Jac);
    s31 = SUM(sigma31*wgts*Jac);
    s32 = SUM(sigma32*wgts*Jac);
    s33 = SUM(sigma33*wgts*Jac);
    
    ! Stress tensor per element
    s_el(1,:) = (/s11, s12, s13/)
    s_el(2,:) = (/s21, s22, s23/)
    s_el(3,:) = (/s31, s32, s33/)

END SUBROUTINE Homogenized_stress_tensor
!===============================================================================!
SUBROUTINE Homogenized_strain_tensor(el,n,e_el)

    INTEGER :: el, no1, no2, no3, i, npts
    REAL(8) :: n(3), u1(3), u2(3), u3(3), X1(3), X2(3), X3(3)

    REAL(8) :: wgts(npoints_Nsi), N1(npoints_Nsi), N2(npoints_Nsi), N3(npoints_Nsi)
    REAL(8) :: dN1dqsi(npoints_Nsi), dN2dqsi(npoints_Nsi), dN3dqsi(npoints_Nsi)
    REAL(8) :: dN1deta(npoints_Nsi), dN2deta(npoints_Nsi), dN3deta(npoints_Nsi)
    REAL(8) :: ux(npoints_Nsi), uy(npoints_Nsi), uz(npoints_Nsi)  
    REAL(8) :: dxdqsi(npoints_Nsi), dydqsi(npoints_Nsi), dzdqsi(npoints_Nsi)
    REAL(8) :: dxdeta(npoints_Nsi), dydeta(npoints_Nsi), dzdeta(npoints_Nsi)
    REAL(8) :: g1(npoints_Nsi), g2(npoints_Nsi), g3(npoints_Nsi), Jac(npoints_Nsi)
    REAL(8) :: u(3)
    REAL(8) :: epsilon11(npoints_Nsi), epsilon12(npoints_Nsi), epsilon13(npoints_Nsi)
    REAL(8) :: epsilon21(npoints_Nsi), epsilon22(npoints_Nsi), epsilon23(npoints_Nsi)
    REAL(8) :: epsilon31(npoints_Nsi), epsilon32(npoints_Nsi), epsilon33(npoints_Nsi)
    REAL(8) :: e11, e12, e13, e21, e22, e23, e31, e32, e33, e_el(3,3)

    ! Node 1 ---------------------------------
    no1 = ELEM_GEO_NODES_global(el,2)
    ! Diplacement and coordinates 
    u1 = Displacement(no1,2:4); 
    X1 = PHY_NODES_global(no1,2:4);
    X1 = X1 + u1         
    ! Node 2 ---------------------------------
    no2 = ELEM_GEO_NODES_global(el,3) 
    ! Diplacement and coordinates 
    u2 = Displacement(no2,2:4);  
    X2 = PHY_NODES_global(no2,2:4);
    X2 = X2 + u2
    ! Node 3 ---------------------------------
    no3 = ELEM_GEO_NODES_global(el,4) 
    ! Diplacement and coordinates 
    u3 = Displacement(no3,2:4);   
    X3 = PHY_NODES_global(no3,2:4);
    X3 = X3 + u3
    ! Weight for numeric integration 
    wgts = weights_ns(:,12)
            
    ! Shape functions of discontinuous element
    N1 = weights_ns(:,3); N2 = weights_ns(:,4); N3 = weights_ns(:,5);
            
    ! Derivatives of the shape functions of discontinous element
    dN1dqsi = weights_ns(:,6);  
    dN2dqsi = weights_ns(:,7);
    dN3dqsi = weights_ns(:,8); 

    dN1deta = weights_ns(:,9); 
    dN2deta = weights_ns(:,10)
    dN3deta = weights_ns(:,11);

    ! Value for displacement in each coordinate
    ux = N1*u1(1) + N2*u2(1) + N3*u3(1)
    uy = N1*u1(2) + N2*u2(2) + N3*u3(2) 
    uz = N1*u1(3) + N2*u2(3) + N3*u3(3) 

    ! Jacobian
    dxdqsi = X1(1)*dN1dqsi + X2(1)*dN2dqsi + X3(1)*dN3dqsi 
    dydqsi = X1(2)*dN1dqsi + X2(2)*dN2dqsi + X3(2)*dN3dqsi 
    dzdqsi = X1(3)*dN1dqsi + X2(3)*dN2dqsi + X3(3)*dN3dqsi 

    dxdeta = X1(1)*dN1deta + X2(1)*dN2deta + X3(1)*dN3deta 
    dydeta = X1(2)*dN1deta + X2(2)*dN2deta + X3(2)*dN3deta 
    dzdeta = X1(3)*dN1deta + X2(3)*dN2deta + X3(3)*dN3deta 

    g1 = dydqsi*dzdeta - dzdqsi*dydeta
    g2 = dzdqsi*dxdeta - dxdqsi*dzdeta
    g3 = dxdqsi*dydeta - dydqsi*dxdeta

    Jac = DSQRT(g1**2 + g2**2 + g3**2)

    npts = npoints_Nsi

    DO i = 1,npts ! Over all points in the element

        u = (/ux(i), uy(i), uz(i)/)

        epsilon11(i) = u(1)*n(1) + n(1)*u(1); 
        epsilon12(i) = u(1)*n(2) + n(1)*u(2); 
        epsilon13(i) = u(1)*n(3) + n(1)*u(3); 
        epsilon21(i) = u(2)*n(1) + n(2)*u(1); 
        epsilon22(i) = u(2)*n(2) + n(2)*u(2); 
        epsilon23(i) = u(2)*n(3) + n(2)*u(3); 
        epsilon31(i) = u(3)*n(1) + n(3)*u(1); 
        epsilon32(i) = u(3)*n(2) + n(3)*u(2); 
        epsilon33(i) = u(3)*n(3) + n(3)*u(3); 

    END DO

    ! Components of strain tensor
    e11 = SUM(epsilon11*wgts*Jac);
    e12 = SUM(epsilon12*wgts*Jac);
    e13 = SUM(epsilon13*wgts*Jac);
    e21 = SUM(epsilon21*wgts*Jac);
    e22 = SUM(epsilon22*wgts*Jac);
    e23 = SUM(epsilon23*wgts*Jac);
    e31 = SUM(epsilon31*wgts*Jac);
    e32 = SUM(epsilon32*wgts*Jac);
    e33 = SUM(epsilon33*wgts*Jac);

    ! Strain tensor per element
    e_el(1,:) = (/e11, e12, e13/)
    e_el(2,:) = (/e21, e22, e23/)
    e_el(3,:) = (/e31, e32, e33/)

END SUBROUTINE Homogenized_strain_tensor
!===============================================================================!

END MODULE Homogenization
