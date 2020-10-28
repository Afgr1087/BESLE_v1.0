!===============================================================================!
!-------------------------------------------------------------------------------!
!            MODULE TO COMPUTE THE FOURIER COEFFICIENT FOR THE                  !
!                            FUNDAMENTAL SOLUTIONS                              !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Fourier_Coefficients
!-------------------------------------------------------------------------------!
USE Global_variables
USE Global_functions
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!
    SUBROUTINE Fourier_Coeff(me,Test_name,reg,RIt_cell, RIc_cell, RI0m_cell, &
                                        RIm0_cell,R00,pos_vet, size_coef_cells, &
                                        C,S,Anglesv)

        INTEGER :: reg, me
        REAL(8) :: max_val
        CHARACTER(LEN=20)::Test_name

        INTEGER :: alpha, size_coef_cells(4)

        REAL(KIND=8), ALLOCATABLE :: Rt_mat(:,:), Rc_mat(:,:)
        REAL(KIND=8), ALLOCATABLE :: It_mat(:,:), Ic_mat(:,:)
        REAL(KIND=8), ALLOCATABLE :: R0m_vet(:,:), Rm0_vet(:,:)
        REAL(KIND=8), ALLOCATABLE :: I0m_vet(:,:), Im0_vet(:,:)
        REAL(KIND=8):: R00(3,3), Anglesv(nreg,4)

        TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIt_cell
        TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIc_cell
        TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RI0m_cell
        TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIm0_cell
        INTEGER, ALLOCATABLE :: pos_vet(:)

        REAL(KIND=8) :: C(6,6) ! Compliance tensor
        REAL(KIND=8) :: S(6,6) ! Stiffness tensor
                              
	

        CALL Read_Fourier_Coeff(me,max_val,Test_name,reg,Rt_mat,Rc_mat, &
                                It_mat,Ic_mat,R0m_vet,I0m_vet,Rm0_vet,Im0_vet, &
                                R00,alpha,C,S,Anglesv)
		
        CALL RI_mat(max_val,Rt_mat,Rc_mat,It_mat,Ic_mat,R0m_vet,I0m_vet,Rm0_vet, &
                    Im0_vet,RIt_cell, RIc_cell, RI0m_cell,RIm0_cell, pos_vet, &
                    alpha, size_coef_cells)

        C_tot(reg)%C_reg = C
                                    

    END SUBROUTINE Fourier_Coeff
!===============================================================================!
    SUBROUTINE Read_Fourier_Coeff(me,max_val,Test_name,reg,Rt_mat,Rc_mat, &
                                 It_mat,Ic_mat,R0m_vet,I0m_vet,Rm0_vet,Im0_vet, &
                                 R00,alpha,C,S,Anglesv)

        REAL(8) :: max_val
        INTEGER :: i, reg,me
        CHARACTER(LEN=20)::num_str, num_str_tot
        CHARACTER(LEN=20)::Test_name

        INTEGER :: M, N, LDA, INFO, LWORK
        INTEGER, ALLOCATABLE :: IPIV(:)
        REAL(8), ALLOCATABLE :: WORK(:)

        INTEGER :: alpha

        REAL(KIND=8), ALLOCATABLE :: Rt_mat(:,:), Rc_mat(:,:)
        REAL(KIND=8), ALLOCATABLE :: It_mat(:,:), Ic_mat(:,:)
        REAL(KIND=8), ALLOCATABLE :: R0m_vet(:,:), Rm0_vet(:,:)
        REAL(KIND=8), ALLOCATABLE :: I0m_vet(:,:), Im0_vet(:,:)
        REAL(KIND=8):: R00(3,3),Anglesv(nreg,4)

        REAL(KIND=8) :: C(6,6) ! Compliance tensor
        REAL(KIND=8) :: S(6,6) ! Stiffness tensor
              
        IF (nreg.EQ.1) THEN ! for one region

            OPEN(2,file=trim(fileplace_material)//trim(Material_coefficients_file)//'.dat',STATUS='OLD')
            IF (me.EQ.0) THEN   
                WRITE (*,*) '    ',trim(fileplace_material),trim(Material_coefficients_file)
            END IF

        ELSE ! for multiregion
         
            ! Isotropic or Anisotropic - > same properties for each region
            IF (n_Materials.EQ.0) THEN

                OPEN(2,file=trim(fileplace_material)//trim(Material_coefficients_file)//'.dat',STATUS='OLD')
                IF (me.EQ.0) THEN   
                    WRITE (*,*) '    ',trim(fileplace_material),trim(Material_coefficients_file),'.dat'
                END IF
            END IF

            ! Isotropic or Anisotropic - > random distrubuition properties
            IF (n_Materials.GE.1) THEN

                OPEN(2,file=trim(fileplace_material)//trim(Material_coefficients_database)//'/'//trim(Test_name),STATUS='OLD') 
                IF (me.EQ.0) THEN  
                    WRITE (*,*) '    ',trim(fileplace_material),trim(Material_coefficients_database),'/',trim(Test_name)
                END IF
            END IF    

        END IF

        REWIND(2)
      
        READ (2,'(10x,I12)') alpha
        
        READ (2,'(3E30.17)') Anglesv(reg,2:4)
        
        WRITE(num_str,"(I10)") alpha*3
        num_str_tot = "("//TRIM(ADJUSTL(num_str))//"E30.17)"
        
        READ (2,num_str_tot) max_val
                
        ALLOCATE(Rt_mat(alpha*3,alpha*3))  
        DO i = 1,alpha*3
            READ (2,num_str_tot) Rt_mat(i,1:alpha*3)
        END DO
        
        ALLOCATE(Rc_mat(alpha*3,alpha*3))  
        DO i = 1,alpha*3
            READ (2,num_str_tot) Rc_mat(i,1:alpha*3)
        END DO

        ALLOCATE(It_mat(alpha*3,alpha*3))  
        DO i = 1,alpha*3
            READ (2,num_str_tot) It_mat(i,1:alpha*3)
        END DO

        ALLOCATE(Ic_mat(alpha*3,alpha*3))  
        DO i = 1,alpha*3
            READ (2,num_str_tot) Ic_mat(i,:)
        END DO
        
        ALLOCATE(R0m_vet(alpha*3,3))  
        DO i = 1,alpha*3
            READ (2,'(3E30.17)') R0m_vet(i,:)
        END DO

        ALLOCATE(I0m_vet(alpha*3,3))  
        DO i = 1,alpha*3
            READ (2,'(3E30.17)') I0m_vet(i,:)
        END DO

        ALLOCATE(Rm0_vet(alpha*3,3))  
        DO i = 1,alpha*3
            READ (2,'(3E30.17)') Rm0_vet(i,:)
        END DO

        ALLOCATE(Im0_vet(alpha*3,3))
        DO i = 1,alpha*3
            READ (2,'(3E30.17)') Im0_vet(i,:)
        END DO

        DO i = 1,3
            READ (2,'(3E30.17)') R00(i,:)
        END DO

        DO i = 1,6
            READ (2,'(6F30.16)') C(i,:)
        END DO
        
        ! Stiffness tensor [S] = C^-1
        C = C*scale_prop_mat
        S = C
        ! compute an LU factorization of a general M-by-N matrix F_reg
        
        M = SIZE(C,1)
        N = SIZE(C,2)
        LDA = M
        ALLOCATE(IPIV(M))
        IPIV = 0
        INFO = 0
        CALL DGETRF(M, N, S, LDA, IPIV, INFO)
        
        ! compute the inverse of a matrix using the LU	factorization
        LWORK = N
        ALLOCATE(WORK(LWORK))
        WORK = 0.d0
        CALL DGETRI(N, S, LDA, IPIV, WORK, LWORK, INFO)
            
                             
    END SUBROUTINE Read_Fourier_Coeff
!===============================================================================!
    SUBROUTINE RI_mat(max_val,Rt_mat,Rc_mat,It_mat,Ic_mat,R0m_vet,I0m_vet,Rm0_vet, &
                      Im0_vet,RIt_cell, RIc_cell, RI0m_cell,RIm0_cell, pos_vet, &
                      alpha, size_coef_cells)
        
        REAL(8) :: max_val, a1, a2, b1, b2, c1, c2, d1, d2, Mat(3,3)
        INTEGER :: m, n, ca, cb, cc, cd

        INTEGER :: alpha, size_coef_cells(4)

        REAL(KIND=8), ALLOCATABLE :: Rt_mat(:,:), Rc_mat(:,:)
        REAL(KIND=8), ALLOCATABLE :: It_mat(:,:), Ic_mat(:,:)
        REAL(KIND=8), ALLOCATABLE :: R0m_vet(:,:), Rm0_vet(:,:)
        REAL(KIND=8), ALLOCATABLE :: I0m_vet(:,:), Im0_vet(:,:)

        TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIt_cell
        TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIc_cell
        TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RI0m_cell
        TYPE(Cell2), DIMENSION(:), ALLOCATABLE :: RIm0_cell
        INTEGER, ALLOCATABLE :: pos_vet(:)
      

        ALLOCATE(RIt_cell(alpha*alpha), RIc_cell(alpha*alpha))
        ALLOCATE(RI0m_cell(alpha*alpha), RIm0_cell(alpha*alpha))
        ALLOCATE(pos_vet(alpha*alpha))
        pos_vet = 0

       ca=0; cb=0; cc=0; cd=0;
       DO m=1,alpha
            DO n=1,alpha
                Mat = Rt_mat(3*m-2:3*m,3*n-2:3*n)
                a1 = E_NORM(Mat)
                Mat = It_mat(3*m-2:3*m,3*n-2:3*n)
                a2 = E_NORM(Mat)
                Mat = Rc_mat(3*m-2:3*m,3*n-2:3*n)
                b1 = E_NORM(Mat)
                Mat = Ic_mat(3*m-2:3*m,3*n-2:3*n)
                b2 = E_NORM(Mat)


                IF ((a1.GT.error*max_val).OR.(a2.GT.error*max_val)) THEN
                    ca = ca + 1
                    IF (a1.GT.error*max_val) THEN
                        RIt_cell(ca)%re = Rt_mat(3*m-2:3*m,3*n-2:3*n)
                    ELSE
                        RIt_cell(ca)%re = 0.d0
                    END IF
                    IF (a2.GT.error*max_val) THEN
                        RIt_cell(ca)%im = It_mat(3*m-2:3*m,3*n-2:3*n)
                    ELSE
                        RIt_cell(ca)%im = 0.d0
                    END IF
                    RIt_cell(ca)%pos_mat = (/m,n/)
                END IF

                IF ((b1.GT.error*max_val).OR.(b2.GT.error*max_val)) THEN
                    cb = cb + 1
                    IF (b1.GT.error*max_val) THEN
                        RIc_cell(cb)%re = Rc_mat(3*m-2:3*m,3*n-2:3*n)
                    ELSE
                        RIc_cell(cb)%re = 0.d0
                    END IF
                    IF (b2.GT.error*max_val) THEN
                        RIc_cell(cb)%im = Ic_mat(3*m-2:3*m,3*n-2:3*n)
                    ELSE
                        RIc_cell(cb)%im = 0.d0
                    END IF
                    RIc_cell(cb)%pos_mat = (/m,n/)
                END IF

            END DO

            Mat = R0m_vet(3*m-2:3*m,:)
            c1 = E_NORM(Mat)
            Mat = I0m_vet(3*m-2:3*m,:)
            c2 = E_NORM(Mat)
            Mat = Rm0_vet(3*m-2:3*m,:)
            d1 = E_NORM(Mat)
            Mat = Im0_vet(3*m-2:3*m,:)
            d2 = E_NORM(Mat)

            IF ((c1.GT.error*max_val).OR.(c2.GT.error*max_val)) THEN
                cc = cc + 1
                IF (c1.GT.error*max_val) THEN
                    RI0m_cell(cc)%re = R0m_vet(3*m-2:3*m,:)
                ELSE
                    RI0m_cell(cc)%re = 0.d0
                END IF
                IF (c2.GT.error*max_val) THEN
                    RI0m_cell(cc)%im = I0m_vet(3*m-2:3*m,:)
                ELSE
                    RI0m_cell(cc)%im = 0.d0
                END IF
                pos_vet(cc) = m
            END IF

            IF ((d1.GT.error*max_val).OR.(d2.GT.error*max_val)) THEN
                cd = cd + 1
                IF (d1.GT.error*max_val) THEN
                    RIm0_cell(cd)%re = Rm0_vet(3*m-2:3*m,:)
                ELSE
                    RIm0_cell(cd)%re = 0.d0
                END IF
                IF (d2.GT.error*max_val) THEN
                    RIm0_cell(cd)%im = Im0_vet(3*m-2:3*m,:)
                ELSE
                    RIm0_cell(cd)%im = 0.d0
                END IF
                pos_vet(cd) = m
            END IF

       END DO
       size_coef_cells = (/ca,cb,cc,cd/)

       DEALLOCATE(Rt_mat, Rc_mat, It_mat, Ic_mat)
       DEALLOCATE(R0m_vet, Rm0_vet, I0m_vet, Im0_vet)
    
    END SUBROUTINE RI_mat
!===============================================================================!
FUNCTION E_NORM(Mat)

    INTEGER :: m, n , i, j
    REAL(8) :: E_NORM, val, s
    REAL(8) :: Mat(3,3)
    
    s = 0.d0
    m = SIZE(Mat,1); n = SIZE(Mat,2)
    DO i=1,m
        DO j=1,n
            val = Mat(i,j)**2
            s = s + val
        END DO
    END DO
    E_NORM = DSQRT(s)
        
END FUNCTION E_NORM
!===============================================================================!
END MODULE Fourier_Coefficients

