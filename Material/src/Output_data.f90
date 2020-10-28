!===============================================================================!
!-------------------------------------------------------------------------------!
!                               MODULE FOR OUTPUT DATA                          !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Output_data
USE Global_variables 
CONTAINS
!-------------------------------------------------------------------------------!
    SUBROUTINE Print_data(ts1)
!-------------------------------------------------------------------------------!   
        INTEGER::t0,t1,rate
        REAL,INTENT(OUT)::ts1
        INTEGER :: alpha=16, i
        CHARACTER(LEN=20)::num_str, num_str_tot
        
        CALL SYSTEM_CLOCK(t0, rate) ! Initializates the time counter
        
        WRITE (*,*) '----------------------------------------------------------------------'
        WRITE (*,*) ' '
        WRITE (*,'(A)',advance='no') ' 2. Writing output file                  ...'
                
        OPEN(2,file=trim(fileplace)//trim(file_name),STATUS='REPLACE')
        REWIND(2)
        
        WRITE(2,*) alpha

        WRITE(2,'(3E30.17)') phi_1, phi, phi_2

        WRITE(num_str,"(I10)") alpha*3

        num_str_tot = "("//TRIM(ADJUSTL(num_str))//"E30.17)"
          
        WRITE(2,num_str_tot) max_val
               
        DO i=1,alpha*3
            WRITE(2,num_str_tot) Rt_matrix(i,:) 
        END DO

        DO i=1,alpha*3
            WRITE(2,num_str_tot) Rc_matrix(i,:) 
        END DO

        DO i=1,alpha*3
            WRITE(2,num_str_tot) It_matrix(i,:) 
        END DO

        DO i=1,alpha*3
            WRITE(2,num_str_tot) Ic_matrix(i,:) 
        END DO

        DO i=1,alpha*3
            WRITE(2,'(3E30.17)') R0m_vector(i,:) 
        END DO

        DO i=1,alpha*3
            WRITE(2,'(3E30.17)') I0m_vector(i,:) 
        END DO

        DO i=1,alpha*3
            WRITE(2,'(3E30.17)') Rm0_vector(i,:) 
        END DO

        DO i=1,alpha*3
            WRITE(2,'(3E30.17)') Im0_vector(i,:) 
        END DO     
        
        DO i=1,3
            WRITE(2,'(3E30.17)') R00_matrix(i,:) 
        END DO 

        DO i=1,6
            WRITE(2,'(6F30.16)') C(i,:) 
        END DO 
                 
        DEALLOCATE(Rt_mat, Rc_mat)
        DEALLOCATE(It_mat, Ic_mat)
        DEALLOCATE(R0m_vet, Rm0_vet)
        DEALLOCATE(I0m_vet, Im0_vet)

        DEALLOCATE(Rt_matrix)
        DEALLOCATE(It_matrix)
        DEALLOCATE(Rc_matrix)
        DEALLOCATE(Ic_matrix)

        DEALLOCATE(R0m_vector, Rm0_vector)
        DEALLOCATE(I0m_vector, Im0_vector)


        CALL SYSTEM_CLOCK(t1)       ! Time counter stopped
        WRITE (*,'(A,F4.2,A)') 'COMPLETED! (Time :', REAL(t1-t0)/rate,'s)'
        WRITE (*,*) ''
        WRITE (*,*) '----------------------------------------------------------------------'
        ts1 = REAL(t1-t0)/rate
        
        CLOSE(2)
                
    END SUBROUTINE Print_data
!-------------------------------------------------------------------------------!
END MODULE Output_data
