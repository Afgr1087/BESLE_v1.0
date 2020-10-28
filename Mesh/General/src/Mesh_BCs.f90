!===============================================================================!
!-------------------------------------------------------------------------------!
!                               MODULE FOR INPUT DATA                           !
!-------------------------------------------------------------------------------!
!===============================================================================!
MODULE Mesh_BCs
!-------------------------------------------------------------------------------!
USE Global_variables
USE Custom_functions
!-------------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE start_greet
        WRITE(*,*) ''
        WRITE(*,*) '*****************************************************************************'
        WRITE(*,*) '*****************************************************************************'
        WRITE(*,*) '*****************************************************************************'
        WRITE(*,*) '*****************************************************************************'
        WRITE(*,*) '                  MESH AND BOUNDARY CONDITIONS GENERATOR                     '
        WRITE(*,*) '																			 '
		WRITE(*,*) 'By Daniel M. Prada and Andres F. Galvis'
		WRITE(*,*) 'School of Mechanical Engineering'
		WRITE(*,*) 'University of Campinas'
		WRITE(*,*) '17/02/2020'
        WRITE(*,*) '*****************************************************************************'
        WRITE(*,*) '*****************************************************************************'
        WRITE(*,*) '*****************************************************************************'
        WRITE(*,*) '*****************************************************************************'
    END SUBROUTINE start_greet


!===============================================================================!
    SUBROUTINE Mesh_BCs_1

	!=======================================================================
	! determination of the size of the class "Region"
	!=======================================================================
	!-----------------------------------------------------------------------
		!Initializations

	WRITE (*,*) '----------------------------------------------------------------------------'
        WRITE (*,*) ' '
        WRITE (*,'(A)',advance='no') ' 1. Initializations                               ...'
	WRITE (*,*) ''
	WRITE (*,*) ''

		stat=0
		nraux=0
		!opening input file to size the dimensions of the matrices
		open(2,file=filename_in,form='formatted')
			do while(stat==0)
				read(2,'(A100)',iostat=stat) line
			   	!reading the number regions
				    if(line(3:9)=='object') then
					nraux=nraux+1
					n_Regions=nraux		
				    endif
			enddo   
		close(2)
		allocate(Region(n_Regions))
	!-----------------------------------------------------------------------
	!=======================================================================
	!determination of the size of the classes "VerticesCoor" and "Elements"
	!of each class "Region"
	!=======================================================================
	!-----------------------------------------------------------------------
		!Initializations
		stat=0
		nraux=0
		!opening input file to size the dimensions of the matrices
		open(2,file=filename_in,form='formatted')
			do while(stat==0)
			read(2,'(A100)',iostat=stat) line
			!reading the number regions
			if(line(3:9)=='object') then
				nraux=nraux+1
				nvaux=0
				neaux=0
			endif
			!reading the number of verices of the region iRegion
			if(line(1:1)=='v') then
				nvaux=nvaux+1
				Region(nraux)%n_Vertices=nvaux
			endif
			!reading the number of elements of the region iRegion
			if(line(1:1)=='f') then
				neaux=neaux+1
				Region(nraux)%n_Elements=neaux
			endif
		enddo          
		close(2)
		do i=1,n_Regions
			nvaux=Region(i)%n_Vertices
	 	neaux=Region(i)%n_Elements
		allocate(Region(i)%VerticesCoor(nvaux))
		do j=1, nvaux
			Region(i)%VerticesCoor(j)%Coord=0.d0
		end do
		allocate(Region(i)%Elements(neaux))
		!write(*,*) 'nreg', i, 'nver', nvaux, 'nele', neaux
		end do
			
	!-----------------------------------------------------------------------
	!=======================================================================
	! Reading values from the mesh file
	!=======================================================================
	!-----------------------------------------------------------------------
	!-----------------------------------------------------------------------

	WRITE (*,*) '----------------------------------------------------------------------------'
        WRITE (*,*) ' '
        WRITE (*,'(A)',advance='no') ' 2. Reading values from the mesh file             ...'
	WRITE (*,*) ''
	WRITE (*,*) ''

		stat=0
		nraux=0
		NT_Vertices=0
		NT_Elements=0
		!Read the mesh input file

		open(2,file=filename_in,form='formatted')
			do while(stat==0)
			read(2,'(A100)',iostat=stat) line
			!reading regions
			if(line(3:9)=='object') then
				nraux=nraux+1
				nvaux=0
				neaux=0
			endif
			!reading verices 
			if(line(1:1)=='v') then
				nvaux=nvaux+1
				NT_Vertices=NT_Vertices+1
				read(line(3:),*)Region(nraux)%VerticesCoor(nvaux)%Coord(1), &
						Region(nraux)%VerticesCoor(nvaux)%Coord(2), &
						Region(nraux)%VerticesCoor(nvaux)%Coord(3)
		     	endif
			!reading elements 
			if(line(1:1)=='f') then
				NT_Elements=NT_Elements+1
				neaux=neaux+1
				line=line(3:)
				i=index(line,' ')
				j=index(line(:len_trim(line)),' ',.true.)
				block(1)=line(1:i-1)
				block(2)=line(i+1:j-1)
				block(3)=line(j+1:)           
				do k=1,3
					read(block(k),*)Region(nraux)%Elements(neaux)%Vertices(k)         
				enddo 
	 		endif
		enddo   
		close(2)
	!-----------------------------------------------------------------------
	!=======================================================================
	!Reading BCS files
	!=======================================================================
	!-----------------------------------------------------------------------

	WRITE (*,*) '----------------------------------------------------------------------------'
        WRITE (*,*) ' '
        WRITE (*,'(A)',advance='no') ' 3. Reading BCS files                             ...'
	WRITE (*,*) ''
	WRITE (*,*) ''

		do m = 1,NumBCS
	 
		write(BCSfilename_in,"(I10)") m
		BCSfilename_in="BCs_"//trim(adjustl(BCSfilename_in))//".obj" 
		nvaux=0
		neaux=0
		stat=0
		!opening input BCS file to measure the dimensions of the objects
		open(2,file=trim(BCSfilesPlace_in)//trim(BCSfilename_in),form='formatted')
			do while(stat==0)
				read(2,'(A100)',iostat=stat) line
				!reading the number of verices of the region iRegion
				if(line(1:1)=='v') then
					nvaux=nvaux+1
					BC(m)%n_BCVertices=nvaux
				endif
				!reading the number of elements of the region iRegion
				if(line(1:1)=='f') then
					neaux=neaux+1
					BC(m)%n_BCElements=neaux
				endif
			enddo          
		close(2)
			
		nvaux=BC(m)%n_BCVertices
	 	neaux=BC(m)%n_BCElements

		allocate(BC(m)%VerticesCoor(nvaux))
		do j=1, nvaux
			BC(m)%VerticesCoor(j)%Coord=0.d0
		end do
		allocate(BC(m)%Elements(neaux))
		enddo
	!--------------------------------------------------------------------------
	!==========================================================================
	!reading bcs
	!==========================================================================
	!--------------------------------------------------------------------------
		do m = 1,NumBCS
			write(BCSfilename_in,"(I10)") m
			BCSfilename_in="BCs_"//trim(adjustl(BCSfilename_in))//".obj"
			stat=0
			nvaux=0
			neaux=0
			!Reading BCS input file
			open(2,file=trim(BCSfilesPlace_in)//trim(BCSfilename_in),form='formatted')
				do while(stat==0)
					read(2,'(A100)',iostat=stat) line
					!read verices of the region iRegion
					if(line(1:1)=='v') then
						nvaux=nvaux+1
						read(line(3:),*)BC(m)%VerticesCoor(nvaux)%Coord(1), &
								BC(m)%VerticesCoor(nvaux)%Coord(2), &
								BC(m)%VerticesCoor(nvaux)%Coord(3)
			     		endif
					!read elements of the region iRegion
					if(line(1:1)=='f') then
						neaux=neaux+1
						line=line(3:)
						i=index(line,' ')
						j=index(line(:len_trim(line)),' ',.true.)
						block(1)=line(1:i-1)
						block(2)=line(i+1:j-1)
						block(3)=line(j+1:)           
						do k=1,3
							read(block(k),*)BC(m)%Elements(neaux)%Vertices(k)         
						enddo
	 				endif
				enddo   
			close(2)
			
		
		enddo
	!-----------------------------------------------------------------------
	!=======================================================================
	!formating data into matrices
	!=======================================================================
	!-----------------------------------------------------------------------

	WRITE (*,*) '----------------------------------------------------------------------------'
        WRITE (*,*) ' '
        WRITE (*,'(A)',advance='no') ' 4. Formating data into matrices                  ...'	
	WRITE (*,*) ''
	WRITE (*,*) ''

		!allocation of matrices
		Allocate(FACES(NT_Elements,2))
		Allocate(POINTS(NT_Vertices,4),POINTS_int(NT_Vertices,4))
		Allocate(ELEM(NT_Elements,5))
		Allocate(NORMAL_VECTORS(NT_Elements,3))
		Allocate(El_reg(n_Regions,2))
			
		!initialization of matrices
		FACES=0
		POINTS=0.d0
		ELEM=0
		NORMAL_VECTORS=0.d0
		El_reg=0
			
		!Matrix of faces
		do i=1,NT_Elements
		FACES(i,:)=(/i,1/)
		end do

		!Matrix of Points
		k=0
		do i=1,n_Regions
		nvaux=Region(i)%n_Vertices
		do j=1, nvaux
			k=k+1	
			POINTS(k,1)=k
			POINTS(k,2:4)=Region(i)%Verticescoor(j)%Coord
		end do
		end do

		!matrix of elements
		k=1
		neaux2=0
		do i=1,n_Regions
		neaux=Region(i)%n_Elements
		do j=1,neaux
			ELEM(K,1)=K
			ELEM(K,2)=K
			ELEM(K,3)=Region(i)%Elements(j)%Vertices(1)
			ELEM(K,4)=Region(i)%Elements(j)%Vertices(2)
			ELEM(K,5)=Region(i)%Elements(j)%Vertices(3)
			k=k+1
		end do
		end do
			
		!matrix of regions
		do i=1, n_Regions
		El_reg(i,1)=i
		El_reg(i,2)=Region(i)%n_Elements
		end do
		!-----------------------------------------------------------------------
		!-----------------------------------------------------------------------
		!=======================================================================
		!Filtering the numerical presition in order to consider only the number
		!of decimals set by "MeshNumPress" 
		!=======================================================================
		!-----------------------------------------------------------------------
		
		scale_mult = (10.d0)**(MeshNumPress)
		scale_div = (10.d0)**(-MeshNumPress)

		POINTS(:,2:4) = POINTS(:,2:4)*scale_mult
		do i=1,size(POINTS,1)
		POINTS_int(i,2) = IDNINT(POINTS(i,2))
		POINTS_int(i,3) = IDNINT(POINTS(i,3))
		POINTS_int(i,4) = IDNINT(POINTS(i,4)) 
			   
		end do

		do i=1,size(POINTS,1)
		POINTS(i,2) = REAL(POINTS_int(i,2),8)
		POINTS(i,3) = REAL(POINTS_int(i,3),8)
		POINTS(i,4) = REAL(POINTS_int(i,4),8)        
		end do
				
		POINTS(:,2:4) = POINTS(:,2:4)*scale_div
		
		do i=1,NumBCS
		
			do j=1,size(BC(i)%VerticesCoor)
			
				BCsVertReal(1)=BC(i)%VerticesCoor(j)%Coord(1)*scale_mult
				BCsVertReal(2)=BC(i)%VerticesCoor(j)%Coord(2)*scale_mult
				BCsVertReal(3)=BC(i)%VerticesCoor(j)%Coord(3)*scale_mult
				
				BCsVertInt(1)=IDNINT(BCsVertReal(1))
				BCsVertInt(2)=IDNINT(BCsVertReal(2))
				BCsVertInt(3)=IDNINT(BCsVertReal(3))
				
				BCsVertReal(1)=REAL(BCsVertInt(1),8)*scale_div
				BCsVertReal(2)=REAL(BCsVertInt(2),8)*scale_div
				BCsVertReal(3)=REAL(BCsVertInt(3),8)*scale_div
				
				BC(i)%VerticesCoor(j)%Coord(1)=BCsVertReal(1)
				BC(i)%VerticesCoor(j)%Coord(2)=BCsVertReal(2)
				BC(i)%VerticesCoor(j)%Coord(3)=BCsVertReal(3)
				
			end do
				
				
		end do
			    
		!-----------------------------------------------------------------------
		!=======================================================================
		! determination of the normals of each element
		!=======================================================================
		!-----------------------------------------------------------------------

		do i=1, NT_Elements

		nv1=ELEM(i,3)
		nv2=ELEM(i,4)
		nv3=ELEM(i,5)

		V1=POINTS(nv1,2:4)
		V2=POINTS(nv2,2:4)
		V3=POINTS(nv3,2:4)

		NA=V1-V3
		NB=V2-V3
				
		NC(1) = NA(2)*NB(3)-NA(3)*NB(2)
		NC(2) = NA(3)*NB(1)-NA(1)*NB(3)
		NC(3) = NA(1)*NB(2)-NA(2)*NB(1)

		Norm = (1.d0/dsqrt(NC(1)**2+NC(2)**2+NC(3)**2))*NC
		NORMAL_VECTORS(i,:)=Norm
				
		end do
		!-----------------------------------------------------------------------
		!=======================================================================
		! determining the global id of the bcs elements
		!=======================================================================
		!-----------------------------------------------------------------------
		Tol=2*(10.d0)**(-MeshNumPress)
		NT_BCElements=0
		do m=1,NumBCS
		neaux=BC(m)%n_BCElements
		write(*,*) 'Number of elements of BC',m,'=',neaux
		allocate(BC(m)%BCElements(neaux))
		do i=1, neaux !Elements of BCS
			V1BC=BC(m)%Elements(i)%Vertices(1)
			V2BC=BC(m)%Elements(i)%Vertices(2)
			V3BC=BC(m)%Elements(i)%Vertices(3)

			C1V1BC=BC(m)%VerticesCoor(V1BC)%Coord(1)
			C2V1BC=BC(m)%VerticesCoor(V1BC)%Coord(2)
			C3V1BC=BC(m)%VerticesCoor(V1BC)%Coord(3)

			C1V2BC=BC(m)%VerticesCoor(V2BC)%Coord(1)
			C2V2BC=BC(m)%VerticesCoor(V2BC)%Coord(2)
			C3V2BC=BC(m)%VerticesCoor(V2BC)%Coord(3)

			C1V3BC=BC(m)%VerticesCoor(V3BC)%Coord(1)
			C2V3BC=BC(m)%VerticesCoor(V3BC)%Coord(2)
			C3V3BC=BC(m)%VerticesCoor(V3BC)%Coord(3)

			do j=1, NT_Elements  !ELments of Mesh

				V1M=ELEM(j,3)
				V2M=ELEM(j,4)
				V3M=ELEM(j,5)
				
				C1V1M=POINTS(V1M,2)
				C2V1M=POINTS(V1M,3)
				C3V1M=POINTS(V1M,4)

				C1V2M=POINTS(V2M,2)
				C2V2M=POINTS(V2M,3)
				C3V2M=POINTS(V2M,4)

				C1V3M=POINTS(V3M,2)
				C2V3M=POINTS(V3M,3)
				C3V3M=POINTS(V3M,4)

				!case 1
				Er1 = ABS(C1V1BC-C1V1M)
				Er2 = ABS(C2V1BC-C2V1M)
				Er3 = ABS(C3V1BC-C3V1M) 
				Er4 = ABS(C1V2BC-C1V2M) 
				Er5 = ABS(C2V2BC-C2V2M) 
				Er6 = ABS(C3V2BC-C3V2M) 
				Er7 = ABS(C1V3BC-C1V3M) 
				Er8 = ABS(C2V3BC-C2V3M)
				Er9 = ABS(C3V3BC-C3V3M)
				
				if(Er1.LT.Tol .and. Er2.LT.Tol .and. Er3.LT.Tol) then
					if(Er4.LT.Tol .and. Er5.LT.Tol .and. Er6.LT.Tol) then
						if(Er7.LT.Tol .and. Er8.LT.Tol .and. Er9.LT.Tol) then
							BC(m)%BCElements(i)=j
							NT_BCElements=NT_BCElements+1
						end if 
		 			end if
				end if 

				!case 2
				Er1 = ABS(C1V1BC-C1V3M)
				Er2 = ABS(C2V1BC-C2V3M)
				Er3 = ABS(C3V1BC-C3V3M) 
				Er4 = ABS(C1V2BC-C1V1M) 
				Er5 = ABS(C2V2BC-C2V1M) 
				Er6 = ABS(C3V2BC-C3V1M) 
				Er7 = ABS(C1V3BC-C1V2M) 
				Er8 = ABS(C2V3BC-C2V2M)
				Er9 = ABS(C3V3BC-C3V2M)  
				if(Er1.LT.Tol .and. Er2.LT.Tol .and. Er3.LT.Tol) then
					if(Er4.LT.Tol .and. Er5.LT.Tol .and. Er6.LT.Tol) then
						if(Er7.LT.Tol .and. Er8.LT.Tol .and. Er9.LT.Tol) then
							BC(m)%BCElements(i)=j
							NT_BCElements=NT_BCElements+1
						end if 
		 			end if
				end if 
				
				!case 3
				Er1 = ABS(C1V1BC-C1V2M)
				Er2 = ABS(C2V1BC-C2V2M)
				Er3 = ABS(C3V1BC-C3V2M) 
				Er4 = ABS(C1V2BC-C1V3M) 
				Er5 = ABS(C2V2BC-C2V3M) 
				Er6 = ABS(C3V2BC-C3V3M) 
				Er7 = ABS(C1V3BC-C1V1M) 
				Er8 = ABS(C2V3BC-C2V1M)
				Er9 = ABS(C3V3BC-C3V1M)  
				if(Er1.LT.Tol .and. Er2.LT.Tol .and. Er3.LT.Tol) then
					if(Er4.LT.Tol .and. Er5.LT.Tol .and. Er6.LT.Tol) then
						if(Er7.LT.Tol .and. Er8.LT.Tol .and. Er9.LT.Tol) then
							BC(m)%BCElements(i)=j
							NT_BCElements=NT_BCElements+1
						end if 
		 			end if
				end if 
			enddo
		enddo
		enddo 
		write(*,*) 'Global number of elements of the BCs =',NT_BCElements
		!-----------------------------------------------------------------------
		!=======================================================================
		!rearrange of the elements of each interface. Two Elements which make an 
		!interface, have to have the same node's order
		!=======================================================================
		!-----------------------------------------------------------------------
		k=2
		do i=1,size(ELEM,1)
		do j=k,size(ELEM,1)

			E1NV1=ELEM(i,3)
			E1NV2=ELEM(i,4)
			E1NV3=ELEM(i,5)
			E2NV1=ELEM(j,3)
			E2NV2=ELEM(j,4)
			E2NV3=ELEM(j,5)

			E1V1X=POINTS(E1NV1,2)
			E1V1Y=POINTS(E1NV1,3)
			E1V1Z=POINTS(E1NV1,4)

			E1V2X=POINTS(E1NV2,2)
			E1V2Y=POINTS(E1NV2,3)
			E1V2Z=POINTS(E1NV2,4)

			E1V3X=POINTS(E1NV3,2)
			E1V3Y=POINTS(E1NV3,3)
			E1V3Z=POINTS(E1NV3,4)	

			E2V1X=POINTS(E2NV1,2)
			E2V1Y=POINTS(E2NV1,3)
			E2V1Z=POINTS(E2NV1,4)

			E2V2X=POINTS(E2NV2,2)
			E2V2Y=POINTS(E2NV2,3)
			E2V2Z=POINTS(E2NV2,4)

			E2V3X=POINTS(E2NV3,2)
			E2V3Y=POINTS(E2NV3,3)
			E2V3Z=POINTS(E2NV3,4)

				
			!case 1

			Er1 = ABS(E1V1X-E2V1X); 
			Er2 = ABS(E1V1Y-E2V1Y); 
			Er3 = ABS(E1V1Z-E2V1Z);
			Er4 = ABS(E1V2X-E2V3X); 
			Er5 = ABS(E1V2Y-E2V3Y); 
			Er6 = ABS(E1V2Z-E2V3Z);
			Er7 = ABS(E1V3X-E2V2X); 
			Er8 = ABS(E1V3Y-E2V2Y); 
			Er9 = ABS(E1V3Z-E2V2Z);

			if(Er1.LT.Tol .and. Er2.LT.Tol .and. Er3.LT.Tol) then
				if(Er4.LT.Tol .and. Er5.LT.Tol .and. Er6.LT.Tol) then
					if(Er7.LT.Tol .and. Er8.LT.Tol .and. Er9.LT.Tol) then
						ELEM(j,3:5)=ELEM(i,3:5)
						m=m+1
					end if 
		 		end if
			end if

			!case 2
			Er1 = ABS(E1V1X-E2V2X); 
			Er2 = ABS(E1V1Y-E2V2Y); 
			Er3 = ABS(E1V1Z-E2V2Z);
			Er4 = ABS(E1V2X-E2V1X); 
			Er5 = ABS(E1V2Y-E2V1Y); 
			Er6 = ABS(E1V2Z-E2V1Z);
			Er7 = ABS(E1V3X-E2V3X); 
			Er8 = ABS(E1V3Y-E2V3Y); 
			Er9 = ABS(E1V3Z-E2V3Z);
				
			if(Er1.LT.Tol .and. Er2.LT.Tol .and. Er3.LT.Tol) then
				if(Er4.LT.Tol .and. Er5.LT.Tol .and. Er6.LT.Tol) then
					if(Er7.LT.Tol .and. Er8.LT.Tol .and. Er9.LT.Tol) then
						ELEM(j,3:5)=ELEM(i,3:5) 
						m=m+1           
					end if 
		 		end if
			end if

			!case 3
			Er1 = ABS(E1V1X-E2V3X); 
			Er2 = ABS(E1V1Y-E2V3Y); 
			Er3 = ABS(E1V1Z-E2V3Z);
			Er4 = ABS(E1V3X-E2V1X); 
			Er5 = ABS(E1V3Y-E2V1Y); 
			Er6 = ABS(E1V3Z-E2V1Z);
			Er7 = ABS(E1V2X-E2V2X); 
			Er8 = ABS(E1V2Y-E2V2Y); 
			Er9 = ABS(E1V2Z-E2V2Z);

				
			if(Er1.LT.Tol .and. Er2.LT.Tol .and. Er3.LT.Tol) then
				if(Er4.LT.Tol .and. Er5.LT.Tol .and. Er6.LT.Tol) then
					if(Er7.LT.Tol .and. Er8.LT.Tol .and. Er9.LT.Tol) then
						ELEM(j,3:5)=ELEM(i,3:5)
						m=m+1
					end if 
		 		end if
			end if

			!case 4
			Er1 = ABS(E1V1X-E2V2X); 
			Er2 = ABS(E1V1Y-E2V2Y); 
			Er3 = ABS(E1V1Z-E2V2Z)
			Er4 = ABS(E1V2X-E2V3X); 
			Er5 = ABS(E1V2Y-E2V3Y); 
			Er6 = ABS(E1V2Z-E2V3Z)
			Er7 = ABS(E1V3X-E2V1X); 
			Er8 = ABS(E1V3Y-E2V1Y); 
			Er9 = ABS(E1V3Z-E2V1Z)
				
			if(Er1.LT.Tol .and. Er2.LT.Tol .and. Er3.LT.Tol) then
				if(Er4.LT.Tol .and. Er5.LT.Tol .and. Er6.LT.Tol) then
					if(Er7.LT.Tol .and. Er8.LT.Tol .and. Er9.LT.Tol) then
						ELEM(j,3:5)=ELEM(i,3:5)
						m=m+1
					end if 
		 		end if
			end if

			!case 5
			Er1 = ABS(E1V1X-E2V3X); 
			Er2 = ABS(E1V1Y-E2V3Y); 
			Er3 = ABS(E1V1Z-E2V3Z)
			Er4 = ABS(E1V2X-E2V1X); 
			Er5 = ABS(E1V2Y-E2V1Y); 
			Er6 = ABS(E1V2Z-E2V1Z)
			Er7 = ABS(E1V3X-E2V2X); 
			Er8 = ABS(E1V3Y-E2V2Y); 
			Er9 = ABS(E1V3Z-E2V2Z)
		 
			if(Er1.LT.Tol .and. Er2.LT.Tol .and. Er3.LT.Tol) then
				if(Er4.LT.Tol .and. Er5.LT.Tol .and. Er6.LT.Tol) then
					if(Er7.LT.Tol .and. Er8.LT.Tol .and. Er9.LT.Tol) then
						ELEM(j,3:5)=ELEM(i,3:5)
						m=m+1
					end if 
		 		end if
			end if
				
		end do
		k=k+1	
		end do
	!------------------------------------------------------------------------------------
	!====================================================================================
	!Boundary conditions Area's Evaluation
	!====================================================================================
	!------------------------------------------------------------------------------------
	do i=1, NumBCS
		EArea=0.d0
		do j=1, BC(i)%n_BCElements

			V1BC=BC(i)%Elements(j)%Vertices(1)
			V2BC=BC(i)%Elements(j)%Vertices(2)
			V3BC=BC(i)%Elements(j)%Vertices(3)

			C1V1BC=BC(i)%VerticesCoor(V1BC)%Coord(1)
			C2V1BC=BC(i)%VerticesCoor(V1BC)%Coord(2)
			C3V1BC=BC(i)%VerticesCoor(V1BC)%Coord(3)

			C1V2BC=BC(i)%VerticesCoor(V2BC)%Coord(1)
			C2V2BC=BC(i)%VerticesCoor(V2BC)%Coord(2)
			C3V2BC=BC(i)%VerticesCoor(V2BC)%Coord(3)

			C1V3BC=BC(i)%VerticesCoor(V3BC)%Coord(1)
			C2V3BC=BC(i)%VerticesCoor(V3BC)%Coord(2)
			C3V3BC=BC(i)%VerticesCoor(V3BC)%Coord(3)
			
			AuxArea1=(C2V2BC-C2V1BC)*(C3V3BC-C3V1BC)-(C3V2BC-C3V1BC)*(C2V3BC-C2V1BC)
			AuxArea2=(C1V2BC-C1V1BC)*(C3V3BC-C3V1BC)-(C3V2BC-C3V1BC)*(C1V3BC-C1V1BC)
			AuxArea3=(C1V2BC-C1V1BC)*(C2V3BC-C2V1BC)-(C2V2BC-C2V1BC)*(C1V3BC-C1V1BC)

			EArea=EArea+dsqrt(AuxArea1**2+AuxArea2**2+AuxArea3**2)/2
				
		end do
		BC(i)%Area=EArea
		write(*,*) 'The area of the BC',i,'is:', BC(i)%Area, '[mesh units]'
	end do
	!------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------
	!====================================================================================
	!Boundary conditions data arrangement
	!====================================================================================
	!------------------------------------------------------------------------------------
	!if the Analysis is static the number of load steps is 1 
	!------------------------------------------------------------------------------------
	if (NumBCS .gt. 0) then
		
	WRITE (*,*) '----------------------------------------------------------------------------'
        WRITE (*,*) ' '
        WRITE (*,'(A)',advance='no') ' 5. Boundary conditions data arrangement          ...'	
	WRITE (*,*) ''
	WRITE (*,*) ''

   	 do i=1,NumBCS
		select case(BC(i)%DirType)
			case("xyz")

				select case(BC(i)%BCxType)
					case("displacement")
						BC(i)%IDtype(1)=0 
					case("traction")
						BC(i)%IDtype(1)=1
					case("free")
						BC(i)%IDtype(1)=1
					case default
						write(*,*) "Error in BC_",i
						write(*,*) "The Boundary Condition (mode=xyz, direction x) was not set correctly."
						write(*,*) "Options are (displacement, traction, free)"
						write(*,*) "Remember it is case sensitive, please check it."
						call abort 
				end select
	
				select case(BC(i)%BCyType)
					case("displacement")
						BC(i)%IDtype(2)=0
					case("traction")
						BC(i)%IDtype(2)=1
					case("free")
						BC(i)%IDtype(2)=1
					case default
						write(*,*) "Error in BC_",i
						write(*,*) "The Boundary Condition (mode=xyz, direction y) was not set correctly."
						write(*,*) "Options are (displacement, traction, free)"
						write(*,*) "Remember it is case sensitive, please check it."
						call abort 
				end select
	
				select case(BC(i)%BCzType)
					case("displacement")
						BC(i)%IDtype(3)=0
					case("traction")
						BC(i)%IDtype(3)=1
					case("free")
						BC(i)%IDtype(3)=1
					case default
						write(*,*) "Error in BC_",i
						write(*,*) "The Boundary Condition (mode=xyz, direction z) was not set correctly."
						write(*,*) "Options are (displacement, traction, free)"
						write(*,*) "Remember it is case sensitive, please check it."
						call abort 
				end select
	
			case("normal")
	
				select case(BC(i)%BCnType)
					case("displacement")
						BC(i)%IDtype(1)=2 
					case("traction")
						BC(i)%IDtype(1)=3
					case default
						write(*,*) "Error in BC_",i
						write(*,*) "The Boundary Condition (mode=normal) was not set correctly."
						write(*,*) "Options are (displacement, traction)"
						write(*,*) "Remember it is case sensitive, please check it."
						call abort
				end select
	
			case("load")
				BC(i)%IDtype(1:3)=1 
	
			case default
				write(*,*) "Error in BC_",i
				write(*,*) "The Direction mode was not set correctly."
				write(*,*) "Remember it is case sensitive, please check it."
				write(*,*) "Options are (xyz, normal, load)"
				call abort 
			
		end select
    	enddo
	!----------------------------------------------------------------------------------------------------
	do i=1,NumBCS
		select case(BC(i)%DirType)
			case("xyz")
			if(AnalysisType=="elastostatic") then
	
				Nsteps=1
				Allocate(BC(i)%BCxVal(1))
				Allocate(BC(i)%BCyVal(1))
				Allocate(BC(i)%BCzVal(1))
	
				if(BC(i)%BCxType.eq."free") then
					BC(i)%BCxVal(1)=0.d0
				else
					BC(i)%BCxVal(1)=BC(i)%EParam%staticBCxVal
				end if
	
				if(BC(i)%BCyType.eq."free") then
					BC(i)%BCyVal(1)=0.d0
				else
					BC(i)%BCyVal(1)=BC(i)%EParam%staticBCyVal
				end if
	
				if(BC(i)%BCzType.eq."free") then
					BC(i)%BCzVal(1)=0.d0
				else
					BC(i)%BCzVal(1)=BC(i)%EParam%staticBCzVal
				end if
	
			else if(AnalysisType=="quasi-elastostatic".or. AnalysisType=="elastodynamic") then
	
				Allocate(BC(i)%BCxVal(Nsteps))
				Allocate(BC(i)%BCyVal(Nsteps))
				Allocate(BC(i)%BCzVal(Nsteps))
				!--------------------------------------------------------------------------------------------------------------------------------------
				if(BC(i)%BCxType.eq."free") then
					do j=1,Nsteps
						BC(i)%BCxVal(j)=0.d0
					end do
				else
					select case(BC(i)%FuncxType)
						case("linear")
							do j=1,Nsteps
								BC(i)%BCxVal(j)=(1.d0/real(Nsteps,8))*real(j,8)*BC(i)%LParam%LinearBCxMaxVal
							end do
						case("quadratic")
							BC(i)%QParam%QuadSAux=BC(i)%QParam%QuadSOx
							do j=1,Nsteps
								BC(i)%QParam%QuadSAux= BC(i)%QParam%QuadSAux + (BC(i)%QParam%QuadSfx - BC(i)%QParam%QuadSOx)/real(Nsteps,8)
								BC(i)%BCxVal(j)= BC(i)%QParam%QuadAx * (BC(i)%QParam%QuadSAux)**(2.d0) + BC(i)%QParam%QuadBx * &
								BC(i)%QParam%QuadSAux + BC(i)%QParam%QuadCx
							end do
						case("sine")
							BC(i)%SParam%SineSAux= BC(i)%SParam%SineSOx
							do j=1,Nsteps
								BC(i)%SParam%SineSAux= BC(i)%SParam%SineSAux + (BC(i)%SParam%SineSfx - BC(i)%SParam%SineSOx)/real(Nsteps,8)
								BC(i)%BCxVal(j)= BC(i)%SParam%SineAx * sin( BC(i)%SParam%SineBx * BC(i)%SParam%SineSAux + BC(i)%SParam%SineCx )
							end do
						case("custom_1")
							BC(i)%BCxVal=CustomFunc1(nsteps)
						case("custom_2")
							BC(i)%BCxVal=CustomFunc2(nsteps)
						case("custom_3")
							BC(i)%BCxVal=CustomFunc3(nsteps)
						case default
							write(*,*) "Error in BC_",i
							write(*,*) "The BC Function (mode=xyz, x direction) was not set correctly."
							write(*,*) "Remember it is case sensitive, please check it."
							write(*,*) "Options are (linear, quadratic, sine, custom_1, custom_2, custom_3)"
							call abort
					end select
				end if
				!--------------------------------------------------------------------------------------------------------------------------------------
				if(BC(i)%BCyType.eq."free") then
					do j=1,Nsteps
						BC(i)%BCyVal(j)=0.d0
					end do
				else
					select case(BC(i)%FuncyType)

						case("linear")
							do j=1,Nsteps
								BC(i)%BCyVal(j)=(1.d0/real(Nsteps,8))*real(j,8)*BC(i)%LParam%LinearBCyMaxVal
							end do
						case("quadratic")
							BC(i)%QParam%QuadSAux=BC(i)%QParam%QuadSOy
							do j=1,Nsteps
								BC(i)%QParam%QuadSAux= BC(i)%QParam%QuadSAux + (BC(i)%QParam%QuadSfy - BC(i)%QParam%QuadSOy)/real(Nsteps,8)
								BC(i)%BCyVal(j)= BC(i)%QParam%QuadAy * (BC(i)%QParam%QuadSAux)**(2.d0) + BC(i)%QParam%QuadBy * &
								BC(i)%QParam%QuadSAux + BC(i)%QParam%QuadCy
							end do
						case("sine")
							BC(i)%SParam%SineSAux= BC(i)%SParam%SineSOy
							do j=1,Nsteps
								BC(i)%SParam%SineSAux= BC(i)%SParam%SineSAux + (BC(i)%SParam%SineSfy - BC(i)%SParam%SineSOy)/real(Nsteps,8)
								BC(i)%BCyVal(j)= BC(i)%SParam%SineAy * sin( BC(i)%SParam%SineBy * BC(i)%SParam%SineSAux + BC(i)%SParam%SineCy )
							end do
						case("custom_1")
							BC(i)%BCyVal=CustomFunc1(nsteps)
						case("custom_2")
							BC(i)%BCyVal=CustomFunc2(nsteps)
						case("custom_3")
							BC(i)%BCyVal=CustomFunc3(nsteps)

						case default
							write(*,*) "Error in BC_",i
							write(*,*) "The BC Function (mode=xyz, y direction) was not set correctly."
							write(*,*) "Remember it is case sensitive, please check it."
							write(*,*) "Options are (linear, quadratic, sine, custom_1, custom_2, custom_3)"
							call abort

					end select
				end if
				!---------------------------------------------------------------------------------------------------------------------------------------
				if(BC(i)%BCzType.eq."free") then
					do j=1,Nsteps
						BC(i)%BCzVal(j)=0.d0
					end do
				else
					select case(BC(i)%FunczType)
						case("linear")
							do j=1,Nsteps
								BC(i)%BCzVal(j)=(1.d0/real(Nsteps,8))*real(j,8)*BC(i)%LParam%LinearBCzMaxVal
							end do
						case("quadratic")
							BC(i)%QParam%QuadSAux=BC(i)%QParam%QuadSOz
							do j=1,Nsteps
								BC(i)%QParam%QuadSAux= BC(i)%QParam%QuadSAux + (BC(i)%QParam%QuadSfz - BC(i)%QParam%QuadSOz)/real(Nsteps,8)
								BC(i)%BCzVal(j)= BC(i)%QParam%QuadAz * (BC(i)%QParam%QuadSAux)**(2.d0) + BC(i)%QParam%QuadBz * &
								BC(i)%QParam%QuadSAux + BC(i)%QParam%QuadCz
							end do
						case("sine")
							BC(i)%SParam%SineSAux= BC(i)%SParam%SineSOz
							do j=1,Nsteps
								BC(i)%SParam%SineSAux= BC(i)%SParam%SineSAux + (BC(i)%SParam%SineSfz - BC(i)%SParam%SineSOz)/real(Nsteps,8)
								BC(i)%BCzVal(j)= BC(i)%SParam%SineAz * sin( BC(i)%SParam%SineBz * BC(i)%SParam%SineSAux + BC(i)%SParam%SineCz )
							end do
						case("custom_1")
							BC(i)%BCzVal=CustomFunc1(nsteps)
						case("custom_2")
							BC(i)%BCzVal=CustomFunc2(nsteps)
						case("custom_3")
							BC(i)%BCzVal=CustomFunc3(nsteps)

						case default
							write(*,*) "Error in BC_",i
							write(*,*) "The BC Function (mode=xyz, z direction) was not set correctly."
							write(*,*) "Remember it is case sensitive, please check it."
							write(*,*) "Options are (linear, quadratic, sine, custom_1, custom_2, custom_3)"
							call abort

					end select
				end if
				!---------------------------------------------------------------------------------------------------------------------------------------			
	
   		 	else
				write(*,*) "Error in BCS"
				write(*,*) "The type of analysis was not set correctly."
				write(*,*) "Remember it is case sensitive, please check it."
				write(*,*) "Options are (elastostatic, quasi-elastostatic, elastodynamic)"
				call abort 
			end if
			case("normal")
			if(AnalysisType=="elastostatic") then
	
				Nsteps=1
				Allocate(BC(i)%BCxVal(1))
				Allocate(BC(i)%BCyVal(1))
				Allocate(BC(i)%BCzVal(1))
	
				BC(i)%BCxVal(1)=BC(i)%EParam%staticBCnVal
				BC(i)%BCyVal(1)=0.d0
				BC(i)%BCzVal(1)=0.d0
				
	
			else if(AnalysisType=="quasi-elastostatic".or. AnalysisType=="elastodynamic") then
	
				Allocate(BC(i)%BCxVal(Nsteps))
				Allocate(BC(i)%BCyVal(Nsteps))
				Allocate(BC(i)%BCzVal(Nsteps))

				select case(BC(i)%FuncnType)
						case("linear")
							do j=1,Nsteps
								BC(i)%BCxVal(j)=(1.d0/real(Nsteps,8))*real(j,8)*BC(i)%LParam%LinearBCnMaxVal
								BC(i)%BCyVal(j)=0.d0
								BC(i)%BCzVal(j)=0.d0
							end do
						case("quadratic")
							BC(i)%QParam%QuadSAux=BC(i)%QParam%QuadSOn
							do j=1,Nsteps
								BC(i)%QParam%QuadSAux= BC(i)%QParam%QuadSAux + (BC(i)%QParam%QuadSfn - BC(i)%QParam%QuadSOn)/real(Nsteps,8)
								BC(i)%BCxVal(j)= BC(i)%QParam%QuadAn * (BC(i)%QParam%QuadSAux)**(2.d0) + BC(i)%QParam%QuadBn * &
								BC(i)%QParam%QuadSAux + BC(i)%QParam%QuadCn
								BC(i)%BCyVal(j)=0.d0
								BC(i)%BCzVal(j)=0.d0
							end do
						case("sine")
							BC(i)%SParam%SineSAux= BC(i)%SParam%SineSOn
							do j=1,Nsteps
								BC(i)%SParam%SineSAux= BC(i)%SParam%SineSAux + (BC(i)%SParam%SineSfn - BC(i)%SParam%SineSOn)/real(Nsteps,8)
								BC(i)%BCxVal(j)= BC(i)%SParam%SineAn * sin( BC(i)%SParam%SineBn * BC(i)%SParam%SineSAux + BC(i)%SParam%SineCn )
								BC(i)%BCyVal(j)=0.d0
								BC(i)%BCzVal(j)=0.d0
							end do
						case("custom_1")
							BC(i)%BCxVal=CustomFunc1(nsteps)
							BC(i)%BCyVal=0.d0
							BC(i)%BCzVal=0.d0
						case("custom_2")
							BC(i)%BCxVal=CustomFunc2(nsteps)
							BC(i)%BCyVal=0.d0
							BC(i)%BCzVal=0.d0
						case("custom_3")
							BC(i)%BCxVal=CustomFunc3(nsteps)
							BC(i)%BCyVal=0.d0
							BC(i)%BCzVal=0.d0
						case default
							write(*,*) "Error in BC_",i
							write(*,*) "The BC Function (mode=normal) was not set correctly."
							write(*,*) "Remember it is case sensitive, please check it."
							write(*,*) "Options are (linear, quadratic, sine, custom_1, custom_2, custom_3)"
							call abort
				end select
			end if
			case("load")
				
				if(AnalysisType=="elastostatic") then
	
				Nsteps=1
				Allocate(BC(i)%BCxVal(1))
				Allocate(BC(i)%BCyVal(1))
				Allocate(BC(i)%BCzVal(1))
	
				BC(i)%BCxVal(1)=(BC(i)%EParam%staticBClVal*BC(i)%LoadDirection(1))/BC(i)%Area
				BC(i)%BCyVal(1)=(BC(i)%EParam%staticBClVal*BC(i)%LoadDirection(2))/BC(i)%Area
				BC(i)%BCzVal(1)=(BC(i)%EParam%staticBClVal*BC(i)%LoadDirection(3))/BC(i)%Area
				
	
				else if(AnalysisType=="quasi-elastostatic".or. AnalysisType=="elastodynamic") then
					
				Allocate(BC(i)%BCxVal(Nsteps))
				Allocate(BC(i)%BCyVal(Nsteps))
				Allocate(BC(i)%BCzVal(Nsteps))
				
				select case(BC(i)%FunclType)

						case("linear")
							do j=1,Nsteps
								BC(i)%BCxVal(j)=((1.d0/real(Nsteps,8))*real(j,8)*BC(i)%LParam%LinearBClMaxVal*BC(i)%LoadDirection(1))/BC(i)%Area
								BC(i)%BCyVal(j)=((1.d0/real(Nsteps,8))*real(j,8)*BC(i)%LParam%LinearBClMaxVal*BC(i)%LoadDirection(2))/BC(i)%Area
								BC(i)%BCzVal(j)=((1.d0/real(Nsteps,8))*real(j,8)*BC(i)%LParam%LinearBClMaxVal*BC(i)%LoadDirection(3))/BC(i)%Area
							end do
						case("quadratic")
							BC(i)%QParam%QuadSAux=BC(i)%QParam%QuadSOl
							do j=1,Nsteps
								BC(i)%QParam%QuadSAux= BC(i)%QParam%QuadSAux + (BC(i)%QParam%QuadSfl - BC(i)%QParam%QuadSOl)/real(Nsteps,8)
								BC(i)%BCxVal(j)= ((BC(i)%QParam%QuadAl * (BC(i)%QParam%QuadSAux)**(2.d0) + BC(i)%QParam%QuadBl * &
								BC(i)%QParam%QuadSAux + BC(i)%QParam%QuadCl)*BC(i)%LoadDirection(1))/BC(i)%Area
								BC(i)%BCyVal(j)=((BC(i)%QParam%QuadAl * (BC(i)%QParam%QuadSAux)**(2.d0) + BC(i)%QParam%QuadBl * &
								BC(i)%QParam%QuadSAux + BC(i)%QParam%QuadCl)*BC(i)%LoadDirection(2))/BC(i)%Area
								BC(i)%BCzVal(j)=((BC(i)%QParam%QuadAl * (BC(i)%QParam%QuadSAux)**(2.d0) + BC(i)%QParam%QuadBl * &
								BC(i)%QParam%QuadSAux + BC(i)%QParam%QuadCl)*BC(i)%LoadDirection(3))/BC(i)%Area
							end do
						case("sine")
							BC(i)%SParam%SineSAux = BC(i)%SParam%SineSOl
							do j=1,Nsteps
								BC(i)%SParam%SineSAux= BC(i)%SParam%SineSAux + (BC(i)%SParam%SineSfl - BC(i)%SParam%SineSOl)/real(Nsteps,8)
								BC(i)%BCxVal(j)= ((BC(i)%SParam%SineAl * sin( BC(i)%SParam%SineBl * BC(i)%SParam%SineSAux + BC(i)%SParam%SineCl )) * &
								BC(i)%LoadDirection(1))/BC(i)%Area
								BC(i)%BCyVal(j)=((BC(i)%SParam%SineAl * sin( BC(i)%SParam%SineBl * BC(i)%SParam%SineSAux + BC(i)%SParam%SineCl )) * &
								BC(i)%LoadDirection(2))/BC(i)%Area
								BC(i)%BCzVal(j)=((BC(i)%SParam%SineAl * sin( BC(i)%SParam%SineBl * BC(i)%SParam%SineSAux + BC(i)%SParam%SineCl )) * &
								BC(i)%LoadDirection(3))/BC(i)%Area
							end do
						case("custom_1")
							BC(i)%BCxVal=(CustomFunc1(nsteps)*BC(i)%LoadDirection(1))/BC(i)%Area
							BC(i)%BCyVal=(CustomFunc1(nsteps)*BC(i)%LoadDirection(2))/BC(i)%Area
							BC(i)%BCzVal=(CustomFunc1(nsteps)*BC(i)%LoadDirection(3))/BC(i)%Area
						case("custom_2")
							BC(i)%BCxVal=(CustomFunc2(nsteps)*BC(i)%LoadDirection(1))/BC(i)%Area
							BC(i)%BCyVal=(CustomFunc2(nsteps)*BC(i)%LoadDirection(2))/BC(i)%Area
							BC(i)%BCzVal=(CustomFunc2(nsteps)*BC(i)%LoadDirection(3))/BC(i)%Area
						case("custom_3")
							BC(i)%BCxVal=(CustomFunc3(nsteps)*BC(i)%LoadDirection(1))/BC(i)%Area
							BC(i)%BCyVal=(CustomFunc3(nsteps)*BC(i)%LoadDirection(2))/BC(i)%Area
							BC(i)%BCzVal=(CustomFunc3(nsteps)*BC(i)%LoadDirection(3))/BC(i)%Area
						case default
							write(*,*) "Error in BC_",i
							write(*,*) "The BC Function (mode=load) was not set correctly."
							write(*,*) "Remember it is case sensitive, please check it."
							write(*,*) "Options are (linear, quadratic, sine, custom_1, custom_2, custom_3)"
							call abort
				end select
				end if

		end select	
	end do
	end if
	WRITE (*,*) '----------------------------------------------------------------------------'
        WRITE (*,*) ' '
        WRITE (*,'(A)',advance='no') ' 6. writing .dat mesh and BCs files               ...'	
	WRITE (*,*) ''
	WRITE (*,*) ''
	WRITE (*,*) '----------------------------------------------------------------------------'
	
	!====================================================================================
	!------------------------------------------------------------------------------------
	!writing mesh .dat file

		OPEN(3,file='Meshes/DAT/'//filename_out,STATUS='REPLACE')
		    
		REWIND(3)

		WRITE(3,'(1I10)') n_Regions 

		WRITE(3,'(2I10)') NT_Vertices, 3

		DO i=1,NT_Vertices
			WRITE(3,'(3F30.15)') POINTS(i,2:4)
		END DO

		WRITE(3,'(2I10)') NT_Elements, 3
			   
		DO i=1,NT_Elements
			WRITE(3,'(3I10)') ELEM(i,3:5)
		END DO

		WRITE(3,'(2I10)') NT_Elements, 3
			   
		DO i=1,NT_Elements
			WRITE(3,'(3F30.15)') NORMAL_VECTORS(i,:)
		END DO

		WRITE(3,'(2I10)') n_Regions, 1

		DO i=1,n_Regions
			WRITE(3,'(1I10)') El_reg(i,2)
		END DO	

		CLOSE(3)

	!====================================================================================
	!------------------------------------------------------------------------------------
	!writing BCS .dat file
	if (NumBCS .gt. 0) then
		do i=1,NumBCS
		write(BCSfilename_aux,"(I10)") i
		BCSfilename_out = "BCs_"//trim(adjustl(BCSfilename_aux))//".dat"
		OPEN(3,file='BCs/BCs_DAT/'//trim(BCSfilename_out),STATUS='REPLACE')


		WRITE(3,'(1I10)') size(BC(i)%BCElements)
		do j=1,size(BC(i)%BCElements)
			WRITE(3,'(1I10)') BC(i)%BCElements(j)
		end do
		WRITE(3,'(3I10)') BC(i)%IDtype(1),BC(i)%IDtype(2),BC(i)%IDtype(3)
		do j=1, Nsteps
			WRITE(3,'(3F27.15)') BC(i)%BCxVal(j),BC(i)%BCyVal(j),BC(i)%BCzVal(j)
		end do
		CLOSE(3)


		end do
	end if

    END SUBROUTINE Mesh_BCs_1
!===============================================================================!



END MODULE Mesh_BCs
