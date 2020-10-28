MODULE Custom_functions
CONTAINS	
	FUNCTION CustomFunc1(nstepsaux) !Cubic Function
	
		REAL(8) :: CustomFunc1(nstepsaux)
		integer :: nstepsaux, i
		real(8) :: A, B, C, D, So, Sf, Sdeltha, S
		So=-0.62081d0
		Sf=2.84326d0
		Sdeltha=(Sf-So)/real(nstepsaux,8)
		A=-0.6d0
		B=1.9d0
		C=-0.2d0
		D=-1.d0
		S=So
		do i=1, nstepsaux
			S=S+Sdeltha
			CustomFunc1(i) = A*(S**3.d0) + B*(S**2.d0) + C*S + D
		end do
	
	END FUNCTION customFunc1
	
	FUNCTION CustomFunc2(nstepsaux) !Piecewise Funtion (y=0.5s, y=1, y=sin(4s-12)+1)
		REAL(8) :: CustomFunc2(nstepsaux)
		integer :: nstepsaux, i
		real(8) :: Sdeltha, S
		real(8) :: A1, So1
		real(8) :: So2
		real(8) :: A3, B3, C3,So3, Sf3
		real(8), parameter :: PI= 3.141592653589793d0
		A1=0.5d0
		So1=0.d0
		So2=2.d0
		So3=3.d0
		Sf3=3.d0+2.d0*pi
		A3=4.d0
		B3=-12.d0
		C3=1.d0
		Sdeltha=(abs(Sf3-So1))/real(nstepsaux,8)
		S=So1
		do i=1, nstepsaux
		S=S+Sdeltha

			if((S.gt.So1) .and. (S.lt.So2)) then
				CustomFunc2(i) = A1*S
			end if

			if( ((S.eq.So2) .or. (S.gt.So2)) .and. (S.lt.So3) ) then
				CustomFunc2(i) = 1.d0
			end if

			if( (S.eq.So3) .or. (S.gt.So3) ) then
				CustomFunc2(i) = sin(A3*S+B3)+C3
			end if
		
		end do
	
	END FUNCTION customFunc2
	
	FUNCTION CustomFunc3(nstepsaux) !Piecewise Funtion (y=2s, y=2, y=-0.05e^S+3)
		REAL(8) :: CustomFunc3(nstepsaux)
		integer :: nstepsaux, i
		real(8) :: Sdeltha, S
		real(8) :: A1, So1
		real(8) :: So2
		real(8) :: A3, B3, So3, Sf3
		real(8), parameter :: euler= 2.71828182845904d0
		A1=2.d0
		So1=0.d0
		So2=1.d0
		So3=3.d0
		Sf3=4.1d0
		A3=-0.05d0
		B3=3.d0
		Sdeltha=(abs(Sf3-So1))/real(nstepsaux,8)
		S=So1
		do i=1, nstepsaux
		S=S+Sdeltha

			if((S.gt.So1) .and. (S.lt.So2)) then
				CustomFunc3(i) = A1*S
			end if

			if( ((S.eq.So2) .or. (S.gt.So2)) .and. (S.lt.So3) ) then
				CustomFunc3(i) = 2.d0
			end if

			if( (S.eq.So3) .or. (S.gt.So3) ) then
				CustomFunc3(i) = A3*euler**S+B3
			end if
		
		end do

	END FUNCTION customFunc3

END MODULE Custom_functions
