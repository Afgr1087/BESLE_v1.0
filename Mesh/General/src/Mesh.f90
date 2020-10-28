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
! Program: 3D_BEM_Aniso                                                              !
! Description: 3D-BEM PROGRAM TO SOLVE LINEAR ELASTICITY PROBLEMS                    !
!               Using the fundamental solution based on Fourier series               !
!               and 6-node triangular Discontinuous elements                         !
!====================================================================================!
!====================================================================================!
PROGRAM Main
!
!====================================================================================!
! MODULES
!------------------------------------------------------------------------------------!
    USE Global_variables
    USE Set_parameters
    USE Custom_functions
    USE Mesh_BCs

!====================================================================================!
! VARIABLES
!------------------------------------------------------------------------------------!
    IMPLICIT NONE
!====================================================================================!
! INPUT DATA
!------------------------------------------------------------------------------------!

	CALL start_greet

	CALL Setup

	CALL Mesh_BCs_1
    
END PROGRAM Main
