!***********************************************************************************************************************************
!
!                                                       M 3 3 I N V _ M A I N
!
!  Program:      M33INV_MAIN
!
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!
!  Date:         July 22, 2005
!
!  Language:     Fortran-90
!
!  Version:      1.00b  (Feburary 7, 2009)
!
!  Description:  This program is a short "driver" to call function M33INV, which inverts a 3x3 matrix.
!
!  Files:        Source files:
!
!                   m33inv.f90                   Main program
!
!***********************************************************************************************************************************

!       PROGRAM M33INV_MAIN

!       IMPLICIT NONE

!       INTEGER :: I, J
!       DOUBLE PRECISION, DIMENSION(3,3) :: MAT, MATINV
!       LOGICAL :: OK_FLAG

!       LOGICAL :: M33INV

! !-----------------------------------------------------------------------------------------------------------------------------------

! !
! !     Get user input.
! !

!       WRITE (UNIT=*, FMT='(/A/)') ' Enter matrix:'

!       DO I = 1, 3
!          DO J = 1, 3
!             WRITE (UNIT=*, FMT='(A,I1,1H,,I1,A)', ADVANCE='NO') ' A(', I, J, ') = '
!             READ (UNIT=*, FMT=*) MAT(I,J)
!          END DO
!       END DO

! !
! !     Invert the input matrix.
! !

!       CALL M33INV (MAT, MATINV, OK_FLAG)

! !
! !     Print the result.
! !

!       IF (OK_FLAG) THEN
!          WRITE (UNIT=*, FMT='(/A/)') ' Inverse:'
!          WRITE (UNIT=*, FMT='(3ES25.15)') ((MATINV(I,J), J=1,3), I=1,3)
!       ELSE
!          WRITE (UNIT=*, FMT='(/A)') ' Singular matrix.'
!       END IF

!       STOP

!       END PROGRAM M33INV_MAIN



!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************

      SUBROUTINE M33INV (A, AINV, OK_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE M33INV

