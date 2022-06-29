SUBROUTINE read_coords(na,cmat)
 USE constants
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: na
 DOUBLE PRECISION, INTENT(OUT) :: cmat(na,3)
 REAL :: x, y, z
 INTEGER :: i

 OPEN(UNIT=100,FILE="coords.in",STATUS="OLD",ACTION="READ")
 DO i = 1, na
  READ(100,*) x, y, z
  cmat(i,1) = atb*DBLE(x); cmat(i,2) = atb*DBLE(y); cmat(i,3) = atb*DBLE(z)
 END DO
 CLOSE(100) 

END SUBROUTINE
