SUBROUTINE read_orb(nbf,orb)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: nbf
 DOUBLE PRECISION, INTENT(OUT) :: orb(nbf,1)
 INTEGER :: i, io_open
 REAL :: temp(nbf,1)

 OPEN(UNIT=100,FILE="pk_guess.in",STATUS="OLD",ACTION="READ",IOSTAT=io_open)
 IF (io_open .NE. 0) THEN
  STOP "FAILED TO OPEN FILE: pk_guess.in"
 ELSE
  DO i=1,nbf
   READ(100,*) temp(i,1)
  END DO
 END IF
 CLOSE(100)
 orb = DBLE(temp)
RETURN
END SUBROUTINE
