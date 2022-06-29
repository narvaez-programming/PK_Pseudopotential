SUBROUTINE read_t(nbf,tmat)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: nbf
 DOUBLE PRECISION, INTENT(OUT) :: tmat(nbf,nbf)
 INTEGER :: i, j, io_open

 OPEN(UNIT=100,FILE="kinetic.in",STATUS="OLD",ACTION="READ",IOSTAT=io_open)
 IF (io_open .NE. 0) THEN
  STOP "COULD NOT READ kinetic.in FILE"
 ELSE 
  DO i=1,nbf
  READ(100,*) (tmat(i,j), j=1,nbf) 
  END DO
 END IF 
 CLOSE(100)

! OPEN(UNIT=200,FILE="kinetic.out",STATUS="REPLACE",ACTION="WRITE")
! DO i=1,nbf
!  WRITE(200,*) (tmat(i,j), j=1,nbf)
! END DO
! CLOSE(200)

RETURN
END SUBROUTINE
