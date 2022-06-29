SUBROUTINE read_s(nbf,smat)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: nbf
 DOUBLE PRECISION, INTENT(OUT) :: smat(nbf,nbf)
 INTEGER :: i, j, io_open 

 OPEN(UNIT=100,FILE="overlap.in",STATUS="OLD",ACTION="READ",IOSTAT=io_open)
 IF (io_open .NE. 0) THEN
  STOP "COULD NOT READ overlap.in FILE"
 ELSE 
  DO i=1,nbf
  READ(100,*) (smat(i,j), j=1,nbf) 
  END DO
 END IF 
 CLOSE(100)

! OPEN(UNIT=200,FILE="overlap.out",STATUS="REPLACE",ACTION="WRITE")
! DO i=1,nbf
!  WRITE(200,*) (smat(i,j), j=1,nbf)
! END DO
! CLOSE(200)
 
RETURN
END SUBROUTINE
