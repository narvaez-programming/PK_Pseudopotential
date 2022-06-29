SUBROUTINE read_mos(nbf,nmos,mos_mat,val)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: nbf, nmos
 DOUBLE PRECISION, INTENT(OUT):: mos_mat(nbf,nmos), val(nbf)
 DOUBLE PRECISION, EXTERNAL :: DDOT
 INTEGER :: io_open, i, j
 REAL :: temp(nbf,nmos)

 OPEN(UNIT=100,FILE='mos.in',STATUS="OLD",ACTION="READ",IOSTAT=io_open)
 IF (io_open .NE. 0) THEN
  STOP "THE FILE mos.in COULD NOT BE OPENED."
 ELSE
  DO i=1,nbf
   READ(100,*) (temp(i,j), j=1,nmos)
  END DO
 END IF
 CLOSE(100)
 mos_mat = DBLE(temp)
 val(:) = mos_mat(:,nmos)
 
! OPEN(UNIT=200,FILE="mos.out",STATUS="REPLACE",ACTION="WRITE")
! DO i=1,nbf
!  WRITE(200,*) (mos_mat(i,j), j=1,nmos)
! END DO
! CLOSE(200)
 
RETURN
END SUBROUTINE read_mos
