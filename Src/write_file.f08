SUBROUTINE write_file(filename,rows,columns,input)
IMPLICIT NONE
CHARACTER(LEN=*) :: filename
INTEGER, INTENT(IN) :: rows, columns
DOUBLE PRECISION, INTENT(IN) :: input(rows,columns)
INTEGER :: i, j, io_open

OPEN(UNIT=100,FILE=filename,STATUS="REPLACE",ACTION="WRITE",IOSTAT=io_open,ERR=41)
41 IF (io_open .NE. 0) THEN
 WRITE(*,*) 'CAN NOT OPEN:', filename
 STOP
ELSE 
 DO i=1,rows
  WRITE(100,*) (input(i,j), j=1,columns)
 END DO
 CLOSE(100) 
END IF

RETURN
END SUBROUTINE
