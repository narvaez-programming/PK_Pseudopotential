SUBROUTINE grid_cube(grid_print,shft,pts_1d,spacng,tpts,g)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! COMPUTES A GRID USING THE FORMULA            !!
!! x(i) = (i*spacing - L/2)                     !!
!! L = (N -1)*spacing                           !!
!! PROVIDE spacng IN ANGSTROMS                  !!
!! ONLY WORKS FOR CUBIC GRIDS                   !!
!! GRID IS CENTERED AROUND 0 BY DEFAULT         !!
!! shft LETS YOU ALTER THE CENTER OF THE GRID   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 USE constants
 IMPLICIT NONE
 CHARACTER(LEN=1), INTENT(IN) :: grid_print
 INTEGER, INTENT(IN) :: pts_1d, tpts
 DOUBLE PRECISION, INTENT(IN) :: shft, spacng
 DOUBLE PRECISION, INTENT(OUT) :: g(tpts,3)
 INTEGER :: i, j, k, n
 DOUBLE PRECISION :: ln

 ln = DBLE(pts_1d - 1)/DBLE(2)
 n = 1
 DO k=1,pts_1d
  DO j=1,pts_1d
   DO i=1,pts_1d
    g(n,1) = (DBLE(i - 1) - ln)*spacng*atb + shft*atb 
    g(n,2) = (DBLE(j - 1) - ln)*spacng*atb + shft*atb
    g(n,3) = (DBLE(k - 1) - ln)*spacng*atb + shft*atb
    n = n + 1
   END DO
  END DO
 END DO

 IF (grid_print .EQ. 'Y') THEN
  OPEN(UNIT=356,FILE='grid.out',ACTION='WRITE',STATUS='REPLACE')
  87 FORMAT(3E14.6)
  DO i=1,tpts
   WRITE(356,87) (g(i,j), j=1,3) 
  END DO
  CLOSE(356)
 ELSE 
  CONTINUE
 END IF

 
 WRITE(*,*) '------------------------'
 WRITE(*,*) 'THE GRID HAS BEEN BUILT.'
 WRITE(*,*) '------------------------'

RETURN
END SUBROUTINE
