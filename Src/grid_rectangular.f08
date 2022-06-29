SUBROUTINE grid(grid_print,shft,xl,xr,yl,yr,zl,zr,xpts,ypts,zpts,tpts,g)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! L_i = (N_i -1)*spacing_i                             !!
!! INPUT GRID INFORMATION IN ANGSTROMS.                 !!
!! THE OUTPUT WILL BE PROVIDED IN BOHRS.                !!
!! THE GRID IS CENTERED AROUND 0.                       !! 
!! shft MOVES THE CENTER OF THE GRID AWAY FROM 0        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 USE constants
 IMPLICIT NONE
 CHARACTER(LEN=1), INTENT(IN) :: grid_print
 DOUBLE PRECISION, INTENT(IN) :: shft, xl, xr, yl, yr, zl, zr
 INTEGER, INTENT(IN) :: xpts, ypts, zpts, tpts
 DOUBLE PRECISION, INTENT(OUT) :: g(tpts,3)
 INTEGER :: n, i, j, k, pts
 DOUBLE PRECISION :: xd, yd, zd, x, y, z

 IF (xpts .EQ. 1) THEN 
  xd = DBLE(0)
 ELSE 
  xd = (xr - xl)/DBLE(xpts-1)
 END IF
 IF (ypts .EQ. 1) THEN
  yd = DBLE(0)
 ELSE
  yd = (yr - yl)/DBLE(ypts-1)
 END IF
 IF (zpts .EQ. 1) THEN
  zd = DBLE(0)
 ELSE 
  zd = (zr - zl)/DBLE(zpts-1)
 END IF

 n = 1
 DO k=1,zpts
  DO j=1,ypts
   DO i=1,xpts
    x = (DBLE(i-1)*xd + xl + shft)*atb
    y = (DBLE(j-1)*yd + yl + shft)*atb
    z = (DBLE(k-1)*zd + zl + shft)*atb
    g(n,1) = x; g(n,2) = y; g(n,3) = z
    n = n + 1
   END DO
  END DO
 END DO

 IF (grid_print .EQ. 'Y') THEN
  OPEN(UNIT=356,FILE='grid.out',ACTION='WRITE',STATUS='REPLACE')
  87 FORMAT(3E14.6)
  DO i=1,tpts 
   WRITE(356,87) ( g(i,j), j=1,3 )
  END DO
  CLOSE(356)
 ELSE
  CONTINUE
 END IF

 243 FORMAT(' ',A1,3X,F6.2,3X,F6.2,3X,I5,3X,F9.4)
 WRITE(*,*) "----------------------------------------------------------------------" 
 WRITE(*,*) "GRID INFORMATION (AXIS, LEFT LIMIT, RIGHT LIMIT, GRID POINTS, SPACING)"
 WRITE(*,243) "X", xl, xr, xpts, xd 
 WRITE(*,243) "Y", yl, yr, ypts, yd
 WRITE(*,243) "Z", zl, zr, zpts, zd
 WRITE(*,*) "----------------------------------------------------------------------" 

RETURN
END SUBROUTINE
