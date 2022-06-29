SUBROUTINE steinhauser(cut_off,scut_off,t3m,t2mi,r,tswit) 
 IMPLICIT NONE
 DOUBLE PRECISION, INTENT(IN) :: cut_off, scut_off
 DOUBLE PRECISION, INTENT(IN) :: t3m, t2mi, r
 DOUBLE PRECISION, INTENT(OUT) :: tswit
 DOUBLE PRECISION :: difft, diffr

 difft = r - scut_off
 diffr = r - cut_off
  
 tswit = 1.0d0 - (difft**2)*(t3m - r)*t2mi
 tswit = MERGE(tswit,1.0d0,difft .GE. 0.0d0)
 tswit = MERGE(0.0d0,tswit,diffr .GE. 0.0d0)

RETURN
END SUBROUTINE
