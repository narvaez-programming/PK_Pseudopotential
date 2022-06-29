SUBROUTINE taper_atom(na,taper_cut_off,taper_range,charges,x,y,z,molc,ueff,output) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CURRENTLY THE COM CAN ONLY BE CALCULATED   !!
!! FOR AN ATOMIC PSEUDOPOTENTIAL              !!
!! MAY 20, 2019                               !!
!! ONLY TAPERS ATOMIC POTENTIALS RIGHT NOW.   !!
!! JUNE 3, 2019                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE constants
IMPLICIT NONE
INTEGER, INTENT(IN) :: na
DOUBLE PRECISION, INTENT(IN) :: charges(na,1), molc(na,3)
DOUBLE PRECISION, INTENT(IN) :: taper_cut_off, taper_range, x, y, z, ueff
DOUBLE PRECISION, INTENT(OUT) :: output
DOUBLE PRECISION :: r, r2, rinv, tswit
DOUBLE PRECISION :: au_cut_off, scut_off, t2mi, t3m
DOUBLE PRECISION :: temporary, taper
DOUBLE PRECISION :: xa, ya, za, ra, ra2, ra_inv
INTEGER :: i

 au_cut_off = taper_cut_off*atb  
 CALL PREP_TAPERING(taper_cut_off,taper_range,scut_off,t2mi,t3m)
 temporary = 0d0; taper = 0d0

 r2 = x**2d0 + y**2d0 + z**2d0; r = DSQRT(r2); 

 IF (r .GE. au_cut_off) THEN
  DO i=1,na
   xa = x - molc(i,1); ya = y - molc(i,2); za = z - molc(i,3)
   ra2 = xa**2d0 + ya**2d0 + za**2d0; ra = DSQRT(ra2); ra_inv = 1d0/ra 
   temporary = temporary - charges(i,1)*ra_inv
  END DO 
 ELSE
  CALL STEINHAUSER(au_cut_off,scut_off,t3m,t2mi,r,tswit)
  DO i=1,na
   xa = x - molc(i,1); ya = y - molc(i,2); za = z - molc(i,3)
   ra2 = xa**2d0 + ya**2d0 + za**2d0; ra = DSQRT(ra2); ra_inv = 1d0/ra 
   taper = taper + (1.0d0 - tswit)*charges(i,1)*ra_inv 
  END DO 
  temporary = ueff*tswit - taper
 END IF

 output = temporary
RETURN
END SUBROUTINE taper_atom
