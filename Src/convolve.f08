SUBROUTINE convolve(alpha,beta,xe,ye,ze,output)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! alpha := EXPONENT OF THE BASIS FUNCTION                !!
!! beta := EXPONENT OF THE CONVOLUTION GAUSSIAN           !!
!! xe := POWER OF THE x POLYNOMAL FOR THE BASIS FUNCTION  !!
!! output := (1) NEW EXPONENT                             !!
!!           (2-3) CONSTANT                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 USE constants
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: xe, ye, ze
 DOUBLE PRECISION, INTENT(IN) :: alpha, beta
 DOUBLE PRECISION, INTENT(OUT) :: output(3)
 INTEGER :: l
 DOUBLE PRECISION :: expo_sum, nrm
  
 l = xe + ye + ze
 expo_sum = alpha+beta
 output(1) = (alpha*beta)/expo_sum
 IF (l .EQ. 0) THEN
  output(2) = (pi/expo_sum)**1.5d0
  output(3) = DBLE(0)
 ELSE IF (l .EQ. 1) THEN 
  nrm = (pi/expo_sum)**1.5d0
  output(2) = nrm*beta/expo_sum
  output(3) = DBLE(0)
 ELSE IF (l .EQ. 2) THEN
  IF (xe .EQ. 2 .OR. ye .EQ. 2 .OR. ze .EQ. 2) THEN
   nrm = pi/expo_sum
   nrm = nrm*DBLE(0.5)*DSQRT(pi)/DSQRT(expo_sum**5) 
   output(2) = nrm*expo_sum
   output(3) = nrm*DBLE(2)*(beta**2d0)
  ELSE 
   nrm = DSQRT(pi/expo_sum)
   output(2) = nrm*(beta**2d0)*pi/(expo_sum**3d0) 
   output(3) = DBLE(0)
  END IF
 END IF

RETURN  
END SUBROUTINE
