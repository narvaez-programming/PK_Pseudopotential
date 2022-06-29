DOUBLE PRECISION FUNCTION neumaier(l,vctr_in)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! COMPUTES THE SUM OF THE ELEMENTS OF AN       !!
!! INPUT VECTOR USING THE NEUMAIER VERSION      !!
!! OF THE KAHAN SUMMATION ALGORITHM.            !!
!!                                              !!
!! VARIABLES:                                   !!
!! l: LENGTH OF THE INPUT VECTOR.               !!
!! vctr_in: 1D ARRAY WITH THE ELEMENTS THAT'LL  !!
!!          BE SUMMED.                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: l
 DOUBLE PRECISION, INTENT(IN) ::  vctr_in(l)
 DOUBLE PRECISION :: c, t, mysum, neu_temp 
 INTEGER :: i 

 mysum = vctr_in(1); c = DBLE(0)
 IF (l .EQ. 1) THEN
  neu_temp = mysum
 ELSE IF (l .EQ. 2) THEN
  neu_temp = mysum + vctr_in(2)
 ELSE
  DO i=2,l
   t = mysum + vctr_in(i)
   IF (DABS(mysum) .GE. DABS(vctr_in(i))) THEN
    c = c + (mysum - t) + vctr_in(i)
   ELSE 
    c = c + (vctr_in(i) - t) + mysum
   END IF
   mysum = t
  END DO
  neu_temp = mysum + c
 END IF

 neumaier = neu_temp
RETURN
END FUNCTION
