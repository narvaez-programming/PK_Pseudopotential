INTEGER FUNCTION fac2(num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! COMPUTES THE DOUBLE FACTORIAL OF AN  !!
!! INTEGER                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 INTEGER, INTENT(IN) :: num
 INTEGER :: num_temp, fac2_temp
  IF (num .LT. 1) then
   fac2 = 1
   RETURN
  ELSE 
   num_temp = num
   fac2_temp = 1
   DO WHILE (1 .LT. num_temp)
    fac2_temp = fac2_temp*num_temp
    num_temp = num_temp - 2
   END DO
   fac2 = fac2_temp
   RETURN
  END IF
END FUNCTION fac2
