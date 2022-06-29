INTEGER FUNCTION fac(num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! COMPUTES THE FACTORIAL OF AN INTEGER !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 INTEGER, INTENT(IN) :: num
 INTEGER :: num_temp, fac_temp 
 num_temp = num
 fac_temp = 1
 DO WHILE (1 .LT. num_temp)
  fac_temp = fac_temp*num_temp
  num_temp = num_temp - 1
 END DO
 fac = fac_temp
 RETURN
END FUNCTION fac
