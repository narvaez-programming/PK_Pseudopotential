DOUBLE PRECISION FUNCTION norm_prim(a,l,m,n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! THIS WILL NORMALIZE A GAUSSIAN PRIMITIVE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 USE constants
 IMPLICIT NONE
 DOUBLE PRECISION, INTENT(IN) :: a
 INTEGER, INTENT(IN) :: l, m, n
 INTEGER, EXTERNAL :: FAC2
 DOUBLE PRECISION :: temp, num, dem

 num = DBLE(4)*a
 num = num**(l+m+n)
 num = DSQRT(num)
 dem = DBLE(FAC2(2*l - 1)*FAC2(2*m - 1)*FAC2(2*n - 1))
 dem = DSQRT(dem)
 temp = DBLE(2)*a/pi
 temp = temp**DBLE(0.75)
 temp = temp*num/dem
 norm_prim = temp

RETURN
END FUNCTION

