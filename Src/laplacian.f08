DOUBLE PRECISION FUNCTION LAP(alpha, ax, ay, az, a_x, a_y, a_z, x, y, z)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! alpha = ORBTIAL EXPONENT                                           !!
  !! ai = POWER OF X, Y, or Z FOR i                                     !!
  !! a_i = X, Y, or Z CENTER OF ORBITAL i                               !!
  !! x,y,z = POINT USED TO EVALUATE THE LAP                             !!
  !! vctr(1,2,3) = LAP OF THE X, Y, Z COMPONENT RESPECTIVELY            !!
  !! vctr1 = COMPONENTS OF THE X COMPONENT LAPLACIAN                    !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION, INTENT(IN) :: alpha, ax, ay, az
  DOUBLE PRECISION, INTENT(IN) :: a_x, a_y, a_z, x, y, z
  DOUBLE PRECISION, EXTERNAL :: NEUMAIER 
  DOUBLE PRECISION :: lapx, lapy, lapz
  DOUBLE PRECISION :: vctr1(4), vctr2(4), vctr3(4), vctr(3)
  
  IF (ax .LE. DBLE(1)) THEN 
   vctr1(1) = DBLE(0)
  ELSE
   vctr1(1) = ax*(ax-1d0)*((x-a_x)**(ax-2d0))
  END IF
  vctr1(2) = DBLE(-1)*2d0*alpha*(ax+1d0)*((x-a_x)**ax)
  IF (ax .EQ. DBLE(0)) THEN
   vctr1(3) = DBLE(0)
  ELSE
   vctr1(3) = DBLE(-1)*2d0*alpha*ax*((x-a_x)**ax)
  END IF
  vctr1(4) = 4d0*(alpha**2d0)*((x-a_x)**(ax+2d0)) 
  vctr(1) = ((y-a_y)**(ay))*((z-a_z)**(az))*NEUMAIER(4,vctr1)
  vctr(1) = vctr(1)*DEXP(-alpha*((x-a_x)**2d0 + (y-a_y)**2d0 + (z-a_z)**2d0))

  IF (ay .LE. DBLE(1)) THEN 
   vctr2(1) = DBLE(0)
  ELSE
   vctr2(1) = ay*(ay-1d0)*((y-a_y)**(ay-2d0))
  END IF
  vctr2(2) = DBLE(-1)*2d0*alpha*(ay+1d0)*((y-a_y)**ay)
  IF (ay .EQ. DBLE(0)) THEN
   vctr2(3) = DBLE(0)
  ELSE
   vctr2(3) = DBLE(-1)*2d0*alpha*ay*((y-a_y)**ay)
  END IF
  vctr2(4) = 4d0*(alpha**2d0)*((y-a_y)**(ay+2d0)) 
  vctr(2) = ((x-a_x)**(ax))*((z-a_z)**(az))*NEUMAIER(4,vctr2)
  vctr(2) = vctr(2)*DEXP(-alpha*((x-a_x)**2d0 + (y-a_y)**2d0 + (z-a_z)**2d0))

  IF (az .LE. DBLE(1)) THEN 
   vctr3(1) = DBLE(0)
  ELSE
   vctr3(1) = az*(az-1d0)*((z-a_z)**(az-2d0))
  END IF
  vctr3(2) = DBLE(-1)*2d0*alpha*(az+1d0)*((z-a_z)**az)
  IF (az .EQ. DBLE(0)) THEN
   vctr3(3) = DBLE(0)
  ELSE
   vctr3(3) = DBLE(-1)*2d0*alpha*az*((z-a_z)**az)
  END IF
  vctr3(4) = 4d0*(alpha**2d0)*((z-a_z)**(az+2d0)) 
  vctr(3) = ((y-a_y)**(ay))*((x-a_x)**(ax))*NEUMAIER(4,vctr3)
  vctr(3) = vctr(3)*DEXP(-alpha*((x-a_x)**2d0 + (y-a_y)**2d0 + (z-a_z)**2d0)) 

  lap = NEUMAIER(3,vctr)
RETURN
END FUNCTION LAP
