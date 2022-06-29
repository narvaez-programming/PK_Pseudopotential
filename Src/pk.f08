SUBROUTINE pp_no_cnvln(norm,nbf,np,na,ng,grdpts,g,bs,cs,po,eigen,pot)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! COMPUTES THE PSEUDO POTENTIAL FOR   !!
!! THE UNSMOOTHED PSEUDO ORBITAL.      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 USE constants
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: nbf, np, na, ng, grdpts, norm
 DOUBLE PRECISION, INTENT(IN) :: bs(np,10), cs(na+ng,3), po(nbf), eigen
 DOUBLE PRECISION, INTENT(IN) :: g(grdpts,3)
 DOUBLE PRECISION, INTENT(OUT) :: pot(grdpts)
 DOUBLE PRECISION, EXTERNAL :: NORM_PRIM, LAP, NEUMAIER
 DOUBLE PRECISION, ALLOCATABLE :: vctr1(:), vctr2(:)
 DOUBLE PRECISION :: orb, kin_orb, bf1(nbf), bf2(nbf)
 DOUBLE PRECISION :: x, y, z
 DOUBLE PRECISION :: alpha, ax, ay, az, nrm, coeff, xd, yd, zd
 INTEGER :: prim, i, j, k, cntr, atom, xe, ye, ze
 
! OPEN(UNIT=831,FILE='pk_orbital_potential.out',STATUS='REPLACE',ACTION='WRITE')
! 100 FORMAT(E14.6,2X,E14.6,2X,E14.6,2X,E14.6E3,2X,E14.6E3)

 OPEN(UNIT=831,FILE='pk_orbital.out',STATUS='REPLACE',ACTION='WRITE')
 OPEN(UNIT=832,FILE='pk_potential.out',STATUS='REPLACE',ACTION='WRITE')
 100 FORMAT(E14.6E3)

 DO k=1,grdpts
  x = g(k,1); y = g(k,2); z = g(k,3)
  bfloop:DO i=1,nbf
   cntr = 0 
   primitives: DO j=1,np
    IF (IDNINT(bs(j,2)) .NE. i) CYCLE primitives
    cntr = cntr + 1 
   END DO primitives
   ALLOCATE(vctr1(cntr),vctr2(cntr)) 
   prim = 1
   ploop:DO j=1,np 
   IF(IDINT(bs(j,2)) .NE. i) CYCLE ploop     
   alpha = bs(j,7); atom = IDINT(bs(j,1))
   ax = cs(atom,1); ay = cs(atom,2); az = cs(atom,3)
   xe = IDINT(bs(j,4)); ye = IDINT(bs(j,5)); ze = IDINT(bs(j,6)) 
   xd = x - ax; yd = y - ay; zd = z - az
   IF (IDINT(bs(j,3)) .EQ. 0) THEN 
    coeff = bs(j,8)
    IF (norm .EQ. 0) THEN
     nrm = DBLE(1)
    ELSE 
     nrm = NORM_PRIM(alpha,xe,ye,ze)
    END IF
    vctr2(prim) = DBLE(-1)*nrm*bs(j,8)*0.5d0*LAP(alpha,DBLE(xe),DBLE(ye),DBLE(ze),ax,ay,az,x,y,z)
    vctr1(prim) = nrm*coeff*DEXP(-alpha*(xd**2 + yd**2 + zd**2))
   ELSE IF (IDINT(bs(j,3)) .EQ. 1) THEN
    coeff = bs(j,9)*(xd**xe)*(yd**ye)*(zd**ze)
    IF (norm .EQ. 0) THEN
     nrm = DBLE(1)
    ELSE
     nrm = NORM_PRIM(alpha,xe,ye,ze)
    END IF
    vctr2(prim) = DBLE(-1)*nrm*bs(j,9)*0.5d0*LAP(alpha,DBLE(xe),DBLE(ye),DBLE(ze),ax,ay,az,x,y,z)
    vctr1(prim) = nrm*coeff*DEXP(-alpha*(xd**2 + yd**2 + zd**2))
   ELSE IF(IDINT(bs(j,3)) .EQ. 2) THEN 
    coeff = bs(j,10)*(xd**xe)*(yd**ye)*(zd**ze)
    IF (norm .EQ. 0) THEN
     nrm = DBLE(1)
    ELSE
     nrm = NORM_PRIM(alpha,xe,ye,ze)
    END IF
    vctr2(prim) = DBLE(-1)*nrm*bs(j,10)*0.5d0*LAP(alpha,DBLE(xe),DBLE(ye),DBLE(ze),ax,ay,az,x,y,z)
    vctr1(prim) = nrm*coeff*DEXP(-alpha*(xd**2 + yd**2 + zd**2))
   ELSE
    WRITE(*,*) "THE CODE CAN'T HANDLE ORBITALS OF L > 2"
   END IF
   prim = prim + 1
   END DO ploop
   bf1(i) = po(i)*NEUMAIER(prim-1,vctr1)
   bf2(i) = po(i)*NEUMAIER(prim-1,vctr2)
   DEALLOCATE(vctr1,vctr2)
  END DO bfloop
  orb = NEUMAIER(nbf,bf1)
  kin_orb = NEUMAIER(nbf,bf2)
  pot(k) = eigen - kin_orb/orb
!  WRITE(831,100) x, y, z, orb, pot(k)
  WRITE(831,100) DABS(orb)
  WRITE(832,100) pot(k)
 END DO
 CLOSE(831)
 CLOSE(831)

 WRITE(*,*) "COMPUTED THE PSEUDOPOTENTIAL ON A GRID."
 WRITE(*,*) "---------------------------------------"  

RETURN
END SUBROUTINE
