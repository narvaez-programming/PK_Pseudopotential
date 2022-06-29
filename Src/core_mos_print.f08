SUBROUTINE print_mos(norm,nbf,np,na,ng,grdpts,g,bs,cs,nmos,tmos,mos)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GONNA PRINT SOME CORE MOS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 USE constants
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: nbf, np, na, ng, grdpts, norm
 INTEGER, INTENT(IN) :: nmos, tmos
 DOUBLE PRECISION, INTENT(IN) :: bs(np,10), cs(na+ng,3), mos(nbf,tmos)
 DOUBLE PRECISION, INTENT(IN) :: g(grdpts,3)
 DOUBLE PRECISION, EXTERNAL :: NORM_PRIM, LAP, NEUMAIER
 DOUBLE PRECISION, ALLOCATABLE :: vctr1(:)
 DOUBLE PRECISION :: core_mo(1,nmos), tcore(nbf), core(nbf,tmos)
 DOUBLE PRECISION :: x, y, z
 DOUBLE PRECISION :: alpha, ax, ay, az, nrm, coeff, xd, yd, zd
 INTEGER :: cloop, prim, i, j, k, cntr, atom, xe, ye, ze
 
 core(:,:) = 0d0
 OPEN(UNIT=712,FILE='core_orbitals.out',STATUS='REPLACE',ACTION='WRITE')
 100 FORMAT(20E20.6E3) 

 DO k=1,grdpts
  x = g(k,1); y = g(k,2); z = g(k,3)
  bfloop:DO i=1,nbf
   cntr = 0 
   primitives: DO j=1,np
    IF (IDNINT(bs(j,2)) .NE. i) CYCLE primitives
    cntr = cntr +1 
   END DO primitives
   ALLOCATE(vctr1(cntr)) 
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
    vctr1(prim) = nrm*coeff*DEXP(-alpha*(xd**2 + yd**2 + zd**2))
   ELSE IF (IDINT(bs(j,3)) .EQ. 1) THEN
    coeff = bs(j,9)*(xd**xe)*(yd**ye)*(zd**ze)
    IF (norm .EQ. 0) THEN
     nrm = DBLE(1)
    ELSE
     nrm = NORM_PRIM(alpha,xe,ye,ze)
    END IF
    vctr1(prim) = nrm*coeff*DEXP(-alpha*(xd**2 + yd**2 + zd**2))
   ELSE IF(IDINT(bs(j,3)) .EQ. 2) THEN 
    coeff = bs(j,10)*(xd**xe)*(yd**ye)*(zd**ze)
    IF (norm .EQ. 0) THEN
     nrm = DBLE(1)
    ELSE
     nrm = NORM_PRIM(alpha,xe,ye,ze)
    END IF
    vctr1(prim) = nrm*coeff*DEXP(-alpha*(xd**2 + yd**2 + zd**2))
   ELSE
    WRITE(*,*) "THE CODE CAN'T HANDLE ORBITALS OF L > 2"
   END IF
   prim = prim + 1
   END DO ploop
   DO cloop=1,nmos
    core(i,cloop) = mos(i,cloop)*NEUMAIER(prim-1,vctr1)
   END DO
   DEALLOCATE(vctr1)
  END DO bfloop
  DO cloop=1,nmos
   tcore(:) = core(:,cloop)  
   core_mo(1,cloop) = NEUMAIER(nbf,tcore)
  END DO 
  WRITE(712,100) (core_mo(1,cloop), cloop=1,nmos)
 END DO

 CLOSE(712)

 WRITE(*,*) "COMPUTED CORE MOs."
 WRITE(*,*) "---------------------------------------"  

RETURN
END SUBROUTINE
