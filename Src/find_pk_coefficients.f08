SUBROUTINE find_orb(guess,cnvg,nbf,nmos,mosmat,smat,tmat,pmat,orb)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROVIDES THE COEFFICIENTS FOR THE PSEUDO ORBITAL  !!
!! IN THE BASIS SET.                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 IMPLICIT NONE
 CHARACTER(LEN=1) :: guess
 INTEGER, INTENT(IN) :: nbf, nmos
 DOUBLE PRECISION, INTENT(IN) :: mosmat(nbf,nmos),smat(nbf,nbf),tmat(nbf,nbf),pmat(nbf,nbf)
 REAL, INTENT(IN) :: cnvg
 DOUBLE PRECISION, INTENT(OUT) :: orb(nbf)
 DOUBLE PRECISION, EXTERNAL :: DDOT
 REAL, EXTERNAL :: SDOT
 DOUBLE PRECISION :: rh(nbf,1), old_orb(nbf), new_orb(nbf), ipiv(nbf)
 DOUBLE PRECISION :: temp1(nbf), temp2(nbf,1), sp(nbf,nbf), spt(nbf,nbf), m(nbf,nbf)
 DOUBLE PRECISION :: t_orb
 REAL :: err, diff_orbs(nbf)
 INTEGER :: i, j, info_f, info_ls

 WRITE(*,*) "COMPUTING PSEUDO ORBITAL...." 
 !!!!!!!!!!!!!!!!!!!!
 !! COMPUTE S.LUMO !!
 !!!!!!!!!!!!!!!!!!!!
 rh(:,1) = DBLE(0)
 CALL DGEMM('N','N',nbf,1,nbf,DBLE(1),smat,nbf,mosmat(:,nmos),nbf,DBLE(0),rh,nbf)

 i = 1
 DO 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! ORBITAL FOR INITIAL GUESS !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (i .EQ. 1) THEN
   IF (guess .EQ. 'Y') THEN
    CALL read_orb(nbf,old_orb) 
   ELSE  
    DO j=1,nbf
     old_orb(j) = mosmat(j,nmos)
    END DO
   END IF
  ELSE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! ORBITAL FOR ITERATION > 1 !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   old_orb = new_orb
  END IF

  !!!!!!!!!!!!!!!!!
  !! COMPUTE <T> !!
  !!!!!!!!!!!!!!!!!
  temp1(:) = DBLE(0)
  CALL DGEMV('N',nbf,nbf,DBLE(1),tmat,nbf,old_orb,1,DBLE(0),temp1,1) 
  t_orb = DDOT(nbf,old_orb,1,temp1,1) 
  t_orb = t_orb/DDOT(nbf,old_orb,1,old_orb,1)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! PREPARE THE PK MATRIX !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  sp(:,:) = DBLE(0); spt(:,:) = DBLE(0)
  CALL DGEMM('N','N',nbf,nbf,nbf,DBLE(1),smat,nbf,pmat,nbf,DBLE(0),sp,nbf)
  CALL DGEMM('N','N',nbf,nbf,nbf,DBLE(1),sp,nbf,tmat,nbf,DBLE(0),spt,nbf) 
  spt = spt/t_orb
  m = smat - spt  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! SOLVE THE PK EQUATION !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  temp2 = rh; ipiv(:) = DBLE(0)
  CALL DGETRF(nbf,nbf,m,nbf,ipiv,info_f)
  CALL DGETRS('N',nbf,1,m,nbf,ipiv,temp2,nbf,info_ls)  
  DO j=1,nbf
  new_orb(j) = temp2(j,1)
  END DO
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! COMPUTE < NEW - OLD | NEW - OLD > !! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  diff_orbs = REAL(new_orb - old_orb)
  err = SQRT(SDOT(nbf,diff_orbs,1,diff_orbs,1)) 
  IF (err .LE. cnvg) THEN
   WRITE(*,*) "ITERATION: ", i, "ERROR: ", err
   orb = new_orb
   EXIT
  ELSE
   WRITE(*,*) "ITERATION: ", i, "ERROR: ", err
   i = i + 1
   CONTINUE
  END IF
 END DO

 !!!!!!!!!!!!!!!!!!!
 !! PRINT ORBITAL !!
 !!!!!!!!!!!!!!!!!!!
 OPEN(UNIT=275,FILE="pk_coefficient.out",STATUS="REPLACE",ACTION="WRITE")
 278 FORMAT(E17.8E3)
 DO i=1,nbf
  WRITE(275,278) REAL(orb(i))
 END DO
 CLOSE(275)
 
 WRITE(*,*) "PK COEFFICIENTS FOUND." 
 WRITE(*,*) "----------------------"  

RETURN
END SUBROUTINE
