SUBROUTINE build_projection(nbf,nmos,mosmat,pmat)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! COMPUTES |i><i| FOR EVERY MO AND SUMMES THE    !!
 !! RESULTS INTO ONE MATRIX.                       !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! nbf = # OF BASIS FUNCTIONS                     !!
 !! nmos = # OF MOS                                !!
 !! mosmat = MATRIX COMPOSED OF ONE MO PER COLUMN  !!
 !! pmat = SUM(|i><i|) FOR i MOS                   !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: nbf, nmos
 DOUBLE PRECISION, INTENT(IN) :: mosmat(nbf,nmos)
 DOUBLE PRECISION, INTENT(OUT) :: pmat(nbf,nbf)
 INTEGER :: i, j, k
 DOUBLE PRECISION :: temp(nbf,nbf)
 DOUBLE PRECISION :: c_mat(nbf,nbf), y_mat(nbf,nbf), t_mat(nbf,nbf)
  


 CALL DGEMM('N','T',nbf,nbf,1,DBLE(1),mosmat(:,1),nbf,mosmat(:,1),nbf,DBLE(0),pmat,nbf)
 c_mat(:,:) = DBLE(0)
 DO i=2,nmos-1
  CALL DGEMM('N','T',nbf,nbf,1,DBLE(1),mosmat(:,i),nbf,mosmat(:,i),nbf,DBLE(0),temp,nbf)
  DO j=1,nbf
   DO k=1,nbf
    t_mat(j,k) = pmat(j,k) + temp(j,k)
    IF (DABS(pmat(j,k)) .GE. DABS(temp(j,k))) THEN
     c_mat(j,k) = c_mat(j,k) + (pmat(j,k) - t_mat(j,k)) + temp(j,k)
    ELSE 
     c_mat(j,k) = c_mat(j,k) + (temp(j,k) - t_mat(j,k)) + pmat(j,k)
    END IF
    pmat(j,k) = t_mat(j,k)
   END DO
  END DO
 END DO
 pmat = pmat + c_mat

! OPEN(UNIT=100,FILE="projection.out",STATUS="REPLACE",ACTION="WRITE")
! DO i=1,nbf
!  WRITE(100,*) (pmat(i,j), j=1,nbf)
! END DO
! CLOSE(100)

RETURN
END SUBROUTINE
