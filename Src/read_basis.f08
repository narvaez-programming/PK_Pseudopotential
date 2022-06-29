SUBROUTINE read_basis(nbf,np,bfmat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! YOU MUST PROVIDE YOUR GTO-TYPE BASIS SET IN  !!
!! PURE CARTESIAN FORM. WE FOLLOW THE FORMAT    !!
!! USED BY GAUSSIAN.                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 IMPLICIT NONE
 INTEGER, INTENT(IN) :: nbf,np
 DOUBLE PRECISION, INTENT(OUT) :: bfmat(np,10)
 INTEGER :: io_open, io_read, atom, tag1, nprims, bf_start, bf_end, prim_type, i, j
 CHARACTER(LEN=4) :: bf_type, atom_end
 INTEGER :: prim_num, atom_num, bf_cnt, prim_cnt, cntr, x_expo, y_expo, z_expo
 DOUBLE PRECISION :: expo, coeff, s_coeff, p_coeff
 REAL :: sclr, sclr2

 OPEN(UNIT=100,FILE="basis.in",STATUS="OLD",ACTION="READ",IOSTAT=io_open) 
 READ(100,*,IOSTAT=io_read) atom, tag1
 READ(100,*) bf_type, nprims, sclr, sclr2
 IF (io_open .NE. 0) THEN
  STOP "COULD NOT OPEN basis.in FILE"
 ELSE 
  prim_num = 1; bf_start = 1
  read_file:DO WHILE (io_read .EQ. 0) 
   atom_loop:DO WHILE (bf_type .NE. "****")
    IF (bf_type .EQ. "S   ") THEN
     bf_end = bf_start
     DO bf_cnt = bf_start, bf_end
      DO prim_cnt = 1, nprims     
       READ(100,*) expo, coeff
       bfmat(prim_num,1) = DBLE(atom); bfmat(prim_num,2) = DBLE(bf_cnt)
       bfmat(prim_num,3) = DBLE(0); bfmat(prim_num,4) = DBLE(0)
       bfmat(prim_num,5) = DBLE(0); bfmat(prim_num,6) = DBLE(0)
       bfmat(prim_num,7) = expo; bfmat(prim_num,8) = coeff
       bfmat(prim_num,9) = DBLE(0); bfmat(prim_num,10) = DBLE(0)
       prim_num = prim_num + 1
      END DO
      IF (bf_cnt .LT. bf_end) THEN
       DO prim_cnt = 1, nprims
        BACKSPACE(100)
       END DO
      ELSE 
       bf_start = bf_end + 1
       CONTINUE
      END IF
     END DO 
    ELSE IF (bf_type .EQ. "P   ") THEN
     cntr = 1; bf_end = bf_start + 2
     DO bf_cnt = bf_start, bf_end     
      IF (cntr .EQ. 1) THEN 
       x_expo = 1; y_expo = 0; z_expo = 0
      ELSE IF (cntr .EQ. 2) THEN
       x_expo = 0; y_expo = 1; z_expo = 0
      ELSE 
       x_expo = 0; y_expo = 0; z_expo = 1
      END IF  
      DO prim_cnt = 1, nprims
       READ(100,*) expo, coeff
       bfmat(prim_num,1) = DBLE(atom); bfmat(prim_num,2) = DBLE(bf_cnt)    
       bfmat(prim_num,3) = DBLE(1); bfmat(prim_num,4) = DBLE(x_expo)
       bfmat(prim_num,5) = DBLE(y_expo); bfmat(prim_num,6) = DBLE(z_expo)
       bfmat(prim_num,7) = expo; bfmat(prim_num,8) = DBLE(0)
       bfmat(prim_num,9) = coeff; bfmat(prim_num,10) = DBLE(0)
       prim_num = prim_num + 1
      END DO
      cntr = cntr + 1
      IF (bf_cnt .LT. bf_end) THEN
       DO prim_cnt = 1, nprims
        BACKSPACE(100)
       END DO
      ELSE
       bf_start = bf_end + 1
       CONTINUE
      END IF
     END DO
    ELSE IF (bf_type .EQ. "D   ") THEN
     cntr = 1; bf_end = bf_start + 5
     DO bf_cnt = bf_start, bf_end
      IF (cntr .EQ. 1) THEN
       x_expo = 2; y_expo = 0; z_expo = 0
      ELSE IF (cntr .EQ. 2) THEN
       x_expo = 0; y_expo = 2; z_expo = 0
      ELSE IF (cntr .EQ. 3) THEN
       x_expo = 0; y_expo = 0; z_expo = 2
      ELSE IF (cntr .EQ. 4) THEN
       x_expo = 1; y_expo = 1; z_expo = 0
      ELSE IF (cntr .EQ. 5) THEN
       x_expo = 1; y_expo = 0; z_expo = 1
      ELSE 
       x_expo = 0; y_expo = 1; z_expo = 1
      END IF
      DO prim_cnt = 1, nprims                                             
       READ(100,*) expo, coeff                                            
       bfmat(prim_num,1) = DBLE(atom); bfmat(prim_num,2) = DBLE(bf_cnt)   
       bfmat(prim_num,3) = DBLE(2); bfmat(prim_num,4) = DBLE(x_expo)      
       bfmat(prim_num,5) = DBLE(y_expo); bfmat(prim_num,6) = DBLE(z_expo) 
       bfmat(prim_num,7) = expo; bfmat(prim_num,8) = DBLE(0)
       bfmat(prim_num,9) = DBLE(0); bfmat(prim_num,10) = coeff           
       prim_num = prim_num + 1
      END DO
      cntr = cntr + 1
      IF (bf_cnt .LT. bf_end) THEN
       DO prim_cnt = 1, nprims
        BACKSPACE(100)
       END DO
      ELSE
       bf_start = bf_end + 1
       CONTINUE
      END IF
     END DO
    ELSE IF (bf_type .EQ. "SP  ") THEN
     cntr = 1; bf_end = bf_start + 3
     DO bf_cnt = bf_start, bf_end
      IF (cntr .EQ. 1) THEN
       prim_type = 0
       x_expo = 0; y_expo = 0; z_expo = 0
      ELSE IF (cntr .EQ. 2) THEN
       prim_type = 1
       x_expo = 1; y_expo = 0; z_expo = 0
      ELSE IF (cntr .EQ. 3) THEN
       prim_type = 1
       x_expo = 0; y_expo = 1; z_expo = 0
      ELSE 
       prim_type = 1
       x_expo = 0; y_expo = 0; z_expo = 1
      END IF
      DO prim_cnt = 1, nprims
       READ(100,*) expo, s_coeff, p_coeff
       bfmat(prim_num,1) = DBLE(atom); bfmat(prim_num,2) = DBLE(bf_cnt)   
       bfmat(prim_num,3) = DBLE(prim_type); bfmat(prim_num,4) = DBLE(x_expo)      
       bfmat(prim_num,5) = DBLE(y_expo); bfmat(prim_num,6) = DBLE(z_expo) 
       bfmat(prim_num,7) = expo; bfmat(prim_num,8) = s_coeff
       bfmat(prim_num,9) = p_coeff; bfmat(prim_num,10) = DBLE(0) 
       prim_num = prim_num + 1
      END DO
      cntr = cntr + 1
      IF (bf_cnt .LT. bf_end) THEN
       DO prim_cnt = 1, nprims
        BACKSPACE(100)
       END DO
      ELSE
       bf_start = bf_end + 1
       CONTINUE
      END IF
     END DO                            
    END IF
     READ(100,*) bf_type
     IF (bf_type .EQ. "****") THEN
      READ(100,*,IOSTAT=io_read) atom, tag1 
      IF (io_read .LT. 0) THEN
       EXIT read_file
      ELSE IF (io_read .GT. 0) THEN
       WRITE(*,*) "FAILED TO READ basis.in FILE"
       EXIT read_file
      ELSE 
       READ(100,*) bf_type, nprims, sclr, sclr2
       EXIT atom_loop
      END IF
     ELSE
      BACKSPACE(100)
      READ(100,*) bf_type, nprims, sclr, sclr2
     END IF
   END DO atom_loop
  END DO read_file
 END IF
 CLOSE(100)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UN-COMMENT IF YOU WANT TO PRINT OUT !!
!! WHAT WAS STORED AS THE BASIS SET.   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! OPEN(UNIT=200,FILE="basis.out",STATUS="REPLACE",ACTION="WRITE")
! DO i = 1, np
!  WRITE(200,*) (bfmat(i,j), j=1,10)
! END DO
! CLOSE(200)

RETURN
END SUBROUTINE
