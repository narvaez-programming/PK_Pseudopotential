PROGRAM find_pseudo_orb
 USE setup_variables
 IMPLICIT NONE

 CALL read_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PREPARE THE SPATIAL GRID POINTS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 IF (grid_type .EQ. 'R') THEN
  grid_pts = xp*yp*zp
  ALLOCATE( grd(grid_pts,3), pot(grid_pts), tapered_pot(grid_pts) )
  CALL grid(grid_print,grid_shift,lx,rx,ly,ry,lz,rz,xp,yp,zp,grid_pts,grd)
 ELSE IF (grid_type .EQ. 'C') THEN
  grid_pts = pts_1d**3
  ALLOCATE( grd(grid_pts,3), pot(grid_pts), tapered_pot(grid_pts) )
  CALL grid_cube(grid_print,grid_shift,pts_1d,spacng,grid_pts,grd)
 END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ALLOCATE NECESSATY MATRICES !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CALL prepare_arrays 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! READ BASIS SET AND MOLECULAR COORDINATES !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CALL read_basis(tbf,tp,b) 
 CALL read_coords(ta+tg,c)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! READ KINETIC, OVERLAP, AND MOs !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 IF (only_pot .EQ. 'N') THEN
  CALL read_t(tbf,t)
  CALL read_s(tbf,s) 
  CALL read_mos(tbf,tmos,mos,v)
  CALL build_projection(tbf,tmos,mos,p)
 ELSE
   CONTINUE
 END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SOLVE THE PK PSEUDO ORBITAL !!
!! EQUATION                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 IF (only_pot .EQ. 'Y') THEN
  CALL read_orb(tbf,o)
 ELSE
  CALL find_orb(orb_guess,cnvrg_thresh,tbf,tmos,mos,s,t,p,o)
 END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PRINT SOME CORE MOS  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!
 IF (core_mos_print .EQ. 0) THEN
  CONTINUE
 ELSE
  CALL print_mos(nrmlz,tbf,tp,ta,tg,grid_pts,grd,b,c,core_mos_print,tmos,mos) 
 END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! COMPUTE PSEUDO POTENTIAL ON A GRID !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 IF (taper_call .EQ. 'Y') THEN
  CALL pp_cnvln_tapered(taper_cut_off,taper_range,cnvltn_alpha,nrmlz,tbf,tp,ta,tg,grid_pts,grd,b,c,o,eval,pot)
 ELSE 
  IF (cnvltn .EQ. 'Y') THEN
   CALL pp_cnvln(cnvltn_alpha,nrmlz,tbf,tp,ta,tg,grid_pts,grd,b,c,o,eval,pot)
  ELSE 
   CALL pp_no_cnvln(nrmlz,tbf,tp,ta,tg,grid_pts,grd,b,c,o,eval,pot)
  END IF
 END IF

END PROGRAM
