MODULE setup_variables
 IMPLICIT NONE
 INTEGER :: tbf, tp, tmos, ta, tg, nrmlz
 CHARACTER(LEN=1) :: grid_type, grid_print 
 INTEGER :: core_mos_print
 INTEGER :: smooth_pts, cnvltn_pts, grid_pts, pts_1d
 REAL  :: cnvrg_thresh
 DOUBLE PRECISION :: eval, grid_shift
 CHARACTER(LEN=1)  :: only_pot, orb_guess, cnvltn, taper_call
 DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: grd, b, c
 DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: t, s, mos, p
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: o, v
 DOUBLE PRECISION :: cnvltn_alpha, taper_cut_off, taper_range
 DOUBLE PRECISION, ALLOCATABLE :: pk_orb(:), smooth(:), gauss_orb(:)
 DOUBLE PRECISION :: lx, rx, ly, ry, lz, rz, spacng
 DOUBLE PRECISION, ALLOCATABLE :: pot(:), tapered_pot(:) 
INTEGER :: xp, yp, zp

CONTAINS
 SUBROUTINE read_parameters 
  IMPLICIT NONE
  NAMELIST/param/tbf,tp,nrmlz,tmos,ta,tg,cnvrg_thresh,&
&cnvltn,cnvltn_alpha,taper_call,taper_cut_off,taper_range,core_mos_print,orb_guess,only_pot,eval,&
&grid_type,grid_print,grid_shift,lx,rx,ly,ry,lz,rz,xp,yp,zp,&
&pts_1d,spacng
  OPEN(UNIT=12,FILE="param.in",STATUS="OLD",ACTION="READ")
  READ(12,param)
  CLOSE(12)

 RETURN 
 END SUBROUTINE read_parameters

 SUBROUTINE prepare_arrays
  IMPLICIT NONE
  ALLOCATE( b(tp,10), c(ta+tg,3) )
  ALLOCATE( t(tbf,tbf), s(tbf,tbf), mos(tbf,tmos), p(tbf,tbf) )
  ALLOCATE( o(tbf), v(tbf) )

 RETURN
 END SUBROUTINE prepare_arrays

END MODULE

