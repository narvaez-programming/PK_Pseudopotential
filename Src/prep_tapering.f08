SUBROUTINE prep_tapering(cut_off,stap,scut_off,t2mi,t3m) 
 USE constants
 IMPLICIT NONE
 DOUBLE PRECISION, INTENT(IN) :: cut_off, stap
 DOUBLE PRECISION, INTENT(OUT) :: scut_off, t2mi, t3m
 DOUBLE PRECISION :: au_cut_off, au_stap 

 au_cut_off = cut_off*atb; au_stap = stap*atb
 scut_off = au_cut_off - au_stap
 t2mi = 2.0d0/((au_cut_off - scut_off)**3d0)
 t3m = 0.5d0*(3.0d0*au_cut_off - scut_off)

RETURN
END SUBROUTINE prep_tapering


