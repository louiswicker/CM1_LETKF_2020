 subroutine cressman(x, y, ob, xg, yg, roi, anal, nobs, nx, ny)

   implicit none
 
   real(8), intent(in),  dimension(nobs) :: x
   real(8), intent(in),  dimension(nobs) :: y
   real(8), intent(in),  dimension(nobs) :: ob
   real(8), intent(in),  dimension(nx)   :: xg
   real(8), intent(in),  dimension(ny)   :: yg
   real(8), intent(out), dimension(nx,ny) :: anal
   real,    intent(in) :: roi
   integer, intent(in) :: nobs, nx, ny
 
   integer n, i, j
   real(kind=8) dis, R2, w_sum, top, wk, rk2
   real, parameter :: hsp = 1.33

   logical, parameter :: debug = .false.
 
   R2 = roi**2.0
 
   IF( debug ) THEN
      print *, 'Maxval of anal before:  ', maxval(anal)
      print *, 'Minval of anal before:  ', minval(anal)
      print *, ''
      print *, 'Maxval of observations:  ', maxval(ob)
      print *, 'Minval of observations:  ', minval(ob)
      print *, 'Min/Max xob:    ', minval(x), maxval(x)
      print *, 'Min/Max yob:    ', minval(y), maxval(y)
      print *, 'Min/Max xgrid:  ', minval(xg), maxval(xg)
      print *, 'Min/Max ygrid:  ', minval(yg), maxval(yg)
      print *, 'Radius of influence:  ', roi
      print *, ''
   ENDIF
 
   DO j = 1,ny
    DO i = 1,nx
     w_sum = 0.0
     top   = 0.0  
     anal(i,j) = 0.0
     DO n = 1,nobs
       dis = sqrt( (xg(i) - x(n))**2 + (yg(j)-y(n))**2 )
       IF (dis .le. roi) THEN
         rk2 = dis**2.0
         wk = (R2-rk2) / (R2+rk2)
         top = top + wk*ob(n)
         w_sum = w_sum + wk
       ENDIF
 
 !     IF (dis .le. 4.*roi) THEN
 !       wk = exp( -((dis/(hsp*roi))**2) )
 !       top = top + wk*ob(n)
 !       w_sum = w_sum + wk
 !     ENDIF
 
     ENDDO
 
     IF (w_sum .ge. 0.01) THEN
      anal(i,j) = anal(i,j) + top/w_sum
     ENDIF
 
    ENDDO
   ENDDO
 
   IF( debug ) THEN
      print *, 'Maxval of anal after:  ', maxval(anal)
      print *, 'Minval of anal after:  ', minval(anal)
   ENDIF
 end subroutine cressman
