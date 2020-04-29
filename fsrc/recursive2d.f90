!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE RECURSIVE_FILTER            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  recursive filter in xy two direction. according to paper
!  published by Lorenc and Purser.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  Jidong Gao, CAPS, July, 2000
!
!-----------------------------------------------------------------------
!

SUBROUTINE recurfilt_2d(pgrd_in,pgrd,ipass_filt,radius,nx,ny)
!
!
  INTEGER              :: nx, ny
  INTEGER, INTENT(IN)  :: ipass_filt, radius
  REAL,    INTENT(IN)  :: pgrd_in(nx,ny)
  REAL,    INTENT(OUT) :: pgrd(nx,ny)
  
  REAL, DIMENSION (:), allocatable :: temx
  REAL, DIMENSION (:), allocatable :: temy
!
!
  allocate (temx(nx))
  allocate (temy(ny))
!
!
  ee    = REAL(ipass_filt) / REAL(radius * radius)
  alpha = 1 + ee - SQRT(ee*(ee+2.))
  
  print *, 'FORTRAN:  ', maxval(pgrd_in), minval(pgrd_in), ipass_filt, radius
!
  pgrd(:,:) = pgrd_in(:,:)
!
  DO n = 1, ipass_filt
!
    DO j = 1, ny
!
      DO i = 1, nx
        temx (i) = pgrd(i,j)
      END DO
!
      CALL recurfilt_1d( temx, nx, alpha, n )
!
      DO i = 1, nx
        pgrd(i,j) = temx(i)
      END DO
!
    END DO
!
!
    DO i = 1, nx
!
      DO j = 1, ny
        temy (j) = pgrd(i,j)
      END DO
!
      CALL recurfilt_1d( temy, ny, alpha, n)
!
      DO j = 1, ny
        pgrd(i,j) = temy (j)
      END DO
!
    END DO
!
  END DO
!
!
  deallocate (temx)
  deallocate (temy)
!
  RETURN
END SUBROUTINE  recurfilt_2d
!
!
SUBROUTINE  recurfilt_1d(vara, nx, alpha, ipass)
!
!
  INTEGER :: nx
  REAL :: vara(nx)
  REAL :: alpha
  INTEGER :: ipass
  REAL,   DIMENSION (:), allocatable :: varb
  REAL,   DIMENSION (:), allocatable :: varc
!
  allocate (varb(nx))
  allocate (varc(nx))
!
  DO n = 1, nx
    varb(n) = 0.
    varc(n) = 0.
  END DO
!
!  varb(1) = (1-alpha) * vara(1)
  if(ipass==1) varb(1) = (1-alpha) * vara(1)
  if(ipass==2) varb(1) = (1-alpha)/(1-alpha*alpha) * vara(1)
  if(ipass==3) then
     temp = (1-alpha)/((1-alpha*alpha)*(1-alpha*alpha))
     temp2 =alpha*alpha*alpha
     varb(1) = temp * (vara(1)-temp2*vara(2))
  ENDIF
  if(ipass>=4) then
     temp2 =alpha*alpha*alpha
     temp = (1-alpha)/(1-3*alpha*alpha+3*temp2*alpha-temp2*temp2)
     varb(1) = temp * (vara(1)-3*temp2*vara(2)+                         &
           temp2*alpha*alpha*vara(2)+temp2*alpha*vara(3))
  ENDIF
!
!
  DO i = 2, nx, 1
    varb(i) = alpha*varb(i-1) + (1.-alpha)*vara(i)
  END DO
!
!
!  varc(nx) =  (1./(1.+alpha)) * varb(nx)
  if(ipass==0) varc(nx) = (1-alpha) * varb(nx)
  if(ipass==1) varc(nx) = (1-alpha)/(1-alpha*alpha) * varb(nx)
  if(ipass==2) then
     temp = (1-alpha)/((1-alpha*alpha)*(1-alpha*alpha))
     temp2 =alpha*alpha*alpha
     varc(nx) = temp * (varb(nx)-temp2*varb(nx-1))
  ENDIF
  if(ipass>=3) then
     temp2 =alpha*alpha*alpha
     temp = (1-alpha)/(1-3*alpha*alpha+3*temp2*alpha-temp2*temp2)
     varc(nx) = temp * (varb(nx)-3*temp2*varb(nx-1)+                    &
           temp2*alpha*alpha*varb(nx-1)+temp2*alpha*varb(nx-2))
  ENDIF
!
  DO i = nx-1, 1, -1
    varc(i) = alpha*varc(i+1) + (1.-alpha)*varb(i)
  END DO
!
  DO i = 1, nx
    vara (i) = varc (i)
  END DO
!
  deallocate (varb)
  deallocate (varc)

  RETURN
END SUBROUTINE  recurfilt_1d

