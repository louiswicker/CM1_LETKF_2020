MODULE common_letkf
!=======================================================================
!
! [PURPOSE:] Local Ensemble Transform Kalman Filtering (LETKF)
!            Model Independent Core Module
!
! [REFERENCES:]
!  [1] Ott et al., 2004: A local ensemble Kalman filter for atmospheric
!    data assimilation. Tellus, 56A, 415-428.
!  [2] Hunt et al., 2007: Efficient Data Assimilation for Spatiotemporal
!    Chaos: A Local Ensemble Transform Kalman Filter. Physica D, 230,
!    112-126.
!
! [HISTORY:]
!  01/21/2009 Takemasa Miyoshi  Created at U. of Maryland, College Park
!
!=======================================================================
  USE common_mtx
  
  IMPLICIT NONE

  PUBLIC

!-----------------------------------------------------------------------
! Variable size definitions
!-----------------------------------------------------------------------

  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_dble=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  
!-----------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------

  REAL*8,PARAMETER :: pi=3.1415926535
  REAL*8,PARAMETER :: gg=9.81
  REAL*8,PARAMETER :: rd=287.0
  REAL*8,PARAMETER :: cp=1005.7
  REAL*8,PARAMETER :: re=6371.3e3
  REAL*8,PARAMETER :: r_omega=7.292e-5
  REAL*8,PARAMETER :: t0c=273.15
  REAL*8,PARAMETER :: undef=-9.99e33

CONTAINS
!=======================================================================
!  Main Subroutine of LETKF Core
!   INPUT
!     nobs             : array size, but only first nobsl elements are used
!     nobsl            : total number of observation assimilated at the point
!     hdxb(nobs,nbv)   : obs operator times fcst ens perturbations
!     rdiag(nobs)      : observation error variance
!     rloc(nobs)       : localization weighting function
!     dep(nobs)        : observation departure (yo-Hxb)
!     parm_infl        : covariance inflation parameter
!   OUTPUT
!     trans(nbv,nbv)   : transformation matrix
!     wmean(nbv)       : weights for mean analysis 
!     inflout          : adaptive inflation value at point
!=======================================================================
SUBROUTINE letkf_core(nbv,nobs,nobsl,var,hdxb,rdiag,rloc,dep,obIndex,infl_in,trans,wmean,infl_out,hf,rcp_return)

  IMPLICIT NONE

! Passed variables

  INTEGER,INTENT(IN)   :: nobs
  INTEGER,INTENT(IN)   :: nobsl
  INTEGER,INTENT(IN)   :: nbv
  REAL*8,INTENT(IN)    :: infl_in(1)      ! adaptive inflation parameter
  REAL*8,INTENT(OUT)   :: trans(nbv,nbv)  ! weight matrix (nobs,nobs)
  REAL*8,INTENT(OUT)   :: wmean(nbv)      ! weights for mean analysis (nbv x 1)
  REAL*8,INTENT(OUT)   :: infl_out(1)     ! adaptive inflation parameter

  REAL*4,INTENT(IN)    :: var(nbv)
  REAL*8,INTENT(IN)    :: hdxb(nobs,nbv)  ! y - Hx
  REAL*8,INTENT(IN)    :: rdiag(nobs)     ! observation error
  REAL*8,INTENT(IN)    :: rloc(nobs)      ! spatial localization weight
  REAL*8,INTENT(IN)    :: dep(nobs)       ! y - Hx_mean
  INTEGER,INTENT(IN)   :: obIndex(nobs)
  INTEGER,INTENT(IN)   :: hf
  REAL*8,INTENT(INOUT) :: rcp_return(1)

  
! Local variables

  REAL*8 :: YbL_rinv(nobsl,nbv)
  REAL*8 :: eivec(nbv,nbv)
  REAL*8 :: eival(nbv)
  REAL*8 :: pa(nbv,nbv)
  REAL*8 :: work1(nbv,nbv)
  REAL*8 :: work2(nbv,nobsl)
! REAL*8 :: work3(nbv)                                               ! LJW replaced with wmean, to return mean analysis weights
  REAL*8 :: Ybl(nobsl,nbv), rdiagL(nobsl), rlocL(nobsl), depL(nobsl) ! LJW changes
  REAL*8 :: rho, obs_error
  REAL*8 :: HxfL(nbv), Hxfmean, sdHxf
  REAL*8 :: parm(4),sigma_o,gain
  REAL*8 :: mean, stddev, rcp
  REAL*8, PARAMETER :: rcp_weight = 0.5
  REAL*8, PARAMETER :: sigma_b = 0.10d0 ! error stdev of parm_infl
  INTEGER :: i,j,k,m
  
  REAL(KIND=4)  :: F_RCP
  EXTERNAL F_RCP

  LOGICAL, PARAMETER :: wgt_regularization  = .false.

  LOGICAL, PARAMETER :: DEBUG  = .false.
  LOGICAL, PARAMETER :: DEBUG2 = .false.
  LOGICAL, PARAMETER :: DEBUG3 = .false.  ! Debug flag for adaptive inflation stats

  IF ( DEBUG ) THEN
    print *, ""
    print *, "Hello from letkf core"
    print *, "NE: ",nbv, "NOBS:  ", nobs, "NVALID OBS:  ", nobsl
    print *, 'YB(2,3) = ', hdxb(2,3),  'YB(3,1) = ', hdxb(3,1)    
    print *, 'INFLATE IN: ', infl_in(1)
  ENDIF
  
  IF(nobsl == 0) THEN
    trans = 0.0d0
    DO i=1,nbv
      trans(i,i) = SQRT(infl_in(1))
      wmean(i)   = SQRT(infl_in(1))
    END DO
    infl_out(1)  = infl_in(1)
    RETURN
  ELSE
  
  DO m = 1,nobsl
    depL(m)   = dep(obIndex(m))
    YbL(m,:)  = hdxb(obIndex(m),:) - SUM(hdxb(obIndex(m),:))/float(nbv)
    rdiagL(m) = rdiag(obindex(m))
    rlocL(m)  = rloc(obindex(m))
  ENDDO

  IF( wgt_regularization ) THEN
    DO m = 1,nobsl
      HxfL(:) = YbL(m,:)
      call mean_stddev(HxfL,nbv,Hxfmean,sdHxf)
      rlocL(m)  = (rlocL(m)*rdiagL(m)**2 / (rdiagL(m)**2 + sdHxf**2)) &
                / (1.0 - (rlocL(m)*sdHxf**2/(rdiagL(m)**2 + sdHxf**2)))
    ENDDO
  ENDIF

  IF( hf > 0 ) THEN
    rcp_return(1) = 0.0
    DO m = 1,nobsl
      HxfL(:)    = hdxb(obIndex(m),:)
      rcp        = F_RCP(var, HxfL, nbv, hf)
      rlocL(m)   = rcp_weight * rcp + (1.0 - rcp_weight) * rlocL(m)
      rcp_return(1) = rcp_return(1) + rcp
    ENDDO
    rcp_return(1) = rcp_return(1) / float(nobsl)
  ENDIF
  
!-----------------------------------------------------------------------
!  hdxb Rinv
!-----------------------------------------------------------------------
  IF( DEBUG2 )  print*, "Step 4"
  DO j=1,nbv
    DO i=1,nobsl
      YbL_rinv(i,j) = YbL(i,j) / rdiagL(i) * rlocL(i)
    END DO
  END DO
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb
!-----------------------------------------------------------------------
  IF( DEBUG2 ) print *, "Step 5A"
  CALL dgemm('t','n',nbv,nbv,nobsl,1.0d0,YbL_rinv,nobsl,YbL,nobsl,0.0d0,work1,nbv)
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb + (m-1) I / rho (covariance inflation)
!-----------------------------------------------------------------------
  IF( DEBUG2 ) print *,   "Step 5B"
  rho = 1.0d0 / infl_in(1)
! rho = 1.0d0 
  DO i=1,nbv
    work1(i,i) = work1(i,i) + REAL(nbv-1,r_size) * rho
  END DO
!-----------------------------------------------------------------------
!  eigenvalues and eigenvectors of [ hdxb^T Rinv hdxb + (m-1) I ]
!-----------------------------------------------------------------------
  IF( DEBUG2 )  print *, "Step 5C"
  CALL mtx_eigen(1,nbv,work1,eival,eivec,i)
!-----------------------------------------------------------------------
!  Pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv
!-----------------------------------------------------------------------
  IF( DEBUG2 ) print *,  "Step 5D"
  DO j=1,nbv
    DO i=1,nbv
      work1(i,j) = eivec(i,j) / eival(j)
    END DO
  END DO
  CALL dgemm('n','t',nbv,nbv,nbv,1.0d0,work1,nbv,eivec,nbv,0.0d0,pa,nbv)
!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T
!-----------------------------------------------------------------------
  IF( DEBUG2 ) print *,  "Step 6"
  CALL dgemm('n','t',nbv,nobsl,nbv,1.0d0,pa,nbv,YbL_rinv,nobsl,0.0d0,work2,nbv)
!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
  DO i=1,nbv
    wmean(i) = work2(i,1) * depL(1)
    DO j=2,nobsl
      wmean(i) = wmean(i) + work2(i,j) * depL(j)
    END DO
  END DO
  IF( DEBUG2 )  print *,"Mean weight computed in common_letkf"
!-----------------------------------------------------------------------
!  T = sqrt[(m-1)Pa]
!-----------------------------------------------------------------------
  IF( DEBUG2 ) print*, "Step 7"
  DO j=1,nbv
    rho = SQRT( REAL(nbv-1,r_size) / eival(j) )
    DO i=1,nbv
      work1(i,j) = eivec(i,j) * rho
    END DO
  END DO
  CALL dgemm('n','t',nbv,nbv,nbv,1.0d0,work1,nbv,eivec,nbv,0.0d0,trans,nbv)
!-----------------------------------------------------------------------
!  T + Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nbv
      trans(i,j) = trans(i,j) + wmean(i)
    END DO
  END DO
  IF( DEBUG2 )  print *,"Weights computed in common_letkf"
!-----------------------------------------------------------------------
!  Inflation estimation
! 
! Did some fixes where the number of observations < 3, was
! generating Nan...from parm(2) == 0
!
!-----------------------------------------------------------------------
  IF ( DEBUG ) THEN  
  ENDIF

  parm(:) = 0.0d0
  DO i=1,nobsl
    parm(1) = parm(1) + depL(i)*depL(i)/rdiagL(i) * rlocL(i)
    parm(3) = parm(3) + rlocL(i)
  END DO

  DO j=1,nbv
    DO i=1,nobsl
      parm(2) = parm(2) + YbL_rinv(i,j) * YbL(i,j)
    END DO
  END DO

  parm(2)   = max( parm(2) / REAL(nbv-1,r_size), 1.0d-8)
! parm(3)   = SUM(rlocL)
  parm(4)   = (parm(1)-parm(3))/parm(2) - infl_in(1)
! sigma_o   = 1.0d0/REAL(nobsl,r_size)/MAXVAL(rlocL(1:nobsl))
  sigma_o   = 2.0d0/parm(3)*((infl_in(1)*parm(2)+parm(3))/parm(2))**2
  gain      = sigma_b**2 / (sigma_o + sigma_b**2)
  infl_out(1) = infl_in(1) + gain * parm(4)

  IF( DEBUG3 ) THEN
    IF( gain * parm(4) > 1.0 ) THEN
      print *, 'Adaptive inflation....'
      print *, ""
      print *, 'N_OBS:  ', nobsl
      print *, 'PARAM(1): ', parm(1), ' PARAM(2): ', parm(2)
      print *, 'PARAM(3): ', parm(3), ' PARAM(4): ', parm(4)
      print *, 'SIGMA_O:  ', sigma_o, ' KGain: ', gain
      print *, 'INPUT INFLATE:  ', infl_in(1), ' INC: ', gain*parm(4)
    ENDIF
  ENDIF

  IF ( DEBUG ) THEN  
    print *, 'TRANS(2,3) = ', trans(2,3),  'TRANS(3,1) = ', trans(3,1)    
    print *, ""
  ENDIF

  RETURN
  END IF
END SUBROUTINE letkf_core

END MODULE common_letkf
