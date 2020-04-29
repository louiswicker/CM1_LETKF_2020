
      subroutine sfc_and_turb(getsfc,nstep,dt,dosfcflx,cloudvar,qbudget,    &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,  &
                   yh,vh,rvh,yf,vf,rvf,                              &
                   rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,zf,mf,rmf,  &
                   pi0s,rth0s,pi0,rho0,prs0,thv0,th0,qv0,            &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv,        &
                   tsk,thflux,qvflux,cd,ch,cq,u1,v1,s1,tlh,          &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,          &
                   divx,rho,rr,rf,prs,                               &
                   t11,t12,t13,t22,t23,t33,                          &
                   u0,ua,v0,va,wa,                                   &
                   ppi,pp3d,ppten,                                   &
                   tha,th3d,thten,thten1,qa,                         &
                   kmh,kmv,khh,khv,tkea,tke3d,                       &
                   nm,defv,defh,dissten,radsw,radswnet,radlwin,      &
                   thpten,qvpten,qcpten,qipten,upten,vpten,          &
                   lu_index,kpbl2d,psfc,u10,v10,s10,hfx,qfx,xland,znt,ust,    &
                   hpbl,wspd,psim,psih,gz1oz0,br,                    &
                   CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                    &
                   MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                 &
                   CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,    &
                   f2d,gsw,glw,chklowq,capg,snowc,dsxy,wstar,delta,fm,fh,  &
                   mznt,smois,taux,tauy,hpbl2d,evap2d,heat2d,rc2d,   &
                   num_soil_layers,slab_zs,slab_dzs,tslb,tmn,        &
                   tml,t0ml,hml,h0ml,huml,hvml,tmoml,                &
                   reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,               &
                   nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                  &
                   pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,                  &
                   vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,                  &
                   uw31,uw32,ue31,ue32,us31,us32,un31,un32,          &
                   kw1,kw2,ke1,ke2,ks1,ks2,kn1,kn2,                  &
                   rtime,ntdiag,tdiag,update_sfc)
      ! end_sfc_and_turb
      use module_sf_sfclay
      use module_sf_sfclayrev
      use module_sf_slab
      use module_sf_oml
      use module_bl_ysu
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'radcst.incl'
      include 'timestat.incl'

!-----------------------------------------------------------------------
! Arrays and variables passed into solve

      logical, intent(in) :: getsfc
      integer, intent(in) :: nstep
      real, intent(inout) :: dt
      logical, intent(in) :: dosfcflx
      logical, intent(in), dimension(maxq) :: cloudvar
      double precision, intent(inout), dimension(nbudget) :: qbudget
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf,ruf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf,rvf
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf,rmf
      real, intent(in), dimension(ib:ie,jb:je) :: pi0s,rth0s
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,prs0,thv0,th0,qv0
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,znt,ust,thflux,qvflux,cd,ch,cq,u1,v1,s1,xland,psfc,tlh
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,divx,rho,rr,rf,prs
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,ppten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tha,th3d,thten,thten1
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tke3d
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh,dissten
      real, intent(inout), dimension(ni,nj) :: radsw,radswnet,radlwin
      real, intent(inout), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: thpten,qvpten,qcpten,qipten,upten,vpten
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: lu_index,kpbl2d
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: u10,v10,s10,hfx,qfx, &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,       &
                                      CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      f2d,gsw,glw,chklowq,capg,snowc,dsxy,wstar,delta,fm,fh
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: mznt,smois,taux,tauy,hpbl2d,evap2d,heat2d,rc2d
      integer, intent(in) :: num_soil_layers
      real, intent(in), dimension(num_soil_layers) :: slab_zs,slab_dzs
      real, intent(inout), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_p
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, intent(inout), dimension(jmp,kmp) :: vw1,vw2,ve1,ve2
      real, intent(inout), dimension(imp,kmp) :: vs1,vs2,vn1,vn2
      real, intent(inout), dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(jmp,kmt,4)     :: kw1,kw2,ke1,ke2
      real, intent(inout), dimension(imp,kmt,4)     :: ks1,ks2,kn1,kn2
      real, intent(in) :: rtime
      integer, intent(in) :: ntdiag
      real, intent(inout) , dimension(ibd:ied,jbd:jed,kbd:ked,ntdiag) :: tdiag
      logical, intent(in) :: update_sfc

!-----------------------------------------------------------------------

      integer :: i,j,k,k2
      integer :: isfflx,ifsnow
      real :: ep1,ep2,rovg,dtmin,dz1
      real :: SVP1,SVP2,SVP3,SVPT0,p1000mb,eomeg,stbolt,tem,tem1,tem2,tem3

      integer :: NTSFLG
      real :: SFENTH
      logical :: flag_qi
      integer :: p_qi,p_first_scalar
      logical :: disheat
      real :: alpha,beta,var_ric,coef_ric_l,coef_ric_s
      real :: qx

!-----------------------------------------------------------------------

      IF( iturb.ge.1 .or. ipbl.ge.1 .or. idiss.eq.1 .or. output_dissten.eq.1 )THEN

        ! cm1r17:  dissten is defined on w (full) levels:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dissten(i,j,k)=0.0
        enddo
        enddo
        enddo

      ENDIF

      do j=1,nj
      do i=1,ni
        psfc(i,j) = cgs1*prs(i,j,1)+cgs2*prs(i,j,2)+cgs3*prs(i,j,3)
      enddo
      enddo

      IF( (sfcmodel.eq.2) .or. (sfcmodel.eq.3) .or. (sfcmodel.eq.4) .or. (ipbl.eq.1) .or. (ipbl.eq.2) .or.  (oceanmodel.eq.2) )THEN

        ! variables for wrf physics:

        if( ipbl.ge.1 )then
          k2 = nk
        else
          k2 = 1
        endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do j=1,nj
          do k=1,k2
          do i=1,ni
            dum1(i,j,k)=0.5*(ua(i,j,k)+ua(i+1,j,k))
            dum2(i,j,k)=0.5*(va(i,j,k)+va(i,j+1,k))
            dum3(i,j,k)=th0(i,j,k)+tha(i,j,k)
            dum7(i,j,k)=pi0(i,j,k)+ppi(i,j,k)
            dum4(i,j,k)=dum3(i,j,k)*dum7(i,j,k)
          enddo
          enddo
          do k=1,max(2,k2)
          do i=1,ni
            dum5(i,j,k) = dz*rmh(i,j,k)
          enddo
          enddo
          do k=2,max(2,k2)
          do i=1,ni
            dum6(i,j,k) = c1(i,j,k)*prs(i,j,k-1)+c2(i,j,k)*prs(i,j,k)
          enddo
          enddo
          ! surface:
          do i=1,ni
            dum6(i,j,1) = psfc(i,j)
          enddo
          ! top of model:
          if( k2.gt.1 )then
          do i=1,ni
            dum6(i,j,nk+1)= cgt1*prs(i,j,nk)+cgt2*prs(i,j,nk-1)+cgt3*prs(i,j,nk-2)
          enddo
          endif
        enddo

        ! dum1 = u at scalars
        ! dum2 = v at scalars
        ! dum3 = th
        ! dum7 = pi
        ! dum4 = t
        ! dum5 = dz8w
        ! dum6 = p3di

        isfflx = 1
        SVP1=0.6112
        SVP2=17.67
        SVP3=29.65
        SVPT0=273.15
        p1000mb      = 100000.
        EOMEG=7.2921E-5
        STBOLT=5.67051E-8
        ep1 = rv/rd - 1.0
        ep2 = rd/rv
        rovg = rd/g

        IF(radopt.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            gsw(i,j)=radswnet(i,j)
            glw(i,j)=radlwin(i,j)
          enddo
          enddo
        ELSE
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            gsw(i,j)=0.0
            glw(i,j)=0.0
          enddo
          enddo
        ENDIF

      ENDIF

!-----------------------------------------------------------------------

    dosfc:  IF( getsfc )THEN

      bbc3:  IF( bbc.eq.3 )THEN

        !-------------------------------
        ! u1 is u at lowest model level
        ! v1 is v at lowest model level
        ! s1 is horizontal wind speed at lowest model level
        ! (all defined at the scalar point of the staggered grid)
        ! for pertflx=1, account for domain (i.e., surface) motion 
        !                in calculation of wind speed
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
          IF( pertflx.eq.1 .or. imove.eq.1 )THEN
            do i=1,ni
              u1(i,j) = 0.5*( ua(i,j,1)+umove + ua(i+1,j,1)+umove )
              v1(i,j) = 0.5*( va(i,j,1)+vmove + va(i,j+1,1)+vmove )
              s1(i,j) = sqrt(u1(i,j)**2+v1(i,j)**2)
            enddo
          ELSE 
            do i=1,ni 
              u1(i,j) = 0.5*( ua(i,j,1) + ua(i+1,j,1) )
              v1(i,j) = 0.5*( va(i,j,1) + va(i,j+1,1) )
              s1(i,j) = sqrt(u1(i,j)**2+v1(i,j)**2)
            enddo
          ENDIF
          IF(imoist.eq.1)THEN
            do k=1,nk
            do i=1,ni
              divx(i,j,k) = qa(i,j,k,nqv)
            enddo
            enddo
          ELSE
            do k=1,nk
            do i=1,ni
              divx(i,j,k) = 0.0
            enddo
            enddo
          ENDIF
          IF( terrain_flag )THEN
            do i=1,ni
              ppten(i,j,1) = zh(i,j,1)-zs(i,j)
            enddo
          ELSE
            do i=1,ni
              ppten(i,j,1) = zh(i,j,1)
            enddo
          ENDIF

        enddo

        !-------------------------------
        ! NOTE:
        ! divx stores qv
        ! ppten stores height of first model level above surface

        IF( ipbl.eq.0 )THEN
          call gethpbl(zh,th0,tha,divx,hpbl)
          IF( terrain_flag )THEN
            do j=1,nj
            do i=1,ni
              hpbl(i,j) = hpbl(i,j)-zs(i,j)
            enddo
            enddo
          ENDIF
        ENDIF

        IF( sfcmodel.eq.1 )THEN

          call getcecd(u0,v0,u1,v1,s1,ua,va,ppten(ib,jb,1),u10,v10,s10,xland,znt,ust,cd,ch,cq)
          if(isfcflx.eq.1)then
            call sfcflux(dt,ruh,xf,rvh,pi0s,ch,cq,pi0,thv0,th0,u0,v0,tsk,thflux,qvflux,mavail, &
                         rho,rf,u1,v1,s1,ua,va,ppi,tha,qa(ibm,jbm,kbm,nqv), &
                         qbudget(8),psfc,u10,v10,s10,qsfc,znt,rtime)
          endif
          ! get sfc diagnostics needed by pbl scheme:
          call sfcdiags(tsk,thflux,qvflux,cd,ch,cq,u1,v1,s1,             &
                        xland,psfc,qsfc,u10,v10,hfx,qfx,cda,znt,gz1oz0,  &
                        psim,psih,br,zol,mol,hpbl,dsxy,th2,t2,q2,fm,fh,  &
                        zs,ppten(ib,jb,1),pi0s,pi0,th0,ppi,tha,rho,rf,divx,ua,va)
          if( ipbl.ge.1 )then
            do j=1,nj
            do i=1,ni
              wspd(i,j) = max( s1(i,j) , 1.0e-6 )
              CPMM(i,j)=CP*(1.0+0.8*divx(i,j,1))                                   
              hfx(i,j) = thflux(i,j)*CPMM(i,j)*rf(i,j,1)
              qvflux(i,j) = qfx(i,j)*rf(i,j,1)
            enddo
            enddo
          endif

        ENDIF


        IF( (sfcmodel.eq.2) .or. (sfcmodel.eq.3) .or. (sfcmodel.eq.4) )THEN
          ! surface layer:
          ! (needed by sfcmodel=2,3,4 and ipbl=1,2)
        if( sfcmodel.eq.2 )then
          call SFCLAY(dum1,dum2,dum4,qa(ib,jb,kb,nqv),prs,dum5,      &
                       CP,G,ROVCP,RD,XLV,lv1,lv2,PSFC,CHS,CHS2,CQS2,CPMM, &
                       ZNT,UST,hpbl,MAVAIL,ZOL,MOL,REGIME,PSIM,PSIH, &
                       FM,FH,                                        &
                       XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,QGH,QSFC,RMOL, &
                       U10,V10,TH2,T2,Q2,rf(ib,jb,1),                &
                       GZ1OZ0,WSPD,BR,ISFFLX,dsxy,                   &
                       SVP1,SVP2,SVP3,SVPT0,EP1,EP2,                 &
                       KARMAN,EOMEG,STBOLT,                          &
                       P1000mb,                                      &
                       1  ,ni+1 , 1  ,nj+1 , 1  ,nk+1 ,              &
                       ib ,ie , jb ,je , kb ,ke ,                    &
                       1  ,ni , 1  ,nj , 1  ,nk ,                    &
                       ustm,ck,cka,cd,cda,isftcflx,iz0tlnd           )
        elseif( sfcmodel.eq.3 )then
          call SFCLAYREV(dum1,dum2,dum4,qa(ib,jb,kb,nqv),prs,dum5,   &
                       CP,G,ROVCP,RD,XLV,lv1,lv2,PSFC,CHS,CHS2,CQS2,CPMM, &
                       ZNT,UST,hpbl,MAVAIL,ZOL,MOL,REGIME,PSIM,PSIH, &
                       FM,FH,                                        &
                       XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,QGH,QSFC,RMOL, &
                       U10,V10,TH2,T2,Q2,rf(ib,jb,1),                &
                       GZ1OZ0,WSPD,BR,ISFFLX,dsxy,                   &
                       SVP1,SVP2,SVP3,SVPT0,EP1,EP2,                 &
                       KARMAN,EOMEG,STBOLT,                          &
                       P1000mb,                                      &
                       1  ,ni+1 , 1  ,nj+1 , 1  ,nk+1 ,              &
                       ib ,ie , jb ,je , kb ,ke ,                    &
                       1  ,ni , 1  ,nj , 1  ,nk ,                    &
                       ustm,ck,cka,cd,cda,isftcflx,iz0tlnd           )
        endif

          ifsnow = 0
          dtmin = dt/60.0

        IF( update_sfc )THEN
          ! slab scheme (MM5/WRF):
          call SLAB(dum4,qa(ib,jb,kb,nqv),prs,FLHC,FLQC,            &
                       PSFC,XLAND,TMN,HFX,QFX,LH,TSK,QSFC,CHKLOWQ,  &
                       GSW,GLW,CAPG,THC,SNOWC,EMISS,MAVAIL,         &
                       DT,ROVCP,XLV,lv1,lv2,DTMIN,IFSNOW,           &
                       SVP1,SVP2,SVP3,SVPT0,EP2,                    &
                       KARMAN,EOMEG,STBOLT,                         &
                       TSLB,slab_ZS,slab_DZS,num_soil_layers, .true. ,       &
                       P1000mb,                                     &
                         1, ni+1,   1, nj+1,   1, nk+1,             &
                        ib, ie,  jb, je,  kb, ke,                   &
                         1, ni,   1, nj,   1, nk                    )
        ELSE
          ! dont update tsk, but diagnose qsfc:
          do j=1,nj
          do i=1,ni
            qx = divx(i,j,1)
            if ( FLQC(i,j) .ne. 0.) then
               QSFC(i,j)=QX+QFX(i,j)/FLQC(i,j)
            else
               QSFC(i,j) = QX
            end if
            CHKLOWQ(i,j)=MAVAIL(i,j)
          enddo
          enddo
        ENDIF

          ! put WRF parameters into CM1 arrays:
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            ch(i,j) = chs2(i,j)
            cq(i,j) = cqs2(i,j)
            s10(i,j) = sqrt( u10(i,j)**2 + v10(i,j)**2 )
          enddo
          enddo
          IF( dosfcflx .or. output_sfcflx.eq.1 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni
              thflux(i,j) = hfx(i,j)/(CPMM(i,j)*rf(i,j,1))
              qvflux(i,j) = qfx(i,j)/rf(i,j,1)
            enddo
            enddo
          ENDIF

        ENDIF

      ENDIF  bbc3
      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()

    ELSE

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  ... skipping sfc stuff ... '
      if(dowr) write(outfile,*)

    ENDIF  dosfc

!---------------------------------------------
! bc/comms (very important):

  IF( bbc.eq.3 )THEN
    !-------------!
    call bc2d(ust)
#ifdef MPI
    call comm_1s2d_start(ust,uw31(1,1,1),uw32(1,1,1),ue31(1,1,1),ue32(1,1,1),  &
                             us31(1,1,1),us32(1,1,1),un31(1,1,1),un32(1,1,1),reqs_s)
#endif
    call bc2d(u1)
#ifdef MPI
    call comm_1s2d_start(u1 ,uw31(1,1,2),uw32(1,1,2),ue31(1,1,2),ue32(1,1,2),  &
                             us31(1,1,2),us32(1,1,2),un31(1,1,2),un32(1,1,2),reqs_u)
#endif
    call bc2d(v1)
#ifdef MPI
    call comm_1s2d_start(v1 ,uw31(1,1,3),uw32(1,1,3),ue31(1,1,3),ue32(1,1,3),  &
                             us31(1,1,3),us32(1,1,3),un31(1,1,3),un32(1,1,3),reqs_v)
#endif
    call bc2d(s1)
#ifdef MPI
    call comm_1s2d_start(s1 ,uw31(1,1,4),uw32(1,1,4),ue31(1,1,4),ue32(1,1,4),  &
                             us31(1,1,4),us32(1,1,4),un31(1,1,4),un32(1,1,4),reqs_w)
#endif
    call bc2d(znt)
#ifdef MPI
    call comm_1s2d_start(znt,uw31(1,1,5),uw32(1,1,5),ue31(1,1,5),ue32(1,1,5),  &
                             us31(1,1,5),us32(1,1,5),un31(1,1,5),un32(1,1,5),reqs_p)
#endif
    !-------------!

#ifdef MPI
    !-------------!
    call comm_1s2d_end(ust,uw31(1,1,1),uw32(1,1,1),ue31(1,1,1),ue32(1,1,1),  &
                           us31(1,1,1),us32(1,1,1),un31(1,1,1),un32(1,1,1),reqs_s)
    call bcs2_2d(ust)

    call comm_1s2d_end(u1 ,uw31(1,1,2),uw32(1,1,2),ue31(1,1,2),ue32(1,1,2),  &
                           us31(1,1,2),us32(1,1,2),un31(1,1,2),un32(1,1,2),reqs_u)
    call bcs2_2d(u1 )

    call comm_1s2d_end(v1 ,uw31(1,1,3),uw32(1,1,3),ue31(1,1,3),ue32(1,1,3),  &
                           us31(1,1,3),us32(1,1,3),un31(1,1,3),un32(1,1,3),reqs_v)
    call bcs2_2d(v1 )

    call comm_1s2d_end(s1 ,uw31(1,1,4),uw32(1,1,4),ue31(1,1,4),ue32(1,1,4),  &
                           us31(1,1,4),us32(1,1,4),un31(1,1,4),un32(1,1,4),reqs_w)
    call bcs2_2d(s1 )

    call comm_1s2d_end(znt,uw31(1,1,5),uw32(1,1,5),ue31(1,1,5),ue32(1,1,5),  &
                           us31(1,1,5),us32(1,1,5),un31(1,1,5),un32(1,1,5),reqs_p)
    call bcs2_2d(znt)
    !-------------!

    !-------------!
    call comm_2d_corner(ust)
    call comm_2d_corner(u1)
    call comm_2d_corner(v1)
    call comm_2d_corner(s1)
    call comm_2d_corner(znt)
    !-------------!
#endif
  ENDIF

!-------------------------------------------------------------------
! simple ocean mixed layer model based Pollard, Rhines and Thompson (1973)
!   (from WRF)

    IF(oceanmodel.eq.2)THEN
    IF( update_sfc )THEN
      if( getsfc )then

        CALL oceanml(tml,t0ml,hml,h0ml,huml,hvml,ust,dum1,dum2, &
                     tmoml,f2d,g,oml_gamma,                     &
                     xland,hfx,lh,tsk,gsw,glw,emiss,            &
                     dt,STBOLT,                                 &
                       1, ni+1,   1, nj+1,   1, nk+1,           &
                      ib, ie,  jb, je,  kb, ke,                 &
                       1, ni,   1, nj,   1, nk                  )

        if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()

      endif
    ENDIF
    ENDIF

!-------------------------------------------------------------------
!  PBL scheme:

      IF(ipbl.ge.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        DO k=1,nk
          ! store qi in dum8:
          if( nqi.ge.1 )then
            do j=1,nj
            do i=1,ni
              dum8(i,j,k) = qa(i,j,k,nqi)
            enddo
            enddo
          else
            do j=1,nj
            do i=1,ni
              dum8(i,j,k) = 0.0
            enddo
            enddo
          endif
          IF(output_km.eq.1.or.output_kh.eq.1)THEN
            do j=1,nj
            do i=1,ni
              thten1(i,j,k)=0.0
              thten(i,j,k)=0.0
            enddo
            enddo
          ENDIF
        ENDDO

        if( iice.eq.1 )then
          flag_qi = .true.
        else
          flag_qi = .false.
        endif

      if(ipbl.eq.1)then
        ! PBL:
        call ysu(u3d=dum1,v3d=dum2,th3d=dum3,t3d=dum4,qv3d=qa(ib,jb,kb,nqv),        &
                  qc3d=qa(ib,jb,kb,nqc),qi3d=dum8,p3d=prs,p3di=dum6,pi3d=dum7,      &
                  rublten=upten,rvblten=vpten,rthblten=thpten,                      &
                  rqvblten=qvpten,rqcblten=qcpten,rqiblten=qipten,flag_qi=flag_qi,  &
                  cp=cp,g=g,rovcp=rovcp,rd=rd,rovg=rovg,ep1=ep1,ep2=ep2,            &
                  karman=karman,xlv=xlv,lv1=lv1,lv2=lv2,rv=rv,                      &
                  dz8w=dum5 ,psfc=psfc,                                             &
                  znt=znt,ust=ust,hpbl=hpbl,psim=fm,psih=fh,                        &
                  xland=xland,hfx=hfx,qfx=qfx,wspd=wspd,br=br,                      &
                  dt=dt,kpbl2d=kpbl2d,                                              &
                  exch_h=thten1,exch_m=thten,                                       &
                  wstar=wstar,delta=delta,                                          &
                  u10=u10,v10=v10,                                                  &
                  ids=1  ,ide=ni+1 , jds= 1 ,jde=nj+1 , kds=1  ,kde=nk+1 ,          &
                  ims=ib ,ime=ie   , jms=jb ,jme=je   , kms=kb ,kme=ke ,            &
                  its=1  ,ite=ni   , jts=1  ,jte=nj   , kts=1  ,kte=nk ,            &
                  regime=regime)
        IF(output_km.eq.1.or.output_kh.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk+1
          do j=1,nj
          do i=1,ni
            khv(i,j,k) = thten1(i,j,k)
            kmv(i,j,k) = thten(i,j,k)
          enddo
          enddo
          enddo
        ENDIF
      endif
        if(timestats.ge.1) time_pbl=time_pbl+mytime()

        call bcs(upten)
#ifdef MPI
        call comm_1s_start(upten,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
#endif
        call bcs(vpten)
#ifdef MPI
        call comm_1s_start(vpten,vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,reqs_p)
#endif
        ! Dissipative heating from ysu scheme:
        IF( idiss.eq.1 .or. output_dissten.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k,tem)
          do j=1,nj
          do i=1,ni
            ! assume t13,t23 are zero at top of domain:
            t13(i,j,nk+1) = 0.0
            t23(i,j,nk+1) = 0.0
            do k=nk,1,-1
              tem = rho(i,j,k)*dz*rmh(i,j,k)/rf(i,j,k)
              t13(i,j,k) = t13(i,j,k+1)-upten(i,j,k)*tem
              t23(i,j,k) = t23(i,j,k+1)-vpten(i,j,k)*tem
            enddo
            do k=2,nk
              tem = rdz*mf(i,j,k)
              dum3(i,j,k)=( dum1(i,j,k)-dum1(i,j,k-1) )*tem
              dum4(i,j,k)=( dum2(i,j,k)-dum2(i,j,k-1) )*tem
            enddo
            dissten(i,j,1) = dissten(i,j,1) + (ust(i,j)**3)/(karman*znt(i,j))
            do k=2,nk
              ! NOTE:  dissten is defined at w points:
              dissten(i,j,k) = dissten(i,j,k)         &
                            +( t13(i,j,k)*dum3(i,j,k) &
                              +t23(i,j,k)*dum4(i,j,k) )
            enddo
          enddo
          enddo
        ENDIF
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
#ifdef MPI
        call comm_1s_end(upten,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
        call comm_1s_end(vpten,vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,reqs_p)
#endif

      ENDIF

!---------------------------------------------

      IF( iturb.ge.1 .or. output_nm.eq.1 )THEN
        ! squared Brunt-Vaisala frequency:
        call calcnm(c1,c2,mf,pi0,thv0,th0,cloudvar,nm,dum1,dum2,dum3,dum4,   &
                    prs,ppi,tha,qa)
      ENDIF

      IF( iturb.ge.1 .or. output_def.eq.1 )THEN
        ! deformation:
        call calcdef(    rds,sigma,rdsf,sigmaf,zs,gz,rgz,gzu,rgzu,gzv,rgzv,                &
                     xh,rxh,arh1,arh2,uh,xf,rxf,arf1,arf2,uf,vh,vf,mh,c1,c2,mf,defv,defh,  &
                     dum1,dum2,ua,va,wa,t11,t12,t13,t22,t23,t33,gx,gy,rho,rr,rf)

      ENDIF
      if(timestats.ge.1) time_turb=time_turb+mytime()

!--------------------------------------
!  next section is for iturb >= 1 only:

    IF(iturb.ge.1)THEN

      IF(iturb.eq.1)THEN
        call tkekm(nstep,dt,ruh,rvh,rmh,zf,mf,rmf,znt,ust,rf,         &
                   nm,defv,defh,dum1,dum2,dum3,dum4,dum5,             &
                   kmh,kmv,khh,khv,tkea,ua,va,dissten,                &
                   nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                   &
                   kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),       &
                   ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),       &
                   kw1(1,1,2),kw2(1,1,2),ke1(1,1,2),ke2(1,1,2),       &
                   ks1(1,1,2),ks2(1,1,2),kn1(1,1,2),kn2(1,1,2),       &
                   kw1(1,1,3),kw2(1,1,3),ke1(1,1,3),ke2(1,1,3),       &
                   ks1(1,1,3),ks2(1,1,3),kn1(1,1,3),kn2(1,1,3),       &
                   kw1(1,1,4),kw2(1,1,4),ke1(1,1,4),ke2(1,1,4),       &
                   ks1(1,1,4),ks2(1,1,4),kn1(1,1,4),kn2(1,1,4))
      ELSEIF(iturb.eq.2)THEN
        call turbsmag(nstep,dt,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,rf, &
                      nm,defv,defh,dum4,dum5,thten1,zf,znt,ust,        &
                      kmh,kmv,khh,khv,ua,va,dissten,                   &
                      nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                      kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),     &
                      ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),     &
                      kw1(1,1,2),kw2(1,1,2),ke1(1,1,2),ke2(1,1,2),     &
                      ks1(1,1,2),ks2(1,1,2),kn1(1,1,2),kn2(1,1,2))
      ELSEIF(iturb.eq.3)THEN
        call turbparam(nstep,zf,dt,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,rf, &
                      nm,defv,defh,dum4,kmh,kmv,khh,khv,ua,va,dissten,zs,znt,ust,xland,psfc,tlh, &
                      nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                      kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),     &
                      ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),     &
                      kw1(1,1,2),kw2(1,1,2),ke1(1,1,2),ke2(1,1,2),     &
                      ks1(1,1,2),ks2(1,1,2),kn1(1,1,2),kn2(1,1,2))
      ENDIF


        !  now, get turbulent stresses:
        call     gettau(xf,rxf,arf1,arf2,ust,u1,v1,s1,rf,  &
                        kmh,kmv,t11,t12,t13,t22,t23,t33,ua)


      !  last step:  Surface dissipation:
      IF(bbc.eq.2)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,tem1)
        do j=1,nj
        do i=1,ni
          tem1 = ( ( t13(i  ,j,1)*ua(i  ,j,1)     &
                    +t13(i+1,j,1)*ua(i+1,j,1) )   &
                 + ( t23(i,j  ,1)*va(i,j  ,1)     &
                    +t23(i,j+1,1)*va(i,j+1,1) )   &
                 )*rdz*mf(i,j,1)/rf(i,j,1)
          dissten(i,j,1) = dissten(i,j,1) + max(tem1,0.0)
        enddo
        enddo
      ENDIF
      IF( bbc.eq.3 .and. ipbl.eq.0 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dissten(i,j,1) = dissten(i,j,1) + (ust(i,j)**3)/(karman*znt(i,j))
        enddo
        enddo
      ENDIF

      IF(iturb.eq.1)THEN
        do j=0,nj+1
        do i=0,ni+1
          tke3d(i,j,1) = tkea(i,j,1)
        enddo
        enddo
      ENDIF

      if(timestats.ge.1) time_turb=time_turb+mytime()

    ENDIF  ! endif section for iturb.ge.1

      end subroutine sfc_and_turb

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine tkekm(nstep,dt,ruh,rvh,rmh,zf,mf,rmf,znt,ust,rf,           &
                         nm,defv,defh,tk,lenscl,lenh,grdscl,rgrdscl,        &
                         kmh,kmv,khh,khv,tkea,ua,va,dissten,                &
                         nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                   &
                         khcw1,khcw2,khce1,khce2,khcs1,khcs2,khcn1,khcn2,   &
                         khdw1,khdw2,khde1,khde2,khds1,khds2,khdn1,khdn2,   &
                         kvcw1,kvcw2,kvce1,kvce2,kvcs1,kvcs2,kvcn1,kvcn2,   &
                         kvdw1,kvdw2,kvde1,kvde2,kvds1,kvds2,kvdn1,kvdn2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      integer, intent(in) :: nstep
      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: ruh
      real, intent(in), dimension(jb:je) :: rvh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, intent(in),dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf,rmf
      real, intent(in), dimension(ib:ie,jb:je) :: znt,ust
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tk,lenscl,lenh,grdscl,rgrdscl
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: dissten
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmt) :: khcw1,khcw2,khce1,khce2
      real, intent(inout), dimension(imp,kmt) :: khcs1,khcs2,khcn1,khcn2
      real, intent(inout), dimension(jmp,kmt) :: khdw1,khdw2,khde1,khde2
      real, intent(inout), dimension(imp,kmt) :: khds1,khds2,khdn1,khdn2
      real, intent(inout), dimension(jmp,kmt) :: kvcw1,kvcw2,kvce1,kvce2
      real, intent(inout), dimension(imp,kmt) :: kvcs1,kvcs2,kvcn1,kvcn2
      real, intent(inout), dimension(jmp,kmt) :: kvdw1,kvdw2,kvde1,kvde2
      real, intent(inout), dimension(imp,kmt) :: kvds1,kvds2,kvdn1,kvdn2

!----------------------------------------

      integer :: i,j,k
      real :: prinv,tem1,tem2

#ifdef MPI
      integer reqs_khc(8)
      integer reqs_kvc(8)
      integer reqs_khd(8)
      integer reqs_kvd(8)
#endif

!------------------------------------------------------------------
!  Get length scales:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=2,nk

    !  get grid scale
    IF(tconfig.eq.1)THEN
      ! single length scale:  appropriate if dx,dy are nearly the same as dz
      do j=1,nj
      do i=1,ni
        grdscl(i,j,k)=( ((dx*ruh(i))*(dy*rvh(j)))*(dz*rmf(i,j,k)) )**0.33333333
        ! cm1r17:  wall condition near surface
        grdscl(i,j,k) = sqrt(1.0/( 1.0/(grdscl(i,j,k)**2)                                  &
                                  +1.0/((karman*((zf(i,j,k)-zf(i,j,1))+znt(i,j))*rcs)**2)  &
                               ) )
        rgrdscl(i,j,k)=1.0/grdscl(i,j,k)
      enddo
      enddo
    ELSEIF(tconfig.eq.2)THEN
      ! two length scales:  one for horizontal, one for vertical
      do j=1,nj
      do i=1,ni
        lenh(i,j,k)=sqrt( (dx*ruh(i))*(dy*rvh(j)) )
      enddo
      enddo
      do j=1,nj
      do i=1,ni
        grdscl(i,j,k)=dz*rmf(i,j,k)
        ! cm1r17:  wall condition near surface
        grdscl(i,j,k) = sqrt(1.0/( 1.0/(grdscl(i,j,k)**2)                                  &
                                  +1.0/((karman*((zf(i,j,k)-zf(i,j,1))+znt(i,j))*rcs)**2)  &
                               ) )
        rgrdscl(i,j,k)=1.0/grdscl(i,j,k)
      enddo
      enddo
    ENDIF

      ! Get turbulence length scale
      do j=1,nj
      do i=1,ni
        tk(i,j,k)=max(tkea(i,j,k),1.0e-6)
        lenscl(i,j,k)=grdscl(i,j,k)
        if(nm(i,j,k).gt.1.0e-6.and.tk(i,j,k).ge.1.0e-3)then
          lenscl(i,j,k)=c_l*sqrt(tk(i,j,k)/nm(i,j,k))
          lenscl(i,j,k)=min(lenscl(i,j,k),grdscl(i,j,k))
          lenscl(i,j,k)=max(lenscl(i,j,k),1.0e-6*grdscl(i,j,k))
        endif 
      enddo
      enddo

    ENDDO

    if( nstep.eq.0 .and. myid.eq.0 )then
      print *
      print *,'  zf,grdscl:'
      i = 1
      j = 1
      do k=2,nk
        print *,k,(zf(i,j,k)-zf(i,j,1)),grdscl(i,j,k)
      enddo
    endif

!------------------------------------------------------------------

      tem1 = 0.125*dx*dx/dt
      tem2 = 0.125*dy*dy/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,k,prinv)
    DO k=2,nk

    !  Get km, kh
    IF(tconfig.eq.1)THEN

      do j=1,nj
      do i=1,ni
        kmh(i,j,k)=c_m*sqrt(tk(i,j,k))*lenscl(i,j,k)
        kmv(i,j,k)=kmh(i,j,k)
        prinv=3.00
        if(nm(i,j,k).gt.1.0e-6)then
          prinv=min(1.0+2.00*lenscl(i,j,k)*rgrdscl(i,j,k),3.00)
        endif
        khh(i,j,k)=kmh(i,j,k)*prinv
        khv(i,j,k)=khh(i,j,k)
      enddo
      enddo

    ELSEIF(tconfig.eq.2)THEN

      do j=1,nj
      do i=1,ni
        kmh(i,j,k)=c_m*sqrt(tk(i,j,k))*lenh(i,j,k)
        kmv(i,j,k)=c_m*sqrt(tk(i,j,k))*lenscl(i,j,k)
        prinv=3.00
        if(nm(i,j,k).gt.1.0e-6)then
          prinv=min(1.0+2.00*lenscl(i,j,k)*rgrdscl(i,j,k),3.00)
        endif
        khh(i,j,k)=kmh(i,j,k)*prinv
        khv(i,j,k)=kmv(i,j,k)*prinv
      enddo
      enddo

    ENDIF

      !  limit for numerical stability:
      do j=1,nj
      do i=1,ni
        kmh(i,j,k) = min( kmh(i,j,k) , tem1*ruh(i)*ruh(i) , tem2*rvh(j)*rvh(j) )
        khh(i,j,k) = min( khh(i,j,k) , tem1*ruh(i)*ruh(i) , tem2*rvh(j)*rvh(j) )
      enddo
      enddo

    ENDDO

      if(timestats.ge.1) time_turb=time_turb+mytime()

!------------------------------------------------------------
! Set values at boundaries, start comms:

      call bcw(kmh,1)
#ifdef MPI
      call comm_1t_start(kmh,khcw1,khcw2,khce1,khce2,   &
                             khcs1,khcs2,khcn1,khcn2,reqs_khc)
#endif

      call bcw(kmv,1)
#ifdef MPI
      call comm_1t_start(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                             kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
#endif

      call bcw(khh,1)
#ifdef MPI
      call comm_1t_start(khh,khdw1,khdw2,khde1,khde2,   &
                             khds1,khds2,khdn1,khdn2,reqs_khd)
#endif

      call bcw(khv,1)
#ifdef MPI
      call comm_1t_start(khv,kvdw1,kvdw2,kvde1,kvde2,   &
                             kvds1,kvds2,kvdn1,kvdn2,reqs_kvd)
#endif

!--------------------------------------------------------------
!  Dissipation:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=2,nk
      do j=1,nj
      do i=1,ni
        dissten(i,j,k) = dissten(i,j,k)                           &
                     +(c_e1+c_e2*lenscl(i,j,k)*rgrdscl(i,j,k))    &
                     *tk(i,j,k)*sqrt(tk(i,j,k))/lenscl(i,j,k)
      enddo
      enddo
    ENDDO
    if(timestats.ge.1) time_turb=time_turb+mytime()

!--------------------------------------------------------------
!  Finish comms:
#ifdef MPI
      call comm_1t_end(kmh,khcw1,khcw2,khce1,khce2,   &
                           khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call bct2(kmh)
      call comm_1t_end(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                           kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call bct2(kmv)
      call comm_1t_end(khh,khdw1,khdw2,khde1,khde2,   &
                           khds1,khds2,khdn1,khdn2,reqs_khd)
      call comm_1t_end(khv,kvdw1,kvdw2,kvde1,kvde2,   &
                           kvds1,kvds2,kvdn1,kvdn2,reqs_kvd)
      call getcornert(kmh,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call getcornert(kmv,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
#endif

!--------------------------------------------------------------
!  cm1r18: surface

      IF( bbc.eq.3 )THEN

        tem1 = (c_s/c_m)**2

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          tkea(i,j,1) = tem1*ust(i,j)*ust(i,j)
          kmh(i,j,1) = karman*znt(i,j)*ust(i,j)
        enddo
        enddo

        !-----
        call bc2d(tkea(ibt,jbt,1))
#ifdef MPI
        call comm_1s2d_start(tkea(ibt,jbt,1),kvcw1(1,1),kvcw2(1,1),kvce1(1,1),kvce2(1,1),   &
                                             kvcs1(1,1),kvcs2(1,1),kvcn1(1,1),kvcn2(1,1),reqs_kvc)
#endif
        call bc2d(kmh(ibc,jbc,1))
#ifdef MPI
        call comm_1s2d_start(kmh(ibc,jbc,1),khcw1(1,1),khcw2(1,1),khce1(1,1),khce2(1,1),   &
                                            khcs1(1,1),khcs2(1,1),khcn1(1,1),khcn2(1,1),reqs_khc)
        !-----
        call comm_1s2d_end(tkea(ibt,jbt,1),kvcw1(1,1),kvcw2(1,1),kvce1(1,1),kvce2(1,1),   &
                                           kvcs1(1,1),kvcs2(1,1),kvcn1(1,1),kvcn2(1,1),reqs_kvc)
        call bcs2_2d(tkea(ibt,jbt,1))
        call comm_1s2d_end(kmh(ibc,jbc,1),khcw1(1,1),khcw2(1,1),khce1(1,1),khce2(1,1),   &
                                          khcs1(1,1),khcs2(1,1),khcn1(1,1),khcn2(1,1),reqs_khc)
        call bcs2_2d(kmh(ibc,jbc,1))
        !-----
        call comm_2d_corner(tkea(ibt,jbt,1))
        call comm_2d_corner(kmh(ibc,jbc,1))
        !-----
#endif

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          kmv(i,j,1) = kmh(i,j,1)
          khv(i,j,1) = kmv(i,j,1)*(khv(i,j,2)/(1.0e-10+kmv(i,j,2)))
          khh(i,j,1) = kmh(i,j,1)*(khh(i,j,2)/(1.0e-10+kmh(i,j,2)))
        enddo
        enddo

      ENDIF

!--------------------------------------------------------------
!  finished
      
      return
      end subroutine tkekm


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbsmag(nstep,dt,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,rf, &
                          nm,defv,defh,lenscl,grdscl,lenh,zf,znt,ust,      &
                          kmh,kmv,khh,khv,ua,va,dissten,                   &
                          nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                          khcw1,khcw2,khce1,khce2,khcs1,khcs2,khcn1,khcn2, &
                          kvcw1,kvcw2,kvce1,kvce2,kvcs1,kvcs2,kvcn1,kvcn2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      integer, intent(in) :: nstep
      real, intent(in) :: dt
      logical, intent(in) :: dosfcflx
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, dimension(ib:ie,jb:je,kb:ke) :: th0
      real, dimension(ib:ie,jb:je) :: thflux,qvflux,rth0s
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf
      real, dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh
      real, dimension(ib:ie,jb:je,kb:ke) :: lenscl,grdscl,lenh
      real, intent(in), dimension(ib:ie,jb:je) :: znt,ust
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, dimension(ib:ie,jb:je,kb:ke+1) :: dissten
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmt) :: khcw1,khcw2,khce1,khce2
      real, intent(inout), dimension(imp,kmt) :: khcs1,khcs2,khcn1,khcn2
      real, intent(inout), dimension(jmp,kmt) :: kvcw1,kvcw2,kvce1,kvce2
      real, intent(inout), dimension(imp,kmt) :: kvcs1,kvcs2,kvcn1,kvcn2

      integer i,j,k
      real :: tem,temx,temy

#ifdef MPI
      integer reqs_khc(8)
      integer reqs_kvc(8)
#endif

      real, parameter :: cs      = 0.18
      real, parameter :: csinv   = 1.0/cs
      real, parameter :: prandtl = 1.0/3.00
      real, parameter :: prinv   = 1.0/prandtl
      real, parameter :: dmin    = 1.0e-10

!-----------------------------------------------------------------------

      temx = 0.125*dx*dx/dt
      temy = 0.125*dy*dy/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    jloop:  DO j=1,nj

    IF(tconfig.eq.1)THEN
      ! single length scale:  appropriate if dx,dy are nearly the same as dz

      do k=2,nk
      do i=1,ni
        grdscl(i,j,k)=( ((dx*ruh(i))*(dy*rvh(j)))*(dz*rmf(i,j,k)) )**0.33333333
        grdscl(i,j,k) = sqrt(1.0/( 1.0/(grdscl(i,j,k)**2)                                    &
                                  +1.0/((karman*((zf(i,j,k)-zf(i,j,1))+znt(i,j))*csinv)**2)  &
                               ) )
      enddo
      enddo

    ELSEIF(tconfig.eq.2)THEN
      ! two length scales:  one for horizontal, one for vertical

      do i=1,ni
        tem=sqrt( (dx*ruh(i))*(dy*rvh(j)) )
        do k=2,nk
          lenh(i,j,k)=tem
        enddo
      enddo

      do k=2,nk
      do i=1,ni
        grdscl(i,j,k)=dz*rmf(i,j,k)
        grdscl(i,j,k) = sqrt(1.0/( 1.0/(grdscl(i,j,k)**2)                                    &
                                  +1.0/((karman*((zf(i,j,k)-zf(i,j,1))+znt(i,j))*csinv)**2)  &
                               ) )
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------------

    IF(tconfig.eq.1)THEN

      do k=2,nk
      do i=1,ni
        kmh(i,j,k)=((cs*grdscl(i,j,k))**2)     &
                 *sqrt( max(defv(i,j,k)+defh(i,j,k)-nm(i,j,k)*prinv,dmin) )
        kmh(i,j,k) = min( kmh(i,j,k) , temx*ruh(i)*ruh(i)   &
                                     , temy*rvh(j)*rvh(j) )
        kmv(i,j,k)=kmh(i,j,k)
      enddo
      enddo

    ELSEIF(tconfig.eq.2)THEN

      do k=2,nk
      do i=1,ni
        kmh(i,j,k)=((cs*lenh(i,j,k))**2)     &
                 *sqrt( max(defh(i,j,k),dmin) )
        kmh(i,j,k) = min( kmh(i,j,k) , temx*ruh(i)*ruh(i)   &
                                     , temy*rvh(j)*rvh(j) )
        kmv(i,j,k)=((cs*grdscl(i,j,k))**2)     &
                 *sqrt( max(defv(i,j,k)-nm(i,j,k)*prinv,dmin) )
      enddo
      enddo

    ENDIF

    ENDDO  jloop

    if( nstep.eq.0 .and. myid.eq.0 )then
      print *
      print *,'  cs,csinv = ',cs,csinv
      print *,'  zf,grdscl:'
      i = 1
      j = 1
      do k=2,nk
        print *,k,(zf(i,j,k)-zf(i,j,1)),grdscl(i,j,k)
      enddo
    endif

!--------------------------------------------------------------

      if(timestats.ge.1) time_turb=time_turb+mytime()
      call bcw(kmh,1)
      call bcw(kmv,1)
#ifdef MPI
      call comm_1t_start(kmh,khcw1,khcw2,khce1,khce2,   &
                             khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call comm_1t_start(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                             kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call comm_1t_end(kmh,khcw1,khcw2,khce1,khce2,   &
                           khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call comm_1t_end(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                           kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call getcornert(kmh,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call getcornert(kmv,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call bct2(kmh)
      call bct2(kmv)
#endif

!--------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    do j=0,nj+1

      do k=1,nk+1
      do i=0,ni+1
        khh(i,j,k)=kmh(i,j,k)*prinv
        ! limit for numerical stability:
        khh(i,j,k) = min( khh(i,j,k) , temx*ruh(i)*ruh(i)   &
                                     , temy*rvh(j)*rvh(j) )
        khv(i,j,k)=kmv(i,j,k)*prinv
      enddo
      enddo

    IF( idiss.eq.1 .or. output_dissten.eq.1 )THEN
    IF( j.ge.1 .and. j.le.nj )THEN
    IF( tconfig.eq.1 )THEN
      do k=2,nk
      do i=1,ni
        dissten(i,j,k) = dissten(i,j,k) + (kmv(i,j,k)**3)/((cs*grdscl(i,j,k))**4)
      enddo
      enddo
    ELSEIF( tconfig.eq.2 )THEN
      do k=2,nk
      do i=1,ni
        dissten(i,j,k) = dissten(i,j,k) + (kmv(i,j,k)**3)/((cs*grdscl(i,j,k))**4)    &
                                        + (kmh(i,j,k)**3)/((cs*lenh(i,j,k))**4)
      enddo
      enddo
    ENDIF
    ENDIF
    ENDIF

    enddo

!--------------------------------------------------------------
!  cm1r18: surface

      IF( bbc.eq.3 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          kmv(i,j,1) = karman*znt(i,j)*ust(i,j)
        enddo
        enddo

        call bc2d(kmv(ibc,jbc,1))
#ifdef MPI
        call comm_1s2d_start(kmv(ibc,jbc,1),khcw1(1,1),khcw2(1,1),khce1(1,1),khce2(1,1),   &
                                            khcs1(1,1),khcs2(1,1),khcn1(1,1),khcn2(1,1),reqs_khc)
        call comm_1s2d_end(kmv(ibc,jbc,1),khcw1(1,1),khcw2(1,1),khce1(1,1),khce2(1,1),   &
                                          khcs1(1,1),khcs2(1,1),khcn1(1,1),khcn2(1,1),reqs_khc)
        call comm_2d_corner(kmv(ibc,jbc,1))
        call bcs2_2d(kmv(ibc,jbc,1))
#endif

      IF( tconfig.eq.1 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          kmh(i,j,1) = kmv(i,j,1)
          khv(i,j,1) = kmv(i,j,1)*prinv
          khh(i,j,1) = khv(i,j,1)
        enddo
        enddo
      ELSEIF( tconfig.eq.2 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          khv(i,j,1) = kmv(i,j,1)*prinv
        enddo
        enddo
      ENDIF

      ENDIF

!--------------------------------------------------------------

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end subroutine turbsmag


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbparam(nstep,zf,dt,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,rf, &
                          nm,defv,defh,lvz,kmh,kmv,khh,khv,ua,va,dissten,zs,znt,ust,xland,psfc,tlh,  &
                          nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                         &
                          khcw1,khcw2,khce1,khce2,khcs1,khcs2,khcn1,khcn2,         &
                          kvcw1,kvcw2,kvce1,kvce2,kvcs1,kvcs2,kvcn1,kvcn2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      integer, intent(in) :: nstep
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in) :: dt
      logical, intent(in) :: dosfcflx
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, dimension(ib:ie,jb:je,kb:ke) :: th0
      real, dimension(ib:ie,jb:je) :: thflux,qvflux,rth0s
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf
      real, dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh
      real, dimension(ib:ie,jb:je,kb:ke) :: lvz
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, dimension(ib:ie,jb:je,kb:ke+1) :: dissten
      real, intent(in), dimension(ib:ie,jb:je) :: zs,znt,ust,xland,psfc
      real, intent(inout), dimension(ib:ie,jb:je) :: tlh
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmt) :: khcw1,khcw2,khce1,khce2
      real, intent(inout), dimension(imp,kmt) :: khcs1,khcs2,khcn1,khcn2
      real, intent(inout), dimension(jmp,kmt) :: kvcw1,kvcw2,kvce1,kvce2
      real, intent(inout), dimension(imp,kmt) :: kvcs1,kvcs2,kvcn1,kvcn2

      integer i,j,k
      real :: rlinf,tem,tem1,temx,temy

#ifdef MPI
      integer reqs_khc(8)
      integer reqs_kvc(8)
#endif

      real, parameter :: prandtl = 1.0
      real, parameter :: prinv   = 1.0/prandtl
      real, parameter :: dmin    = 1.0e-10

!--------------------------------------------------------------
!  Smagorinsky-type scheme for parameterized turbulence:
!--------------------------------------------------------------
!  Interior:

!!!    tem = 1.0/(1.0e-6+l_inf)
    rlinf = (1.0e-6+l_inf)**(-2)
    if(ny.eq.1)then
      temx =  0.250*dx*dx/dt
      temy = 1000.0*dy*dy/dt
    elseif(nx.eq.1)then
      temx = 1000.0*dx*dx/dt
      temy =  0.250*dy*dy/dt
    else
      temx =  0.125*dx*dx/dt
      temy =  0.125*dy*dy/dt
    endif

  IF( l_h.gt.1.0e-12 .or. lhref1.gt.1.0e-12 .or. lhref2.gt.1.0e-12 )THEN
    ! cm1r18:
    ! Over water, make tlh a function of surface pressure.
    !   (designed for hurricanes)
    ! Over land, simply set to tlh to l_h.
!$omp parallel do default(shared)   &
!$omp private(i,j)
    do j=1,nj
    do i=1,ni
      IF( (xland(i,j).gt.1.5) .and. (zs(i,j).lt.1.0) )THEN
        ! over water (sea level only):
        tlh(i,j) = lhref2+(lhref1-lhref2)   &
                      *(psfc(i,j)-90000.0)  &
                      /( 101500.0-90000.0)
      ELSE
        ! all other cases:
        tlh(i,j) = l_h
      ENDIF
    enddo
    enddo
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      kmh(i,j,k)=(tlh(i,j)**2)*sqrt( max(defh(i,j,k),dmin) )
      kmh(i,j,k) = min( kmh(i,j,k) , temx*ruh(i)*ruh(i) , temy*rvh(j)*rvh(j) )
    enddo
    enddo
    enddo
  ENDIF
  IF( l_inf.gt.1.0e-12 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      lvz(i,j,k)=sqrt( ( rlinf + (karman*((zf(i,j,k)-zf(i,j,1))+znt(i,j)))**(-2) )**(-1) )
      kmv(i,j,k)=(lvz(i,j,k)**2)*sqrt( max(defv(i,j,k)-nm(i,j,k)*prinv,dmin) )
    enddo
    enddo
    enddo
  ENDIF

!--------------------------------------------------------------
! boundary conditions:

      if(timestats.ge.1) time_turb=time_turb+mytime()

    IF( l_h.gt.1.0e-12 .or. lhref1.gt.1.0e-12 .or. lhref2.gt.1.0e-12 )THEN
      call bcw(kmh,1)
#ifdef MPI
      call comm_1t_start(kmh,khcw1,khcw2,khce1,khce2,   &
                             khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call comm_1t_end(kmh,khcw1,khcw2,khce1,khce2,   &
                           khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call getcornert(kmh,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call bct2(kmh)
#endif
    ENDIF
    IF( l_inf.gt.1.0e-12 )THEN
      call bcw(kmv,1)
#ifdef MPI
      call comm_1t_start(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                             kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call comm_1t_end(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                           kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call getcornert(kmv,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call bct2(kmv)
#endif
    ENDIF

        do j=0,nj+1
        do i=0,ni+1
          kmh(i,j,1) = 2.0*kmh(i,j,2) - kmh(i,j,3)
        enddo
        enddo

!--------------------------------------------------------------
!  calculate Kh
!  and also limit horizontal coeffs for numerical stability:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem1)
    do j=0,nj+1

    IF( l_h.gt.1.0e-12 .or. lhref1.gt.1.0e-12 .or. lhref2.gt.1.0e-12 )THEN
      do k=1,nk+1
      do i=0,ni+1
        khh(i,j,k)=kmh(i,j,k)*prinv
        khh(i,j,k) = min( khh(i,j,k) , temx*ruh(i)*ruh(i) , temy*rvh(j)*rvh(j) )
      enddo
      enddo
    ENDIF
    IF( l_inf.gt.1.0e-12 )THEN
      do k=1,nk+1
      do i=0,ni+1
        khv(i,j,k)=kmv(i,j,k)*prinv
      enddo
      enddo
    ENDIF

    IF( idiss.eq.1 .or. output_dissten.eq.1 )THEN
    IF( j.ge.1 .and. j.le.nj )THEN
    IF( l_h.gt.1.0e-12 .or. lhref1.gt.1.0e-12 .or. lhref2.gt.1.0e-12 )THEN
      do k=2,nk
      do i=1,ni
        dissten(i,j,k) = dissten(i,j,k) + (kmh(i,j,k)**3)/(tlh(i,j)**4)
      enddo
      enddo
    ENDIF
    IF( l_inf.gt.1.0e-12 )THEN
      do k=2,nk
      do i=1,ni
        dissten(i,j,k) = dissten(i,j,k) + (kmv(i,j,k)**3)/(lvz(i,j,k)**4)
      enddo
      enddo
    ENDIF
    ENDIF
    ENDIF

    enddo

!--------------------------------------------------------------
!  cm1r18: surface

      IF( bbc.eq.3 .and. l_inf.gt.1.0e-12 )THEN

        do j=1,nj
        do i=1,ni
          kmv(i,j,1) = karman*znt(i,j)*ust(i,j)
        enddo
        enddo

        call bc2d(kmv(ibc,jbc,1))
#ifdef MPI
        call comm_1s2d_start(kmv(ibc,jbc,1),khcw1(1,1),khcw2(1,1),khce1(1,1),khce2(1,1),   &
                                            khcs1(1,1),khcs2(1,1),khcn1(1,1),khcn2(1,1),reqs_khc)
        call comm_1s2d_end(kmv(ibc,jbc,1),khcw1(1,1),khcw2(1,1),khce1(1,1),khce2(1,1),   &
                                          khcs1(1,1),khcs2(1,1),khcn1(1,1),khcn2(1,1),reqs_khc)
        call comm_2d_corner(kmv(ibc,jbc,1))
        call bcs2_2d(kmv(ibc,jbc,1))
#endif

        do j=0,nj+1
        do i=0,ni+1
          khv(i,j,1) = kmv(i,j,1)*prinv
        enddo
        enddo

      ENDIF

!--------------------------------------------------------------

      if(nstep.eq.1.and.myid.eq.0)then
        print *,'  k,zf,lvz:  znt = ',znt(1,1)
        do k=2,nk
          print *,k,(zf(1,1,k)-zf(1,1,1)),lvz(1,1,k)
        enddo
      endif

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end subroutine turbparam


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      ! gettau
      subroutine gettau(xf,rxf,arf1,arf2,ust,u1,v1,s1,rf,  &
                        kmh,kmv,t11,t12,t13,t22,t23,t33,ua)
      implicit none
      
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2
      real, intent(in), dimension(ib:ie,jb:je) :: ust,u1,v1,s1
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv
      real, dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
        
      integer i,j,k
      real :: tem

!----------------------------------------------------------------------
!
!  This subroutine calculates the subgrid stress terms.
!
!    t_ij  =  2 * rho * K * S_ij
!
!  NOTE:  upon entering this subroutine, the t_ij arrays must already 
!         contain rho * S_ij  (see calcdef subroutine)
!
!  Since cm1r18, surface stress (ie, surface drag) is incorporated into
!  the stress arrays here.
!
!  Note:  Turbulent viscosities are defined on w points.
!
!  Note:  For axisymmetric simulations, t11 and t12 herein are 
!         actually not stresses:  the actual stresses are
!         combined in a convienent form for the sake of flux-form
!         calculations in the turbu and turbv subroutines.
!         Also note that t22 is never calculated.
!         So, if you need the actual stress components for something, 
!         beware that you will need to re-calculate t11,t12,t22.
!
!----------------------------------------------------------------------

  IF(axisymm.eq.0)THEN

    ! Cartesian grid:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem)
    do k=1,nk

      do j=0,nj+1
      do i=0,ni+1
        !  2.0 * 0.5 = 1.0
        tem = (kmh(i,j,k)+kmh(i,j,k+1))
        t11(i,j,k)=t11(i,j,k)*tem
        t22(i,j,k)=t22(i,j,k)*tem
        t33(i,j,k)=t33(i,j,k)*tem
      enddo
      enddo

      do j=1,nj+1
      do i=1,ni+1
        !  2.0 * 0.125 = 0.25
        t12(i,j,k)=t12(i,j,k)*0.25                                            &
     *( ( (kmh(i-1,j-1,k  )+kmh(i,j,k  ))+(kmh(i-1,j,k  )+kmh(i,j-1,k  )) )   &
       +( (kmh(i-1,j-1,k+1)+kmh(i,j,k+1))+(kmh(i-1,j,k+1)+kmh(i,j-1,k+1)) ) )
      enddo
      enddo
          !-----
          ! lateral boundary conditions:
          if(wbc.eq.3.and.ibw.eq.1)then
            ! free slip b.c.
            do j=1,nj+1
              t12(1,j,k) = t12(2,j,k)
            enddo
          endif
          if(ebc.eq.3.and.ibe.eq.1)then
            ! free slip b.c.
            do j=1,nj+1
              t12(ni+1,j,k) = t12(ni,j,k)
            enddo
          endif
          !-----
          !-----
          if(sbc.eq.3.and.ibs.eq.1)then
            ! free slip b.c.
            do i=1,ni+1
              t12(i,1,k) = t12(i,2,k)
            enddo
          endif
          if(nbc.eq.3.and.ibn.eq.1)then
            ! free slip b.c.
            do i=1,ni+1
              t12(i,nj+1,k) = t12(i,nj,k)
            enddo
          endif
          !-----

    IF(k.ge.2)THEN
      do j=1,nj+1
      do i=1,ni+1
        !  2.0 x 0.5 = 1.0
        t13(i,j,k)=t13(i,j,k)*( kmv(i-1,j,k)+kmv(i,j,k) )
        t23(i,j,k)=t23(i,j,k)*( kmv(i,j-1,k)+kmv(i,j,k) )
      enddo
      enddo
            !-----
            ! lateral boundary conditions:
            if(wbc.eq.3.and.ibw.eq.1)then
              ! free slip b.c.
              do j=1,nj
                t13(1,j,k) = t13(2,j,k)
              enddo
            endif
            if(ebc.eq.3.and.ibe.eq.1)then
              ! free slip b.c.
              do j=1,nj
                t13(ni+1,j,k) = t13(ni,j,k)
              enddo
            endif
            !-----
            !-----
            if(sbc.eq.3.and.ibs.eq.1)then
              ! free slip b.c.
              do i=1,ni
                t23(i,1,k) = t23(i,2,k)
              enddo
            endif
            if(nbc.eq.3.and.ibn.eq.1)then
              ! free slip b.c.
              do i=1,ni
                t23(i,nj+1,k) = t23(i,nj,k)
              enddo
            endif
            !-----
    ENDIF

    enddo

!------------------------------------

  ELSE

    ! axisymmetric grid:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem)
    DO k=1,nk

      do j=1,nj
      do i=1,ni+1
        !  2.0 * 0.5 = 1.0
        tem = (kmh(i,j,k)+kmh(i,j,k+1))
        t11(i,j,k)=t11(i,j,k)*tem
        t33(i,j,k)=t33(i,j,k)*tem
        !  2.0 * 0.25  =  0.5
        t12(i,j,k)=0.5*t12(i,j,k)*( arf2(i)*(kmh(i  ,j,k+1)+kmh(i  ,j,k)) &
                                   +arf1(i)*(kmh(i-1,j,k+1)+kmh(i-1,j,k)) )
      enddo
      enddo
          !-----
          ! lateral boundary conditions:
          j = 1
          if(wbc.eq.3)then
            ! free slip b.c.
!!!            t12(1,j,k) = t12(2,j,k)
            t12(1,j,k) = 0.0
          endif
          if(ebc.eq.3)then
            ! free slip b.c.
            t12(ni+1,j,k) = t12(ni,j,k)
          endif
          !-----
    IF(k.ge.2)THEN
      do j=1,nj
      do i=1,ni+1
        !  2.0 * 0.5  =  1.0
        t13(i,j,k)=t13(i,j,k)*( arf1(i)*kmv(i-1,j,k)+arf2(i)*kmv(i,j,k) )
        t23(i,j,k)=2.0*t23(i,j,k)*kmv(i,j,k)
      enddo
      enddo
    ENDIF

    ENDDO

  ENDIF

!------------------------------------------------------------------
!  open boundary conditions:

    IF( wbc.eq.2 .or. ebc.eq.2 .or. sbc.eq.2 .or. nbc.eq.2 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
        !-----
        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=0,nj+1
            t11(0,j,k) = t11(1,j,k)
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=0,nj+1
            t11(ni+1,j,k) = t11(ni,j,k)
          enddo
        ENDIF
        !-----
        !ccccc
        !-----
        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=0,ni+1
            t22(i,0,k) = t22(i,1,k)
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=0,ni+1
            t22(i,nj+1,k) = t22(i,nj,k)
          enddo
        ENDIF
        !-----
        !ccccc
        !-----
        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj+1
            t12(1,j,k) = t12(2,j,k)
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj+1
            t12(ni+1,j,k) = t12(ni,j,k)
          enddo
        ENDIF
        !-----
        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni+1
            t12(i,1,k) = t12(i,2,k)
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni+1
            t12(i,nj+1,k) = t12(i,nj,k)
          enddo
        ENDIF
        !-----
        ! corner points:
        !-----
        IF( sbc.eq.2 .and. ibs.eq.1 .and. &
            wbc.eq.2 .and. ibw.eq.1 )THEN
          t12(1,1,k) = t12(2,2,k)
        ENDIF
        IF( sbc.eq.2 .and. ibs.eq.1 .and. &
            ebc.eq.2 .and. ibe.eq.1 )THEN
          t12(ni+1,1,k) = t12(ni,2,k)
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 .and. &
            wbc.eq.2 .and. ibw.eq.1 )THEN
          t12(1,nj+1,k) = t12(2,nj,k)
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 .and. &
            ebc.eq.2 .and. ibe.eq.1 )THEN
          t12(ni+1,nj+1,k) = t12(ni,nj,k)
        ENDIF
        !-----
      ENDDO
    ENDIF

!--------------------------------------------------------------
!  lower boundary conditions

    IF(bbc.eq.1)THEN
      ! free slip:

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni+1
        t13(i,j,1)=t13(i,j,2)
        t23(i,j,1)=t23(i,j,2)
      enddo
      enddo

    ELSEIF(bbc.eq.2)THEN
      ! no slip:

      IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni+1
          t13(i,j,1)=t13(i,j,1)*0.5*( kmv(i-1,j,2)+kmv(i,j,2) )
          t23(i,j,1)=t23(i,j,1)*0.5*( kmv(i,j-1,2)+kmv(i,j,2) )
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni+1
          t13(i,j,1)=t13(i,j,1)*0.5*( arf1(i)*kmv(i-1,j,2)+arf2(i)*kmv(i,j,2) )
          t23(i,j,1)=t23(i,j,1)*kmv(i,j,2)
        enddo
        enddo

      ENDIF

    ELSEIF(bbc.eq.3)THEN
      !--------------------------------------------------------!
      !--------  surface stress for semi-slip lower bc --------!
      !-------- (this is where "drag" is set for bbc=3) -------!

      IF(axisymm.eq.0)THEN

        ! Cartesian grid:
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni+1
        t13(i,j,1) = 0.25*( (ust(i-1,j)**2)*(u1(i-1,j)/max(s1(i-1,j),0.01))    &
                           +(ust(i  ,j)**2)*(u1(i  ,j)/max(s1(i  ,j),0.01)) )  &
                         *( rf(i-1,j,1)+rf(i,j,1) )
        t23(i,j,1) = 0.25*( (ust(i,j-1)**2)*(v1(i,j-1)/max(s1(i,j-1),0.01))    &
                           +(ust(i,j  )**2)*(v1(i,j  )/max(s1(i,j  ),0.01)) )  &
                         *( rf(i,j-1,1)+rf(i,j,1) )
      enddo
      enddo

      ELSE

        ! axisymmetric grid:
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni+1
        t13(i,j,1) = 0.25*( arf1(i)*(ust(i-1,j)**2)*(u1(i-1,j)/max(s1(i-1,j),0.01))    &
                           +arf2(i)*(ust(i  ,j)**2)*(u1(i  ,j)/max(s1(i  ,j),0.01)) )  &
                         *(arf1(i)*rf(i-1,j,1)+arf2(i)*rf(i,j,1))
        t23(i,j,1) = rf(i,j,1)*(ust(i,j)**2)*(v1(i,j)/max(s1(i,j),0.01))
      enddo
      enddo

      ENDIF

    ENDIF

!--------------------------------------------------------------
!  upper boundary conditions

    IF(tbc.eq.1)THEN
      ! free slip:

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni+1
        t13(i,j,nk+1)=t13(i,j,nk)
        t23(i,j,nk+1)=t23(i,j,nk)
      enddo
      enddo

    ELSEIF(tbc.eq.2)THEN
      ! no slip:

      IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni+1
          t13(i,j,nk+1)=t13(i,j,nk+1)*0.5*( kmv(i-1,j,nk)+kmv(i,j,nk) )
          t23(i,j,nk+1)=t23(i,j,nk+1)*0.5*( kmv(i,j-1,nk)+kmv(i,j,nk) )
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni+1
          t13(i,j,nk+1)=t13(i,j,nk+1)*0.5*( arf1(i)*kmv(i-1,j,nk)+arf2(i)*kmv(i,j,nk) )
          t23(i,j,nk+1)=t23(i,j,nk+1)*kmv(i,j,nk)
        enddo
        enddo

      ENDIF

    ENDIF

!--------------------------------------------------------------

    IF( axisymm.eq.1 )THEN
      ! lateral boundary condition:
!$omp parallel do default(shared)   &
!$omp private(k)
      do k=0,nk+1
        t13(1,1,k)=0.0
      enddo
    ENDIF

!--------------------------------------------------------------
!  finished

      if(timestats.ge.1) time_turb=time_turb+mytime()
 
      return
      end subroutine gettau
      ! gettau


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      ! calcdef
      subroutine calcdef(rds,sigma,rdsf,sigmaf,zs,gz,rgz,gzu,rgzu,gzv,rgzv,                &
                     xh,rxh,arh1,arh2,uh,xf,rxf,arf1,arf2,uf,vh,vf,mh,c1,c2,mf,defv,defh,  &
                     dum1,dum2,ua,va,wa,s11,s12,s13,s22,s23,s33,gx,gy,rho,rr,rf)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: defv,defh
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: s11,s12,s13,s22,s23,s33
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gy
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf
        
      integer :: i,j,k
      real :: r1,r2,r3,r4
      real :: tmp11,tmp22,tmp33,tmp12,tmp13,tmp23,rrf
      real :: temz

!----------------------------------------------------------------------
!
!  This subroutine calculates the strain rate terms
!
!    S_ij  =  0.5 * ( d(u_i)/d(x_j) + d(u_j)/d(x_i) )
!
!  (note: multiplied by density herein)
!  and then uses these variables to calculate deformation.
!
!  Note:
!  Since cm1r18, surface stress (ie, surface drag) is no longer 
!  calculated in this subroutine.  See gettau subroutine instead.
!
!  Note:  For axisymmetric simulations, s11 and s12 herein are 
!         actually not rate-of-strain components:  the actual 
!         components have been combined mathematically in a 
!         way to be consistent with the flux-form calculations 
!         in the turbu and turbv subroutines.
!         Also note that s22 is never calculated.
!         So, if you need the actual strain components for something, 
!         beware that you will need to re-calculate s11,s12,s22.
!
!----------------------------------------------------------------------

  IF(.not.terrain_flag)THEN

  IF( axisymm.eq.0 )THEN
    ! Cartesian without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,temz)
    DO k=1,nk

      temz = rdz*mh(1,1,k)

      do j=0,nj+1
      do i=0,ni+1 
        s11(i,j,k)=rho(i,j,k)*(ua(i+1,j,k)-ua(i,j,k))*rdx*uh(i)
        s22(i,j,k)=rho(i,j,k)*(va(i,j+1,k)-va(i,j,k))*rdy*vh(j)
        s33(i,j,k)=rho(i,j,k)*(wa(i,j,k+1)-wa(i,j,k))*temz
      enddo
      enddo
      do j=1,nj+1 
      do i=1,ni+1
        s12(i,j,k)=0.5*( (ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)   &
                        +(va(i,j,k)-va(i-1,j,k))*rdx*uf(i) ) &
              *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
      enddo
      enddo       
          !-----
          ! lateral boundary conditions:
          if(wbc.eq.3.and.ibw.eq.1)then
            ! free slip b.c.
            do j=1,nj+1
              s12(1,j,k) = s12(2,j,k)
            enddo
          elseif(wbc.eq.4.and.ibw.eq.1)then
            ! no slip b.c.
            i = 1
            do j=1,nj+1
              s12(1,j,k) = 2.0*va(1,j,k)*rdx*uf(1)   &
                   *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
            enddo
          endif
          if(ebc.eq.3.and.ibe.eq.1)then
            ! free slip b.c.
            do j=1,nj+1
              s12(ni+1,j,k) = s12(ni,j,k)
            enddo
          elseif(ebc.eq.4.and.ibe.eq.1)then
            ! no slip b.c.
            i = ni+1
            do j=1,nj+1
              s12(ni+1,j,k) = -2.0*va(ni,j,k)*rdx*uf(ni+1)   &
                   *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
            enddo
          endif
          !-----
          !-----
          if(sbc.eq.3.and.ibs.eq.1)then
            ! free slip b.c.
            do i=1,ni+1
              s12(i,1,k) = s12(i,2,k)
            enddo
          elseif(sbc.eq.4.and.ibs.eq.1)then
            ! no slip b.c.
            j = 1
            do i=1,ni+1
              s12(i,1,k) = 2.0*ua(i,1,k)*rdy*vf(1)   &
                   *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
            enddo
          endif
          if(nbc.eq.3.and.ibn.eq.1)then
            ! free slip b.c.
            do i=1,ni+1
              s12(i,nj+1,k) = s12(i,nj,k)
            enddo
          elseif(nbc.eq.4.and.ibn.eq.1)then
            ! no slip b.c.
            j = nj+1
            do i=1,ni+1
              s12(i,nj+1,k) = -2.0*ua(i,nj,k)*rdy*vf(nj+1)   &
                   *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
            enddo
          endif
          !-----
    IF(k.ge.2)THEN
      do j=1,nj
      do i=1,ni+1
        s13(i,j,k)=0.5*( (wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i)   &
                        +(ua(i,j,k)-ua(i,j,k-1))*rdz*0.5*(mf(i-1,j,k)+mf(i,j,k))  &
                       )*0.5*( rf(i-1,j,k)+rf(i,j,k) )
      enddo
      enddo
            !-----
            ! lateral boundary conditions:
            if(wbc.eq.3.and.ibw.eq.1)then
              ! free slip b.c.
              do j=1,nj
                s13(1,j,k) = s13(2,j,k)
              enddo
            elseif(wbc.eq.4.and.ibw.eq.1)then
              ! no slip b.c.
              do j=1,nj
                s13(1,j,k) = 2.0*wa(1,j,k)*rdx*uf(1)
              enddo
            endif
            if(ebc.eq.3.and.ibe.eq.1)then
              ! free slip b.c.
              do j=1,nj
                s13(ni+1,j,k) = s13(ni,j,k)
              enddo
            elseif(ebc.eq.4.and.ibe.eq.1)then
              ! no slip b.c.
              do j=1,nj
                s13(ni+1,j,k) = -2.0*wa(ni,j,k)*rdx*uf(ni+1)
              enddo
            endif
            !-----
      do j=1,nj+1   
      do i=1,ni
        s23(i,j,k)=0.5*( (wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j)   &
                        +(va(i,j,k)-va(i,j,k-1))*rdz*0.5*(mf(i,j-1,k)+mf(i,j,k))  &
                       )*0.5*( rf(i,j-1,k)+rf(i,j,k) )
      enddo
      enddo
            !-----
            if(sbc.eq.3.and.ibs.eq.1)then
              ! free slip b.c.
              do i=1,ni
                s23(i,1,k) = s23(i,2,k)
              enddo
            elseif(sbc.eq.4.and.ibs.eq.1)then
              ! no slip b.c.
              do i=1,ni
                s23(i,1,k) = 2.0*wa(i,1,k)*rdy*vf(1)
              enddo
            endif
            if(nbc.eq.3.and.ibn.eq.1)then
              ! free slip b.c.
              do i=1,ni
                s23(i,nj+1,k) = s23(i,nj,k)
              enddo
            elseif(nbc.eq.4.and.ibn.eq.1)then
              ! no slip b.c.
              do i=1,ni
                s23(i,nj+1,k) = -2.0*wa(i,nj,k)*rdy*vf(nj+1)
              enddo
            endif
            !-----
    ENDIF  ! endif for k.ge.2

    ENDDO  ! endif for k-loop

!-------------------------------------------------------------------------------

  ELSE
    ! axisymmetric:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk

      do j=1,nj
      do i=1,ni+1
        s11(i,j,k)=rho(i,j,k)*(ua(i+1,j,k)*arf1(i+1)-ua(i,j,k)*arf2(i))*rdx*uh(i)
        s33(i,j,k)=rho(i,j,k)*(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)
        !  0.5 * 0.5  =  0.25
        s12(i,j,k)=0.25*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k))   &
                       *(arh1(i)*va(i,j,k)-arh2(i-1)*va(i-1,j,k))*rdx*uf(i)
      enddo
      enddo
          !-----
          ! lateral boundary conditions:
          j = 1
          if(wbc.eq.3)then
            ! free slip b.c.
!!!            s12(1,j,k) = s12(2,j,k)
            s12(1,j,k) = 0.0
          elseif(wbc.eq.4)then
            ! no slip b.c.
            i = 1
            s12(1,j,k) = 2.0*va(1,j,k)*rdx*uf(1)   &
                      *0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k))
          endif
          if(ebc.eq.3)then
            ! free slip b.c.
            s12(ni+1,j,k) = s12(ni,j,k)
          elseif(ebc.eq.4)then
            ! no slip b.c.
            i = ni+1
            s12(ni+1,j,k) = -2.0*va(ni,j,k)*rdx*uf(ni+1)   &
                      *0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k))
          endif
          !-----
    IF(k.ge.2)THEN
      do j=1,nj
      do i=1,ni+1
        !  0.5 * 0.5  =  0.25
        s13(i,j,k)=0.25*(arf1(i)*rf(i-1,j,k)+arf2(i)*rf(i,j,k))  &
                       *( (ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k)  &
                         +(wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i) )
        s23(i,j,k)=0.5*rf(i,j,k)*(va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)
      enddo
      enddo
    ENDIF  ! endif for k.ge.2

    ENDDO  ! endif for k-loop

  ENDIF

!-------------------------------------------------------------------------------
!  Cartesian with terrain:

  ELSE

    ! dum1 stores u at w-pts:
    ! dum2 stores v at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
    do j=0,nj+2
      ! lowest model level:
      do i=0,ni+2
        dum1(i,j,1) = cgs1*ua(i,j,1)+cgs2*ua(i,j,2)+cgs3*ua(i,j,3)
        dum2(i,j,1) = cgs1*va(i,j,1)+cgs2*va(i,j,2)+cgs3*va(i,j,3)
      enddo

      ! upper-most model level:
      do i=0,ni+2
        dum1(i,j,nk+1) = cgt1*ua(i,j,nk)+cgt2*ua(i,j,nk-1)+cgt3*ua(i,j,nk-2)
        dum2(i,j,nk+1) = cgt1*va(i,j,nk)+cgt2*va(i,j,nk-1)+cgt3*va(i,j,nk-2)
      enddo

      ! interior:
      do k=2,nk
      r2 = (sigmaf(k)-sigma(k-1))*rds(k)
      r1 = 1.0-r2
      do i=0,ni+2
        dum1(i,j,k) = r1*ua(i,j,k-1)+r2*ua(i,j,k)
        dum2(i,j,k) = r1*va(i,j,k-1)+r2*va(i,j,k)
      enddo
      enddo
    enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1)
    DO k=1,nk
      do j=0,nj+1
      do i=0,ni+1 
        s11(i,j,k)=gz(i,j)*(ua(i+1,j,k)*rgzu(i+1,j)-ua(i,j,k)*rgzu(i,j))*rdx*uh(i) &
                  +( gx(i,j,k+1)*(dum1(i,j,k+1)+dum1(i+1,j,k+1))      &
                    -gx(i,j,k  )*(dum1(i,j,k  )+dum1(i+1,j,k  ))      &
                   )*0.5*rdsf(k)
        s11(i,j,k)=s11(i,j,k)*rho(i,j,k)
        s22(i,j,k)=gz(i,j)*(va(i,j+1,k)*rgzv(i,j+1)-va(i,j,k)*rgzv(i,j))*rdy*vh(j) &
                  +( gy(i,j,k+1)*(dum2(i,j,k+1)+dum2(i,j+1,k+1))      &
                    -gy(i,j,k  )*(dum2(i,j,k  )+dum2(i,j+1,k  ))      &
                   )*0.5*rdsf(k)
        s22(i,j,k)=s22(i,j,k)*rho(i,j,k)
        s33(i,j,k)=(wa(i,j,k+1)-wa(i,j,k))*rdsf(k)*gz(i,j)
        s33(i,j,k)=s33(i,j,k)*rho(i,j,k)
      enddo
      enddo
      do j=1,nj+1 
      do i=1,ni+1
        r1 = 0.25*( ( rho(i-1,j-1,k)*gz(i-1,j-1)   &
                     +rho(i  ,j  ,k)*gz(i  ,j  ) ) &
                   +( rho(i-1,j  ,k)*gz(i-1,j  )   &
                     +rho(i  ,j-1,k)*gz(i  ,j-1) ) )
        s12(i,j,k)=0.5*(                                                         &
                   ( r1*(ua(i,j,k)*rgzu(i,j)-ua(i,j-1,k)*rgzu(i,j-1))*rdy*vf(j)  &
                    +0.5*( (zt-sigmaf(k+1))*(dum1(i,j-1,k+1)+dum1(i,j,k+1))      &
                          -(zt-sigmaf(k  ))*(dum1(i,j-1,k  )+dum1(i,j,k  ))      &
                         )*rdsf(k)*r1*(rgzu(i,j)-rgzu(i,j-1))*rdy*vf(j) )        &
                  +( r1*(va(i,j,k)*rgzv(i,j)-va(i-1,j,k)*rgzv(i-1,j))*rdx*uf(i)  &
                    +0.5*( (zt-sigmaf(k+1))*(dum2(i-1,j,k+1)+dum2(i,j,k+1))      &
                          -(zt-sigmaf(k  ))*(dum2(i-1,j,k  )+dum2(i,j,k  ))      &
                         )*rdsf(k)*r1*(rgzv(i,j)-rgzv(i-1,j))*rdx*uf(i) )    )
      enddo
      enddo       
    ENDDO

    ! now, dum1 stores w at scalar-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=0,nj+1
      do i=0,ni+1
        dum1(i,j,k)=0.5*(wa(i,j,k)+wa(i,j,k+1))
      enddo
      enddo
    ENDDO
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=2,nk
      do j=1,nj
      do i=1,ni+1
        s13(i,j,k)=0.5*(                                                              &
                   (ua(i,j,k)-ua(i,j,k-1))*rds(k)                                     &
                  +(wa(i,j,k)*rgz(i,j)-wa(i-1,j,k)*rgz(i-1,j))*rdx*uf(i)              &
                  +0.5*rds(k)*( (zt-sigma(k  ))*(dum1(i,j,k  )+dum1(i-1,j,k  ))       &
                               -(zt-sigma(k-1))*(dum1(i,j,k-1)+dum1(i-1,j,k-1)) )     &
                             *(rgz(i,j)-rgz(i-1,j))*rdx*uf(i)                         )
        s13(i,j,k)=s13(i,j,k)*0.5*( gz(i-1,j)*rf(i-1,j,k)+gz(i,j)*rf(i,j,k) )
      enddo
      enddo
      do j=1,nj+1   
      do i=1,ni
        s23(i,j,k)=0.5*(                                                              &
                   (va(i,j,k)-va(i,j,k-1))*rds(k)                                     &
                  +(wa(i,j,k)*rgz(i,j)-wa(i,j-1,k)*rgz(i,j-1))*rdy*vf(j)              &
                  +0.5*rds(k)*( (zt-sigma(k  ))*(dum1(i,j,k  )+dum1(i,j-1,k  ))       &
                               -(zt-sigma(k-1))*(dum1(i,j,k-1)+dum1(i,j-1,k-1)) )     &
                             *(rgz(i,j)-rgz(i,j-1))*rdy*vf(j)                         )
        s23(i,j,k)=s23(i,j,k)*0.5*( gz(i,j-1)*rf(i,j-1,k)+gz(i,j)*rf(i,j,k) )
      enddo
      enddo
    ENDDO

  ENDIF

!  end of calculations for terrain
!-------------------------------------------------------------------------------
!  open boundary conditions:

    IF( wbc.eq.2 .or. ebc.eq.2 .or. sbc.eq.2 .or. nbc.eq.2 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
        !-----
        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=0,nj+1
            s11(0,j,k) = s11(1,j,k)
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=0,nj+1
            s11(ni+1,j,k) = s11(ni,j,k)
          enddo
        ENDIF
        !-----
        !ccccc
        !-----
        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=0,ni+1
            s22(i,0,k) = s22(i,1,k)
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=0,ni+1
            s22(i,nj+1,k) = s22(i,nj,k)
          enddo
        ENDIF
        !-----
        !ccccc
        !-----
        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj+1
            s12(1,j,k) = s12(2,j,k)
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj+1
            s12(ni+1,j,k) = s12(ni,j,k)
          enddo
        ENDIF
        !-----
        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni+1
            s12(i,1,k) = s12(i,2,k)
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni+1
            s12(i,nj+1,k) = s12(i,nj,k)
          enddo
        ENDIF
        !-----
        ! corner points:
        !-----
        IF( sbc.eq.2 .and. ibs.eq.1 .and. &
            wbc.eq.2 .and. ibw.eq.1 )THEN
          s12(1,1,k) = s12(2,2,k)
        ENDIF
        IF( sbc.eq.2 .and. ibs.eq.1 .and. &
            ebc.eq.2 .and. ibe.eq.1 )THEN
          s12(ni+1,1,k) = s12(ni,2,k)
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 .and. &
            wbc.eq.2 .and. ibw.eq.1 )THEN
          s12(1,nj+1,k) = s12(2,nj,k)
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 .and. &
            ebc.eq.2 .and. ibe.eq.1 )THEN
          s12(ni+1,nj+1,k) = s12(ni,nj,k)
        ENDIF
        !-----
      ENDDO
    ENDIF

!----------------------------------------------------------------------
!  if l_h or l_v is zero, set appropriate terms to zero:
!    (just to be sure)

    IF( iturb.eq.3 .and. l_h*lhref1*lhref2.lt.1.0e-12 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        s11(i,j,k) = 0.0
        s12(i,j,k) = 0.0
        s33(i,j,k) = 0.0
        s22(i,j,k) = 0.0
      enddo
      enddo
      enddo
    ENDIF

    IF( iturb.eq.3 .and. l_inf.lt.tsmall )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        s13(i,j,k) = 0.0
        s23(i,j,k) = 0.0
      enddo
      enddo
      enddo
    ENDIF

!--------------------------------------------------------------

    IF( axisymm.eq.1 )THEN
      ! lateral boundary condition:
!$omp parallel do default(shared)   &
!$omp private(k)
      do k=0,nk+1
        s13(1,1,k)=0.0
      enddo
    ENDIF

!----------------------------------------------------------------------
!  calculate deformation:
!  Note:  deformation is defined at w points.

    IF(axisymm.eq.0)THEN
      ! Cartesian domain:

      ! Def = 2.0 * S_ij * S_ij
      !
      !     = 2.0 * (  S11*S11 + S12*S12 + S13*S13 
      !              + S21*S21 + S22*S22 + S23*S23 
      !              + S31*S31 + S32*S32 + S33*S33 )
      !
      !     =   2.0*( S11*S11 + S22*S22 + S33*S33 )
      !       + 4.0*( S12*S12 + S13*S13 + S23*S23 )

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tmp11,tmp22,tmp33,tmp12,tmp13,tmp23,rrf)
      do k=2,nk
      do j=1,nj
      do i=1,ni

        tmp11=( c1(i,j,k)*s11(i,j,k-1)**2 + c2(i,j,k)*s11(i,j,k)**2 )
        tmp22=( c1(i,j,k)*s22(i,j,k-1)**2 + c2(i,j,k)*s22(i,j,k)**2 )
        tmp33=( c1(i,j,k)*s33(i,j,k-1)**2 + c2(i,j,k)*s33(i,j,k)**2 )

        tmp12=0.25*( c1(i,j,k)*( ( s12(i,j  ,k-1)**2 + s12(i+1,j+1,k-1)**2 )     &
                               + ( s12(i,j+1,k-1)**2 + s12(i+1,j  ,k-1)**2 ) )   &
                    +c2(i,j,k)*( ( s12(i,j  ,k  )**2 + s12(i+1,j+1,k  )**2 )     &
                               + ( s12(i,j+1,k  )**2 + s12(i+1,j  ,k  )**2 ) ) )

        tmp13=0.5*( s13(i,j,k)**2 + s13(i+1,j,k)**2 )

        tmp23=0.5*( s23(i,j,k)**2 + s23(i,j+1,k)**2 )

        rrf = 1.0/(rf(i,j,k)**2)

        defv(i,j,k)= 4.0*( tmp13 + tmp23 )*rrf

        defh(i,j,k) = ( 2.0*( ( tmp11 + tmp22 ) + tmp33 ) + 4.0*tmp12 )*rrf

      enddo
      enddo
      enddo

!--------------------------------------------
    ELSE
      ! axisymmetric domain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tmp11,tmp22,tmp33,tmp12,tmp13,tmp23,rrf,r1,r2,r3,r4)
      do k=2,nk
      do j=1,nj
      do i=1,ni

        tmp11=( c1(1,1,k)*(s11(i,j,k-1)**2) + c2(1,1,k)*(s11(i,j,k)**2) )
        tmp33=( c1(1,1,k)*(s33(i,j,k-1)**2) + c2(1,1,k)*(s33(i,j,k)**2) )

        tmp12=0.5*(  c1(1,1,k)*( s12(i,j  ,k-1)**2 + s12(i+1,j  ,k-1)**2 )     &
                   + c2(1,1,k)*( s12(i,j  ,k  )**2 + s12(i+1,j  ,k  )**2 ) )

        tmp13=0.5*( s13(i,j,k)**2 + s13(i+1,j,k)**2 )

        tmp23=      s23(i,j,k)**2

        rrf = 1.0/(rf(i,j,k)**2)

        defv(i,j,k)= 4.0*( tmp13 + tmp23 )*rrf

        defh(i,j,k) = ( 2.0*( tmp11 + tmp33 ) + 4.0*tmp12 )*rrf

      enddo
      enddo
      enddo

    ENDIF  ! endif for axisymm

!--------------------------------------------------------------
!  finished

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end subroutine calcdef
      ! calcdef


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine calcnm(c1,c2,mf,pi0,thv0,th0,cloudvar,nm,t,qt,thv,cloud,   &
                        prs,pp,th,qa)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'goddard.incl'

      logical, dimension(maxq) :: cloudvar
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0,th0
      real, dimension(ib:ie,jb:je,kb:ke+1) :: nm
      real, dimension(ib:ie,jb:je,kb:ke) :: t,qt,thv,cloud,prs
      real, dimension(ib:ie,jb:je,kb:ke) :: pp,th
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa

      integer i,j,k,n
      real pavg,tavg,qtavg,esl,qvs,lhv,cpml,gamma,qiavg,qsavg,qgavg,drdt
      real qlavg,qvl,qvi,fliq,fice,nmtmp

!----------------------------------------------------------------------
!  Dry nm

    IF(imoist.eq.0)then

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
    do j=1,nj

      do k=2,nk
      do i=1,ni
        nm(i,j,k)=alog( (th0(i,j,k)+th(i,j,k))/(th0(i,j,k-1)+th(i,j,k-1)) ) &
                    *g*rdz*mf(i,j,k)
      enddo
      enddo
      do i=1,ni
        nm(i,j,   1)=0.0
        nm(i,j,nk+1)=0.0
      enddo

    enddo

!-----------------------------------------------------------------------
!  Moist nm

    ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,pavg,tavg,qtavg,esl,qvs,lhv,cpml,drdt,gamma,nmtmp)
    DO j=1,nj

      do k=1,nk
      do i=1,ni
        t(i,j,k)=(th0(i,j,k)+th(i,j,k))*(pi0(i,j,k)+pp(i,j,k))
      enddo
      enddo

      do k=1,nk
      do i=1,ni
        qt(i,j,k)=0.0
      enddo
      enddo

      DO n=1,numq
        IF( (n.eq.nqv) .or.                                 &
            (n.ge.nql1.and.n.le.nql2) .or.                  &
            (n.ge.nqs1.and.n.le.nqs2.and.iice.eq.1) )THEN
          do k=1,nk
          do i=1,ni
            qt(i,j,k)=qt(i,j,k)+qa(i,j,k,n)
          enddo
          enddo
        ENDIF
      ENDDO

      do k=1,nk
      do i=1,ni
        thv(i,j,k)=(th0(i,j,k)+th(i,j,k))*(1.0+reps*qa(i,j,k,nqv))   &
                                         /(1.0+qt(i,j,k))
      enddo
      enddo

      do k=2,nk
      do i=1,ni
        nm(i,j,k)=g*alog(thv(i,j,k)/thv(i,j,k-1))*rdz*mf(i,j,k)
      enddo
      enddo

      do i=1,ni
        nm(i,j,   1)=0.0
        nm(i,j,nk+1)=0.0
      enddo

      do k=1,nk
      do i=1,ni
        cloud(i,j,k)=0.0
      enddo
      enddo
      do n=1,numq
        if(cloudvar(n))then
          do k=1,nk
          do i=1,ni
            cloud(i,j,k)=cloud(i,j,k)+qa(i,j,k,n)
          enddo
          enddo
        endif
      enddo

      do k=2,nk
      do i=1,ni
        IF( (cloud(i,j,k).ge.clwsat) .or. (cloud(i,j,k-1).ge.clwsat) )THEN
          pavg = c1(i,j,k)*prs(i,j,k-1)+c2(i,j,k)*prs(i,j,k)
          tavg =   c1(i,j,k)*t(i,j,k-1)+  c2(i,j,k)*t(i,j,k)
          qtavg=  c1(i,j,k)*qt(i,j,k-1)+ c2(i,j,k)*qt(i,j,k)
          esl = 611.2*exp( 17.67 * ( tavg - 273.15 ) / ( tavg - 29.65 ) )
          qvs = eps*esl/(pavg-esl)
          lhv=lv1-lv2*tavg
          cpml=cp+cpv*qvs+cpl*(qtavg-qvs)

          drdt=17.67*(273.15-29.65)*qvs/((tavg-29.65)**2)
          gamma=g*(1.0+qtavg)*(1.0+lhv*qvs/(rd*tavg))/(cpml+lhv*drdt)
          nmtmp=g*( ( alog(t(i,j,k)/t(i,j,k-1))*rdz*mf(i,j,k)      &
                            +gamma/tavg )*(1.0+tavg*drdt/(eps+qvs))   &
                         -alog((1.0+qt(i,j,k))/(1.0+qt(i,j,k-1)))*rdz*mf(i,j,k) )
        IF( (cloud(i,j,k).ge.clwsat) .and. (cloud(i,j,k-1).ge.clwsat) )THEN
          nm(i,j,k)=nmtmp
        ELSE
          nm(i,j,k)=0.5*(nm(i,j,k)+nmtmp)
        ENDIF
        ENDIF
      enddo
      enddo

    ENDDO

    ENDIF    ! endif for imoist

!----------------------------------------------------------------------

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end subroutine calcnm


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbs(iflux,dt,dosfcflx,xh,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,sflux,  &
                       rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                       turbx,turby,turbz,dumx,dumy,dumz,rho,rr,rf,s,sten,khh,khv,dum7,dum8)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      integer iflux
      real :: dt
      logical, intent(in) :: dosfcflx
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf,arf1,arf2,uf
      real, dimension(jb:je) :: vh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je) :: sflux
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dumx,dumy,dumz,rho,rr,rf,s,sten
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: khh,khv
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: dum7,dum8

      integer :: i,j,k,k1,k2
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2,cfa,cfb,cfc,cfd

!---------------------------------------------------------------

  IF(.not.terrain_flag)THEN

    IF(axisymm.eq.0)THEN
      ! Cartesian without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        !  x-direction
        do j=1,nj
        do i=1,ni+1
          dumx(i,j,k)= -0.125*( rho(i,j,k)+rho(i-1,j,k) )           &
                             *(  (khh(i,j,k  )+ khh(i-1,j,k  ))     &
                                +(khh(i,j,k+1)+ khh(i-1,j,k+1)) )   &
                             *(    s(i,j,k)-   s(i-1,j,k) )*rdx*uf(i)
        enddo
        enddo
        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj
            dumx(1,j,k) = dumx(2,j,k)
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj
            dumx(ni+1,j,k) = dumx(ni,j,k)
          enddo
        ENDIF
        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=-(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
        enddo
        enddo

        !  y-direction
        do j=1,nj+1
        do i=1,ni
          dumy(i,j,k)= -0.125*( rho(i,j,k)+rho(i,j-1,k) )           &
                             *(  (khh(i,j,k  )+ khh(i,j-1,k  ))     &
                                +(khh(i,j,k+1)+ khh(i,j-1,k+1)) )   &
                            *(    s(i,j,k)-   s(i,j-1,k) )*rdy*vf(j)
        enddo
        enddo
        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni
            dumy(i,1,k) = dumy(i,2,k)
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni
            dumy(i,nj+1,k) = dumy(i,nj,k)
          enddo
        ENDIF
        do j=1,nj
        do i=1,ni
          turby(i,j,k)=-(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
        enddo
        enddo

      enddo

    ELSE
      ! axisymmetric:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj

        do i=1,ni+1
          dumx(i,j,k)= -0.125*( rho(i,j,k)+rho(i-1,j,k) )           &
                             *(  (khh(i,j,k  )+ khh(i-1,j,k  ))     &
                                +(khh(i,j,k+1)+ khh(i-1,j,k+1)) )   &
                             *(    s(i,j,k)-   s(i-1,j,k) )*rdx*uf(i)
        enddo
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          dumx(ni+1,j,k) = arh1(ni)*dumx(ni,j,k)/arh2(ni)
        ENDIF
        !-----
        if(wbc.eq.3.or.wbc.eq.4)then
          ! assume zero flux:
          dumx(1,j,k) = 0.0
        endif
        if(ebc.eq.3.or.ebc.eq.4)then
        ! assume zero flux:
          dumx(ni+1,j,k) = 0.0
        endif
        !-----
        do i=1,ni
          turbx(i,j,k)=-(arh2(i)*dumx(i+1,j,k)-arh1(i)*dumx(i,j,k))*rdx*uh(i)
        enddo
        do i=1,ni
          turby(i,j,k)=0.0
        enddo

      enddo
      enddo

    ENDIF   ! endif for axisymm check

!---------------------------------------------------------------

  ELSE
      ! Cartesian with terrain:

      ! use turbz as a temporary array for s at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do j=0,nj+1

        ! lowest model level:
        do i=0,ni+1
          turbz(i,j,1) = cgs1*s(i,j,1)+cgs2*s(i,j,2)+cgs3*s(i,j,3)
        enddo

        ! upper-most model level:
        do i=0,ni+1
          turbz(i,j,nk+1) = cgt1*s(i,j,nk)+cgt2*s(i,j,nk-1)+cgt3*s(i,j,nk-2)
        enddo

        ! interior:
        do k=2,nk
        r2 = (sigmaf(k)-sigma(k-1))*rds(k)
        r1 = 1.0-r2
        do i=0,ni+1
          turbz(i,j,k) = r1*s(i,j,k-1)+r2*s(i,j,k)
        enddo
        enddo

      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        ! x-flux
        do j=1,nj
        do i=1,ni+1
          dumx(i,j,k)= -0.125*( gz(i,j)*rho(i,j,k)+gz(i-1,j)*rho(i-1,j,k) )  &
                             *(  (khh(i,j,k  )+ khh(i-1,j,k  ))     &
                                +(khh(i,j,k+1)+ khh(i-1,j,k+1)) )*( &
                  (s(i,j,k)*rgz(i,j)-s(i-1,j,k)*rgz(i-1,j))         &
                   *rdx*uf(i)                                       &
              +0.5*( gxu(i,j,k+1)*(turbz(i,j,k+1)+turbz(i-1,j,k+1)) &
                    -gxu(i,j,k  )*(turbz(i,j,k  )+turbz(i-1,j,k  )) &
                   )*rdsf(k)*rgzu(i,j) )
        enddo
        enddo

        ! y-flux
        do j=1,nj+1
        do i=1,ni
          dumy(i,j,k)= -0.125*( gz(i,j)*rho(i,j,k)+gz(i,j-1)*rho(i,j-1,k) )  &
                             *(  (khh(i,j,k  )+ khh(i,j-1,k  ))     &
                                +(khh(i,j,k+1)+ khh(i,j-1,k+1)) )*( &
                  (s(i,j,k)*rgz(i,j)-s(i,j-1,k)*rgz(i,j-1))         &
                   *rdy*vf(j)                                       &
              +0.5*( gyv(i,j,k+1)*(turbz(i,j,k+1)+turbz(i,j-1,k+1)) &
                    -gyv(i,j,k  )*(turbz(i,j,k  )+turbz(i,j-1,k  )) &
                   )*rdsf(k)*rgzv(i,j) )
        enddo
        enddo

      enddo

      ! use turbz,dumz as temporary arrays for fluxes at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do j=1,nj+1
        ! lowest model level:
        do i=1,ni+1
          turbz(i,j,1) = cgs1*dumx(i,j,1)+cgs2*dumx(i,j,2)+cgs3*dumx(i,j,3)
           dumz(i,j,1) = cgs1*dumy(i,j,1)+cgs2*dumy(i,j,2)+cgs3*dumy(i,j,3)
        enddo

        ! upper-most model level:
        do i=1,ni+1
          turbz(i,j,nk+1) = cgt1*dumx(i,j,nk)+cgt2*dumx(i,j,nk-1)+cgt3*dumx(i,j,nk-2)
           dumz(i,j,nk+1) = cgt1*dumy(i,j,nk)+cgt2*dumy(i,j,nk-1)+cgt3*dumy(i,j,nk-2)
        enddo

        ! interior:
        do k=2,nk
        r2 = (sigmaf(k)-sigma(k-1))*rds(k)
        r1 = 1.0-r2
        do i=1,ni+1
          turbz(i,j,k) = r1*dumx(i,j,k-1)+r2*dumx(i,j,k)
           dumz(i,j,k) = r1*dumy(i,j,k-1)+r2*dumy(i,j,k)
        enddo
        enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        ! x-tendency
        do j=1,nj
        do i=1,ni
          turbx(i,j,k) = -gz(i,j)*( dumx(i+1,j,k)*rgzu(i+1,j)             &
                                   -dumx(i  ,j,k)*rgzu(i  ,j) )*rdx*uh(i) &
                -( ( gx(i,j,k+1)*(turbz(i,j,k+1)+turbz(i+1,j,k+1))        &
                    -gx(i,j,k  )*(turbz(i,j,k  )+turbz(i+1,j,k  )) )      &
                 )*0.5*rdsf(k)
        enddo
        enddo

        ! y-tendency
        do j=1,nj
        do i=1,ni
          turby(i,j,k) = -gz(i,j)*( dumy(i,j+1,k)*rgzv(i,j+1)             &
                                   -dumy(i,j  ,k)*rgzv(i,j  ) )*rdy*vh(j) &
                -( ( gy(i,j,k+1)*( dumz(i,j,k+1)+ dumz(i,j+1,k+1))        &
                    -gy(i,j,k  )*( dumz(i,j,k  )+ dumz(i,j+1,k  )) )      &
                 )*0.5*rdsf(k)
        enddo
        enddo

      enddo

    IF( wbc.eq.2 .or. ebc.eq.2 .or. sbc.eq.2 .or. nbc.eq.2 )THEN
      !  open boundary conditions:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk

        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj
            turbx(1,j,k) = 0.0
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj
            turbx(ni,j,k) = 0.0
          enddo
        ENDIF

        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni
            turby(i,1,k) = 0.0
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni
            turby(i,nj,k) = 0.0
          enddo
        ENDIF

      ENDDO
    ENDIF

  ENDIF  ! endif for terrain check

!---------------------------------------------------------------------
!  z-direction

    IF( iturb.eq.3 .and. l_inf.lt.tsmall )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

    ELSE

      rdt = 1.0/dt
      tema = -1.0*dt*vialpha*rdz*rdz
      temb = dt*vibeta*rdz*rdz
      temc = dt*rdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2,cfa,cfb,cfc,cfd,tem)
    DO j=1,nj
      k = 1
      DO i=1,ni
          r1 = 0.0
          r2 = dum8(i,j,k)
          cfa = 0.0
          cfc = tema*r2
          cfb = 1.0 - cfc
          cfd = s(i,j,k) + temb*( -r2*s(i,j,k)+r2*s(i,j,k+1) )
          IF(bcturbs.eq.1)THEN
            dumz(i,j,1)=0.0
          ELSEIF(bcturbs.eq.2)THEN
            dumz(i,j,1) = -khv(i,j,2)*(s(i,j,2)-s(i,j,1))*rdz*mf(i,j,2)*rf(i,j,2)
          ENDIF
          if(iflux.eq.1 .and. dosfcflx)then
            dumz(i,j,1)=sflux(i,j)*rf(i,j,1)
          endif
          cfd = cfd + temc*dumz(i,j,1)*mh(i,j,1)*rr(i,j,1)
        tem = 1.0/cfb
        dumx(i,j,1)=-cfc*tem
        dumy(i,j,1)= cfd*tem
      ENDDO
        !--------
        do k=2,nk-1
        do i=1,ni
          r1 = dum7(i,j,k)
          r2 = dum8(i,j,k)
          cfa = tema*r1
          cfc = tema*r2
          cfb = 1.0 - cfa - cfc
          cfd = s(i,j,k) + temb*(r1*s(i,j,k-1)-(r1+r2)*s(i,j,k)+r2*s(i,j,k+1) )
          tem = 1.0/(cfa*dumx(i,j,k-1)+cfb)
          dumx(i,j,k)=-cfc*tem
          dumy(i,j,k)=(cfd-cfa*dumy(i,j,k-1))*tem
        enddo
        enddo
        !--------
        k = nk
        do i=1,ni
          r1 = dum7(i,j,k)
          r2 = 0.0
          cfa = tema*r1
          cfc = 0.0
          cfb = 1.0 - cfa
          cfd = s(i,j,k) + temb*( r1*s(i,j,k-1)-r1*s(i,j,k) )
          IF(bcturbs.eq.1)THEN
            dumz(i,j,nk+1)=0.0
          ELSEIF(bcturbs.eq.2)THEN
            dumz(i,j,nk+1) = -khv(i,j,nk)*(s(i,j,nk)-s(i,j,nk-1))*rdz*mf(i,j,nk)*rf(i,j,nk)
          ENDIF
          cfd = cfd - temc*dumz(i,j,nk+1)*mh(i,j,nk)*rr(i,j,nk)
          tem = 1.0/(cfa*dumx(i,j,k-1)+cfb)
!!!          dumx(i,j,k)=-cfc*tem
          dumy(i,j,k)=(cfd-cfa*dumy(i,j,k-1))*tem
          dumz(i,j,k)=dumy(i,j,k)
          turbz(i,j,k) = rho(i,j,k)*(dumz(i,j,k)-s(i,j,k))*rdt
        enddo
        !--------

      do k=nk-1,1,-1
      DO i=1,ni
          dumz(i,j,k)=dumx(i,j,k)*dumz(i,j,k+1)+dumy(i,j,k)
          turbz(i,j,k) = rho(i,j,k)*(dumz(i,j,k)-s(i,j,k))*rdt
      ENDDO
      enddo

    ENDDO

  ENDIF

!---------------------------------------------------------------------

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=sten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))*rr(i,j,k)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=sten(i,j,k)+(turbx(i,j,k)+turbz(i,j,k))*rr(i,j,k)
      enddo
      enddo
      enddo

    ENDIF

!---------------------------------------------------------------------

      if(timestats.ge.1) time_tmix=time_tmix+mytime()

      return
      end subroutine turbs


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbt(dt,xh,rxh,uh,xf,uf,vh,vf,mh,mf,rho,rr,rf,          &
                       rds,sigma,gz,rgz,gzu,rgzu,gzv,rgzv,                &
                       turbx,turby,turbz,dumx,dumy,dumz,t,tten,kmh,kmv)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,rxh,uh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: rf
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dumx,dumy,dumz
      real, intent(in), dimension(ibt:iet,jbt:jet,kbt:ket) :: t
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tten
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv

      integer :: i,j,k
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2,rrf
      real :: cfa,cfb,cfc,cfd

!---------------------------------------------------------------

    IF(.not.terrain_flag)THEN
      ! Cartesian without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        !  x-direction
        do j=1,nj
        do i=1,ni+1
          ! note:  K is multiplied by 2:
          dumx(i,j,k)= -0.25*( rf(i,j,k)+rf(i-1,j,k) )   &
                       *2.0*( kmh(i,j,k)+kmh(i-1,j,k) )   &
                           *(   t(i,j,k)-  t(i-1,j,k) )*rdx*uf(i)
        enddo
        enddo
        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj
            dumx(1,j,k) = dumx(2,j,k)
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj
            dumx(ni+1,j,k) = dumx(ni,j,k)
          enddo
        ENDIF
        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=-(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
        enddo
        enddo

        !  y-direction
        do j=1,nj+1
        do i=1,ni
          ! note:  K is multiplied by 2:
          dumy(i,j,k)= -0.25*( rf(i,j,k)+rf(i,j-1,k) )   &
                       *2.0*( kmh(i,j,k)+kmh(i,j-1,k) )   &
                           *(   t(i,j,k)-  t(i,j-1,k) )*rdy*vf(j)
        enddo
        enddo
        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni
            dumy(i,1,k) = dumy(i,2,k)
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni
            dumy(i,nj+1,k) = dumy(i,nj,k)
          enddo
        ENDIF
        do j=1,nj
        do i=1,ni
          turby(i,j,k)=-(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
        enddo
        enddo

      enddo

!---------------------------------------------------------------
!  Cartesian with terrain:

    ELSE

      ! turbz stores t at s-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          turbz(i,j,k) = 0.5*(t(i,j,k)+t(i,j,k+1))
        enddo
        enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        ! x-flux:
        do j=1,nj
        do i=1,ni+1
          ! note:  K is multiplied by 2:
          dumx(i,j,k)= -0.25*( gz(i,j)*rf(i,j,k)+gz(i-1,j)*rf(i-1,j,k) )                 &
                       *2.0*( kmh(i,j,k)+kmh(i-1,j,k) )*(                                &
                            (t(i,j,k)*rgz(i,j)-t(i-1,j,k)*rgz(i-1,j))*rdx*uf(i)          &
                     +0.5*( (zt-sigma(k  ))*(turbz(i-1,j,k  )+turbz(i,j,k  ))            &
                           -(zt-sigma(k-1))*(turbz(i-1,j,k-1)+turbz(i,j,k-1))            &
                          )*rds(k)*(rgz(i,j)-rgz(i-1,j))*rdx*uf(i)                       &
                                                        )
        enddo
        enddo

        ! y-flux:
        do j=1,nj+1
        do i=1,ni
          ! note:  K is multiplied by 2:
          dumy(i,j,k)= -0.25*( gz(i,j)*rf(i,j,k)+gz(i,j-1)*rf(i,j-1,k) )                 &
                       *2.0*( kmh(i,j,k)+kmh(i,j-1,k) )*(                                &
                            (t(i,j,k)*rgz(i,j)-t(i,j-1,k)*rgz(i,j-1))*rdy*vf(j)          &
                     +0.5*( (zt-sigma(k  ))*(turbz(i,j-1,k  )+turbz(i,j,k  ))            &
                           -(zt-sigma(k-1))*(turbz(i,j-1,k-1)+turbz(i,j,k-1))            &
                          )*rds(k)*(rgz(i,j)-rgz(i,j-1))*rdy*vf(j)                       &
                                                        )
        enddo
        enddo

      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do j=1,nj+1
        do i=1,ni+1
          dumx(i,j,   1)=0.0
          dumx(i,j,nk+1)=0.0
          dumy(i,j,   1)=0.0
          dumy(i,j,nk+1)=0.0
        enddo
        enddo

      ! turbz stores dumx at s-pts:
      !  dumz stores dumy at s-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          turbz(i,j,k)=0.5*(dumx(i,j,k)+dumx(i,j,k+1))
           dumz(i,j,k)=0.5*(dumy(i,j,k)+dumy(i,j,k+1))
        enddo
        enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        ! x-tendency:
        do j=1,nj
        do i=1,ni
          turbx(i,j,k) = -(dumx(i+1,j,k)*rgzu(i+1,j)-dumx(i,j,k)*rgzu(i,j))*gz(i,j)*rdx*uh(i) &
                         -0.5*( (zt-sigma(k  ))*(turbz(i,j,k  )+turbz(i+1,j,k  ))             &
                               -(zt-sigma(k-1))*(turbz(i,j,k-1)+turbz(i+1,j,k-1))             &
                              )*rds(k)*(rgzu(i+1,j)-rgzu(i,j))*gz(i,j)*rdx*uh(i)
        enddo
        enddo

        ! y-tendency:
        do j=1,nj
        do i=1,ni
          turby(i,j,k) = -(dumy(i,j+1,k)*rgzv(i,j+1)-dumy(i,j,k)*rgzv(i,j))*gz(i,j)*rdy*vh(j) &
                         -0.5*( (zt-sigma(k  ))*( dumz(i,j,k  )+ dumz(i,j+1,k  ))             &
                               -(zt-sigma(k-1))*( dumz(i,j,k-1)+ dumz(i,j+1,k-1))             &
                              )*rds(k)*(rgzv(i,j+1)-rgzv(i,j))*gz(i,j)*rdy*vh(j)
        enddo
        enddo

      enddo

    IF( wbc.eq.2 .or. ebc.eq.2 .or. sbc.eq.2 .or. nbc.eq.2 )THEN
      !  open boundary conditions:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=2,nk

        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj
            turbx(1,j,k) = 0.0
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj
            turbx(ni,j,k) = 0.0
          enddo
        ENDIF

        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni
            turby(i,1,k) = 0.0
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni
            turby(i,nj,k) = 0.0
          enddo
        ENDIF

      ENDDO
    ENDIF

    ENDIF  ! endif for terrain check

!---------------------------------------------------------------------
!  z-direction

      rdt = 1.0/dt
      tema = -1.0*dt*vialpha*rdz*rdz
      temb =      dt*vibeta*rdz*rdz
      temc = dt*rdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2,cfa,cfb,cfc,cfd,tem,rrf)
      do j=1,nj

        k=2
        do i=1,ni
          rrf = mf(i,j,k)/rf(i,j,k)
          r2 = (kmh(i,j,k  )+kmh(i,j,k+1))*mh(i,j,k  )*rho(i,j,k  )*rrf
          cfc = tema*r2
          cfb = 1.0 - cfc
          cfd = t(i,j,k) + temb*( r2*t(i,j,k+1)-r2*t(i,j,k) )
          tem = -(kmv(i,j,k-1)+kmv(i,j,k))*(t(i,j,k)-t(i,j,k-1))*rdz*mh(i,j,k-1)*rho(i,j,k-1)
          cfd = cfd + temc*tem*rrf
          tem = 1.0/cfb
          dumx(i,j,k) = -cfc*tem
          dumy(i,j,k) =  cfd*tem
        enddo

        do k=3,(nk-1)
        do i=1,ni
          rrf = mf(i,j,k)/rf(i,j,k)
          r1 = (kmh(i,j,k-1)+kmh(i,j,k  ))*mh(i,j,k-1)*rho(i,j,k-1)*rrf
          r2 = (kmh(i,j,k  )+kmh(i,j,k+1))*mh(i,j,k  )*rho(i,j,k  )*rrf
          cfa = tema*r1
          cfc = tema*r2
          cfb = 1.0 - cfa - cfc
          cfd = t(i,j,k) + temb*(r2*t(i,j,k+1)-(r1+r2)*t(i,j,k)+r1*t(i,j,k-1))
          tem = 1.0/(cfa*dumx(i,j,k-1)+cfb)
          dumx(i,j,k) = -cfc*tem
          dumy(i,j,k) = (cfd-cfa*dumy(i,j,k-1))*tem
        enddo
        enddo

        k = nk
        do i=1,ni
          rrf = mf(i,j,k)/rf(i,j,k)
          r1 = (kmh(i,j,k-1)+kmh(i,j,k  ))*mh(i,j,k-1)*rho(i,j,k-1)*rrf
          cfa = tema*r1
          cfb = 1.0 - cfa
          cfd = t(i,j,k) + temb*( -r1*t(i,j,k)+r1*t(i,j,k-1) )
          tem = -(kmv(i,j,k)+kmv(i,j,k+1))*(t(i,j,k+1)-t(i,j,k))*rdz*mh(i,j,k)*rho(i,j,k)
          cfd = cfd - temc*tem*rrf
          tem = 1.0/(cfa*dumx(i,j,k-1)+cfb)
          dumy(i,j,k) = (cfd-cfa*dumy(i,j,k-1))*tem
          !---
          dumz(i,j,k) = dumy(i,j,k)
          turbz(i,j,k) = rf(i,j,k)*(dumz(i,j,k)-t(i,j,k))*rdt
        enddo

        do k=(nk-1),2,-1
        do i=1,ni
          dumz(i,j,k) = dumx(i,j,k)*dumz(i,j,k+1)+dumy(i,j,k)
          turbz(i,j,k) = rf(i,j,k)*(dumz(i,j,k)-t(i,j,k))*rdt
        enddo
        enddo

      enddo

!---------------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        tten(i,j,k)=tten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))/rf(i,j,k)
      enddo
      enddo
      enddo

!---------------------------------------------------------------------

      if(timestats.ge.1) time_tmix=time_tmix+mytime()

      return
      end subroutine turbt


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbu(dt,xh,ruh,xf,rxf,arf1,arf2,uf,vh,mh,mf,rmf,rho,rf,  &
                       zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gxu,     &
                       turbx,turby,turbz,dum1,dum2,dum3,u,uten,w,t11,t12,t13,t22,kmv)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,gzv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gxu
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dum1,dum2,dum3
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: uten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmv

      integer :: i,j,k,ip
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2,rru0
      real :: cfa,cfb,cfc,cfd

!---------------------------------------------------------------

  IF(.not.terrain_flag)THEN

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        !  x-direction
        do j=1,nj
        do i=1,ni+1
          turbx(i,j,k)=(t11(i,j,k)-t11(i-1,j,k))*rdx*uf(i)
        enddo
        enddo

        !  y-direction
        do j=1,nj
        do i=1,ni+1
          turby(i,j,k)=(t12(i,j+1,k)-t12(i,j,k))*rdy*vh(j)
        enddo
        enddo

      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(j,k)
      do k=1,nk

        do j=1,nj
        turbx(1,j,k)=0.0
        do i=2,ni+1
          turbx(i,j,k) = ( arf2(i)*arf2(i)*t11(i,j,k) - arf1(i)*arf1(i)*t11(i-1,j,k) )*rdx*uf(i)
        enddo
        IF(ebc.eq.3.or.ebc.eq.4)THEN
          turbx(ni+1,j,k)=0.0
        ENDIF
        enddo

      enddo

    ENDIF

!---------------------------------------------------------------
!  Terrain:

  ELSE

      ! dum1 stores t11 at w-pts:
      ! dum2 stores t12 at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do j=1,nj+1

          ! lowest model level:
          do i=0,ni+1
            dum1(i,j,1) = cgs1*t11(i,j,1)+cgs2*t11(i,j,2)+cgs3*t11(i,j,3)
            dum2(i,j,1) = cgs1*t12(i,j,1)+cgs2*t12(i,j,2)+cgs3*t12(i,j,3)
          enddo

          ! upper-most model level:
          do i=0,ni+1
            dum1(i,j,nk+1) = cgt1*t11(i,j,nk)+cgt2*t11(i,j,nk-1)+cgt3*t11(i,j,nk-2)
            dum2(i,j,nk+1) = cgt1*t12(i,j,nk)+cgt2*t12(i,j,nk-1)+cgt3*t12(i,j,nk-2)
          enddo

          ! interior:
          do k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do i=0,ni+1
            dum1(i,j,k) = r1*t11(i,j,k-1)+r2*t11(i,j,k)
            dum2(i,j,k) = r1*t12(i,j,k-1)+r2*t12(i,j,k)
          enddo
          enddo

      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do k=1,nk

        !  x-direction
        do j=1,nj
        do i=1,ni+1
          turbx(i,j,k)=gzu(i,j)*(t11(i,j,k)*rgz(i,j)-t11(i-1,j,k)*rgz(i-1,j))*rdx*uf(i)  &
                      +0.5*( gxu(i,j,k+1)*(dum1(i-1,j,k+1)+dum1(i,j,k+1))                &
                            -gxu(i,j,k  )*(dum1(i-1,j,k  )+dum1(i,j,k  )) )*rdsf(k)
        enddo
        enddo

        !  y-direction
        do j=1,nj
        do i=1,ni+1
          r1 = 0.25*((rgz(i-1,j-1)+rgz(i,j))+(rgz(i-1,j)+rgz(i,j-1)))
          r2 = 0.25*((rgz(i-1,j+1)+rgz(i,j))+(rgz(i-1,j)+rgz(i,j+1)))
          turby(i,j,k)=gzu(i,j)*(t12(i,j+1,k)*r2-t12(i,j,k)*r1)*rdy*vh(j)      &
                      +0.5*( (zt-sigmaf(k+1))*(dum2(i,j,k+1)+dum2(i,j+1,k+1))  &
                            -(zt-sigmaf(k  ))*(dum2(i,j,k  )+dum2(i,j+1,k  ))  &
                           )*gzu(i,j)*(r2-r1)*rdy*vh(j)*rdsf(k)
        enddo
        enddo

      enddo

  ENDIF  ! endif for terrain check

!-----------------------------------------------------------------
!  open boundary conditions:

    IF( wbc.eq.2 .or. ebc.eq.2 .or. sbc.eq.2 .or. nbc.eq.2 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk

        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj
            turbx(1,j,k) = 0.0
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj
            turbx(ni+1,j,k) = 0.0
          enddo
        ENDIF

        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni+1
            turby(i,1,k) = 0.0
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni+1
            turby(i,nj,k) = 0.0
          enddo
        ENDIF

      ENDDO
    ENDIF

!-----------------------------------------------------------------
!  z-direction

    IF( iturb.eq.3 .and. l_inf.lt.tsmall )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

    ELSE

    IF(axisymm.eq.0)THEN
      rdt = 0.5/dt
      tema = -0.0625*dt*vialpha*rdz*rdz
      temb =  0.0625*dt*vibeta*rdz*rdz
      temc =  0.5*dt*rdz
    ELSE
      rdt = 0.5/dt
      tema = -0.25*dt*vialpha*rdz*rdz
      temb =  0.25*dt*vibeta*rdz*rdz
      temc =  dt*rdz
    ENDIF

      ip = 0
      if( axisymm.eq.1 ) ip = 1

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2,cfa,cfb,cfc,cfd,tem,rru0)
    DO j=1,nj

      check_grid:  IF(axisymm.eq.0)THEN
        ! Cartesian grid:

        !--------
        k = 1
        do i=1+ip,ni+1
          rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
          tem = (mh(i-1,j,k)+mh(i,j,k))*rru0
          r1 = 0.0
          r2 = (kmv(i-1,j,k+1)+kmv(i,j,k+1))*(mf(i-1,j,k+1)+mf(i,j,k+1))   &
              *(rf(i-1,j,k+1)+rf(i,j,k+1))*tem
          cfa = 0.0
          cfc = tema*r2
          cfb = 1.0 - cfc
          cfd = u(i,j,k) + temb*( r2*u(i,j,k+1)-r2*u(i,j,k) )
          tem = temc*t13(i,j,1)*(mh(i-1,j,1)+mh(i,j,1))*rru0
          cfd = cfd - tem
          tem = 1.0/cfb
          dum1(i,j,1)=-cfc*tem
          dum2(i,j,1)= cfd*tem
        enddo
        !--------
        do k=2,nk-1
        do i=1+ip,ni+1
          rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
          tem = (mh(i-1,j,k)+mh(i,j,k))*rru0
          r1 = (kmv(i-1,j,k  )+kmv(i,j,k  ))*(mf(i-1,j,k  )+mf(i,j,k  ))   &
              *(rf(i-1,j,k  )+rf(i,j,k  ))*tem
          r2 = (kmv(i-1,j,k+1)+kmv(i,j,k+1))*(mf(i-1,j,k+1)+mf(i,j,k+1))   &
              *(rf(i-1,j,k+1)+rf(i,j,k+1))*tem
          cfa = tema*r1
          cfc = tema*r2
          cfb = 1.0 - cfa - cfc
          cfd = u(i,j,k) + temb*( r2*u(i,j,k+1)-(r1+r2)*u(i,j,k)+r1*u(i,j,k-1) )
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum1(i,j,k)=-cfc*tem
          dum2(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
        enddo
        enddo
        !--------
        k = nk
        do i=1+ip,ni+1
          rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
          tem = (mh(i-1,j,k)+mh(i,j,k))*rru0
          r1 = (kmv(i-1,j,k  )+kmv(i,j,k  ))*(mf(i-1,j,k  )+mf(i,j,k  ))   &
              *(rf(i-1,j,k  )+rf(i,j,k  ))*tem
          r2 = 0.0
          cfa = tema*r1
          cfc = 0.0
          cfb = 1.0 - cfa
          cfd = u(i,j,k) + temb*( -r1*u(i,j,k)+r1*u(i,j,k-1) )
          tem = temc*t13(i,j,nk+1)*(mh(i-1,j,nk)+mh(i,j,nk))*rru0
          cfd = cfd + tem
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
!!!          dum1(i,j,k)=-cfc*tem
          dum2(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
          dum3(i,j,k)=dum2(i,j,k)
          turbz(i,j,k) = (rho(i-1,j,k)+rho(i,j,k))*(dum3(i,j,k)-u(i,j,k))*rdt
        enddo
        !--------

        do k=nk-1,1,-1
        do i=1+ip,ni+1
          dum3(i,j,k)=dum1(i,j,k)*dum3(i,j,k+1)+dum2(i,j,k)
          turbz(i,j,k) = (rho(i-1,j,k)+rho(i,j,k))*(dum3(i,j,k)-u(i,j,k))*rdt
        enddo
        enddo

      !------------------------------------------------------------
      !------------------------------------------------------------
      !------------------------------------------------------------

      ELSEIF(axisymm.eq.1)THEN
        ! axisymmetric grid:

        !--------
        k = 1
        do i=1+ip,ni+1
          rru0 = 1.0/(0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k)))
          tem = mh(1,1,k)*rru0
          r1 = 0.0
          r2 = (kmv(i-1,j,k+1)+kmv(i,j,k+1))*mf(1,1,k+1)   &
              *(arf1(i)*rf(i-1,j,k+1)+arf2(i)*rf(i,j,k+1))*tem
          cfa = 0.0
          cfc = tema*r2
          cfb = 1.0 - cfc
          cfd = u(i,j,k) + temb*( r2*u(i,j,k+1)-r2*u(i,j,k) )
          tem = temc*t13(i,j,1)*mh(1,1,1)*rru0
          cfd = cfd - tem
          tem = 1.0/cfb
          dum1(i,j,1)=-cfc*tem
          dum2(i,j,1)= cfd*tem
        enddo
        !--------
        do k=2,nk-1
        do i=1+ip,ni+1
          rru0 = 1.0/(0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k)))
          tem = mh(1,1,k)*rru0
          r1 = (kmv(i-1,j,k  )+kmv(i,j,k  ))*mf(1,1,k  )   &
              *(arf1(i)*rf(i-1,j,k  )+arf2(i)*rf(i,j,k  ))*tem
          r2 = (kmv(i-1,j,k+1)+kmv(i,j,k+1))*mf(1,1,k+1)   &
              *(arf1(i)*rf(i-1,j,k+1)+arf2(i)*rf(i,j,k+1))*tem
          cfa = tema*r1
          cfc = tema*r2
          cfb = 1.0 - cfa - cfc
          cfd = u(i,j,k) + temb*( r2*u(i,j,k+1)-(r1+r2)*u(i,j,k)+r1*u(i,j,k-1) )
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum1(i,j,k)=-cfc*tem
          dum2(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
        enddo
        enddo
        !--------
        k = nk
        do i=1+ip,ni+1
          rru0 = 1.0/(0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k)))
          tem = mh(1,1,k)*rru0
          r1 = (kmv(i-1,j,k  )+kmv(i,j,k  ))*mf(1,1,k  )   &
              *(arf1(i)*rf(i-1,j,k  )+arf2(i)*rf(i,j,k  ))*tem
          r2 = 0.0
          cfa = tema*r1
          cfc = 0.0
          cfb = 1.0 - cfa
          cfd = u(i,j,k) + temb*( -r1*u(i,j,k)+r1*u(i,j,k-1) )
          tem = temc*t13(i,j,nk+1)*mh(1,1,nk)*rru0
          cfd = cfd + tem
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
!!!          dum1(i,j,k)=-cfc*tem
          dum2(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
          dum3(i,j,k)=dum2(i,j,k)
          turbz(i,j,k) = (arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k))*(dum3(i,j,k)-u(i,j,k))*rdt
        enddo
        !--------

        do k=nk-1,1,-1
        do i=1+ip,ni+1
          dum3(i,j,k)=dum1(i,j,k)*dum3(i,j,k+1)+dum2(i,j,k)
          turbz(i,j,k) = (arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k))*(dum3(i,j,k)-u(i,j,k))*rdt
        enddo
        enddo

      ENDIF  check_grid

    ENDDO

  IF( terrain_flag )THEN
    ! dum1 stores w at scalar-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=0,nj+1
      do i=0,ni+1
        dum1(i,j,k)=0.5*(w(i,j,k)+w(i,j,k+1))
      enddo
      enddo
    ENDDO
  ENDIF

      ! explicit piece ... dwdx term

    DO j=1,nj
      do i=1+ip,ni+1
        dum2(i,j,1) = 0.0
        dum2(i,j,nk+1) = 0.0
      enddo
      IF(.not.terrain_flag)THEN
        IF( axisymm.eq.0 )THEN
          do k=2,nk
          do i=1+ip,ni+1
            dum2(i,j,k)=(w(i,j,k)-w(i-1,j,k))*rdx*uf(i)
            dum2(i,j,k)=dum2(i,j,k)*0.25*( kmv(i-1,j,k)+kmv(i,j,k) )                 &
                                        *( rf(i-1,j,k)+rf(i,j,k) )
          enddo
          enddo
        ELSEIF( axisymm.eq.1 )THEN
          do k=2,nk
          do i=1+ip,ni+1
            dum2(i,j,k)=(w(i,j,k)-w(i-1,j,k))*rdx*uf(i)
            dum2(i,j,k)=dum2(i,j,k)*0.25*( arf1(i)*kmv(i-1,j,k)+arf2(i)*kmv(i,j,k) )  &
                                        *(arf1(i)*rf(i-1,j,k)+arf2(i)*rf(i,j,k))
          enddo
          enddo
        ENDIF
      ELSE
        do k=2,nk
        do i=1+ip,ni+1
          dum2(i,j,k)=(w(i,j,k)*rgz(i,j)-w(i-1,j,k)*rgz(i-1,j))*rdx*uf(i)          &
                  +0.5*rds(k)*( (zt-sigma(k  ))*(dum1(i,j,k  )+dum1(i-1,j,k  ))    &
                               -(zt-sigma(k-1))*(dum1(i,j,k-1)+dum1(i-1,j,k-1)) )  &
                             *(rgz(i,j)-rgz(i-1,j))*rdx*uf(i)
          dum2(i,j,k)=dum2(i,j,k)*0.25*( kmv(i-1,j,k)+kmv(i,j,k) )                 &
                                      *( gz(i-1,j)*rf(i-1,j,k)+gz(i,j)*rf(i,j,k) )
        enddo
        enddo
      ENDIF
      do k=1,nk
      do i=1+ip,ni+1
        turbz(i,j,k)=turbz(i,j,k)+(dum2(i,j,k+1)-dum2(i,j,k))*rdz*0.5*(mh(i-1,j,k)+mh(i,j,k))
      enddo
      enddo
    ENDDO


      IF(axisymm.eq.1)THEN
        DO k=1,nk
          turbz(1,1,k) = 0.0
        ENDDO
        IF( ebc.eq.3 .or. ebc.eq.4 )THEN
          do k=1,nk
            turbz(ni+1,1,k)=0.0
          enddo
        ENDIF
      ENDIF

  ENDIF

!-----------------------------------------------------------------

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,rru0)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
        uten(i,j,k)=uten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))*rru0
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k,rru0)
      do k=1,nk
      do j=1,nj
      do i=2,ni+1
        rru0 = 1.0/(0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k)))
        uten(i,j,k)=uten(i,j,k)+(turbx(i,j,k)+turbz(i,j,k))*rru0
      enddo
      enddo
      enddo

    ENDIF

!-------------------------------------------------------------------
!  All done

      if(timestats.ge.1) time_tmix=time_tmix+mytime()

      return
      end subroutine turbu


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbv(dt,xh,rxh,arh1,arh2,uh,xf,rvh,vf,mh,mf,rho,rr,rf,   &
                       zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gyv,  &
                       turbx,turby,turbz,dum1,dum2,dum3,v,vten,w,t12,t22,t23,kmv)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: rvh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,gzv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gyv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dum1,dum2,dum3
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: t12,t22,t23
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmv
 
      integer :: i,j,k,ip
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2,rrv0
      real :: cfa,cfb,cfc,cfd
      real :: foo1,foo2

!---------------------------------------------------------------

  IF(.not.terrain_flag)THEN

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        !  x-direction
        do j=1,nj+1
        do i=1,ni
          turbx(i,j,k)=(t12(i+1,j,k)-t12(i,j,k))*rdx*uh(i)
        enddo
        enddo

        !  y-direction
        do j=1,nj+1
        do i=1,ni
          turby(i,j,k)=(t22(i,j,k)-t22(i,j-1,k))*rdy*vf(j)
        enddo
        enddo

      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=(arh2(i)*arh2(i)*t12(i+1,j,k)-arh1(i)*arh1(i)*t12(i,j,k))*rdx*uh(i)
        enddo
        enddo

        do j=1,nj
        do i=1,ni
          turby(i,j,k)=0.0
        enddo
        enddo

      enddo

    ENDIF

!---------------------------------------------------------------
!  Terrain:

  ELSE

      ! dum1 stores t12 at w-pts:
      ! dum2 stores t22 at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do j=0,nj+1

          ! lowest model level:
          do i=1,ni+1
            dum1(i,j,1) = cgs1*t12(i,j,1)+cgs2*t12(i,j,2)+cgs3*t12(i,j,3)
            dum2(i,j,1) = cgs1*t22(i,j,1)+cgs2*t22(i,j,2)+cgs3*t22(i,j,3)
          enddo

          ! upper-most model level:
          do i=1,ni+1
            dum1(i,j,nk+1) = cgt1*t12(i,j,nk)+cgt2*t12(i,j,nk-1)+cgt3*t12(i,j,nk-2)
            dum2(i,j,nk+1) = cgt1*t22(i,j,nk)+cgt2*t22(i,j,nk-1)+cgt3*t22(i,j,nk-2)
          enddo

          ! interior:
          do k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do i=1,ni+1
            dum1(i,j,k) = r1*t12(i,j,k-1)+r2*t12(i,j,k)
            dum2(i,j,k) = r1*t22(i,j,k-1)+r2*t22(i,j,k)
          enddo
          enddo

      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do k=1,nk

        !  x-direction
        do j=1,nj+1
        do i=1,ni
          r1 = 0.25*((rgz(i-1,j-1)+rgz(i,j))+(rgz(i-1,j)+rgz(i,j-1)))
          r2 = 0.25*((rgz(i+1,j-1)+rgz(i,j))+(rgz(i+1,j)+rgz(i,j-1)))
          turbx(i,j,k)=gzv(i,j)*(t12(i+1,j,k)*r2-t12(i,j,k)*r1)*rdx*uh(i)      &
                      +0.5*( (zt-sigmaf(k+1))*(dum1(i,j,k+1)+dum1(i+1,j,k+1))  &
                            -(zt-sigmaf(k  ))*(dum1(i,j,k  )+dum1(i+1,j,k  ))  &
                           )*gzv(i,j)*(r2-r1)*rdx*uh(i)*rdsf(k)
        enddo
        enddo

        !  y-direction
        do j=1,nj+1
        do i=1,ni
          turby(i,j,k)=gzv(i,j)*(t22(i,j,k)*rgz(i,j)-t22(i,j-1,k)*rgz(i,j-1))*rdy*vf(j)  &
                      +0.5*( gyv(i,j,k+1)*(dum2(i,j-1,k+1)+dum2(i,j,k+1))                &
                            -gyv(i,j,k  )*(dum2(i,j-1,k  )+dum2(i,j,k  )) )*rdsf(k)
        enddo
        enddo

      enddo

  ENDIF  ! endif for terrain check

!-----------------------------------------------------------------
!  open boundary conditions:

    IF( wbc.eq.2 .or. ebc.eq.2 .or. sbc.eq.2 .or. nbc.eq.2 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk

        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj+1
            turbx(1,j,k) = 0.0
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj+1
            turbx(ni,j,k) = 0.0
          enddo
        ENDIF

        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni
            turby(i,1,k) = 0.0
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni
            turby(i,nj+1,k) = 0.0
          enddo
        ENDIF

      ENDDO
    ENDIF

!-----------------------------------------------------------------
!  z-direction

    IF( iturb.eq.3 .and. l_inf.lt.tsmall )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

    ELSE

      rdt = 0.5/dt
      tema = -0.0625*dt*vialpha*rdz*rdz
      temb =  0.0625*dt*vibeta*rdz*rdz
      temc =  0.5*dt*rdz

      ip = 1
      if( axisymm.eq.1 ) ip = 0

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2,cfa,cfb,cfc,cfd,tem,rrv0)
    do j=1,nj+ip

        !--------
        k = 1
        do i=1,ni
          rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
          tem = (mh(i,j-1,k)+mh(i,j,k))*rrv0
          r1 = 0.0
          r2 = (kmv(i,j-1,k+1)+kmv(i,j,k+1))*(mf(i,j-1,k+1)+mf(i,j,k+1))   &
              *(rf(i,j-1,k+1)+rf(i,j,k+1))*tem
          cfa = 0.0
          cfc = tema*r2
          cfb = 1.0 - cfc
          cfd = v(i,j,k) + temb*( r2*v(i,j,k+1)-r2*v(i,j,k) )
          tem = temc*t23(i,j,1)*(mh(i,j-1,1)+mh(i,j,1))*rrv0
          cfd = cfd - tem
          tem = 1.0/cfb
          dum1(i,j,1)=-cfc*tem
          dum2(i,j,1)= cfd*tem
        enddo
        !--------
        do k=2,nk-1
        do i=1,ni
          rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
          tem = (mh(i,j-1,k)+mh(i,j,k))*rrv0
          r1 = (kmv(i,j-1,k  )+kmv(i,j,k  ))*(mf(i,j-1,k  )+mf(i,j,k  ))   &
              *(rf(i,j-1,k  )+rf(i,j,k  ))*tem
          r2 = (kmv(i,j-1,k+1)+kmv(i,j,k+1))*(mf(i,j-1,k+1)+mf(i,j,k+1))   &
              *(rf(i,j-1,k+1)+rf(i,j,k+1))*tem
          cfa = tema*r1
          cfc = tema*r2
          cfb = 1.0 - cfa - cfc
          cfd = v(i,j,k) + temb*( r2*v(i,j,k+1)-(r1+r2)*v(i,j,k)+r1*v(i,j,k-1) )
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum1(i,j,k)=-cfc*tem
          dum2(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
        enddo
        enddo
        !--------
        k = nk
        do i=1,ni
          rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
          tem = (mh(i,j-1,k)+mh(i,j,k))*rrv0
          r1 = (kmv(i,j-1,k  )+kmv(i,j,k  ))*(mf(i,j-1,k  )+mf(i,j,k  ))   &
              *(rf(i,j-1,k  )+rf(i,j,k  ))*tem
          r2 = 0.0
          cfa = tema*r1
          cfc = 0.0
          cfb = 1.0 - cfa
          cfd = v(i,j,k) + temb*( -r1*v(i,j,k)+r1*v(i,j,k-1) )
          tem = temc*t23(i,j,nk+1)*(mh(i,j-1,nk)+mh(i,j,nk))*rrv0
          cfd = cfd + tem
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
!!!          dum1(i,j,k)=-cfc*tem
          dum2(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
          dum3(i,j,k)=dum2(i,j,k)
          turbz(i,j,k) = (rho(i,j-1,k)+rho(i,j,k))*(dum3(i,j,k)-v(i,j,k))*rdt
        enddo
        !--------

        do k=nk-1,1,-1
        do i=1,ni
          dum3(i,j,k)=dum1(i,j,k)*dum3(i,j,k+1)+dum2(i,j,k)
          turbz(i,j,k) = (rho(i,j-1,k)+rho(i,j,k))*(dum3(i,j,k)-v(i,j,k))*rdt
        enddo
        enddo

    enddo

  IF( terrain_flag )THEN
    ! dum1 stores w at scalar-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=0,nj+1
      do i=0,ni+1
        dum1(i,j,k)=0.5*(w(i,j,k)+w(i,j,k+1))
      enddo
      enddo
    ENDDO
  ENDIF

    IF(axisymm.eq.0)THEN
      ! explicit piece ... dwdy term
      ! Cartesian grid only:
    DO j=1,nj+ip
      do i=1,ni
        dum2(i,j,1) = 0.0
        dum2(i,j,nk+1) = 0.0
      enddo
      IF(.not.terrain_flag)THEN
        do k=2,nk
        do i=1,ni
          dum2(i,j,k)=(w(i,j,k)-w(i,j-1,k))*rdy*vf(j)
          dum2(i,j,k)=dum2(i,j,k)*0.25*( kmv(i,j-1,k)+kmv(i,j,k) )                 &
                                      *( rf(i,j-1,k)+rf(i,j,k) )
        enddo
        enddo
      ELSE
        do k=2,nk
        do i=1,ni
          dum2(i,j,k)=(w(i,j,k)*rgz(i,j)-w(i,j-1,k)*rgz(i,j-1))*rdy*vf(j)          &
                  +0.5*rds(k)*( (zt-sigma(k  ))*(dum1(i,j,k  )+dum1(i,j-1,k  ))    &
                               -(zt-sigma(k-1))*(dum1(i,j,k-1)+dum1(i,j-1,k-1)) )  &
                             *(rgz(i,j)-rgz(i,j-1))*rdy*vf(j)
          dum2(i,j,k)=dum2(i,j,k)*0.25*( kmv(i,j-1,k)+kmv(i,j,k) )                 &
                                      *( gz(i,j-1)*rf(i,j-1,k)+gz(i,j)*rf(i,j,k) )
        enddo
        enddo
      ENDIF
      do k=1,nk
      do i=1,ni
        turbz(i,j,k)=turbz(i,j,k)+(dum2(i,j,k+1)-dum2(i,j,k))*rdz*0.5*(mh(i,j-1,k)+mh(i,j,k))
      enddo
      enddo
    ENDDO
    ENDIF

  ENDIF

!-----------------------------------------------------------------

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,rrv0)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
        vten(i,j,k)=vten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))*rrv0
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        vten(i,j,k)=vten(i,j,k)+(turbx(i,j,k)+turbz(i,j,k))*rr(i,j,k)
      enddo
      enddo
      enddo

    ENDIF

!-------------------------------------------------------------------
!  All done
 
      if(timestats.ge.1) time_tmix=time_tmix+mytime()
 
      return
      end subroutine turbv


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
      subroutine turbw(dt,xh,rxh,arh1,arh2,uh,xf,vh,mh,mf,rho,rf,gz,rgzu,rgzv,rds,sigma,   &
                       turbx,turby,turbz,dum1,dum2,dum3,w,wten,t13,t23,t33,t22,kmh)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real :: dt
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, dimension(jb:je) :: vh
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: rho,rf
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgzu,rgzv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dum1,dum2,dum3
      real, dimension(ib:ie,jb:je,kb:ke+1) :: w,wten
      real, dimension(ib:ie,jb:je,kb:ke) :: t13,t23,t33,t22
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh
 
      integer :: i,j,k
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2,rrf
      real :: cfa,cfb,cfc,cfd
      real :: foo1,foo2

!----------------------------------------------------------------

  IF(.not.terrain_flag)THEN

    IF(axisymm.eq.0)THEN
      ! Cartesian without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        !  x-direction
        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=(t13(i+1,j,k)-t13(i,j,k))*rdx*uh(i)
        enddo
        enddo

        !  y-direction
        do j=1,nj
        do i=1,ni
          turby(i,j,k)=(t23(i,j+1,k)-t23(i,j,k))*rdy*vh(j)
        enddo
        enddo

      enddo

    ELSE
      ! axisymmetric:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=(arh2(i)*t13(i+1,j,k)-arh1(i)*t13(i,j,k))*rdx*uh(i)
        enddo
        enddo

        !  y-direction
        do j=1,nj
        do i=1,ni
          turby(i,j,k)=0.0
        enddo
        enddo

      enddo

    ENDIF

!----------------------------------------------------------------

  ELSE
      ! Cartesian with terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = 0.25*( (t13(i,j,k+1)+t13(i+1,j,k+1)) &
                              +(t13(i,j,k  )+t13(i+1,j,k  )) )
        enddo
        enddo
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = 0.25*( (t23(i,j,k+1)+t23(i,j+1,k+1)) &
                              +(t23(i,j,k  )+t23(i,j+1,k  )) )
        enddo
        enddo
      enddo


!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=gz(i,j)*( t13(i+1,j,k)*rgzu(i+1,j)             &
                                -t13(i  ,j,k)*rgzu(i  ,j) )*rdx*uh(i) &
              +( (zt-sigma(k  ))*dum1(i,j,k  )                        &
                -(zt-sigma(k-1))*dum1(i,j,k-1)                        &
               )*gz(i,j)*(rgzu(i+1,j)-rgzu(i,j))*rdx*uh(i)*rds(k)
        enddo
        enddo

        do j=1,nj
        do i=1,ni
          turby(i,j,k)=gz(i,j)*( t23(i,j+1,k)*rgzv(i,j+1)             &
                                -t23(i,j  ,k)*rgzv(i,j  ) )*rdy*vh(j) &
              +( (zt-sigma(k  ))*dum2(i,j,k  )                        &
                -(zt-sigma(k-1))*dum2(i,j,k-1)                        &
               )*gz(i,j)*(rgzv(i,j+1)-rgzv(i,j))*rdy*vh(j)*rds(k)
        enddo
        enddo

      enddo

  ENDIF

!-----------------------------------------------------------------
!  open boundary conditions:

    IF( wbc.eq.2 .or. ebc.eq.2 .or. sbc.eq.2 .or. nbc.eq.2 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=2,nk

        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj
            turbx(1,j,k) = 0.0
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj
            turbx(ni,j,k) = 0.0
          enddo
        ENDIF

        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni
            turby(i,1,k) = 0.0
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni
            turby(i,nj,k) = 0.0
          enddo
        ENDIF

      ENDDO
    ENDIF

!-----------------------------------------------------------------
!  z-direction

    IF( iturb.eq.3 .and. l_inf.lt.tsmall )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

    ELSE

      rdt = 1.0/dt
      tema = -1.0*dt*vialpha*rdz*rdz
      temb =      dt*vibeta*rdz*rdz
      temc = dt*rdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2,cfa,cfb,cfc,cfd,tem,rrf)
      do j=1,nj

        k=2
        do i=1,ni
          rrf = mf(i,j,k)/rf(i,j,k)
          r2 = (kmh(i,j,k  )+kmh(i,j,k+1))*mh(i,j,k  )*rho(i,j,k  )*rrf
          cfc = tema*r2
          cfb = 1.0 - cfc
          cfd = w(i,j,k) + temb*( r2*w(i,j,k+1)-r2*w(i,j,k) )
          cfd = cfd - temc*t33(i,j,k-1)*rrf
          tem = 1.0/cfb
          dum1(i,j,k) = -cfc*tem
          dum2(i,j,k) =  cfd*tem
        enddo

        do k=3,(nk-1)
        do i=1,ni
          rrf = mf(i,j,k)/rf(i,j,k)
          r1 = (kmh(i,j,k-1)+kmh(i,j,k  ))*mh(i,j,k-1)*rho(i,j,k-1)*rrf
          r2 = (kmh(i,j,k  )+kmh(i,j,k+1))*mh(i,j,k  )*rho(i,j,k  )*rrf
          cfa = tema*r1
          cfc = tema*r2
          cfb = 1.0 - cfa - cfc
          cfd = w(i,j,k) + temb*(r2*w(i,j,k+1)-(r1+r2)*w(i,j,k)+r1*w(i,j,k-1))
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum1(i,j,k) = -cfc*tem
          dum2(i,j,k) = (cfd-cfa*dum2(i,j,k-1))*tem
        enddo
        enddo

        k = nk
        do i=1,ni
          rrf = mf(i,j,k)/rf(i,j,k)
          r1 = (kmh(i,j,k-1)+kmh(i,j,k  ))*mh(i,j,k-1)*rho(i,j,k-1)*rrf
          cfa = tema*r1
          cfb = 1.0 - cfa
          cfd = w(i,j,k) + temb*( -r1*w(i,j,k)+r1*w(i,j,k-1) )
          cfd = cfd + temc*t33(i,j,k)*rrf
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum2(i,j,k) = (cfd-cfa*dum2(i,j,k-1))*tem
          !---
          dum3(i,j,k) = dum2(i,j,k)
          turbz(i,j,k) = rf(i,j,k)*(dum3(i,j,k)-w(i,j,k))*rdt
        enddo

        do k=(nk-1),2,-1
        do i=1,ni
          dum3(i,j,k) = dum1(i,j,k)*dum3(i,j,k+1)+dum2(i,j,k)
          turbz(i,j,k) = rf(i,j,k)*(dum3(i,j,k)-w(i,j,k))*rdt
        enddo
        enddo

      enddo

    ENDIF

!-----------------------------------------------------------------

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        wten(i,j,k)=wten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))/rf(i,j,k)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        wten(i,j,k)=wten(i,j,k)+(turbx(i,j,k)+turbz(i,j,k))/rf(i,j,k)
      enddo
      enddo
      enddo

    ENDIF

!-------------------------------------------------------------------
!  All done

      if(timestats.ge.1) time_tmix=time_tmix+mytime()
 
      return
      end subroutine turbw

