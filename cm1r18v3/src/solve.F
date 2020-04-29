!-----------------------------------------------------------------------
!
!  CM1 Numerical Model, Release 18.3  (cm1r18.3)
!  7 October 2015
!  http://www2.mmm.ucar.edu/people/bryan/cm1/
!
!  Copyright (c) 2000-2015 by George H. Bryan, National Center for 
!  Atmospheric Research, Boulder, Colorado, USA. 
!
!-----------------------------------------------------------------------
!  Quick Index:
!    ua/u3d     = velocity in x-direction (m/s)
!    va/v3d     = velocity in y-direction (m/s)
!    wa/w3d     = velocity in z-direction (m/s)
!    tha/th3d   = perturbation potential temperature (K)
!    ppi/pp3d   = perturbation nondimensional pressure ("Exner function")
!    qa/q3d     = mixing ratios of moisture (kg/kg)
!    tkea/tke3d = SUBGRID turbulence kinetic energy (m^2/s^2)
!    kmh/kmv    = turbulent diffusion coefficients for momentum (m^2/s)
!    khh/khv    = turbulent diffusion coefficients for scalars (m^2/s)
!                 (h = horizontal, v = vertical)
!    prs        = pressure (Pa)
!    rho        = density (kg/m^3)
!
!    th0,pi0,prs0,etc = base-state arrays
!
!    xh         = x (m) at scalar points
!    xf         = x (m) at u points
!    yh         = y (m) at scalar points
!    yf         = y (m) at v points
!    zh         = z (m above sea level) of scalar points (aka, "half levels")
!    zf         = z (m above sea level) of w points (aka, "full levels")
!
!    For the axisymmetric model (axisymm=1), xh and xf are radius (m).
!
!  See "The governing equations for CM1" for more details:
!        http://www2.mmm.ucar.edu/people/bryan/cm1/cm1_equations.pdf
!-----------------------------------------------------------------------
!  Some notes:
!
!  - Upon entering solve, the arrays ending in "a" (eg, ua,wa,tha,qa,etc)
!    are equivalent to the arrays ending in "3d" (eg, u3d,w3d,th3d,q3d,etc).
!  - The purpose of solve is to update the variables from time "t" to time
!    "t+dt".  Values at time "t+dt" are stored in the "3d" arrays.
!  - The "ghost zones" (boundaries beyond the computational subdomain) are
!    filled out completely (3 rows/columns) for the "3d" arrays.  To save 
!    unnecessary computations, starting with cm1r15 the "ghost zones" of 
!    the "a" arrays are only filled out to 1 row/column.  Hence, if you 
!    need to do calculations that use a large stencil, you must use the 
!    "3d" arrays (not the "a") arrays.
!  - Arrays named "ten" store tendencies.  Those ending "ten1" store
!    pre-RK tendencies that are calculated once and then held fixed during
!    the RK (Runge-Kutta) sub-steps. 
!  - CM1 uses a low-storage three-step Runge-Kutta scheme.  See Wicker
!    and Skamarock (2002, MWR, p 2088) for more information.
!  - CM1 uses a staggered C grid.  Hence, u arrays have one more grid point
!    in the i direction, v arrays have one more grid point in the j 
!    direction, and w arrays have one more grid point in the k direction
!    (compared to scalar arrays).
!  - CM1 assumes the subgrid turbulence parameters (tke,km,kh) are located
!    at the w points. 
!-----------------------------------------------------------------------

      subroutine solve(nstep,nrec,prec,nwrite,nrst,rbufsz,num_soil_layers,ndt,ntdiag,nqdiag, &
                   dt,dtlast,mtime,stattim,taptim,rsttim,radtim,prcltim,adt,acfl,dbldt,mass1,mass2, &
                   dosfcflx,cloudvar,rhovar,qname,budname,bud,bud2,qbudget,asq,bsq, &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,yh,vh,rvh,yf,vf,rvf,   &
                   xfref,yfref,dumk1,dumk2,rds,sigma,rdsf,sigmaf,    &
                   tauh,taus,zh,mh,rmh,c1,c2,tauf,zf,mf,rmf,         &
                   rstat,rho0s,pi0s,prs0s,rth0s,pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0, &
                   qi0,rr0,rf0,rrf0,                                 &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv, &
                   rain,sws,svs,sps,srs,sgs,sus,shs,                 &
                   tsk,thflux,qvflux,cd,ch,cq,u1,v1,s1,tlh,          &
                   radbcw,radbce,radbcs,radbcn,                      &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,          &
                   divx,rho,rr,rf,prs,                               &
                   t11,t12,t13,t22,t23,t33,                          &
                   u0,rru,us,ua,u3d,uten,uten1,                      &
                   v0,rrv,vs,va,v3d,vten,vten1,                      &
                   rrw,ws,wa,w3d,wten,wten1,                         &
                   ppi,pp3d,piadv,ppten,sten,ppx,                    &
                   tha,th3d,thadv,thten,thten1,thterm,               &
                   qpten,qtten,qvten,qcten,qa,q3d,qten,              &
                   kmh,kmv,khh,khv,tkea,tke3d,tketen,                &
                   nm,defv,defh,dissten,                             &
                   thpten,qvpten,qcpten,qipten,upten,vpten,          &
                   swten,lwten,o30,radsw,rnflx,radswnet,radlwin,dsr,olr,rad2d,  &
                   effc,effi,effs,effr,effg,effis,                   &
                   lu_index,kpbl2d,psfc,u10,v10,s10,hfx,qfx,xland,znt,ust,  &
                   hpbl,wspd,psim,psih,gz1oz0,br,                    &
                   CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                    &
                   MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                 &
                   CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,    &
                   f2d,gsw,glw,chklowq,capg,snowc,dsxy,wstar,delta,fm,fh,  &
                   mznt,smois,taux,tauy,hpbl2d,evap2d,heat2d,rc2d,   &
                   slab_zs,slab_dzs,tslb,tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,       &
                   pta,pt3d,ptten,pdata,packet,ploc,                 &
                   cfb,cfa,cfc,ad1,ad2,pdt,deft,rhs,trans,flag,      &
                   reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_tk,reqs_q,reqs_t, &
                   nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                  &
                   n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,          &
                   ww1,ww2,we1,we2,ws1,ws2,wn1,wn2,                  &
                   pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,                  &
                   vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,                  &
                   uw31,uw32,ue31,ue32,us31,us32,un31,un32,          &
                   vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,          &
                   ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,          &
                   sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,          &
                   rw31,rw32,re31,re32,rs31,rs32,rn31,rn32,          &
                   qw31,qw32,qe31,qe32,qs31,qs32,qn31,qn32,          &
                   tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,          &
                   kw1,kw2,ke1,ke2,ks1,ks2,kn1,kn2,                  &
                   tw1,tw2,te1,te2,ts1,ts2,tn1,tn2,                  &
                   dat1,dat2,dat3,reqt,tdiag,qdiag)
        ! end_solve
      use module_mp_thompson
      use module_mp_graupel
      use module_mp_nssl_2mom, only : zscale, nssl_2mom_driver
#ifdef MPI
      use mpi
#endif
      use module_restart
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'radcst.incl'
      include 'timestat.incl'

!-----------------------------------------------------------------------
! Arrays and variables passed into solve

      integer, intent(in) :: nstep
      integer, intent(inout) :: nrec,prec,nwrite,nrst
      integer, intent(in) :: rbufsz,num_soil_layers
      integer, intent(inout) :: ndt
      integer, intent(in) :: ntdiag,nqdiag
      real, intent(inout) :: dt,dtlast
      double precision, intent(inout) :: mtime,stattim,taptim,rsttim,radtim,prcltim
      double precision, intent(inout) :: adt,acfl,dbldt
      double precision, intent(in   ) :: mass1
      double precision, intent(inout) :: mass2
      logical, intent(in) :: dosfcflx
      logical, intent(in), dimension(maxq) :: cloudvar,rhovar
      character*3, intent(in), dimension(maxq) :: qname
      character*6, intent(in), dimension(maxq) :: budname
      double precision, intent(inout), dimension(nk) :: bud
      double precision, intent(inout), dimension(nj) :: bud2
      double precision, intent(inout), dimension(nbudget) :: qbudget
      double precision, intent(inout), dimension(numq) :: asq,bsq
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf,ruf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf,rvf
      real, intent(in), dimension(-2:nx+4) :: xfref
      real, intent(in), dimension(-2:ny+4) :: yfref
      double precision, intent(inout), dimension(kb:ke) :: dumk1,dumk2
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: tauh,taus,zh,mh,rmh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: tauf,zf,mf,rmf
      real, intent(inout), dimension(stat_out) :: rstat
      real, intent(in), dimension(ib:ie,jb:je) :: rho0s,pi0s,prs0s,rth0s
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: qi0,rr0,rf0,rrf0
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(inout), dimension(ib:ie,jb:je,nrain) :: rain,sws,svs,sps,srs,sgs,sus,shs
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,znt,ust,thflux,qvflux,cd,ch,cq,u1,v1,s1,xland,psfc,tlh
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: divx,rho,rr,rf,prs
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: rru,us,ua,u3d,uten,uten1
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: rrv,vs,va,v3d,vten,vten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw,ws,wa,w3d,wten,wten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,piadv,ppten,sten,ppx
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tha,th3d,thadv,thten,thten1,thterm
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem) :: qpten,qtten,qvten,qcten
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa,q3d,qten
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tke3d,tketen
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh,dissten
      real, intent(inout), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: thpten,qvpten,qcpten,qipten,upten,vpten
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten
      real, intent(in), dimension(ibr:ier,jbr:jer,kbr:ker) :: o30
      real, intent(inout), dimension(ni,nj) :: radsw,rnflx,radswnet,radlwin,dsr,olr
      real, intent(inout), dimension(ni,nj,nrad2d) :: rad2d
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: effc,effi,effs,effr,effg,effis
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: lu_index
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: kpbl2d
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: u10,v10,s10,hfx,qfx, &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,       &
                                      CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      f2d,gsw,glw,chklowq,capg,snowc,dsxy,wstar,delta,fm,fh
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: mznt,smois,taux,tauy,hpbl2d,evap2d,heat2d,rc2d
      real, intent(in), dimension(num_soil_layers) :: slab_zs,slab_dzs
      real, intent(inout), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, intent(inout), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta,pt3d,ptten
      real, intent(inout), dimension(npvals,nparcels) :: pdata,packet
      real, intent(inout), dimension(3,nparcels) :: ploc
      real, intent(in), dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: cfb
      real, intent(in), dimension(kpb:kpe) :: cfa,cfc,ad1,ad2
      complex, intent(inout), dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: pdt,deft
      complex, intent(inout), dimension(ipb:ipe,jpb:jpe) :: rhs,trans
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_tk
      integer, intent(inout), dimension(rmp,numq) :: reqs_q
      integer, intent(inout), dimension(rmp,npt) :: reqs_t
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(cmp,cmp,kmt+1) :: n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2
      real, intent(inout), dimension(jmp,kmp-1) :: ww1,ww2,we1,we2
      real, intent(inout), dimension(imp,kmp-1) :: ws1,ws2,wn1,wn2
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, intent(inout), dimension(jmp,kmp) :: vw1,vw2,ve1,ve2
      real, intent(inout), dimension(imp,kmp) :: vs1,vs2,vn1,vn2
      real, intent(inout), dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, intent(inout), dimension(imp,cmp,kmp)   :: vs31,vs32,vn31,vn32
      real, intent(inout), dimension(cmp,jmp,kmp-1) :: ww31,ww32,we31,we32
      real, intent(inout), dimension(imp,cmp,kmp-1) :: ws31,ws32,wn31,wn32
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, intent(inout), dimension(cmp,jmp,kmp,2) :: rw31,rw32,re31,re32
      real, intent(inout), dimension(imp,cmp,kmp,2) :: rs31,rs32,rn31,rn32
      real, intent(inout), dimension(cmp,jmp,kmp,numq) :: qw31,qw32,qe31,qe32
      real, intent(inout), dimension(imp,cmp,kmp,numq) :: qs31,qs32,qn31,qn32
      real, intent(inout), dimension(cmp,jmp,kmt)   :: tkw1,tkw2,tke1,tke2
      real, intent(inout), dimension(imp,cmp,kmt)   :: tks1,tks2,tkn1,tkn2
      real, intent(inout), dimension(jmp,kmt,4)     :: kw1,kw2,ke1,ke2
      real, intent(inout), dimension(imp,kmt,4)     :: ks1,ks2,kn1,kn2
      real, intent(inout), dimension(cmp,jmp,kmp,npt) :: tw1,tw2,te1,te2
      real, intent(inout), dimension(imp,cmp,kmp,npt) :: ts1,ts2,tn1,tn2
      real, intent(inout), dimension(ni+1,nj+1) :: dat1
      real, intent(inout), dimension(d2i,d2j) :: dat2
      real, intent(inout), dimension(d3i,d3j,d3n) :: dat3
      integer, intent(inout), dimension(d3t) :: reqt
      real, intent(inout) , dimension(ibd:ied,jbd:jed,kbd:ked,ntdiag) :: tdiag
      real, intent(inout) , dimension(ibd:ied,jbd:jed,kbd:ked,nqdiag) :: qdiag

!-----------------------------------------------------------------------
! Arrays and variables defined inside solve

      integer :: i,j,k,n,nrk,bflag,pdef,nn,fnum,diffit,k1

      real :: delqv,delpi,delth,delt,fac,epsd,dheat,dz1,xs
      real :: foo1,foo2
      real :: dttmp,rtime,rdt,tem,tem0,thrad,prad,ql
      real :: cpm,cvm
      real :: r1,r2,tnew,pnew,pinew,thnew,qvnew
      real :: gamm,aiu

      double precision :: weps,afoo,bfoo,p0,p2,tout

      logical :: getdbz,getvt
      logical :: dorad,reset
      logical :: doirrad,dosorad
      logical :: get_time_avg,dostat,dowrite,dorestart,doprclout

      real :: saltitude,sazimuth,zen
      real, dimension(2) :: x1
      real, dimension(2) :: y1

      ! 1d arrays for radiation scheme:
      real, dimension(rbufsz) :: radbuf
      real, dimension(nkr) :: swtmp,lwtmp
      real, dimension(nkr) :: tem1,tem2,tem3,tem4,tem5,   &
                              tem6,tem7,tem8,tem9,tem10,   &
                              tem11,tem12,tem13,tem14,tem15,   &
                              tem16,tem17
      real, dimension(nkr) :: teffc,teffi,teffs,teffr,teffg,teffis
      real, dimension(nkr) :: ptprt,pprt,qv,qc,qr,qi,qs,qh,cvr,   &
                              ptbar,pbar,appi,rhostr,zpp,o31

      real :: rad2deg,albedo,albedoz,tema,temb,frac_snowcover

      real, parameter  ::  cfl_limit   =  0.95    ! maximum CFL allowed  (actually a "target" value)
      real, parameter  ::  max_change  =  0.10    ! maximum (percentage) change in timestep

!--------------------------------------------------------------------

      if( adapt_dt.eq.1 .and. myid.eq.0 ) print *,'cflmax,dt,nsound:',cflmax,dt,nsound

#ifdef MPI
      nf=0
      nu=0
      nv=0
      nw=0
#endif

      afoo=0.0d0
      bfoo=0.0d0

      dowrite = .false.

      ! check if writeout will be called at end of timestep:
      if( (mtime+dt).ge.(taptim-0.1*dt) ) dowrite = .true.
      IF( myid.eq.0 .and. dowrite ) print *,'  dowrite = ',dowrite

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   radiation  ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IF( radopt.ge.1 )THEN

        ! just to be sure:
        call setradwrk(nir,njr,nkr)

        ! time at beginning of timestep:
        rtime=sngl(mtime)
        dorad = .false.
        IF( mtime.ge.(radtim-0.1*dt) ) dorad = .true.
        dtrad = max( dtrad , dt )

        IF( dorad )THEN
          i = 1
          j = 1
          rtime=sngl(mtime+dt)
          if(dowr) write(outfile,*) '  Calculating radiation tendency:'
          if(timestats.ge.1) time_rad=time_rad+mytime()
          call bcs(prs)
          CALL zenangl( ni,nj,      zf(1,1,1),    &
                rad2d(1,1,ncosz), rad2d(1,1,ncosss), radsw,              &
                dum1(1,1,1),dum1(1,1,2),dum1(1,1,3),dum1(1,1,4),        &
                dum2(1,1,1),dum2(1,1,2),dum2(1,1,3),dum2(1,1,4),        &
                saltitude,sazimuth,dx,dy,dt,rtime,                     &
                ctrlat,ctrlon,year,month,day,hour,minute,second,jday )
          if(myid.eq.0)then
            print *,'    solar zenith angle  (degrees) = ',   &
                                   acos(rad2d(ni,nj,ncosz))*degdpi
            print *,'    solar azimuth angle (degrees) = ',sazimuth*degdpi
          endif
!-----------------------------------------------------------------------
!
!  Calculate surface albedo which is dependent on solar zenith angle
!  and soil moisture. Set the albedo for different types of solar
!  flux to be same.
!
!    rsirbm   Solar IR surface albedo for beam radiation
!    rsirdf   Solar IR surface albedo for diffuse radiation
!    rsuvbm   Solar UV surface albedo for beam radiation
!    rsuvdf   Solar UV surface albedo for diffuse radiation
!
!-----------------------------------------------------------------------
!
  rad2deg = 180.0/3.141592654

          radbuf = 0.0
          tem1 = 0.0
          tem2 = 0.0
          tem3 = 0.0
          tem4 = 0.0
          tem5 = 0.0
          tem6 = 0.0
          tem7 = 0.0
          tem8 = 0.0
          tem9 = 0.0
          tem10 = 0.0
          tem11 = 0.0
          tem12 = 0.0
          tem13 = 0.0
          tem14 = 0.0
          tem15 = 0.0
          tem16 = 0.0
          tem17 = 0.0
          ptprt = 0.0
          pprt = 0.0
          qv = 0.0
          qc = 0.0
          qr = 0.0
          qi = 0.0
          qs = 0.0
          qh = 0.0
          cvr = 0.0
          ptbar = 0.0
          pbar = 0.0
          appi = 0.0
          rhostr = 0.0
          zpp = 0.0
          o31 = 0.0

!$omp parallel do default(shared)  &
!$omp private(i,j,albedo,albedoz,frac_snowcover,tema)
  DO j=1,nj
    DO i=1,ni

      ! let's just use MM5/WRF value, instead:
      albedo = albd(i,j)

      ! arps code for albedo:
      ! (not sure I trust this.....)

!      albedoz = 0.01 * ( EXP( 0.003286         & ! zenith dependent albedo
!          * SQRT( ( ACOS(rad2d(i,j,ncosz))*rad2deg ) ** 3 ) ) - 1.0 )
!
!      IF ( soilmodel == 0 ) THEN             ! soil type not defined
!!!!        stop 12321
!        tema = 0
!      ELSE
!        tema = qsoil(i,j,1)/wsat(soiltyp(i,j))
!      END IF
!
!      frac_snowcover = MIN(snowdpth(i,j)/snowdepth_crit, 1.0)
!
!      IF ( tema > 0.5 ) THEN
!        albedo = albedoz + (1.-frac_snowcover)*0.14                     &
!                         + frac_snowcover*snow_albedo
!      ELSE
!        albedo = albedoz + (1.-frac_snowcover)*(0.31 - 0.34 * tema)     &
!                         + frac_snowcover*snow_albedo
!      END IF
!        albedo = albedoz

      rad2d(i,j,nrsirbm) = albedo
      rad2d(i,j,nrsirdf) = albedo
      rad2d(i,j,nrsuvbm) = albedo
      rad2d(i,j,nrsuvdf) = albedo

    END DO
  END DO
          ! big OpenMP parallelization loop:
!$omp parallel do default(shared)  &
!$omp private(i,j,k,ptprt,pprt,qv,qc,qr,qi,qs,qh,cvr,appi,o31,        &
!$omp tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,tem10,        &
!$omp tem11,tem12,tem13,tem14,tem15,tem16,tem17,radbuf,swtmp,lwtmp,   &
!$omp doirrad,dosorad,zpp,ptbar,pbar,rhostr,x1,y1,  &
!$omp teffc,teffi,teffs,teffr,teffg,teffis)
        do j=1,nj
        do i=1,ni
          swtmp = 0.0
          lwtmp = 0.0
          do k=1,nk+2
            ptprt(k) =  tha(i,j,k-1)
             pprt(k) =  prs(i,j,k-1) - prs0(i,j,k-1)
               qv(k) =   qa(i,j,k-1,nqv)
               qc(k) =   0.0
               if( nqc.ge.1 ) qc(k) = qa(i,j,k-1,nqc)
               qr(k) =   0.0
               if( nqr.ge.1 ) qr(k) = qa(i,j,k-1,nqr)
               qi(k) =   0.0
               if( nqi.ge.1 ) qi(k) = qa(i,j,k-1,nqi)
               qs(k) =   0.0
               if( nqs.ge.1 ) qs(k) = qa(i,j,k-1,nqs)
               qh(k) =   0.0
               if( nqg.ge.1 ) qh(k) = qa(i,j,k-1,nqg)
              cvr(k) = cv+cvv*qv(k)+cpl*(qc(k)+qr(k))+cpi*(qi(k)+qs(k)+qh(k))
             appi(k) =  pi0(i,j,k-1) + ppi(i,j,k-1)
              o31(k) =  o30(i,j,k-1)
            teffc(k) = effc(i,j,k-1)
            teffi(k) = effi(i,j,k-1)
            teffs(k) = effs(i,j,k-1)
            teffr(k) = effr(i,j,k-1)
            teffg(k) = effg(i,j,k-1)
           teffis(k) = effis(i,j,k-1)
          enddo
          ptprt(1) = ptprt(2)
           pprt(1) =  pprt(2)
          ptprt(nk+2) = ptprt(nk+1)
           pprt(nk+2) =  pprt(nk+1)
          x1(1) = xf(i)
          x1(2) = xf(i+1)
          y1(1) = yf(j)
          y1(2) = yf(j+1)
          do k=1,nk+3
            zpp(k) =   zf(i,j,k-1)
          enddo
          do k=1,nk+2
            ptbar(k) =  th0(i,j,k-1)
             pbar(k) = prs0(i,j,k-1)
           rhostr(k) =  rho(i,j,k-1)
          enddo
            ptbar(1) = rth0s(i,j)**(-1)
             pbar(1) = prs0s(i,j)
           rhostr(1) = rho0s(i,j)
            doirrad = .true.
            dosorad = .true.
          CALL radtrns(nir,njr,nkr, rbufsz, 0,myid,dx,dy,            &
                 ib,ie,jb,je,kb,ke,xh,yh,prs0s(i,j),olr(i,j),dsr(i,j),  &  ! MS add olr,dsr
                 ptprt,pprt,qv,qc,qr,qi,qs,qh,cvr,                      &
                 ptbar,pbar,appi,o31,rhostr, tsk(i,j), zpp ,                                 &
                 radsw(i,j),rnflx(i,j),radswnet(i,j),radlwin(i,j), rad2d(i,j,ncosss),            &
                 rad2d(i,j,nrsirbm),rad2d(i,j,nrsirdf),rad2d(i,j,nrsuvbm),                       &
                 rad2d(i,j,nrsuvdf), rad2d(i,j,ncosz),sazimuth,                                  &
                 rad2d(i,j,nfdirir),rad2d(i,j,nfdifir),rad2d(i,j,nfdirpar),rad2d(i,j,nfdifpar),  &
                 tem1, tem2, tem3, tem4, tem5,                &
                 tem6, tem7, tem8, tem9, tem10,               &
                 tem11,tem12,tem13,tem14,tem15,tem16,         &
                 radbuf(1), tem17,swtmp,lwtmp,doirrad,dosorad, &
                 teffc,teffi,teffs,teffr,teffg,teffis,        &
                 cgs1,cgs2,cgs3,cgt1,cgt2,cgt3,ptype,g,cp)
          do k=1,nk
            swten(i,j,k) = swtmp(k+1)
            lwten(i,j,k) = lwtmp(k+1)
          enddo
        enddo
        enddo
          radtim=radtim+dtrad
        ENDIF
        if(timestats.ge.1) time_rad=time_rad+mytime()

      ENDIF


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   subgrid turbulence schemes  cccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!--------------------------------------------------------------------
!  get RHS for tke scheme:

      IF(iturb.eq.1)THEN
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        DO k=2,nk
          !  Buoyancy, Dissipation, and Shear terms:
          do j=1,nj
          do i=1,ni
            tketen(i,j,k) = -khv(i,j,k)*nm(i,j,k)  &
                            -dissten(i,j,k)
          enddo
          enddo
          ! Shear term 
          IF(tconfig.eq.1)THEN
            do j=1,nj
            do i=1,ni
              tketen(i,j,k)=tketen(i,j,k)+kmv(i,j,k)*max(0.0,(defv(i,j,k)+defh(i,j,k)))
            enddo
            enddo
          ELSEIF(tconfig.eq.2)THEN
            do j=1,nj
            do i=1,ni
              tketen(i,j,k)=tketen(i,j,k)+kmv(i,j,k)*max(0.0,defv(i,j,k))   &
                                         +kmh(i,j,k)*max(0.0,defh(i,j,k))
            enddo
            enddo
          ENDIF
        ENDDO
        if(timestats.ge.1) time_turb=time_turb+mytime()

        call turbt(dt,xh,rxh,uh,xf,uf,vh,vf,mh,mf,rho,rr,rf,  &
                   rds,sigma,gz,rgz,gzu,rgzu,gzv,rgzv,        &
                   dum1,dum2,dum3,dum4,dum5,sten,tkea,tketen,kmh,kmv)

      ENDIF


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   Pre-RK calculations   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!--------------------------------------------------------------------
!  radbc
 
      if(irbc.eq.1)then

        if(ibw.eq.1 .or. ibe.eq.1) call radbcew(radbcw,radbce,ua)
 
        if(ibs.eq.1 .or. ibn.eq.1) call radbcns(radbcs,radbcn,va)

      endif

!--------------------------------------------------------------------
!  U-equation

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
!!!        uten1(i,j,k)=0.
        uten1(i,j,k)=-rdalpha*0.5*(tauh(i-1,j,k)+tauh(i,j,k))*(ua(i,j,k)-u0(i,j,k))
      enddo
      enddo
      enddo

      IF( iinit.eq.10 .and. mtime.lt.t2_uforce )THEN
        ! u-forcing for squall-line initialization:
        ! (Morrison et al, 2015, JAS, pg 315)
        gamm = 1.0
        if(mtime.ge.t1_uforce)THEN
          gamm = 1.0+(0.0-1.0)*(mtime-t1_uforce)/(t2_uforce-t1_uforce)
        endif
        if(myid.eq.0) print *,'  mtime,gamm = ',mtime,gamm
!$omp parallel do default(shared)  &
!$omp private(i,j,k,aiu)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          if( abs(xf(i)-xc_uforce).lt.xr_uforce .and. abs(zf(i,j,k)-zs(i,j)).lt.zr_uforce )then
            aiu = alpha_uforce*cos(0.5*pi*(xf(i)-xc_uforce)/xr_uforce)   &
                              *((cosh(2.5*(zf(i,j,k)-zs(i,j))/zr_uforce))**(-2))
            uten1(i,j,k)=uten1(i,j,k)+gamm*aiu
          endif
        enddo
        enddo
        enddo
      ENDIF
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

      if(idiff.ge.1)then
        if(difforder.eq.2)then
          call diff2u(1,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                      dum1,dum2,dum3,dum4,uten1,ust,rho,rr,rf,divx,t11,t12,t13)
        endif
      endif

      if(dns.eq.1)then
        call diff2u(2,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                    dum1,dum2,dum3,dum4,uten1,ust,rho,rr,rf,divx,t11,t12,t13)
      endif
 
      if(iturb.ge.1)then
        call turbu(dt,xh,ruh,xf,rxf,arf1,arf2,uf,vh,mh,mf,rmf,rho,rf,  &
                   zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gxu,     &
                   dum1,dum2,dum3,dum4,dum5,dum6,ua,uten1,wa,t11,t12,t13,t22,kmv)
      endif

      if(ipbl.ge.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          uten1(i,j,k) = uten1(i,j,k) + 0.5*( upten(i-1,j,k)+ upten(i,j,k))
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
      endif

!--------------------------------------------------------------------
!  V-equation
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
!!!        vten1(i,j,k)=0.
        vten1(i,j,k)=-rdalpha*0.5*(tauh(i,j-1,k)+tauh(i,j,k))*(va(i,j,k)-v0(i,j,k))
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

      if(idiff.ge.1)then
        if(difforder.eq.2)then
          call diff2v(1,xh,arh1,arh2,uh,rxf,arf1,arf2,uf,vh,vf,mh,mf,  &
                      dum1,dum2,dum3,dum4,vten1,ust,rho,rr,rf,divx,t22,t12,t23)
        endif
      endif

      if(dns.eq.1)then
        call diff2v(2,xh,arh1,arh2,uh,rxf,arf1,arf2,uf,vh,vf,mh,mf,  &
                    dum1,dum2,dum3,dum4,vten1,ust,rho,rr,rf,divx,t22,t12,t23)
      endif
 
      if(iturb.ge.1)then
        call turbv(dt,xh,rxh,arh1,arh2,uh,xf,rvh,vf,mh,mf,rho,rr,rf,   &
                   zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gyv,  &
                   dum1,dum2,dum3,dum4,dum5,dum6,va,vten1,wa,t12,t22,t23,kmv)
      endif

      if(ipbl.ge.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          vten1(i,j,k) = vten1(i,j,k) + 0.5*( vpten(i,j-1,k)+ vpten(i,j,k))
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
      endif
 
!--------------------------------------------------------------------
!  W-equation

!$omp parallel do default(shared)  &
!$omp private(i,j,k,xs)
      do k=1,nk+1
      do j=1,nj
      do i=1,ni
!!!        wten1(i,j,k)=0.0
        wten1(i,j,k)=-rdalpha*tauf(i,j,k)*wa(i,j,k)
!!!        ! forcing term from Nolan (2005, JAS):
!!!        xs = sqrt( ((zf(i,j,k)-3000.0)**2)/(2000.0**2) &
!!!                  +(xh(i)**2)/(1000.0**2) )
!!!        if( xs.lt.1.0 ) wten1(i,j,k)=wten1(i,j,k)+1.26*cos(0.5*pi*xs)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

      if(idiff.ge.1)then
        if(difforder.eq.2)then
          call diff2w(1,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                      dum1,dum2,dum3,dum4,wten1,rho,rr,rf,divx,t33,t13,t23)
        endif
      endif

      if(dns.eq.1)then
        call diff2w(2,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                    dum1,dum2,dum3,dum4,wten1,rho,rr,rf,divx,t33,t13,t23)
      endif
 
      if(iturb.ge.1)then
        call turbw(dt,xh,rxh,arh1,arh2,uh,xf,vh,mh,mf,rho,rf,gz,rgzu,rgzv,rds,sigma,   &
                   dum1,dum2,dum3,dum4,dum5,dum6,wa,wten1,t13,t23,t33,t22,kmh)
      endif

!--------------------------------------------------------------------
!  Arrays for vimpl turbs:
!    NOTE:  do not change dum7,dum8 from here to RK loop

      if(iturb.ge.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum7(i,j,k) = khv(i,j,k  )*mf(i,j,k  )*rf(i,j,k  )*mh(i,j,k)*rr(i,j,k)
          dum8(i,j,k) = khv(i,j,k+1)*mf(i,j,k+1)*rf(i,j,k+1)*mh(i,j,k)*rr(i,j,k)
        enddo
        enddo
        enddo
      endif

!--------------------------------------------------------------------
!  THETA-equation

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
!!!        thten1(i,j,k)=0.0
        thten1(i,j,k)=-rdalpha*taus(i,j,k)*tha(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

      if(idiff.eq.1)then
        if(difforder.eq.2)then
          call diff2s(1,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                      dum1,dum2,dum3,dum4,tha,thten1,rho,rr,rf)
        endif
      endif

      !----- cvm (if needed) -----!

      IF( eqtset.eq.2 .and. imoist.eq.1 .and. (idiss.eq.1.or.rterm.eq.1) )THEN
        ! for energy-conserving moist thermodynamics:
        ! store cvm in dum1:
        ! store ql  in dum2:
        ! store qi  in dum3:
!$omp parallel do default(shared)  &
!$omp private(i,j,k,n)
        DO k=1,nk
          do j=1,nj
          do i=1,ni
            dum2(i,j,k)=qa(i,j,k,nql1)
          enddo
          enddo
          do n=nql1+1,nql2
            do j=1,nj
            do i=1,ni
              dum2(i,j,k)=dum2(i,j,k)+qa(i,j,k,n)
            enddo
            enddo
          enddo
          IF(iice.eq.1)THEN
            do j=1,nj
            do i=1,ni
              dum3(i,j,k)=qa(i,j,k,nqs1)
            enddo
            enddo
            do n=nqs1+1,nqs2
              do j=1,nj
              do i=1,ni
                dum3(i,j,k)=dum3(i,j,k)+qa(i,j,k,n)
              enddo
              enddo
            enddo
          ELSE
            do j=1,nj
            do i=1,ni
              dum3(i,j,k)=0.0
            enddo
            enddo
          ENDIF
          do j=1,nj
          do i=1,ni
            dum1(i,j,k)=cv+cvv*qa(i,j,k,nqv)+cpl*dum2(i,j,k)+cpi*dum3(i,j,k)
          enddo
          enddo
        ENDDO
      ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k,n)
        DO k=1,nk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k)=cv
          enddo
          enddo
        ENDDO
      ENDIF

      !----- store appropriate rho for budget calculations in dum2 -----!

      IF(axisymm.eq.1)THEN
       ! for axisymmetric grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = rho(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
        enddo
        enddo
        enddo
      ELSE
       ! for Cartesian grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = rho(i,j,k)
        enddo
        enddo
        enddo
      ENDIF

      !-------------------------------------------------------------

      !  budget calculations:
      if(dosfcflx.and.imoist.eq.1)then
        tem0 = dt*dx*dy*dz
!$omp parallel do default(shared)  &
!$omp private(i,j,k,delpi,delth,delqv,delt,n)
        do j=1,nj
        bud2(j) = 0.0d0
        do i=1,ni
          k = 1
          delth = rf0(i,j,1)*rr0(i,j,1)*rdz*mh(i,j,1)*thflux(i,j)
          delqv = rf0(i,j,1)*rr0(i,j,1)*rdz*mh(i,j,1)*qvflux(i,j)
          delpi = rddcv*(pi0(i,j,1)+ppi(i,j,1))*(           &
                                delqv/(eps+qa(i,j,1,nqv))   &
                               +delth/(th0(i,j,1)+tha(i,j,1))  )
          delt = (pi0(i,j,k)+ppi(i,j,k))*delth   &
                +(th0(i,j,k)+tha(i,j,k))*delpi
          bud2(j) = bud2(j) + dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*(        &
                  cv*delt                                                   &
                + cvv*qa(i,j,k,nqv)*delt                                    &
                + cvv*(pi0(i,j,k)+ppi(i,j,k))*(th0(i,j,k)+tha(i,j,k))*delqv &
                + g*zh(i,j,k)*delqv   )
          do n=nql1,nql2
            bud2(j) = bud2(j) + dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*cpl*qa(i,j,k,n)*delt
          enddo
          if(iice.eq.1)then
            do n=nqs1,nqs2
              bud2(j) = bud2(j) + dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*cpi*qa(i,j,k,n)*delt
            enddo
          endif
        enddo
        enddo
        do j=1,nj
          qbudget(9) = qbudget(9) + tem0*bud2(j)
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()
      endif

      !---- Dissipative heating term:

      IF(idiss.eq.1)THEN
        IF( output_dissheat.eq.1 .and. dowrite  )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            tdiag(i,j,k,td_diss) = thten1(i,j,k)
          enddo
          enddo
          enddo
        ENDIF
        ! note:  dissten array stores epsilon (dissipation rate) at w points
        if( bbc.eq.3 )then
          k1 = 2
        else
          k1 = 1
        endif
        if(imoist.eq.1.and.eqtset.eq.2)then
          ! moist, new equations:
!$omp parallel do default(shared)  &
!$omp private(i,j,k,epsd,dheat)
          do k=k1,nk
          do j=1,nj
          do i=1,ni
            epsd = 0.5*(dissten(i,j,k)+dissten(i,j,k+1))
            dheat=epsd/( cpdcv*dum1(i,j,k)*(pi0(i,j,k)+ppi(i,j,k)) )
            thten1(i,j,k)=thten1(i,j,k)+dheat
          enddo
          enddo
          enddo
          if( bbc.eq.3 )then
            k = 1
!$omp parallel do default(shared)  &
!$omp private(i,j,dz1,epsd,dheat)
            do j=1,nj
            do i=1,ni
              dz1 = zf(i,j,2)-zf(i,j,1)
              epsd = (ust(i,j)**3)*alog((dz1+znt(i,j))/znt(i,j))/(karman*dz1)
              dheat=epsd/( cpdcv*dum1(i,j,k)*(pi0(i,j,k)+ppi(i,j,k)) )
              thten1(i,j,k)=thten1(i,j,k)+dheat
            enddo
            enddo
          endif
        else
          ! traditional cloud-modeling equations (also dry equations):
!$omp parallel do default(shared)  &
!$omp private(i,j,k,epsd,dheat)
          do k=k1,nk
          do j=1,nj
          do i=1,ni
            epsd = 0.5*(dissten(i,j,k)+dissten(i,j,k+1))
            dheat=epsd/( cp*(pi0(i,j,k)+ppi(i,j,k)) )
            thten1(i,j,k)=thten1(i,j,k)+dheat
          enddo
          enddo
          enddo
          if( bbc.eq.3 )then
            k = 1
!$omp parallel do default(shared)  &
!$omp private(i,j,dz1,epsd,dheat)
            do j=1,nj
            do i=1,ni
              dz1 = zf(i,j,2)-zf(i,j,1)
              epsd = (ust(i,j)**3)*alog((dz1+znt(i,j))/znt(i,j))/(karman*dz1)
              dheat=epsd/( cp*(pi0(i,j,k)+ppi(i,j,k)) )
              thten1(i,j,k)=thten1(i,j,k)+dheat
            enddo
            enddo
          endif
        endif
        IF( output_dissheat.eq.1 .and. dowrite )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            tdiag(i,j,k,td_diss) = thten1(i,j,k)-tdiag(i,j,k,td_diss)
            if( abs(tdiag(i,j,k,td_diss)).lt.1.0e-14 ) tdiag(i,j,k,td_diss) = 0.0
          enddo
          enddo
          enddo
        ENDIF
      ENDIF

      !---- Rotunno-Emanuel "radiation" term
      !---- (currently capped at 2 K/day ... see RE87 p 546)

      IF(rterm.eq.1)THEN
        tem0 = dt*dx*dy*dz
!$omp parallel do default(shared)  &
!$omp private(i,j,k,thrad,prad)
        do k=1,nk
        bud(k)=0.0d0
        do j=1,nj
        do i=1,ni
          ! NOTE:  thrad is a POTENTIAL TEMPERATURE tendency
          thrad = -tha(i,j,k)/(12.0*3600.0)
          if( tha(i,j,k).gt. 1.0 ) thrad = -1.0/(12.0*3600.0)
          if( tha(i,j,k).lt.-1.0 ) thrad =  1.0/(12.0*3600.0)
          thten1(i,j,k)=thten1(i,j,k)+thrad
          ! associated pressure tendency:
          prad = (pi0(i,j,k)+ppi(i,j,k))*rddcv*thrad/(th0(i,j,k)+tha(i,j,k))
          ! budget:
          bud(k) = bud(k) + dum1(i,j,k)*dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*( &
                            thrad*(pi0(i,j,k)+ppi(i,j,k))    &
                           + prad*(th0(i,j,k)+tha(i,j,k)) )
        enddo
        enddo
        enddo
        do k=1,nk
          qbudget(10) = qbudget(10) + tem0*bud(k)
        enddo
      ENDIF
      if(timestats.ge.1) time_misc=time_misc+mytime()

      if( (dns.eq.1).or.(radopt.ge.1) )then
        ! use thadv to store total potential temperature:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          thadv(i,j,k)=th0(i,j,k)+tha(i,j,k)
        enddo
        enddo
        enddo
      endif

      if(dns.eq.1)then
        call diff2s(2,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                    dum1,dum2,dum3,dum4,thadv,thten1,rho,rr,rf)
      endif

      IF( radopt.ge.1 )THEN
        ! tendency from radiation scheme:
        rdt = 1.0/dt
!$omp parallel do default(shared)  &
!$omp private(i,j,k,tnew,pnew,thnew)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ! cm1r17:  swten and lwten now store TEMPERATURE tendencies:
          ! NOTE:  thadv stores theta (see above)
          tnew = thadv(i,j,k)*(pi0(i,j,k)+ppi(i,j,k)) + dt*(swten(i,j,k)+lwten(i,j,k))
          pnew = rho(i,j,k)*(rd+rv*qa(i,j,k,nqv))*tnew
          thnew = tnew/((pnew*rp00)**rovcp)
          thten1(i,j,k) = thten1(i,j,k) + (thnew-thadv(i,j,k))*rdt
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_rad=time_rad+mytime()
      ENDIF

      IF( ipbl.ge.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          thten1(i,j,k) = thten1(i,j,k) + thpten(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
      ENDIF

      if(iturb.ge.1)then
        ! cm1r18: subtract th0r from theta (as in advection scheme)
        !         (reduces roundoff error)
        IF(.not.terrain_flag)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k,tem)
          do k=1,nk
          tem = th0(1,1,k)-th0r
          do j=0,nj+1
          do i=0,ni+1
            thadv(i,j,k)=tem+tha(i,j,k)
          enddo
          enddo
          enddo
        ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k,tem)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            thadv(i,j,k)=(th0(i,j,k)-th0r)+tha(i,j,k)
          enddo
          enddo
          enddo
        ENDIF
        call turbs(1,dt,dosfcflx,xh,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,thflux,   &
                   rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                   dum1,dum2,dum3,dum4,dum5,sten,rho,rr,rf,thadv,thten1,khh,khv,dum7,dum8)
      endif

!-------------------------------------------------------------------
!  Passive Tracers

      if(iptra.eq.1)then
        do n=1,npt
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ptten(i,j,k,n)=0.0
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_misc=time_misc+mytime()
          if(idiff.eq.1)then
            if(difforder.eq.2)then
              call diff2s(1,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                          dum1,dum2,dum3,dum4,pta(ib,jb,kb,n),ptten(ib,jb,kb,n),rho,rr,rf)
            endif
          endif
          if(iturb.ge.1)then
            call turbs(0,dt,dosfcflx,xh,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,qvflux,   &
                       rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                       dum1,dum2,dum3,dum4,dum5,sten,rho,rr,rf,pta(ib,jb,kb,n),ptten(ib,jb,kb,n),khh,khv,dum7,dum8)
          endif
        enddo
      endif

!-------------------------------------------------------------------
!  Moisture

      if(imoist.eq.1)then
        DO n=1,numq
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,n)=0.0
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_misc=time_misc+mytime()
!---------------------------
          ! qv:
          if(n.eq.nqv)then
            if(idiff.eq.1)then
              if(difforder.eq.2)then
                call diff2s(1,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                            dum1,dum2,dum3,dum4,qa(ib,jb,kb,n),qten(ib,jb,kb,n),rho,rr,rf)
              endif
            endif
            if(iturb.ge.1)then
              call turbs(1,dt,dosfcflx,xh,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,qvflux,   &
                         rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                         dum1,dum2,dum3,dum4,dum5,sten,rho,rr,rf,qa(ib,jb,kb,n),qten(ib,jb,kb,n),khh,khv,dum7,dum8)
            endif
            if(ipbl.ge.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
              do k=1,nk
              do j=1,nj
              do i=1,ni
                qten(i,j,k,nqv) = qten(i,j,k,nqv) + qvpten(i,j,k)
              enddo
              enddo
              enddo
              if(timestats.ge.1) time_pbl=time_pbl+mytime()
            endif
!---------------------------
          ! not qv:
          else
            if(idiff.eq.1)then
              if(difforder.eq.2)then
                call diff2s(1,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                            dum1,dum2,dum3,dum4,qa(ib,jb,kb,n),qten(ib,jb,kb,n),rho,rr,rf)
              endif
            endif
            if(iturb.ge.1)then
              call turbs(0,dt,dosfcflx,xh,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,qvflux,   &
                         rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                         dum1,dum2,dum3,dum4,dum5,sten,rho,rr,rf,qa(ib,jb,kb,n),qten(ib,jb,kb,n),khh,khv,dum7,dum8)
            endif
          endif
!---------------------------
        ENDDO
        IF(ipbl.ge.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if(nqc.ne.0)   &
            qten(i,j,k,nqc) = qten(i,j,k,nqc) + qcpten(i,j,k)
            if(nqi.ne.0)   &
            qten(i,j,k,nqi) = qten(i,j,k,nqi) + qipten(i,j,k)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_pbl=time_pbl+mytime()
        ENDIF
      endif

!-------------------------------------------------------------------
!    NOTE:  now ok to change dum7,dum8
!-------------------------------------------------------------------
!  contribution to pressure tendency from potential temperature:
!  (for mass conservation)
!  plus, some other stuff:

      IF(eqtset.eq.1)THEN
        ! traditional cloud modeling:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ppten(i,j,k)=0.0
        enddo
        enddo
        enddo
      ELSE
        ! mass-conserving pressure eqt:  different sections for moist/dry cases:
        rdt = 1.0/dt
        tem = 0.0001*tsmall
        IF(imoist.eq.1)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k,tnew,pnew,pinew,thnew,qvnew)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            !-----
            ! cm1r17:
            ! note:  nothing in pre-RK section should modify rho
            IF( abs(dt*thten1(i,j,k)).gt.tem .or.  &
                abs(dt*qten(i,j,k,nqv)).gt.qsmall )THEN
              thnew = tha(i,j,k)+dt*thten1(i,j,k)
              qvnew = qa(i,j,k,nqv)+dt*qten(i,j,k,nqv)
              pinew = (rho(i,j,k)*(th0(i,j,k)+thnew)*(rd+rv*qvnew)*rp00)**rddcv - pi0(i,j,k)
              ppten(i,j,k) = (pinew-ppi(i,j,k))*rdt
            ELSE
              ppten(i,j,k) = 0.0
            ENDIF
            !-----
            ! use diabatic tendencies from last timestep as a good estimate:
            ppten(i,j,k)=ppten(i,j,k)+qpten(i,j,k)
            thten1(i,j,k)=thten1(i,j,k)+qtten(i,j,k)
            qten(i,j,k,nqv)=qten(i,j,k,nqv)+qvten(i,j,k)
            qten(i,j,k,nqc)=qten(i,j,k,nqc)+qcten(i,j,k)
          enddo
          enddo
          enddo
        ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k,tnew,pnew,pinew,thnew,qvnew)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            !-----
            ! cm1r17:
            ! note:  nothing in pre-RK section should modify rho
            IF( abs(dt*thten1(i,j,k)).gt.tem )THEN
              thnew = tha(i,j,k)+dt*thten1(i,j,k)
              pinew = (rho(i,j,k)*(th0(i,j,k)+thnew)*rd*rp00)**rddcv - pi0(i,j,k)
              ppten(i,j,k) = (pinew-ppi(i,j,k))*rdt
            ELSE
              ppten(i,j,k)=0.0
            ENDIF
            !-----
          enddo
          enddo
          enddo
        ENDIF  ! endif for moist/dry
      ENDIF    ! endif for eqtset 1/2

        if(timestats.ge.1) time_integ=time_integ+mytime()


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   Begin RK section   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ! time at end of full timestep:
      rtime=sngl(mtime+dt)

!--------------------------------------------------------------------
! RK3 begin

      DO NRK=1,3

        dttmp=dt/float(4-nrk)

!--------------------------------------------------------------------
        IF(nrk.ge.2)THEN
#ifdef MPI
          call comm_3u_end(u3d,uw31,uw32,ue31,ue32,   &
                               us31,us32,un31,un32,reqs_u)
          call comm_3v_end(v3d,vw31,vw32,ve31,ve32,   &
                               vs31,vs32,vn31,vn32,reqs_v)
          call comm_3w_end(w3d,ww31,ww32,we31,we32,   &
                               ws31,ws32,wn31,wn32,reqs_w)
#endif
          if(terrain_flag)then
            call bcwsfc(gz,dzdx,dzdy,u3d,v3d,w3d)
            call bc2d(w3d(ib,jb,1))
          endif
        ENDIF
!--------------------------------------------------------------------
!  Get rru,rrv,rrw,divx
!  (NOTE:  do not change these arrays until after RK steps)

    IF(.not.terrain_flag)THEN
      ! without terrain:

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      DO k=1,nk
        do j=0,nj+1
        do i=0,ni+2
          rru(i,j,k)=rho0(1,1,k)*u3d(i,j,k)
        enddo
        enddo
        do j=0,nj+2
        do i=0,ni+1
          rrv(i,j,k)=rho0(1,1,k)*v3d(i,j,k)
        enddo
        enddo
        IF(k.eq.1)THEN
          do j=0,nj+1
          do i=0,ni+1
            rrw(i,j,   1) = 0.0
            rrw(i,j,nk+1) = 0.0
          enddo
          enddo
        ELSE
          do j=0,nj+1
          do i=0,ni+1
            rrw(i,j,k)=rf0(1,1,k)*w3d(i,j,k)
          enddo
          enddo
        ENDIF
      ENDDO

    ELSE
      ! with terrain:

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      DO k=1,nk
        do j=0,nj+1
        do i=0,ni+2
          rru(i,j,k)=0.5*(rho0(i-1,j,k)+rho0(i,j,k))*u3d(i,j,k)*rgzu(i,j)
        enddo
        enddo
        do j=0,nj+2
        do i=0,ni+1
          rrv(i,j,k)=0.5*(rho0(i,j-1,k)+rho0(i,j,k))*v3d(i,j,k)*rgzv(i,j)
        enddo
        enddo
      ENDDO

!$omp parallel do default(shared)  &
!$omp private(i,j,k,r1,r2)
      DO k=1,nk
        IF(k.eq.1)THEN
          do j=0,nj+1
          do i=0,ni+1
            rrw(i,j,   1) = 0.0
            rrw(i,j,nk+1) = 0.0
          enddo
          enddo
        ELSE
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do j=0,nj+1
          do i=0,ni+1
            rrw(i,j,k)=rf0(i,j,k)*w3d(i,j,k)                                  &
                      +0.5*( ( r2*(rru(i,j,k  )+rru(i+1,j,k  ))               &
                              +r1*(rru(i,j,k-1)+rru(i+1,j,k-1)) )*dzdx(i,j)   &
                            +( r2*(rrv(i,j,k  )+rrv(i,j+1,k  ))               &
                              +r1*(rrv(i,j,k-1)+rrv(i,j+1,k-1)) )*dzdy(i,j)   &
                           )*(sigmaf(k)-zt)*gz(i,j)*rzt
          enddo
          enddo
        ENDIF
      ENDDO

    ENDIF
    if(timestats.ge.1) time_advs=time_advs+mytime()

        IF(terrain_flag)THEN
          call bcw(rrw,0)
#ifdef MPI
          call comm_1w_start(rrw,ww1,ww2,we1,we2,   &
                                 ws1,ws2,wn1,wn2,reqs_w)
          call comm_1w_end(rrw,ww1,ww2,we1,we2,   &
                               ws1,ws2,wn1,wn2,reqs_w)
#endif
        ENDIF

      IF(.not.terrain_flag)THEN
        IF(axisymm.eq.0)THEN
          ! Cartesian without terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            divx(i,j,k)=(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)        &
                       +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)        &
                       +(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k)
            if(abs(divx(i,j,k)).lt.smeps) divx(i,j,k)=0.0
          enddo
          enddo
          enddo
        ELSE
          ! axisymmetric:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            divx(i,j,k)=(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i)   &
                       +(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k)
            if(abs(divx(i,j,k)).lt.smeps) divx(i,j,k)=0.0
          enddo
          enddo
          enddo
        ENDIF
      ELSE
          ! Cartesian with terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            divx(i,j,k)=(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)        &
                       +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)        &
                       +(rrw(i,j,k+1)-rrw(i,j,k))*rdsf(k)
            if(abs(divx(i,j,k)).lt.smeps) divx(i,j,k)=0.0
          enddo
          enddo
          enddo
      ENDIF
      if(timestats.ge.1) time_divx=time_divx+mytime()

!--------------------------------------------------------------------
#ifdef MPI
        IF(nrk.ge.2)THEN
          call comm_1s_end(rho,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
          call getcorner(rho,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
          call bcs2(rho)
        ENDIF
#endif
!--------------------------------------------------------------------
!  U-equation

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          uten(i,j,k)=uten1(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()


        ! Coriolis acceleration:
        if(icor.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
          IF(axisymm.eq.0)THEN
            ! for Cartesian grid:
            if(pertcor.eq.1)then
              do j=1,nj+1
              do i=0,ni+1
                dum1(i,j,k)=v3d(i,j,k)-v0(i,j,k)
              enddo
              enddo
            else
              do j=1,nj+1
              do i=0,ni+1
                dum1(i,j,k)=v3d(i,j,k)
              enddo
              enddo
            endif
            do j=1,nj
            do i=1,ni+1
              uten(i,j,k)=uten(i,j,k)+fcor*             &
               0.25*( (dum1(i  ,j,k)+dum1(i  ,j+1,k))   &
                     +(dum1(i-1,j,k)+dum1(i-1,j+1,k)) )
            enddo
            enddo
          ELSE
            ! for axisymmetric grid:
            do j=1,nj
            do i=2,ni+1
              uten(i,j,k)=uten(i,j,k)+fcor*0.5*(v3d(i,j,k)+v3d(i-1,j,k))
            enddo
            enddo
          ENDIF
        enddo
        if(timestats.ge.1) time_cor=time_cor+mytime()
        endif


        ! inertial term for axisymmetric grid:
        if(axisymm.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
            do i=1,ni+1
              dum1(i,j,k)=(v3d(i,j,k)**2)*rxh(i)
            enddo
            if(ebc.eq.3)then
              dum1(ni+1,j,k) = -dum1(ni,j,k)
            endif
            do i=2,ni+1
              uten(i,j,k)=uten(i,j,k)+0.5*(dum1(i-1,j,k)+dum1(i,j,k))
            enddo
          enddo
          enddo
        endif


          call advu(nrk,arh1,arh2,xf,rxf,arf1,arf2,uf,vh,gz,rgz,gzu,mh,rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,  &
                     rru,u3d,uten,rrv,rrw,rdsf,c1,c2,rho,dttmp)

!--------------------------------------------------------------------
!  V-equation
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          vten(i,j,k)=vten1(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()


        ! Coriolis acceleration:
        ! note for axisymmetric grid: since cm1r18, this term is included in advvaxi
        if( icor.eq.1 .and. axisymm.eq.0 )then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
            ! for Cartesian grid:
            if(pertcor.eq.1)then
              do j=0,nj+1
              do i=1,ni+1
                dum1(i,j,k)=u3d(i,j,k)-u0(i,j,k)
              enddo
              enddo
            else
              do j=0,nj+1
              do i=1,ni+1
                dum1(i,j,k)=u3d(i,j,k)
              enddo
              enddo
            endif
            do j=1,nj+1
            do i=1,ni
              vten(i,j,k)=vten(i,j,k)-fcor*             &
               0.25*( (dum1(i,j  ,k)+dum1(i+1,j  ,k))   &
                     +(dum1(i,j-1,k)+dum1(i+1,j-1,k)) )
            enddo
            enddo
        enddo
        if(timestats.ge.1) time_cor=time_cor+mytime()
        endif


!!!        ! since cm1r17, this term is included in advvaxi
!!!        if(axisymm.eq.1)then
!!!          ! for axisymmetric grid:
!!!
!!!!$omp parallel do default(shared)  &
!!!!$omp private(i,j,k)
!!!          do k=1,nk
!!!          do j=1,nj
!!!          do i=1,ni
!!!            vten(i,j,k)=vten(i,j,k)-(v3d(i,j,k)*rxh(i))*0.5*(xf(i)*u3d(i,j,k)+xf(i+1)*u3d(i+1,j,k))*rxh(i)
!!!          enddo
!!!          enddo
!!!          enddo
!!!
!!!        endif


          call advv(nrk,xh,rxh,arh1,arh2,uh,xf,vf,gz,rgz,gzv,mh,rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,  &
                     rru,rrv,v3d,vten,rrw,rdsf,c1,c2,rho,dttmp)


!--------------------------------------------------------------------
!  finish comms for q/theta:
#ifdef MPI
        IF(nrk.ge.2)THEN
          call comm_3r_end(th3d,pp3d,rw31,rw32,re31,re32,   &
                                     rs31,rs32,rn31,rn32,reqs_p)
          if(imoist.eq.1)then
            call comm_3q_end(q3d,qw31,qw32,qe31,qe32,   &
                                 qs31,qs32,qn31,qn32,reqs_q(1,1))
          endif
        ENDIF
#endif
!--------------------------------------------------------------------
!  Calculate misc. variables
!
!    These arrays store variables that are used later in the
!    SOUND subroutine.  Do not modify t11 or t22 until after sound!
!
!    dum1 = vapor
!    dum2 = all liquid
!    dum3 = all solid
!    t11 = theta_rho
!    t22 = ppterm

        IF(imoist.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,cpm,cvm)
          do k=1,nk

            do j=1,nj
            do i=1,ni
              dum2(i,j,k)=q3d(i,j,k,nql1)
            enddo
            enddo
            do n=nql1+1,nql2
              do j=1,nj
              do i=1,ni
                dum2(i,j,k)=dum2(i,j,k)+q3d(i,j,k,n)
              enddo
              enddo
            enddo
            IF(iice.eq.1)THEN
              do j=1,nj
              do i=1,ni
                dum3(i,j,k)=q3d(i,j,k,nqs1)
              enddo
              enddo
              do n=nqs1+1,nqs2
                do j=1,nj
                do i=1,ni
                  dum3(i,j,k)=dum3(i,j,k)+q3d(i,j,k,n)
                enddo
                enddo
              enddo
            ELSE
              do j=1,nj
              do i=1,ni
                dum3(i,j,k)=0.0
              enddo
              enddo
            ENDIF
            ! save qv,ql,qi for buoyancy calculation:
          IF(eqtset.eq.2)THEN
            do j=1,nj
            do i=1,ni
              t12(i,j,k)=max(q3d(i,j,k,nqv),0.0)
              t13(i,j,k)=max(0.0,dum2(i,j,k)+dum3(i,j,k))
              t11(i,j,k)=(th0(i,j,k)+th3d(i,j,k))*(1.0+reps*t12(i,j,k))     &
                         /(1.0+t12(i,j,k)+t13(i,j,k))
      ! terms in theta and pi equations for proper mass/energy conservation
      ! Reference:  Bryan and Fritsch (2002, MWR), Bryan and Morrison (2012, MWR)
              dum4(i,j,k)=cpl*max(0.0,dum2(i,j,k))+cpi*max(0.0,dum3(i,j,k))
              cpm=cp+cpv*t12(i,j,k)+dum4(i,j,k)
              cvm=1.0/(cv+cvv*t12(i,j,k)+dum4(i,j,k))
              thterm(i,j,k)=(th0(i,j,k)+th3d(i,j,k))*( rd+rv*t12(i,j,k)-rovcp*cpm )*cvm
              t22(i,j,k)=(pi0(i,j,k)+pp3d(i,j,k))*rovcp*cpm*cvm
            enddo
            enddo
          ELSEIF(eqtset.eq.1)THEN
            do j=1,nj
            do i=1,ni
              t12(i,j,k)=max(q3d(i,j,k,nqv),0.0)
              t13(i,j,k)=max(0.0,dum2(i,j,k)+dum3(i,j,k))
              t11(i,j,k)=(th0(i,j,k)+th3d(i,j,k))*(1.0+reps*t12(i,j,k))     &
                         /(1.0+t12(i,j,k)+t13(i,j,k))
              t22(i,j,k)=(pi0(i,j,k)+pp3d(i,j,k))*rddcv
            enddo
            enddo
          ENDIF

          enddo

        ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            t11(i,j,k)=th0(i,j,k)+th3d(i,j,k)
            t22(i,j,k)=(pi0(i,j,k)+pp3d(i,j,k))*rddcv
          enddo
          enddo
          enddo

        ENDIF

        if(timestats.ge.1) time_buoyan=time_buoyan+mytime()

!--------------------------------------------------------------------
        call bcs(t11)
#ifdef MPI
        call comm_1s_start(t11,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
#endif
!--------------------------------------------------------------------
!  W-equation

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          wten(i,j,k)=wten1(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()
 
        if( imoist.eq.1 )then
          ! buoyancy (moisture terms):
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do j=1,nj
            do k=1,nk
            do i=1,ni
              dum2(i,j,k) = repsm1*(t12(i,j,k)-qv0(i,j,k)) - (t13(i,j,k)-qc0(i,j,k)-qi0(i,j,k))
            enddo
            enddo
            do k=2,nk
            do i=1,ni
              wten(i,j,k)=wten(i,j,k)+g*(c1(i,j,k)*dum2(i,j,k-1)+c2(i,j,k)*dum2(i,j,k))
            enddo
            enddo
          enddo
          if(timestats.ge.1) time_buoyan=time_buoyan+mytime()
        endif

        if(psolver.eq.1.or.psolver.eq.4.or.psolver.eq.5)then
          ! buoyancy for non-time-split solvers
          ! (i.e.,truly-compressible/anelastic/incompressible solvers:)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do j=1,nj
            do k=1,nk
            do i=1,ni
              dum2(i,j,k) = th3d(i,j,k)*rth0(i,j,k)
            enddo
            enddo
            do k=2,nk
            do i=1,ni
              wten(i,j,k)=wten(i,j,k)+g*(c1(i,j,k)*dum2(i,j,k-1)+c2(i,j,k)*dum2(i,j,k))
            enddo
            enddo
          enddo
          if(timestats.ge.1) time_buoyan=time_buoyan+mytime()
        endif

          call advw(nrk,xh,rxh,arh1,arh2,uh,xf,vh,gz,rgz,mf,rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,  &
                     rru,rrv,rrw,w3d,wten,rds,c1,c2,rho,dttmp)

!--------------------------------------------------------------------
!  THETA-equation

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem)
      do k=1,nk
        do j=1,nj
        do i=1,ni
          thten(i,j,k)=thten1(i,j,k)
        enddo
        enddo
        IF(.not.terrain_flag)THEN
          tem = th0(1,1,k)-th0r
          do j=jb,je
          do i=ib,ie
            thadv(i,j,k)=tem+th3d(i,j,k)
          enddo
          enddo
        ELSE
          do j=jb,je
          do i=ib,ie
            thadv(i,j,k)=(th0(i,j,k)-th0r)+th3d(i,j,k)
          enddo
          enddo
        ENDIF
      enddo
      if(timestats.ge.1) time_misc=time_misc+mytime()


        weps = 10.0*epsilon
        diffit = 0
        if( idiff.eq.1 .and. difforder.eq.6 ) diffit = 1
        call advs(nrk,1,0,bfoo,xh,rxh,arh1,arh2,uh,ruh,xf,vh,rvh,gz,rgz,mh,rmh, &
                   rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,dum5,dum6,   &
                   rru,rrv,rrw,tha,thadv,thten,0,dttmp,weps,           &
                   flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,rdsf,c1,c2,rho,rr,diffit)

!--------------------------------------------------------------------
!  Pressure equation

      IF( psolver.eq.1 .or. psolver.eq.2 .or. psolver.eq.3 )THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k,tem)
      do k=1,nk
        do j=1,nj
        do i=1,ni
          sten(i,j,k)=ppten(i,j,k)
        enddo
        enddo
        IF(.not.terrain_flag)THEN
          tem = pi0(1,1,k)
          do j=jb,je
          do i=ib,ie
            piadv(i,j,k)=tem+pp3d(i,j,k)
          enddo
          enddo
        ELSE
          do j=jb,je
          do i=ib,ie
            piadv(i,j,k)=pi0(i,j,k)+pp3d(i,j,k)
          enddo
          enddo
        ENDIF
      enddo
      if(timestats.ge.1) time_misc=time_misc+mytime()


          weps = epsilon
          diffit = 0
          call advs(nrk,0,0,bfoo,xh,rxh,arh1,arh2,uh,ruh,xf,vh,rvh,gz,rgz,mh,rmh, &
                     rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,dum5,dum6,   &
                     rru,rrv,rrw,ppi,piadv,sten,0,dttmp,weps,            &
                     flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,rdsf,c1,c2,rho,rr,diffit)

      ENDIF

!--------------------------------------------------------------------
#ifdef MPI
        call comm_1s_end(  t11,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
#endif
!--------------------------------------------------------------------
!  call sound

        get_time_avg = .false.
        IF( psolver.eq.2 .or. psolver.eq.3 .or. psolver.eq.6 )THEN
          IF( imoist.eq.1 .and. numq.ge.1    ) get_time_avg = .true.
          IF( iptra.eq.1 .and. npt.ge.1      ) get_time_avg = .true.
          IF( iturb.eq.1                     ) get_time_avg = .true.
          IF( iprcl.eq.1                     ) get_time_avg = .true.
        ENDIF


        IF(psolver.eq.1)THEN

          call   soundns(xh,rxh,arh1,arh2,uh,xf,uf,yh,vh,yf,vf,           &
                         zh,mh,c1,c2,mf,pi0,thv0,rr0,rf0,                 &
                         rds,sigma,rdsf,sigmaf,                           &
                         zs,gz,rgz,gzu,rgzu,gzv,rgzv,                     &
                         dzdx,dzdy,gx,gxu,gy,gyv,                         &
                         radbcw,radbce,radbcs,radbcn,                     &
                         dum1,dum2,dum3,dum4,                             &
                         u0,ua,u3d,uten,                                  &
                         v0,va,v3d,vten,                                  &
                         wa,w3d,wten,                                     &
                         ppi,pp3d, sten,t11,   t22,dttmp,nrk,rtime,       &
                         th0,tha,th3d,thten,thterm)

        ELSEIF(psolver.eq.2)THEN

          ! NOTE:  rr,rf,prs are used as dummy arrays by sounde

          call   sounde(dt,xh,arh1,arh2,uh,ruh,xf,uf,yh,vh,rvh,yf,vf,     &
                        rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,mf,         &
                        pi0,rho0,rr0,rf0,rrf0,th0,rth0,zs,                &
                        gz,rgz,gzu,rgzu,gzv,rgzv,                         &
                        dzdx,dzdy,gx,gxu,gy,gyv,                          &
                        radbcw,radbce,radbcs,radbcn,                      &
                        dum1,dum2,dum3,dum4,dum5,dum6,                    &
                        dum7,dum8,rr ,rf ,prs,t12,t13,t23,t33,            &
                        u0,rru,us,ua,u3d,uten,                            &
                        v0,rrv,vs,va,v3d,vten,                            &
                        rrw,ws,wa,w3d,wten,                               &
                        ppi,pp3d,piadv,sten ,ppx,                         &
                        tha,th3d,thadv,thten,thterm,                      &
                        t11,t22   ,nrk,dttmp,rtime,get_time_avg,          &
                        pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)

        ELSEIF(psolver.eq.3)THEN

          ! NOTE:  rr,rf,prs are used as dummy arrays by sound

          call   sound( dt,xh,arh1,arh2,uh,ruh,xf,uf,yh,vh,rvh,yf,vf,     &
                        rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,mf,         &
                        pi0,rho0,rr0,rf0,rrf0,th0,rth0,zs,                &
                        gz,rgz,gzu,rgzu,gzv,rgzv,                         &
                        dzdx,dzdy,gx,gxu,gy,gyv,                          &
                        radbcw,radbce,radbcs,radbcn,                      &
                        dum1,dum2,dum3,dum4,dum5,dum6,                    &
                        dum7,dum8,rr ,rf ,prs,t12,t13,t23,t33,            &
                        u0,rru,us,ua,u3d,uten,                            &
                        v0,rrv,vs,va,v3d,vten,                            &
                        rrw,ws,wa,w3d,wten,                               &
                        ppi,pp3d,piadv,sten ,ppx,                         &
                        tha,th3d,thadv,thten,thterm,                      &
                        t11,t22   ,nrk,dttmp,rtime,get_time_avg,          &
                        pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)

        ELSEIF(psolver.eq.4.or.psolver.eq.5)THEN
          ! anelastic/incompressible solver:

          call   anelp(xh,uh,xf,uf,yh,vh,yf,vf,                     &
                       zh,mh,rmh,mf,rmf,pi0,thv0,rho0,prs0,rf0,     &
                       radbcw,radbce,radbcs,radbcn,dum1,divx,       &
                       u0,us,ua,u3d,uten,                           &
                       v0,vs,va,v3d,vten,                           &
                       ws,wa,w3d,wten,                              &
                       ppi,pp3d,tha,th3d,thten,t11,cfb,cfa,cfc,     &
                       ad1,ad2,pdt,deft,rhs,trans,dttmp,nrk,rtime)

        ELSEIF(psolver.eq.6)THEN

          ! NOTE:  rr,rf,prs are used as dummy arrays

          call   soundcb(dt,xh,arh1,arh2,uh,ruh,xf,uf,yh,vh,rvh,yf,vf,    &
                        rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,mf,         &
                        pi0,rho0,rr0,rf0,rrf0,th0,rth0,zs,                &
                        gz,rgz,gzu,rgzu,gzv,rgzv,                         &
                        dzdx,dzdy,gx,gxu,gy,gyv,                          &
                        radbcw,radbce,radbcs,radbcn,                      &
                        dum1,dum2,dum3,dum4,dum5,dum6,                    &
                        dum7,dum8,rr ,rf ,prs,t12,t13,t23,t33,            &
                        u0,rru,us,ua,u3d,uten,                            &
                        v0,rrv,vs,va,v3d,vten,                            &
                        rrw,ws,wa,w3d,wten,                               &
                        ppi,pp3d,piadv,sten ,ppx,                         &
                        tha,th3d,thadv,thten,thterm,                      &
                        t11,t22   ,nrk,dttmp,rtime,get_time_avg,          &
                        pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)

        ENDIF


!--------------------------------------------------------------------
!  re-compute divx if using time-averaged velocities:

    IF( get_time_avg )THEN
      IF(.not.terrain_flag)THEN
        IF(axisymm.eq.0)THEN
          ! Cartesian without terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            divx(i,j,k)=(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)        &
                       +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)        &
                       +(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k)
            if(abs(divx(i,j,k)).lt.smeps) divx(i,j,k)=0.0
          enddo
          enddo
          enddo
        ELSE
          ! axisymmetric:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            divx(i,j,k)=(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i)   &
                       +(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k)
            if(abs(divx(i,j,k)).lt.smeps) divx(i,j,k)=0.0
          enddo
          enddo
          enddo
        ENDIF
      ELSE
          ! Cartesian with terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            divx(i,j,k)=(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)        &
                       +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)        &
                       +(rrw(i,j,k+1)-rrw(i,j,k))*rdsf(k)
            if(abs(divx(i,j,k)).lt.smeps) divx(i,j,k)=0.0
          enddo
          enddo
          enddo
      ENDIF
      if(timestats.ge.1) time_divx=time_divx+mytime()
    ENDIF

!--------------------------------------------------------------------
!  Update v for axisymmetric model simulations:

        IF(axisymm.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            v3d(i,j,k)=va(i,j,k)+dttmp*vten(i,j,k)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_misc=time_misc+mytime()

          call bcv(v3d)

        ENDIF

!--------------------------------------------------------------------

        IF( nrk.eq.3 )THEN
          call calccflquick(dt,uh,vh,mh,u3d,v3d,w3d)
        ENDIF

!--------------------------------------------------------------------
!  radbc

        if(irbc.eq.4)then

          if(ibw.eq.1 .or. ibe.eq.1)then
            call radbcew4(ruf,radbcw,radbce,ua,u3d,dttmp)
          endif

          if(ibs.eq.1 .or. ibn.eq.1)then
            call radbcns4(rvf,radbcs,radbcn,va,v3d,dttmp)
          endif

        endif

!--------------------------------------------------------------------
!  Moisture:

  IF(imoist.eq.1)THEN

    DO n=1,numq

      ! t33 = dummy

      bflag=0
      if(stat_qsrc.eq.1 .and. nrk.eq.3) bflag=1

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=qten(i,j,k,n)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_misc=time_misc+mytime()

      if(nrk.eq.3)then
        pdef = 1
      else
        pdef = 0
      endif


      ! Note: epsilon = 1.0e-18
      weps = 0.01*epsilon
      IF( idm.eq.1 .and. n.ge.nnc1 .and. n .le. nnc2 ) weps = 1.0e5*epsilon
      IF( idmplus.eq.1 .and. n.ge.nzl1 .and. n .le. nzl2 ) weps = 1.d-30/zscale
      diffit = 0
      if( idiff.eq.1 .and. difforder.eq.6 ) diffit = 1
      call advs(nrk,1,bflag,bsq(n),xh,rxh,arh1,arh2,uh,ruh,xf,vh,rvh,gz,rgz,mh,rmh,    &
                 rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,dum5,dum6,   &
                 rru,rrv,rrw,qa(ib,jb,kb,n),q3d(ib,jb,kb,n),sten,pdef,dttmp,weps, &
                 flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,rdsf,c1,c2,rho,rr,diffit)

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        q3d(i,j,k,n)=qa(i,j,k,n)+dttmp*sten(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()

      IF(nrk.lt.3)THEN
        call bcs(q3d(ib,jb,kb,n))
      ENDIF

    ENDDO   ! enddo for n loop

  ENDIF    ! endif for imoist=1

!--------------------------------------------------------------------
!  Get pressure
!  Get density

    pscheck:  IF(psolver.eq.4.or.psolver.eq.5.or.psolver.eq.6)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        prs(i,j,k)=prs0(i,j,k)
        rho(i,j,k)=rho0(i,j,k)
        rr(i,j,k)=rr0(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_prsrho=time_prsrho+mytime()

    ELSE

      IF(imoist.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk

          IF(nrk.eq.3.and.eqtset.eq.2)THEN
          do j=1,nj
          do i=1,ni
            ! subtract-off estimated diabatic terms used during RK steps:
            pp3d(i,j,k)=pp3d(i,j,k)-dt*qpten(i,j,k)
            th3d(i,j,k)=th3d(i,j,k)-dt*qtten(i,j,k)
            q3d(i,j,k,nqv)=q3d(i,j,k,nqv)-dt*qvten(i,j,k)
            q3d(i,j,k,nqc)=q3d(i,j,k,nqc)-dt*qcten(i,j,k)
            ! save values before calculating microphysics:
            qpten(i,j,k)=pp3d(i,j,k)
            qtten(i,j,k)=th3d(i,j,k)
            qvten(i,j,k)=q3d(i,j,k,nqv)
            qcten(i,j,k)=q3d(i,j,k,nqc)
          enddo
          enddo
          ENDIF

          do j=1,nj
          do i=1,ni
            prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
            rho(i,j,k)=prs(i,j,k)                         &
               /( (th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))     &
                 *(rd+max(0.0,q3d(i,j,k,nqv))*rv) )
            rr(i,j,k) = 1.0/rho(i,j,k)
          enddo
          enddo

        enddo

      ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
          rho(i,j,k)=prs(i,j,k)   &
             /(rd*(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k)))
          rr(i,j,k) = 1.0/rho(i,j,k)
        enddo
        enddo
        enddo

      ENDIF


      !-----------------------------------------------
      pmod:  IF( apmasscon.eq.1 .and. nrk.eq.3 )THEN
        ! cm1r18:  adjust average pressure perturbation to ensure 
        !          conservation of total dry-air mass

        dumk1 = 0.0
        dumk2 = 0.0

        IF( axisymm.eq.0 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dumk1(k) = dumk1(k) + rho(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)
            dumk2(k) = dumk2(k) + (pi0(i,j,k)+pp3d(i,j,k))
          enddo
          enddo
          enddo
        ELSEIF( axisymm.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dumk1(k) = dumk1(k) + rho(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)
            dumk2(k) = dumk2(k) + (pi0(i,j,k)+pp3d(i,j,k))
          enddo
          enddo
          enddo
        ENDIF

        mass2 = 0.0
        p2 = 0.0

        do k=1,nk
          mass2 = mass2 + dumk1(k) 
          p2 = p2 + dumk2(k) 
        enddo

#ifdef MPI
        p0=0.0d0
        call MPI_ALLREDUCE(mass2,p0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        mass2 = p0*(dx*dy*dz)
        !---
        p0=0.0d0
        call MPI_ALLREDUCE(p2,p0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        p2 = p0/dble(nx*ny*nz)
#else
        mass2 = mass2*(dx*dy*dz)
        p2 = p2/dble(nx*ny*nz)
#endif

        tem = ( (mass1/mass2)**(dble(rd)/dble(cv)) - 1.0d0 )*p2

        IF( imoist.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            pp3d(i,j,k) = pp3d(i,j,k) + tem
            prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
            rho(i,j,k)=prs(i,j,k)                         &
               /( (th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))     &
                 *(rd+max(0.0,q3d(i,j,k,nqv))*rv) )
            rr(i,j,k) = 1.0/rho(i,j,k)
          enddo
          enddo
          enddo
        ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            pp3d(i,j,k) = pp3d(i,j,k) + tem
            prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
            rho(i,j,k)=prs(i,j,k)   &
               /(rd*(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k)))
            rr(i,j,k) = 1.0/rho(i,j,k)
          enddo
          enddo
          enddo
        ENDIF

      ENDIF  pmod
      !-----------------------------------------------

      if(timestats.ge.1) time_prsrho=time_prsrho+mytime()

    ENDIF  pscheck

      call bcs(rho)

!--------------------------------------------------------------------
#ifdef MPI
      call comm_3u_start(u3d,uw31,uw32,ue31,ue32,   &
                             us31,us32,un31,un32,reqs_u)
      call comm_3v_start(v3d,vw31,vw32,ve31,ve32,   &
                             vs31,vs32,vn31,vn32,reqs_v)
      call comm_3w_start(w3d,ww31,ww32,we31,we32,   &
                             ws31,ws32,wn31,wn32,reqs_w)
      call comm_1s_start(rho,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
      if(nrk.lt.3)then
        call comm_3r_start(th3d,pp3d,rw31,rw32,re31,re32,   &
                                     rs31,rs32,rn31,rn32,reqs_p)
      endif
      IF(imoist.eq.1)THEN
        !  start comms for q:
        ! dont communicate on last rk step
        if(nrk.lt.3)then
          call comm_3q_start(q3d,qw31,qw32,qe31,qe32,   &
                                 qs31,qs32,qn31,qn32,reqs_q(1,1))
        endif
      ENDIF
#endif
!--------------------------------------------------------------------
!  TKE advection
 
        IF(iturb.eq.1)THEN

          ! use wten for tke tendency, step tke forward:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            wten(i,j,k)=tketen(i,j,k)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_misc=time_misc+mytime()

#ifdef MPI
        IF(nrk.ge.2)THEN
          call comm_3t_end(tke3d,tkw1,tkw2,tke1,tke2,   &
                                 tks1,tks2,tkn1,tkn2,reqs_tk)
        ENDIF
#endif

            call advw(nrk,xh,rxh,arh1,arh2,uh,xf,vh,gz,rgz,mf,rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,  &
                       rru,rrv,rrw,tke3d,wten,rds,c1,c2,rho,dttmp)

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
          do j=1,nj
          do i=1,ni
            tke3d(i,j,k)=tkea(i,j,k)+dttmp*wten(i,j,k)
            if(tke3d(i,j,k).lt.1.0e-6) tke3d(i,j,k)=0.0
          enddo
          enddo
        enddo
        if(timestats.ge.1) time_integ=time_integ+mytime()


          call bcw(tke3d,1)
#ifdef MPI
          call comm_3t_start(tke3d,tkw1,tkw2,tke1,tke2,   &
                                   tks1,tks2,tkn1,tkn2,reqs_tk)
#endif

        ENDIF

!--------------------------------------------------------------------
!  Passive Tracers

    if(iptra.eq.1)then

      if( nrk.eq.3 .and. pdtra.eq.1 )then
        pdef = 1
      else
        pdef = 0
      endif

    DO n=1,npt

      ! t33 = dummy

      bflag=0
      if(stat_qsrc.eq.1 .and. nrk.eq.3) bflag=1

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=ptten(i,j,k,n)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_misc=time_misc+mytime()


#ifdef MPI
          IF(nrk.ge.2)THEN
            call comm_3s_end(pt3d(ib,jb,kb,n),                           &
                  tw1(1,1,1,n),tw2(1,1,1,n),te1(1,1,1,n),te2(1,1,1,n),   &
                  ts1(1,1,1,n),ts2(1,1,1,n),tn1(1,1,1,n),tn2(1,1,1,n),   &
                  reqs_t(1,n))
          ENDIF
#endif

      weps = 1.0*epsilon
      diffit = 0
      if( idiff.eq.1 .and. difforder.eq.6 ) diffit = 1
      call advs(nrk,1,bflag,bfoo,xh,rxh,arh1,arh2,uh,ruh,xf,vh,rvh,gz,rgz,mh,rmh,        &
                 rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,dum5,dum6,     &
                 rru,rrv,rrw,pta(ib,jb,kb,n),pt3d(ib,jb,kb,n),sten,pdef,dttmp,weps, &
                 flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,rdsf,c1,c2,rho,rr,diffit)

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        pt3d(i,j,k,n)=pta(i,j,k,n)+dttmp*sten(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()

      IF(nrk.le.2)THEN
        call bcs(pt3d(ib,jb,kb,n))
#ifdef MPI
        call comm_3s_start(pt3d(ib,jb,kb,n)   &
                     ,tw1(1,1,1,n),tw2(1,1,1,n),te1(1,1,1,n),te2(1,1,1,n)     &
                     ,ts1(1,1,1,n),ts2(1,1,1,n),tn1(1,1,1,n),tn2(1,1,1,n)     &
                     ,reqs_t(1,n) )
#endif
      ENDIF

    ENDDO
    endif

!--------------------------------------------------------------------
! RK loop end

      ENDDO


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   End of RK section   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!--------------------------------------------------------------------
!  Final step for Passive Tracers
!  (using final value of rho)

    if(iptra.eq.1)then
      DO n=1,npt
        if( pdtra.eq.1 ) call pdefq(0.0,afoo,ruh,rvh,rmh,rho,pt3d(ib,jb,kb,n))
        call bcs(pt3d(ib,jb,kb,n))
#ifdef MPI
        call comm_3s_start(pt3d(ib,jb,kb,n)   &
                     ,tw1(1,1,1,n),tw2(1,1,1,n),te1(1,1,1,n),te2(1,1,1,n)     &
                     ,ts1(1,1,1,n),ts2(1,1,1,n),tn1(1,1,1,n),tn2(1,1,1,n)     &
                     ,reqs_t(1,n) )
#endif
      ENDDO
    endif


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   BEGIN microphysics   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF(imoist.eq.1)THEN

        getdbz = .false.
        IF(output_dbz.eq.1)THEN
          if( ((mtime+dbldt).ge.(taptim-0.1*dt)) .or. stopit )then
            getdbz = .true.
          endif
          if( ((mtime+dbldt).ge.(rsttim-0.1*dt)) .and. rstfrq.gt.0.0001 )then
            getdbz = .true.
          endif
          if( iprcl.eq.1 )then
            if( ((mtime+dbldt).ge.(prcltim-0.1*dt)) .or. prclfrq.le.0.0 )then
              getdbz = .true.
            endif
          endif
          if(getdbz)then
            if(dowr) write(outfile,*) '  Getting dbz ... '
          endif
        ENDIF

        getvt = .false.
        IF( efall.eq.1 ) getvt = .true.
        IF( dowrite .and. output_fallvel.eq.1 ) getvt = .true.

        ! sten = dbz
        ! dum1 = T
        ! dum3 = appropriate rho for budget calculations
        ! store copy of T in thten array:

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n)
        DO k=1,nk

          do j=1,nj
          do i=1,ni
            dum1(i,j,k)=(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
            thten(i,j,k)=dum1(i,j,k)
            qten(i,j,k,nqv)=q3d(i,j,k,nqv)
          enddo
          enddo

          if(getdbz)then
            ! store dbz in sten array:
            do j=1,nj
            do i=1,ni
              sten(i,j,k)=0.0
            enddo
            enddo
          endif

          IF(axisymm.eq.0)THEN
            ! for Cartesian grid:
            do j=1,nj
            do i=1,ni
              dum3(i,j,k)=rho(i,j,k)
            enddo
            enddo
          ELSE
            ! for axisymmetric grid:
            do j=1,nj
            do i=1,ni
              dum3(i,j,k) = rho(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
            enddo
            enddo
          ENDIF

          IF( output_mptend.eq.1 .and. dowrite  )THEN
            do j=1,nj
            do i=1,ni
              tdiag(i,j,k,td_mptend) = th3d(i,j,k)
            enddo
            enddo
          ENDIF

        ENDDO


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  NOTES:
!           sten       is used for     dbz
!
!           dum1   is   T
!           dum3   is   rho for budget calculations
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Kessler scheme   cccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        IF(ptype.eq.1)THEN
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq( qsmall,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call k_fallout(rho,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call geterain(dt,cpl,lv1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          if(efall.ge.1)then
            call getcvm(dum2,q3d)
            call getefall(1,cpl,mf,dum1,dum2,dum4,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk-1
            do j=1,nj
            do i=1,ni
              if( abs(dt*dum4(i,j,k)).ge.tsmall )then
                dum1(i,j,k) = dum1(i,j,k) + dt*dum4(i,j,k)
                prs(i,j,k)=rho(i,j,k)*rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps)
                pp3d(i,j,k)=(prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
              endif
            enddo
            enddo
            enddo
          endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call kessler(dt,qbudget(3),qbudget(4),qbudget(5),ruh,rvh,rmh,pi0,th0,dum1,   &
                       rho,dum3,pp3d,th3d,prs,                            &
                       q3d(ib,jb,kb,nqv),q3d(ib,jb,kb,2),q3d(ib,jb,kb,3))
          call satadj(4,dt,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,   &
                      rho,dum3,pp3d,prs,th3d,q3d)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Goddard LFO scheme   cccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ELSEIF(ptype.eq.2)THEN
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq( qsmall,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq( qsmall,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq( qsmall,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq( qsmall,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
          call goddard(dt,qbudget(3),qbudget(4),qbudget(5),ruh,rvh,rmh,pi0,th0,             &
                       rho,dum3,prs,pp3d,th3d,                            &
     q3d(ib,jb,kb,1), q3d(ib,jb,kb,2),q3d(ib,jb,kb,3),qten(ib,jb,kb,3),   &
     q3d(ib,jb,kb,4),qten(ib,jb,kb,4),q3d(ib,jb,kb,5),qten(ib,jb,kb,5),   &
     q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
          call satadj_ice(4,dt,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,     &
                          rho,dum3,pp3d,prs,th3d,                     &
              q3d(ib,jb,kb,1),q3d(ib,jb,kb,2),q3d(ib,jb,kb,3),   &
              q3d(ib,jb,kb,4),q3d(ib,jb,kb,5),q3d(ib,jb,kb,6))
          call geterain(dt,cpl,lv1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call geterain(dt,cpi,ls1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,4),qten(ib,jb,kb,4))
          call geterain(dt,cpi,ls1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,5),qten(ib,jb,kb,5))
          call geterain(dt,cpi,ls1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
          if(efall.ge.1)then
            call getcvm(dum2,q3d)
            call getefall(1,cpl,mf,dum1,dum2,dum4,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
            call getefall(0,cpi,mf,dum1,dum2,dum4,q3d(ib,jb,kb,4),qten(ib,jb,kb,4))
            call getefall(0,cpi,mf,dum1,dum2,dum4,q3d(ib,jb,kb,5),qten(ib,jb,kb,5))
            call getefall(0,cpi,mf,dum1,dum2,dum4,q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk-1
            do j=1,nj
            do i=1,ni
              if( abs(dt*dum4(i,j,k)).ge.tsmall )then
                dum1(i,j,k) = dum1(i,j,k) + dt*dum4(i,j,k)
                prs(i,j,k)=rho(i,j,k)*rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps)
                pp3d(i,j,k)=(prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
              endif
            enddo
            enddo
            enddo
          endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,4),qten(ib,jb,kb,4))
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,5),qten(ib,jb,kb,5))
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
          if(getdbz) call calcdbz(rho,q3d(ib,jb,kb,3),q3d(ib,jb,kb,5),q3d(ib,jb,kb,6),sten)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Thompson scheme   ccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ELSEIF(ptype.eq.3)THEN
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq( qsmall,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq( qsmall,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq( qsmall,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq( qsmall,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
!!!          call pdefq(    1.0,asq(7),ruh,rvh,rmh,rho,q3d(ib,jb,kb,7))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! cm1r17:  to make things easier to understand, use same arrays 
            !          that are used for morrison code:
            ! dum1 = T  (this should have been calculated already)
            ! dum2 = pi (nondimensional pressure)
            ! dum4 = dz
            ! thten = copy of T  (this should have been calculated already)
            dum2(i,j,k)=pi0(i,j,k)+pp3d(i,j,k)
            dum4(i,j,k)=dz*rmh(i,j,k)
          enddo
          enddo
          enddo
          call mp_gt_driver(q3d(ib,jb,kb,1),q3d(ib,jb,kb,2),q3d(ib,jb,kb,3), &
                            q3d(ib,jb,kb,4),q3d(ib,jb,kb,5),q3d(ib,jb,kb,6), &
                            q3d(ib,jb,kb,7),q3d(ib,jb,kb,8),                 &
                            th0,dum1,dum2,prs,dum4,dt,rain,                  &
                            qbudget(1),qbudget(2),qbudget(5),qbudget(6),     &
                            ruh,rvh,rmh,rho,dum3,sten,getdbz)
          ! Get final values for th3d,pp3d,prs:
          ! Note:  dum1 stores temperature, thten stores old temperature:
          IF( eqtset.eq.2 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              if( abs(dum1(i,j,k)-thten(i,j,k)).ge.tsmall .or.  &
                  abs(q3d(i,j,k,nqv)-qten(i,j,k,nqv)).ge.qsmall )then
                prs(i,j,k)=rho(i,j,k)*(rd+rv*q3d(i,j,k,nqv))*dum1(i,j,k)
                pp3d(i,j,k)=(prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
              endif
            enddo
            enddo
            enddo
          ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              if( abs(dum1(i,j,k)-thten(i,j,k)).ge.tsmall .or.  &
                  abs(q3d(i,j,k,nqv)-qten(i,j,k,nqv)).ge.qsmall )then
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
                rho(i,j,k)=prs(i,j,k)/(rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps))
              endif
            enddo
            enddo
            enddo
          ENDIF
          if(timestats.ge.1) time_microphy=time_microphy+mytime()

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   GSR LFO scheme   cccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.4)THEN
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq( qsmall,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq( qsmall,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq( qsmall,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq( qsmall,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
          call lfo_ice_drive(dt, mf, pi0, prs0, pp3d, prs, th0, th3d,    &
                             qv0, rho0, q3d, qten, dum1)
          do n=2,numq
            call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                         q3d(ib,jb,kb,n),qten(ib,jb,kb,n))
          enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Morrison scheme   cccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.5)THEN
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq( qsmall,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq( qsmall,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq( qsmall,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq( qsmall,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
!!!          call pdefq(    1.0,asq(7),ruh,rvh,rmh,rho,q3d(ib,jb,kb,7))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! dum1 = T  (this should have been calculated already)
            ! dum4 = dz
            ! thten = copy of T  (this should have been calculated already)
            dum4(i,j,k)=dz*rmh(i,j,k)
          enddo
          enddo
          enddo
          IF(numq.eq.11)THEN
            ! ppten stores ncc:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              ppten(i,j,k) = q3d(i,j,k,11)
            enddo
            enddo
            enddo
          ENDIF
          ! cm1r17:  get fall velocities (store in qten array)
          call MP_GRAUPEL(nstep,dum1,                                 &
                          q3d(ib,jb,kb, 1),q3d(ib,jb,kb, 2),q3d(ib,jb,kb, 3), &
                          q3d(ib,jb,kb, 4),q3d(ib,jb,kb, 5),q3d(ib,jb,kb, 6), &
                          q3d(ib,jb,kb, 7),q3d(ib,jb,kb, 8),q3d(ib,jb,kb, 9), &
                          q3d(ib,jb,kb,10),ppten,                             &
                               prs,rho,dt,dum4,w3d,rain,                      &
                          effc,effi,effs,effr,effg,effis,                     &
                          qbudget(1),qbudget(2),qbudget(5),qbudget(6),        &
                          ruh,rvh,rmh,dum3,sten,getdbz,                       &
                          qten(ib,jb,kb,nqc),qten(ib,jb,kb,nqr),qten(ib,jb,kb,nqi),  &
                          qten(ib,jb,kb,nqs),qten(ib,jb,kb,nqg),getvt)
          IF(numq.eq.11)THEN
            ! ppten stores ncc:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              q3d(i,j,k,11) = ppten(i,j,k)
            enddo
            enddo
            enddo
          ENDIF
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
          IF(efall.eq.1)THEN
            ! dum1 = T
            ! dum2 = cvm
            ! dum4 = T tendency
            call getcvm(dum2,q3d)
            call getefall(1,cpl,mf,dum1,dum2,dum4,q3d(ib,jb,kb,nqc),qten(ib,jb,kb,nqc))
            call getefall(0,cpl,mf,dum1,dum2,dum4,q3d(ib,jb,kb,nqr),qten(ib,jb,kb,nqr))
            call getefall(0,cpi,mf,dum1,dum2,dum4,q3d(ib,jb,kb,nqi),qten(ib,jb,kb,nqi))
            call getefall(0,cpi,mf,dum1,dum2,dum4,q3d(ib,jb,kb,nqs),qten(ib,jb,kb,nqs))
            call getefall(0,cpi,mf,dum1,dum2,dum4,q3d(ib,jb,kb,nqg),qten(ib,jb,kb,nqg))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk-1
            do j=1,nj
            do i=1,ni
              dum1(i,j,k) = dum1(i,j,k) + dt*dum4(i,j,k)
            enddo
            enddo
            enddo
          ENDIF
          ! Get final values for th3d,pp3d,prs:
          ! Note:  dum1 stores temperature, thten stores old temperature:
          IF( eqtset.eq.2 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              if( abs(dum1(i,j,k)-thten(i,j,k)).ge.tsmall .or.  &
                  abs(q3d(i,j,k,nqv)-qten(i,j,k,nqv)).ge.qsmall )then
                prs(i,j,k)=rho(i,j,k)*(rd+rv*q3d(i,j,k,nqv))*dum1(i,j,k)
                pp3d(i,j,k)=(prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
              endif
            enddo
            enddo
            enddo
          ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              if( abs(dum1(i,j,k)-thten(i,j,k)).ge.tsmall .or.  &
                  abs(q3d(i,j,k,nqv)-qten(i,j,k,nqv)).ge.qsmall )then
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
                rho(i,j,k)=prs(i,j,k)/(rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps))
              endif
            enddo
            enddo
            enddo
          ENDIF
          IF( output_fallvel.eq.1 .and. dowrite  )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              qdiag(i,j,k,qd_vtc) = qten(i,j,k,nqc)
              qdiag(i,j,k,qd_vtr) = qten(i,j,k,nqr)
              qdiag(i,j,k,qd_vts) = qten(i,j,k,nqs)
              qdiag(i,j,k,qd_vtg) = qten(i,j,k,nqg)
              qdiag(i,j,k,qd_vti) = qten(i,j,k,nqi)
            enddo
            enddo
            enddo
          ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   RE87 scheme   ccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.6)THEN
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if(q3d(i,j,k,2).gt.0.001)then
              qten(i,j,k,2) = v_t
            else
              qten(i,j,k,2) = 0.0
            endif
          enddo
          enddo
          enddo
          call geterain(dt,cpl,lv1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,2),qten(ib,jb,kb,2))
          if(efall.ge.1)then
            call getcvm(dum2,q3d)
            call getefall(1,cpl,mf,dum1,dum2,dum4,q3d(ib,jb,kb,2),qten(ib,jb,kb,2))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk-1
            do j=1,nj
            do i=1,ni
              if( abs(dt*dum4(i,j,k)).ge.tsmall )then
                dum1(i,j,k) = dum1(i,j,k) + dt*dum4(i,j,k)
                prs(i,j,k)=rho(i,j,k)*rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps)
                pp3d(i,j,k)=(prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
              endif
            enddo
            enddo
            enddo
          endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,2),qten(ib,jb,kb,2))
          call satadj(4,dt,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,   &
                      rho,dum3,pp3d,prs,th3d,q3d)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Ziegler/Mansell (NSSL) two-moment scheme
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        ELSEIF(ptype.ge.26)THEN
          IF ( ptype .eq. 26 ) THEN
            j = 13
          ELSEIF ( ptype .eq. 27 ) THEN
            j = 16
          ELSEIF ( ptype .eq. 28) THEN ! single moment
            j = 6
          ELSEIF ( ptype .eq. 29 ) THEN
            j = 19
          ENDIF
          DO i = 1,j
            call pdefq(0.0,asq(i),ruh,rvh,rmh,rho,q3d(ib,jb,kb,i))
          ENDDO

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = pi0(i,j,k)+pp3d(i,j,k)
            dum2(i,j,k) = dz*rmh(i,j,k)
            dum4(i,j,k) = th0(i,j,k)+th3d(i,j,k)
            ! store old theta in thten array:
            thten(i,j,k)=dum4(i,j,k)
          enddo
          enddo
          enddo
          
          IF ( ptype .eq. 26 ) THEN  ! graupel only
             call nssl_2mom_driver(                          &
                               th  = dum4,                   &
                               qv  = q3d(ib,jb,kb, 1),       &
                               qc  = q3d(ib,jb,kb, 2),       &
                               qr  = q3d(ib,jb,kb, 3),       &
                               qi  = q3d(ib,jb,kb, 4),       &
                               qs  = q3d(ib,jb,kb, 5),       &
                               qh  = q3d(ib,jb,kb, 6),       &
                               cn  = q3d(ib,jb,kb, 7),       &
                               ccw = q3d(ib,jb,kb, 8),       &
                               crw = q3d(ib,jb,kb, 9),       &
                               cci = q3d(ib,jb,kb, 10),      &
                               csw = q3d(ib,jb,kb, 11),      &
                               chw = q3d(ib,jb,kb, 12),      &
                               vhw = q3d(ib,jb,kb, 13),      &
                               pii = dum1,                   &
                               p   =  prs,                   &
                               w   =  w3d,                   &
                               dn  =  rho,                   &
                               dz  =  dum2,                  &
                               dtp = dt,                     &
                               itimestep = nstep,            &
                              RAIN = rain,                   &
                              nrain = nrain,                 &
                              dbz = sten,                    &
                              ruh = ruh, rvh = rvh, rmh = rmh, &
                              dx = dx, dy = dy,              &
                              tcond = qbudget(1),            &
                              tevac = qbudget(2),            &
                              tevar = qbudget(5),            &
                              train = qbudget(6),            &
                              rr    = dum3,                  &
                              diagflag = getdbz,                  &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
         ELSEIF ( ptype .eq. 27 ) THEN
             call nssl_2mom_driver(                          &
                               th  = dum4,                   &
                               qv  = q3d(ib,jb,kb, 1),       &
                               qc  = q3d(ib,jb,kb, 2),       &
                               qr  = q3d(ib,jb,kb, 3),       &
                               qi  = q3d(ib,jb,kb, 4),       &
                               qs  = q3d(ib,jb,kb, 5),       &
                               qh  = q3d(ib,jb,kb, 6),       &
                               qhl = q3d(ib,jb,kb, 7),       &
                               cn  = q3d(ib,jb,kb, 8),       &
                               ccw = q3d(ib,jb,kb, 9),       &
                               crw = q3d(ib,jb,kb,10),       &
                               cci = q3d(ib,jb,kb, 11),      &
                               csw = q3d(ib,jb,kb, 12),      &
                               chw = q3d(ib,jb,kb, 13),      &
                               chl = q3d(ib,jb,kb, 14),      &
                               vhw = q3d(ib,jb,kb, 15),      &
                               vhl = q3d(ib,jb,kb, 16),      &
                               pii = dum1,                   &
                               p   =  prs,                   &
                               w   =  w3d,                   &
                               dn  =  rho,                   &
                               dz  =  dum2,                  &
                               dtp = dt,                     &
                               itimestep = nstep,            &
                              RAIN = rain,                   &
                              nrain = nrain,                 &
                              dbz = sten,                    &
                              ruh = ruh, rvh = rvh, rmh = rmh, &
                              dx = dx, dy = dy,              &
                              tcond = qbudget(1),            &
                              tevac = qbudget(2),            &
                              tevar = qbudget(5),            &
                              train = qbudget(6),            &
                              rr    = dum3,                  &
                              diagflag = getdbz,             &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
          ELSEIF ( ptype .eq. 28 ) THEN  ! single moment
             call nssl_2mom_driver(                          &
                               th  = dum4,                   &
                               qv  = q3d(ib,jb,kb, 1),       &
                               qc  = q3d(ib,jb,kb, 2),       &
                               qr  = q3d(ib,jb,kb, 3),       &
                               qi  = q3d(ib,jb,kb, 4),       &
                               qs  = q3d(ib,jb,kb, 5),       &
                               qh  = q3d(ib,jb,kb, 6),       &
                               pii = dum1,                   &
                               p   =  prs,                   &
                               w   =  w3d,                   &
                               dn  =  rho,                   &
                               dz  =  dum2,                  &
                               dtp = dt,                     &
                               itimestep = nstep,            &
                              RAIN = rain,                   &
                              nrain = nrain,                 &
                              dbz = sten,                    &
                              ruh = ruh, rvh = rvh, rmh = rmh, &
                              dx = dx, dy = dy,              &
                              tcond = qbudget(1),            &
                              tevac = qbudget(2),            &
                              tevar = qbudget(5),            &
                              train = qbudget(6),            &
                              rr    = dum3,                  &
                              diagflag = getdbz,                  &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)

!         ELSEIF ( ptype .eq. 29 ) THEN ! 3-moment
!             call nssl_2mom_driver(                          &
!                               th  = dum4,                   &
!                               qv  = q3d(ib,jb,kb, 1),       &
!                               qc  = q3d(ib,jb,kb, 2),       &
!                               qr  = q3d(ib,jb,kb, 3),       &
!                               qi  = q3d(ib,jb,kb, 4),       &
!                               qs  = q3d(ib,jb,kb, 5),       &
!                               qh  = q3d(ib,jb,kb, 6),       &
!                               qhl = q3d(ib,jb,kb, 7),       &
!                               cn  = q3d(ib,jb,kb, 8),       &
!                               ccw = q3d(ib,jb,kb, 9),       &
!                               crw = q3d(ib,jb,kb,10),       &
!                               cci = q3d(ib,jb,kb, 11),      &
!                               csw = q3d(ib,jb,kb, 12),      &
!                               chw = q3d(ib,jb,kb, 13),      &
!                               chl = q3d(ib,jb,kb, 14),      &
!                               zrw = q3d(ib,jb,kb, 15),      &
!                               zhw = q3d(ib,jb,kb, 16),      &
!                               zhl = q3d(ib,jb,kb, 17),      &
!                               vhw = q3d(ib,jb,kb, 18),      &
!                               vhl = q3d(ib,jb,kb, 19),      &
!                               pii = dum1,                   &
!                               p   =  prs,                   &
!                               w   =  w3d,                   &
!                               dn  =  rho,                   &
!                               dz  =  dum2,                  &
!                               dtp = dt,                     &
!                               itimestep = nstep,            &
!                              RAIN = rain,                   &
!                              nrain = nrain,                 &
!                              dbz = sten,                    &
!                              ruh = ruh, rvh = rvh, rmh = rmh, &
!                              dx = dx, dy = dy,              &
!                              tcond = qbudget(1),            &
!                              tevac = qbudget(2),            &
!                              tevar = qbudget(5),            &
!                              train = qbudget(6),            &
!                              rr    = dum3,                  &
!                              diagflag = getdbz,             &
!                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
!                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
!          
          ENDIF

        IF(eqtset.eq.2)THEN
          ! for mass conservation:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( abs(dum4(i,j,k)-thten(i,j,k)).ge.tsmall .or.  &
                abs(q3d(i,j,k,nqv)-qten(i,j,k,nqv)).ge.qsmall )then
              prs(i,j,k) = rho(i,j,k)*rd*dum4(i,j,k)*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps)
              pp3d(i,j,k) = (prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
              th3d(i,j,k) = dum4(i,j,k)*dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
            endif
          enddo
          enddo
          enddo
        ELSE
          ! traditional thermodynamics:  p,pi remain unchanged
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( abs(dum4(i,j,k)-thten(i,j,k)).ge.tsmall .or.  &
                abs(q3d(i,j,k,nqv)-qten(i,j,k,nqv)).ge.qsmall )then
              th3d(i,j,k)= dum4(i,j,k) - th0(i,j,k)
              rho(i,j,k)=prs(i,j,k)/(rd*dum4(i,j,k)*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps))
            endif
          enddo
          enddo
          enddo
        ENDIF

          if(timestats.ge.1) time_microphy=time_microphy+mytime()

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  insert new microphysics schemes here
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        ELSEIF(ptype.eq.8)THEN
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! otherwise, stop for undefined ptype
        ELSE
          print *,'  Undefined ptype!'
          call stopcm1
        ENDIF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   END microphysics   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          IF( output_mptend.eq.1 .and. dowrite  )THEN
            rdt = 1.0/dt
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              tdiag(i,j,k,td_mptend) = (th3d(i,j,k)-tdiag(i,j,k,td_mptend))*rdt
            enddo
            enddo
            enddo
          ENDIF

        if(timestats.ge.1) time_microphy=time_microphy+mytime()
      ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Begin:  message passing

          call bcs(th3d)
          call bcs(pp3d)
#ifdef MPI
          call comm_3r_start(th3d,pp3d,rw31,rw32,re31,re32,   &
                                       rs31,rs32,rn31,rn32,reqs_p)
#endif

      IF( imoist.eq.1 )THEN
          DO n=1,numq
            call bcs(q3d(ib,jb,kb,n))
          ENDDO
#ifdef MPI
          call comm_3q_start(q3d,qw31,qw32,qe31,qe32,   &
                                 qs31,qs32,qn31,qn32,reqs_q(1,1))
#endif
      ENDIF

!Done:  message passing
!-----------------------------------------------------------------
!  cm1r17:  diabatic tendencies for next timestep:

        IF(imoist.eq.1.and.eqtset.eq.2)THEN
          ! get diabatic tendencies (will be used in next timestep):
          rdt = 1.0/dt
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qpten(i,j,k)=(pp3d(i,j,k)-qpten(i,j,k))*rdt
            qtten(i,j,k)=(th3d(i,j,k)-qtten(i,j,k))*rdt
            qvten(i,j,k)=(q3d(i,j,k,nqv)-qvten(i,j,k))*rdt
            qcten(i,j,k)=(q3d(i,j,k,nqc)-qcten(i,j,k))*rdt
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
        ENDIF

!-----------------------------------------------------------------
!  Equate the two arrays

#ifdef MPI
      call comm_3u_end(u3d,uw31,uw32,ue31,ue32,   &
                           us31,us32,un31,un32,reqs_u)
      call comm_3v_end(v3d,vw31,vw32,ve31,ve32,   &
                           vs31,vs32,vn31,vn32,reqs_v)
      call comm_3w_end(w3d,ww31,ww32,we31,we32,   &
                           ws31,ws32,wn31,wn32,reqs_w)
#endif

      if(terrain_flag)then
        call bcwsfc(gz,dzdx,dzdy,u3d,v3d,w3d)
        call bc2d(w3d(ib,jb,1))
      endif

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
        do j=0,nj+1
        do i=0,ni+2
          ua(i,j,k)=u3d(i,j,k)
        enddo
        enddo
        do j=0,nj+2
        do i=0,ni+1
          va(i,j,k)=v3d(i,j,k)
        enddo
        enddo
        do j=0,nj+1
        do i=0,ni+1
          wa(i,j,k)=w3d(i,j,k)
        enddo
        enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()

!----------

      if(iturb.eq.1)then
#ifdef MPI
        call comm_3t_end(tke3d,tkw1,tkw2,tke1,tke2,   &
                               tks1,tks2,tkn1,tkn2,reqs_tk)
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=0,nj+1
        do i=0,ni+1
          tkea(i,j,k)=tke3d(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_integ=time_integ+mytime()
      endif

!----------

      if(iptra.eq.1)then
        do n=1,npt
#ifdef MPI
          call comm_3s_end(pt3d(ib,jb,kb,n),                           &
                tw1(1,1,1,n),tw2(1,1,1,n),te1(1,1,1,n),te2(1,1,1,n),   &
                ts1(1,1,1,n),ts2(1,1,1,n),tn1(1,1,1,n),tn2(1,1,1,n),   &
                reqs_t(1,n))
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            pta(i,j,k,n)=pt3d(i,j,k,n)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_integ=time_integ+mytime()
        enddo
      endif

!----------

#ifdef MPI
      call comm_3r_end(th3d,pp3d,rw31,rw32,re31,re32,   &
                                 rs31,rs32,rn31,rn32,reqs_p)
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          ppi(i,j,k)=pp3d(i,j,k)
        enddo
        enddo
        do j=0,nj+1
        do i=0,ni+1
          tha(i,j,k)=th3d(i,j,k)
        enddo
        enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()

!----------

      if(imoist.eq.1)then
#ifdef MPI
        call comm_3q_end(q3d,qw31,qw32,qe31,qe32,   &
                             qs31,qs32,qn31,qn32,reqs_q(1,1))
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k,n)
        do k=1,nk
          do n=1,numq
          do j=0,nj+1
          do i=0,ni+1
            qa(i,j,k,n)=q3d(i,j,k,n)
          enddo
          enddo
          enddo
        enddo
        if(timestats.ge.1) time_integ=time_integ+mytime()
      endif

!----------

#ifdef MPI
        call comm_1s_end(rho,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
        call getcorner(rho,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
        call bcs2(rho)
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=0,nj+1
          k = 1
          do i=0,ni+1
            ! cm1r17, 2nd-order extrapolation:
            rf(i,j,1) = cgs1*rho(i,j,1)+cgs2*rho(i,j,2)+cgs3*rho(i,j,3)
            rr(i,j,1) = 1.0/rho(i,j,1)
          enddo
          do k=2,nk
          do i=0,ni+1
            rf(i,j,k) = (c1(i,j,k)*rho(i,j,k-1)+c2(i,j,k)*rho(i,j,k))
            rr(i,j,k) = 1.0/rho(i,j,k)
          enddo
          enddo
          do i=0,ni+1
            ! cm1r17, 2nd-order extrapolation:
            rf(i,j,nk+1) = cgt1*rho(i,j,nk)+cgt2*rho(i,j,nk-1)+cgt3*rho(i,j,nk-2)
          enddo
        enddo
        if(timestats.ge.1) time_prsrho=time_prsrho+mytime()

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc  Update parcel locations  ccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IF(iprcl.eq.1)THEN
        ! cm1r18:  use velocities averaged over small time steps (for psolver=2,3,6)
        !          (note:  parcel_driver assumes no terrain)
        IF( psolver.eq.2 .or. psolver.eq.3 .or. psolver.eq.6 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
            do j=1,nj
            do i=1,ni+1
              us(i,j,k)=rru(i,j,k)*rr0(1,1,k)
            enddo
            enddo
            IF(axisymm.eq.0)THEN
              ! Cartesian grid:
              do j=1,nj+1
              do i=1,ni
                vs(i,j,k)=rrv(i,j,k)*rr0(1,1,k)
              enddo
              enddo
            ENDIF
            IF(k.gt.1)THEN
              do j=1,nj
              do i=1,ni
                ws(i,j,k)=rrw(i,j,k)*rrf0(1,1,k)
              enddo
              enddo
            ENDIF
          enddo
        ELSE
          ! psolver=1,4,5:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
            do j=1,nj
            do i=1,ni+1
              us(i,j,k)=u3d(i,j,k)
            enddo
            enddo
            IF(axisymm.eq.0)THEN
              ! Cartesian grid:
              do j=1,nj+1
              do i=1,ni
                vs(i,j,k)=v3d(i,j,k)
              enddo
              enddo
            ENDIF
            IF(k.gt.1)THEN
              do j=1,nj
              do i=1,ni
                ws(i,j,k)=w3d(i,j,k)
              enddo
              enddo
            ENDIF
          enddo
        ENDIF
        if(timestats.ge.1) time_parcels=time_parcels+mytime()
        call     parcel_driver(dt,xh,uh,ruh,xf,yh,vh,rvh,yf,zh,mh,rmh,zf,mf,    &
                               znt,rho,us,vs,ws,pdata,packet(1,1),ploc,         &
                               reqs_p,reqs_u,reqs_v,reqs_w,                     &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,                 &
                               nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                               n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,         &
                               uw31,uw32,ue31,ue32,us31,us32,un31,un32,         &
                               vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,         &
                               ww31,ww32,we31,we32,ws31,ws32,wn31,wn32)
        if(timestats.ge.1) time_parcels=time_parcels+mytime()
      ENDIF


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   All done   cccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!!!#ifdef MPI
!!!      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!!!      if(timestats.ge.1) time_mpb=time_mpb+mytime()
!!!#endif

!--------------------------------------------------------------------
!  Calculate surface "swaths."  Move surface (if necessary). 
!--------------------------------------------------------------------

    IF( output_rain.eq.1 )THEN

      if(imove.eq.1.and.imoist.eq.1)then
        weps = 10.0*epsilon
        call movesfc(0.0,dt,weps,uh,vh,rain(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

    ENDIF

!--------------------------------------------------------------------
! Maximum horizontal wind speed at lowest model level: 
! (include domain movement in calculation)

    IF( output_sws.eq.1 )THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = sqrt( (umove+0.5*(ua(i,j,1)+ua(i+1,j,1)))**2    &
                   +(vmove+0.5*(va(i,j,1)+va(i,j+1,1)))**2 ) 
        do n=1,nrain
          sws(i,j,n)=max(sws(i,j,n),tem)
        enddo
      enddo
      enddo
      if(timestats.ge.1) time_swath=time_swath+mytime()

      if(imove.eq.1)then
        weps = 10.0*epsilon
        call movesfc(0.0,dt,weps,uh,vh,sws(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3),  &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

    ENDIF

!--------------------------------------------------------------------
!  Maximum vertical vorticity at lowest model level:

  IF( output_svs.eq.1 )THEN

  IF(axisymm.eq.0)THEN
    IF(.not.terrain_flag)THEN
      ! Cartesian grid, without terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj+1
      do i=1,ni+1
        tem = (va(i,j,1)-va(i-1,j,1))*rdx*uf(i)   &
             -(ua(i,j,1)-ua(i,j-1,1))*rdy*vf(j)
        do n=1,nrain
          svs(i,j,n)=max(svs(i,j,n),tem)
        enddo
      enddo
      enddo
    ELSE
      ! Cartesian grid, with terrain:
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
        ! interior:
        do k=2,2
        r2 = (sigmaf(k)-sigma(k-1))*rds(k)
        r1 = 1.0-r2
        do i=0,ni+2
          dum1(i,j,k) = r1*ua(i,j,k-1)+r2*ua(i,j,k)
          dum2(i,j,k) = r1*va(i,j,k-1)+r2*va(i,j,k)
        enddo
        enddo
      enddo
      k = 1
!$omp parallel do default(shared)  &
!$omp private(i,j,n,r1,tem)
      do j=1,nj+1
      do i=1,ni+1
        r1 = zt/(zt-0.25*((zs(i-1,j-1)+zs(i,j))+(zs(i-1,j)+zs(i,j-1))))
        tem = ( r1*(va(i,j,k)*rgzv(i,j)-va(i-1,j,k)*rgzv(i-1,j))*rdx*uf(i)  &
               +0.5*( (zt-sigmaf(k+1))*(dum2(i-1,j,k+1)+dum2(i,j,k+1))      &
                     -(zt-sigmaf(k  ))*(dum2(i-1,j,k  )+dum2(i,j,k  ))      &
                    )*rdsf(k)*r1*(rgzv(i,j)-rgzv(i-1,j))*rdx*uf(i) )        &
             -( r1*(ua(i,j,k)*rgzu(i,j)-ua(i,j-1,k)*rgzu(i,j-1))*rdy*vf(j)  &
               +0.5*( (zt-sigmaf(k+1))*(dum1(i,j-1,k+1)+dum1(i,j,k+1))      &
                     -(zt-sigmaf(k  ))*(dum1(i,j-1,k  )+dum1(i,j,k  ))      &
                    )*rdsf(k)*r1*(rgzu(i,j)-rgzu(i,j-1))*rdy*vf(j) )
        do n=1,nrain
          svs(i,j,n)=max(svs(i,j,n),tem)
        enddo
      enddo
      enddo
    ENDIF
  ELSE
      ! Axisymmetric grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj+1
      do i=1,ni+1
        tem = (va(i,j,1)*xh(i)-va(i-1,j,1)*xh(i-1))*rdx*uf(i)*rxf(i)
        do n=1,nrain
          svs(i,j,n)=max(svs(i,j,n),tem)
        enddo
      enddo
      enddo
  ENDIF
      if(timestats.ge.1) time_swath=time_swath+mytime()

      if(imove.eq.1)then
        weps = 1.0*epsilon
        call movesfc(-1000.0,dt,weps,uh,vh,svs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

  ENDIF

!--------------------------------------------------------------------
!  Minimum pressure perturbation at lowest model level:

  IF( output_sps.eq.1 )THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = prs(i,j,1)-prs0(i,j,1)
        do n=1,nrain
          sps(i,j,n)=min(sps(i,j,n),tem)
        enddo
      enddo
      enddo
      if(timestats.ge.1) time_swath=time_swath+mytime()

      if(imove.eq.1)then
        weps = 1000.0*epsilon
        call movesfc(-200000.0,dt,weps,uh,vh,sps(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

  ENDIF

!--------------------------------------------------------------------
!  Maximum rainwater mixing ratio (qr) at lowest model level:

  IF( output_srs.eq.1 )THEN

    IF(imoist.eq.1.and.nqr.ne.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = qa(i,j,1,nqr)
        do n=1,nrain
          srs(i,j,n)=max(srs(i,j,n),tem)
        enddo
      enddo
      enddo
      if(timestats.ge.1) time_swath=time_swath+mytime()

      if(imove.eq.1)then
        weps = 0.01*epsilon
        call movesfc(0.0,dt,weps,uh,vh,srs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3),  &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif
    ENDIF

  ENDIF

!--------------------------------------------------------------------
!  Maximum graupel/hail mixing ratio (qg) at lowest model level:

  IF( output_sgs.eq.1 )THEN

    IF(imoist.eq.1.and.nqg.ne.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = qa(i,j,1,nqg)
        do n=1,nrain
          sgs(i,j,n)=max(sgs(i,j,n),tem)
        enddo
      enddo
      enddo
      if(timestats.ge.1) time_swath=time_swath+mytime()

      if(imove.eq.1)then
        weps = 0.01*epsilon
        call movesfc(0.0,dt,weps,uh,vh,sgs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3),  &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif
    ENDIF

  ENDIF

!--------------------------------------------------------------------

  IF( output_sus.eq.1 )THEN

      ! get height AGL:
      if( .not. terrain_flag )then
        ! without terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dum3(i,j,k) = zh(i,j,k)
          wten(i,j,k) = zf(i,j,k)
        enddo
        enddo
        enddo
      else
        ! get height AGL:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dum3(i,j,k) = zh(i,j,k)-zs(i,j)
          wten(i,j,k) = zf(i,j,k)-zs(i,j)
        enddo
        enddo
        enddo
      endif

!--------------------------------------------------------------------
!  Maximum updraft velocity (w) at 5 km AGL:

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,tem)
      do j=1,nj
      do i=1,ni
        k = 2
        ! wten is height AGL:
        do while( wten(i,j,k).lt.5000.0 .and. k.lt.nk )
          k = k + 1
        enddo
        tem = w3d(i,j,k)
        do n=1,nrain
          sus(i,j,n)=max(sus(i,j,n),tem)
        enddo
      enddo
      enddo
      if(timestats.ge.1) time_swath=time_swath+mytime()

      if(imove.eq.1)then
        weps = 10.0*epsilon
        call movesfc(-1000.0,dt,weps,uh,vh,sus(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

    ENDIF

!--------------------------------------------------------------------
!  Maximum integrated updraft helicity:

    IF( output_shs.eq.1 )THEN

      ! dum3 is zh (agl), wten is zf (agl)
      call calcuh(uf,vf,dum3,wten,ua,va,wa,dum1(ib,jb,1),dum2,dum5,dum6, &
                  zs,rgzu,rgzv,rds,sigma,rdsf,sigmaf)
!$omp parallel do default(shared)  &
!$omp private(i,j,n)
      do j=1,nj
      do i=1,ni
        do n=1,nrain
          shs(i,j,n)=max(shs(i,j,n),dum1(i,j,1))
        enddo
      enddo
      enddo
      if(timestats.ge.1) time_swath=time_swath+mytime()

      if(imove.eq.1)then
        weps = 100.0*epsilon
        call movesfc(0.0,dt,weps,uh,vh,shs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3),  &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

    ENDIF

!  Done with "swaths"
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Step time forward, take care of a few odds and ends .... 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      mtime = mtime + dbldt

      if( convinit.eq.1 )then
        if( mtime.gt.convtime ) convinit = 0
      endif

      rtime=sngl(mtime)

      dostat = .false.
      dowrite = .false.
      dorestart = .false.
      doprclout = .false.

      if( mtime.ge.(stattim-0.1*dt) .or. statfrq.lt.0.0 .or. stopit ) dostat = .true.
      if( mtime.ge.(taptim-0.1*dt) .or. stopit ) dowrite = .true.
      if( mtime.ge.(rsttim-0.1*dt) .and. rstfrq.gt.0.0 ) dorestart = .true.
      if(iprcl.eq.1)then
        IF( mtime.ge.(prcltim-0.1*dt) .or. prclfrq.lt.0.0 ) doprclout = .true.
      endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Adaptive timestepping:
!   (assumes cflmax has already been calculated)
!  cm1r18:  get new timestep before calling sfc_and_turb
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IF( adapt_dt.eq.1 )THEN
        dtlast = dt

      IF( myid.eq.0 )THEN
        ! only processor 0 does this:

        cflmax = max(cflmax,1.0e-10)

        IF( cflmax.gt.cfl_limit )THEN
          ! decrease timestep:
          dbldt = dt*cfl_limit/cflmax
        ELSE
          ! increase timestep:
          dbldt = min( 1.0+max_change , cfl_limit/cflmax )*dt
        ENDIF

        ! don't allow dt to exceed twice initial timestep
        dbldt = min( dbldt , dble(2.0*dtl) )

      IF( taptim.gt.0.0 )THEN
        ! ramp-down timestep when approaching output time
        if( dowrite )then
          tout = ( ( taptim - mtime ) + tapfrq )
        else
          tout = ( taptim - mtime )
        endif
        if( tout.gt.(2.0*dbldt) .and. tout.le.(3.0*dbldt)  )then
          dbldt = tout/3.0
        elseif( tout.gt.dbldt .and. tout.le.(2.0*dbldt)  )then
          dbldt = tout/2.0
        elseif( tout.le.dbldt )then
          dbldt = tout
        endif
      ENDIF

      IF( rsttim.gt.0.0 )THEN
        ! ramp-down timestep when approaching restart time
        if( dorestart )then
          tout = ( ( rsttim - mtime ) + rstfrq )
        else
          tout = ( rsttim - mtime )
        endif
        if( tout.gt.(2.0*dbldt) .and. tout.le.(3.0*dbldt)  )then
          dbldt = tout/3.0
        elseif( tout.gt.dbldt .and. tout.le.(2.0*dbldt)  )then
          dbldt = tout/2.0
        elseif( tout.le.dbldt )then
          dbldt = tout
        endif
      ENDIF

      IF( stattim.gt.0.0 )THEN
        ! ramp-down timestep when approaching stat time
        if( dostat )then
          tout = ( ( stattim - mtime ) + statfrq )
        else
          tout = ( stattim - mtime )
        endif
        if( tout.gt.(2.0*dbldt) .and. tout.le.(3.0*dbldt)  )then
          dbldt = tout/3.0
        elseif( tout.gt.dbldt .and. tout.le.(2.0*dbldt)  )then
          dbldt = tout/2.0
        elseif( tout.le.dbldt )then
          dbldt = tout
        endif
      ENDIF

        ! end of processor 0 stuff
      ENDIF

        ! all processors:
#ifdef MPI
        call MPI_BCAST(dbldt ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!!        call MPI_BCAST(dt    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
!!!        call MPI_BCAST(nsound,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
        dt = dbldt
        call dtsmall(dt)
        ndt = ndt + 1
        adt = adt + dt
        acfl = acfl + cflmax
        if(timestats.ge.1) time_misc=time_misc+mytime()
        IF( dt.ne.dtlast )THEN
          IF( (imoist.eq.1).and.(ptype.eq.2) )then
            call consat2(dt)
            if(timestats.ge.1) time_microphy=time_microphy+mytime()
          ENDIF
          IF( (imoist.eq.1).and.(ptype.eq.4) )then
            call lfoice_init(dt)
            if(timestats.ge.1) time_microphy=time_microphy+mytime()
          ENDIF
        ENDIF
      ENDIF

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   Prepare turbulence vars for next time step   cccccccccccccccccc
!cc     (new since cm1r17)                         cccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sfc_and_turb(.true.,nstep,dt,dosfcflx,cloudvar,qbudget,   &
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
                   rtime,ntdiag,tdiag,.true.)
        ! end_sfc_and_turb

!--------------------------------------------------------------------
!  get stress terms for explicit diffusion scheme:

      IF( ((idiff.ge.1).and.(difforder.eq.2)) .or. (dns.eq.1) )THEN
        call diff2def(uh,arh1,arh2,uf,arf1,arf2,vh,vf,mh,c1,c2,mf,ust,znt,u1,v1,s1,  &
                      divx,rho,rr,rf,t11,t12,t13,t22,t23,t33,ua,va,wa,dissten)
      ENDIF

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Get statistics
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      reset = .false.

      if( dostat )then
        IF(axisymm.eq.0)THEN
          ! for Cartesian grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ppten(i,j,k)=rho(i,j,k)
          enddo
          enddo
          enddo
        ELSE
          ! for axisymmetric grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ppten(i,j,k) = rho(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
          enddo
          enddo
          enddo
        ENDIF
        call     statpack(nrec,ndt,dt,dtlast,rtime,adt,acfl,cloudvar,qname,budname,qbudget,asq,bsq, &
                          xh,rxh,uh,ruh,xf,uf,yh,vh,rvh,vf,zh,mh,rmh,mf,    &
                          zs,rgzu,rgzv,rds,sigma,rdsf,sigmaf,               &
                          rstat,pi0,rho0,thv0,th0,qv0,u0,v0,                &
                          dum1,dum2,dum3,dum4,dum5,ppten,prs,               &
                          ua,va,wa,ppi,tha,qa,qten,kmh,kmv,khh,khv,tkea,pta,u10,v10,reset)
        if(statfrq.gt.0.0) stattim=stattim+statfrq
      endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Writeout and stuff
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(myid.eq.0)then
        if(timeformat.eq.1)then
          write(6,110) nstep,rtime,' sec '
        elseif(timeformat.eq.2)then
          write(6,110) nstep,rtime/60.0,' min '
        elseif(timeformat.eq.3)then
          write(6,110) nstep,rtime/3600.0,' hour'
        elseif(timeformat.eq.4)then
          write(6,110) nstep,rtime/86400.0,' day '
        else
          write(6,110) nstep,rtime,' sec'
        endif
110     format(2x,i12,4x,f18.6,a5)
      endif
      if(timestats.ge.1) time_misc=time_misc+mytime()

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if( dowrite )then
        nwrite=nwrite+1
      IF(output_format.eq.1.or.output_format.eq.2)THEN
        nn = 1
        if(terrain_flag .and. output_interp.eq.1) nn = 2
        if(output_format.eq.2) nn = 1
        DO n=1,nn
          if(n.eq.1)then
            fnum = 51
          else
            fnum = 71
          endif
          IF( stopit )THEN
            ! diag code does not account for stopit, so just set arrays to zero
            ! to avoid confusion
            tdiag = 0.0
            qdiag = 0.0
          ENDIF
          call writeout(rtime,dt,fnum,nwrite,qname,xh,xf,uf,yh,yf,vf,xfref,yfref,              &
                        rds,sigma,rdsf,sigmaf,zh,zf,mf,gx,gy,                                  &
                        pi0,prs0,rho0,rr0,rf0,rrf0,th0,qv0,u0,v0,                              &
                        zs,rgzu,rgzv,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,psfc,      &
                        rxh,arh1,arh2,uh,ruh,rxf,arf1,arf2,vh,rvh,mh,rmf,rr,rf,                &
                        gz,rgz,gzu,gzv,gxu,gyv,dzdx,dzdy,c1,c2,                                &
                        cd,ch,cq,tlh,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,                  &
                        t11,t12,t13,t22,t23,t33,rho,prs,sten,                                  &
                        rru,ua,u3d,uten,rrv,va,v3d,vten,rrw,wa,w3d,wten,ppi,tha,               &
                        us,vs,ws,thadv,thten,nm,defv,defh,dissten,                             &
                        thpten,qvpten,qcpten,qipten,upten,vpten,                               &
                        lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,  &
                        qa,kmh,kmv,khh,khv,tkea,swten,lwten,                                   &
                        radsw,rnflx,radswnet,radlwin,dsr,olr,pta,                              &
                        num_soil_layers,u10,v10,t2,q2,znt,ust,u1,v1,s1,                        &
                        hpbl,zol,mol,br,psim,psih,qsfc,                                        &
                        dat1,dat2,dat3,reqt,ntdiag,nqdiag,tdiag,qdiag,                         &
                        nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
        ENDDO
#ifdef HDFOUT
      ELSEIF(output_format.ge.3)THEN
        call writeout_mult_hdf5(rtime,qname,rds,sigma,rdsf,sigmaf,xh,xf,uf,yh,yf,vf,mh,zh,mf,zf, &
                      pi0,prs0,rho0,th0,qv0,u0,v0,                     &
                      zs,gz,rgz,rgzu,rgzv,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,cd,ch,cq,dum1,dum2,dum3,dum4,  &
                      dum5,dum6,rho,prs,sten,ua,uten,va,vten,wa,wten,ppi,tha,    &
                      qa,kmh,kmv,khh,khv,tkea,pta,num_soil_layers,   &
                      lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,   &
                      radsw,rnflx,radswnet,radlwin,u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br,   &
                      dissten,thpten,qvpten,qcpten,qipten,upten,vpten,swten,lwten)
! See above note about sten being used for dbz calculation
#endif
      ENDIF
        if(tapfrq.gt.0.0) taptim=taptim+tapfrq
        if(timestats.ge.1) time_write=time_write+mytime()
      endif

!--------------------------------------------------------------------
!  Write parcel data:
!--------------------------------------------------------------------

      if(iprcl.eq.1)then
      IF( doprclout )THEN
        call     parcel_interp(dt,xh,uh,ruh,xf,uf,yh,vh,rvh,yf,vf,     &
                               zh,mh,rmh,zf,mf,znt,ust,c1,c2,          &
                               pi0,th0,thv0,qv0,qc0,qi0,rth0,          &
                               dum1,dum2,dum3,dum4,dum5,dum6,prs,rho,  &
                               sten,dum7,dum8,wten,wten1,              &
                               u3d,v3d,w3d,pp3d,thten,thten1,th3d,q3d, &
                               kmh,kmv,khh,khv,tke3d,pt3d,pdata,       &
                               packet,reqs_p,                          &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,        &
                               nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
        call parcel_write(prec,rtime,qname,pdata,ploc)
        if(prclfrq.gt.0.0) prcltim = prcltim + prclfrq
        if(timestats.ge.1) time_parcels=time_parcels+mytime()
      ENDIF
      endif

!-------------------------------------------------------------------

      if( dorestart )then
        nrst=nrst+1
        if(rstfrq.gt.0.0) rsttim=rsttim+rstfrq
        call     write_restart(nstep,nrec,prec,nwrite,nrst,num_soil_layers,         &
                               dt,dtlast,mtime,ndt,adt,acfl,dbldt,mass1,            &
                               stattim,taptim,rsttim,radtim,prcltim,                &
                               qbudget,asq,bsq,qname,                               &
                               xfref,yfref,zh,zf,sigma,sigmaf,zs,                   &
                               th0,prs0,pi0,rho0,qv0,u0,v0,                         &
                               rain,sws,svs,sps,srs,sgs,sus,shs,                    &
                               tsk,znt,ust,cd,ch,cq,u1,v1,s1,thflux,qvflux,         &
                               radbcw,radbce,radbcs,radbcn,                         &
                               rho,prs,ua,va,wa,ppi,tha,qa,tkea,                    &
                               swten,lwten,radsw,rnflx,radswnet,radlwin,rad2d,      &
                               effc,effi,effs,effr,effg,effis,                      &
                               lu_index,kpbl2d,psfc,u10,v10,s10,hfx,qfx,xland,      &
                               hpbl,wspd,psim,psih,gz1oz0,br,                       &
                               CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                       &
                               MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                    &
                               CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,       &
                               f2d,gsw,glw,chklowq,capg,snowc,fm,fh,tslb,           &
                               tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,               &
                               qpten,qtten,qvten,qcten,pta,pdata,ploc,ppx,sten,     &
                               ntdiag,nqdiag,tdiag,qdiag,                           &
                               dum1,dat1,dat2,dat3,reqt)
        if(timestats.ge.1) time_restart=time_restart+mytime()
      endif

!-------------------------------------------------------------------

      IF( adapt_dt.eq.1 .and. reset )THEN
        ndt  = 0
        adt  = 0.0
        acfl = 0.0
      ENDIF

      if(stopit)then
        if(myid.eq.0)then
          print *
          print *,' Courant number has exceeded 1.5 '
          print *
          print *,' Stopping model .... '
          print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif

      return
      end subroutine solve

