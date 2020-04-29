
      subroutine param(dt,dtlast,stattim,taptim,rsttim,radtim,prcltim,  &
                       cloudvar,rhovar,qname,budname,                   &
                       xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf, &
                       yh,vh,rvh,yf,vf,rvf,xfref,yfref,                 &
                       rds,sigma,rdsf,sigmaf,tauh,taus,zh,mh,rmh,cc1,cc2,tauf,zf,mf,rmf, &
                       zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,  &
                       ntdiag,nqdiag,                                   &
                       reqs_u,reqs_v,reqs_s,reqs_p,                     &
                       nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                       n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,         &
                       sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,         &
                       uw31,uw32,ue31,ue32,us31,us32,un31,un32,         &
                       vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,         &
                       ww31,ww32,we31,we32,ws31,ws32,wn31,wn32)



      use module_mp_thompson
      use module_mp_graupel
      use module_mp_nssl_2mom, only:    &
                        nssl_2mom_init, &
                        rho_qr,         &
                        cnor,           &
                        rho_qs,         &
                        cnos,           &
                        rho_qh,         &
                        cnoh,           &
                        ccn,            &
                        infall,         &
                        alphah,         &
                        alphahl,        &
                        imurain,        &
                        icdx,           &
                        icdxhl,         &
                        dfrz,           &
                        hldnmn,         &
                        iferwisventr,   &
                        iehw,iehlw,     &
                        ehw0,ehlw0,     &
                        dmrauto,        &
                        ioldlimiter
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real :: dt,dtlast
      double precision :: stattim,taptim,rsttim,radtim,prcltim
      logical, dimension(maxq) :: cloudvar,rhovar
      character*3, dimension(maxq) :: qname
      character*6, dimension(maxq) :: budname
      real, dimension(ib:ie) :: xh,rxh,arh1,arh2,uh,ruh
      real, dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf,ruf
      real, dimension(jb:je) :: yh,vh,rvh
      real, dimension(jb:je+1) :: yf,vf,rvf
      real, dimension(-2:nx+4) :: xfref
      real, dimension(-2:ny+4) :: yfref
      real, dimension(kb:ke) :: rds,sigma
      real, dimension(kb:ke+1) :: rdsf,sigmaf
      real, dimension(ib:ie,jb:je,kb:ke) :: tauh,taus,zh,mh,rmh,cc1,cc2
      real, dimension(ib:ie,jb:je,kb:ke+1) :: tauf,zf,mf,rmf
      real, dimension(ib:ie,jb:je) :: zs
      real, dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      integer, intent(inout) ::ntdiag,nqdiag
      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_s,reqs_p
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(cmp,cmp,kmt+1) :: n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, intent(inout), dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, intent(inout), dimension(imp,cmp,kmp)   :: vs31,vs32,vn31,vn32
      real, intent(inout), dimension(cmp,jmp,kmp-1) :: ww31,ww32,we31,we32
      real, intent(inout), dimension(imp,cmp,kmp-1) :: ws31,ws32,wn31,wn32

!-----------------------------------------------------------------------

      integer i,j,k,n,m,nn,kst,ni1,ni2,ni3,nj1,nj2,nj3,nk1,nk2,nk3
      integer ival,jval
      integer iterrain
      integer :: inum
      real :: var
      real :: zfw1d(kb:ke+1), zfs1d(kb:ke+1)
      double precision :: gzc(kb:ke+1), gze(kb:ke+1)
      integer :: nbndlyr    = 0
      real    :: rtop       = 1.0     ! upper level stretch factor
      real    :: ztopstr    = 100000. ! height to start upper level stretching
      real    :: dzmaxtop   = 700.    ! max upper level dz




      real c1,c2,nominal_dx,nominal_dy,nominal_dz,z1,z2,z3,mult
      real x1,x2,y1,y2
      logical :: doit
      double precision, dimension(:), allocatable :: xfdp,yfdp

      integer, parameter :: bigm = 4   ! highest deriv
      integer, parameter :: bign = 3   ! number of grid points minus 1
      double precision :: x0,b1,b2,b3
      double precision, dimension(0:bign) :: alpha
      ! delta(n,m,nu):
      double precision, dimension(0:bign,-1:bigm,0:bign) :: delta

      namelist /param1/ dx,dy,dz,dtl,timax,run_time,                    &
          tapfrq,rstfrq,statfrq,prclfrq
      namelist /param2/                                                 &
          adapt_dt,irst,rstnum,iconly,                                  &
          hadvordrs,vadvordrs,hadvordrv,vadvordrv,pdscheme,apmasscon,   &
          advwenos,advwenov,idiff,mdiff,difforder,imoist,iturb,         &
          tconfig,bcturbs,dns,                                          &
          irdamp,hrdamp,psolver,nsound,ptype,ihail,iautoc,              &
          icor,pertcor,eqtset,idiss,efall,rterm,                        &
          wbc,ebc,sbc,nbc,bbc,tbc,irbc,roflux,isnd,iwnd,itern,iinit,    &
          irandp,ibalance,iorigin,axisymm,imove,iptra,npt,pdtra,        &
          iprcl,nparcels
      namelist /param3/ kdiff2,kdiff6,fcor,kdiv,alph,rdalpha,zd,xhd,    &
                        umove,vmove,v_t,l_h,lhref1,lhref2,l_inf,ndcnst
      namelist /param4/ stretch_x,dx_inner,dx_outer,nos_x_len,tot_x_len
      namelist /param5/ stretch_y,dy_inner,dy_outer,nos_y_len,tot_y_len
      namelist /param6/ stretch_z,ztop,str_bot,str_top,dz_bot,dz_top
      namelist /param7/ bc_temp,ptc_top,ptc_bot,viscosity,pr_num
      namelist /param8/ var1,var2,var3,var4,var5,var6,var7,var8,var9,var10
      namelist /param9/                                                       &
              output_path,output_basename,output_format,output_filetype,      &
              output_interp,                                                  &
              output_rain,output_sws,output_svs,output_sps,output_srs,        &
              output_sgs,output_sus,output_shs,output_coldpool,               &
              output_sfcflx,output_sfcparams,output_sfcdiags,                 &
              output_psfc,output_zs,output_zh,output_basestate,               &
              output_th,output_thpert,output_prs,output_prspert,              &
              output_pi,output_pipert,output_rho,output_rhopert,output_tke,   &
              output_km,output_kh,                                            &
              output_qv,output_qvpert,output_q,output_dbz,output_buoyancy,    &
              output_u,output_upert,output_uinterp,                           &
              output_v,output_vpert,output_vinterp,output_w,output_winterp,   &
              output_vort,output_pv,output_uh,output_pblten,                  &
              output_dissten,output_dissheat,output_mptend,output_fallvel,    &
              output_nm,output_def,output_turbten,output_impdiften,           &
              output_radten,                                                  &
              restart_format,restart_filetype,                                &
              restart_file_theta,restart_file_dbz,restart_file_th0,           &
              restart_file_prs0,restart_file_pi0,restart_file_rho0,           &
              restart_file_qv0,restart_file_u0,restart_file_v0,               &
              restart_file_zs,restart_file_zh,restart_file_zf,                &
              restart_file_diags,restart_use_theta,restart_reset_frqtim
      namelist /param10/                                                      &
              stat_w,stat_u,stat_v,stat_rmw,stat_pipert,stat_prspert,         &
              stat_thpert,stat_q,                                             &
              stat_tke,stat_km,stat_kh,stat_div,stat_rh,stat_rhi,stat_the,    &
              stat_cloud,stat_sfcprs,stat_wsp,stat_cfl,stat_vort,             &
              stat_tmass,stat_tmois,stat_qmass,stat_tenerg,stat_mo,stat_tmf,  &
              stat_pcn,stat_qsrc
      namelist /param11/                                                      &
              radopt,dtrad,ctrlat,ctrlon,year,month,day,hour,minute,second
      namelist /param12/                                                      &
              isfcflx,sfcmodel,oceanmodel,ipbl,initsfc,                       &
              tsk0,tmn0,xland0,lu0,season,cecd,pertflx,cnstce,cnstcd,         &
              isftcflx,iz0tlnd,oml_hml0,oml_gamma
      namelist /param13/                                                      &
              prcl_th,prcl_t,prcl_prs,prcl_ptra,prcl_q,prcl_nc,               &
              prcl_km,prcl_kh,prcl_tke,prcl_dbz,prcl_b,prcl_vpg,prcl_vort,    &
              prcl_rho,prcl_qsat,prcl_sfc

      NAMELIST /nssl2mom_params/            &
                        rho_qr,         &
                        cnor,           &
                        rho_qs,         &
                        cnos,           &
                        rho_qh,         &
                        cnoh,           &
                        ccn,            &
                        infall,         &
                        alphah,         &
                        alphahl,        &
                        imurain,        &
                        icdx,           &
                        icdxhl,         &
                        dfrz,           &
                        hldnmn,         &
                        iferwisventr,   &
                        iehw,iehlw,     &
                        ehw0,ehlw0,     &
                        dmrauto,        &
                        ioldlimiter

!--------------------------------------------------------------










      if(dowr) write(outfile,*) 'Inside PARAM'


!--------------------------------------------------------------

      if(nodex.ne.1 .or. nodey.ne.1)then
        print *
        print *,'  For non-MPI runs, nodex and nodey must be = 1 !'
        print *
        call stopcm1
      endif


!--------------------------------------------------------------





      open(unit=20,file='namelist.input',form='formatted',status='old',    &
           access='sequential',err=8000)
      read(20,nml=param1,err=8001,end=8001)
      read(20,nml=param2,err=8002,end=8002)
      read(20,nml=param3,err=8003,end=8003)
      read(20,nml=param11,err=8011,end=8011)
      read(20,nml=param12,err=8012,end=8012)
      read(20,nml=param4,err=8004,end=8004)
      read(20,nml=param5,err=8005,end=8005)
      read(20,nml=param6,err=8006,end=8006)
      read(20,nml=param7,err=8007,end=8007)
      read(20,nml=param8,err=8008,end=8008)
      read(20,nml=param9,err=8009,end=8009)
      read(20,nml=param10,err=8010,end=8010)
      if( iprcl.eq.1 )then
        read(20,nml=param13,err=8013,end=8013)
      endif
      IF ( ptype .ge. 26 ) THEN
         read(20,nml=nssl2mom_params,err=8051,end=8051)
      ENDIF
      close(unit=20)



!-----------------------------------------------------------------------
!  Some "dummy" checks:

      eqtset = max( eqtset , 1 )
      eqtset = min( eqtset , 2 )
      if(imoist.ne.1) efall=0
      if(imove.eq.0) umove=0.0
      if(imove.eq.0) vmove=0.0
      irst = max( irst , 0 )
      irst = min( irst , 1 )

      IF( nk.lt.5 )THEN



        if(myid.eq.0)then
        print *,'  nk = ',nk
        print *,'  nk must be >= 5 '
        endif
        call stopcm1
      ENDIF

!------------------------------------------------
!  begin non-fatal checks:

      IF( hadvordrs.lt.5 .or. hadvordrs.gt.6 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  hadvordrs = ',hadvordrs
        print *
        print *,'  This value is invalid ... setting to 5 '
        print *
        print *,'  -------------------------------- '
        endif
        hadvordrs = 5
      ENDIF
      IF( hadvordrv.lt.5 .or. hadvordrv.gt.6 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  hadvordrv = ',hadvordrv
        print *
        print *,'  This value is invalid ... setting to 5 '
        print *
        print *,'  -------------------------------- '
        endif
        hadvordrv = 5
      ENDIF
      IF( vadvordrs.lt.5 .or. vadvordrs.gt.6 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  vadvordrs = ',vadvordrs
        print *
        print *,'  This value is invalid ... setting to 5 '
        print *
        print *,'  -------------------------------- '
        endif
        vadvordrs = 5
      ENDIF
      IF( vadvordrv.lt.5 .or. vadvordrv.gt.6 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  vadvordrv = ',vadvordrv
        print *
        print *,'  This value is invalid ... setting to 5 '
        print *
        print *,'  -------------------------------- '
        endif
        vadvordrv = 5
      ENDIF
      IF( sfcmodel.ge.1 )THEN
        IF( oceanmodel.le.0 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  sfcmodel = ',sfcmodel
        print *
        print *,'  but oceanmodel = ',oceanmodel
        print *
        print *,'  setting oceanmodel to 1  (just in case its needed) '
        print *
        print *,'  -------------------------------- '
        endif
        oceanmodel = 1
        ENDIF
      ENDIF
      IF( (sfcmodel.ge.1) .or. (oceanmodel.eq.2) .or. (ipbl.ge.1) )then
        IF( bbc.ne.3 .and. dns.eq.0 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  sfcmodel   = ',sfcmodel
        print *,'  oceanmodel = ',oceanmodel
        print *,'  ipbl       = ',ipbl
        print *
        print *,'  at least one of these options requires bbc = 3 '
        print *,'  ... so, setting bbc to 3 '
        print *
        print *,'  -------------------------------- '
        endif
        bbc = 3
        ENDIF
      ENDIF
      IF( irst.eq.1 .and. iinit.ne.0 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  irst       = ',irst
        print *,'  iinit      = ',iinit
        print *
        print *,'  This is a restart. '
        print *,'  so, setting iinit to 0 '
        print *
        print *,'  -------------------------------- '
        endif
        iinit = 0
      ENDIF
      IF( irst.eq.1 .and. irandp.ne.0 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  irst       = ',irst
        print *,'  irandp     = ',irandp
        print *
        print *,'  This is a restart. '
        print *,'  so, setting irandp to 0 '
        print *
        print *,'  -------------------------------- '
        endif
        irandp = 0
      ENDIF
      IF( (sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4) .and. isfcflx.eq.0 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  sfcmodel   = ',sfcmodel
        print *,'  isfcflx    = ',isfcflx
        print *
        print *,'  sfcmodel=2,3,4 requires isfcflx=1 '
        print *,'  so, setting isfcflx to 1 '
        print *
        print *,'  -------------------------------- '
        endif
        isfcflx = 1
      ENDIF
      IF( (psolver.eq.4.or.psolver.eq.5.or.psolver.eq.6) .and. eqtset.eq.2 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  psolver    = ',psolver
        print *,'  eqtset     = ',eqtset
        print *
        print *,'  psolver=4,5,6 requires eqtset=1 '
        print *,'  ... setting eqtset to 1 ... '
        print *
        print *,'  -------------------------------- '
        endif
        eqtset = 1
      ENDIF
      IF( (psolver.eq.4.or.psolver.eq.5.or.psolver.eq.6) .and. apmasscon.eq.1 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  psolver    = ',psolver
        print *,'  apmasscon  = ',apmasscon
        print *
        print *,'  psolver=4,5,6 requires apmasscon=0 '
        print *,'  ... setting apmasscon to 0 ... '
        print *
        print *,'  -------------------------------- '
        endif
        apmasscon = 0
      ENDIF
      IF( restart_use_theta )THEN
      IF( .not. restart_file_theta )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  restart_use_theta = ',restart_use_theta
        print *
        print *,'  ... setting restart_file_theta to true ... '
        print *
        print *,'  -------------------------------- '
        endif
        restart_file_theta = .true.
      ENDIF
      ENDIF

      IF( restart_filetype.ge.3 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  Single processor run: '
        print *
        print *,'  restart_filetype >= 3 not available '
        print *
        print *,'  ... setting restart_filetype to 2 ... '
        print *
        print *,'  -------------------------------- '
        endif
        restart_filetype = 2
      ENDIF

      IF( restart_format.eq.1 .and. restart_filetype.eq.1 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  binary-format restart file: '
        print *
        print *,'  restart_filetype = 1 not available '
        print *
        print *,'  ... setting restart_filetype to 2 ... '
        print *
        print *,'  -------------------------------- '
        endif
        restart_filetype = 2
      ENDIF

      IF( restart_format.eq.2 .and. restart_filetype.ge.3 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  netcdf-format restart file: '
        print *
        print *,'  restart_filetype >= 3 not available '
        print *
        print *,'  ... setting restart_filetype to 2 ... '
        print *
        print *,'  -------------------------------- '
        endif
        restart_filetype = 2
      ENDIF

      IF( (iturb.eq.3) .and. (tconfig.ne.2) )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  iturb=3  requires  tconfig=2 '
        print *
        print *,'  ... setting tconfig to 2 ... '
        print *
        print *,'  -------------------------------- '
        endif
        tconfig = 2
      ENDIF
      IF( ipbl.ge.1 .and. iturb.eq.3 .and. abs(l_inf).gt.1.0e-6 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  ipbl >= 1  requires  l_inf = 0 '
        print *
        print *,'  ... setting l_inf to 0 ... '
        print *
        print *,'  -------------------------------- '
        endif
        l_inf = 0.0
      ENDIF

!  end non-fatal checks:
!------------------------------------------------
!  begin fatal checks  (ie, model stops)

      IF( psolver.lt.1 .or. psolver.gt.6 )THEN
        if(myid.eq.0)then
        print *
        print *,'  psolver  = ',psolver
        print *
        print *,'  invalid value for psolver '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( psolver.eq.1 .and. adapt_dt.eq.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  psolver  = ',psolver
        print *,'  adapt_dt = ',adapt_dt
        print *
        print *,'  Cannot use adapt_dt with psolver=1 '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( psolver.eq.2 .or. psolver.eq.3 )THEN
        IF( nsound.lt.4 )THEN
          if(myid.eq.0)then
          print *
          print *,'  nsound = ',nsound
          print *
          print *,'  nsound must be >= 4 '
          print *
          print *,'   stopping model .... '
          print *
          endif



          call stopcm1
        ENDIF
        IF( mod(nsound,2).ne.0 )THEN
          if(myid.eq.0)then
          print *
          print *,'  nsound = ',nsound
          print *
          print *,'  nsound must be an even integer '
          print *
          print *,'   stopping model .... '
          print *
          endif



          call stopcm1
        ENDIF
      ENDIF
      IF(dns.gt.1.or.dns.lt.0)THEN
        if(myid.eq.0)then
        print *
        print *,'  dns   = ',dns
        print *
        print *,'  dns must be either 0 or 1'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(sfcmodel.eq.1)THEN
      IF(cecd.lt.1.or.cecd.gt.3)THEN
        if(myid.eq.0)then
        print *
        print *,'  cecd  = ',cecd
        print *
        print *,'  cecd must be 1,2,3'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      ENDIF
      IF(iturb.gt.3.or.iturb.lt.0)THEN
        if(myid.eq.0)then
        print *
        print *,'  iturb   = ',iturb
        print *
        print *,'  iturb must be either 0, 1, 2, or 3'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(iturb.ge.1 .and. dns.ge.1)THEN
        if(myid.eq.0)then
        print *
        print *,'  iturb = ',iturb
        print *,'  dns   = ',dns
        print *
        print *,'  For dns = 1, iturb must be 0'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(iturb.ge.1 .and. (idiff.eq.1.and.difforder.eq.2))THEN
        if(myid.eq.0)then
        print *
        print *,'  iturb     = ',iturb
        print *,'  idiff     = ',idiff
        print *,'  difforder = ',difforder
        print *
        print *,'  For idiff=1 with difforder=2, iturb > 0 cannot be used '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(bcturbs.lt.1.or.bcturbs.gt.2)THEN
        if(myid.eq.0)then
        print *
        print *,'  bcturbs = ',bcturbs
        print *
        print *,'  bcturbs must be 1 or 2'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(dns.ge.1 .and. imoist.ge.1)THEN
        if(myid.eq.0)then
        print *
        print *,'  imoist = ',imoist
        print *,'  dns    = ',dns
        print *
        print *,'  For dns = 1, imoist must be 0'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(imoist.eq.1 .and. (             isnd.eq.2   &
                        .or.isnd.eq.3.or.isnd.eq.8.) )THEN
        if(myid.eq.0)then
        print *
        print *,'  imoist = ',imoist
        print *,'  isnd   = ',isnd
        print *
        print *,'  For this value of isnd, imoist must be 0'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( bbc.lt.1 .or. bbc.gt.3 )THEN
        if(myid.eq.0)then
        print *
        print *,'  bbc = ',bbc
        print *
        print *,'  bbc must be 1, 2, or 3'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( tbc.lt.1 .or. tbc.gt.2 )THEN
        if(myid.eq.0)then
        print *
        print *,'  tbc = ',tbc
        print *
        print *,'  tbc must be 1 or 2'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(dns.eq.1 .and. (bc_temp.le.0 .or. bc_temp.ge.3))THEN
        if(myid.eq.0)then
        print *
        print *,'  dns     = ',dns
        print *,'  bc_temp = ',bc_temp
        print *
        print *,'  for dns = 1, bc_temp must be either 1 or 2'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(ihail.lt.0.or.ihail.gt.1)THEN
        if(myid.eq.0)then
        print *
        print *,'  ihail   = ',ihail
        print *
        print *,'  ihail must be 0 or 1'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(imoist.eq.1.and.output_dbz.eq.1.and.ptype.ne.2.and.ptype.ne.3.and.ptype.ne.5  &
          .and. (.not. ptype.ge.26))then
        if(myid.eq.0)then
        print *
        print *,'  ptype      = ',ptype
        print *,'  output_dbz = ',output_dbz
        print *
        print *,'  output_dbz is only available for ptype=2,3,5,26,27,28'
        print *
        endif
        IF(ptype.eq.4)THEN
          print *,'   stopping model .... '
          print *



          call stopcm1
        ELSE
          output_dbz = 0
        ENDIF
      ENDIF
      IF(imoist.eq.1 .and. eqtset.ge.2 .and. ptype.eq.4)THEN
        if(myid.eq.0)then
        print *
        print *,'  eqtset  = ',eqtset
        print *,'  ptype   = ',ptype
        print *
        print *,'  eqtset = 2 is not available for ptype = 4'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(imoist.eq.1 .and. efall.eq.1)THEN
      IF(ptype.ne.1.and.ptype.ne.2.and.ptype.ne.5.and.ptype.ne.6)THEN
        if(myid.eq.0)then
        print *
        print *,'  efall   = ',efall
        print *,'  ptype   = ',ptype
        print *
        print *,'  efall = 1 is only supported with ptype = 1,2,5,6'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      ENDIF
      IF((imoist.eq.1).and.(ptype.eq.4).and.terrain_flag)THEN
        if(myid.eq.0)then
        print *
        print *,'  ptype   = ',ptype
        print *,'  terrain_flag = ',terrain_flag
        print *
        print *,'  ptype = 4 does not work with terrain '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(imoist.eq.1 .and. ndcnst.le.1.0e-6)THEN
        if(myid.eq.0)then
        print *
        print *,'  imoist  = ',imoist
        print *,'  ndcnst  = ',ndcnst
        print *
        print *,'  ndcnst is too small.  Please enter a larger value '
        print *,'  in the param3 section of namelist.input '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(terrain_flag .and. (iprcl.ne.0) )THEN
        if(myid.eq.0)then
        print *
        print *,'  terrain_flag = ',terrain_flag
        print *,'  iprcl        = ',iprcl
        print *
        print *,'  cannot use parcels with terrain (for now) '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(terrain_flag .and. (psolver.eq.4.or.psolver.eq.5) )THEN
        if(myid.eq.0)then
        print *
        print *,'  terrain_flag = ',terrain_flag
        print *,'  psolver      = ',psolver
        print *
        print *,'  for psolver = 4,5 terrain_flag must be .false.'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (psolver.eq.4.or.psolver.eq.5) .and.    &
          (wbc.eq.2.or.ebc.eq.2.or.sbc.eq.2.or.nbc.eq.2) )THEN
        if(myid.eq.0)then
        print *
        print *,'  psolver = ',psolver
        print *
        print *,'  cannot use open boundary conditions for psolver = 4 and 5 (at the moment)'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(terrain_flag .and. ibalance.eq.2)THEN
        if(myid.eq.0)then
        print *
        print *,'  terrain_flag = ',terrain_flag
        print *,'  ibalance     = ',ibalance
        print *
        print *,'  for ibalance.eq.2, terrain_flag must be .false.'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF(iinit.eq.6)THEN
        if(myid.eq.0)then
        print *
        print *,'  iinit        = ',iinit
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (output_format.le.0) .or. (output_format.ge.6) )THEN
        if(myid.eq.0)then
        print *
        print *,'  output_format = ',output_format
        print *
        print *,'  only output_format = 1,2,3,4,5 are currently supported'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (axisymm.eq.1) .and. (iorigin.ne.1) )THEN
        if(myid.eq.0)then
        print *
        print *,'  iorigin = ',iorigin
        print *
        print *,'  axisymm=1 requires iorigin=1'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (axisymm.eq.1) .and. (imove.ne.0) )THEN
        if(myid.eq.0)then
        print *
        print *,'  imove = ',imove
        print *
        print *,'  axisymm=1 requires imove=0'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( terrain_flag .and. (imove.ne.0) )THEN
        if(myid.eq.0)then
        print *
        print *,'  imove = ',imove
        print *
        print *,'  imove must be 0 when using terrain '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (axisymm.eq.1) .and. terrain_flag )THEN
        if(myid.eq.0)then
        print *
        print *,'  terrain_flag = ',terrain_flag
        print *
        print *,'  axisymm=1 cannot be used with terrain '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( icor.eq.0 ) fcor = 0.0
      IF( (axisymm.eq.1) .and. (wbc.ne.3) )THEN
        if(myid.eq.0)then
        print *
        print *,'  wbc = ',wbc
        print *
        print *,'  axisymm=1 requires wbc=3 '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (axisymm.eq.1) .and. ( (sbc.ne.1).or.(nbc.ne.1) ) )THEN
        if(myid.eq.0)then
        print *
        print *,'  sbc = ',sbc
        print *,'  nbc = ',nbc
        print *
        print *,'  axisymm=1 requires sbc=nbc=1 '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (axisymm.eq.1).and.(ny.gt.1) )THEN
        if(myid.eq.0)then
        print *
        print *,'  ny = ',ny
        print *
        print *,'  axisymm=1 requires ny=1'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (axisymm.eq.1.and.iturb.ge.1).and.iturb.ne.3 )THEN
        if(myid.eq.0)then
        print *
        print *,'  iturb    = ',iturb
        print *
        print *,'  axisymm=1 is only available with iturb=3'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( axisymm.eq.1.and.(psolver.eq.1.or.psolver.eq.4.or.psolver.eq.5) )THEN
        if(myid.eq.0)then
        print *
        print *,'  psolver    = ',psolver
        print *
        print *,'  axisymm=1 is only available with psolver=2,3,6'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( axisymm.eq.1 .and. idiff.eq.1 .and. difforder.eq.6 )THEN
        idiff = 0
        difforder = 0
      ENDIF
      IF( (bbc.eq.3) .and. (sfcmodel.le.0) )THEN
        if(myid.eq.0)then
        print *
        print *,'  bbc      = ',bbc
        print *,'  sfcmodel = ',sfcmodel
        print *
        print *,'  bbc=3 requires a setting for sfcmodel '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (bbc.eq.3).or.(isfcflx.eq.1) )THEN
      IF( iturb.eq.0 .and. ipbl.eq.0 .and. dns.ne.1 .and. idiff.ne.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  bbc      = ',bbc
        print *,'  isfcflx  = ',isfcflx
        print *
        print *,'  these options require the use of a vertical diffusion/turbulence scheme'
        print *
        print *,'  iturb    = ',iturb
        print *,'  dns      = ',dns
        print *,'  ipbl     = ',ipbl
        print *,'  idiff    = ',idiff
        print *
        print *,'  Use iturb = 1,2,3 or ipbl = 1,2 or idiff = 1'
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      ENDIF
      IF( ipbl.ge.1 .and. (iturb.eq.1.or.iturb.eq.2) )THEN
        if(myid.eq.0)then
        print *
        print *,'  ipbl  = ',ipbl
        print *,'  iturb = ',iturb
        print *
        print *,'  cannot use PBL scheme and LES subgrid turbulence scheme at same time '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( isfcflx.ne.0 )THEN
      IF( sfcmodel.lt.1 .or. sfcmodel.gt.4 )THEN
        if(myid.eq.0)then
        print *
        print *,'  sfcmodel   = ',sfcmodel
        print *
        print *,'  sfcmodel must be 1,2,3,4 '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      ENDIF
      IF( (sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4).and.imove.ne.0 )THEN
        if(myid.eq.0)then
        print *
        print *,'  sfcmodel  = ',sfcmodel
        print *,'  imove     = ',imove
        print *
        print *,'  domain translation is now allowed with sfcmodel = 2,3,4 '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4).and.(season.le.0.or.season.ge.3) )THEN
        if(myid.eq.0)then
        print *
        print *,'  sfcmodel = ',sfcmodel
        print *,'  season   = ',season
        print *
        print *,'  season must have a value of 1 or 2 '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( pertflx.eq.1 .and. sfcmodel.ge.2 )THEN
        if(myid.eq.0)then
        print *
        print *,'  pertflx  = ',pertflx
        print *,'  sfcmodel = ',sfcmodel
        print *
        print *,'  pertflx can only be used with sfcmodel = 1  '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( sfcmodel.eq.1 .and. oceanmodel.ne.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  sfcmodel   = ',sfcmodel
        print *,'  oceanmodel = ',oceanmodel
        print *
        print *,'  sfcmodel = 1 requires oceanmodel = 1 '
        print *,'  (oceanmodel = 2 requires sfcmodel = 2 ) '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( radopt.lt.0 .or. radopt.gt.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  radopt   = ',radopt
        print *
        print *,'  radopt must be 0,1 '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (radopt.eq.1.or.radopt.eq.2) .and. imoist.eq.0 )THEN
        if(myid.eq.0)then
        print *
        print *,'  radopt   = ',radopt
        print *,'  imoist   = ',imoist
        print *
        print *,'  radopt=1 requires imoist=1 (for now) '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( ipbl.ge.1 .and. imoist.eq.0 )THEN
        if(myid.eq.0)then
        print *
        print *,'  ipbl     = ',ipbl
        print *,'  imoist   = ',imoist
        print *
        print *,'  ipbl=1,2 requires imoist=1 (for now) '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (radopt.eq.1.or.radopt.eq.2) .and. rterm.eq.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  radopt   = ',radopt
        print *,'  rterm    = ',rterm
        print *
        print *,'  cannot use radopt and rterm at the same time '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (radopt.eq.1.or.radopt.eq.2) .and. (ptype.eq.1.or.ptype.eq.6)  )THEN
        if(myid.eq.0)then
        print *
        print *,'  radopt   = ',radopt
        print *,'  ptype    = ',ptype
        print *
        print *,'  radopt=1 requires an ice microphysics scheme (for now) '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (radopt.eq.1.or.radopt.eq.2) .and. sfcmodel.eq.0 )THEN
        if(myid.eq.0)then
        print *
        print *,'  radopt   = ',radopt
        print *,'  sfcmodel = ',sfcmodel
        print *
        print *,'  radopt=1 requires a surface model '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF
      IF( (sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4) .and. imoist.eq.0 )THEN
        if(myid.eq.0)then
        print *
        print *,'  sfcmodel = ',sfcmodel
        print *,'  imoist   = ',imoist
        print *
        print *,'  sfcmodel=2,3,4 requires imoist=1 (for now) '
        print *
        print *,'   stopping model .... '
        print *
        endif



        call stopcm1
      ENDIF

      IF( terrain_flag .and. output_interp.ne.0 .and. output_format.eq.2 )THEN
        if(myid.eq.0)then
        print *
        print *,'  output_interp = ',output_interp
        print *
        print *,'  output_interp=1 is not currently available for netcdf output'
        print *
        endif



        call stopcm1
      ENDIF


      IF(output_format.eq.3.or.output_format.eq.4.or.output_format.eq.5)THEN
        if(myid.eq.0)then
        print *
        print *,'  output_format = ',output_format
        print *
        print *,'  You have requested hdf output, but you have not'
        print *,'  compiled the code with hdf capability.  Modify the'
        print *,'  Makefile, clean, and recompile'
        print *
        endif



        call stopcm1
      ENDIF



!--------------------------------------------------------------
!  Check domain size (MPI only)

!--------------------------------------------------------------
!  Check that lateral bc combinations make sense:

      if(ebc.eq.1 .and. wbc.ne.1)then
        print *,"Can not have periodic b.c.'s on one side only!"



        call stopcm1
      endif
      if(wbc.eq.1 .and. ebc.ne.1)then
        print *,"Can not have periodic b.c.'s on one side only!"



        call stopcm1
      endif
      if(nbc.eq.1 .and. sbc.ne.1)then
        print *,"Can not have periodic b.c.'s on one side only!"



        call stopcm1
      endif
      if(sbc.eq.1 .and. nbc.ne.1)then
        print *,"Can not have periodic b.c.'s on one side only!"



        call stopcm1
      endif

!  end fatal checks  (ie, model stops)
!--------------------------------------------------------------
!  Some basic checks:

      ! for passive tracers:
      iptra    = max(0,min(1,iptra))
      if(iptra.eq.1)then
        npt      = max(1,npt)
      else
        npt      = 1
      endif

      ! for parcels:
      nparcels = max(1,nparcels)

!-----

      if(stretch_z.lt.1) ztop = dz*float(nk)
      IF ( stretch_z == 2 ) dz = ztop/float(nk) ! nk is the number of scalar levels

      IF( advwenos.eq.2 )THEN
        hadvordrs = 5
        vadvordrs = 5
      ENDIF
      IF( advwenos.lt.0 .or. advwenos.gt.2 )THEN
        print *
        print *,'  advwenos = ',advwenos
        print *
        print *,'  unrecognized value for advwenos '
        print *
        print *,'   stopping model .... '
        print *



        call stopcm1
      ENDIF

      !-----

      IF( advwenov.eq.2 )THEN
        hadvordrv = 5
        vadvordrv = 5
      ENDIF
      IF( advwenov.lt.0 .or. advwenov.gt.2 )THEN
        print *
        print *,'  advwenov = ',advwenov
        print *
        print *,'  unrecognized value for advwenov '
        print *
        print *,'   stopping model .... '
        print *



        call stopcm1
      ENDIF

!--------------------------------------------------------------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'dx        =',dx
      if(dowr) write(outfile,*) 'dy        =',dy
      if(dowr) write(outfile,*) 'dz        =',dz
      if(dowr) write(outfile,*) 'dtl       =',dtl
      if(dowr) write(outfile,*) 'timax     =',timax
      if(dowr) write(outfile,*) 'run_time  =',run_time
      if(dowr) write(outfile,*) 'tapfrq    =',tapfrq
      if(dowr) write(outfile,*) 'rstfrq    =',rstfrq
      if(dowr) write(outfile,*) 'statfrq   =',statfrq
      if(dowr) write(outfile,*) 'prclfrq   =',prclfrq
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'adapt_dt  =',adapt_dt
      if(dowr) write(outfile,*) 'irst      =',irst
      if(dowr) write(outfile,*) 'rstnum    =',rstnum
      if(dowr) write(outfile,*) 'iconly    =',iconly
      if(dowr) write(outfile,*) 'hadvordrs =',hadvordrs
      if(dowr) write(outfile,*) 'vadvordrs =',vadvordrs
      if(dowr) write(outfile,*) 'hadvordrv =',hadvordrv
      if(dowr) write(outfile,*) 'vadvordrv =',vadvordrv
      if(dowr) write(outfile,*) 'pdscheme  =',pdscheme
      if(dowr) write(outfile,*) 'apmasscon =',apmasscon
      if(dowr) write(outfile,*) 'advwenos  =',advwenos
      if(dowr) write(outfile,*) 'advwenov  =',advwenov
      if(dowr) write(outfile,*) 'idiff     =',idiff
      if(dowr) write(outfile,*) 'mdiff     =',mdiff
      if(dowr) write(outfile,*) 'difforder =',difforder
      if(dowr) write(outfile,*) 'imoist    =',imoist
      if(dowr) write(outfile,*) 'iturb     =',iturb
      if(dowr) write(outfile,*) 'tconfig   =',tconfig
      if(dowr) write(outfile,*) 'bcturbs   =',bcturbs
      if(dowr) write(outfile,*) 'dns       =',dns
      if(dowr) write(outfile,*) 'irdamp    =',irdamp
      if(dowr) write(outfile,*) 'hrdamp    =',hrdamp
      if(dowr) write(outfile,*) 'psolver   =',psolver
      if(dowr) write(outfile,*) 'nsound    =',nsound
      if(dowr) write(outfile,*) 'ptype     =',ptype
      if(dowr) write(outfile,*) 'ihail     =',ihail
      if(dowr) write(outfile,*) 'iautoc    =',iautoc
      if(dowr) write(outfile,*) 'icor      =',icor
      if(dowr) write(outfile,*) 'pertcor   =',pertcor
      if(dowr) write(outfile,*) 'eqtset    =',eqtset
      if(dowr) write(outfile,*) 'idiss     =',idiss
      if(dowr) write(outfile,*) 'efall     =',efall
      if(dowr) write(outfile,*) 'rterm     =',rterm
      if(dowr) write(outfile,*) 'wbc       =',wbc
      if(dowr) write(outfile,*) 'ebc       =',ebc
      if(dowr) write(outfile,*) 'sbc       =',sbc
      if(dowr) write(outfile,*) 'nbc       =',nbc
      if(dowr) write(outfile,*) 'bbc       =',bbc
      if(dowr) write(outfile,*) 'tbc       =',tbc
      if(dowr) write(outfile,*) 'irbc      =',irbc
      if(dowr) write(outfile,*) 'roflux    =',roflux
      if(dowr) write(outfile,*) 'isnd      =',isnd
      if(dowr) write(outfile,*) 'iwnd      =',iwnd
      if(dowr) write(outfile,*) 'itern     =',itern
      if(dowr) write(outfile,*) 'iinit     =',iinit
      if(dowr) write(outfile,*) 'irandp    =',irandp
      if(dowr) write(outfile,*) 'ibalance  =',ibalance
      if(dowr) write(outfile,*) 'iorigin   =',iorigin
      if(dowr) write(outfile,*) 'axisymm   =',axisymm
      if(dowr) write(outfile,*) 'imove     =',imove
      if(dowr) write(outfile,*) 'iptra     =',iptra
      if(dowr) write(outfile,*) 'npt       =',npt
      if(dowr) write(outfile,*) 'pdtra     =',pdtra
      if(dowr) write(outfile,*) 'iprcl     =',iprcl
      if(dowr) write(outfile,*) 'nparcels  =',nparcels
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'kdiff2    =',kdiff2
      if(dowr) write(outfile,*) 'kdiff6    =',kdiff6
      if(dowr) write(outfile,*) 'fcor      =',fcor
      if(dowr) write(outfile,*) 'kdiv      =',kdiv
      if(dowr) write(outfile,*) 'alph      =',alph
      if(dowr) write(outfile,*) 'rdalpha   =',rdalpha
      if(dowr) write(outfile,*) 'zd        =',zd
      if(dowr) write(outfile,*) 'xhd       =',xhd
      if(dowr) write(outfile,*) 'umove     =',umove
      if(dowr) write(outfile,*) 'vmove     =',vmove
      if(dowr) write(outfile,*) 'v_t       =',v_t
      if(dowr) write(outfile,*) 'l_h       =',l_h
      if(dowr) write(outfile,*) 'lhref1    =',lhref1
      if(dowr) write(outfile,*) 'lhref2    =',lhref2
      if(dowr) write(outfile,*) 'l_inf     =',l_inf
      if(dowr) write(outfile,*) 'ndcnst    =',ndcnst
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'radopt    =',radopt
      if(dowr) write(outfile,*) 'dtrad     =',dtrad
      if(dowr) write(outfile,*) 'ctrlat    =',ctrlat
      if(dowr) write(outfile,*) 'ctrlon    =',ctrlon
      if(dowr) write(outfile,*) 'year      =',year
      if(dowr) write(outfile,*) 'month     =',month
      if(dowr) write(outfile,*) 'day       =',day
      if(dowr) write(outfile,*) 'hour      =',hour
      if(dowr) write(outfile,*) 'minute    =',minute
      if(dowr) write(outfile,*) 'second    =',second
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'isfcflx   =',isfcflx
      if(dowr) write(outfile,*) 'sfcmodel  =',sfcmodel
      if(dowr) write(outfile,*) 'oceanmodel=',oceanmodel
      if(dowr) write(outfile,*) 'ipbl      =',ipbl
      if(dowr) write(outfile,*) 'initsfc   =',initsfc
      if(dowr) write(outfile,*) 'tsk0      =',tsk0
      if(dowr) write(outfile,*) 'tmn0      =',tmn0
      if(dowr) write(outfile,*) 'xland0    =',xland0
      if(dowr) write(outfile,*) 'lu0       =',lu0
      if(dowr) write(outfile,*) 'season    =',season
      if(dowr) write(outfile,*) 'cecd      =',cecd
      if(dowr) write(outfile,*) 'pertflx   =',pertflx
      if(dowr) write(outfile,*) 'cnstce    =',cnstce
      if(dowr) write(outfile,*) 'cnstcd    =',cnstcd
      if(dowr) write(outfile,*) 'isftcflx  =',isftcflx
      if(dowr) write(outfile,*) 'iz0tlnd   =',iz0tlnd
      if(dowr) write(outfile,*) 'oml_hml0  =',oml_hml0
      if(dowr) write(outfile,*) 'oml_gamma =',oml_gamma
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'stretch_x =',stretch_x
      if(dowr) write(outfile,*) 'dx_inner  =',dx_inner
      if(dowr) write(outfile,*) 'dx_outer  =',dx_outer
      if(dowr) write(outfile,*) 'nos_x_len =',nos_x_len
      if(dowr) write(outfile,*) 'tot_x_len =',tot_x_len
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'stretch_y =',stretch_y
      if(dowr) write(outfile,*) 'dy_inner  =',dy_inner
      if(dowr) write(outfile,*) 'dy_outer  =',dy_outer
      if(dowr) write(outfile,*) 'nos_y_len =',nos_y_len
      if(dowr) write(outfile,*) 'tot_y_len =',tot_y_len
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'stretch_z =',stretch_z
      if(dowr) write(outfile,*) 'ztop      =',ztop
      if(dowr) write(outfile,*) 'str_bot   =',str_bot
      if(dowr) write(outfile,*) 'str_top   =',str_top
      if(dowr) write(outfile,*) 'dz_bot    =',dz_bot
      if(dowr) write(outfile,*) 'dz_top    =',dz_top
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'bc_temp   =',bc_temp
      if(dowr) write(outfile,*) 'ptc_top   =',ptc_top
      if(dowr) write(outfile,*) 'ptc_bot   =',ptc_bot
      if(dowr) write(outfile,*) 'viscosity =',viscosity
      if(dowr) write(outfile,*) 'pr_num    =',pr_num
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'var1      =',var1
      if(dowr) write(outfile,*) 'var2      =',var2
      if(dowr) write(outfile,*) 'var3      =',var3
      if(dowr) write(outfile,*) 'var4      =',var4
      if(dowr) write(outfile,*) 'var5      =',var5
      if(dowr) write(outfile,*) 'var6      =',var6
      if(dowr) write(outfile,*) 'var7      =',var7
      if(dowr) write(outfile,*) 'var8      =',var8
      if(dowr) write(outfile,*) 'var9      =',var9
      if(dowr) write(outfile,*) 'var10     =',var10
      if(dowr) write(outfile,*)


      IF ( ptype >= 26 .and. dowr ) THEN
        write(outfile,NML=nssl2mom_params)
!        write(outfile,*) 'alphah    =',alphah
!        write(outfile,*) 'alphahl   =',alphahl
!        write(outfile,*) 'dfrz      =',dfrz
!        write(outfile,*) 'hldnmn    =',hldnmn
!        write(outfile,*) 'imurain   =',imurain
!        write(outfile,*) 'ccn       =',ccn
!        write(outfile,*) 'icdx      =',icdx
!        write(outfile,*) 'icdxhl    =',icdxhl
!        write(outfile,*) 'iferwisventr =',iferwisventr
!        write(outfile,*) 'iehw      =',iehw
!        write(outfile,*) 'iehlw     =',iehlw
!        write(outfile,*) 'ehw0      =',ehw0
!        write(outfile,*) 'ehlw0     =',ehlw0
!        write(outfile,*) 'dmrauto   =',dmrauto
!        write(outfile,*) 'ioldlimiter=',ioldlimiter
      ENDIF


!--------------------------------------------------------------
!  Configuration for simulations with moisture
!

      !--- begin: define defaults (please do not change) ---------
      iice     = 0
      idm      = 0
      idmplus  = 0
      numq     = 1
      nqv      = 1
      nql1     = 1
      nql2     = 1
      nqs1     = 1
      nqs2     = 1
      nnc1     = 1
      nnc2     = 1
      nzl1     = 1
      nzl2     = 1
      nvl1     = 1
      nvl2     = 1
      nbudget  = 10
      budrain  = 1
      cloudvar = .false.
      rhovar   = .false.
      !--- end: define defaults ----------------------------------

      IF(imoist.eq.1)THEN

!-----------------------------------------------------------------------
!-------   BEGIN:  modify stuff below here -----------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
        IF(ptype.eq.1)THEN        ! Kessler scheme

          numq = 3    ! there are 3 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array

          cloudvar(1) = .false.
          cloudvar(2) = .true.
          cloudvar(3) = .false.

          qname(1) = 'qv '
          qname(2) = 'qc '
          qname(3) = 'qr '

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... using Kessler microphysics scheme ... '
          if(dowr) write(outfile,*) '         numq   = ',numq
          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
        ELSEIF((ptype.eq.2).or.(ptype.eq.4))THEN    ! Goddard-LFO or 
                                                    ! GSR-LFO scheme

          iice = 1    ! this means that ptype=2,4 are ice schemes

          numq = 6    ! there are 6 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array

          cloudvar(1) = .false.
          cloudvar(2) = .true.
          cloudvar(3) = .false.
          cloudvar(4) = .true.
          cloudvar(5) = .false.
          cloudvar(6) = .false.

          qname(1) = 'qv '
          qname(2) = 'qc '
          qname(3) = 'qr '
          qname(4) = 'qi '
          qname(5) = 'qs '
          qname(6) = 'qg '

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the Goddard or GSR LFO scheme -----

          if(ptype.eq.2)THEN

            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) 'Calling CONSAT'
            if(dowr) write(outfile,*)

            call consat
            call consat2(dtl)

            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) '  ... using Goddard LFO microphysics scheme ... '
            if(dowr) write(outfile,*) '         numq   = ',numq
            if(dowr) write(outfile,*) '         ihail  = ',ihail
            if(dowr) write(outfile,*)

          endif

          if(ptype.eq.4)then

            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) 'Calling lfoice_init'
            if(dowr) write(outfile,*)

            call lfoice_init(dtl)

            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) '  ... using GSR LFO microphysics scheme ... '
            if(dowr) write(outfile,*) '         numq   = ',numq
            if(dowr) write(outfile,*)

          endif

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
        ELSEIF(ptype.eq.3)THEN    ! Thompson scheme

          iice = 1    ! this means that ptype=3 is an ice scheme
          idm  = 1    ! this means that ptype=3 has at least one double moment

          numq = 8    ! there are 8 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array
          nnc1 = 7    ! the first number concentration var is the seventh array
          nnc2 = 8    ! the last number concentration var is the eighth array

          cloudvar(1) = .false.
          cloudvar(2) = .true.
          cloudvar(3) = .false.
          cloudvar(4) = .true.
          cloudvar(5) = .false.
          cloudvar(6) = .false.
          cloudvar(7) = .false.
          cloudvar(8) = .false.

          qname(1) = 'qv '
          qname(2) = 'qc '
          qname(3) = 'qr '
          qname(4) = 'qi '
          qname(5) = 'qs '
          qname(6) = 'qg '
          qname(7) = 'nci'
          qname(8) = 'ncr'

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the Thompson scheme -----

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) 'Calling thompson_init'
          if(dowr) write(outfile,*) '(this can take several minutes ... please be patient)'

          call thompson_init

          if(dowr) write(outfile,*) 'Done with thompson_init'
          if(dowr) write(outfile,*)

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... using Thompson microphysics scheme ... '
          if(dowr) write(outfile,*) '         numq   = ',numq
          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        ELSEIF(ptype.eq.5)THEN    ! Morrison scheme

          !----- initialize the Morrison scheme -----

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) 'Calling GRAUPEL_INIT'

          call graupel_init(ihail,inum,ndcnst)

          if(dowr) write(outfile,*) 'Returned from GRAUPEL_INIT'
          if(dowr) write(outfile,*)

          !------------------------------------------

          iice = 1    ! this means that ptype=5 is an ice scheme
          idm  = 1    ! this means that ptype=5 has at least one double moment

        if(inum.eq.1)then
          ! constant cloud-drop concentration

          numq = 10   ! there are 10 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array
          nnc1 = 7    ! the first number concentration var is the seventh array
          nnc2 = 10   ! the last number concentration var is the tenth array

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .false.
          cloudvar( 6) = .false.
          cloudvar( 7) = .false.
          cloudvar( 8) = .false.
          cloudvar( 9) = .false.
          cloudvar(10) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qs '
          qname( 6) = 'qg '
          qname( 7) = 'nci'
          qname( 8) = 'ncs'
          qname( 9) = 'ncr'
          qname(10) = 'ncg'

        elseif(inum.eq.0)then
          ! cloud-droplet concentration is a predicted variable

          numq = 11   ! there are 11 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array
          nnc1 = 7    ! the first number concentration var is the seventh array
          nnc2 = 11   ! the last number concentration var is the eleventh array

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .false.
          cloudvar( 6) = .false.
          cloudvar( 7) = .false.
          cloudvar( 8) = .false.
          cloudvar( 9) = .false.
          cloudvar(10) = .false.
          cloudvar(11) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qs '
          qname( 6) = 'qg '
          qname( 7) = 'nci'
          qname( 8) = 'ncs'
          qname( 9) = 'ncr'
          qname(10) = 'ncg'
          qname(11) = 'ncc'

        else

          print *,'  unrecognized value for inum '



          call stopcm1

        endif

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... using Morrison microphysics scheme ... '
          if(dowr) write(outfile,*) '         numq   = ',numq
          if(dowr) write(outfile,*) '         ihail  = ',ihail
        if(inum.eq.1)then
          if(dowr) write(outfile,*) '         assuming constant cloud droplet concentration' 
          if(dowr) write(outfile,*) '         ndcnst = ',ndcnst
        else
          if(dowr) write(outfile,*) '         predicting cloud droplet concentration'
        endif
          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        ELSEIF(ptype.eq.6)THEN        ! Rotunno-Emanuel scheme

          numq = 2    ! there are 2 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 2    ! the last liquid variable is the second array

          cloudvar(1) = .false.
          cloudvar(2) = .true.

          qname(1) = 'qv '
          qname(2) = 'ql '

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... using Rotunno-Emanuel microphysics scheme ... '
          if(dowr) write(outfile,*) '         numq   = ',numq
          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        ELSEIF ( ptype .eq. 26 ) THEN    ! ZVD scheme (with no hail category)

          iice = 1    ! this means that ptype=26 is an ice scheme
          idm  = 1    ! this means that ptype=26 has at least one double moment
          idmplus = 1

          numq = 13   ! number of variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array
          nnc1 = 7    ! the first number concentration var is the seventh array
          nnc2 = 12   ! the last number concentration var is the eleventh array
          nvl1 = 13   ! the first partical volume var is the seventh array
          nvl2 = 13   ! the last partical volume var is the eleventh array

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .false.
          cloudvar( 6) = .false.
          cloudvar( 7) = .false.
          cloudvar( 8) = .false.
          cloudvar( 9) = .false.
          cloudvar(10) = .false.
          cloudvar(11) = .false.
          cloudvar(12) = .false.
          cloudvar(13) = .false.
!          cloudvar(14) = .false.
!          cloudvar(15) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qs '
          qname( 6) = 'qg '
!          qname( 7) = 'qhl'
          qname( 7) = 'ccn' ! CCN concentration
          qname( 8) = 'ccw' ! droplet conc
          qname( 9) = 'crw' ! rain conc
          qname(10) = 'cci' ! ice crystal conc
          qname(11) = 'csw' ! snow conc
          qname(12) = 'chw' ! graupel conc
!          qname(14) = 'chl' ! hail conc
          qname(13) = 'vhw' ! graupel volume

          rhovar( 1) = .false.
          rhovar( 2) = .false.
          rhovar( 3) = .false.
          rhovar( 4) = .false.
          rhovar( 5) = .false.
          rhovar( 6) = .false.
          rhovar( 7) = .true.
          rhovar( 8) = .true.
          rhovar( 9) = .true.
          rhovar(10) = .true.
          rhovar(11) = .true.
          rhovar(12) = .true.
!          rhovar(13) = .false.
          rhovar(13) = .true.
!          rhovar(14) = .true.
!          rhovar(15) = .true.

!          ipconc = 5
!          lr = 4
!          li = 5
!          ls = 6
!          lh = 7
!          lg = lh
!          lhab = lh
!          lhl = 0
!          lqe  = lhab
!
!          lccn = 8
!          lnc  = 9
!          lnr  = 10
!          lni  = 11
!          lns  = 12
!          lnh  = 13
!          lnhl = 0
!          lss  = 14
!          lvh  = 15
!
!          lsch = 0
!          lschab = 0
!          lscw = 0
!          lscb = lscw
!          lscni = 0
!          lscpi = 0
!          lsce = lscni
!          lsceq= lschab
!
!          lsw  = 0
!          lhw  = 0
!          lhlw = 0

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the ZVD scheme -----

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) 'Calling index_module_init'
          if(dowr) write(outfile,*)

!          write(0,*) 'option 26 currently not available'
!          STOP
         CALL nssl_2mom_init(ipctmp=5,mixphase=0,ihvol=-1,eqtset_tmp=eqtset)
!          call INDEX_MODULE_INIT(ptype)

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) 'Returned from index_module_init'
          if(dowr) write(outfile,*)

        ELSEIF ( ptype .eq. 27 ) THEN    ! ZVDH scheme (with hail)

          iice = 1    ! this means that ptype=27 is an ice scheme
          idm  = 1    ! this means that ptype=27 has at least one double moment
          idmplus = 1

          numq = 16   ! number of variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 7    ! the last solid variable is the sixth array
          nnc1 = 8    ! the first number concentration var is the seventh array
          nnc2 = 14   ! the last number concentration var is the eleventh array
          nvl1 = 15   ! the first partical volume var is the seventh array
          nvl2 = 16   ! the last partical volume var is the eleventh array

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .false.
          cloudvar( 6) = .false.
          cloudvar( 7) = .false.
          cloudvar( 8) = .false.
          cloudvar( 9) = .false.
          cloudvar(10) = .false.
          cloudvar(11) = .false.
          cloudvar(12) = .false.
          cloudvar(13) = .false.
          cloudvar(14) = .false.
          cloudvar(15) = .false.
          cloudvar(16) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qs '
          qname( 6) = 'qg '
          qname( 7) = 'qhl'
          qname( 8) = 'ccn' ! CCN concentration
          qname( 9) = 'ccw' ! droplet conc
          qname(10) = 'crw' ! rain conc
          qname(11) = 'cci' ! ice crystal conc
          qname(12) = 'csw' ! snow conc
          qname(13) = 'chw' ! graupel conc
          qname(14) = 'chl' ! hail conc
!          qname(15) = 'ss ' ! max supersaturation
          qname(15) = 'vhw' ! graupel volume
          qname(16) = 'vhl' ! hail volume

          rhovar( 1) = .false.
          rhovar( 2) = .false.
          rhovar( 3) = .false.
          rhovar( 4) = .false.
          rhovar( 5) = .false.
          rhovar( 6) = .false.
          rhovar( 7) = .false.
          rhovar( 8) = .true.
          rhovar( 9) = .true.
          rhovar(10) = .true.
          rhovar(11) = .true.
          rhovar(12) = .true.
          rhovar(13) = .true.
          rhovar(14) = .true.
!          rhovar(15) = .false.
          rhovar(15) = .true.
          rhovar(16) = .true.

!          ipconc = 5
!          lr = 4
!          li = 5
!          ls = 6
!          lh = 7
!          lhl = 8
!          lg = lh
!          lhab = lhl
!          lqe  = lhab
!
!          lccn = 9
!          lnc  = 10
!          lnr  = 11
!          lni  = 12
!          lns  = 13
!          lnh  = 14
!          lnhl = 15
!          lss  = 16
!          lvh  = 17
!          lvhl = 18
!
!          lsch = 0
!          lschab = 0
!          lscw = 0
!          lscb = lscw
!          lscni = 0
!          lscpi = 0
!          lsce = lscni
!          lsceq= lschab
!
!          lsw  = 0
!          lhw  = 0
!          lhlw = 0

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the ZVD scheme -----

!          if(dowr) write(outfile,*)
!          if(dowr) write(outfile,*) 'Calling graupel_init'
!          if(dowr) write(outfile,*)

!          call INDEX_MODULE_INIT(ptype)
         CALL nssl_2mom_init(ipctmp=5,mixphase=0,ihvol=1,eqtset_tmp=eqtset)

!          if(dowr) write(outfile,*)
!          if(dowr) write(outfile,*) 'Returned from graupel_init'
!          if(dowr) write(outfile,*)

        ELSEIF ( ptype .eq. 28 ) THEN    ! single moment ZIEG scheme (without hail)

          iice = 1    ! this means that ptype=28 is an ice scheme
          idm  = 0    ! this means that ptype=28 has at least one double moment
          idmplus = 1

          numq = 6   ! number of variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .false.
          cloudvar( 6) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qs '
          qname( 6) = 'qg '

          rhovar( 1) = .false.
          rhovar( 2) = .false.
          rhovar( 3) = .false.
          rhovar( 4) = .false.
          rhovar( 5) = .false.
          rhovar( 6) = .false.

!          ipconc = 0
!          lr = 4
!          li = 5
!          ls = 6
!          lh = 7
!          lg = lh
!          lhab = lh
!          lhl = 0
!          lqe  = lhab
!
!          lccn = 0
!          lnc  = 0
!          lnr  = 0
!          lni  = 0
!          lns  = 0
!          lnh  = 0
!          lnhl = 0
!          lss  = 0
!          lvh  = 0
!
!          lsch = 0
!          lschab = 0
!          lscw = 0
!          lscb = lscw
!          lscni = 0
!          lscpi = 0
!          lsce = lscni
!          lsceq= lschab
!
!          lsw  = 0
!          lhw  = 0
!          lhlw = 0

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the ZVD scheme -----

!          if(dowr) write(outfile,*)
!          if(dowr) write(outfile,*) 'Calling graupel_init'
!          if(dowr) write(outfile,*)

         CALL nssl_2mom_init(ipctmp=0,mixphase=0,ihvol=-1,eqtset_tmp=eqtset)


!          if(dowr) write(outfile,*)
!          if(dowr) write(outfile,*) 'Returned from graupel_init'
!          if(dowr) write(outfile,*)

        ELSEIF ( ptype .eq. 29 ) THEN    ! 3-moment ZVDH scheme (with hail)

          iice = 1    ! this means that ptype=27 is an ice scheme
          idm  = 1    ! this means that ptype=27 has at least one double moment
          idmplus = 1

          numq = 19   ! number of variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 7    ! the last solid variable is the sixth array
          nnc1 = 8    ! the first number concentration var is the seventh array
          nnc2 = 14   ! the last number concentration var is the eleventh array
          nzl1 = 15   ! the first reflectivity var is the seventh array
          nzl2 = 17   ! the last reflectivity var is the eleventh array
          nvl1 = 18   ! the first partical volume var is the seventh array
          nvl2 = 19   ! the last partical volume var is the eleventh array

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .false.
          cloudvar( 6) = .false.
          cloudvar( 7) = .false.
          cloudvar( 8) = .false.
          cloudvar( 9) = .false.
          cloudvar(10) = .false.
          cloudvar(11) = .false.
          cloudvar(12) = .false.
          cloudvar(13) = .false.
          cloudvar(14) = .false.
          cloudvar(15) = .false.
          cloudvar(16) = .false.
          cloudvar(17) = .false.
          cloudvar(18) = .false.
          cloudvar(19) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qs '
          qname( 6) = 'qg '
          qname( 7) = 'qhl'
          qname( 8) = 'ccn' ! CCN concentration
          qname( 9) = 'ccw' ! droplet conc
          qname(10) = 'crw' ! rain conc
          qname(11) = 'cci' ! ice crystal conc
          qname(12) = 'csw' ! snow conc
          qname(13) = 'chw' ! graupel conc
          qname(14) = 'chl' ! hail conc
!          qname(15) = 'ss ' ! max supersaturation
          qname(15) = 'zrw' ! rain reflectivity moment
          qname(16) = 'zhw' ! graupel reflectivity moment
          qname(17) = 'zhl' ! hail reflectivity moment
          qname(18) = 'vhw' ! graupel volume
          qname(19) = 'vhl' ! hail volume

          rhovar( 1) = .false.
          rhovar( 2) = .false.
          rhovar( 3) = .false.
          rhovar( 4) = .false.
          rhovar( 5) = .false.
          rhovar( 6) = .false.
          rhovar( 7) = .false.
          rhovar( 8) = .true.
          rhovar( 9) = .true.
          rhovar(10) = .true.
          rhovar(11) = .true.
          rhovar(12) = .true.
          rhovar(13) = .true.
          rhovar(14) = .true.
!          rhovar(15) = .false.
          rhovar(15) = .true.
          rhovar(16) = .true.
          rhovar(17) = .true.
          rhovar(18) = .true.
          rhovar(19) = .true.


          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the ZVD scheme -----

!          if(dowr) write(outfile,*)
!          if(dowr) write(outfile,*) 'Calling graupel_init'
!          if(dowr) write(outfile,*)

!          call INDEX_MODULE_INIT(ptype)
         CALL nssl_2mom_init(ipctmp=8,mixphase=0,ihvol=1,eqtset_tmp=eqtset)

!          if(dowr) write(outfile,*)
!          if(dowr) write(outfile,*) 'Returned from graupel_init'
!          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  insert new ptype here

!!!        ELSEIF(ptype.eq.8)THEN    ! new microphysics scheme

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        ELSE

          IF(myid.eq.0)THEN
            print *
            print *,'  ptype = ',ptype
            print *
            print *,'  Unrecognized value for ptype '
            print *
            print *,'  ... stopping cm1 ... '
            print *
          ENDIF




          call stopcm1

        ENDIF    ! endif for ptype

      ENDIF    ! endif for imoist=1

!-----------------------------------------------------------------------
!-------   END:  modify stuff above here -------------------------------
!-----------------------------------------------------------------------

      IF( (radopt.eq.1.or.radopt.eq.2) .and. iice.ne.1 )THEN
        print *
        print *,'  radopt   = ',radopt
        print *,'  iice     = ',iice
        print *
        print *,'  radopt=1 requires an ice microphysics scheme '
        print *
        print *,'   stopping model .... '
        print *



        call stopcm1
      ENDIF

!-----------------------------------------------------------------------

      nqc = 0
      nqr = 0
      nqi = 0
      nqs = 0
      nqg = 0

      do n=1,numq
        if( qname(n).eq.'qc ' .or. qname(n).eq.'ql ' ) nqc = n
        if( qname(n).eq.'qr ' ) nqr = n
        if( qname(n).eq.'qi ' ) nqi = n
        if( qname(n).eq.'qs ' ) nqs = n
        if( qname(n).eq.'qg ' ) nqg = n
      enddo

      if(numq .gt. maxq)then
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  WARNING!   numq > maxq'
        if(dowr) write(outfile,*) '  You need to increase maxq in input.incl and recompile'
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Stopping model ....'
        if(dowr) write(outfile,*)



        call stopcm1
      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'iice      =',iice
      if(dowr) write(outfile,*) 'idm       =',idm
      if(dowr) write(outfile,*) 'idmplus   =',idmplus
      if(dowr) write(outfile,*) 'numq      =',numq
      if(dowr) write(outfile,*) 'nqv       =',nqv
      if(dowr) write(outfile,*) 'nqc       =',nqc
      if(dowr) write(outfile,*) 'nqr       =',nqr
      if(dowr) write(outfile,*) 'nqi       =',nqi
      if(dowr) write(outfile,*) 'nqs       =',nqs
      if(dowr) write(outfile,*) 'nqg       =',nqg
      if(dowr) write(outfile,*) 'nql1      =',nql1
      if(dowr) write(outfile,*) 'nql2      =',nql2
      if(dowr) write(outfile,*) 'nqs1      =',nqs1
      if(dowr) write(outfile,*) 'nqs2      =',nqs2
      if(dowr) write(outfile,*) 'nnc1      =',nnc1
      if(dowr) write(outfile,*) 'nnc2      =',nnc2
      if(dowr) write(outfile,*) 'nvl1      =',nvl1
      if(dowr) write(outfile,*) 'nvl2      =',nvl2
      if(dowr) write(outfile,*) 'nzl1      =',nzl1
      if(dowr) write(outfile,*) 'nzl2      =',nzl2
      if(dowr) write(outfile,*)

!--------------------------------------------------------------

      iterrain = 0
      if(terrain_flag) iterrain = 1

      output_interp  = max(0,min(1,output_interp))*iterrain
      output_rain    = max(0,min(1,output_rain))
      output_sws     = max(0,min(1,output_sws))
      output_svs     = max(0,min(1,output_svs))
      output_sps     = max(0,min(1,output_sps))
      output_srs     = max(0,min(1,output_srs))
      output_sgs     = max(0,min(1,output_sgs))
      output_sus     = max(0,min(1,output_sus))
      output_shs     = max(0,min(1,output_shs))
      output_coldpool= max(0,min(1,output_coldpool))
      output_sfcflx  = max(0,min(1,output_sfcflx))
      output_sfcparams = max(0,min(1,output_sfcparams))
      output_sfcdiags = max(0,min(1,output_sfcdiags))
      output_psfc    = max(0,min(1,output_psfc))
      output_zs      = max(0,min(1,output_zs))*iterrain
      output_zh      = max(0,min(1,output_zh))
      output_basestate = max(0,min(1,output_basestate))
      output_th      = max(0,min(1,output_th))
      output_thpert  = max(0,min(1,output_thpert))
      output_prs     = max(0,min(1,output_prs))
      output_prspert = max(0,min(1,output_prspert))
      output_pi      = max(0,min(1,output_pi))
      output_pipert  = max(0,min(1,output_pipert))
      output_rho     = max(0,min(1,output_rho))
      output_rhopert = max(0,min(1,output_rhopert))
      output_tke     = max(0,min(1,output_tke))
      output_km      = max(0,min(1,output_km))
      output_kh      = max(0,min(1,output_kh))
      output_qv      = max(0,min(1,output_qv))
      output_qvpert  = max(0,min(1,output_qvpert))
      output_q       = max(0,min(1,output_q))
      output_dbz     = max(0,min(1,output_dbz))
      output_buoyancy= max(0,min(1,output_buoyancy))
      output_u       = max(0,min(1,output_u))
      output_upert   = max(0,min(1,output_upert))
      output_uinterp = max(0,min(1,output_uinterp))
      output_v       = max(0,min(1,output_v))
      output_vpert   = max(0,min(1,output_vpert))
      output_vinterp = max(0,min(1,output_vinterp))
      output_w       = max(0,min(1,output_w))
      output_winterp = max(0,min(1,output_winterp))
      output_vort    = max(0,min(1,output_vort))
      output_pv      = max(0,min(1,output_pv))
      output_uh      = max(0,min(1,output_uh))
      output_pblten  = max(0,min(1,output_pblten))
      output_dissten = max(0,min(1,output_dissten))
      output_dissheat = max(0,min(1,output_dissheat))
      output_mptend   = max(0,min(1,output_mptend))
      output_fallvel  = max(0,min(1,output_fallvel))
      output_nm      = max(0,min(1,output_nm))
      output_def     = max(0,min(1,output_def))
      output_turbten = max(0,min(1,output_turbten))
      output_impdiften = max(0,min(1,output_impdiften))
      output_radten  = max(0,min(1,output_radten))


      nrain = 1
      if(imove.eq.1) nrain = 2

      if(dowr) write(outfile,*) 'nrain     =',nrain
      if(dowr) write(outfile,*)

      if(imoist.eq.0)then
        output_rain=0
        output_srs=0
        output_sgs=0
        output_qv=0
        output_qvpert=0
        output_q=0
        output_dbz=0
      endif
      if( nqr.le.0 ) output_srs = 0
      if( nqg.le.0 ) output_sgs = 0
      if( (iturb.eq.0.or.dns.eq.1) )then
        output_tke=0
      endif
      if( (iturb.eq.0.or.dns.eq.1).and.(ipbl.eq.0) )then
        output_km=0
        output_kh=0
      endif
      if(iturb.eq.2.or.iturb.eq.3)then
        output_tke=0
      endif
      if(ipbl.lt.1)then
        output_pblten=0
      endif
      if(radopt.eq.0)then
        output_radten=0
      endif
      if( iturb.eq.0 )then
        output_turbten = 0
      endif
      if( output_impdiften.eq.1 )then
        if( hadvordrv.ne.5 .and. vadvordrv.ne.5 ) output_impdiften = 0
      endif
      if( bbc.ne.3 .and. isfcflx.eq.0 )then
!!!        output_sfcflx = 0
        output_sfcparams = 0
        output_sfcdiags = 0
      endif
      if(terrain_flag)then
        output_zs = 1
        output_zh = 1
      endif
      if( imoist.eq.0 )then
        output_mptend = 0
      endif
      if( imoist.eq.0 .or. ptype.ne.5 )then
        output_fallvel = 0
      endif
      if( idiss.eq.0 )then
        output_dissheat = 0
      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'output_path      = ',output_path
      if(dowr) write(outfile,*) 'output_basename  = ',output_basename
      if(dowr) write(outfile,*) 'output_format    =',output_format
      if(dowr) write(outfile,*) 'output_filetype  =',output_filetype
      if(dowr) write(outfile,*) 'output_interp    =',output_interp
      if(dowr) write(outfile,*) 'output_rain      =',output_rain
      if(dowr) write(outfile,*) 'output_sws       =',output_sws
      if(dowr) write(outfile,*) 'output_svs       =',output_svs
      if(dowr) write(outfile,*) 'output_sps       =',output_sps
      if(dowr) write(outfile,*) 'output_srs       =',output_srs
      if(dowr) write(outfile,*) 'output_sgs       =',output_sgs
      if(dowr) write(outfile,*) 'output_sus       =',output_sus
      if(dowr) write(outfile,*) 'output_shs       =',output_shs
      if(dowr) write(outfile,*) 'output_coldpool  =',output_coldpool
      if(dowr) write(outfile,*) 'output_sfcflx    =',output_sfcflx
      if(dowr) write(outfile,*) 'output_sfcparams =',output_sfcparams
      if(dowr) write(outfile,*) 'output_sfcdiags  =',output_sfcdiags
      if(dowr) write(outfile,*) 'output_psfc      =',output_psfc
      if(dowr) write(outfile,*) 'output_zs        =',output_zs
      if(dowr) write(outfile,*) 'output_zh        =',output_zh
      if(dowr) write(outfile,*) 'output_basestate =',output_basestate
      if(dowr) write(outfile,*) 'output_th        =',output_th
      if(dowr) write(outfile,*) 'output_thpert    =',output_thpert
      if(dowr) write(outfile,*) 'output_prs       =',output_prs
      if(dowr) write(outfile,*) 'output_prspert   =',output_prspert
      if(dowr) write(outfile,*) 'output_pi        =',output_pi
      if(dowr) write(outfile,*) 'output_pipert    =',output_pipert
      if(dowr) write(outfile,*) 'output_rho       =',output_rho
      if(dowr) write(outfile,*) 'output_rhopert   =',output_rhopert
      if(dowr) write(outfile,*) 'output_tke       =',output_tke
      if(dowr) write(outfile,*) 'output_km        =',output_km
      if(dowr) write(outfile,*) 'output_kh        =',output_kh
      if(dowr) write(outfile,*) 'output_qv        =',output_qv
      if(dowr) write(outfile,*) 'output_qvpert    =',output_qvpert
      if(dowr) write(outfile,*) 'output_q         =',output_q
      if(dowr) write(outfile,*) 'output_dbz       =',output_dbz
      if(dowr) write(outfile,*) 'output_buoyancy  =',output_buoyancy
      if(dowr) write(outfile,*) 'output_u         =',output_u
      if(dowr) write(outfile,*) 'output_upert     =',output_upert
      if(dowr) write(outfile,*) 'output_uinterp   =',output_uinterp
      if(dowr) write(outfile,*) 'output_v         =',output_v
      if(dowr) write(outfile,*) 'output_vpert     =',output_vpert
      if(dowr) write(outfile,*) 'output_vinterp   =',output_vinterp
      if(dowr) write(outfile,*) 'output_w         =',output_w
      if(dowr) write(outfile,*) 'output_winterp   =',output_winterp
      if(dowr) write(outfile,*) 'output_vort      =',output_vort
      if(dowr) write(outfile,*) 'output_pv        =',output_pv
      if(dowr) write(outfile,*) 'output_uh        =',output_uh
      if(dowr) write(outfile,*) 'output_pblten    =',output_pblten
      if(dowr) write(outfile,*) 'output_dissten   =',output_dissten
      if(dowr) write(outfile,*) 'output_dissheat  =',output_dissheat
      if(dowr) write(outfile,*) 'output_mptend    =',output_mptend
      if(dowr) write(outfile,*) 'output_fallvel   =',output_fallvel
      if(dowr) write(outfile,*) 'output_nm        =',output_nm
      if(dowr) write(outfile,*) 'output_def       =',output_def
      if(dowr) write(outfile,*) 'output_turbten   =',output_turbten
      if(dowr) write(outfile,*) 'output_impdiften =',output_impdiften
      if(dowr) write(outfile,*) 'output_radten    =',output_radten
      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'restart_format       = ',restart_format
      if(dowr) write(outfile,*) 'restart_filetype     = ',restart_filetype
      if(dowr) write(outfile,*) 'restart_file_theta   = ',restart_file_theta
      if(dowr) write(outfile,*) 'restart_file_dbz     = ',restart_file_dbz
      if(dowr) write(outfile,*) 'restart_file_th0     = ',restart_file_th0
      if(dowr) write(outfile,*) 'restart_file_prs0    = ',restart_file_prs0
      if(dowr) write(outfile,*) 'restart_file_pi0     = ',restart_file_pi0
      if(dowr) write(outfile,*) 'restart_file_rho0    = ',restart_file_rho0
      if(dowr) write(outfile,*) 'restart_file_qv0     = ',restart_file_qv0
      if(dowr) write(outfile,*) 'restart_file_u0      = ',restart_file_u0
      if(dowr) write(outfile,*) 'restart_file_v0      = ',restart_file_v0
      if(dowr) write(outfile,*) 'restart_file_zs      = ',restart_file_zs
      if(dowr) write(outfile,*) 'restart_file_zh      = ',restart_file_zh
      if(dowr) write(outfile,*) 'restart_file_zf      = ',restart_file_zf
      if(dowr) write(outfile,*) 'restart_file_diags   = ',restart_file_diags
      if(dowr) write(outfile,*) 'restart_use_theta    = ',restart_use_theta
      if(dowr) write(outfile,*) 'restart_reset_frqtim = ',restart_reset_frqtim
      if(dowr) write(outfile,*)

      maxk = nk

      if(dowr) write(outfile,*) '  maxk = ',maxk
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  Define dimensions for allocatable arrays

      if(imoist.eq.1)then
        ibm=ib
        iem=ie
        jbm=jb
        jem=je
        kbm=kb
        kem=ke
        if(ptype.ge.26)then
          ibzvd=ib
          iezvd=ie
          jbzvd=jb
          jezvd=je
          kbzvd=kb
          kezvd=ke
          nqzvd = numq + 1
        else
          ibzvd=1
          iezvd=1
          jbzvd=1
          jezvd=1
          kbzvd=1
          kezvd=1
          nqzvd=1
        endif
      else
        ibm=1
        iem=1
        jbm=1
        jem=1
        kbm=1
        kem=1
        ibzvd=1
        iezvd=1
        jbzvd=1
        jezvd=1
        kbzvd=1
        kezvd=1
        nqzvd=1
      endif

      if(iice.eq.1)then
        ibi=ib
        iei=ie
        jbi=jb
        jei=je
        kbi=kb
        kei=ke
      else
        ibi=1
        iei=1
        jbi=1
        jei=1
        kbi=1
        kei=1
      endif

      if(radopt.ge.1)then
        ibr=ib
        ier=ie
        jbr=jb
        jer=je
        kbr=kb
        ker=ke
      else
        ibr=1
        ier=1
        jbr=1
        jer=1
        kbr=1
        ker=1
      endif

      if(ipbl.ge.1)then
        ibb=ib
        ieb=ie
        jbb=jb
        jeb=je
        kbb=kb
        keb=ke
      else
        ibb=1
        ieb=1
        jbb=1
        jeb=1
        kbb=1
        keb=1
      endif

      if( (sfcmodel.ge.1) .or. (oceanmodel.eq.2) .or. (ipbl.ge.1) .or. (bbc.eq.3) )then
        ibl=ib
        iel=ie
        jbl=jb
        jel=je
      else
        ibl=1
        iel=1
        jbl=1
        jel=1
      endif

      if((iturb.ge.1).or.(ipbl.ge.1))then
        ibc=ib
        iec=ie
        jbc=jb
        jec=je
        kbc=kb
        kec=ke+1
      else
        ibc=1
        iec=1
        jbc=1
        jec=1
        kbc=1
        kec=1
      endif

      if(iturb.eq.1)then
        ibt=ib
        iet=ie
        jbt=jb
        jet=je
        kbt=kb
        ket=ke+1
      else
        ibt=1
        iet=1
        jbt=1
        jet=1
        kbt=1
        ket=1
      endif

      if(iptra.eq.1)then
        ibp=ib
        iep=ie
        jbp=jb
        jep=je
        kbp=kb
        kep=ke
      else
        ibp=1
        iep=1
        jbp=1
        jep=1
        kbp=1
        kep=1
      endif

      if(psolver.eq.4.or.psolver.eq.5.or.ibalance.eq.2)then

        imirror = 0
        jmirror = 0

        ipb=1
        ipe=ni

        jpb=1
        jpe=nj

        if( (wbc.eq.2.or.wbc.eq.3).or.(ebc.eq.2.or.ebc.eq.3) )then

          imirror = 1
          ipe = ni*2

        endif

        if( (sbc.eq.2.or.sbc.eq.3).or.(nbc.eq.2.or.nbc.eq.3) )then

          jmirror = 1
          jpe = nj*2

        endif

        kpb=0
        kpe=nk+1

      else

        ipb=1
        ipe=1
        jpb=1
        jpe=1
        kpb=1
        kpe=1

      endif

      ! for diagnostics:
      ibd = 1
      ied = 1
      jbd = 1
      jed = 1
      kbd = 1
      ked = 1
      ntdiag = 0
      nqdiag = 0
      td_diss = 0
      td_mptend = 0
      qd_vtc = 0
      qd_vtr = 0
      qd_vts = 0
      qd_vtg = 0
      qd_vti = 0

      doit = .false.
      IF( output_dissheat.eq.1 .or. ipbl.eq.2 )THEN
        doit = .true.
        ntdiag = ntdiag + 1
        td_diss = ntdiag
      ENDIF
      IF( output_mptend.eq.1 )THEN
        doit = .true.
        ntdiag = ntdiag + 1
        td_mptend = ntdiag
      ENDIF

      IF( output_fallvel.eq.1 )THEN
        if( ptype.eq.5 )then
          doit = .true.

          nqdiag = nqdiag+1
          qd_vtc = nqdiag

          nqdiag = nqdiag+1
          qd_vtr = nqdiag

          nqdiag = nqdiag+1
          qd_vts = nqdiag

          nqdiag = nqdiag+1
          qd_vtg = nqdiag

          nqdiag = nqdiag+1
          qd_vti = nqdiag

        endif
      ENDIF

      ntdiag = max( ntdiag , 1 )
      nqdiag = max( nqdiag , 1 )

      IF( doit )THEN
        ibd = ib
        ied = ie
        jbd = jb
        jed = je
        kbd = kb
        ked = ke
      ENDIF






      mynode = 0
      nodemaster = 0
      nodes = 1
      ppnode = 1


      d2i  = 1
      d2is = 1
      d2iu = 1
      d2iv = 1

      d2j  = 1
      d2js = 1
      d2ju = 1
      d2jv = 1

      d3i  = 1
      d3is = 1
      d3iu = 1
      d3iv = 1

      d3j  = 1
      d3js = 1
      d3ju = 1
      d3jv = 1

      d3n = 1
      d3t = 1

      IF( myid.eq.nodemaster )THEN

        d2i   = (ni+1)*nodex
        d2is  = nx
        d2iu  = nx+1
        d2iv  = nx

        d2j   = (nj+1)*nodey
        d2js  = ny
        d2ju  = ny
        d2jv  = ny+1

        d3i   = ni+1
        d3is  = ni
        d3iu  = ni+1
        d3iv  = ni

        d3j   = nj+1
        d3js  = nj
        d3ju  = nj
        d3jv  = nj+1

        d3n = max( numprocs , 1 )

        d3t = max( ppnode-1 + nodes-1 , 1 )

      ENDIF

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  myid,mynode             = ',myid,mynode
      if(dowr) write(outfile,*) '  nodemaster,nodes,ppnode = ',nodemaster,nodes,ppnode
      if(dowr) write(outfile,*) '  nodex,nodey             = ',nodex,nodey
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  d2i   = ',d2i
      if(dowr) write(outfile,*) '  d2is  = ',d2is
      if(dowr) write(outfile,*) '  d2iu  = ',d2iu
      if(dowr) write(outfile,*) '  d2iv  = ',d2iv
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  d2j   = ',d2j
      if(dowr) write(outfile,*) '  d2js  = ',d2js
      if(dowr) write(outfile,*) '  d2ju  = ',d2ju
      if(dowr) write(outfile,*) '  d2jv  = ',d2jv
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  d3i   = ',d3i
      if(dowr) write(outfile,*) '  d3is  = ',d3is
      if(dowr) write(outfile,*) '  d3iu  = ',d3iu
      if(dowr) write(outfile,*) '  d3iv  = ',d3iv
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  d3j   = ',d3j
      if(dowr) write(outfile,*) '  d3js  = ',d3js
      if(dowr) write(outfile,*) '  d3ju  = ',d3ju
      if(dowr) write(outfile,*) '  d3jv  = ',d3jv
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  d3n   = ',d3n
      if(dowr) write(outfile,*) '  d3t   = ',d3t
      if(dowr) write(outfile,*)

!--------------------------------------------------------------

      stat_w       = max(0,min(1,stat_w))
      stat_u       = max(0,min(1,stat_u))
      stat_v       = max(0,min(1,stat_v))
      stat_rmw     = max(0,min(1,stat_rmw))
      IF(axisymm.ne.1) stat_rmw = 0
      stat_pipert  = max(0,min(1,stat_pipert))
      stat_prspert = max(0,min(1,stat_prspert))
      stat_thpert  = max(0,min(1,stat_thpert))
      stat_q       = max(0,min(1,stat_q))
      stat_tke     = max(0,min(1,stat_tke))
      stat_km      = max(0,min(1,stat_km))
      stat_kh      = max(0,min(1,stat_kh))
      stat_div     = max(0,min(1,stat_div))
      stat_rh      = max(0,min(1,stat_rh))
      stat_rhi     = max(0,min(1,stat_rhi))
      stat_the     = max(0,min(1,stat_the))
      stat_cloud   = max(0,min(1,stat_cloud))
      stat_sfcprs  = max(0,min(1,stat_sfcprs))
      stat_wsp     = max(0,min(1,stat_wsp))
      stat_cfl     = max(0,min(1,stat_cfl))
      stat_vort    = max(0,min(1,stat_vort))
      stat_tmass   = max(0,min(1,stat_tmass))
      stat_tmois   = max(0,min(1,stat_tmois))
      stat_qmass   = max(0,min(1,stat_qmass))
      stat_tenerg  = max(0,min(1,stat_tenerg))
      stat_mo      = max(0,min(1,stat_mo))
      stat_tmf     = max(0,min(1,stat_tmf))
      stat_pcn     = max(0,min(1,stat_pcn))
      stat_qsrc    = max(0,min(1,stat_qsrc))


      if(imoist.eq.0)then
        stat_q=0
        stat_rh=0
        stat_rhi=0
        stat_the=0
        stat_cloud=0
        stat_tmois=0
        stat_qmass=0
        stat_pcn=0
        stat_qsrc=0
      endif
      if(iice.eq.0)then
        stat_rhi=0
      endif
      if(iturb.eq.0.or.dns.eq.1)then
        stat_tke=0
        stat_km=0
        stat_kh=0
      endif 
      if(iturb.eq.2.or.iturb.eq.3)then
        stat_tke=0
      endif


      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'stat_w       = ',stat_w
      if(dowr) write(outfile,*) 'stat_u       = ',stat_u
      if(dowr) write(outfile,*) 'stat_v       = ',stat_v
      if(dowr) write(outfile,*) 'stat_rmw     = ',stat_rmw
      if(dowr) write(outfile,*) 'stat_pipert  = ',stat_pipert
      if(dowr) write(outfile,*) 'stat_prspert = ',stat_prspert
      if(dowr) write(outfile,*) 'stat_thpert  = ',stat_thpert
      if(dowr) write(outfile,*) 'stat_q       = ',stat_q
      if(dowr) write(outfile,*) 'stat_tke     = ',stat_tke
      if(dowr) write(outfile,*) 'stat_km      = ',stat_km
      if(dowr) write(outfile,*) 'stat_kh      = ',stat_kh
      if(dowr) write(outfile,*) 'stat_div     = ',stat_div
      if(dowr) write(outfile,*) 'stat_rh      = ',stat_rh
      if(dowr) write(outfile,*) 'stat_rhi     = ',stat_rhi
      if(dowr) write(outfile,*) 'stat_the     = ',stat_the
      if(dowr) write(outfile,*) 'stat_cloud   = ',stat_cloud
      if(dowr) write(outfile,*) 'stat_sfcprs  = ',stat_sfcprs
      if(dowr) write(outfile,*) 'stat_wsp     = ',stat_wsp
      if(dowr) write(outfile,*) 'stat_cfl     = ',stat_cfl
      if(dowr) write(outfile,*) 'stat_vort    = ',stat_vort
      if(dowr) write(outfile,*) 'stat_tmass   = ',stat_tmass
      if(dowr) write(outfile,*) 'stat_tmois   = ',stat_tmois
      if(dowr) write(outfile,*) 'stat_qmass   = ',stat_qmass
      if(dowr) write(outfile,*) 'stat_tenerg  = ',stat_tenerg
      if(dowr) write(outfile,*) 'stat_mo      = ',stat_mo
      if(dowr) write(outfile,*) 'stat_tmf     = ',stat_tmf
      if(dowr) write(outfile,*) 'stat_pcn     = ',stat_pcn
      if(dowr) write(outfile,*) 'stat_qsrc    = ',stat_qsrc
      if(dowr) write(outfile,*)

      stat_out=2*(stat_w+stat_pipert+stat_prspert+numq*stat_q+              &
              stat_tke+2*stat_km+2*stat_kh+stat_div+stat_rh+stat_rhi+       &
              stat_cloud+2*stat_sfcprs+2*stat_wsp)  +                       &
              4*(stat_thpert+stat_u+stat_v)  + 2*stat_rmw +                 &
              3*stat_cfl  +  6*stat_vort  +  stat_tmass  +  stat_tmois  +   &
              (1+(1+nql2-nql1)+iice*(1+nqs2-nqs1))*stat_qmass +             &
              5*stat_tenerg  +  3*stat_mo  +                                &
              nbudget*stat_pcn  + numq*2*stat_qsrc +                        &
              4*stat_the  +  2*stat_tmf + 2*iptra*npt
      IF( adapt_dt.eq.1 ) stat_out = stat_out + 1
      IF( stat_wsp.eq.1 .and. bbc.eq.3 ) stat_out = stat_out + 2
      IF( stat_cfl.ge.1 .and. iturb.lt.1 ) stat_out = stat_out - 2
!!!      stat_out = max(1,stat_out)

      if(dowr) write(outfile,*) 'stat_out = ',stat_out
      if(dowr) write(outfile,*)


      prcl_th       = max(0,min(1,prcl_th))
      prcl_t        = max(0,min(1,prcl_t))
      prcl_prs      = max(0,min(1,prcl_prs))
      prcl_ptra     = max(0,min(1,prcl_ptra))
      prcl_q        = max(0,min(1,prcl_q))
      prcl_nc       = max(0,min(1,prcl_nc))
      prcl_km       = max(0,min(1,prcl_km))
      prcl_kh       = max(0,min(1,prcl_kh))
      prcl_tke      = max(0,min(1,prcl_tke))
      prcl_dbz      = max(0,min(1,prcl_dbz))
      prcl_b        = max(0,min(1,prcl_b))
      prcl_vpg      = max(0,min(1,prcl_vpg))
      prcl_vort     = max(0,min(1,prcl_vort))
      prcl_rho      = max(0,min(1,prcl_rho))
      prcl_qsat     = max(0,min(1,prcl_qsat))
      prcl_sfc      = max(0,min(1,prcl_sfc))


      if( iptra.eq.0 )then
        prcl_ptra = 0
      endif
      if( imoist.eq.0 )then
        prcl_q  = 0
        prcl_nc = 0
        prcl_dbz = 0
        prcl_qsat = 0
      else
        if( ptype.eq.1 .or. ptype.eq.4 )then
          prcl_dbz = 0
        endif
        if( idm.eq.0 )then
          prcl_nc = 0
        endif
      endif
      if( iturb.eq.0 )then
        prcl_km = 0
        prcl_kh = 0
      endif
      if( iturb.ne.1 )then
        prcl_tke = 0
      endif
      if( bbc.ne.3 )then
        prcl_sfc = 0
      endif


      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'prcl_th      = ',prcl_th
      if(dowr) write(outfile,*) 'prcl_t       = ',prcl_t
      if(dowr) write(outfile,*) 'prcl_prs     = ',prcl_prs
      if(dowr) write(outfile,*) 'prcl_ptra    = ',prcl_ptra
      if(dowr) write(outfile,*) 'prcl_q       = ',prcl_q
      if(dowr) write(outfile,*) 'prcl_nc      = ',prcl_nc
      if(dowr) write(outfile,*) 'prcl_km      = ',prcl_km
      if(dowr) write(outfile,*) 'prcl_kh      = ',prcl_kh
      if(dowr) write(outfile,*) 'prcl_tke     = ',prcl_tke
      if(dowr) write(outfile,*) 'prcl_dbz     = ',prcl_dbz
      if(dowr) write(outfile,*) 'prcl_b       = ',prcl_b
      if(dowr) write(outfile,*) 'prcl_vpg     = ',prcl_vpg
      if(dowr) write(outfile,*) 'prcl_vort    = ',prcl_vort
      if(dowr) write(outfile,*) 'prcl_rho     = ',prcl_rho
      if(dowr) write(outfile,*) 'prcl_qsat    = ',prcl_qsat
      if(dowr) write(outfile,*) 'prcl_sfc     = ',prcl_sfc
      if(dowr) write(outfile,*)

!------------------------------------------------------------------
!  constants in subgrid turbulence schemes:

      c_e1 = c_m * c_l * c_l * ( 1.0 / ri_c - 1.0 )
      c_e2 = max( 0.0 , c_m * pi * pi - c_e1 )

      c_s = ( c_m * c_m * c_m / ( c_e1 + c_e2 ) )**0.25
      rcs = 1.0/c_s

      IF( dowr .and. ( iturb.eq.1 .or. iturb.eq.2 ) )THEN
        write(outfile,*)
        write(outfile,*) '  c_m,c_l,ri_c  = ',c_m,c_l,ri_c
        write(outfile,*) '  c_e1,c_e2,sum = ',c_e1,c_e2,c_e1+c_e2
        write(outfile,*) '  c_s,rcs       = ',c_s,rcs
      ENDIF

!--------------------------------------------------------------

      rdx=1.0/dx
      rdy=1.0/dy
      rdz=1.0/dz
      rdx2=1.0/(2.0*dx)
      rdy2=1.0/(2.0*dy)
      rdz2=1.0/(2.0*dz)
      rdx4=1.0/(4.0*dx)
      rdy4=1.0/(4.0*dy)
      rdz4=1.0/(4.0*dz)

      thec_mb=0.0
      qt_mb=0.0

      stattim=statfrq
      taptim=tapfrq
      rsttim=rstfrq
      ! GHB:  now setting prcltim in cm1.F
!!!      prcltim=prclfrq
      if( iprcl.ne.1 ) prcltim = 1.0d30
      radtim=0.0

!--------------------------------------------------------------

      npvals = 1

      prx = 0
      pry = 0
      prz = 0
      pru = 0
      prv = 0
      prw = 0
      prth = 0
      prt = 0
      prprs = 0
      prpt1 = 0
      prpt2 = 0
      prqv = 0
      prq1 = 0
      prq2 = 0
      prnc1 = 0
      prnc2 = 0
      prkm = 0
      prkh = 0
      prtke = 0
      prdbz = 0
      prb = 0
      prvpg = 0
      przv = 0
      prrho = 0
      prqsl = 0
      prqsi = 0
      prznt = 0
      prust = 0

      ! for parcels:
      if(iprcl.eq.1)then

        ! 6 basic variables for all simulations:
        ! (x,y,z,u,v,w)
        !  1 2 3 4 5 6
        npvals = 6

        prx = 1
        pry = 2
        prz = 3
        pru = 4
        prv = 5
        prw = 6

        if( prcl_th.eq.1 )then
          npvals = npvals+1
          prth = npvals
        endif
        if( prcl_t.eq.1 )then
          npvals = npvals+1
          prt = npvals
        endif
        if( prcl_prs.eq.1 )then
          npvals = npvals+1
          prprs = npvals
        endif

        ! passive tracers:
        if( prcl_ptra.eq.1 )then
          prpt1 = npvals+1
          npvals = npvals + npt
          prpt2 = npvals
        endif

        ! moisture variables:
        if( prcl_q.eq.1 )then
          !---
          npvals = npvals+1
          prqv = npvals
          prq1 = npvals+1
          !---
          npvals = npvals+(nql2-nql1+1)
          if( iice.eq.1 )then
            npvals = npvals+(nqs2-nqs1+1)
          endif
          !---
          prq2 = npvals
        endif
        if( prcl_nc.eq.1 .and. idm.eq.1 )then
          prnc1 = npvals+1
          npvals = npvals+(nnc2-nnc1+1)
          prnc2 = npvals
        endif


        ! turbulence parameters:
        if( iturb.ge.1 )then
          if( prcl_km.eq.1 )then
            prkm = npvals+1
            npvals = npvals+2
          endif
          if( prcl_kh.eq.1 )then
            prkh = npvals+1
            npvals = npvals+2
          endif
          if( prcl_tke.eq.1 )then
            npvals = npvals+1
            prtke = npvals
          endif
        endif

        if( prcl_dbz.eq.1 )then
          npvals = npvals+1
          prdbz = npvals
        endif
        if( prcl_b.eq.1 )then
          npvals = npvals+1
          prb = npvals
        endif
        if( prcl_vpg.eq.1 )then
          npvals = npvals+1
          prvpg = npvals
        endif
        if( prcl_vort.eq.1 )then
          npvals = npvals+1
          przv = npvals
        endif
        if( prcl_rho.eq.1 )then
          npvals = npvals+1
          prrho = npvals
        endif
        if( prcl_qsat.eq.1 )then
          npvals = npvals+1
          prqsl = npvals
          if( iice.eq.1 )then
            npvals = npvals+1
            prqsi = npvals
          endif
        endif
        if( prcl_sfc.eq.1 )then
          npvals = npvals+1
          prznt = npvals
          npvals = npvals+1
          prust = npvals
        endif

      else

        nparcels = 1

      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  nparcels = ',nparcels
      if(dowr) write(outfile,*) '  npvals   = ',npvals
      if(dowr) write(outfile,*) '    prx      = ',prx
      if(dowr) write(outfile,*) '    pry      = ',pry
      if(dowr) write(outfile,*) '    prz      = ',prz
      if(dowr) write(outfile,*) '    pru      = ',pru
      if(dowr) write(outfile,*) '    prv      = ',prv
      if(dowr) write(outfile,*) '    prw      = ',prw
      if(dowr) write(outfile,*) '    prth     = ',prth
      if(dowr) write(outfile,*) '    prt      = ',prt
      if(dowr) write(outfile,*) '    prprs    = ',prprs
      if(dowr) write(outfile,*) '    prpt1    = ',prpt1
      if(dowr) write(outfile,*) '    prpt2    = ',prpt2
      if(dowr) write(outfile,*) '    prqv     = ',prqv
      if(dowr) write(outfile,*) '    prq1     = ',prq1
      if(dowr) write(outfile,*) '    prq2     = ',prq2
      if(dowr) write(outfile,*) '    prnc1    = ',prnc1
      if(dowr) write(outfile,*) '    prnc2    = ',prnc2
      if(dowr) write(outfile,*) '    prkm     = ',prkm
      if(dowr) write(outfile,*) '    prkh     = ',prkh
      if(dowr) write(outfile,*) '    prtke    = ',prtke
      if(dowr) write(outfile,*) '    prdbz    = ',prdbz
      if(dowr) write(outfile,*) '    prb      = ',prb
      if(dowr) write(outfile,*) '    prvpg    = ',prvpg
      if(dowr) write(outfile,*) '    przv     = ',przv
      if(dowr) write(outfile,*) '    prrho    = ',prrho
      if(dowr) write(outfile,*) '    prqsl    = ',prqsl
      if(dowr) write(outfile,*) '    prqsi    = ',prqsi
      if(dowr) write(outfile,*) '    prznt    = ',prznt
      if(dowr) write(outfile,*) '    prust    = ',prust
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  Get identity

      ibw=0
      ibe=0
      ibs=0
      ibn=0

      patchsws = .false.
      patchsww = .false.
      patchses = .false.
      patchsee = .false.
      patchnwn = .false.
      patchnww = .false.
      patchnen = .false.
      patchnee = .false.

      p2tchsws = .false.
      p2tchsww = .false.
      p2tchses = .false.
      p2tchsee = .false.
      p2tchnwn = .false.
      p2tchnww = .false.
      p2tchnen = .false.
      p2tchnee = .false.

      myi=1
      myj=1



    allocate( xfdp(-2:nx+4) )
    allocate( yfdp(-2:ny+4) )

    IF(iorigin.eq.1)THEN

      do i=-2,nx+4
        xfdp(i)=dble(dx)*(i-1)
      enddo
      do i=ib,ie+1
        xf(i)=xfdp(i+(myi-1)*ni)
      enddo
      do i=ib,ie
        xh(i)=0.5d0*(xfdp(i+1+(myi-1)*ni)+xfdp(i+(myi-1)*ni))
      enddo

        i = ib
        if(i+(myi-1)*nx/nodex.eq.(-2) .and. wbc.ne.1) ibw=1
        i = ie
        if(i+(myi-1)*nx/nodex.eq.(nx+3) .and. ebc.ne.1) ibe=1

      do j=-2,ny+4
        yfdp(j)=dble(dy)*(j-1)
      enddo
      do j=jb,je+1
        yf(j)=yfdp(j+(myj-1)*nj)
      enddo
      do j=jb,je
        yh(j)=0.5d0*(yfdp(j+1+(myj-1)*nj)+yfdp(j+(myj-1)*nj))
      enddo

        j = jb
        if(j+(myj-1)*ny/nodey.eq.(-2) .and. sbc.ne.1) ibs=1
        j = je
        if(j+(myj-1)*ny/nodey.eq.(ny+3) .and. nbc.ne.1) ibn=1

    ELSEIF(iorigin.eq.2)THEN

      do i=-2,nx+4
        xfdp(i)=dble(dx)*(i-1)-0.5d0*dble(dx)*nx
      enddo
      do i=ib,ie+1
        xf(i)=xfdp(i+(myi-1)*ni)
      enddo
      do i=ib,ie
        xh(i)=0.5d0*(xfdp(i+1+(myi-1)*ni)+xfdp(i+(myi-1)*ni))
      enddo

        i = ib
        if(i+(myi-1)*nx/nodex.eq.(-2) .and. wbc.ne.1) ibw=1
        i = ie
        if(i+(myi-1)*nx/nodex.eq.(nx+3) .and. ebc.ne.1) ibe=1

      do j=-2,ny+4
        yfdp(j)=dble(dy)*(j-1)-0.5d0*dble(dy)*ny
      enddo
      do j=jb,je+1
        yf(j)=yfdp(j+(myj-1)*nj)
      enddo
      do j=jb,je
        yh(j)=0.5d0*(yfdp(j+1+(myj-1)*nj)+yfdp(j+(myj-1)*nj))
      enddo

        j = jb
        if(j+(myj-1)*ny/nodey.eq.(-2) .and. sbc.ne.1) ibs=1
        j = je
        if(j+(myj-1)*ny/nodey.eq.(ny+3) .and. nbc.ne.1) ibn=1

    ELSE

      print *,'  invalid option for iorigin'



      call stopcm1

    ENDIF


      if(wbc.eq.2)then
        ibw=1
      endif

      if(ebc.eq.2)then
        ibe=1
      endif

      if(sbc.eq.2)then
        ibs=1
      endif

      if(nbc.eq.2)then
        ibn=1
      endif


!--------------------------------------------------------------

      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) 'g     =',g
      if(dowr) write(outfile,*) 'to    =',to
      if(dowr) write(outfile,*) 'rd    =',rd
      if(dowr) write(outfile,*) 'rv    =',rv
      if(dowr) write(outfile,*) 'cp    =',cp
      if(dowr) write(outfile,*) 'cv    =',cv
      if(dowr) write(outfile,*) 'cpv   =',cpv
      if(dowr) write(outfile,*) 'cvv   =',cvv
      if(dowr) write(outfile,*) 'p00   =',p00
      if(dowr) write(outfile,*) 'rp00  =',rp00
      if(dowr) write(outfile,*) 'th0r  =',th0r
      if(dowr) write(outfile,*) 'rcp   =',rcp
      if(dowr) write(outfile,*) 'pi    =',pi

      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) 'cpdcv =',cpdcv
      if(dowr) write(outfile,*) 'rovcp =',rovcp
      if(dowr) write(outfile,*) 'rddcv =',rddcv
      if(dowr) write(outfile,*) 'cvdrd =',cvdrd
      if(dowr) write(outfile,*) 'cpdrd =',cpdrd
      if(dowr) write(outfile,*) 'eps   =',eps
      if(dowr) write(outfile,*) 'reps  =',reps
      if(dowr) write(outfile,*) 'repsm1=',repsm1
      if(dowr) write(outfile,*) 'cpt   =',cpt
      if(dowr) write(outfile,*) 'cvt   =',cvt
      if(dowr) write(outfile,*) 'pnum  =',pnum
      if(dowr) write(outfile,*) 'xlv   =',xlv
      if(dowr) write(outfile,*) 'xls   =',xls
      if(dowr) write(outfile,*) 'lvdcp =',lvdcp
      if(dowr) write(outfile,*) 'condc =',condc
      if(dowr) write(outfile,*) 'cpl   =',cpl
      if(dowr) write(outfile,*) 'cpi   =',cpi
      if(dowr) write(outfile,*) 'lv1   =',lv1
      if(dowr) write(outfile,*) 'lv2   =',lv2
      if(dowr) write(outfile,*) 'ls1   =',ls1
      if(dowr) write(outfile,*) 'ls2   =',ls2
      if(dowr) write(outfile,*) 'karman=',karman

      if( psolver.eq.6 )then
        !----------------
        ! speed of sound for compressible-Boussinesq equations (psolver=6)
        csound = 300.0
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) 'csound=',csound
        if(dowr) write(outfile,*)
      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'timeformat   =',timeformat
      if(dowr) write(outfile,*) 'timestats    =',timestats
      if(dowr) write(outfile,*) 'terrain_flag =',terrain_flag
      if(dowr) write(outfile,*) 'procfiles    =',procfiles


      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'nx    =',nx
      if(dowr) write(outfile,*) 'ny    =',ny
      if(dowr) write(outfile,*) 'nz    =',nz










      if(dowr) write(outfile,*)
 
      if(dowr) write(outfile,*) 'ni    =',ni
      if(dowr) write(outfile,*) 'nj    =',nj
      if(dowr) write(outfile,*) 'nk    =',nk
      if(dowr) write(outfile,*) 'nkp1  =',nkp1

      if(dowr) write(outfile,*)
 
      if(dowr) write(outfile,130) 'ib,ibm,ibi,ibc,ibt=',ib,ibm,ibi,ibc,ibt
      if(dowr) write(outfile,130) 'ie,iem,iei,iec,iet=',ie,iem,iei,iec,iet
      if(dowr) write(outfile,130) 'jb,jbm,jbi,jbc,jbt=',jb,jbm,jbi,jbc,jbt
      if(dowr) write(outfile,130) 'je,jem,jei,jec,jet=',je,jem,jei,jec,jet
      if(dowr) write(outfile,130) 'kb,kbm,kbi,kbc,kbt=',kb,kbm,kbi,kbc,kbt
      if(dowr) write(outfile,130) 'ke,kem,kei,kec,ket=',ke,kem,kei,kec,ket

130   format(1x,a19,5(4x,i5))

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'imirror,jmirror = ',imirror,jmirror
      if(dowr) write(outfile,*)

      if(dowr) write(outfile,131) 'ibp,itb,ipb,ibr,ibb=',ibp,itb,ipb,ibr,ibb
      if(dowr) write(outfile,131) 'iep,ite,ipe,ier,ieb=',iep,ite,ipe,ier,ieb
      if(dowr) write(outfile,131) 'jbp,jtb,jpb,jbr,jbb=',jbp,jtb,jpb,jbr,jbb
      if(dowr) write(outfile,131) 'jep,jte,jpe,jer,jeb=',jep,jte,jpe,jer,jeb
      if(dowr) write(outfile,131) 'kbp,ktb,kpb,kbr,kbb=',kbp,ktb,kpb,kbr,kbb
      if(dowr) write(outfile,131) 'kep,kte,kpe,ker,keb=',kep,kte,kpe,ker,keb

131   format(1x,a20,5(4x,i5))

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,132) 'ibd               =',ibd
      if(dowr) write(outfile,132) 'ied               =',ied
      if(dowr) write(outfile,132) 'jbd               =',jbd
      if(dowr) write(outfile,132) 'jed               =',jed
      if(dowr) write(outfile,132) 'kbd               =',kbd
      if(dowr) write(outfile,132) 'ked               =',ked

      if( dowr )then
        if( td_diss   .gt.0 ) write(outfile,*),'  td_diss   = ',td_diss
        if( td_mptend .gt.0 ) write(outfile,*),'  td_mptend = ',td_mptend
        if( qd_vtc    .gt.0 ) write(outfile,*),'  qd_vtc    = ',qd_vtc
        if( qd_vtr    .gt.0 ) write(outfile,*),'  qd_vtr    = ',qd_vtr
        if( qd_vts    .gt.0 ) write(outfile,*),'  qd_vts    = ',qd_vts
        if( qd_vtg    .gt.0 ) write(outfile,*),'  qd_vtg    = ',qd_vtg
        if( qd_vti    .gt.0 ) write(outfile,*),'  qd_vti    = ',qd_vti
      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,132) 'ibl               =',ibl
      if(dowr) write(outfile,132) 'iel               =',iel
      if(dowr) write(outfile,132) 'jbl               =',jbl
      if(dowr) write(outfile,132) 'jel               =',jel

132   format(1x,a19,1(4x,i5))

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,133) imp,jmp,kmp,kmt,rmp,cmp

133   format(' imp,jmp,kmp,kmt,rmp,cmp =',6(2x,i5))

!----------

      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) 'rdx    =',rdx
      if(dowr) write(outfile,*) 'rdy    =',rdy
      if(dowr) write(outfile,*) 'rdz    =',rdz
      if(dowr) write(outfile,*) 'rdx2   =',rdx2
      if(dowr) write(outfile,*) 'rdy2   =',rdy2
      if(dowr) write(outfile,*) 'rdz2   =',rdz2
      if(dowr) write(outfile,*) 'rdx4   =',rdx4
      if(dowr) write(outfile,*) 'rdy4   =',rdy4
      if(dowr) write(outfile,*) 'rdz4   =',rdz4
      if(dowr) write(outfile,*) 'govtwo =',govtwo
      if(dowr) write(outfile,*) 'clwsat =',clwsat
      if(dowr) write(outfile,*) 'smeps  =',smeps
      if(dowr) write(outfile,*) 'tsmall =',tsmall
      if(dowr) write(outfile,*) 'qsmall =',qsmall
      if(dowr) write(outfile,*) 'cstar  =',cstar
      if(dowr) write(outfile,*) 'csmax  =',csmax
      if(dowr) write(outfile,*) 'epsilon=',epsilon

      if(dowr) write(outfile,*)

!--------------------------------------------------------------

      do i=ib,ie
        uh(i)=1.0
      enddo

      do i=ib,ie+1
        uf(i)=1.0
      enddo

      strx:  IF(stretch_x.ge.1)THEN

!!!        ibw=0
!!!        ibe=0

        ni1 = 0
        ni2 = 0
        ni3 = 0

!-----------------------------------------------------------------------
!  Begin hard-wired analytic stretching function

        nominal_dx = 0.5*( dx_inner + dx_outer )

      IF(stretch_x.eq.1)THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_x = 1 ... stretching on both west and east sides of domain:'
        if(dowr) write(outfile,*)
        ni1 = nint( (tot_x_len-nos_x_len)*0.5/nominal_dx )
        ni2 = nint( nos_x_len/dx_inner )
        ni3 = ni1
        if(dowr) write(outfile,*) '  ni1,ni2,ni3 = ',(tot_x_len-nos_x_len)*0.5/nominal_dx,   &
                         nos_x_len/dx_inner,(tot_x_len-nos_x_len)*0.5/nominal_dx
        if(dowr) write(outfile,*) '    (note:  ni1,ni2,ni3 need to be exact integers for this to work correctly)'
      ELSEIF(stretch_x.eq.2)THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_x = 2 ... stretching on east side of domain only:'
        if(dowr) write(outfile,*)
        ni1 = 0
        ni2 = nint( nos_x_len/dx_inner )
        ni3 = nint( (tot_x_len-nos_x_len)/nominal_dx )
        if(dowr) write(outfile,*) '  ni1,ni2,ni3 = ',0.0,nos_x_len/dx_inner,(tot_x_len-nos_x_len)/nominal_dx
        if(dowr) write(outfile,*) '    (note:  ni1,ni2,ni3 need to be exact integers for this to work correctly)'
      ELSE
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_x must be either 1 or 2'
        if(dowr) write(outfile,*)
        call stopcm1
      ENDIF

        c2=(nominal_dx-dx_inner)/(nominal_dx*nominal_dx*float(ni3-1))
        c1=(dx_inner/nominal_dx)-c2*nominal_dx

        if(dowr) write(outfile,*) '  nominal_dx  = ',nominal_dx
        if(dowr) write(outfile,*) '  c1,c2       = ',c1,c2
        if(dowr) write(outfile,*)

        ! Test to see if nx is kosher.
      IF(stretch_x.eq.1)THEN
        if( nx.ne.(ni1+ni2+ni3) .or. ni1.lt.0 .or. ni2.lt.0 .or. ni3.lt.0 )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  There is a problem with the settings for horizontal grid stretching'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of nx = ',nx
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '       ni1,ni2,ni3 = ',ni1,ni2,ni3
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value needed for these settings ...'
          if(dowr) write(outfile,*) '       dx_inner  = ',dx_inner
          if(dowr) write(outfile,*) '       dx_outer  = ',dx_outer
          if(dowr) write(outfile,*) '       nos_x_len = ',nos_x_len
          if(dowr) write(outfile,*) '       tot_x_len = ',tot_x_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... would be nx = ',(nos_x_len/dx_inner)+(tot_x_len-nos_x_len)/(0.5*(dx_inner+dx_outer))
          if(dowr) write(outfile,*) '  (if this number is an integer) '
          if(dowr) write(outfile,*) '  (and if ni1,ni2,ni3 are all integers) '
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
          call stopcm1
        endif
      ELSEIF(stretch_x.eq.2)THEN
        if( nx.ne.(ni1+ni2+ni3) .or. ni1.lt.0 .or. ni2.lt.0 .or. ni3.lt.0 )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  There is a problem with the settings for horizontal grid stretching'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of nx = ',nx
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '       ni1,ni2,ni3 = ',ni1,ni2,ni3
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value for these settings ...'
          if(dowr) write(outfile,*) '       dx_inner  = ',dx_inner
          if(dowr) write(outfile,*) '       dx_outer  = ',dx_outer
          if(dowr) write(outfile,*) '       nos_x_len = ',nos_x_len
          if(dowr) write(outfile,*) '       tot_x_len = ',tot_x_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... would be nx = ',(nos_x_len/dx_inner)+(tot_x_len-nos_x_len)/(0.5*(dx_inner+dx_outer))
          if(dowr) write(outfile,*) '  (if this number is an integer) '
          if(dowr) write(outfile,*) '  (and if ni1,ni2,ni3 are all integers) '
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
          call stopcm1
        endif
      ENDIF

        mult = 0.0
        if(iorigin.eq.2) mult = 0.5

      IF(stretch_x.eq.1)THEN

        do i=ni1+1,ni1+ni2+1
            xfdp(i)=ni1*nominal_dx+(i-ni1-1)*dx_inner - mult*tot_x_len
        enddo
        do i=ni1+ni2+2,ni1+ni2+ni3+4
            xfdp(i)=ni1*nominal_dx+(ni1+ni2+1-ni1-1)*dble(dx_inner)   &
                 +(c1+c2*dble(i-1-ni1-ni2)*nominal_dx)   &
                 *dble(i-1-ni1-ni2)*nominal_dx - mult*tot_x_len
        enddo
        do i=-2,ni1
            xfdp(i)=ni1*nominal_dx+(ni1+1-ni1-1)*dble(dx_inner)    &
                 -(c1+c2*dble(ni1+1-i)*nominal_dx)   &
                 *dble(ni1+1-i)*nominal_dx - mult*tot_x_len
        enddo

      ELSEIF(stretch_x.eq.2)THEN

        do i=ni1+1,ni1+ni2+1
            xfdp(i)=ni1*nominal_dx+(i-ni1-1)*dx_inner - mult*tot_x_len
        enddo
        do i=ni1+ni2+2,ni1+ni2+ni3+3
            xfdp(i)=ni1*nominal_dx+(ni1+ni2+1-ni1-1)*dble(dx_inner)   &
                 +(c1+c2*dble(i-1-ni1-ni2)*nominal_dx)   &
                 *dble(i-1-ni1-ni2)*nominal_dx - mult*tot_x_len
        enddo
        do i=-2,ni1
            xfdp(i)=ni1*nominal_dx+(ni1+1-ni1-1)*dble(dx_inner)    &
                 -(c1+c2*dble(ni1+1-i)*nominal_dx)   &
                 *dble(ni1+1-i)*nominal_dx - mult*tot_x_len
        enddo

      ENDIF

!!!        if( xf(ib).lt.0.0  .and. wbc.ne.1 ) ibw=1
!!!        if( xf(ie).gt.maxx .and. ebc.ne.1 ) ibe=1

        IF(stretch_x.eq.1)THEN
          xfdp( 0)=xfdp(1)-1*dx_outer
          xfdp(-1)=xfdp(1)-2*dx_outer
          xfdp(-2)=xfdp(1)-3*dx_outer
        ELSEIF(stretch_x.eq.2)THEN
          xfdp( 0)=xfdp(1)-1*dx_inner
          xfdp(-1)=xfdp(1)-2*dx_inner
          xfdp(-2)=xfdp(1)-3*dx_inner
        ENDIF

          xfdp(nx+2)=xfdp(nx+1)+1*dx_outer
          xfdp(nx+3)=xfdp(nx+1)+2*dx_outer
          xfdp(nx+4)=xfdp(nx+1)+3*dx_outer

!  End hard-wired analytic stretching function
!-----------------------------------------------------------------------
!
!  Optional:  to use a different stretching function, or to use 
!  arbitrarily located grid points, simply comment out the 
!  "hard-wired" section above, and then specify values for xfref
!  here.  Do not change anything below here!
!
!  Note:  xfref stores the location of the staggered u points for
!  the entire domain (from x=-2 to x=nx+4) (note: this includes
!  the boundary points that extend 3 gridpoints beyond the
!  computational domain.
!
!-----------------------------------------------------------------------

      ENDIF  strx

        do i=ib,ie+1
          xf(i)=xfdp(i+(myi-1)*ni)
        enddo

        do i=ib,ie
          xh(i)=0.5d0*(xfdp(i+1+(myi-1)*ni)+xfdp(i+(myi-1)*ni))
          uh(i)=dble(dx)/(xfdp(i+1+(myi-1)*ni)-xfdp(i+(myi-1)*ni))
        enddo

        arh1 = 1.0
        arh2 = 1.0
        arf1 = 1.0
        arf2 = 1.0

      IF(axisymm.eq.1)THEN

        print *
        do i=ib,ie
          arh1(i) = ( xfdp(i  )/( 0.5d0*(xfdp(i+1)+xfdp(i)) ) )
          arh2(i) = ( xfdp(i+1)/( 0.5d0*(xfdp(i+1)+xfdp(i)) ) )
          print *,'  arh1,arh2 = ',i,arh1(i),arh2(i),0.5*(arh1(i)+arh2(i))
        enddo
        print *
        print *
        do i=ib+1,ie
          if( abs(xfdp(i)).le.smeps )then
            arf1(i) = 1.0
            arf2(i) = 1.0
          else
            arf1(i) = ( 0.5d0*(xfdp(i-1)+xfdp(i)) / xfdp(i) )
            arf2(i) = ( 0.5d0*(xfdp(i+1)+xfdp(i)) / xfdp(i) )
          endif
          print *,'  arf1,arf2 = ',i,arf1(i),arf2(i),0.5*(arf1(i)+arf2(i))
        enddo
        print *

      ENDIF

        do i=ib+1,ie
          uf(i)=dble(dx)/( 0.5d0*(xfdp(i+1+(myi-1)*ni)+xfdp(i+(myi-1)*ni)) &
                          -0.5d0*(xfdp(i-1+(myi-1)*ni)+xfdp(i+(myi-1)*ni)) )
        enddo

        if(ibw.eq.1)then
          uf( 0)=uf(1)
          uf(-1)=uf(1)
          uf(-2)=uf(1)
        endif

        if(ibe.eq.1)then
          uf(ni+2)=uf(ni+1)
          uf(ni+3)=uf(ni+1)
          uf(ni+4)=uf(ni+1)
        endif

      do i=ib,ie
        rxh(i)=1.0/(smeps+xh(i))
        ruh(i)=1.0/uh(i)
      enddo

      do i=ib,ie+1
        rxf(i)=1.0/(smeps+xf(i))
        ruf(i)=1.0/uf(i)
      enddo

      do i=-2,nx+4
        xfref(i) = xfdp(i)
      enddo

      minx = xfref(1)
      maxx = xfref(nx+1)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'x:'
      if(dowr) write(outfile,124)
      if(dowr) write(outfile,125)

      do i=ib,ib+2
        if(dowr) write(outfile,122) i,xf(i),xh(i),xf(i+1)-xf(i),uf(i),uh(i),'   x'
      enddo
      do i=ib+3,ie-3
        if(dowr) write(outfile,122) i,xf(i),xh(i),xf(i+1)-xf(i),uf(i),uh(i),'    '
      enddo
      do i=ie-2,ie
        if(dowr) write(outfile,122) i,xf(i),xh(i),xf(i+1)-xf(i),uf(i),uh(i),'   x'
      enddo
      if(dowr) write(outfile,123) ie+1,xf(ie+1),uf(ie+1)



      if(dowr) write(outfile,*)

122   format(3x,i5,3x,f11.2,3x,f11.2,3x,f9.2,3x,f8.4,3x,f8.4,a4)
123   format(3x,i5,3x,f11.2,29x,f8.4)
124   format('      i         xf           xh         dx         uf         uh')
125   format(' ---------------------------------------------------------------')

!--------------------------------------------------------------

      do j=jb,je
        vh(j)=1.0
      enddo

      do j=jb,je+1
        vf(j)=1.0
      enddo

      IF(stretch_y.ge.1)THEN

!!!        ibs=0
!!!        ibn=0

!-----------------------------------------------------------------------
!  Begin hard-wired analytic stretching function

        nominal_dy = 0.5*( dy_inner + dy_outer )

      IF(stretch_y.eq.1)THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_y = 1 ... stretching on both south and north sides of domain:'
        if(dowr) write(outfile,*)
        nj1 = nint( (tot_y_len-nos_y_len)*0.5/nominal_dy )
        nj2 = nint( nos_y_len/dy_inner )
        nj3 = nj1
        if(dowr) write(outfile,*) '  nj1,nj2,nj3 = ',(tot_y_len-nos_y_len)*0.5/nominal_dy,   &
                         nos_y_len/dy_inner,(tot_y_len-nos_y_len)*0.5/nominal_dy
        if(dowr) write(outfile,*) '    (note:  nj1,nj2,nj3 need to be exact integers for this to work correctly)'
      ELSEIF(stretch_y.eq.2)THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_y = 2 ... stretching on north side of domain only:'
        if(dowr) write(outfile,*)
        nj1 = 0
        nj2 = nint( nos_y_len/dy_inner )
        nj3 = nint( (tot_y_len-nos_y_len)/nominal_dy )
        if(dowr) write(outfile,*) '  nj1,nj2,nj3 = ',0.0,nos_y_len/dy_inner,(tot_y_len-nos_y_len)/nominal_dy
        if(dowr) write(outfile,*) '    (note:  nj1,nj2,nj3 need to be exact integers for this to work correctly)'
      ELSE
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_y must be either 1 or 2'
        if(dowr) write(outfile,*)
        call stopcm1
      ENDIF

        c2=(nominal_dy-dy_inner)/(nominal_dy*nominal_dy*float(nj3-1))
        c1=(dy_inner/nominal_dy)-c2*nominal_dy

        if(dowr) write(outfile,*) '  nominal_dy  = ',nominal_dy
        if(dowr) write(outfile,*) '  c1,c2       = ',c1,c2
        if(dowr) write(outfile,*)

        ! Test to see if ny is kosher.
      IF(stretch_y.eq.1)THEN
        if( ny.ne.(nj1+nj2+nj3) .or. nj1.lt.0 .or. nj2.lt.0 .or. nj3.lt.0 )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  There is a problem with the settings for horizontal grid stretching'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of ny = ',ny
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '       nj1,nj2,nj3 = ',nj1,nj2,nj3
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value needed for these settings ...'
          if(dowr) write(outfile,*) '       dy_inner  = ',dy_inner
          if(dowr) write(outfile,*) '       dy_outer  = ',dy_outer
          if(dowr) write(outfile,*) '       nos_y_len = ',nos_y_len
          if(dowr) write(outfile,*) '       tot_y_len = ',tot_y_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... would be ny = ',(nos_y_len/dy_inner)+(tot_y_len-nos_y_len)/(0.5*(dy_inner+dy_outer))
          if(dowr) write(outfile,*) '  (if this number is an integer) '
          if(dowr) write(outfile,*) '  (and if nj1,nj2,nj3 are all integers) '
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
          call stopcm1
        endif
      ELSEIF(stretch_y.eq.2)THEN
        if( ny.ne.(nj1+nj2+nj3) .or. nj1.lt.0 .or. nj2.lt.0 .or. nj3.lt.0 )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  There is a problem with the settings for horizontal grid stretching'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of ny = ',ny
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '       nj1,nj2,nj3 = ',nj1,nj2,nj3
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value for these settings ...'
          if(dowr) write(outfile,*) '       dy_inner  = ',dy_inner
          if(dowr) write(outfile,*) '       dy_outer  = ',dy_outer
          if(dowr) write(outfile,*) '       nos_y_len = ',nos_y_len
          if(dowr) write(outfile,*) '       tot_y_len = ',tot_y_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... would be ny = ',(nos_y_len/dy_inner)+(tot_y_len-nos_y_len)/(0.5*(dy_inner+dy_outer))
          if(dowr) write(outfile,*) '  (if this number is an integer) '
          if(dowr) write(outfile,*) '  (and if nj1,nj2,nj3 are all integers) '
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
          call stopcm1
        endif
      ENDIF

        mult = 0.0
        if(iorigin.eq.2) mult = 0.5

      IF(stretch_y.eq.1)THEN

        do j=nj1+1,nj1+nj2+1
            yfdp(j)=nj1*nominal_dy+(j-nj1-1)*dy_inner - mult*tot_y_len
        enddo
        do j=nj1+nj2+2,nj1+nj2+nj3+4
            yfdp(j)=nj1*nominal_dy+(nj1+nj2+1-nj1-1)*dble(dy_inner)   &
                 +(c1+c2*dble(j-1-nj1-nj2)*nominal_dy)   &
                 *dble(j-1-nj1-nj2)*nominal_dy - mult*tot_y_len
        enddo
        do j=-2,nj1
            yfdp(j)=nj1*nominal_dy+(nj1+1-nj1-1)*dble(dy_inner)    &
                 -(c1+c2*dble(nj1+1-j)*nominal_dy)   &
                 *dble(nj1+1-j)*nominal_dy - mult*tot_y_len
        enddo

      ELSEIF(stretch_y.eq.2)THEN

        do j=nj1+1,nj1+nj2+1
            yfdp(j)=nj1*nominal_dy+(j-nj1-1)*dy_inner - mult*tot_y_len
        enddo
        do j=nj1+nj2+2,nj1+nj2+nj3+3
            yfdp(j)=nj1*nominal_dy+(nj1+nj2+1-nj1-1)*dble(dy_inner)   &
                 +(c1+c2*dble(j-1-nj1-nj2)*nominal_dy)   &
                 *dble(j-1-nj1-nj2)*nominal_dy - mult*tot_y_len
        enddo
        do j=-2,nj1
            yfdp(j)=nj1*nominal_dy+(nj1+1-nj1-1)*dble(dy_inner)    &
                 -(c1+c2*dble(nj1+1-j)*nominal_dy)   &
                 *dble(nj1+1-j)*nominal_dy - mult*tot_y_len
        enddo

      ENDIF

!!!        if( yf(jb).lt.0.0  .and. sbc.ne.1 ) ibs=1
!!!        if( yf(je).gt.maxy .and. nbc.ne.1 ) ibn=1

        IF(stretch_y.eq.1)THEN
          yfdp( 0)=yfdp(1)-1*dy_outer
          yfdp(-1)=yfdp(1)-2*dy_outer
          yfdp(-2)=yfdp(1)-3*dy_outer
        ELSEIF(stretch_y.eq.2)THEN
          yfdp( 0)=yfdp(1)-1*dy_inner
          yfdp(-1)=yfdp(1)-2*dy_inner
          yfdp(-2)=yfdp(1)-3*dy_inner
        ENDIF

          yfdp(ny+2)=yfdp(ny+1)+1*dy_outer
          yfdp(ny+3)=yfdp(ny+1)+2*dy_outer
          yfdp(ny+4)=yfdp(ny+1)+3*dy_outer

!  End hard-wired analytic stretching function
!-----------------------------------------------------------------------
!
!  Optional:  to use a different stretching function, or to use 
!  arbitrarily located grid points, simply comment out the 
!  "hard-wired" section above, and then specify values for yfref
!  here.  Do not change anything below here!
!
!  Note:  yfref stores the location of the staggered v points for
!  the entire domain (from y=-2 to y=ny+4) (note: this includes
!  the boundary points that extend 3 gridpoints beyond the
!  computational domain.
!
!-----------------------------------------------------------------------

        do j=jb,je+1
          yf(j)=yfdp(j+(myj-1)*nj)
        enddo

        do j=jb,je
          yh(j)=0.5d0*(yfdp(j+1+(myj-1)*nj)+yfdp(j+(myj-1)*nj))
          vh(j)=dble(dy)/(yfdp(j+1+(myj-1)*nj)-yfdp(j+(myj-1)*nj))
        enddo

        do j=jb+1,je
          vf(j)=dble(dy)/( 0.5d0*(yfdp(j+1+(myj-1)*nj)+yfdp(j+(myj-1)*nj)) &
                          -0.5d0*(yfdp(j-1+(myj-1)*nj)+yfdp(j+(myj-1)*nj)) )
        enddo

        if(ibs.eq.1)then
          vf( 0)=vf(1)
          vf(-1)=vf(1)
          vf(-2)=vf(1)
        endif

        if(ibn.eq.1)then
          vf(nj+2)=vf(nj+1)
          vf(nj+3)=vf(nj+1)
          vf(nj+4)=vf(nj+1)
        endif

      ENDIF

      do j=jb,je
        rvh(j)=1.0/vh(j)
      enddo

      do j=jb,je+1
        rvf(j)=1.0/vf(j)
      enddo

      do j=-2,ny+4
        yfref(j) = yfdp(j)
      enddo

      miny = yfref(1)
      maxy = yfref(ny+1)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'y:'
      if(dowr) write(outfile,134)
      if(dowr) write(outfile,125)

      do j=jb,jb+2
        if(dowr) write(outfile,122) j,yf(j),yh(j),yf(j+1)-yf(j),vf(j),vh(j),'   x'
      enddo
      do j=jb+3,je-3
        if(dowr) write(outfile,122) j,yf(j),yh(j),yf(j+1)-yf(j),vf(j),vh(j),'    '
      enddo
      do j=je-2,je
        if(dowr) write(outfile,122) j,yf(j),yh(j),yf(j+1)-yf(j),vf(j),vh(j),'   x'
      enddo
      if(dowr) write(outfile,123) je+1,yf(je+1),vf(je+1)



      if(dowr) write(outfile,*)

134   format('      j         yf           yh         dy         vf         vh')

!--------------------------------------------------------------

      do k=kb,ke+1
      do j=jb,je
      do i=ib,ie
        zf(i,j,k)=dz*(k-1)
        mf(i,j,k)=1.0
      enddo
      enddo
      enddo

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        zh(i,j,k)=0.5*(zf(i,j,k)+zf(i,j,k+1))
        mh(i,j,k)=1.0
      enddo
      enddo
      enddo

      do k=kb,ke+1
        sigmaf(k)=zf(1,1,k)
      enddo
      do k=kb,ke
        sigma(k)=0.5*(sigmaf(k)+sigmaf(k+1))
      enddo

    IF(stretch_z.ge.1)THEN

!-----------------------------------------------------------------------
!  Begin hard-wired analytic stretching function

      strz:  IF ( stretch_z == 1 ) THEN

        nominal_dz = 0.5*(dz_bot+dz_top)

        nk1 = nint( str_bot/dz_bot )
        nk3 = nint( (ztop-str_top)/dz_top )
        nk2 = nk-(nk1+nk3)

        ! dummy checks:
        if(dowr) write(outfile,*) '  bot: ',nk1*dz_bot,str_bot,nk1*dz_bot-str_bot
        if( abs(nk1*dz_bot-str_bot).gt.0.01 )then
          if(dowr) write(outfile,*) '  depth of bottom layer does not exactly divide by dz_bot! '
          if(dowr) write(outfile,*) '  nk1*dz_bot = ',nk1*dz_bot
          if(dowr) write(outfile,*) '  str_bot    = ',str_bot
          if(dowr) write(outfile,*) '  diff       = ',nk1*dz_bot-str_bot
          if(dowr) write(outfile,*) '  stopping cm1 ... '
          call stopcm1
        endif
        if(dowr) write(outfile,*) '  mid: ',nk2*nominal_dz,(str_top-str_bot),nk2*nominal_dz-(str_top-str_bot)
        if( abs(nk2*nominal_dz-(str_top-str_bot)).ge.0.01 )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) amod(str_top-str_bot,nominal_dz),1.0e-6*(str_top-str_bot)
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  depth of middle layer does not exactly divide by nominal_dz! '
          if(dowr) write(outfile,*) '  nk2*nominal_dz  = ',nk2*nominal_dz
          if(dowr) write(outfile,*) '  str_top-str_bot = ',str_top-str_bot
          if(dowr) write(outfile,*) '  diff            = ',nk2*nominal_dz-(str_top-str_bot)
          if(dowr) write(outfile,*) '  stopping cm1 ... '
          call stopcm1
        endif
        if(dowr) write(outfile,*) '  top: ',nk3*dz_top,(ztop-str_top),nk3*dz_top-(ztop-str_top)
        if( abs(nk3*dz_top-(ztop-str_top)).gt.0.01 )then
          if(dowr) write(outfile,*) '  depth of top layer does not exactly divide by dz_top! '
          if(dowr) write(outfile,*) '  nk3*dz_top   = ',nk3*dz_top
          if(dowr) write(outfile,*) '  ztop-str_top = ',ztop-str_top
          if(dowr) write(outfile,*) '  diff         = ',nk3*dz_top-(ztop-str_top)
          if(dowr) write(outfile,*) '  stopping cm1 ... '
          call stopcm1
        endif
        if( (nk1+nk2+nk3)-nk .ne. 0 )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  (nk1+nk2+nk3) does not equal nk '
          if(dowr) write(outfile,*) '  (nk1+nk2+nk3) = ',(nk1+nk2+nk3)
          if(dowr) write(outfile,*) '   nk           = ',nk
          if(dowr) write(outfile,*)
          call stopcm1
        endif

        nominal_dz=(str_top-str_bot)/nk2

        c2=(nominal_dz-dz_bot)/(nominal_dz*nominal_dz*float(nk2-1))
        c1=(dz_bot/nominal_dz)-c2*nominal_dz

        ! Test to see if nk is kosher.
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  actual-nk,test-nk:',float(nk),float(nk1+nk3)+(str_top-str_bot)/(0.5*(dz_bot+dz_top))
        if( abs(float(nk)-(float(nk1+nk3)+(str_top-str_bot)/(0.5*(dz_bot+dz_top)))).gt.1.0e-3  &
              .or. nk1.lt.0 .or. nk2.lt.0 .or. nk3.lt.0 )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of nz = ',nz
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '       nk1,nk2,nk3 = ',nk1,nk2,nk3
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value needed for these settings:'
          if(dowr) write(outfile,*) '       ztop      = ',ztop
          if(dowr) write(outfile,*) '       str_bot   = ',str_bot
          if(dowr) write(outfile,*) '       str_top   = ',str_top
          if(dowr) write(outfile,*) '       dz_bot    = ',dz_bot
          if(dowr) write(outfile,*) '       dz_top    = ',dz_top
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  would be nz = ',nk1+nk3+(str_top-str_bot)/(0.5*(dz_bot+dz_top))
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
          call stopcm1
        endif

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  nk1,nk2,nk3,ntot=',nk1,nk2,nk3,(nk1+nk2+nk3)
        if(dowr) write(outfile,*) '  nominal_dz =',nominal_dz
        if(dowr) write(outfile,*) '  c1,c2 = ',c1,c2
        if(dowr) write(outfile,*)

      do j=jb,je
      do i=ib,ie

        do k=1,nk1+1
          zf(i,j,k)=(k-1)*dz_bot
        enddo
        do k=(nk1+1),(nk1+nk2+1)
          zf(i,j,k)=zf(i,j,nk1+1)+(c1+c2*float(k-1-nk1)*nominal_dz)   &
                         *float(k-1-nk1)*nominal_dz
        enddo
        do k=(nk1+nk2+2),(nk1+nk2+nk3+1)
          zf(i,j,k)=zf(i,j,k-1)+dz_top
        enddo

      enddo
      enddo

!!!      if(terrain_flag)then

        do k=1,nk1+1
          sigmaf(k)=(k-1)*dz_bot
        enddo
        do k=(nk1+1),(nk1+nk2+1)
          sigmaf(k)=sigmaf(nk1+1)+(c1+c2*float(k-1-nk1)*nominal_dz)   &
                         *float(k-1-nk1)*nominal_dz
        enddo
        do k=(nk1+nk2+2),(nk1+nk2+nk3+1)
          sigmaf(k)=sigmaf(k-1)+dz_top
        enddo

        sigmaf(0)=-sigmaf(2)
        sigmaf(nk+2)=sigmaf(nk+1)+(sigmaf(nk+1)-sigmaf(nk))

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ELSEIF ( stretch_z == 2 ) THEN ! geometric stretching from COMMAS

! nz is ke (number of w points)
! kb = 0
! so ng = 1 and zfw1d(-ng+1:nz+ng)
! SUBROUTINE ZGRID(dz, dz_stretch, nz, nbndlyr, gzc, gze, ng, dzmax,ztopstr,rtop,dzmaxtop)
       
       nbndlyr = Int( str_bot/dz_bot + 0.01) - 1
       CALL ZGRID(dz, dz_bot, ke, nbndlyr, gzc, gze, 1, dz_top,ztopstr,rtop,dzmaxtop)
       
         IF ( myid == 0 ) THEN
           DO k = 1,ke
             write(6,*) 'k,gzc,gze = ',k,gzc(k),gze(k)
           ENDDO
         ENDIF
         gze(0) = -gze(2)
         gze(ke+1) = gze(ke) + (gze(ke) - gze(ke-1))

        DO k = kb,ke+1
         DO j = jb,je
          DO i = ib,ie
           zf(i,j,k) = gze(k)
          ENDDO
         ENDDO
        ENDDO

        if(terrain_flag)then
          write(0,*) 'terrain not yet compatible with stretch_z == 2'
          STOP
        ENDIF

        do k=kb,ke+1
          sigmaf(k)=zf(1,1,k)
        enddo
        do k=kb,ke
          sigma(k)=0.5*(sigmaf(k)+sigmaf(k+1))
        enddo

      ENDIF  strz

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!  End hard-wired analytic stretching function
!-----------------------------------------------------------------------
!
!  Optional:  to use a different stretching function, or to use 
!  arbitrarily located grid points, simply comment out the 
!  "hard-wired" section above, and then specify values for zf
!  here.  Do not change anything below here!
!
!  Note:  zf stores the location of the staggered w points. 
!
!  Note:  if you are using terrain, you need to also specify the nominal 
!  locations of the zf points in the sigmaf array.
!
!-----------------------------------------------------------------------

      do j=jb,je
      do i=ib,ie

        zf(i,j,0)=-zf(i,j,2)
        zf(i,j,nk+2)=zf(i,j,nk+1)+(zf(i,j,nk+1)-zf(i,j,nk))

        do k=0,nk+1
          zh(i,j,k)=0.5*(zf(i,j,k+1)+zf(i,j,k))
          mh(i,j,k)=dz/(zf(i,j,k+1)-zf(i,j,k))
        enddo
        zh(i,j,0)=-zh(i,j,1)
        zh(i,j,nk+1)=zh(i,j,nk)+2.0*(zf(i,j,nk+1)-zh(i,j,nk))

        do k=1,nk+1
          mf(i,j,k)=dz/(zh(i,j,k)-zh(i,j,k-1))
        enddo
        mf(i,j,0)=mf(i,j,1)
        mf(i,j,nk+2)=mf(i,j,nk+1)

      enddo
      enddo

    ENDIF

! end vertical stretching section
!-----------------------------------------------------------------------

      do k=kb,ke
        sigma(k)=0.5*(sigmaf(k)+sigmaf(k+1))
      enddo

      maxz = sigmaf(nk+1)

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        rmh(i,j,k)=1.0/mh(i,j,k)
      enddo
      enddo
      enddo

      do k=kb,ke+1
      do j=jb,je
      do i=ib,ie
        rmf(i,j,k)=1.0/mf(i,j,k)
      enddo
      enddo
      enddo

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'model heights:'
      if(dowr) write(outfile,104)
104   format('     k       zf         zh         dz         mf         mh')
      if(dowr) write(outfile,105)
105   format(' ---------------------------------------------------------------')
      do k=1,nk
        if(dowr) write(outfile,102) k,zf(1,1,k),zh(1,1,k),zf(1,1,k+1)-zf(1,1,k),mf(1,1,k),mh(1,1,k)
102     format(3x,i4,3x,f8.2,3x,f8.2,3x,f8.2,3x,f8.4,3x,f8.4)
      enddo
      if(dowr) write(outfile,103) nk+1,zf(1,1,nk+1),mf(1,1,nk+1)
103   format(3x,i4,3x,f8.2,25x,f8.4)
      if(dowr) write(outfile,*)

!-----------------------------------------------------------------------



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                  BEGIN TERRAIN !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      rds  = 1.0
      rdsf = 1.0

      zs=0.0

      gz=1.0
      rgz=1.0
      gzu=1.0
      rgzu=1.0
      gzv=1.0
      rgzv=1.0
      dzdx=0.0
      dzdy=0.0

      gx=0.0
      gxu=0.0
      gy=0.0
      gyv=0.0

      IF(terrain_flag)THEN

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Terrain included!'
        if(dowr) write(outfile,*)

        ! moved this section of code to init_terrain in cm1r15:
        call init_terrain(xh,uh,xf,uf,yh,vh,yf,vf,rds,sigma,rdsf,sigmaf,  &
                          zh,zf,zs,gz,rgz,gzu,rgzu,gzv,rgzv,         &
                          dzdx,dzdy,gx,gxu,gy,gyv,                   &
                          reqs_u,reqs_v,reqs_s,reqs_p,               &
                          nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,           &
                          sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,   &
                          uw31,uw32,ue31,ue32,us31,us32,un31,un32,   &
                          vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,   &
                          ww31(1,1,1),ww32(1,1,1),we31(1,1,1),we32(1,1,1), &
                          ws31(1,1,1),ws32(1,1,1),wn31(1,1,1),wn32(1,1,1))

      ENDIF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                  END   TERRAIN !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        mh(i,j,k)=dz/(zf(i,j,k+1)-zf(i,j,k))
        rmh(i,j,k)=1.0/mh(i,j,k)
      enddo
      enddo
      enddo

      do k=kb+1,ke
      do j=jb,je
      do i=ib,ie
        mf(i,j,k)=dz/(zh(i,j,k)-zh(i,j,k-1))
      enddo
      enddo
      enddo

      do j=jb,je
      do i=ib,ie
        mf(i,j,0)=mf(i,j,1)
        mf(i,j,nk+2)=mf(i,j,nk+1)
      enddo
      enddo

      do k=kb,ke+1
      do j=jb,je
      do i=ib,ie
        rmf(i,j,k)=1.0/mf(i,j,k)
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------

      if(dowr) write(outfile,*) '  minx = ',minx
      if(dowr) write(outfile,*) '  maxx = ',maxx
      if(dowr) write(outfile,*) '  miny = ',miny
      if(dowr) write(outfile,*) '  maxy = ',maxy
      if(dowr) write(outfile,*) '  maxz = ',maxz
      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) '  ibw =',ibw
      if(dowr) write(outfile,*) '  ibe =',ibe
      if(dowr) write(outfile,*) '  ibs =',ibs
      if(dowr) write(outfile,*) '  ibn =',ibn
      if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!  Get min/max dx,dy,dz on grid
!  (needed for adapt_dt ... but interesting to report, nontheless)

      min_dx = 1.0e20
      min_dy = 1.0e20
      min_dz = 1.0e20

      max_dx = 0.0
      max_dy = 0.0
      max_dz = 0.0

      do i=1,ni
        min_dx = min( min_dx , xf(i+1)-xf(i) )
        max_dx = max( max_dx , xf(i+1)-xf(i) )
      enddo

      do j=1,nj
        min_dy = min( min_dy , yf(j+1)-yf(j) )
        max_dy = max( max_dy , yf(j+1)-yf(j) )
      enddo

      do k=1,nk
      do j=1,nj
      do i=1,ni
        min_dz = min( min_dz , zf(i,j,k+1)-zf(i,j,k) )
        max_dz = max( max_dz , zf(i,j,k+1)-zf(i,j,k) )
      enddo
      enddo
      enddo



      if(dowr) write(outfile,*) '  min_dx = ',min_dx
      if(dowr) write(outfile,*) '  max_dx = ',max_dx
      if(dowr) write(outfile,*) '  min_dy = ',min_dy
      if(dowr) write(outfile,*) '  max_dy = ',max_dy
      if(dowr) write(outfile,*) '  min_dz = ',min_dz
      if(dowr) write(outfile,*) '  max_dz = ',max_dz
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  new (cm1r16) arrays for vertical interpolation:

      do k=1,nk+1
      do j=jb,je
      do i=ib,ie
        cc2(i,j,k)=(zf(i,j,k)-zh(i,j,k-1))/(zh(i,j,k)-zh(i,j,k-1))
      enddo
      enddo
      enddo

      call bcs(cc2)





      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  k,c1,c2,zhm1,zf,zh:'
      do k=1,nk+1
      do j=jb,je
      do i=ib,ie
        cc1(i,j,k)=1.0-cc2(i,j,k)
        if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,141) k,cc1(i,j,k),cc2(i,j,k),zh(i,j,k-1),zf(i,j,k),zh(i,j,k)
141     format(3x,i4,2(3x,f7.4),3(3x,f8.2))
      enddo
      enddo
      enddo
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  Specify coefficient for Rayleigh damper in vertical

      if( (irdamp.eq.1).and.(zd.lt.maxz) )then

        IF( zd.lt.(0.5*maxz) )THEN
          if(myid.eq.0)then
          print *
          print *,'  Warning:  with these settings, Rayleigh damping would  '
          print *,'  be applied over MORE than half the domain '
          print *
          print *,'  zd,maxz = ',zd,maxz
          print *
          print *,'   stopping model .... '
          print *
          endif



          call stopcm1
        ENDIF

        do j=jb,je
        do i=ib,ie
          do k=1,nk
            if(zh(i,j,k).gt.zd)then
            tauh(i,j,k)=0.5*(1.0-cos(pi*(zh(i,j,k)-zd)/(zf(i,j,nk+1)-zd)))
            taus(i,j,k)=tauh(i,j,k)
            endif
          enddo
          enddo
        enddo
 
        do j=jb,je
        do i=ib,ie
          do k=1,nk+1
            if(zf(i,j,k).gt.zd)then
            tauf(i,j,k)=0.5*(1.0-cos(pi*(zf(i,j,k)-zd)/(zf(i,j,nk+1)-zd)))
            endif
          enddo
          enddo
        enddo

      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  ------ tauf, tauh -----'
      do k=1,nk
        if(dowr) write(outfile,*) k,tauf(1,1,k),tauh(1,1,k)
      enddo
      if(dowr) write(outfile,*) nk+1,tauf(1,1,nk+1)
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  Rayleigh damping near lateral boundaries:

      IF(hrdamp.eq.1)THEN

        IF( nx.gt.1 )THEN
        IF( xhd.gt.(0.5*(maxx-minx)) )THEN
          if(myid.eq.0)then
          print *
          print *,'  Warning:  with these settings, Rayleigh damping would  '
          print *,'  be applied over MORE than half the domain '
          print *
          print *,'  xhd,minx,maxx = ',xhd,minx,maxx
          print *
          print *,'   stopping model .... '
          print *
          endif



          call stopcm1
        ENDIF
        ENDIF
        IF( ny.gt.1 )THEN
        IF( xhd.gt.(0.5*(maxy-miny)) )THEN
          if(myid.eq.0)then
          print *
          print *,'  Warning:  with these settings, Rayleigh damping would  '
          print *,'  be applied over MORE than half the domain '
          print *
          print *,'  xhd,miny,maxy = ',xhd,miny,maxy
          print *
          print *,'   stopping model .... '
          print *
          endif



          call stopcm1
        ENDIF
        ENDIF

        do k=1,nk
        do j=jb,je
        do i=ib,ie
          ! skip this section of code for 2d simulations:
          IF(nx.gt.1)THEN
            ! west boundary:
            IF( axisymm.ne.1 )THEN
              x1 = (xhd+minx)-xh(i)
              if( x1.gt.0.0 )then
                tauh(i,j,k) = max( tauh(i,j,k) , 0.5*(1.0-cos(pi*x1/xhd)) )
                tauf(i,j,k) = max( tauf(i,j,k) , 0.5*(1.0-cos(pi*x1/xhd)) )
              endif
            ENDIF
            ! east boundary:
            x2 = xh(i)-(maxx-xhd)
            if( x2.gt.0.0 )then
              tauh(i,j,k) = max( tauh(i,j,k) , 0.5*(1.0-cos(pi*x2/xhd)) )
              tauf(i,j,k) = max( tauf(i,j,k) , 0.5*(1.0-cos(pi*x2/xhd)) )
            endif
          ENDIF
          ! skip this section of code for 2d simulations:
          IF(ny.gt.1)THEN
            ! south boundary:
            y1 = (xhd+miny)-yh(j)
            if( y1.gt.0.0 )then
              tauh(i,j,k) = max( tauh(i,j,k) , 0.5*(1.0-cos(pi*y1/xhd)) )
              tauf(i,j,k) = max( tauf(i,j,k) , 0.5*(1.0-cos(pi*y1/xhd)) )
            endif
            ! north boundary:
            y2 = yh(j)-(maxy-xhd)
            if( y2.gt.0.0 )then
              tauh(i,j,k) = max( tauh(i,j,k) , 0.5*(1.0-cos(pi*y2/xhd)) )
              tauf(i,j,k) = max( tauf(i,j,k) , 0.5*(1.0-cos(pi*y2/xhd)) )
            endif
          ENDIF
        enddo
        enddo
        enddo

        IF( nx.gt.1 )THEN
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ------ tauh for horizontal Rayleigh damping (W-to-E) -----'
          j = nint(0.5*float(nj))
          j = max( j , 1 )
          j = min( j , nj )
          if(dowr) write(outfile,*) '  i,tauh:'
          do i=0,ni+1
            if(dowr) write(outfile,*) i,tauh(i,j,1)
          enddo
          if(dowr) write(outfile,*)
        ENDIF

        IF( ny.gt.1 )THEN
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ------ tauh for horizontal Rayleigh damping (S-to-N) -----'
          i = nint(0.5*float(ni))
          i = max( i , 1 )
          i = min( i , ni )
          if(dowr) write(outfile,*) '  j,tauh:'
          do j=0,nj+1
            if(dowr) write(outfile,*) j,tauh(i,j,1)
          enddo
          if(dowr) write(outfile,*)
        ENDIF

      ENDIF

!--------------------------------------------------------------
!  vertically implicit turbulent diffusion:

      ! Set vialpha:
      !      0.0 = forward-in-time (unstable if K dt / (dz^2) > 0.5)
      !      0.5 = centered-in-time (Crank-Nicholson) (stable but oscillatory)
      !      1.0 = backward-in-time (stable)
!      vialpha = 1.0

      ! Do not change this:
!      vibeta  = 1.0 - vialpha

!      NOTE:  these are now set in constants.incl files

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  vialpha,vibeta = ',vialpha,vibeta
        if(dowr) write(outfile,*)

!--------------------------------------------------------------

      dt = dtl
      dtlast = dt

      deallocate( xfdp )
      deallocate( yfdp )

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Get 2nd-order extrapolation coefficients (Fornberg 1988)   ccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  DO nn=1,2

    IF(nn.eq.1)THEN
      x0 = sigmaf(1)/(sigmaf(2)-sigmaf(1))
      alpha(0) = sigma(1)/(sigmaf(2)-sigmaf(1))
      alpha(1) = sigma(2)/(sigmaf(2)-sigmaf(1))
      alpha(2) = sigma(3)/(sigmaf(2)-sigmaf(1))
      alpha(3) = sigma(4)/(sigmaf(2)-sigmaf(1))
    ELSE
      x0 = sigmaf(nk+1)/(sigmaf(nk+1)-sigmaf(nk))
      alpha(0) = sigma(nk  )/(sigmaf(nk+1)-sigmaf(nk))
      alpha(1) = sigma(nk-1)/(sigmaf(nk+1)-sigmaf(nk))
      alpha(2) = sigma(nk-2)/(sigmaf(nk+1)-sigmaf(nk))
      alpha(3) = sigma(nk-3)/(sigmaf(nk+1)-sigmaf(nk))
    ENDIF

      delta = 0.0

      delta(0,0,0) = 1.0
      b1 = 1.0

      do n = 1,bign
        b2 = 1.0
        do nu = 0,n-1
          b3 = alpha(n)-alpha(nu)
          b2 = b2*b3
          if( n.le.bigm ) delta(n-1,n,nu) = 0.0
          do m = 0,min(n,bigm)
            delta(n,m,nu) = ( (alpha(n)-x0)*delta(n-1,m,nu) - m*delta(n-1,m-1,nu) )/b3
          enddo
        enddo
        do m = 0,min(n,bigm)
          delta(n,m,n) = (b1/b2)*( m*delta(n-1,m-1,n-1) - (alpha(n-1)-x0)*delta(n-1,m,n-1) )
        enddo
        b1 = b2
      enddo

    IF(nn.eq.1)THEN
      cgs1 = delta(2,0,0)
      cgs2 = delta(2,0,1)
      cgs3 = delta(2,0,2)
      var = cgs1*sigma(1)+cgs2*sigma(2)+cgs3*sigma(3)
      dgs1 = delta(2,1,0)
      dgs2 = delta(2,1,1)
      dgs3 = delta(2,1,2)
    ELSE
      cgt1 = delta(2,0,0)
      cgt2 = delta(2,0,1)
      cgt3 = delta(2,0,2)
      var = cgt1*sigma(nk)+cgt2*sigma(nk-1)+cgt3*sigma(nk-2)
      dgt1 = delta(2,1,0)
      dgt2 = delta(2,1,1)
      dgt3 = delta(2,1,2)
    ENDIF

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  x0,alpha,delta,sum,predict,maxz:'
      if(dowr) write(outfile,*) sngl(x0)
      if(dowr) write(outfile,*) sngl(alpha(0)),sngl(alpha(1)),sngl(alpha(2)),sngl(alpha(3))
      if(dowr) write(outfile,*) sngl(delta(2,0,0)),sngl(delta(2,0,1)),sngl(delta(2,0,2)),sngl(delta(2,0,3))
      if(dowr) write(outfile,*) sngl(delta(2,0,0)+delta(2,0,1)+delta(2,0,2)+delta(2,0,3)),var,maxz
      if(dowr) write(outfile,*)

      IF( nn.eq.1 )then
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  cgs1,cgs2,cgs3 = ',cgs1,cgs2,cgs3
        if(dowr) write(outfile,*) '  sum            = ',cgs1+cgs2+cgs3
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  dgs1,dgs2,dgs3 = ',dgs1,dgs2,dgs3
        if(dowr) write(outfile,*) '  sum            = ',dgs1+dgs2+dgs3
        if(dowr) write(outfile,*)
      ELSEIF( nn.eq.2 )THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  cgt1,cgt2,cgt3 = ',cgt1,cgt2,cgt3
        if(dowr) write(outfile,*) '  sum            = ',cgt1+cgt2+cgt3
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  dgt1,dgt2,dgt3 = ',dgt1,dgt2,dgt3
        if(dowr) write(outfile,*) '  sum            = ',dgt1+dgt2+dgt3
        if(dowr) write(outfile,*)
      ENDIF

  ENDDO

!--------------------------------------------------------------
!  cm1r18:  Set ghost points for zh
!  [Note:  since cm1r17, the array index (i,j,0) means the surface]
!     (upper/lower ghost points are used by parcel subroutines only)

    DO j=jb,je
    DO i=ib,ie
      zh(i,j,0) = zf(i,j,1)
      zh(i,j,nk+1) = zf(i,j,nk+1)
    ENDDO
    ENDDO

!--------------------------------------------------------------

      if(dowr) write(outfile,*) 'Leaving PARAM'

      return

8000  print *
      print *,'  8000: error opening namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8001  print *
      print *,'  8001: error reading param1 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8002  print *
      print *,'  8002: error reading param2 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8003  print *
      print *,'  8003: error reading param3 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8004  print *
      print *,'  8004: error reading param4 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8005  print *
      print *,'  8005: error reading param5 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8006  print *
      print *,'  8006: error reading param6 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8007  print *
      print *,'  8007: error reading param7 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8008  print *
      print *,'  8008: error reading param8 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8009  print *
      print *,'  8009: error reading param9 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8010  print *
      print *,'  8010: error reading param10 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8011  print *
      print *,'  8011: error reading param11 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8012  print *
      print *,'  8012: error reading param12 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8013  print *
      print *,'  8013: error reading param13 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8051  print *
      print *,'  8051: error reading nssl2mom_params section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

      end subroutine param


!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE ZGRID  <<<<<<<<<<<<<<<<<<<<<<<<<<< !
!
!-------------------------------------------------------------------------------
! Creates the vertical grid using a geometric stretch 
!
! Option - set nbndlyr > 0 input deck.  A suggested value is ~ nz/4 to create a 
! layer having constant resolution near the surface.
!
!-------------------------------------------------------------------------------

 SUBROUTINE ZGRID(dz, dz_stretch, nz, nbndlyr, gzc, gze, ng, dzmax,ztopstr,rtop,dzmaxtop)

   implicit none


   real dz, dz_stretch
   real  dzmax,ztopstr,rtop,dzmaxtop
   integer, intent(in) :: nbndlyr, nz, ng
   double precision gzc(-ng+1:nz+ng), gze(-ng+1:nz+ng)  ! stay   real gzc(nz), gze(nz)  ??

!   real dzmax
!   parameter( dzmax = 700. )
   integer n, k
   double precision stretch, zx, xmid, fmid, ztop

   real zheight
!   external zheight

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

! Use Newton interation to find grid coefficients

   ztop    = dz * (nz-1)
   zx      = 1.0d0
   stretch = 1.0d0

   IF( dz .gt. dz_stretch ) THEN

    DO n = 1,50

     IF( abs(zx) .gt. 1.0e-12 ) THEN
      zx   = zx * 0.5
      xmid = stretch + zx
      fmid = ZHEIGHT(dz_stretch,xmid,nz-1,dzmax,nbndlyr,ztopstr,rtop,dzmaxtop) - ztop
!      write(6,*) 'Stretch: ',n,zx,xmid,fmid,fmid + ztop
      IF( fmid .le. 0.0 ) stretch = xmid
      IF ( fmid .eq. 0.0d0 ) EXIT
     ENDIF

    ENDDO

   ENDIF

   write(6,*)
   IF( stretch .gt. 1.1 ) THEN
    write(6,*) 'STRETCH FAC TOO BIG! - NUMERICAL ERRORS WILL BE LARGE'
    write(6,*) 'STRETCH FAC  = ',stretch
    write(6,*) 'INCREASE NZ in namelist or increase dz_bot'
    STOP
   ELSE
    write(6,*) 'ZGRID:  STRETCH FAC  = ',stretch
    write(6,*) 'ZGRID:  DOMAIN  HGT  = ',ZHEIGHT(dz_stretch,stretch,nz-1,dzmax,nbndlyr,ztopstr,rtop,dzmaxtop)
   ENDIF

   gze(1) = 0.0
   DO k = 1,nz-1
    gze(k+1) = ZHEIGHT(dz_stretch,stretch,k,dzmax,nbndlyr,ztopstr,rtop,dzmaxtop)
    gzc(k)   = 0.5 * ( gze(k) + gze(k+1) )
!    write(6,*) 'ZGRID: gze,gzc = ',k,gze(k+1),gzc(k)
   ENDDO

   gzc(nz) = 2.*gzc(nz-1) - gzc(nz-2) 


  END SUBROUTINE ZGRID

!--------------------------------------------------------------------------
! FUNCTION ZHEIGHT:  Computes the height of a geometrically stretched grid
!                    with a few wrinkles:  It can have a layer of constant
!                    dz at the bottom 'n1' layers thick, it also limits
!                    the size of dz at the top of the model to be 'dzmax'.

 REAL FUNCTION ZHEIGHT(dzbot,r,nz,dzmax,n1,zctop,ztopr,dzmax2)

  implicit none
  integer nz, n1, k, k2
  integer n2
  real dzbot
  double precision r
  double precision sum
  double precision dznew, dzmaxdp, dzmax2dp
  real dzmax
  real zctop  ! height for upper level stretch
  real ztopr  ! upper level stretch factor
  real dzmax2 ! maximum upper dz
  real dzm

  sum = 0.0d0
  dzmaxdp = dzmax
  dzmax2dp = dzmax2
  
!  zctop = 10000.
!  ztopr = 1.09
!  dzmax2 = 1000. ! 2*dzmax
  
  n2 = 0

  DO k = 1,nz

   IF( k .le. n1 ) THEN
    dznew=dzbot
   ELSE
    k2=k-n1
    dznew = Min(dzbot * r**(k2-1),dzmaxdp)
      IF ( sum .ge. zctop ) THEN
       IF ( n2 .eq. 0 ) dzm = Min(dznew,dzmaxdp)
       n2 = n2 + 1
       dznew = Min(dzm * ztopr**n2, dzmax2)
      ENDIF
   ENDIF
   sum = sum + dznew

  ENDDO

  ZHEIGHT = sum

 RETURN
 END FUNCTION ZHEIGHT


