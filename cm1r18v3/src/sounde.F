

      subroutine sounde(dt,xh,arh1,arh2,uh,ruh,xf,uf,yh,vh,rvh,yf,vf,     &
                        rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,mf,         &
                        pi0,rho0,rr0,rf0,rrf0,th0,rth0,zs,                &
                        gz,rgz,gzu,rgzu,gzv,rgzv,                         &
                        dzdx,dzdy,gx,gxu,gy,gyv,                          &
                        radbcw,radbce,radbcs,radbcn,                      &
                        dum1,dum2,dum3,dum4,dum5,dum6,                    &
                        ppd ,fpk ,qk ,pk1,pk2,ftk,sk ,tk1,tk2,            &
                        u0,rru,us,ua,u3d,uten,                            &
                        v0,rrv,vs,va,v3d,vten,                            &
                        rrw,ws,wa,w3d,wten,                               &
                        ppi,pp3d,piadv,ppten,ppx,                         &
                        tha,th3d,thadv,thten,thterm,                      &
                        thv,ppterm,nrk,dttmp,rtime,get_time_avg,          &
                        pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,arh1,arh2,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,rr0,rf0,rrf0,th0,rth0
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppd,fpk,qk,pk1,pk2,ftk,sk,tk1,tk2
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: rru,us,ua,u3d,uten
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: rrv,vs,va,v3d,vten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw,ws,wa,w3d,wten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,piadv,ppten,ppx
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tha,th3d,thadv,thten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: thterm,thv,ppterm
      integer, intent(in) :: nrk
      real, intent(in)  :: dttmp,rtime
      logical, intent(in) :: get_time_avg
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      integer, intent(inout), dimension(rmp) :: reqs_p

!-----

      integer :: i,j,k,n,nloop
      real :: tem,tem1,tem2,tem3,tem4,r1,r2,dts

      real :: temx,temy,u1,u2,v1,v2,w1,w2,ww,tavg,div

!---------------------------------------------------------------------

      if(nrk.eq.1)then
!!!        nloop=1
!!!        dts=dt/3.
        nloop=nint(float(nsound)/3.0)
        dts=dt/(nloop*3.0)
        if( dts.gt.(dt/nsound) )then
          nloop=nloop+1
          dts=dt/(nloop*3.0)
        endif
      elseif(nrk.eq.2)then
        nloop=0.5*nsound
        dts=dt/nsound
      elseif(nrk.eq.3)then
        nloop=nsound
        dts=dt/nsound
      endif

!!!      print *,'  nloop,dts,dttmp = ',nloop,dts,nloop*dts

!-----------------------------------------------------------------
!  define ppd first, then start comm:

      if( nrk.eq.1 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ppd(i,j,k)=ppi(i,j,k)+ppx(i,j,k)
        enddo
        enddo
        enddo

      else

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ppd(i,j,k)=ppi(i,j,k)+ppx(i,j,k)
          pp3d(i,j,k)=ppi(i,j,k)
        enddo
        enddo
        enddo

      endif
      if(timestats.ge.1) time_sound=time_sound+mytime()

      call bcs(ppd)
#ifdef MPI
      call comm_1s_start(ppd,pw1,pw2,pe1,pe2,   &
                             ps1,ps2,pn1,pn2,reqs_p)
#endif

!---------------------------------------------------------------------
!  Arrays for vadv:

      IF(.not.terrain_flag)THEN

        ! without terrain:
        ! "s" velocities ARE NOT coupled with reference density
!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem,tem1,r1,r2)
        do k=1,nk
          r2 = dts*rdz*mh(1,1,k)*rr0(1,1,k)*rf0(1,1,k+1)
          r1 = dts*rdz*mh(1,1,k)*rr0(1,1,k)*rf0(1,1,k)
          do j=1,nj
          do i=1,ni
            pk2(i,j,k) = r2*( -c2(1,1,k+1)*piadv(i,j,k+1)+(1.0-c1(1,1,k+1))*piadv(i,j,k) )
            pk1(i,j,k) = r1*( +c1(1,1,k  )*piadv(i,j,k-1)+(c2(1,1,k  )-1.0)*piadv(i,j,k) )
            tk2(i,j,k) = r2*( -c2(1,1,k+1)*thadv(i,j,k+1)+(1.0-c1(1,1,k+1))*thadv(i,j,k) )
            tk1(i,j,k) = r1*( +c1(1,1,k  )*thadv(i,j,k-1)+(c2(1,1,k  )-1.0)*thadv(i,j,k) )
          enddo
          enddo
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

      ELSE

        ! with terrain:
        ! "s" velocities ARE coupled with reference density
!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem,tem1,r1,r2)
        do k=1,nk
          do j=1,nj
          do i=1,ni
            r2 = dts*gz(i,j)*rdsf(k)*rr0(i,j,k)*rf0(i,j,k+1)
            r1 = dts*gz(i,j)*rdsf(k)*rr0(i,j,k)*rf0(i,j,k)
            pk2(i,j,k) = r2*( -c2(i,j,k+1)*piadv(i,j,k+1)+(1.0-c1(i,j,k+1))*piadv(i,j,k) )
            pk1(i,j,k) = r1*( +c1(i,j,k  )*piadv(i,j,k-1)+(c2(i,j,k  )-1.0)*piadv(i,j,k) )
            tk2(i,j,k) = r2*( -c2(i,j,k+1)*thadv(i,j,k+1)+(1.0-c1(i,j,k+1))*thadv(i,j,k) )
            tk1(i,j,k) = r1*( +c1(i,j,k  )*thadv(i,j,k-1)+(c2(i,j,k  )-1.0)*thadv(i,j,k) )
          enddo
          enddo
          do j=1,nj
          do i=1,ni+1
            us(i,j,k)=rru(i,j,k)
          enddo
          enddo
          IF(axisymm.eq.0)THEN
            ! Cartesian grid:
            do j=1,nj+1
            do i=1,ni
              vs(i,j,k)=rrv(i,j,k)
            enddo
            enddo
          ENDIF
          IF( k.ge.2 )THEN
            do j=1,nj
            do i=1,ni
              ws(i,j,k)=rrw(i,j,k)
            enddo
            enddo
          ELSE
            do j=1,nj
            do i=1,ni
              dum3(i,j,1)=0.0
              dum3(i,j,nk+1)=0.0
            enddo
            enddo
          ENDIF
        enddo

      ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        pk1(i,j,1) = 0.0
        tk1(i,j,1) = 0.0
        pk2(i,j,nk) = 0.0
        tk2(i,j,nk) = 0.0
      enddo
      enddo

!---------------------------------------------------------------------
!  Prepare for acoustic steps

      if( nrk.ge.2 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          do j=1,nj
          do i=1,ni+1
            u3d(i,j,k)=ua(i,j,k)
          enddo
          enddo
          IF(axisymm.eq.0)THEN
            ! Cartesian grid:
            do j=1,nj+1
            do i=1,ni
              v3d(i,j,k)=va(i,j,k)
            enddo
            enddo
          ENDIF
          IF(k.ge.2)THEN
            do j=1,nj
            do i=1,ni
              w3d(i,j,k)=wa(i,j,k)
            enddo
            enddo
          ENDIF
          do j=1,nj
          do i=1,ni
            th3d(i,j,k)=tha(i,j,k)
          enddo
          enddo
        enddo

      endif

!---------------------------------------------------------------------

      IF( get_time_avg )THEN
        tavg = 1.0/float(nloop)
!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem)
          do k=1,nk
            do j=1,nj
            do i=1,ni+1
              rru(i,j,k)=0.0
            enddo
            enddo
            IF(axisymm.eq.0)THEN
              ! Cartesian grid:
              do j=1,nj+1
              do i=1,ni
                rrv(i,j,k)=0.0
              enddo
              enddo
            ENDIF
            IF(k.ge.2)THEN
              do j=1,nj
              do i=1,ni
                rrw(i,j,k)=0.0
              enddo
              enddo
            ENDIF
          enddo
    ENDIF      ! endif for get_time_avg check


      if(timestats.ge.1) time_sound=time_sound+mytime()

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  Begin small steps:

      small_step_loop:  DO N=1,NLOOP

!-----

        if(irbc.eq.2)then
 
          if(ibw.eq.1 .or. ibe.eq.1) call radbcew(radbcw,radbce,u3d)
 
          if(ibs.eq.1 .or. ibn.eq.1) call radbcns(radbcs,radbcn,v3d)
 
        endif

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  Open boundary conditions:

        IF(wbc.eq.2.and.ibw.eq.1)THEN
          ! west open bc tendency:
          call   ssopenbcw(uh,rds,sigma,rdsf,sigmaf,gz,rgzu,gx,radbcw,dum1,u3d,uten,dts)
        ENDIF
        IF(ebc.eq.2.and.ibe.eq.1)THEN
          ! east open bc tendency:
          call   ssopenbce(uh,rds,sigma,rdsf,sigmaf,gz,rgzu,gx,radbce,dum1,u3d,uten,dts)
        ENDIF

        IF(roflux.eq.1)THEN
          call restrict_openbc_we(rvh,rmh,rho0,u3d)
        ENDIF

!-----

      IF(axisymm.eq.0)THEN
        IF(sbc.eq.2.and.ibs.eq.1)THEN
          ! south open bc tendency:
          call   ssopenbcs(vh,rds,sigma,rdsf,sigmaf,gz,rgzv,gy,radbcs,dum1,v3d,vten,dts)
        ENDIF
        IF(nbc.eq.2.and.ibn.eq.1)THEN
          ! north open bc tendency:
          call   ssopenbcn(vh,rds,sigma,rdsf,sigmaf,gz,rgzv,gy,radbcn,dum1,v3d,vten,dts)
        ENDIF

        IF(roflux.eq.1)THEN
          call restrict_openbc_sn(ruh,rmh,rho0,v3d)
        ENDIF
      ENDIF

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  integrate u,v forward in time:

#ifdef MPI
          call comm_1s_end(ppd,pw1,pw2,pe1,pe2,   &
                               ps1,ps2,pn1,pn2,reqs_p)
#endif

!-----

    IF(.not.terrain_flag)THEN

      IF(axisymm.eq.0)THEN
        ! Cartesian grid without terrain:

        tem1 = rdx*cp*0.5
        tem2 = rdy*cp*0.5
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=u3d(i,j,k)+dts*( uten(i,j,k)         &
                   -tem1*(ppd(i,j,k)-ppd(i-1,j,k))*uf(i)    &
                        *(thv(i,j,k)+thv(i-1,j,k)) )
          enddo
          enddo
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            v3d(i,j,k)=v3d(i,j,k)+dts*( vten(i,j,k)         &
                   -tem2*(ppd(i,j,k)-ppd(i,j-1,k))*vf(j)    &
                        *(thv(i,j,k)+thv(i,j-1,k)) )
          enddo
          enddo
        enddo

      ELSE
        ! axisymmetric grid:

        tem1 = rdx*cp*0.5
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=u3d(i,j,k)+dts*( uten(i,j,k)         &
                   -tem1*(ppd(i,j,k)-ppd(i-1,j,k))*uf(i)    &
                        *(thv(i,j,k)+thv(i-1,j,k)) )
          enddo
          enddo
        enddo

      ENDIF

    ELSE

        ! Cartesian grid with terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
        do j=0,nj+1
          ! dum1 stores ppd at w-pts:
          ! lowest model level:
          do i=0,ni+1
            dum1(i,j,1) = cgs1*ppd(i,j,1)+cgs2*ppd(i,j,2)+cgs3*ppd(i,j,3)
          enddo
          ! upper-most model level:
          do i=0,ni+1
            dum1(i,j,nk+1) = cgt1*ppd(i,j,nk)+cgt2*ppd(i,j,nk-1)+cgt3*ppd(i,j,nk-2)
          enddo
          ! interior:
          do k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do i=0,ni+1
            dum1(i,j,k) = r1*ppd(i,j,k-1)+r2*ppd(i,j,k)
          enddo
          enddo
        enddo

        tem1 = rdx*cp*0.5
        tem2 = rdy*cp*0.5

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          ! x-dir
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=u3d(i,j,k)+dts*( uten(i,j,k)               &
                   -cp*0.5*(thv(i,j,k)+thv(i-1,j,k))*(            &
                   ( ppd(i  ,j,k)*rgz(i  ,j)                      &
                    -ppd(i-1,j,k)*rgz(i-1,j)                      &
                   )*gzu(i,j)*rdx*uf(i)                           &
              +0.5*( gxu(i,j,k+1)*(dum1(i,j,k+1)+dum1(i-1,j,k+1)) &
                    -gxu(i,j,k  )*(dum1(i,j,k  )+dum1(i-1,j,k  )) &
                   )*rdsf(k) ) )
          enddo
          enddo
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            v3d(i,j,k)=v3d(i,j,k)+dts*( vten(i,j,k)               &
                   -cp*0.5*(thv(i,j,k)+thv(i,j-1,k))*(            &
                   ( ppd(i,j  ,k)*rgz(i,j  )                      &
                    -ppd(i,j-1,k)*rgz(i,j-1)                      &
                   )*gzv(i,j)*rdy*vf(j)                           &
              +0.5*( gyv(i,j,k+1)*(dum1(i,j,k+1)+dum1(i,j-1,k+1)) &
                    -gyv(i,j,k  )*(dum1(i,j,k  )+dum1(i,j-1,k  )) &
                   )*rdsf(k) ) )
          enddo
          enddo
        enddo

    ENDIF

!----------------------------------------------
!  convergence forcing:

        IF( convinit.eq.1 )THEN
          IF( rtime.le.convtime .and. nx.gt.1 )THEN
            call convinitu(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibw,ibe,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xf,yh,zh,u0,u3d)
          ENDIF
        ENDIF

!----------------------------------------------

      IF(axisymm.eq.0)THEN
        ! Cartesian grid:

!----------------------------------------------
!  convergence forcing:

        IF( convinit.eq.1 )THEN
          IF( rtime.le.convtime .and. ny.gt.1 )THEN
            call convinitv(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibs,ibn,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xh,yf,zh,v0,v3d)
          ENDIF
        ENDIF

!----------------------------------------------

      ENDIF

      if(timestats.ge.1) time_sound=time_sound+mytime()

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  integrate w forward in time:

      IF(.not.terrain_flag)THEN
        ! without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem1)
        do k=2,nk
        tem1 = rdz*cp*mf(1,1,k)
        do j=1,nj
        do i=1,ni
          w3d(i,j,k)=w3d(i,j,k)+dts*( wten(i,j,k)                     &
                  -tem1*(ppd(i,j,k)-ppd(i,j,k-1))                     &
                       *(c2(1,1,k)*thv(i,j,k)+c1(1,1,k)*thv(i,j,k-1)) &
                +g*( c2(1,1,k)*th3d(i,j,k)*rth0(1,1,k)                &
                    +c1(1,1,k)*th3d(i,j,k-1)*rth0(1,1,k-1) ) )
        enddo
        enddo
        enddo

      ELSE
        ! with terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem1)
        do k=2,nk
        tem1 = rds(k)*cp
        do j=1,nj
        do i=1,ni
          w3d(i,j,k)=w3d(i,j,k)+dts*( wten(i,j,k)                     &
                  -tem1*(ppd(i,j,k)-ppd(i,j,k-1))*gz(i,j)             &
                       *(c2(i,j,k)*thv(i,j,k)+c1(i,j,k)*thv(i,j,k-1)) &
                +g*( c2(i,j,k)*th3d(i,j,k)*rth0(i,j,k)                &
                    +c1(i,j,k)*th3d(i,j,k-1)*rth0(i,j,k-1) ) )
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcwsfc(gz,dzdx,dzdy,u3d,v3d,w3d)

      ENDIF

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  get terms for div,vadv (terrain only):

      IF(terrain_flag)THEN
        ! Cartesian grid with terrain:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        DO k=1,nk
          do j=1,nj
          do i=1,ni+1
            dum1(i,j,k)=u3d(i,j,k)*rgzu(i,j)
            dum4(i,j,k)=0.5*(rho0(i-1,j,k)+rho0(i,j,k))*dum1(i,j,k)
          enddo
          enddo
          do j=1,nj+1
          do i=1,ni
            dum2(i,j,k)=v3d(i,j,k)*rgzv(i,j)
            dum5(i,j,k)=0.5*(rho0(i,j-1,k)+rho0(i,j,k))*dum2(i,j,k)
          enddo
          enddo
        ENDDO
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
        DO k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do j=1,nj
          do i=1,ni
            dum3(i,j,k)=w3d(i,j,k)                                               &
                       +0.5*( ( r2*(dum1(i,j,k  )+dum1(i+1,j,k  ))               &
                               +r1*(dum1(i,j,k-1)+dum1(i+1,j,k-1)) )*dzdx(i,j)   &
                             +( r2*(dum2(i,j,k  )+dum2(i,j+1,k  ))               &
                               +r1*(dum2(i,j,k-1)+dum2(i,j+1,k-1)) )*dzdy(i,j)   &
                            )*(sigmaf(k)-zt)*gz(i,j)*rzt
            ! NOTE:  dum6 is NOT coupled with density
            dum6(i,j,k)=w3d(i,j,k)                                               &
                       +0.5*( ( r2*(dum4(i,j,k  )+dum4(i+1,j,k  ))               &
                               +r1*(dum4(i,j,k-1)+dum4(i+1,j,k-1)) )*dzdx(i,j)   &
                             +( r2*(dum5(i,j,k  )+dum5(i,j+1,k  ))               &
                               +r1*(dum5(i,j,k-1)+dum5(i,j+1,k-1)) )*dzdy(i,j)   &
                            )*(sigmaf(k)-zt)*gz(i,j)*rzt * rrf0(i,j,k)
          enddo
          enddo
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  get new pp,th

      temx = dts*0.5*rdx
      temy = dts*0.5*rdy

!$omp parallel do default(shared)   &
!$omp private(i,j,k,div,u1,u2,v1,v2,w1,w2,tem)
      DO k=1,nk

    IF( axisymm.eq.0 )THEN
      IF(.not.terrain_flag)THEN
        ! Cartesian grid, without terrain:
        do j=1,nj
        do i=1,ni
          div=(u3d(i+1,j,k)-u3d(i,j,k))*rdx*uh(i)    &
             +(v3d(i,j+1,k)-v3d(i,j,k))*rdy*vh(j)    &
             +(w3d(i,j,k+1)-w3d(i,j,k))*rdz*mh(1,1,k)
          if(abs(div).lt.smeps) div=0.0
          u2 = temx*(u3d(i+1,j,k)-us(i+1,j,k))*uh(i)
          u1 = temx*(u3d(i  ,j,k)-us(i  ,j,k))*uh(i)
          v2 = temy*(v3d(i,j+1,k)-vs(i,j+1,k))*vh(j)
          v1 = temy*(v3d(i,j  ,k)-vs(i,j  ,k))*vh(j)
          w2 = w3d(i,j,k+1)-ws(i,j,k+1)
          w1 = w3d(i,j,k  )-ws(i,j,k  )
          !-----
          ppd(i,j,k)=pp3d(i,j,k)
          pp3d(i,j,k)=pp3d(i,j,k)+dts*( ppten(i,j,k)-ppterm(i,j,k)*div )  &
                 +( -( u2*(piadv(i+1,j,k)-piadv(i  ,j,k))                 &
                      +u1*(piadv(i  ,j,k)-piadv(i-1,j,k)) )               &
                    -( v2*(piadv(i,j+1,k)-piadv(i,j  ,k))                 &
                      +v1*(piadv(i,j  ,k)-piadv(i,j-1,k)) ) )             &
                    +( w1*pk1(i,j,k)+w2*pk2(i,j,k) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          dum1(i,j,k)=kdiv*( pp3d(i,j,k)-ppd(i,j,k) )
          ppd(i,j,k)=pp3d(i,j,k)+dum1(i,j,k)
          !-----
          th3d(i,j,k)=th3d(i,j,k)+dts*( thten(i,j,k)-thterm(i,j,k)*div )  &
                 +( -( u2*(thadv(i+1,j,k)-thadv(i  ,j,k))                 &
                      +u1*(thadv(i  ,j,k)-thadv(i-1,j,k)) )               &
                    -( v2*(thadv(i,j+1,k)-thadv(i,j  ,k))                 &
                      +v1*(thadv(i,j  ,k)-thadv(i,j-1,k)) ) )             &
                    +( w1*tk1(i,j,k)+w2*tk2(i,j,k) )
          if(abs(th3d(i,j,k)).lt.smeps) th3d(i,j,k)=0.0
          !-----
        enddo
        enddo
      ELSE
        ! Cartesian grid, with terrain:
        do j=1,nj
        do i=1,ni
          div = gz(i,j)*( (dum1(i+1,j,k)-dum1(i,j,k))*rdx*uh(i)  &
                         +(dum2(i,j+1,k)-dum2(i,j,k))*rdy*vh(j)  &
                         +(dum3(i,j,k+1)-dum3(i,j,k))*rdsf(k) )
          if(abs(div).lt.smeps) div=0.0
          u2 = temx*(dum4(i+1,j,k)-us(i+1,j,k))*uh(i)
          u1 = temx*(dum4(i  ,j,k)-us(i  ,j,k))*uh(i)
          v2 = temy*(dum5(i,j+1,k)-vs(i,j+1,k))*vh(j)
          v1 = temy*(dum5(i,j  ,k)-vs(i,j  ,k))*vh(j)
          w2 = dum6(i,j,k+1)-ws(i,j,k+1)*rrf0(i,j,k+1)
          w1 = dum6(i,j,k  )-ws(i,j,k  )*rrf0(i,j,k  )
          !-----
          ppd(i,j,k)=pp3d(i,j,k)
          pp3d(i,j,k)=pp3d(i,j,k)+dts*( ppten(i,j,k)-ppterm(i,j,k)*div )  &
                 +( -( u2*(piadv(i+1,j,k)-piadv(i  ,j,k))                 &
                      +u1*(piadv(i  ,j,k)-piadv(i-1,j,k)) )               &
                    -( v2*(piadv(i,j+1,k)-piadv(i,j  ,k))                 &
                      +v1*(piadv(i,j  ,k)-piadv(i,j-1,k)) ) )*rr0(i,j,k)*gz(i,j) &
                    +( w1*pk1(i,j,k)+w2*pk2(i,j,k) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          dum1(i,j,k)=kdiv*( pp3d(i,j,k)-ppd(i,j,k) )
          ppd(i,j,k)=pp3d(i,j,k)+dum1(i,j,k)
          !-----
          th3d(i,j,k)=th3d(i,j,k)+dts*( thten(i,j,k)-thterm(i,j,k)*div )  &
                 +( -( u2*(thadv(i+1,j,k)-thadv(i  ,j,k))                 &
                      +u1*(thadv(i  ,j,k)-thadv(i-1,j,k)) )               &
                    -( v2*(thadv(i,j+1,k)-thadv(i,j  ,k))                 &
                      +v1*(thadv(i,j  ,k)-thadv(i,j-1,k)) ) )*rr0(i,j,k)*gz(i,j) &
                    +( w1*tk1(i,j,k)+w2*tk2(i,j,k) )
          if(abs(th3d(i,j,k)).lt.smeps) th3d(i,j,k)=0.0
          !-----
        enddo
        enddo
      ENDIF
    ELSE
        ! axisymmetric grid:
        do j=1,nj
        do i=1,ni
          div=(arh2(i)*u3d(i+1,j,k)-arh1(i)*u3d(i,j,k))*rdx*uh(i)   &
             +(w3d(i,j,k+1)-w3d(i,j,k))*rdz*mh(1,1,k)
          if(abs(div).lt.smeps) div=0.0
          u2 = temx*(u3d(i+1,j,k)-us(i+1,j,k))*uh(i)*arh2(i)
          u1 = temx*(u3d(i  ,j,k)-us(i  ,j,k))*uh(i)*arh1(i)
          w2 = w3d(i,j,k+1)-ws(i,j,k+1)
          w1 = w3d(i,j,k  )-ws(i,j,k  )
          !-----
          ppd(i,j,k)=pp3d(i,j,k)
          pp3d(i,j,k)=pp3d(i,j,k)+dts*( ppten(i,j,k)-ppterm(i,j,k)*div )  &
                    -( u2*(piadv(i+1,j,k)-piadv(i  ,j,k))                 &
                      +u1*(piadv(i  ,j,k)-piadv(i-1,j,k)) )               &
                    +( w1*pk1(i,j,k)+w2*pk2(i,j,k) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          dum1(i,j,k)=kdiv*( pp3d(i,j,k)-ppd(i,j,k) )
          ppd(i,j,k)=pp3d(i,j,k)+dum1(i,j,k)
          !-----
          th3d(i,j,k)=th3d(i,j,k)+dts*( thten(i,j,k)-thterm(i,j,k)*div )  &
                    -( u2*(thadv(i+1,j,k)-thadv(i  ,j,k))                 &
                      +u1*(thadv(i  ,j,k)-thadv(i-1,j,k)) )               &
                    +( w1*tk1(i,j,k)+w2*tk2(i,j,k) )
          if(abs(th3d(i,j,k)).lt.smeps) th3d(i,j,k)=0.0
          !-----
        enddo
        enddo
    ENDIF

      ENDDO
      if(timestats.ge.1) time_sound=time_sound+mytime()

        IF( n.lt.nloop )THEN
          call bcs(ppd)
#ifdef MPI
          call comm_1s_start(ppd,pw1,pw2,pe1,pe2,   &
                                 ps1,ps2,pn1,pn2,reqs_p)
#endif
        ENDIF

!--------------------------------------------------------------------
!  time-averaged velocities:

    IF( get_time_avg )THEN
      IF(.not.terrain_flag)THEN
        !-----
        ! without terrain:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem)
        DO k=1,nk
          tem = rho0(1,1,k)*tavg
          do j=1,nj
          do i=1,ni+1
            rru(i,j,k)=rru(i,j,k)+u3d(i,j,k)*tem
          enddo
          enddo
          IF( axisymm.eq.0 )THEN
            do j=1,nj+1
            do i=1,ni
              rrv(i,j,k)=rrv(i,j,k)+v3d(i,j,k)*tem
            enddo
            enddo
          ENDIF
          IF( k.ge.2 )THEN
            tem = rf0(1,1,k)*tavg
            do j=1,nj
            do i=1,ni
              rrw(i,j,k)=rrw(i,j,k)+w3d(i,j,k)*tem
            enddo
            enddo
          ENDIF
        ENDDO
        !-----
      ELSE
        !-----
        ! with terrain:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem)
        DO k=1,nk
          do j=1,nj
          do i=1,ni+1
            rru(i,j,k)=rru(i,j,k)+dum4(i,j,k)*tavg
          enddo
          enddo
          IF( axisymm.eq.0 )THEN
            do j=1,nj+1
            do i=1,ni
              rrv(i,j,k)=rrv(i,j,k)+dum5(i,j,k)*tavg
            enddo
            enddo
          ENDIF
          IF( k.ge.2 )THEN
            do j=1,nj
            do i=1,ni
              rrw(i,j,k)=rrw(i,j,k)+dum6(i,j,k)*rf0(i,j,k)*tavg
            enddo
            enddo
          ENDIF
        ENDDO
        !-----
      ENDIF    ! endif for terrain check
    ENDIF      ! endif for get_time_avg check
      if(timestats.ge.1) time_sound=time_sound+mytime()

!--------------------------------------------------------------------

      ENDDO  small_step_loop

!  end of small steps
!--------------------------------------------------------------------

      IF( nrk.eq.3 )THEN
        ! pressure tendency term: save for next timestep:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ppx(i,j,k)=dum1(i,j,k)
        enddo
        enddo
        enddo
      ENDIF
      if(timestats.ge.1) time_sound=time_sound+mytime()


      call bcu(u3d)
      IF(axisymm.eq.0)THEN
        call bcv(v3d)
      ENDIF
      call bcw(w3d,1)
      if(terrain_flag) call bcwsfc(gz,dzdx,dzdy,u3d,v3d,w3d)
      if(nrk.lt.3)then
        call bcs(th3d)
        call bcs(pp3d)
      endif


      end subroutine sounde


