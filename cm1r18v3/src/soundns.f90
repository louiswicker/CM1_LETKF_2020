

      subroutine soundns(xh,rxh,arh1,arh2,uh,xf,uf,yh,vh,yf,vf,           &
                         zh,mh,c1,c2,mf,pi0,thv0,rr0,rf0,                 &
                         rds,sigma,rdsf,sigmaf,                           &
                         zs,gz,rgz,gzu,rgzu,gzv,rgzv,                     &
                         dzdx,dzdy,gx,gxu,gy,gyv,                         &
                         radbcw,radbce,radbcs,radbcn,                     &
                         dum1,dum2,dum3,div ,                             &
                         u0,ua,u3d,uten,                                  &
                         v0,va,v3d,vten,                                  &
                         wa,w3d,wten,                                     &
                         ppi,pp3d,ppten,thv,ppterm,dttmp,nrk,rtime,       &
                         th0,tha,th3d,thten,thterm)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, dimension(ib:ie+1) :: xf,uf
      real, dimension(jb:je) :: yh,vh
      real, dimension(jb:je+1) :: yf,vf
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,mh,c1,c2
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0,rr0,rf0
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, dimension(jb:je,kb:ke) :: radbcw,radbce
      real, dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,div
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u0,ua,u3d,uten
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v0,va,v3d,vten
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa,w3d,wten
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,ppten
      real, dimension(ib:ie,jb:je,kb:ke) :: thv,ppterm
      real, intent(in) :: dttmp
      integer, intent(in) :: nrk
      real, intent(in) :: rtime
      real, dimension(ib:ie,jb:je,kb:ke) :: th0,tha,th3d,thten,thterm

!-----

      integer :: i,j,k
      real :: tem,tem1,tem2,r1,r2

!---------------------------------------------------------------------

        if(irbc.eq.2)then
 
          if(ibw.eq.1 .or. ibe.eq.1) call radbcew(radbcw,radbce,ua)
 
          if(ibs.eq.1 .or. ibn.eq.1) call radbcns(radbcs,radbcn,va)
 
        endif

!-----

        if(wbc.eq.2.and.ibw.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj
            u3d(1,j,k)=ua(1,j,k)
          enddo
          enddo
          call   ssopenbcw(uh,rds,sigma,rdsf,sigmaf,gz,rgzu,gx,radbcw,dum1,u3d,uten,dttmp)
        endif

        if(ebc.eq.2.and.ibe.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj
            u3d(ni+1,j,k)=ua(ni+1,j,k)
          enddo
          enddo
          call   ssopenbce(uh,rds,sigma,rdsf,sigmaf,gz,rgzu,gx,radbce,dum1,u3d,uten,dttmp)
        endif

!-----

      IF(axisymm.eq.0)THEN

        if(sbc.eq.2.and.ibs.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            v3d(i,1,k)=va(i,1,k)
          enddo
          enddo
          call   ssopenbcs(vh,rds,sigma,rdsf,sigmaf,gz,rgzv,gy,radbcs,dum1,v3d,vten,dttmp)
        endif
 
        if(nbc.eq.2.and.ibn.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            v3d(i,nj+1,k)=va(i,nj+1,k)
          enddo
          enddo
          call   ssopenbcn(vh,rds,sigma,rdsf,sigmaf,gz,rgzv,gy,radbcn,dum1,v3d,vten,dttmp)
        endif

      ENDIF

!-----

    IF(.not.terrain_flag)THEN

      IF(axisymm.eq.0)THEN

        tem1=dttmp*rdx*cp*0.5
        tem2=dttmp*rdy*cp*0.5

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=ua(i,j,k)+dttmp*uten(i,j,k)             &
                   -(tem1*(pp3d(i,j,k)-pp3d(i-1,j,k))*uf(i)    &
                         *(thv(i,j,k)+thv(i-1,j,k)))
          enddo
          enddo
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            v3d(i,j,k)=va(i,j,k)+dttmp*vten(i,j,k)             &
                   -(tem2*(pp3d(i,j,k)-pp3d(i,j-1,k))*vf(j)    &
                         *(thv(i,j,k)+thv(i,j-1,k)))
          enddo
          enddo
        enddo

      ELSE

        tem1=dttmp*rdx*cp*0.5

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1+ibw,ni+1-ibe
          u3d(i,j,k)=ua(i,j,k)+dttmp*uten(i,j,k)             &
                 -(tem1*(pp3d(i,j,k)-pp3d(i-1,j,k))*uf(i)    &
                       *(thv(i,j,k)+thv(i-1,j,k)))
        enddo
        enddo
        enddo

      ENDIF

    ELSE

        ! Cartesian grid with terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
        do j=0,nj+1
          ! dum1 stores pp3d at w-pts:
          ! lowest model level:
          do i=0,ni+1
            dum1(i,j,1) = cgs1*pp3d(i,j,1)+cgs2*pp3d(i,j,2)+cgs3*pp3d(i,j,3)
          enddo
          ! upper-most model level:
          do i=0,ni+1
            dum1(i,j,nk+1) = cgt1*pp3d(i,j,nk)+cgt2*pp3d(i,j,nk-1)+cgt3*pp3d(i,j,nk-2)
          enddo
          ! interior:
          do k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do i=0,ni+1
            dum1(i,j,k) = r1*pp3d(i,j,k-1)+r2*pp3d(i,j,k)
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
            u3d(i,j,k)=ua(i,j,k)+dttmp*( uten(i,j,k)              &
                   -cp*0.5*(thv(i,j,k)+thv(i-1,j,k))*(            &
                   ( pp3d(i  ,j,k)*rgz(i  ,j)                     &
                    -pp3d(i-1,j,k)*rgz(i-1,j)                     &
                   )*gzu(i,j)*rdx*uf(i)                           &
              +0.5*( gxu(i,j,k+1)*(dum1(i,j,k+1)+dum1(i-1,j,k+1)) &
                    -gxu(i,j,k  )*(dum1(i,j,k  )+dum1(i-1,j,k  )) &
                   )*rdsf(k) ) )
          enddo
          enddo
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            v3d(i,j,k)=va(i,j,k)+dttmp*( vten(i,j,k)              &
                   -cp*0.5*(thv(i,j,k)+thv(i,j-1,k))*(            &
                   ( pp3d(i,j  ,k)*rgz(i,j  )                     &
                    -pp3d(i,j-1,k)*rgz(i,j-1)                     &
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

        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcu(u3d)

!----------------------------------------------
!  convergence forcing:

      IF(axisymm.eq.0)THEN

        IF( convinit.eq.1 )THEN
          IF( rtime.le.convtime .and. ny.gt.1 )THEN
            call convinitv(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibs,ibn,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xh,yf,zh,v0,v3d)
          ENDIF
        ENDIF

!----------------------------------------------

        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcv(v3d)

      ENDIF

!-----

    IF(.not.terrain_flag)THEN
        ! without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem1)
        do k=2,nk
        tem1 = rdz*cp*mf(1,1,k)
        do j=1,nj
        do i=1,ni
          w3d(i,j,k)=wa(i,j,k)+dttmp*( wten(i,j,k)                    &
                  -tem1*(pp3d(i,j,k)-pp3d(i,j,k-1))                   &
                       *(c2(1,1,k)*thv(i,j,k)+c1(1,1,k)*thv(i,j,k-1)) )
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
          w3d(i,j,k)=wa(i,j,k)+dttmp*( wten(i,j,k)                    &
                  -tem1*(pp3d(i,j,k)-pp3d(i,j,k-1))*gz(i,j)           &
                       *(c2(i,j,k)*thv(i,j,k)+c1(i,j,k)*thv(i,j,k-1)) )
        enddo
        enddo
        enddo

      ENDIF

        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcw(w3d,1)
        if( terrain_flag ) call bcwsfc(gz,dzdx,dzdy,u3d,v3d,w3d)

!-----

  IF(.not.terrain_flag)THEN

    IF(axisymm.eq.0)THEN
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        div(i,j,k)=(u3d(i+1,j,k)-u3d(i,j,k))*rdx*uh(i)    &
                  +(v3d(i,j+1,k)-v3d(i,j,k))*rdy*vh(j)    &
                  +(w3d(i,j,k+1)-w3d(i,j,k))*rdz*mh(1,1,k)
      enddo
      enddo
      enddo

    ELSE
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        div(i,j,k)=(arh2(i)*u3d(i+1,j,k)-arh1(i)*u3d(i,j,k))*rdx*uh(i)   &
                  +(w3d(i,j,k+1)-w3d(i,j,k))*rdz*mh(1,1,k)
      enddo
      enddo
      enddo

    ENDIF

  ELSE

        ! Cartesian grid with terrain:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        DO k=1,nk
          do j=1,nj
          do i=1,ni+1
            dum1(i,j,k)=u3d(i,j,k)*rgzu(i,j)
          enddo
          enddo
          do j=1,nj+1
          do i=1,ni
            dum2(i,j,k)=v3d(i,j,k)*rgzv(i,j)
          enddo
          enddo
        ENDDO
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
        DO k=1,nk
          IF(k.eq.1)THEN
            do j=1,nj
            do i=1,ni
              dum3(i,j,1)=0.0
              dum3(i,j,nk+1)=0.0
            enddo
            enddo
          ELSE
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
            enddo
            enddo
          ENDIF
        ENDDO
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1)
        do k=1,nk
          do j=1,nj
          do i=1,ni
            div(i,j,k) = gz(i,j)*( (dum1(i+1,j,k)-dum1(i,j,k))*rdx*uh(i)  &
                                  +(dum2(i,j+1,k)-dum2(i,j,k))*rdy*vh(j)  &
                                  +(dum3(i,j,k+1)-dum3(i,j,k))*rdsf(k) )
            if(abs(div(i,j,k)).lt.smeps) div(i,j,k)=0.0
          enddo
          enddo
        enddo

  ENDIF

!-----

      IF( imoist.eq.1 .and. eqtset.eq.2 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          pp3d(i,j,k)=ppi(i,j,k)+dttmp*( ppten(i,j,k)-ppterm(i,j,k)*div(i,j,k) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          th3d(i,j,k)=tha(i,j,k)+dttmp*( thten(i,j,k)-thterm(i,j,k)*div(i,j,k) )
          if(abs(th3d(i,j,k)).lt.smeps) th3d(i,j,k)=0.0
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          pp3d(i,j,k)=ppi(i,j,k)+dttmp*( ppten(i,j,k)-ppterm(i,j,k)*div(i,j,k) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          th3d(i,j,k)=tha(i,j,k)+dttmp*thten(i,j,k)
          if(abs(th3d(i,j,k)).lt.smeps) th3d(i,j,k)=0.0
        enddo
        enddo
        enddo

      ENDIF

        if(timestats.ge.1) time_sound=time_sound+mytime()
 
          if(nrk.lt.3)then
            call bcs(th3d)
            call bcs(pp3d)
          endif

!-------------------------------------------------------------------- 

      return
      end



