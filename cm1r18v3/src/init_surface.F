
      subroutine init_surface(num_soil_layers,dosfcflx,xh,ruh,xf,yh,rvh,yf,   &
           lu_index,xland,tsk,slab_zs,slab_dzs,tslb, &
           emiss,thc,albd,znt,mavail,dsxy,prs0s,prs0,   &
           tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml)
      use module_sf_oml
      implicit none

      include 'input.incl'
      include 'constants.incl'

      integer, intent(in) :: num_soil_layers
      logical, intent(inout) :: dosfcflx
      real, intent(in), dimension(ib:ie) :: xh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh,rvh
      real, intent(in), dimension(jb:je+1) :: yf
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: lu_index
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,znt,xland
      real, intent(inout), dimension(num_soil_layers) :: slab_zs,slab_dzs
      real, intent(inout), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: emiss,thc,albd,mavail,dsxy
      real, intent(in), dimension(ib:ie,jb:je) :: prs0s
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: prs0
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml

      integer :: i,j,k,l
      real :: x1,x2

!-----------------------------------------------------------------------
!  BEGIN user settings

    IF( (sfcmodel.ge.1) .or. (bbc.eq.3) )THEN

      ! Initialize the surface:

      ! Key:
      ! tsk      = "skin temperature" (K) of soil/water  (~1 cm deep)
      ! tmn      = deep-layer temperature (K) of soil  (sfcmodel=2,3 only)
      ! xland    = land/water flag:   1 for land,   2 for water
      ! lu_index = land use index   (1-33) (see LANDUSE.TBL file)

!----------
!  Deafult ... fill in arrays with values from Namelist:

      do j=jb,je
      do i=ib,ie
        tsk(i,j)      = tsk0
      enddo
      enddo

      do j=jbl,jel
      do i=ibl,iel
        tmn(i,j)      = tmn0
        xland(i,j)    = xland0
        lu_index(i,j) = lu0
      enddo
      enddo

!----------
!  Initialize the sea breeze test case from WRF:

      IF( initsfc.eq.2 )THEN 

        ! sea breeze test case from WRF:
        ! (modified on 150830 to match WRF3.7 setup)

        x1 = minx + (3.0/8.0)*(maxx-minx)
        x2 = minx + (5.0/8.0)*(maxx-minx)

        do j=jb,je
        do i=ib,ie
          if( xh(i).ge.x1 .and. xh(i).le.x2 )then
            ! land:
            xland(i,j) = 1.0
            lu_index(i,j) = 18
            tsk(i,j) = 280.0
          else
            ! water:
            xland(i,j) = 2.0
            lu_index(i,j) = 16
            tsk(i,j) = 287.0
          endif
          tmn(i,j) = tsk(i,j)-0.5
        enddo
        enddo

!----------
!  if initsfc is not 1,2:

      ELSEIF( initsfc.ne.1 .and. initsfc.ne.2 )THEN

        ! build your own initial conditions here:

!!!        do j=jb,je
!!!        do i=ib,ie
!!!          tsk(i,j)      = ?
!!!          tmn(i,j)      = ?
!!!          xland(i,j)    = ?
!!!          lu_index(i,j) = ?
!!!        enddo
!!!        enddo

!----------

      ENDIF     ! endif for initsfc options

!     END of user settings
!-----------------------------------------------------------------------
!     NO NOT CHANGE ANYTHING BELOW HERE

    IF( sfcmodel.eq.2 .or. sfcmodel.eq.3 .or. sfcmodel.eq.4 )THEN

      IF (num_soil_layers.NE.1)THEN
      slab_dzs(1)=.01
      slab_zs(1)=.5*slab_dzs(1)
      DO l=2,num_soil_layers
         slab_dzs(l)=2*slab_dzs(l-1)
         slab_zs(l)=slab_zs(l-1)+.5*slab_dzs(l-1)+.5*slab_dzs(l)
      ENDDO
      DO l=1,num_soil_layers
         DO j=jb,je
            DO i=ib,ie
              tslb(i,j,l)=( tsk(i,j)*(slab_zs(num_soil_layers)-slab_zs(l))   &
                          + tmn(i,j)*(slab_zs(l)-slab_zs(1)) ) / &
                          ( slab_zs(num_soil_layers)-slab_zs(1) )
            ENDDO
         ENDDO
      ENDDO
      ENDIF

    ENDIF

    IF( sfcmodel.ge.1 )THEN

      IF( axisymm.eq.1 .or. ny.eq.1 )THEN

        do j=jb,je
        do i=ib,ie
          dsxy(i,j) = dx*ruh(i)
        enddo
        enddo

      ELSEIF( nx.eq.1 )THEN

        do j=jb,je
        do i=ib,ie
          dsxy(i,j) = dy*rvh(j)
        enddo
        enddo

      ELSE

        do j=jb,je
        do i=ib,ie
          dsxy(i,j) = sqrt( (dx*ruh(i))**2 + (dy*rvh(j))**2 )
        enddo
        enddo

      ENDIF

    ENDIF

  ENDIF    ! endif for sfcmodel >= 1

!-----------------------------------------------------------------------
! ensure interoperability of surface schemes and turbulence schemes:
!   NOTE:  DO NOT modify anything here unless you really, really
!          know what you are doing

      ! parameters needed in turb subroutines:
      dosfcflx = .false.

      IF(iturb.ge.1)THEN
        if( isfcflx.eq.1 ) dosfcflx = .true.
        if( sfcmodel.eq.2 .or. sfcmodel.eq.3 .or. sfcmodel.eq.4 )then
          dosfcflx = .true.
        endif
      ENDIF
      IF( ipbl.ge.1 )THEN
        ! sfc flux/drag is handled by pbl subroutine:
        dosfcflx = .false.
      ENDIF

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Settings used in this simulation:'
      if(dowr) write(outfile,*) '  isfcflx  = ',isfcflx
      if(dowr) write(outfile,*) '  bbc      = ',bbc
      if(dowr) write(outfile,*) '  dosfcflx = ',dosfcflx
      if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!  Initialize land-surface parameters:

      if( sfcmodel.ge.1 .or. bbc.eq.3 )then
        ! initialize MM5/WRF landuse categories:
        call getlanduse(season,myid,ib,ie,jb,je,ibl,iel,jbl,jel,   &
                        lu_index,xland,emiss,thc,albd,znt,mavail)
      endif

!-----------------------------------------------------------------------
!  reality checks:

      IF( sfcmodel.ge.1 .or. bbc.eq.3 )THEN

        do j=jb,je
        do i=ib,ie
          IF( nint(xland(i,j)).le.0 .or. nint(xland(i,j)).ge.3 )THEN
            print *
            print *,'  Invalid value for xland '
            print *
            print *,'    xland must be 1 (for land) or 2 (for water) '
            print *
            print *,'  Stopping CM1 .... '
            print *
            call stopcm1
          ENDIF
          if( sfcmodel.eq.2 .or. sfcmodel.eq.3 )then
          IF( lu_index(i,j).le.0 .or. lu_index(i,j).ge.34 )THEN
            print *
            print *,'  Invalid value for lu_index '
            print *
            print *,'    lu_index must be between 1 and 33 '
            print *,'    (see LANDUSE.TBL for more info) '
            print *
            print *,'  Stopping CM1 .... '
            print *
            call stopcm1
          ENDIF
          endif
          IF( tsk(i,j).lt.100.0 )THEN
            print *
            print *,'  WARNING:  tsk is less than 100 K '
            print *
            print *,'  Stopping CM1 .... '
            print *
            call stopcm1
          ENDIF
          if( sfcmodel.eq.2 .or. sfcmodel.eq.3 )then
          IF( tmn(i,j).lt.100.0 )THEN
            print *
            print *,'  WARNING:  tmn is less than 100 K '
            print *
            print *,'  Stopping CM1 .... '
            print *
            call stopcm1
          ENDIF
          endif
        enddo
        enddo

      ENDIF

!-----------------------------------------------------------------------
!  Initialize the ocean mixed layer model:

      if(oceanmodel.eq.2)then
        call omlinit(outfile,oml_hml0, tsk,                     &
                       tml,t0ml,hml,h0ml,huml,hvml,tmoml,       &
                       .false.,         .true. ,                &
                         1, ni+1,   1, nj+1,   1, nk+1,         &
                        ib, ie,  jb, je,  kb, ke,               &
                         1, ni,   1, nj,   1, nk                )
      endif

!-----------------------------------------------------------------------

      end subroutine init_surface

