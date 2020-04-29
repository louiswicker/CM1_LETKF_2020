
!-------------------------------------------------------------
!
!  This subroutine writes data in NetCDF files.
!
!  Code originally written by Daniel Kirshbaum
!  Code converted to netcdf4 (f90) by George Bryan, May 2013
!  Code last modified by George Bryan, 130910
!
!-------------------------------------------------------------


      subroutine netcdf_prelim(rtime,nwrite,ncid,time_index,qname,xh,xf,yh,yf, &
                               xfref,yfref,sigma,sigmaf,zs,zh,zf,        &
                               d2d,ds,du,dv,                             &
                               dumz,dumx,dumy)

      use netcdf

      implicit none
      include 'input.incl'
      include 'constants.incl'

      real, intent(in) :: rtime
      integer, intent(in) :: nwrite
      integer, intent(inout) :: ncid,time_index
      character*3, dimension(maxq), intent(in) :: qname
      real, dimension(ib:ie),   intent(in) :: xh
      real, dimension(ib:ie+1), intent(in) :: xf
      real, dimension(jb:je),   intent(in) :: yh
      real, dimension(jb:je+1), intent(in) :: yf
      real, intent(in), dimension(-2:nx+4) :: xfref
      real, intent(in), dimension(-2:ny+4) :: yfref
      real, dimension(kb:ke)  , intent(in) :: sigma
      real, dimension(kb:ke+1), intent(in) :: sigmaf
      real, dimension(ib:ie,jb:je), intent(in) :: zs
      real, dimension(ib:ie,jb:je,kb:ke),   intent(in) :: zh
      real, dimension(ib:ie,jb:je,kb:ke+1), intent(in) :: zf
      real, dimension(ni,nj) :: d2d
      real, dimension(ni,nj,nk) :: ds
      real, dimension(ni+1,nj,nk) :: du
      real, dimension(ni,nj+1,nk) :: dv
      real, intent(inout), dimension(nz+1) :: dumz
      real, intent(inout), dimension(nx+1) :: dumx
      real, intent(inout), dimension(ny+1) :: dumy


      integer i,j,k,n,irec

      ! Users of GrADS might want to set coards to .true.
      logical, parameter :: coards = .false.

      integer :: cdfid    ! ID for the netCDF file to be created
      integer :: status,dimid,varid
      integer :: niid,njid,nkid,nip1id,njp1id,nkp1id,timeid,oneid,tfile,num_write
      character*8 chid

!-------------------------------------------------------------
! Declare and set integer values for the netCDF dimensions 
!-------------------------------------------------------------

      integer :: ival,jval,kval,ivalp1,jvalp1,kvalp1
      real :: actual_time

      logical :: allinfo

      integer, parameter :: shuffle       = 1
      integer, parameter :: deflate       = 1
      integer, parameter :: deflate_level = 2

      integer :: chkx,chky,chkxp1,chkyp1

      logical, parameter :: stop_on_error = .true.

!--------------------------------------------------------------
! Initializing some things
!--------------------------------------------------------------

    if(coards)then
      if(tapfrq.lt.60.0)then
        print *
        print *,'  Output frequency cannot be less than 60 s for coards format'
        print *
        call stopcm1
      endif
      actual_time = rtime/60.0
    else
      actual_time = rtime
    endif

!--------------------------------------------------------------
!  Write data to cdf file
!--------------------------------------------------------------

    IF(output_filetype.eq.1)THEN
      string(totlen+1:totlen+22) = '.nc                   '
    ELSEIF(output_filetype.eq.2)THEN
      string(totlen+1:totlen+22) = '_XXXXXX.nc            '
      write(string(totlen+2:totlen+7),100) nwrite
    ELSEIF(output_filetype.eq.3)THEN
      string(totlen+1:totlen+17) = '_XXXXXX_YYYYYY.nc     '
      write(string(totlen+ 2:totlen+ 7),100) myid
      write(string(totlen+ 9:totlen+14),100) nwrite
    ELSE
      if(dowr) write(outfile,*) '  for netcdf output, output_filetype must be either 1,2, or 3 '
      call stopcm1
    ENDIF

      if(myid.eq.0) print *,string

100   format(i6.6)

!--------------------------------------------------------------
!  Dimensions of data:

    IF( output_filetype.eq.1 .or. output_filetype.eq.2 )THEN
      ival = nx
      jval = ny
      chkx = nx
      chky = ny
      chkxp1 = nx+1
      chkyp1 = ny+1
    ELSEIF( output_filetype.eq.3 )THEN
      ival = ni
      jval = nj
      chkx = ni
      chky = nj
      chkxp1 = ni+1
      chkyp1 = nj+1
    ELSE
      print *,'  unrecognized value for output_filetype '
      call stopcm1
    ENDIF

    ivalp1 = ival+1
    jvalp1 = jval+1

    kval = min(maxk,nk)
    kvalp1 = min(maxk+1,nk+1)

!--------------------------------------------------------------
!  if this is the start of a file, then do this stuff:

    num_write = nwrite

    allinfo = .false.
    IF(num_write.eq.1) allinfo=.true.
    IF(output_filetype.ge.2)THEN
      allinfo=.true.
      num_write=1
    ENDIF

    IF( num_write.ne.1 )THEN
      ! cm1r18:  Try to open file.
      !          If error, set num_write to 1 and write all info.
      status = nf90_open(string,nf90_write,ncid)
      if( status.eq.nf90_noerr )then
        ! no error, file exists.  Get number of time levels in file:
        call disp_err( nf90_inq_dimid(ncid,'time',timeid) , .true. )
        call disp_err( nf90_inquire_dimension(ncid=ncid,dimid=timeid,len=tfile), .true. )
        if( (tfile+1).lt.num_write )then
          if(myid.eq.0) print *,'  tfile,num_write = ',tfile,num_write
          num_write = tfile+1
        endif
      else
        ! if error opening file, then write all info:
        if(myid.eq.0) print *,'  status = ',status
!!!        if(myid.eq.0) print *,nf90_strerror(status)
        allinfo = .true.
        num_write = 1
      endif
    ENDIF

    time_index = num_write

    ifallinfo: IF(allinfo)THEN
!!!      print *,'  allinfo ... '

!-----------------------------------------------------------------------
!  BEGIN NEW:

      if( myid.eq.0 ) print *,'  calling nf90_create '

!--- works with netcdf 4.2, but not 4.0 (grumble)
!!     call disp_err( nf90_create(string,IOR(nf90_netcdf4, nf90_classic_model),ncid) , .true. )
       call disp_err( nf90_create(string,nf90_netcdf4,ncid) , .true. )





      status = nf90_def_dim(ncid,'ni',ival,niid)
      status = nf90_def_dim(ncid,'nj',jval,njid)
      status = nf90_def_dim(ncid,'nk',kval,nkid)
      status = nf90_def_dim(ncid,'nip1',ivalp1,nip1id)
      status = nf90_def_dim(ncid,'njp1',jvalp1,njp1id)
      status = nf90_def_dim(ncid,'nkp1',kvalp1,nkp1id)
      status = nf90_def_dim(ncid,'time',nf90_unlimited,timeid)
      status = nf90_def_dim(ncid,'one',1,oneid)

    IF(icor.eq.1)THEN
      status = nf90_def_var(ncid,"f_cor",nf90_float,oneid,varid)
      status = nf90_put_att(ncid,varid,"long_name","Coriolis parameter")
      status = nf90_put_att(ncid,varid,"units","1/s")
    ENDIF

    if(.not.coards)then
      status = nf90_def_var(ncid,"ztop",nf90_float,oneid,varid)
      status = nf90_put_att(ncid,varid,"units","km")
    endif

    IF(coards)THEN

      status = nf90_def_var(ncid,"time",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","time since beginning of simulation")
      status = nf90_put_att(ncid,varid,"units","minutes since 2000-07-03 00:00:00")

      status = nf90_def_var(ncid,"ni",nf90_float,niid,varid)
      status = nf90_put_att(ncid,varid,"long_name","west-east location of scalar grid points")
      status = nf90_put_att(ncid,varid,"units","degree_east")

      status = nf90_def_var(ncid,"nip1",nf90_float,nip1id,varid)
      status = nf90_put_att(ncid,varid,"long_name","west-east location of staggered u grid points")
      status = nf90_put_att(ncid,varid,"units","degree_east")

      status = nf90_def_var(ncid,"nj",nf90_float,njid,varid)
      status = nf90_put_att(ncid,varid,"long_name","south-north location of scalar grid points")
      status = nf90_put_att(ncid,varid,"units","degree_north")

      status = nf90_def_var(ncid,"njp1",nf90_float,njp1id,varid)
      status = nf90_put_att(ncid,varid,"long_name","south-north location of staggered v grid points")
      status = nf90_put_att(ncid,varid,"units","degree_north")

      status = nf90_def_var(ncid,"nk",nf90_float,nkid,varid)
      status = nf90_put_att(ncid,varid,"long_name","nominal height of scalar grid points")
      status = nf90_put_att(ncid,varid,"units","km")

      status = nf90_def_var(ncid,"nkp1",nf90_float,nkp1id,varid)
      status = nf90_put_att(ncid,varid,"long_name","nominal height of staggered w grid points")
      status = nf90_put_att(ncid,varid,"units","km")

    ELSE

      status = nf90_def_var(ncid,"time",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","time since beginning of simulation")
      status = nf90_put_att(ncid,varid,"units","seconds since 2003-05-08 21:40:00")

      status = nf90_def_var(ncid,"xh",nf90_float,niid,varid)
      status = nf90_put_att(ncid,varid,"long_name","west-east location of scalar grid points")
      status = nf90_put_att(ncid,varid,"units","km")

      status = nf90_def_var(ncid,"xf",nf90_float,nip1id,varid)
      status = nf90_put_att(ncid,varid,"long_name","west-east location of staggered u grid points")
      status = nf90_put_att(ncid,varid,"units","km")

      status = nf90_def_var(ncid,"yh",nf90_float,njid,varid)
      status = nf90_put_att(ncid,varid,"long_name","south-north location of scalar grid points")
      status = nf90_put_att(ncid,varid,"units","km")

      status = nf90_def_var(ncid,"yf",nf90_float,njp1id,varid)
      status = nf90_put_att(ncid,varid,"long_name","south-north location of staggered v grid points")
      status = nf90_put_att(ncid,varid,"units","km")

      status = nf90_def_var(ncid,"z",nf90_float,nkid,varid)
      status = nf90_put_att(ncid,varid,"long_name","nominal height of scalar grid points")
      status = nf90_put_att(ncid,varid,"units","km")

      status = nf90_def_var(ncid,"zf",nf90_float,nkp1id,varid)
      status = nf90_put_att(ncid,varid,"long_name","nominal height of staggered w grid points")
      status = nf90_put_att(ncid,varid,"units","km")

    ENDIF

!--------------------------------------------------------
!  Just to be sure:

        status = nf90_inq_dimid(ncid,'time',timeid)
        status = nf90_inq_dimid(ncid,'ni',niid)
        status = nf90_inq_dimid(ncid,'nj',njid)
        status = nf90_inq_dimid(ncid,'nk',nkid)
        status = nf90_inq_dimid(ncid,'nip1',nip1id)
        status = nf90_inq_dimid(ncid,'njp1',njp1id)
        status = nf90_inq_dimid(ncid,'nkp1',nkp1id)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--- 2D vars:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(output_rain.eq.1)then
        status = nf90_def_var(ncid,"rain",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","accumulated surface rainfall")
        status = nf90_put_att(ncid,varid,"units","cm")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_sws.eq.1) then
        status = nf90_def_var(ncid,"sws",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","max windspeed at lowest level")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_svs.eq.1) then
        status = nf90_def_var(ncid,"svs",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","max vert vorticity at lowest level")
        status = nf90_put_att(ncid,varid,"units","s-1")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_sps.eq.1) then
        status = nf90_def_var(ncid,"sps",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","min pressure at lowest level")
        status = nf90_put_att(ncid,varid,"units","Pa")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_srs.eq.1) then
        status = nf90_def_var(ncid,"srs",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","max surface rainwater")
        status = nf90_put_att(ncid,varid,"units","kg/kg")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_sgs.eq.1) then
        status = nf90_def_var(ncid,"sgs",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","max surface graupel/hail")
        status = nf90_put_att(ncid,varid,"units","kg/kg")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_sus.eq.1) then
        status = nf90_def_var(ncid,"sus",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","max w at 5 km AGL")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_shs.eq.1) then
        status = nf90_def_var(ncid,"shs",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","max integrated updraft helicity")
        status = nf90_put_att(ncid,varid,"units","m2/s2")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      IF(nrain.eq.2)THEN
        if(output_rain.eq.1)then
          status = nf90_def_var(ncid,"rain2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name","accumulated surface rainfall, translated with moving domain")
          status = nf90_put_att(ncid,varid,"units","cm")

          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        endif
        if(output_sws.eq.1) then
          status = nf90_def_var(ncid,"sws2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name","max windspeed at lowest level, translated with moving domain")
          status = nf90_put_att(ncid,varid,"units","m/s")

          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        endif
        if(output_svs.eq.1) then
          status = nf90_def_var(ncid,"svs2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name","max vorticity at lowest level, translated with moving domain")
          status = nf90_put_att(ncid,varid,"units","s-1")

          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        endif
        if(output_sps.eq.1) then
          status = nf90_def_var(ncid,"sps2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name","min pressure at lowest level, translated with moving domain")
          status = nf90_put_att(ncid,varid,"units","Pa")

          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        endif
        if(output_srs.eq.1) then
          status = nf90_def_var(ncid,"srs2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name","max surface rainwater, translated with moving domain")
          status = nf90_put_att(ncid,varid,"units","kg/kg")

          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        endif
        if(output_sgs.eq.1) then
          status = nf90_def_var(ncid,"sgs2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name","max surface graupel/hail, translated with moving domain")
          status = nf90_put_att(ncid,varid,"units","kg/kg")

          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        endif
        if(output_sus.eq.1) then
          status = nf90_def_var(ncid,"sus2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name","max w at 5 km AGL, translated with moving domain")
          status = nf90_put_att(ncid,varid,"units","m/s")

          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        endif
        if(output_shs.eq.1) then
          status = nf90_def_var(ncid,"shs2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name","translated max integrated updraft helicity")
          status = nf90_put_att(ncid,varid,"units","m2/s2")

          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        endif
      ENDIF

      IF(output_uh.eq.1)THEN
        status = nf90_def_var(ncid,"uh",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","integrated updraft helicity")
        status = nf90_put_att(ncid,varid,"units","m2/s2")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      ENDIF

      IF (output_coldpool.eq.1) THEN
        status = nf90_def_var(ncid,"cpc",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","cold pool intensity C")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        status = nf90_def_var(ncid,"cph",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","cold pool depth h")
        status = nf90_put_att(ncid,varid,"units","m AGL")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      ENDIF

      if (output_sfcflx.eq.1) then
        status = nf90_def_var(ncid,"thflux",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","surface potential temperature flux")
        status = nf90_put_att(ncid,varid,"units","K m s^{-1}")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"qvflux",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","surface water vapor flux")
        status = nf90_put_att(ncid,varid,"units","kg kg^{-1} m s^{-1}")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"tsk",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","soil/ocean temperature")
        status = nf90_put_att(ncid,varid,"units","K")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if(output_sfcparams.eq.1)then
        status = nf90_def_var(ncid,"cd",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","cd")
        status = nf90_put_att(ncid,varid,"units","nondimensional")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"ch",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","ch")
        status = nf90_put_att(ncid,varid,"units","nondimensional")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"cq",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","cq")
        status = nf90_put_att(ncid,varid,"units","nondimensional")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"tlh",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","horizontal turbulence lengthscale for iturb=3")
        status = nf90_put_att(ncid,varid,"units","m")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if (output_psfc.eq.1) then
        status = nf90_def_var(ncid,"psfc",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","surface pressure")
        status = nf90_put_att(ncid,varid,"units","Pa")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if (output_zs.eq.1) then
        status = nf90_def_var(ncid,"zs",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","terrain height")
        status = nf90_put_att(ncid,varid,"units","m")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      IF (output_dbz.eq.1) THEN
        status = nf90_def_var(ncid,"cref",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","composite reflectivity (dBZ)")
        status = nf90_put_att(ncid,varid,"units","dBZ")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      ENDIF

      if(output_sfcparams.eq.1)then
        status = nf90_def_var(ncid,"xland",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","land/water flag (1=land,2=water)")
        status = nf90_put_att(ncid,varid,"units","integer flag")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"lu",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","land use index")
        status = nf90_put_att(ncid,varid,"units","integer flag")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"mavail",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","surface moisture availability")
        status = nf90_put_att(ncid,varid,"units","integer flag")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.sfcmodel.eq.3.or.oceanmodel.eq.2))then
        status = nf90_def_var(ncid,"tmn",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","deep-layer soil temperature")
        status = nf90_put_att(ncid,varid,"units","K")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"hfx",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","heat flux at surface")
        status = nf90_put_att(ncid,varid,"units","W/m^2")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"qfx",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","surface moisture flux")
        status = nf90_put_att(ncid,varid,"units","W/m^2")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"gsw",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","downward SW flux at surface")
        status = nf90_put_att(ncid,varid,"units","W/m^2")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"glw",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","downward LW flux at surface")
        status = nf90_put_att(ncid,varid,"units","W/m^2")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.sfcmodel.eq.3))then
        status = nf90_def_var(ncid,"tslb1",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","soil temp, layer 1")
        status = nf90_put_att(ncid,varid,"units","K")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"tslb2",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","soil temp, layer 2")
        status = nf90_put_att(ncid,varid,"units","K")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"tslb3",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","soil temp, layer 3")
        status = nf90_put_att(ncid,varid,"units","K")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"tslb4",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","soil temp, layer 4")
        status = nf90_put_att(ncid,varid,"units","K")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"tslb5",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","soil temp, layer 5")
        status = nf90_put_att(ncid,varid,"units","K")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if(output_sfcparams.eq.1.and.oceanmodel.eq.2)then
        status = nf90_def_var(ncid,"tml",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","ocean mixed layer temperature")
        status = nf90_put_att(ncid,varid,"units","K")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"hml",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","ocean mixed layer depth")
        status = nf90_put_att(ncid,varid,"units","m")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"huml",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","ocean mixed layer u vel.")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )


        status = nf90_def_var(ncid,"hvml",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","ocean mixed layer v vel.")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if(output_radten.eq.1)then
        status = nf90_def_var(ncid,"radsw",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","solar radiation at surface")
        status = nf90_put_att(ncid,varid,"units","w/m^2")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"rnflx",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","net radiation absorbed by surface")
        status = nf90_put_att(ncid,varid,"units","W/m^2")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"radswnet",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","net solar radiation")
        status = nf90_put_att(ncid,varid,"units","W/m^2")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"radlwin",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","incoming longwave radiation")
        status = nf90_put_att(ncid,varid,"units","W/m^2")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        ! MS add
        status = nf90_def_var(ncid,"olr",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","TOA net outgoing longwave radiation")
        status = nf90_put_att(ncid,varid,"units","W/m^2")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        status = nf90_def_var(ncid,"dsr",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","TOA net incoming solar radiation")
        status = nf90_put_att(ncid,varid,"units","W/m^2")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      IF(output_sfcdiags.eq.1)THEN
        status = nf90_def_var(ncid,"u10",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","diagnostic 10 m u wind")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"v10",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","diagnostic 10 m v wind")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"t2",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","diagnostic 2 m temperature")
        status = nf90_put_att(ncid,varid,"units","K")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"q2",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","diagnostic 2 m mixing ratio")
        status = nf90_put_att(ncid,varid,"units","kg/kg")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"znt",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","roughness length")
        status = nf90_put_att(ncid,varid,"units","m")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"ust",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","friction velocity")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"hpbl",nf90_float,(/niid,njid,timeid/),varid)
      if(ipbl.eq.1)then
        status = nf90_put_att(ncid,varid,"long_name","PBL height (from PBL scheme)")
      else
        status = nf90_put_att(ncid,varid,"long_name","rough estimate of PBL height")
      endif
        status = nf90_put_att(ncid,varid,"units","m")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"zol",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","z/L (z over Monin-Obukhov length)")
        status = nf90_put_att(ncid,varid,"units","   ")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"mol",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","T* (similarity theory)")
        status = nf90_put_att(ncid,varid,"units","K")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )



        status = nf90_def_var(ncid,"br",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","bulk Richardson number in surface layer")
        status = nf90_put_att(ncid,varid,"units","dimensionless")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"psim",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","similarity stability function for momentum")
        status = nf90_put_att(ncid,varid,"units","dimensionless")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"psih",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","similarity stability function for heat")
        status = nf90_put_att(ncid,varid,"units","dimensionless")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--- 3D vars:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (output_zh.eq.1) then
        status = nf90_def_var(ncid,"zh",nf90_float,(/niid,njid,nkid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","height (above sea level) of scalar grid points")
        status = nf90_put_att(ncid,varid,"units","m")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if(output_th.eq.1)then
        status = nf90_def_var(ncid,"th",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","potential temperature")
        status = nf90_put_att(ncid,varid,"units","K")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_thpert.eq.1)then
        status = nf90_def_var(ncid,"thpert",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","perturbation potential temperature")
        status = nf90_put_att(ncid,varid,"units","K")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if(output_prs.eq.1)then
        status = nf90_def_var(ncid,"prs",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","pressure")
        status = nf90_put_att(ncid,varid,"units","Pa")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_prspert.eq.1)then
        status = nf90_def_var(ncid,"prspert",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","perturbation pressure")
        status = nf90_put_att(ncid,varid,"units","Pa")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if(output_pi.eq.1)then
        status = nf90_def_var(ncid,"pi",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","nondimensional pressure")
        status = nf90_put_att(ncid,varid,"units","dimensionless")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_pipert.eq.1)then
        status = nf90_def_var(ncid,"pipert",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","perturbation nondimensional pressure")
        status = nf90_put_att(ncid,varid,"units","dimensionless")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_rho.eq.1)then
        status = nf90_def_var(ncid,"rho",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","density of dry air")
        status = nf90_put_att(ncid,varid,"units","kg/m^3")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_rhopert.eq.1)then
        status = nf90_def_var(ncid,"rhopert",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","perturbation density of dry air")
        status = nf90_put_att(ncid,varid,"units","kg/m^3")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
    IF(iptra.eq.1)THEN
      do n=1,npt
        chid = 'pt      '
        if( n.le.9 )then
          write(chid(3:3),111) n
111       format(i1)
        else
          write(chid(3:4),141) n
141       format(i2)
        endif
        status = nf90_def_var(ncid,chid,nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","mixing ratio of passive tracer")
        status = nf90_put_att(ncid,varid,"units","kg/kg")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      enddo
    ENDIF
    IF(imoist.eq.1)THEN
      if(output_qv.eq.1)then
        status = nf90_def_var(ncid,"qv",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","water vapor mixing ratio")
        status = nf90_put_att(ncid,varid,"units","kg/kg")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_qvpert.eq.1)then
        status = nf90_def_var(ncid,"qvpert",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","perturbation water vapor mixing ratio")
        status = nf90_put_att(ncid,varid,"units","kg/kg")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
        if(output_q.eq.1)then
        do n=1,numq
          if(n.ne.nqv)then
            status = nf90_def_var(ncid,qname(n),nf90_float,(/niid,njid,nkid,timeid/),varid)
            if(idm.eq.1.and.n.ge.nnc1.and.n.le.nnc2)then
              status = nf90_put_att(ncid,varid,"long_name","number concentration")
              status = nf90_put_att(ncid,varid,"units","kg^{-1}")
            elseif(idm .eq. 1 .and. n.ge.nzl1 .and. n .le. nzl2)then
              status = nf90_put_att(ncid,varid,"long_name","reflectivity moment")
              status = nf90_put_att(ncid,varid,"units","Z m^{-3}kg^{-1}")
            elseif(idm .eq. 1 .and. n.ge.nvl1 .and. n .le. nvl2)then
              status = nf90_put_att(ncid,varid,"long_name","particle volume")
              status = nf90_put_att(ncid,varid,"units","m^{3}kg^{-1}")
            else
              status = nf90_put_att(ncid,varid,"long_name","mixing ratio")
              status = nf90_put_att(ncid,varid,"units","kg/kg")
            endif

            call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
            call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

          endif
        enddo
      endif
      if(output_dbz.eq.1)then
        status = nf90_def_var(ncid,"dbz",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","reflectivity")
        status = nf90_put_att(ncid,varid,"units","dBZ")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
    ENDIF

      if(output_buoyancy.eq.1)then
        status = nf90_def_var(ncid,"buoyancy",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","buoyancy")
        status = nf90_put_att(ncid,varid,"units","m s^-2")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if(output_uinterp.eq.1)then
        status = nf90_def_var(ncid,"uinterp",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","velocity in x-direction, interpolated to scalar points")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_vinterp.eq.1)then
        status = nf90_def_var(ncid,"vinterp",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","velocity in y-direction, interpolated to scalar points")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_winterp.eq.1)then
        status = nf90_def_var(ncid,"winterp",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","velocity in z-direction, interpolated to scalar points")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if(output_vort.eq.1)then
        status = nf90_def_var(ncid,"xvort",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","horizontal vorticity (x)")
        status = nf90_put_att(ncid,varid,"units","s^-1")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"yvort",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","horizontal vorticity (y)")
        status = nf90_put_att(ncid,varid,"units","s^-1")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"zvort",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","vertical vorticity")
        status = nf90_put_att(ncid,varid,"units","s^-1")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if(output_pv.eq.1)then
        status = nf90_def_var(ncid,"pv",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","potential vorticity")
        status = nf90_put_att(ncid,varid,"units","K m^2 kg^-1 s^-1")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if (output_basestate.eq.1) then

        status = nf90_def_var(ncid,"pi0",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","base-state nondimensional pressure")
        status = nf90_put_att(ncid,varid,"units","dimensionless")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"th0",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","base-state potential temperature")
        status = nf90_put_att(ncid,varid,"units","K")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"prs0",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","base-state pressure")
        status = nf90_put_att(ncid,varid,"units","Pa")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"qv0",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","base-state water vapor mixing ratio")
        status = nf90_put_att(ncid,varid,"units","kg/kg")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      endif

      if(output_pblten.eq.1)then
        status = nf90_def_var(ncid,"thpten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","pbl tendency: theta")
        status = nf90_put_att(ncid,varid,"units","    ")

        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )


        status = nf90_def_var(ncid,"qvpten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","pbl tendency: qv")
        status = nf90_put_att(ncid,varid,"units","    ")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"qcpten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","pbl tendency: qc")
        status = nf90_put_att(ncid,varid,"units","    ")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"qipten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","pbl tendency: qi")
        status = nf90_put_att(ncid,varid,"units","    ")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"upten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","pbl tendency: u")
        status = nf90_put_att(ncid,varid,"units","    ")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"vpten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","pbl tendency: v")
        status = nf90_put_att(ncid,varid,"units","    ")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if(output_radten.eq.1)then
        status = nf90_def_var(ncid,"swten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","pot temp tendency, sw rad")
        status = nf90_put_att(ncid,varid,"units","    ")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


        status = nf90_def_var(ncid,"lwten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","pot temp tendency, lw rad")
        status = nf90_put_att(ncid,varid,"units","K/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif


      IF( output_turbten.eq.1 )THEN
        status = nf90_def_var(ncid,"ftt",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","theta tendency: turbulence scheme")
        status = nf90_put_att(ncid,varid,"units","K/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        IF(imoist.eq.1)THEN
          status = nf90_def_var(ncid,"ftq",nf90_float,(/niid,njid,nkid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name","qv tendency: turbulence scheme")
          status = nf90_put_att(ncid,varid,"units","kg/kg/s")

          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        ENDIF
      ENDIF

      IF( output_dissheat.eq.1 )THEN
        status = nf90_def_var(ncid,"dissheat",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","dissipative heating (potential temperature tendency)")
        status = nf90_put_att(ncid,varid,"units","K/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      ENDIF


      IF( output_mptend.eq.1 )THEN
        status = nf90_def_var(ncid,"mptend",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","potential temperature tendency from microphysics scheme")
        status = nf90_put_att(ncid,varid,"units","K/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      ENDIF

      IF( output_fallvel.eq.1 )THEN
        if( qd_vtc.gt.0 )then
          status = nf90_def_var(ncid,"vtc",nf90_float,(/niid,njid,nkid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qc")
          status = nf90_put_att(ncid,varid,"units","m/s")

          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        endif
        if( qd_vtr.gt.0 )then
          status = nf90_def_var(ncid,"vtr",nf90_float,(/niid,njid,nkid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qr")
          status = nf90_put_att(ncid,varid,"units","m/s")

          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        endif
        if( qd_vts.gt.0 )then
          status = nf90_def_var(ncid,"vts",nf90_float,(/niid,njid,nkid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qs")
          status = nf90_put_att(ncid,varid,"units","m/s")

          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        endif
        if( qd_vtg.gt.0 )then
          status = nf90_def_var(ncid,"vtg",nf90_float,(/niid,njid,nkid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qg")
          status = nf90_put_att(ncid,varid,"units","m/s")

          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        endif
        if( qd_vti.gt.0 )then
          status = nf90_def_var(ncid,"vti",nf90_float,(/niid,njid,nkid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qi")
          status = nf90_put_att(ncid,varid,"units","m/s")

          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        endif
      ENDIF

      if(output_u.eq.1)then
        status = nf90_def_var(ncid,"u",nf90_float,(/nip1id,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","velocity in x-direction")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkxp1,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_upert.eq.1)then
        status = nf90_def_var(ncid,"upert",nf90_float,(/nip1id,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","perturbation velocity in x-direction")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkxp1,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if (output_basestate.eq.1) then
        status = nf90_def_var(ncid,"u0",nf90_float,(/nip1id,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","base-state x-component of velocity")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkxp1,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_turbten.eq.1)then
        status = nf90_def_var(ncid,"ftu",nf90_float,(/nip1id,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","u tendency: turbulence scheme")
        status = nf90_put_att(ncid,varid,"units","m/s/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkxp1,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_impdiften.eq.1)then
        status = nf90_def_var(ncid,"fdu",nf90_float,(/nip1id,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","u tendency: implicit diffusion")
        status = nf90_put_att(ncid,varid,"units","m/s/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkxp1,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if(output_v.eq.1)then
        status = nf90_def_var(ncid,"v",nf90_float,(/niid,njp1id,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","velocity in y-direction")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chkyp1,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_vpert.eq.1)then
        status = nf90_def_var(ncid,"vpert",nf90_float,(/niid,njp1id,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","perturbation velocity in y-direction")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chkyp1,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if (output_basestate.eq.1) then
        status = nf90_def_var(ncid,"v0",nf90_float,(/niid,njp1id,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","base-state y-component of velocity")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chkyp1,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_turbten.eq.1)then
        status = nf90_def_var(ncid,"ftv",nf90_float,(/niid,njp1id,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","v tendency: turbulence scheme")
        status = nf90_put_att(ncid,varid,"units","m/s/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chkyp1,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_impdiften.eq.1)then
        status = nf90_def_var(ncid,"fdv",nf90_float,(/niid,njp1id,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","v tendency: implicit diffusion")
        status = nf90_put_att(ncid,varid,"units","m/s/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chkyp1,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

      if(output_w.eq.1)then
        status = nf90_def_var(ncid,"w",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","velocity in z-direction")
        status = nf90_put_att(ncid,varid,"units","m/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      IF((iturb.eq.1).and.(output_tke.eq.1))THEN
        status = nf90_def_var(ncid,"tke",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","subgrid turbulence kinetic energy")
        status = nf90_put_att(ncid,varid,"units","m^2/s^2")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      ENDIF
      IF(output_km.eq.1)THEN
        !----
        status = nf90_def_var(ncid,"kmh",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","eddy mixing coefficient for momentum in the horizontal direction")
        status = nf90_put_att(ncid,varid,"units","m^2/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        !----
        status = nf90_def_var(ncid,"kmv",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","eddy mixing coefficient for momentum in the vertical direction")
        status = nf90_put_att(ncid,varid,"units","m^2/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        !----
      ENDIF
      IF(output_kh.eq.1)THEN
        !----
        status = nf90_def_var(ncid,"khh",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","eddy mixing coefficient for scalars in the horizontal direction")
        status = nf90_put_att(ncid,varid,"units","m^2/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        !----
        status = nf90_def_var(ncid,"khv",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","eddy mixing coefficient for scalars in the vertical direction")
        status = nf90_put_att(ncid,varid,"units","m^2/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        !----
      ENDIF

      if(output_dissten.eq.1)then
        status = nf90_def_var(ncid,"dissten",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","dissipation rate")
        status = nf90_put_att(ncid,varid,"units","m^2/s^3")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_nm.eq.1)then
        status = nf90_def_var(ncid,"nm",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","squared Brunt-Vaisala freq")
        status = nf90_put_att(ncid,varid,"units","s^{-2}")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_def.eq.1)then
        status = nf90_def_var(ncid,"defv",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","vertical deformation")
        status = nf90_put_att(ncid,varid,"units","s^{-2}")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

        status = nf90_def_var(ncid,"defh",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","horizontal deformation")
        status = nf90_put_att(ncid,varid,"units","s^{-2}")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_turbten.eq.1)then
        status = nf90_def_var(ncid,"ftw",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","w tendency: turbulence scheme")
        status = nf90_put_att(ncid,varid,"units","m/s/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if(output_impdiften.eq.1)then
        status = nf90_def_var(ncid,"fdw",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","w tendency: implicit diffusion")
        status = nf90_put_att(ncid,varid,"units","m/s/s")

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif

!--------------------------------------------------

    if(coards)then
      status = nf90_put_att(ncid,NF90_GLOBAL,'Conventions','COARDS')
    endif
      status = nf90_put_att(ncid,NF90_GLOBAL,'x_units','km')
      status = nf90_put_att(ncid,NF90_GLOBAL,'x_label','x')
      status = nf90_put_att(ncid,NF90_GLOBAL,'y_units','km')
      status = nf90_put_att(ncid,NF90_GLOBAL,'y_label','y')
      status = nf90_put_att(ncid,NF90_GLOBAL,'z_units','km')
      status = nf90_put_att(ncid,NF90_GLOBAL,'z_label','z')

      status = nf90_enddef(ncid)

! ... end of defs
!--------------------------------------------------
! begin data ... initial time ...

    IF(icor.eq.1)THEN
      status = nf90_inq_varid(ncid,'f_cor',varid)
      status = nf90_put_var(ncid,varid,fcor)
    ENDIF

    if(.not.coards)then
      status = nf90_inq_varid(ncid,'ztop',varid)
      status = nf90_put_var(ncid,varid,0.001*ztop)
    endif

      if(coards)then
        status = nf90_inq_varid(ncid,'ni',varid)
      else
        status = nf90_inq_varid(ncid,'xh',varid)
      endif
      do i=1,nx
        dumx(i) = 0.001*0.5*(xfref(i)+xfref(i+1))
      enddo
      status = nf90_put_var(ncid,varid,dumx,(/1/),(/nx/))

      if(coards)then
        status = nf90_inq_varid(ncid,'nip1',varid)
      else
        status = nf90_inq_varid(ncid,'xf',varid)
      endif
      do i=1,nx+1
        dumx(i) = 0.001*xfref(i)
      enddo
      status = nf90_put_var(ncid,varid,dumx,(/1/),(/nx+1/))

      if(coards)then
        status = nf90_inq_varid(ncid,'nj',varid)
      else
        status = nf90_inq_varid(ncid,'yh',varid)
      endif
      do j=1,ny
        dumy(j) = 0.001*0.5*(yfref(j)+yfref(j+1))
      enddo
      status = nf90_put_var(ncid,varid,dumy,(/1/),(/ny/))

      if(coards)then
        status = nf90_inq_varid(ncid,'njp1',varid)
      else
        status = nf90_inq_varid(ncid,'yf',varid)
      endif
      do j=1,ny+1
        dumy(j) = 0.001*yfref(j)
      enddo
      status = nf90_put_var(ncid,varid,dumy,(/1/),(/ny+1/))

      if(coards)then
        status = nf90_inq_varid(ncid,'nk',varid)
      else
        status = nf90_inq_varid(ncid,'z',varid)
      endif
      if(terrain_flag)then
        do k=1,kval
          dumz(k) = 0.001*sigma(k)
        enddo
        status = nf90_put_var(ncid,varid,dumz,(/1/),(/kval/))
      else
        do k=1,kval
          dumz(k) = 0.001*zh(1,1,k)
        enddo
        status = nf90_put_var(ncid,varid,dumz,(/1/),(/kval/))
      endif

      if(coards)then
        status = nf90_inq_varid(ncid,'nkp1',varid)
      else
        status = nf90_inq_varid(ncid,'zf',varid)
      endif
      if(terrain_flag)then
        do k=1,kvalp1
          dumz(k) = 0.001*sigmaf(k)
        enddo
        status = nf90_put_var(ncid,varid,dumz,(/1/),(/kvalp1/))
      else
        do k=1,kvalp1
          dumz(k) = 0.001*zf(1,1,k)
        enddo
        status = nf90_put_var(ncid,varid,dumz,(/1/),(/kvalp1/))
      endif

      ! ... end if info at initial time only

!----------------------------------------------------------

    ENDIF ifallinfo

      call disp_err( nf90_inq_varid(ncid,'time',varid) , stop_on_error )
      call disp_err( nf90_put_var(ncid,varid,actual_time,(/time_index/)) , stop_on_error )



    return
    end subroutine netcdf_prelim


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



    ! cm1r18:  netcdf restart files
    subroutine restart_prelim(nrst,ncid,mtime,xfref,yfref,zh,zf,sigma,sigmaf,  &
                              qname,num_soil_layers,nrad2d,dumx,dumy,dumz)

    use netcdf
    use module_restart
    implicit none

    include 'input.incl'
    include 'constants.incl'

    integer, intent(in) :: nrst
    integer, intent(inout) :: ncid
    double precision, intent(in) :: mtime
    real, intent(in), dimension(-2:nx+4) :: xfref
    real, intent(in), dimension(-2:ny+4) :: yfref
    real, intent(in), dimension(ib:ie,jb:je,kb:ke)   :: zh
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
    real, intent(in), dimension(kb:ke)   :: sigma
    real, intent(in), dimension(kb:ke+1) :: sigmaf
    character*3, intent(in), dimension(maxq) :: qname
    integer, intent(in) :: num_soil_layers,nrad2d
    real, intent(inout), dimension(nx+1) :: dumx
    real, intent(inout), dimension(ny+1) :: dumy
    real, intent(inout), dimension(nz+1) :: dumz

    integer :: i,j,k,n
    integer :: ival,jval,kval,ivalp1,jvalp1,kvalp1
    integer :: chkx,chky,chkxp1,chkyp1
    integer :: time_index,status,varid
    real :: actual_time
    integer :: niid,njid,nkid,nip1id,njp1id,nkp1id,timeid
    integer :: nbudgetid,numqid,plocid,nparcelsid
    character*8 :: text1

    integer, parameter :: shuffle       = 1
    integer, parameter :: deflate       = 1
    integer, parameter :: deflate_level = 2

    chkx = nx
    chky = ny

    chkxp1 = nx+1
    chkyp1 = ny+1

    time_index = 1

    actual_time = mtime

    string(totlen+1:totlen+22) = '_rst_XXXXXX.nc        '
    write(string(totlen+6:totlen+11),100) nrst
100 format(i6.6)

    if(myid.eq.0) print *,'  string = ',string

!--------------------------------------------------------------
!  Dimensions of data:

    ival = nx
    jval = ny
    kval = nk

    ivalp1 = ival+1
    jvalp1 = jval+1
    kvalp1 = kval+1

!--------------------------------------------------------------
!  if this is the start of a file, then do this stuff:



!--- works with netcdf 4.2, but not 4.0 (grumble)
      call disp_err( nf90_create(string,IOR(nf90_netcdf4, nf90_classic_model),ncid) , .true. )
      IF(myid.eq.0) print *,' USING COMPRESSED netCDF4 '






!-----------------------------------------------------------------------
!  BEGIN

      !------------------
      ! define dims:

      call disp_err( nf90_def_dim(ncid,'ni',ival,niid) , .true. )
      call disp_err( nf90_def_dim(ncid,'nj',jval,njid) , .true. )
      call disp_err( nf90_def_dim(ncid,'nk',kval,nkid) , .true. )
      call disp_err( nf90_def_dim(ncid,'nip1',ivalp1,nip1id) , .true. )
      call disp_err( nf90_def_dim(ncid,'njp1',jvalp1,njp1id) , .true. )
      call disp_err( nf90_def_dim(ncid,'nkp1',kvalp1,nkp1id) , .true. )
      call disp_err( nf90_def_dim(ncid,'time',1,timeid) , .true. )
      call disp_err( nf90_def_dim(ncid,'nbudget',nbudget,nbudgetid) , .true. )
      call disp_err( nf90_def_dim(ncid,'numq',numq,numqid) , .true. )
    if( iprcl.eq.1 )then
      n = 3
      call disp_err( nf90_def_dim(ncid,'nploc',n,plocid) , .true. )
      n = max( 1 , nparcels )
      call disp_err( nf90_def_dim(ncid,'nparcels',n,nparcelsid) , .true. )
    endif

      !------------------
      ! define vars:

      call disp_err( nf90_def_var(ncid,"time",nf90_float,timeid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","time since beginning of simulation") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","seconds since 2003-05-08 21:40:00") , .true. )

      call disp_err( nf90_def_var(ncid,"xh",nf90_float,niid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","west-east location of scalar grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

      call disp_err( nf90_def_var(ncid,"xf",nf90_float,nip1id,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","west-east location of staggered u grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

      call disp_err( nf90_def_var(ncid,"yh",nf90_float,njid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","south-north location of scalar grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

      call disp_err( nf90_def_var(ncid,"yf",nf90_float,njp1id,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","south-north location of staggered v grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

      call disp_err( nf90_def_var(ncid,"zh",nf90_float,nkid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","nominal height of scalar grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

      call disp_err( nf90_def_var(ncid,"zf",nf90_float,nkp1id,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","nominal height of staggered w grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

  !------------------
  ! vars:

      call disp_err( nf90_def_var(ncid,"nstep"   ,nf90_int,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"nrec"    ,nf90_int,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"prec"    ,nf90_int,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"nwrite"  ,nf90_int,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"nrst"    ,nf90_int,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"ndt"     ,nf90_int,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"old_format",nf90_int,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"npt"     ,nf90_int,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"numparcels",nf90_int,varid) , .true. )

      call disp_err( nf90_def_var(ncid,"dt"      ,nf90_float,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","timestep") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","s") , .true. )

      call disp_err( nf90_def_var(ncid,"dtlast"  ,nf90_float,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","previous timestep") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","s") , .true. )

      call disp_err( nf90_def_var(ncid,"cflmax"  ,nf90_float,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max Courant number from previous timestep") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )

      call disp_err( nf90_def_var(ncid,"mtime"   ,nf90_double,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","model time (i.e., time since beginning of simulation)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","s") , .true. )

      call disp_err( nf90_def_var(ncid,"stattim" ,nf90_double,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"taptim"  ,nf90_double,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"rsttim"  ,nf90_double,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"radtim"  ,nf90_double,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"prcltim" ,nf90_double,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"adt"     ,nf90_double,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"acfl"    ,nf90_double,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"dbldt"   ,nf90_double,varid) , .true. )
      call disp_err( nf90_def_var(ncid,"mass1"   ,nf90_double,varid) , .true. )

      call disp_err( nf90_def_var(ncid,"qbudget" ,nf90_double,(/nbudgetid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"asq"     ,nf90_double,(/numqid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"bsq"     ,nf90_double,(/numqid/),varid) , .true. )

  !------------------
  ! 2d/3d vars:

      call disp_err( nf90_def_var(ncid,"rain"    ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","accumulated surface rainfall") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","cm") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"sws"     ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max windspeed at lowest level") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"svs"     ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max vert vorticity at lowest level") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","s^(-1)") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"sps"     ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","min pressure at lowest level") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","Pa") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"srs"     ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max surface rainwater") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"sgs"     ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max surface graupel/hail") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"sus"     ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max w at 5 km AGL") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"shs"     ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max integrated updraft helicity") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m2/s2") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


    if( nrain.eq.2 )then

      call disp_err( nf90_def_var(ncid,"rain2"   ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","accumulated surface rainfall (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","cm") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"sws2"    ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max windspeed at lowest level (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"svs2"    ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max vert vorticity at lowest level (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","s^(-1)") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"sps2"    ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","min pressure at lowest level (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","Pa") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"srs2"    ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max surface rainwater (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"sgs2"    ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max surface graupel/hail (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"sus2"    ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max w at 5 km AGL (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"shs2"    ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max integrated updraft helicity (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m2/s2") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


    endif

      call disp_err( nf90_def_var(ncid,"tsk"     ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","soil/ocean temperature") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","K") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


! 3D VARS

      call disp_err( nf90_def_var(ncid,"rho"     ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","dry-air density") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/m^3") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"prs"     ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","pressure") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","Pa") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"ua"      ,nf90_float,(/nip1id,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","west-east velocity (at u points)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkxp1,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"va"      ,nf90_float,(/niid,njp1id,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","south-north velocity (at v points)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chkyp1,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"wa"      ,nf90_float,(/niid,njid,nkp1id/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","vertical velocity (at w points)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"ppi"     ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","perturbation non-dimensional pressure") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"tha"     ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","perturbation potential temperature") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","K") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"ppx"     ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","change in nondimensional pressure used for forward-time-weighting on small steps") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


    IF( imoist.eq.1 )THEN
      do n=1,numq
        text1 = '        '
        write(text1(1:3),156) qname(n)
156     format(a3)
        call disp_err( nf90_def_var(ncid,text1     ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
        if( n.eq.nqv )then
          call disp_err( nf90_put_att(ncid,varid,"long_name","water vapor mixing ratio") , .true. )
          call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )
        elseif(idm.eq.1.and.n.ge.nnc1.and.n.le.nnc2)then
          call disp_err( nf90_put_att(ncid,varid,"long_name","number concentration") , .true. )
          call disp_err( nf90_put_att(ncid,varid,"units","kg^{-1}") , .true. )
        elseif(idm .eq. 1 .and. n.ge.nzl1 .and. n .le. nzl2)then
          call disp_err( nf90_put_att(ncid,varid,"long_name","reflectivity moment") , .true. )
          call disp_err( nf90_put_att(ncid,varid,"units","Z m^{-3}kg^{-1}") , .true. )
        elseif(idm .eq. 1 .and. n.ge.nvl1 .and. n .le. nvl2)then
          call disp_err( nf90_put_att(ncid,varid,"long_name","particle volume") , .true. )
          call disp_err( nf90_put_att(ncid,varid,"units","m^{3}kg^{-1}") , .true. )
        else
          call disp_err( nf90_put_att(ncid,varid,"long_name","mixing ratio") , .true. )
          call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )
        endif

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      enddo
    ENDIF
    if(imoist.eq.1.and.eqtset.eq.2)then
      call disp_err( nf90_def_var(ncid,"qpten"   ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","ppi tendency from microphysics on previous timestep (related to h_diabatic in WRF)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","s^(-1)") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      !---
      call disp_err( nf90_def_var(ncid,"qtten"   ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","theta tendency from microphysics on previous timestep (related to h_diabatic in WRF)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","K/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      !---
      call disp_err( nf90_def_var(ncid,"qvten"   ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","qv tendency from microphysics on previous timestep (related to qv_diabatic in WRF)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/kg/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      !---
      call disp_err( nf90_def_var(ncid,"qcten"   ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","qc tendency from microphysics on previous timestep (related to qc_diabatic in WRF)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/kg/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

    endif
    if(iturb.eq.1)then
      call disp_err( nf90_def_var(ncid,"tkea"    ,nf90_float,(/niid,njid,nkp1id/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","subgrid turbulence kinetic energy") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m^2/s^2") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

    endif

    if( radopt.eq.1 )then

      call disp_err( nf90_def_var(ncid,"lwten"   ,nf90_float,(/niid,njid,nkid/),varid) , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"swten"   ,nf90_float,(/niid,njid,nkid/),varid) , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"radsw"   ,nf90_float,(/niid,njid/),varid) , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"rnflx"   ,nf90_float,(/niid,njid/),varid) , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )

      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
      call disp_err( nf90_def_var(ncid,"radswnet",nf90_float,(/niid,njid/),varid) , .true. )


      call disp_err( nf90_def_var(ncid,"radlwin" ,nf90_float,(/niid,njid/),varid) , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      do n=1,nrad2d
        if( n.lt.10 )then
          text1 = 'radX    '
          write(text1(4:4),181) n
181       format(i1.1)
        elseif( n.lt.100 )then
          text1 = 'radXX   '
          write(text1(4:5),182) n
182       format(i2.2)
        elseif( n.lt.1000 )then
          text1 = 'radXXX  '
          write(text1(4:6),183) n
183       format(i3.3)
        else
          stop 11611
        endif
        call disp_err( nf90_def_var(ncid,text1     ,nf90_float,(/niid,njid/),varid) , .true. )
      enddo
     endif

    if( radopt.ge.1 .and. ptype.eq.5 )then

      call disp_err( nf90_def_var(ncid,"effc"    ,nf90_float,(/niid,njid,nkid/),varid) , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"effi"    ,nf90_float,(/niid,njid,nkid/),varid) , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"effs"    ,nf90_float,(/niid,njid,nkid/),varid) , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"effr"    ,nf90_float,(/niid,njid,nkid/),varid) , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"effg"    ,nf90_float,(/niid,njid,nkid/),varid) , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


      call disp_err( nf90_def_var(ncid,"effis"   ,nf90_float,(/niid,njid,nkid/),varid) , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )


    endif

    if((oceanmodel.eq.2).or.(ipbl.eq.1).or.(sfcmodel.ge.1))then
      if(sfcmodel.ge.1)then
        call disp_err( nf90_def_var(ncid,"ust"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","friction velocity") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"znt"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","roughness length") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

        call disp_err( nf90_def_var(ncid,"cd"      ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","surface drag coefficient") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )

        call disp_err( nf90_def_var(ncid,"ch"      ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","surface exchange coefficient for sensible heat") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )

        call disp_err( nf90_def_var(ncid,"cq"      ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","surface exchange coefficient for moisture") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )

        call disp_err( nf90_def_var(ncid,"u1"      ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","u component of velocity at lowest model level") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"v1"      ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","v component of velocity at lowest model level") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"s1"      ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","horizontal windspeed at lowest model level") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"u10"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","u component of windspeed at z=10 m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"v10"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","v component of windspeed at z=10 m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"s10"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","horizontal windspeed at z=10 m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"xland"   ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","land/water flag (1=land,2=water)") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","integer flag") , .true. )

        call disp_err( nf90_def_var(ncid,"thflux"  ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","surface potential temperature flux") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","K m s^{-1}") , .true. )

        call disp_err( nf90_def_var(ncid,"qvflux"  ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","surface water vapor flux") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","kg kg^{-1} m s^{-1}") , .true. )

        call disp_err( nf90_def_var(ncid,"psfc"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","surface pressure") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","Pa") , .true. )
      endif
      if(sfcmodel.eq.2.or.sfcmodel.eq.3)then
        call disp_err( nf90_def_var(ncid,"lu_index",nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"kpbl2d"  ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"hfx"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"qfx"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"hpbl"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"wspd"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"psim"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"psih"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"gz1oz0"  ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"br"      ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"chs"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"chs2"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"cqs2"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"cpmm"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"zol"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"mavail"  ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"mol"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"rmol"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"regime"  ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"lh"      ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"tmn"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"flhc"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"flqc"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"qgh"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"ck"      ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"cka"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"cda"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"ustm"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"qsfc"    ,nf90_float,(/niid,njid/),varid) , .true. )

        call disp_err( nf90_def_var(ncid,"t2"      ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","temperature at z=2m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","K") , .true. )

        call disp_err( nf90_def_var(ncid,"q2"      ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","water vapor mixing ratio at z=2m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )

        call disp_err( nf90_def_var(ncid,"th2"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","potential temperature at z=2m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","K") , .true. )

        call disp_err( nf90_def_var(ncid,"emiss"   ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"thc"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"albd"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"gsw"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"glw"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"chklowq" ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"capg"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"snowc"   ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"fm"      ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"fh"      ,nf90_float,(/niid,njid/),varid) , .true. )
        do n=1,num_soil_layers
          if( n.lt.10 )then
            text1 = 'tslbX   '
            write(text1(5:5),171) n
171         format(i1.1)
          elseif( n.lt.100 )then
            text1 = 'tslbXX  '
            write(text1(5:6),172) n
172         format(i2.2)
          else
            stop 22122
          endif
          call disp_err( nf90_def_var(ncid,text1     ,nf90_float,(/niid,njid/),varid) , .true. )
        enddo
      endif
      if(oceanmodel.eq.2)then
        call disp_err( nf90_def_var(ncid,"tml"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"t0ml"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"hml"     ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"h0ml"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"huml"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"hvml"    ,nf90_float,(/niid,njid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"tmoml"   ,nf90_float,(/niid,njid/),varid) , .true. )
      endif
    endif

    if( iptra.eq.1 )then
      do n=1,npt
        if( n.le.9 )then
          text1 = 'ptX     '
          write(text1(3:3),161) n
161       format(i1.1)
        elseif( n.lt.100 )then
          text1 = 'ptXX    '
          write(text1(3:4),162) n
162       format(i2.2)
        else
          stop 11511
        endif
        call disp_err( nf90_def_var(ncid,text1,nf90_float,(/niid,njid,nkid/),varid) , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      enddo
    endif

    if( iprcl.eq.1 )then
      call disp_err( nf90_def_var(ncid,"ploc"    ,nf90_float,(/plocid,nparcelsid/),varid) , .true. )
    endif

    if(irbc.eq.4)then
      if( wbc.eq.2 )then
        call disp_err( nf90_def_var(ncid,"radbcw"  ,nf90_float,(/njid,nkid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","estimated gravity wave phase speed on west boundary") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/ny,nz/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if( ebc.eq.2 )then
        call disp_err( nf90_def_var(ncid,"radbce"  ,nf90_float,(/njid,nkid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","estimated gravity wave phase speed on east boundary") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/ny,nz/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if( sbc.eq.2 )then
        call disp_err( nf90_def_var(ncid,"radbcs"  ,nf90_float,(/niid,nkid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","estimated gravity wave phase speed on south boundary") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/nx,nz/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if( nbc.eq.2 )then
        call disp_err( nf90_def_var(ncid,"radbcn"  ,nf90_float,(/niid,nkid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","estimated gravity wave phase speed on north boundary") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/nx,nz/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
    endif

  !------------------
  ! 150820:  optionals

    !-----
    IF( restart_file_theta )THEN
      call disp_err( nf90_def_var(ncid,"theta"   ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","potential temperature") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","K") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

    ENDIF
    IF( restart_file_dbz )THEN
      call disp_err( nf90_def_var(ncid,"dbz"     ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","reflectivity") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","dBZ") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

    ENDIF
    !-----
    IF( restart_file_th0 )THEN
      call disp_err( nf90_def_var(ncid,"th0"     ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","base-state potential temperature") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","K") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

    ENDIF
    IF( restart_file_prs0 )THEN
      call disp_err( nf90_def_var(ncid,"prs0"    ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","base-state pressure") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","Pa") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

    ENDIF
    IF( restart_file_pi0 )THEN
      call disp_err( nf90_def_var(ncid,"pi0"     ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","base-state nondimensional pressure") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

    ENDIF
    IF( restart_file_rho0 )THEN
      call disp_err( nf90_def_var(ncid,"rho0"    ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","base-state dry-air density") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/m^3") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

    ENDIF
    IF( restart_file_qv0 )THEN
      call disp_err( nf90_def_var(ncid,"qv0"     ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","base-state water vapor mixing ratio") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

    ENDIF
    IF( restart_file_u0 )THEN
      call disp_err( nf90_def_var(ncid,"u0"      ,nf90_float,(/nip1id,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","base-state x-component velocity") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkxp1,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

    ENDIF
    IF( restart_file_v0 )THEN
      call disp_err( nf90_def_var(ncid,"v0"      ,nf90_float,(/niid,njp1id,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","base-state y-component velocity") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chkyp1,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

    ENDIF
    !-----
    IF( restart_file_zs )THEN
      call disp_err( nf90_def_var(ncid,"zs"      ,nf90_float,(/niid,njid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","terrain height") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

    ENDIF
    IF( restart_file_zh )THEN
      call disp_err( nf90_def_var(ncid,"zhalf"   ,nf90_float,(/niid,njid,nkid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","height of half (scalar) grid points (3d array)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

    ENDIF
    IF( restart_file_zf )THEN
      call disp_err( nf90_def_var(ncid,"zfull"   ,nf90_float,(/niid,njid,nkp1id/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","height of full (w) grid points (3d array)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

      call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
      call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

    ENDIF
	
!-----
    IF( restart_file_diags )THEN
      if( td_diss.gt.0 )then
        call disp_err( nf90_def_var(ncid,"dissheat",nf90_float,(/niid,njid,nkid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","dissipative heating (potential temperature tendency)") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","K/s") , .true. )

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if( td_mptend.gt.0 )then
        call disp_err( nf90_def_var(ncid,"mptend",nf90_float,(/niid,njid,nkid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","potential temperature tendency from microphysics scheme") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","K/s") , .true. )

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if( qd_vtc.gt.0 )then
        call disp_err( nf90_def_var(ncid,"vtc",nf90_float,(/niid,njid,nkid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qc") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if( qd_vtr.gt.0 )then
        call disp_err( nf90_def_var(ncid,"vtr",nf90_float,(/niid,njid,nkid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qr") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if( qd_vts.gt.0 )then
        call disp_err( nf90_def_var(ncid,"vts",nf90_float,(/niid,njid,nkid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qs") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if( qd_vtg.gt.0 )then
        call disp_err( nf90_def_var(ncid,"vtg",nf90_float,(/niid,njid,nkid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qg") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
      if( qd_vti.gt.0 )then
        call disp_err( nf90_def_var(ncid,"vti",nf90_float,(/niid,njid,nkid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qi") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )

      endif
    ENDIF
    !-----

  !------------------
  ! end definitions:

  call disp_err( nf90_enddef(ncid) , .true. )

  !---------------------------------------------
  !ccccccccccccccccccccccccccccccccccccccccccccc
  !---------------------------------------------


      status = nf90_inq_varid(ncid,'xh',varid)
      do i=1,nx
        dumx(i) = 0.5*(xfref(i)+xfref(i+1))
      enddo
      status = nf90_put_var(ncid,varid,dumx,(/1/),(/nx/))


      status = nf90_inq_varid(ncid,'xf',varid)
      do i=1,nx+1
        dumx(i) = xfref(i)
      enddo
      status = nf90_put_var(ncid,varid,dumx,(/1/),(/nx+1/))


      status = nf90_inq_varid(ncid,'yh',varid)
      do j=1,ny
        dumy(j) = 0.5*(yfref(j)+yfref(j+1))
      enddo
      status = nf90_put_var(ncid,varid,dumy,(/1/),(/ny/))


      status = nf90_inq_varid(ncid,'yf',varid)
      do j=1,ny+1
        dumy(j) = yfref(j)
      enddo
      status = nf90_put_var(ncid,varid,dumy,(/1/),(/ny+1/))


      status = nf90_inq_varid(ncid,'zh',varid)
      if(terrain_flag)then
        do k=1,kval
          dumz(k) = sigma(k)
        enddo
        status = nf90_put_var(ncid,varid,dumz,(/1/),(/kval/))
      else
        do k=1,kval
          dumz(k) = zh(1,1,k)
        enddo
        status = nf90_put_var(ncid,varid,dumz,(/1/),(/kval/))
      endif


      status = nf90_inq_varid(ncid,'zf',varid)
      if(terrain_flag)then
        do k=1,kvalp1
          dumz(k) = sigmaf(k)
        enddo
        status = nf90_put_var(ncid,varid,dumz,(/1/),(/kvalp1/))
      else
        do k=1,kvalp1
          dumz(k) = zf(1,1,k)
        enddo
        status = nf90_put_var(ncid,varid,dumz,(/1/),(/kvalp1/))
      endif


!  END
!-----------------------------------------------------------------------

      call disp_err( nf90_inq_varid(ncid,'time',varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,actual_time,(/time_index/)) , .true. )

      print *,'  leaving restart_prelim '

    return
    end subroutine restart_prelim


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine checkstatus(status)
      use netcdf
      implicit none

      integer :: status

      if(status.ne.nf90_noerr)then
        print *,'  Error ... '
        print *,nf90_strerror(status)
        call stopcm1
      endif

      return
      end subroutine checkstatus

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine write2d_nc(chid,ncid,time_index,ni,nj,d2d)
      use netcdf
      implicit none

      character*8, intent(in) :: chid
      integer, intent(in) :: ncid,time_index,ni,nj
      real, dimension(ni,nj), intent(in) :: d2d

      integer :: varid,status

!----------------------------------

      status = nf90_inq_varid(ncid,chid,varid)
      if(status.ne.nf90_noerr)then
        print *,'  Error1 in write2d_nc, chid = ',chid
        print *,nf90_strerror(status)
        call stopcm1
      endif

      status = nf90_put_var(ncid,varid,d2d,(/1,1,time_index/),(/ni,nj,1/))
      if(status.ne.nf90_noerr)then
        print *,'  Error2 in write2d_nc, chid = ',chid
        print *,nf90_strerror(status)
        call stopcm1
      endif

!----------------------------------

      return
      end subroutine write2d_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine write3d_nc(chid,k,ncid,time_index,ni,nj,ds)
      use netcdf
      implicit none

      character*8, intent(in) :: chid
      integer, intent(in) :: k,ncid,time_index,ni,nj
      real, dimension(ni,nj), intent(in) :: ds

      integer :: varid,status

!----------------------------------

      status = nf90_inq_varid(ncid,chid,varid)
      if(status.ne.nf90_noerr)then
        print *,'  Error1 in write3d_nc, chid = ',chid
        print *,'  status = ',status
        print *,nf90_strerror(status)
        call stopcm1
      endif

      status = nf90_put_var(ncid,varid,ds,(/1,1,k,time_index/),(/ni,nj,1,1/))
      if(status.ne.nf90_noerr)then
        print *,'  Error2 in write3d_nc, chid = ',chid
        print *,'  status = ',status
        print *,nf90_strerror(status)
        call stopcm1
      endif

!----------------------------------

      return
      end subroutine write3d_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine writestat_nc(nrec,rtime,nstat,rstat,qname,budname)
      use netcdf
      implicit none

      include 'input.incl'

      integer, intent(inout) :: nrec
      real,    intent(in)    :: rtime
      integer, intent(in)    :: nstat
      real, dimension(stat_out), intent(in) :: rstat
      character*3, dimension(maxq), intent(in) :: qname
      character*6, dimension(maxq), intent(in) :: budname

      integer :: n,ncid,status,dimid,varid,time_index,timeid,tfile
      integer :: xhid,yhid,zhid
      character*8  :: text1
      character*30 :: text2
      logical :: allinfo

    string(totlen+1:totlen+22) = '_stats.nc             '

    print *,string

    allinfo = .false.
    IF(nrec.eq.1) allinfo=.true.

    IF( nrec.ne.1 )THEN
      ! cm1r18:  Try to open file.
      !          If error, set nrec to 1 and write all info.
      status = nf90_open(string,nf90_write,ncid)
      if( status.eq.nf90_noerr )then
        ! no error, file exists.  Get number of time levels in file:
        call disp_err( nf90_inq_dimid(ncid,'time',timeid) , .true. )
        call disp_err( nf90_inquire_dimension(ncid=ncid,dimid=timeid,len=tfile), .true. )
        if( (tfile+1).lt.nrec )then
          if(myid.eq.0) print *,'  tfile,nrec = ',tfile,nrec
          nrec = tfile+1
        endif
      else
        ! if error opening file, then write all info:
        if(myid.eq.0) print *,'  status = ',status
!!!        if(myid.eq.0) print *,nf90_strerror(status)
        allinfo = .true.
        nrec = 1
      endif
    ENDIF

    if( myid.eq.0 ) print *,'  nrec = ',nrec


  allinfo2:  IF( allinfo )THEN
    ! Definitions/descriptions:

      if( myid.eq.0 ) print *,'  calling nf90_create '

!--- works with netcdf 4.2, but not 4.0 (grumble)
      call disp_err( nf90_create(string,IOR(nf90_netcdf4, nf90_classic_model),ncid) , .true. )





    status = nf90_def_dim(ncid,"xh",1,xhid)
    status = nf90_def_dim(ncid,"yh",1,yhid)
    status = nf90_def_dim(ncid,"zh",1,zhid)
    status = nf90_def_dim(ncid,"time",nf90_unlimited,timeid)

    status = nf90_def_var(ncid,"xh",nf90_float,(/xhid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","west-east location")
    status = nf90_put_att(ncid,varid,"units","degree_east")

    status = nf90_def_var(ncid,"yh",nf90_float,(/yhid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","south-north location")
    status = nf90_put_att(ncid,varid,"units","degree_north")

    status = nf90_def_var(ncid,"zh",nf90_float,(/zhid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","height")
    status = nf90_put_att(ncid,varid,"units","m")

    status = nf90_def_var(ncid,"time",nf90_float,(/timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","time")
    status = nf90_put_att(ncid,varid,"units","seconds")

    IF(adapt_dt.eq.1)THEN
      status = nf90_def_var(ncid,"dt",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","average timestep dt           ")
      status = nf90_put_att(ncid,varid,"units","seconds")
    ENDIF
    IF(stat_w.eq.1)THEN
      status = nf90_def_var(ncid,"wmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","maximum vertical velocity     ")
      status = nf90_put_att(ncid,varid,"units","m/s")
      status = nf90_def_var(ncid,"wmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","minimum vertical velocity     ")
      status = nf90_put_att(ncid,varid,"units","m/s")
    ENDIF
    IF(stat_u.eq.1)THEN
      status = nf90_def_var(ncid,"umax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max E-W velocity              ")
      status = nf90_put_att(ncid,varid,"units","m/s")
      status = nf90_def_var(ncid,"umin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min E-W velocity              ")
      status = nf90_put_att(ncid,varid,"units","m/s")
      status = nf90_def_var(ncid,"sumax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max E-W velocity at lowest lvl")
      status = nf90_put_att(ncid,varid,"units","m/s")
      status = nf90_def_var(ncid,"sumin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min E-W velocity at lowest lvl")
      status = nf90_put_att(ncid,varid,"units","m/s")
    ENDIF
    IF(stat_v.eq.1)THEN
      status = nf90_def_var(ncid,"vmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max N-S velocity              ")
      status = nf90_put_att(ncid,varid,"units","m/s")
      status = nf90_def_var(ncid,"vmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min N-S velocity              ")
      status = nf90_put_att(ncid,varid,"units","m/s")
      status = nf90_def_var(ncid,"svmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max N-S velocity at lowest lvl")
      status = nf90_put_att(ncid,varid,"units","m/s")
      status = nf90_def_var(ncid,"svmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min N-S velocity at lowest lvl")
      status = nf90_put_att(ncid,varid,"units","m/s")
    ENDIF
    IF(stat_rmw.eq.1)THEN
      status = nf90_def_var(ncid,"rmw",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","radius of maximum windspeed")
      status = nf90_put_att(ncid,varid,"units","m")
      status = nf90_def_var(ncid,"zmw",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","height (ASL) of maximum windspeed")
      status = nf90_put_att(ncid,varid,"units","m")
    ENDIF
    IF(stat_pipert.eq.1)THEN
      status = nf90_def_var(ncid,"ppimax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max pi pert.                  ")
      status = nf90_put_att(ncid,varid,"units","nondimensional")
      status = nf90_def_var(ncid,"ppimin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min pi pert.                  ")
      status = nf90_put_att(ncid,varid,"units","nondimensional")
    ENDIF
    IF(stat_prspert.eq.1)THEN
      status = nf90_def_var(ncid,"ppmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max pressure pert.            ")
      status = nf90_put_att(ncid,varid,"units","Pa")
      status = nf90_def_var(ncid,"ppmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min pressure pert.            ")
      status = nf90_put_att(ncid,varid,"units","Pa")
    ENDIF
    IF(stat_thpert.eq.1)THEN
      status = nf90_def_var(ncid,"thpmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max potential temp. pert.     ")
      status = nf90_put_att(ncid,varid,"units","K")
      status = nf90_def_var(ncid,"thpmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min potential temp. pert.     ")
      status = nf90_put_att(ncid,varid,"units","K")
      status = nf90_def_var(ncid,"sthpmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max pot temp pert lowest level")
      status = nf90_put_att(ncid,varid,"units","K")
      status = nf90_def_var(ncid,"sthpmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min pot temp pert lowest level")
      status = nf90_put_att(ncid,varid,"units","K")
    ENDIF
    IF(stat_q.eq.1)THEN
      do n=1,numq
        text1='max     '
        text2='max                           '
        write(text1(4:6),156) qname(n)
        write(text2(5:7),156) qname(n)
156     format(a3)
        status = nf90_def_var(ncid,text1,nf90_float,timeid,varid)
        status = nf90_put_att(ncid,varid,"long_name",text2)
        status = nf90_put_att(ncid,varid,"units","kg/kg")
        text1='min     '
        text2='min                           '
        write(text1(4:6),156) qname(n)
        write(text2(5:7),156) qname(n)
        status = nf90_def_var(ncid,text1,nf90_float,timeid,varid)
        status = nf90_put_att(ncid,varid,"long_name",text2)
        status = nf90_put_att(ncid,varid,"units","kg/kg")
      enddo
    ENDIF
    IF(stat_tke.eq.1)THEN
      status = nf90_def_var(ncid,"tkemax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max tke                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s^2")
      status = nf90_def_var(ncid,"tkemin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min tke                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s^2")
    ENDIF
    IF(stat_km.eq.1)THEN
      status = nf90_def_var(ncid,"kmhmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max kmh                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
      status = nf90_def_var(ncid,"kmhmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min kmh                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
      status = nf90_def_var(ncid,"kmvmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max kmv                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
      status = nf90_def_var(ncid,"kmvmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min kmv                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
    ENDIF
    IF(stat_kh.eq.1)THEN
      status = nf90_def_var(ncid,"khhmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max khh                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
      status = nf90_def_var(ncid,"khhmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min khh                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
      status = nf90_def_var(ncid,"khvmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max khv                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
      status = nf90_def_var(ncid,"khvmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min khv                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
    ENDIF
    IF(stat_div.eq.1)THEN
      status = nf90_def_var(ncid,"divmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max 3d divergence             ")
      status = nf90_def_var(ncid,"divmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min 3d divergence             ")
    ENDIF
    IF(stat_rh.eq.1)THEN
      status = nf90_def_var(ncid,"rhmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max relative humidity         ")
      status = nf90_def_var(ncid,"rhmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min relative humidity         ")
    ENDIF
    IF(stat_rhi.eq.1)THEN
      status = nf90_def_var(ncid,"rhimax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max relative humidity wrt ice ")
      status = nf90_def_var(ncid,"rhimin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min relative humidity wrt ice ")
    ENDIF
    IF(iptra.eq.1)then
      do n=1,npt
        text1='maxpt   '
        text2='max pt                        '
        if( n.le.9 )then
          write(text1(6:6),157) n
          write(text2(7:7),157) n
157       format(i1)
        else
          write(text1(6:7),167) n
          write(text2(7:8),167) n
167       format(i2)
        endif
        status = nf90_def_var(ncid,text1,nf90_float,timeid,varid)
        status = nf90_put_att(ncid,varid,"long_name",text2)
        text1='minpt   '
        text2='min pt                        '
        if( n.le.9 )then
          write(text1(6:6),157) n
          write(text2(7:7),157) n
        else
          write(text1(6:7),167) n
          write(text2(7:8),167) n
        endif
        status = nf90_def_var(ncid,text1,nf90_float,timeid,varid)
        status = nf90_put_att(ncid,varid,"long_name",text2)
      enddo
    endif
    IF(stat_the.eq.1)THEN
      status = nf90_def_var(ncid,"themax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max theta-e below 10 km       ")
      status = nf90_put_att(ncid,varid,"units","K")
      status = nf90_def_var(ncid,"themin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min theta-e below 10 km       ")
      status = nf90_put_att(ncid,varid,"units","K")
      status = nf90_def_var(ncid,"sthemax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max theta-e at lowest level   ")
      status = nf90_put_att(ncid,varid,"units","K")
      status = nf90_def_var(ncid,"sthemin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min theta-e at lowest level   ")
      status = nf90_put_att(ncid,varid,"units","K")
    ENDIF
    IF(stat_cloud.eq.1)THEN
      status = nf90_def_var(ncid,"qctop",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max cloud top height          ")
      status = nf90_put_att(ncid,varid,"units","m")
      status = nf90_def_var(ncid,"qcbot",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min cloud base height         ")
      status = nf90_put_att(ncid,varid,"units","m")
    ENDIF
    IF(stat_sfcprs.eq.1)THEN
      status = nf90_def_var(ncid,"sprsmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max pressure at lowest level  ")
      status = nf90_put_att(ncid,varid,"units","Pa")
      status = nf90_def_var(ncid,"sprsmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min pressure at lowest level  ")
      status = nf90_put_att(ncid,varid,"units","Pa")
      status = nf90_def_var(ncid,"psfcmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max surface pressure")
      status = nf90_put_att(ncid,varid,"units","Pa")
      status = nf90_def_var(ncid,"psfcmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min surface pressure")
      status = nf90_put_att(ncid,varid,"units","Pa")
    ENDIF
    IF(stat_wsp.eq.1)THEN
      status = nf90_def_var(ncid,"wspmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max wind speed                ")
      status = nf90_put_att(ncid,varid,"units","m/s")

      status = nf90_def_var(ncid,"wspmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min wind speed                ")
      status = nf90_put_att(ncid,varid,"units","m/s")

      status = nf90_def_var(ncid,"swspmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max wind speed at lowest level")
      status = nf90_put_att(ncid,varid,"units","m/s")

      status = nf90_def_var(ncid,"swspmin",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min wind speed at lowest level")
      status = nf90_put_att(ncid,varid,"units","m/s")

    IF(bbc.eq.3)THEN
      status = nf90_def_var(ncid,"wsp10max",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max 10 m wind speed           ")
      status = nf90_put_att(ncid,varid,"units","m/s")

      status = nf90_def_var(ncid,"wsp10min",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","min 10 m wind speed           ")
      status = nf90_put_att(ncid,varid,"units","m/s")
    ENDIF
    ENDIF
    IF(stat_cfl.eq.1)THEN
    IF(adapt_dt.eq.1)THEN
      status = nf90_def_var(ncid,"cflmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max Courant number (average)  ")
      status = nf90_put_att(ncid,varid,"units","nondimensional")
    ELSE
      status = nf90_def_var(ncid,"cflmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max Courant number            ")
      status = nf90_put_att(ncid,varid,"units","nondimensional")
    ENDIF
      status = nf90_def_var(ncid,"kshmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max horiz K stability factor  ")
      status = nf90_put_att(ncid,varid,"units","nondimensional")
      status = nf90_def_var(ncid,"ksvmax",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max vert K stability factor   ")
      status = nf90_put_att(ncid,varid,"units","nondimensional")
    ENDIF
    IF(stat_vort.eq.1)THEN
      status = nf90_def_var(ncid,"vortsfc",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max vert. vort. at lowest lvl ")
      status = nf90_put_att(ncid,varid,"units","1/s")
      status = nf90_def_var(ncid,"vort1km",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max vert. vort. at z = 1 km   ")
      status = nf90_put_att(ncid,varid,"units","1/s")
      status = nf90_def_var(ncid,"vort2km",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max vert. vort. at z = 2 km   ")
      status = nf90_put_att(ncid,varid,"units","1/s")
      status = nf90_def_var(ncid,"vort3km",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max vert. vort. at z = 3 km   ")
      status = nf90_put_att(ncid,varid,"units","1/s")
      status = nf90_def_var(ncid,"vort4km",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max vert. vort. at z = 4 km   ")
      status = nf90_put_att(ncid,varid,"units","1/s")
      status = nf90_def_var(ncid,"vort5km",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","max vert. vort. at z = 5 km   ")
      status = nf90_put_att(ncid,varid,"units","1/s")
    ENDIF
    IF(stat_tmass.eq.1)THEN
      status = nf90_def_var(ncid,"tmass",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","total mass of (dry) air       ")
      status = nf90_put_att(ncid,varid,"units","kg")
    ENDIF
    IF(stat_tmois.eq.1)THEN
      status = nf90_def_var(ncid,"tmois",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","total moisture                ")
      status = nf90_put_att(ncid,varid,"units","kg")
    ENDIF
    IF(stat_qmass.eq.1)THEN
      do n=1,numq
        IF( (n.eq.nqv) .or.                                 &
            (n.ge.nql1.and.n.le.nql2) .or.                  &
            (n.ge.nqs1.and.n.le.nqs2.and.iice.eq.1) )THEN
          text1='mass    '
          text2='total mass of                 '
          write(text1( 5: 7),156) qname(n)
          write(text2(15:17),156) qname(n)
          status = nf90_def_var(ncid,text1,nf90_float,timeid,varid)
          status = nf90_put_att(ncid,varid,"long_name",text2)
          status = nf90_put_att(ncid,varid,"units","kg")
        ENDIF
      enddo
    ENDIF
    IF(stat_tenerg.eq.1)THEN
      status = nf90_def_var(ncid,"ek",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","total kinetic energy          ")
      status = nf90_def_var(ncid,"ei",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","total internal energy         ")
      status = nf90_def_var(ncid,"ep",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","total potential energy        ")
      status = nf90_def_var(ncid,"le",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","total latent energy           ")
      status = nf90_def_var(ncid,"et",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","total energy                  ")
    ENDIF
    IF(stat_mo.eq.1)THEN
      status = nf90_def_var(ncid,"tmu",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","total E-W momentum            ")
      status = nf90_def_var(ncid,"tmv",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","total N-S momentum            ")
      status = nf90_def_var(ncid,"tmw",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","total vertical momentum       ")
    ENDIF
    IF(stat_tmf.eq.1)THEN
      status = nf90_def_var(ncid,"tmfu",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","total upward mass flux        ")
      status = nf90_def_var(ncid,"tmfd",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"long_name","total downward mass flux      ")
    ENDIF
    IF(stat_pcn.eq.1)THEN
      do n=1,nbudget
        text1='        '
        text2='                              '
        write(text1(1:6),158) budname(n)
        write(text2(1:6),158) budname(n)
158     format(a6)
        status = nf90_def_var(ncid,text1,nf90_float,timeid,varid)
        status = nf90_put_att(ncid,varid,"long_name",text2)
      enddo
    ENDIF
    IF(stat_qsrc.eq.1)THEN
      do n=1,numq
        text1='as      '
        text2='artificial source of          '
        write(text1( 3: 5),156) qname(n)
        write(text2(22:24),156) qname(n)
        status = nf90_def_var(ncid,text1,nf90_float,timeid,varid)
        status = nf90_put_att(ncid,varid,"long_name",text2)
      enddo
      do n=1,numq
        text1='bs      '
        text2='bndry source/sink of          '
        write(text1( 3: 5),156) qname(n)
        write(text2(22:24),156) qname(n)
        status = nf90_def_var(ncid,text1,nf90_float,timeid,varid)
        status = nf90_put_att(ncid,varid,"long_name",text2)
      enddo
    ENDIF

    status = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions','COARDS')

    status = nf90_enddef(ncid)

    status = nf90_put_var(ncid,1,0.0)
    status = nf90_put_var(ncid,2,0.0)
    status = nf90_put_var(ncid,3,0.0)

  ELSE

    ! open file:

    call disp_err( nf90_open(string,nf90_write,ncid), .true. )

  ENDIF  allinfo2

    ! Write data:

    time_index = nrec

    call disp_err( nf90_inq_varid(ncid,'time',timeid) , .true. )
    call disp_err( nf90_put_var(ncid,timeid,rtime,(/time_index/)) , .true. )

    DO n=1,nstat
      varid = timeid + n
      status = nf90_put_var(ncid,varid,rstat(n),(/time_index/))
    ENDDO

    ! close file

    call disp_err( nf90_close(ncid) , .true. )

    nrec = nrec + 1

    ! all done

      return
      end subroutine writestat_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine disp_err( status , stop_on_error )
      use netcdf
      implicit none

      integer, intent(in) :: status
      logical, intent(in) :: stop_on_error

      IF( status.ne.nf90_noerr )THEN
        IF( stop_on_error )THEN
          print *,'  netcdf status returned an error: ', status,' ... stopping program'
          print *
          print *,nf90_strerror(status)
          print *
          call stopcm1
        ENDIF
      ENDIF

      return
      end subroutine disp_err

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine writepdata_nc(prec,rtime,qname,pdata,pdata2)
      use netcdf
      implicit none

      include 'input.incl'

      integer, intent(inout) :: prec
      real, intent(in) :: rtime
      character*3, intent(in), dimension(maxq) :: qname
      real, intent(in), dimension(npvals,nparcels) :: pdata
      real, intent(inout), dimension(nparcels) :: pdata2

      integer :: ncid,status,dimid,varid,time_index,n,n2,np,timeid,tfile,xid
      integer :: npid,yhid,zhid
      character*8 text1
      character*30 text2
      logical :: allinfo

!-----------------------------------------------------------------------

    string(totlen+1:totlen+22) = '_pdata.nc             '

    allinfo = .false.
    IF(prec.eq.1) allinfo=.true.

    IF( prec.ne.1 )THEN
      ! cm1r18:  Try to open file.
      !          If error, set prec to 1 and write all info.
      status = nf90_open(string,nf90_write,ncid)
      if( status.eq.nf90_noerr )then
        ! no error, file exists.  Get number of time levels in file:
        call disp_err( nf90_inq_dimid(ncid,'time',timeid) , .true. )
        call disp_err( nf90_inquire_dimension(ncid=ncid,dimid=timeid,len=tfile), .true. )
        if( (tfile+1).lt.prec )then
          if(myid.eq.0) print *,'  tfile,prec = ',tfile,prec
          prec = tfile+1
        endif
      else
        ! if error opening file, then write all info:
        if(myid.eq.0) print *,'  status = ',status
!!!        if(myid.eq.0) print *,nf90_strerror(status)
        allinfo = .true.
        prec = 1
      endif
    ENDIF

    if( myid.eq.0 ) print *,'  pdata prec = ',prec


  allinfo3:  IF( allinfo )THEN
    ! Definitions/descriptions:


!--- works with netcdf 4.2, but not 4.0 (grumble)
      call disp_err( nf90_create('cm1out_pdata.nc',IOR(nf90_netcdf4, nf90_classic_model),ncid) , .true. )





    status = nf90_def_dim(ncid,"xh",nparcels,npid)
    status = nf90_def_dim(ncid,"yh",1,yhid)
    status = nf90_def_dim(ncid,"zh",1,zhid)
    status = nf90_def_dim(ncid,"time",nf90_unlimited,timeid)

    status = nf90_def_var(ncid,"xh",nf90_float,(/npid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","west-east location ... actually, really parcel ID number")
    status = nf90_put_att(ncid,varid,"units","degree_east")

    status = nf90_def_var(ncid,"yh",nf90_float,(/yhid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","south-north location")
    status = nf90_put_att(ncid,varid,"units","degree_north")

    status = nf90_def_var(ncid,"zh",nf90_float,(/zhid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","height")
    status = nf90_put_att(ncid,varid,"units","m")

    status = nf90_def_var(ncid,"time",nf90_float,(/timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","time")
    status = nf90_put_att(ncid,varid,"units","seconds")

!------------------------

    status = nf90_def_var(ncid,"x",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","x                             ")
    status = nf90_put_att(ncid,varid,"units","m")

    status = nf90_def_var(ncid,"y",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","y                             ")
    status = nf90_put_att(ncid,varid,"units","m")

    status = nf90_def_var(ncid,"z",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","z                             ")
    status = nf90_put_att(ncid,varid,"units","m")

    status = nf90_def_var(ncid,"u",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","u                             ")
    status = nf90_put_att(ncid,varid,"units","m/s")

    status = nf90_def_var(ncid,"v",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","v                             ")
    status = nf90_put_att(ncid,varid,"units","m/s")

    status = nf90_def_var(ncid,"w",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","w                             ")
    status = nf90_put_att(ncid,varid,"units","m/s")

  if(prth.ge.1)then
    status = nf90_def_var(ncid,"th",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","potential temperature         ")
    status = nf90_put_att(ncid,varid,"units","K")
  endif
  if(prt.ge.1)then
    status = nf90_def_var(ncid,"t",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","temperature")
    status = nf90_put_att(ncid,varid,"units","K")
  endif
  if(prprs.ge.1)then
    status = nf90_def_var(ncid,"prs",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","pressure                      ")
    status = nf90_put_att(ncid,varid,"units","Pa")
  endif

    if(prpt1.ge.1)then
      do n=1,npt
        text1='pt      '
        if(n.le.9)then
          write(text1(3:3),155) n
        else
          write(text1(3:4),154) n
        endif
154     format(i2.2)
155     format(i1.1)
        status = nf90_def_var(ncid,text1,nf90_float,(/npid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"long_name","passive tracer conc.")
        status = nf90_put_att(ncid,varid,"units","kg/kg")
      enddo
    endif

  if(prqv.ge.1)then
    status = nf90_def_var(ncid,"qv",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","water vapor mixing ratio")
    status = nf90_put_att(ncid,varid,"units","kg/kg")
  endif

        if(prq1.ge.1)then
          n2 = nql2
          if( iice.eq.1 ) n2 = nqs2
          do n=nql1,n2
          text1='        '
          text2='                              '
          write(text1(1:3),156) qname(n)
          write(text2(1:3),156) qname(n)
156       format(a3)
          status = nf90_def_var(ncid,text1,nf90_float,(/npid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name",text2)
          enddo
        endif

        if(prnc1.ge.1)then
          do n=nnc1,nnc2
          text1='        '
          text2='                              '
          write(text1(1:3),156) qname(n)
          write(text2(1:3),156) qname(n)
          status = nf90_def_var(ncid,text1,nf90_float,(/npid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"long_name",text2)
          enddo
        endif

  if(prkm.ge.1)then
    status = nf90_def_var(ncid,"kmh",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","turb. coef. for momentum (horiz)")
    status = nf90_put_att(ncid,varid,"units","m^2/s")
  endif
  if(prkm.ge.1)then
    status = nf90_def_var(ncid,"kmv",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","turb. coef. for momentum (vert)")
    status = nf90_put_att(ncid,varid,"units","m^2/s")
  endif
  if(prkh.ge.1)then
    status = nf90_def_var(ncid,"khh",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","turb. coef. for scalar (horiz)")
    status = nf90_put_att(ncid,varid,"units","m^2/s")
  endif
  if(prkh.ge.1)then
    status = nf90_def_var(ncid,"khv",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","turb. coef. for scalar (vert)")
    status = nf90_put_att(ncid,varid,"units","m^2/s")
  endif

  if(prtke.ge.1)then
    status = nf90_def_var(ncid,"tke",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","subgrid turbulence kinetic energy")
    status = nf90_put_att(ncid,varid,"units","m^2/s^2")
  endif

  if(prdbz.ge.1)then
    status = nf90_def_var(ncid,"dbz",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","reflectivity")
    status = nf90_put_att(ncid,varid,"units","dbz")
  endif
  if(prb.ge.1)then
    status = nf90_def_var(ncid,"b",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","buoyancy")
    status = nf90_put_att(ncid,varid,"units","m/s^2")
  endif
  if(prvpg.ge.1)then
    status = nf90_def_var(ncid,"vpg",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","vertical perturbation pressure gradient")
    status = nf90_put_att(ncid,varid,"units","m/s^2")
  endif
  if(przv.ge.1)then
    status = nf90_def_var(ncid,"zv",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","vertical vorticity")
    status = nf90_put_att(ncid,varid,"units","s^-1")
  endif
  if(prrho.ge.1)then
    status = nf90_def_var(ncid,"rho",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","dry-air density")
    status = nf90_put_att(ncid,varid,"units","kg/m^3")
  endif
  if(prqsl.ge.1)then
    status = nf90_def_var(ncid,"qsl",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","saturation mixing ratio wrt liquid")
    status = nf90_put_att(ncid,varid,"units","kg/kg")
  endif
  if(prqsi.ge.1)then
    status = nf90_def_var(ncid,"qsi",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","saturation mixing ratio wrt ice")
    status = nf90_put_att(ncid,varid,"units","kg/kg")
  endif
  if(prznt.ge.1)then
    status = nf90_def_var(ncid,"znt",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","sfc roughness length")
    status = nf90_put_att(ncid,varid,"units","m")
  endif
  if(prust.ge.1)then
    status = nf90_def_var(ncid,"ust",nf90_float,(/npid,timeid/),varid)
    status = nf90_put_att(ncid,varid,"long_name","sfc friction velocity")
    status = nf90_put_att(ncid,varid,"units","m/s")
  endif

!------------------------

    status = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions','COARDS')

    status = nf90_enddef(ncid)

  do np=1,nparcels
    status = nf90_put_var(ncid,npid,float(np),(/np/))
  enddo
    status = nf90_put_var(ncid,yhid,0.0)
    status = nf90_put_var(ncid,zhid,0.0)

!------------------------

  ENDIF  allinfo3

      ! Write data:

      time_index = prec

      call disp_err( nf90_inq_varid(ncid,'time',timeid) , .true. )
      call disp_err( nf90_put_var(ncid,timeid,rtime,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,'x',xid) , .true. )

      DO n=1,npvals
        varid = (xid-1) + n
        do np=1,nparcels
          pdata2(np) = pdata(n,np)
        enddo
        call disp_err( nf90_put_var(ncid,varid,pdata2,(/1,time_index/),(/nparcels,1/)) , .true. )
      ENDDO

      ! close file

      call disp_err( nf90_close(ncid) , .true. )

      prec = prec + 1

      ! all done

      return
      end subroutine writepdata_nc

