!############################################################################
!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######        SUBROUTINE READ_LETKF_WEIGHTS                 ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!
!
!############################################################################
!
!     PURPOSE:
!
!     Writes out the 3D mean and 4D transformation weights from LETKF
!
!     Author:  Lou Wicker
!
!     Creation Date:  April 2014
!
!############################################################################

      SUBROUTINE READ_LETKF_WEIGHTS(mean, trans, nx, ny, nz, ne)

      USE netcdf

      implicit none

!---- input parameters

      integer nx, ny, nz, ne                  ! dimensions of the arrays
      real(kind=4) mean(ne, nz, ny, nx)       ! weights for field mean
	real(kind=4) trans(ne, ne, nz, ny, nx)  ! weights for full ensemble

!---- local variables

      integer ncid                    ! netCDF file ID
	  
!---- variable IDs

      integer mean_id
      integer trans_id

!---- Interface block for SUBROUTINE CHECK

      INTERFACE
        SUBROUTINE CHECK(status, message)
          integer, intent (in) :: status
          character(len=*), optional :: message
        END SUBROUTINE CHECK
      END INTERFACE

!---- open file

      call check( nf90_open("weightsA.nc", NF90_NOWRITE, ncid) )

!---- Get the varid of the data variable, based on its name.

      call check( nf90_inq_varid(ncid, "trans", trans_id) )
      call check( nf90_inq_varid(ncid, "mean", mean_id) )

!---- Read the data.

      call check( nf90_get_var(ncid, trans_id, trans) )
      call check( nf90_get_var(ncid, mean_id, mean) )

!---- Close the file, freeing all resources.

      call check( nf90_close(ncid) )

      write(6,*) 
      write(6,*) "*** SUCCESS reading file weightsA.nc"
      write(6,*) 

      return
      end

!############################################################################
!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######        SUBROUTINE WRITE_LETKF_WEIGHTS                ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!
!
!############################################################################
!
!     PURPOSE:
!
!     Writes out the 3D mean and 4D transformation weights from LETKF
!
!     Author:  Lou Wicker
!
!     Creation Date:  April 2014
!
!############################################################################

      SUBROUTINE WRITE_LETKF_WEIGHTS(mean, trans, xc, yc, zc, time, nx, ny, nz, ne)

      USE netcdf

      implicit none

!---- input parameters

      integer nx, ny, nz, ne                  ! dimensions of the arrays
      real(kind=4) mean(ne, nz, ny, nx)       ! weights for field mean
	real(kind=4) trans(ne, ne, nz, ny, nx)  ! weights for full ensemble
      real(kind=4) xc(nx), yc(ny), zc(nz)
      real(kind=8) time

!---- local variables

      integer ncid                    ! netCDF file ID
	  
!---- dimension IDs

      integer nx_dim      
      integer ny_dim
      integer nz_dim
      integer ne1_dim      
      integer ne2_dim

!---- variable IDs

      integer mean_id
      integer trans_id
      integer xc_id, yc_id, zc_id

!---- variable shapes

      integer four_dims(4), five_dims(5)

! DEFLATE LEVEL

      integer, parameter :: DEFLATE_LEVEL = 4
	  
!---- Interface block for SUBROUTINE CHECK

      INTERFACE
        SUBROUTINE CHECK(status, message)
          integer, intent (in) :: status
          character(len=*), optional :: message
        END SUBROUTINE CHECK
      END INTERFACE

!---- enter define mode

      call check( nf90_create("weights.nc", NF90_NETCDF4, ncid) )

!---- define dimensions

      call check( nf90_def_dim(ncid, 'NX',  nx,   nx_dim) )
      call check( nf90_def_dim(ncid, 'NY',  ny,   ny_dim) )
      call check( nf90_def_dim(ncid, 'NZ',  nz,   nz_dim) )
      call check( nf90_def_dim(ncid, 'NE1',  ne,  ne1_dim) )
      call check( nf90_def_dim(ncid, 'NE2',  ne,  ne2_dim) )
	 
!---- define variables and attributes

      four_dims(1) = ne1_dim
      four_dims(2) = nz_dim
      four_dims(3) = ny_dim
      four_dims(4) = nx_dim
			
      call check( nf90_def_var(ncid, "xc", nf90_float, four_dims(4), xc_id, &
                 shuffle = .TRUE., fletcher32 = .TRUE., deflate_level = DEFLATE_LEVEL), &
                  message='NC_DEFVAR xc failed' )
	 
      call check( nf90_def_var(ncid, "yc", nf90_float, four_dims(3), yc_id, &
                 shuffle = .TRUE., fletcher32 = .TRUE., deflate_level = DEFLATE_LEVEL), &
                  message='NC_DEFVAR yc failed' )
	 
      call check( nf90_def_var(ncid, "zc", nf90_float, four_dims(2), zc_id, &
                 shuffle = .TRUE., fletcher32 = .TRUE., deflate_level = DEFLATE_LEVEL), &
                  message='NC_DEFVAR zc failed' )
	 
      call check( nf90_def_var(ncid, "mean", nf90_float, four_dims, mean_id, &
                 shuffle = .TRUE., fletcher32 = .TRUE., deflate_level = DEFLATE_LEVEL), &
                  message='NC_DEFVAR mean failed' )
	 
      five_dims(1) = ne1_dim
      five_dims(2) = ne2_dim
      five_dims(3) = nz_dim
      five_dims(4) = ny_dim
	five_dims(5) = nx_dim
	 
	call check( nf90_def_var(ncid, "trans", nf90_float, five_dims, trans_id, &
                 shuffle = .TRUE., fletcher32 = .TRUE., deflate_level = DEFLATE_LEVEL), &
                  message='NC_DEFVAR trans failed' )

!---- leave define mode

	call check(nf90_enddef(ncid), message='NC_ENDDEF failed') 

!---- write global attribute for time

      call check( nf90_put_att(ncid, NF90_GLOBAL, "Time", time), message='NC_PUTATTR time failed')

!---- write variables

      call check( nf90_put_var(ncid, xc_id, xc), message='NC_PUTVAR xc failed' )
      call check( nf90_put_var(ncid, yc_id, yc), message='NC_PUTVAR yc failed' )
      call check( nf90_put_var(ncid, xc_id, zc), message='NC_PUTVAR zc failed' )

      call check( nf90_put_var(ncid, mean_id, mean), message='NC_PUTVAR mean failed' )
      call check( nf90_put_var(ncid, trans_id, trans), message='NC_PUTVAR trans failed' )

	call check( nf90_close(ncid), message='NC_CLOSE failed')

  RETURN

  END SUBROUTINE WRITE_LETKF_WEIGHTS
 
  !############################################################################
  !
  !     ##################################################################
  !     ##################################################################
  !     ######                                                      ######
  !     ######                SUBROUTINE CHECK                      ######
  !     ######                                                      ######
  !     ##################################################################
  !     ##################################################################
  !
  !
  !############################################################################
  SUBROUTINE CHECK(status, message)

    USE NETCDF
    integer, intent (in) :: status
    character(len=*), optional :: message

    IF( status /= nf90_noerr ) THEN
      IF( present(message) ) THEN
        write(6,*) message//"  "//trim(nf90_strerror(status))
        stop "Stopped by Subroutine CHECK"
      ELSE
        write(6,*) trim(nf90_strerror(status))
        stop "Stopped by Subroutine CHECK"
      ENDIF
    END IF

  RETURN
  END SUBROUTINE CHECK
