! On Snellius compile this file with the command:
! gfortran -O3 -fconvert=big-endian -I/sw/arch/RHEL8/EB_production/2022/software/netCDF-Fortran/4.6.0-gompi-2022a/include -o calculate_msf_1deg calculate_msf_1deg.f -lnetcdf -lnetcdff

      program calc_msf

!-------------------------------------------------------------------
!  calculate meridional overturning streamfunction and write
!    out in netcdf format.
!-------------------------------------------------------------------

      implicit none

      integer, parameter ::  num_2d_fields = 3
     &  , num_1d_fields = 3
     &  , num_tries_max = 5

      integer :: imt, jmt, km, i, j, k, m, ncid, ierr, dim, ny_mht
     &  , u_first_record, v_first_record, nrec, rec_length
     &  , uet_first_record, vnt_first_record, num_tries, ip1
     &  , status

      integer, dimension(2) :: dims, start, count
      integer, dimension(1) :: latid, latvid, depid, depvid

      integer, dimension(num_2d_fields) :: fieldid_2d
      integer, dimension(num_1d_fields) :: fieldid_1d

      integer, allocatable, dimension(:,:) :: KMT

      real, parameter ::  special_value = -1.e34

      real, allocatable, dimension(:) :: lat_mht, mht_temp
      real, allocatable, dimension(:,:) :: mht
      real, allocatable, dimension(:,:) :: PSI_TEMP, PSI_TEMP_OLD
      real, allocatable, dimension(:,:,:) :: PSI

      real, allocatable, dimension(:,:,:) :: U, V, W, UET, VNT
     &, DZT, DZU

      real :: ysouth_mht, ynorth_mht, dy_mht, pi, dzw

      real, allocatable, dimension(:) :: dz, depth

      real, allocatable, dimension(:,:) :: 
     &  TLAT, TLONG, DXU, DYU, TAREA

      double precision, allocatable, dimension(:,:) :: 
     &  ULAT, ULONG, DZBC

      character (len=200) attr_name, attr_valc, dim_name, var_name
     &  ,kmt_file, kmt_atl_file, kmt_indopac_file, in_depths, grid_file
     &  ,u_file, v_file, uet_file, vnt_file, out_file, pbcfile

      logical :: verbose, do_msf, do_mht

      namelist /msf_stuff/ imt, jmt, km
     &  ,ysouth_mht, ynorth_mht, dy_mht
     &  ,kmt_file, kmt_atl_file, kmt_indopac_file, in_depths, grid_file
     &  ,u_file, v_file, out_file, u_first_record, v_first_record
     &  ,uet_file, vnt_file, uet_first_record, vnt_first_record
     &  ,verbose, do_msf, do_mht, pbcfile

      include 'netcdf.inc'

!-------------------------------------------------------------------
!-------------------------------------------------------------------

      pi = 4.0*atan(1.)

!-------------------------------------------------------------------
!  read in namelist
!-------------------------------------------------------------------

      call nml_defaults (
     &   imt, jmt, km, ysouth_mht, ynorth_mht, dy_mht
     &  ,kmt_file, kmt_atl_file, kmt_indopac_file, in_depths, grid_file
     &  ,u_file, v_file, u_first_record, v_first_record
     &  ,uet_file, vnt_file, uet_first_record, vnt_first_record
     &  , out_file, verbose, do_msf, do_mht, pbcfile)

      !open(1, file = 'in_msf', status = 'old')
      !read(1, nml=msf_stuff)
      !close(1)
      read(unit=5, nml=msf_stuff) !, iostat=status, err=10, end=26)

      if (verbose) write(*,nml=msf_stuff)

!-------------------------------------------------------------------
!  read in horizontal grid
!-------------------------------------------------------------------
      allocate(ULAT(imt,jmt), ULONG(imt,jmt),
     &         TLAT(imt,jmt), TLONG(imt,jmt),
     &         TAREA(imt,jmt), DXU(imt,jmt)
     &       , DYU(imt,jmt), KMT(imt,jmt))

      call grid_stuff(grid_file, ULAT, ULONG, DXU, DYU, TAREA, imt, jmt)

      ULAT  = ULAT *180.0/pi
      ULONG = ULONG*180.0/pi

      call calc_tpoints(TLONG,TLAT, ULONG, ULAT, imt, jmt, pi)

      deallocate(ULAT, ULONG, TLONG)  !  not needed anymore

!-------------------------------------------------------------------
!  set up vertical grid
!  note that depth is for W points
!-------------------------------------------------------------------

      allocate (dz(km), depth(km))
      open(1,file=in_depths,status='old')
      do k = 1,km
        read(1,*) dz(k)
      enddo
      close(1)
!     depth(1) = dz(1)
      depth(1) = 0.
      do k = 2, km
        depth(k) = depth(k-1) + dz(k-1)
      enddo
!  try vertical T grid
!     depth(1) = 0.5*dz(1)
!     do k = 1, km-1
!       dzw = 0.5*(dz(k) + dz(k+1))
!       depth(k+1) = depth(k) + dzw
!     enddo
      depth = depth*0.01  ! convert to m

!-------------------------------------------------------------------
!  read in KMT
!-------------------------------------------------------------------

      inquire (iolength=rec_length) KMT
      open(1,file=kmt_file,access='direct',form='unformatted',
     &          recl=rec_length,status='unknown')
      read(1,rec=1)KMT
      close(1)
      write(*,*)' done reading in KMT'

!-----------------------------------------------------------------------
!  read PBC file
!  28 feb 2014: Michael Kliphuis set DZBC = 0 for 1 degree cases because
!  there are no partial bottom cells (PBC = .false.)
!-----------------------------------------------------------------------
      allocate(DZBC(imt,jmt),DZT(imt,jmt,km),DZU(imt,jmt,km))

! inquire(iolength = rec_length) DZBC
! MK  open(1,file=pbcfile,access='direct',form='unformatted',   
! MK &          recl=rec_length,status='unknown')
! MK  read(1,rec=1) DZBC
! MK  close(1)
! MK  write(*,*)' read file: ',pbcfile

      DZBC = 0
      do k=1,km
         do j=1,jmt
         do i=1,imt
            if (KMT(i,j) == k) then
               DZT(i,j,k) = DZBC(i,j)
            else
               DZT(i,j,k) = dz(k)
            end if
         end do
         end do

         !*** DZU = min of surrounding DZTs

         do j=1,jmt-1
         do i=1,imt
            ip1 = merge(1,i+1,i==imt)
            DZU(i,j,k) = min(DZT(i  ,j  ,k),  
     &                        DZT(ip1,j  ,k),  
     &                        DZT(i  ,j+1,k),  
     &                        DZT(ip1,j+1,k))
         end do
         end do

! assume top row is land
         DZU(:,jmt,k) = 0.

      end do

      deallocate(DZBC,DZT)

!-------------------------------------------------------------------
!  if doing MSF, read in U, V, then calculate W
!-------------------------------------------------------------------

      if( do_msf ) then

        allocate (U(imt,jmt,km))
        write(*,*)' allocated U'
        allocate (V(imt,jmt,km))
        write(*,*)' allocated V'
        allocate (W(imt,jmt,km))
        write(*,*)' allocated W'

        inquire (iolength=rec_length) U(:,:,1)

        num_tries = 0
  101   nrec = u_first_record
        num_tries = num_tries + 1
        if (num_tries == num_tries_max) then
          write(*,*)' Giving up on reading U;  QUITTING'
          go to 999
        endif
        if (num_tries /= 1) then
          write(*,*)' I/O error reading U, trying again'
          close(1)
        endif
        open(1,file=u_file,access='direct',form='unformatted',
     &          recl=rec_length,status='unknown')
        do k = 1, km
          read(1,rec=nrec,err=101)U(:,:,k)
          nrec = nrec + 1
        enddo
        close(1)
        write(*,*)' done reading in U'

        num_tries = 0
  102   nrec = v_first_record
        num_tries = num_tries + 1
        if (num_tries == num_tries_max) then
          write(*,*)' Giving up on reading V;  QUITTING'
          go to 999
        endif
        if (num_tries /= 1) then
          write(*,*)' I/O error reading V, trying again'
          close(1)
        endif
        open(1,file=v_file,access='direct',form='unformatted',
     &          recl=rec_length,status='unknown')
        do k = 1, km
          read(1,rec=nrec,err=102)V(:,:,k)
          nrec = nrec + 1
        enddo
        close(1)
        write(*,*)' done reading in V'

        call wcalc(W, U, V, DXU, DYU, TAREA, KMT, DZU, dz, imt, jmt, km)

        deallocate (U)

      endif

!-------------------------------------------------------------------
!  if doing MHT, read in UET, VNT
!-------------------------------------------------------------------

      if( do_mht ) then

        allocate (UET(imt,jmt,km), VNT(imt,jmt,km))

        inquire (iolength=rec_length) UET(:,:,1)

        num_tries = 0
  103   nrec = uet_first_record
        num_tries = num_tries + 1
        if (num_tries == num_tries_max) then
          write(*,*)' Giving up on reading UET;  QUITTING'
          go to 999
        endif
        if (num_tries /= 1) then
          write(*,*)' I/O error reading UET, trying again'
          close(1)
        endif
        open(1,file=uet_file,access='direct',form='unformatted',
     &          recl=rec_length,status='unknown')
        do k = 1, km
          read(1,rec=nrec,err=103)UET(:,:,k)
          nrec = nrec + 1
        enddo
        close(1)
        write(*,*)' done reading in UET'

        num_tries = 0
  104   nrec = vnt_first_record
        num_tries = num_tries + 1
        if (num_tries == num_tries_max) then
          write(*,*)' Giving up on reading VNT;  QUITTING'
          go to 999
        endif
        if (num_tries /= 1) then
          write(*,*)' I/O error reading VNT, trying again'
          close(1)
        endif
        open(1,file=vnt_file,access='direct',form='unformatted',
     &          recl=rec_length,status='unknown')
        do k = 1, km
          read(1,rec=nrec,err=104)VNT(:,:,k)
          nrec = nrec + 1
        enddo
        close(1)
        write(*,*)' done reading in VNT'

      endif

!-------------------------------------------------------------------
!  set up latitude grid to be used and allocate arrays
!-------------------------------------------------------------------

      ny_mht = nint( (ynorth_mht - ysouth_mht)/dy_mht ) + 1
      write(*,*)'  using ',ny_mht,' lat points for merid streamfunction'

      allocate(lat_mht(ny_mht))

      do j = 1, ny_mht
        lat_mht(j) = ysouth_mht + (j-1)*dy_mht
      enddo

      if(do_mht) allocate(mht(ny_mht,num_1d_fields), mht_temp(ny_mht))

      if(do_msf) allocate(PSI_TEMP_OLD(km,ny_mht), PSI_TEMP(ny_mht,km)
     &       , PSI(ny_mht,km,num_2d_fields))

!-------------------------------------------------------------------
!  calculate overturning streamfunctions and heat transports
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  first do global
!-------------------------------------------------------------------

      inquire (iolength=rec_length) KMT
      open(44,file=kmt_file,access='direct',form='unformatted',
     &          recl=rec_length,status='unknown')
      read(44,rec=1)KMT
      close(44)
      
      if(do_msf) then
        call moc(V,W,KMT,TLAT,DXU,TAREA,DZU,dz
     &        ,imt,jmt,km,lat_mht,ny_mht,PSI_TEMP_OLD)
        do k = 1, km
          PSI(:,k,1) = PSI_TEMP_OLD(k,:)
        enddo
      endif

      if(do_mht) then
        call meridional_heat(UET,VNT,KMT,TLAT,TAREA,dz
     &  ,imt,jmt,km,lat_mht,ny_mht,mht_temp)
        mht(:,1) = mht_temp
      endif

!-------------------------------------------------------------------
!  next do atlantic
!-------------------------------------------------------------------
      open(44,file=kmt_atl_file,access='direct',form='unformatted',
     &          recl=rec_length,status='unknown')
      read(44,rec=1)KMT
      close(44)

      if(do_msf) then
        call moc(V,W,KMT,TLAT,DXU,TAREA,DZU,dz
     &        ,imt,jmt,km,lat_mht,ny_mht,PSI_TEMP_OLD)
        do k = 1, km
          PSI(:,k,2) = PSI_TEMP_OLD(k,:)
        enddo
      endif

      if(do_mht) then
        call meridional_heat(UET,VNT,KMT,TLAT,TAREA,dz
     &  ,imt,jmt,km,lat_mht,ny_mht,mht_temp)
        mht(:,2) = mht_temp
      endif

!-------------------------------------------------------------------
!  now do indo-pacific
!-------------------------------------------------------------------
      open(44,file=kmt_indopac_file,access='direct',form='unformatted',
     &          recl=rec_length,status='unknown')
      read(44,rec=1)KMT
      close(44)

      if(do_msf) then
        call moc(V,W,KMT,TLAT,DXU,TAREA,DZU,dz
     &        ,imt,jmt,km,lat_mht,ny_mht,PSI_TEMP_OLD)
        do k = 1, km
          PSI(:,k,3) = PSI_TEMP_OLD(k,:)
        enddo
      endif

      if(do_mht) then
        call meridional_heat(UET,VNT,KMT,TLAT,TAREA,dz
     &  ,imt,jmt,km,lat_mht,ny_mht,mht_temp)
        mht(:,3) = mht_temp
      endif

!-------------------------------------------------------------------
!  start netcdf stuff
!-------------------------------------------------------------------

      ncid = NCCRE(out_file, NCCLOB, ierr)

!-------------------------------------------------------------------
!  define depth and latitude grids
!-------------------------------------------------------------------

      dim_name = 'depth_w'
      depid(1) = NCDDEF(ncid,dim_name,km,ierr)

      dim_name = 'lat_mht'
      latid(1) = NCDDEF(ncid,dim_name,ny_mht,ierr)

!-------------------------------------------------------------------
!  define grid variables
!-------------------------------------------------------------------

      var_name = 'depth_w'
      dim = 1
      dims(1) = depid(1)
      depvid(1) = NCVDEF(ncid,var_name,NCFLOAT,dim,dims,ierr)
      attr_name = 'units'
      attr_valc = 'm'
      call NCAPTC(ncid,depvid(1),attr_name,NCCHAR,100,attr_valc,ierr)

      var_name = 'lat_mht'
      dim = 1
      dims(1) = latid(1)
      latvid(1) = NCVDEF(ncid,var_name,NCFLOAT,dim,dims,ierr)
      attr_name = 'units'
      attr_valc = 'latitude'
      call NCAPTC(ncid,latvid(1),attr_name,NCCHAR,100,attr_valc,ierr)

!-------------------------------------------------------------------
!  define 1d variables
!-------------------------------------------------------------------

      if(do_mht) then
        dim = 1
        dims(1) = latid(1)
        do m = 1,num_1d_fields
           write(*,*)' m = ',m
           if(m == 1) then
              var_name = 'MHTG'
           elseif(m == 2) then
              var_name = 'MHTA'
           elseif(m == 3) then
              var_name = 'MHTIP'
           endif
           fieldid_1d(m) = NCVDEF(ncid,var_name,NCFLOAT,dim,dims,ierr)
           attr_name = 'units'
           attr_valc = 'PW'
           call NCAPTC(ncid,fieldid_1d(m),attr_name
     &              ,NCCHAR,100,attr_valc,ierr)
        enddo
        write(*,*)' done defining 1d variables'
      endif

!-------------------------------------------------------------------
!  define 2d variables
!-------------------------------------------------------------------

      if(do_msf) then
        dim = 2
        dims(1) = latid(1)
        dims(2) = depid(1)
        do m = 1,num_2d_fields
           write(*,*)' m = ',m
           if(m == 1) then
              var_name = 'TMTG'
           elseif(m == 2) then
              var_name = 'TMTA'
           elseif(m == 3) then
              var_name = 'TMTIP'
           endif
           fieldid_2d(m) = NCVDEF(ncid,var_name,NCFLOAT,dim,dims,ierr)
           attr_name = 'units'
           attr_valc = 'Sv'
           call NCAPTC(ncid,fieldid_2d(m),attr_name
     &              ,NCCHAR,100,attr_valc,ierr)
           attr_name = 'missing_value'
           call NCAPT(ncid,fieldid_2d(m),attr_name
     &              ,NCFLOAT,1,special_value,ierr)
        enddo
        write(*,*)' done defining 2d variables'
      endif

      call NCENDF(ncid,ierr)

!-------------------------------------------------------------------
!  start writing the file
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  write the grid
!-------------------------------------------------------------------

      start(1) = 1
      count(1) = km
      call NCVPT(ncid,depvid(1),start,count,depth,ierr)
      start(1) = 1
      count(1) = ny_mht
      call NCVPT(ncid,latvid(1),start,count,lat_mht,ierr)

!-------------------------------------------------------------------
!  write 1d variables
!-------------------------------------------------------------------

      if(do_mht) then

        start(1) = 1
        count(1) = ny_mht

        do m = 1,num_1d_fields
          mht_temp = mht(:,m)
          call NCVPT(ncid,fieldid_1d(m),start,count,mht_temp,ierr)
          write(*,*)'...done writing 1D field #',m
        enddo

      endif

!-------------------------------------------------------------------
!  write 2d variables
!-------------------------------------------------------------------

      if(do_msf) then

        start(1) = 1
        start(2) = 1
        count(1) = ny_mht
        count(2) = km

        do m = 1,num_2d_fields
          do k = 1, km
            PSI_TEMP(:,k) = PSI(:,k,m)
            where (PSI_TEMP(:,k) == 0.) PSI_TEMP(:,k) = special_value
          enddo
          call NCVPT(ncid,fieldid_2d(m),start,count,PSI_TEMP,ierr)
          write(*,*)'...done writing 2D field #',m
        enddo

      endif

      call NCCLOS(ncid,ierr)

      !10 print *, 'error reading main namelist, iostat = ', status
      !26 stop 'error reading main namelist'
  999 continue
      end

!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
      subroutine grid_stuff(grid_file, ULAT, ULONG, DXU, DYU, TAREA
     & , imt, jmt)

      implicit none

      integer :: imt, jmt, i, j, rec_length

      real, dimension(imt,jmt) :: DXU, DYU, TAREA
      real, allocatable, dimension(:,:) :: DXT, DYT

      double precision, dimension(imt,jmt) :: ULAT, ULONG
      double precision, allocatable, dimension(:,:) :: HTN, HTE

      character (len=100) grid_file

!---------------------------------------------------------------
!  read in grid and create grid-related quantities need for
!   future calculations.
!---------------------------------------------------------------

      allocate ( DXT(imt,jmt), DYT(imt,jmt)
     &         , HTN(imt,jmt), HTE(imt,jmt) )

      inquire (iolength=rec_length) HTN
      write(*,*)'opening file: ',grid_file
      open(1,file=grid_file,access='direct',form='unformatted',
     &          recl=rec_length,status='unknown')
      read(1,rec=1)ULAT
      read(1,rec=2)ULONG
      read(1,rec=3)HTN
      read(1,rec=4)HTE
      close(1)

      do j = 2, jmt - 1
        do i = 2, imt - 1
          DXU(i,j) = 0.5*(HTN(i,j) + HTN(i+1,j))
          DYU(i,j) = 0.5*(HTE(i,j) + HTE(i,j+1))
          DXT(i,j) = 0.5*(HTN(i,j) + HTN(i,j-1))
          DYT(i,j) = 0.5*(HTE(i,j) + HTE(i-1,j))
        enddo
      enddo

      do j = 2, jmt - 1
        i = 1
        DXU(i,j) = 0.5*(HTN(i,j) + HTN(i+1,j))
        DYU(i,j) = 0.5*(HTE(i,j) + HTE(i,j+1))
        DXT(i,j) = 0.5*(HTN(i,j) + HTN(i,j-1))
        DYT(i,j) = 0.5*(HTE(i,j) + HTE(imt,j))
        i = imt
        DXU(i,j) = 0.5*(HTN(i,j) + HTN(1  ,j))
        DYU(i,j) = 0.5*(HTE(i,j) + HTE(i,j+1))
        DXT(i,j) = 0.5*(HTN(i,j) + HTN(i,j-1))
        DYT(i,j) = 0.5*(HTE(i,j) + HTE(i-1,j))
      enddo

      do i = 2, imt - 1
        j = 1
        DXU(i,j) = 0.5*(HTN(i,j) + HTN(i+1,j))
        DYU(i,j) = 0.5*(HTE(i,j) + HTE(i,j+1))
        DXT(i,j) = 0.5*(HTN(i,j) + HTN(i,1  ))
        DYT(i,j) = 0.5*(HTE(i,j) + HTE(i-1,j))
        j = jmt
        DXU(i,j) = 0.5*(HTN(i,j) + HTN(i+1,j))
        DYU(i,j) = 0.5*(HTE(i,j) + HTE(i,jmt))
        DXT(i,j) = 0.5*(HTN(i,j) + HTN(i,j-1))
        DYT(i,j) = 0.5*(HTE(i,j) + HTE(i-1,j))
      enddo

      i = 1
      j = 1
      DXU(i,j) = 0.5*(HTN(i,j) + HTN(i+1,j))
      DYU(i,j) = 0.5*(HTE(i,j) + HTE(i,j+1))
      DXT(i,j) = 0.5*(HTN(i,j) + HTN(i,1  ))
      DYT(i,j) = 0.5*(HTE(i,j) + HTE(imt,j))

      i = 1
      j = jmt
      DXU(i,j) = 0.5*(HTN(i,j) + HTN(i+1,j))
      DYU(i,j) = 0.5*(HTE(i,j) + HTE(i,jmt))
      DXT(i,j) = 0.5*(HTN(i,j) + HTN(i,j-1))
      DYT(i,j) = 0.5*(HTE(i,j) + HTE(imt,j))

      i = imt
      j = 1
      DXU(i,j) = 0.5*(HTN(i,j) + HTN(1  ,j))
      DYU(i,j) = 0.5*(HTE(i,j) + HTE(i,j+1))
      DXT(i,j) = 0.5*(HTN(i,j) + HTN(i,1  ))
      DYT(i,j) = 0.5*(HTE(i,j) + HTE(i-1,j))

      i = imt
      j = jmt
      DXU(i,j) = 0.5*(HTN(i,j) + HTN(1  ,j))
      DYU(i,j) = 0.5*(HTE(i,j) + HTE(i,jmt))
      DXT(i,j) = 0.5*(HTN(i,j) + HTN(i,j-1))
      DYT(i,j) = 0.5*(HTE(i,j) + HTE(i-1,j))

      TAREA = DXT*DYT

      deallocate ( DXT, DYT, HTN, HTE )

      end

!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
      subroutine wcalc(W, U, V, DXU, DYU, TAREA, KMT, DZU, dz, 
     &  imt, jmt, km)

      implicit none

      integer :: imt, jmt, km, i, j, k, im1, jm1

      integer, dimension(imt,jmt) :: KMT

      real, dimension(imt,jmt) :: DXU, DYU, TAREA
      real, dimension(imt,jmt,km) :: U, V, W, DZU
      real, dimension(km) :: dz

      double precision :: fvn, fvs, fue, fuw, wtkb, c0=0., p5=0.5

!-------------------------------------------------------------------
!-------------------------------------------------------------------

      write(*,*)' starting to calculate W...'

      do j=1,jmt
      jm1 = j - 1
      if (j == 1) jm1 = 1
      do i=1,imt
      im1 = i - 1
      if (i == 1) im1 = imt

         wtkb = c0  ! vertical velocity zero at bottom of bottom cell

         do k=KMT(i,j),1,-1

            !***
            !*** advection fluxes
            !***

            fue = p5*(U(i  ,j  ,k)*DYU(i  ,j    )*      
     &                               DZU(i  ,j  ,k) +     
     &                U(i  ,jm1,k)*DYU(i  ,jm1  )*      
     &                               DZU(i  ,jm1,k))
            fuw = p5*(U(im1,j  ,k)*DYU(im1,j    )*      
     &                               DZU(im1,j  ,k) +     
     &                U(im1,jm1,k)*DYU(im1,jm1  )*      
     &                               DZU(im1,jm1,k))
            fvn = p5*(V(i  ,j  ,k)*DXU(i  ,j    )*      
     &                               DZU(i  ,j  ,k) +     
     &                V(im1,j  ,k)*DXU(im1,j    )*      
     &                               DZU(im1,j  ,k))
            fvs = p5*(V(i  ,jm1,k)*DXU(i  ,jm1  )*      
     &                               DZU(i  ,jm1,k) +     
     &                V(im1,jm1,k)*DXU(im1,jm1  )*      
     &                               DZU(im1,jm1,k))

            !***
            !*** vertical velocity at top of box from continuity eq.
            !***

            W(i,j,k) = wtkb -  
     &              (fvn - fvs + fue - fuw)/TAREA(i,j)

            wtkb = W(i,j,k) ! top value becomes bottom for next pass

         enddo

      end do ! horizontal loops
      end do

      write(*,*)'...done'

      end  ! wcalc

!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
      subroutine calc_tpoints(TLONG,TLAT, ULONG, ULAT, imt, jmt, pi)

      implicit none

      integer :: imt, jmt

      real, dimension(imt,jmt) :: TLONG,TLAT
      real, allocatable, dimension(:,:) :: AT0,ATS,ATW,ATSW,WORKX
      real :: pi
      real dlong, left_long, right_long, cut_long

      double precision, dimension(imt,jmt) :: ULONG,ULAT

      logical all_pos_longs
      logical, allocatable, dimension(:,:) :: MASK

!---------------------------------------------------------------
!  The following is just a simple average of neighboring points
!---------------------------------------------------------------

      allocate ( MASK(imt,jmt), AT0(imt,jmt), ATS(imt,jmt)
     &, ATW(imt,jmt), ATSW(imt,jmt), WORKX(imt,jmt) )

      AT0  = 0.25
      ATS  = 0.25
      ATW  = 0.25
      ATSW = 0.25

      all_pos_longs = .true.
      if(ANY(ULONG .lt. 0.)) all_pos_longs = .false.
c     dlong = pi/36.  !  5  degrees
      dlong = pi/9.   !  20 degrees

      if(all_pos_longs) then
         left_long  = 2.*pi - dlong
         right_long = dlong
         cut_long   = 2.*pi
      else
         left_long  = pi - dlong
         right_long = pi + dlong
         cut_long   = pi
      endif

      MASK = .false.
      where( (ULONG .gt. left_long) .or. (ULONG .lt. right_long))
     &   MASK = .true.
      WORKX = ULONG
      call sw_4pt(TLAT ,AT0,ATS,ATW,ATSW,WORKX,imt,jmt)  !  use as temp
      right_long = right_long + dlong
      where( ULONG .lt. right_long ) WORKX = WORKX + 2.*pi

      call sw_4pt(TLONG,AT0,ATS,ATW,ATSW,WORKX,imt,jmt)
      where(.not. MASK) TLONG = TLAT
      WORKX = ULAT
      call sw_4pt(TLAT ,AT0,ATS,ATW,ATSW,WORKX,imt,jmt)

      where(TLONG .gt. cut_long) TLONG = TLONG - 2.*pi

      deallocate ( MASK,AT0,ATS,ATW,ATSW,WORKX )
      end

!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------

      subroutine sw_4pt(XOUT,CC,CS,CW,CSW,X,imt,jmt)

      implicit none

      integer i,j,imt,jmt

      real, dimension(imt,jmt) :: X,XOUT,CC,CS,CW,CSW

      do j=2,jmt
        do i=2,imt
          XOUT(i,j) = CC (i,j)*X(i,j) +
     &                CS (i,j)*X(i  ,j-1) + CW(i,j)*X(i-1,j) +
     &                CSW(i,j)*X(i-1,j-1)
        end do
      end do

      j=1
        do i=2,imt
          XOUT(i,j) = XOUT(i,j+1)
        end do

      i = 1
c     do j=2,jmt
c         XOUT(i,j) = XOUT(i+1,j)
c       end do
c  periodic
      do j=2,jmt
          XOUT(i,j) = CC (i,j)*X(i,j) +
     &                CS (i,j)*X(i  ,j-1) + CW(i,j)*X(imt,j) +
     &                CSW(i,j)*X(imt,j-1)
      end do
       XOUT(1,1) = XOUT(2,2)

      end  ! sw_4pt

!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------

      subroutine moc(V,W,KMTB,TLAT,DXU,TAREA,DZU,dz
     &              ,imt,jmt,km,lats,ny,PSI)

!---------------------------------------------------------------------
!     calculate meridional overturning streamfunction
!---------------------------------------------------------------------

      implicit none

      integer imt,jmt,km,j,k,jj,ny

      integer, dimension(imt,jmt) :: KMTB

      real :: southern_lat

      real, dimension(imt,jmt) :: TLAT, DXU, TAREA
      real, allocatable, dimension(:,:) :: WORK

      real, dimension(imt,jmt,km) :: W, V, DZU
      real, allocatable, dimension(:,:,:) :: WORK3D

      real, dimension(km,ny) :: PSI
      real, allocatable, dimension(:,:) :: AYZ
      real, allocatable, dimension(:) :: psi0

      real, dimension(ny) :: lats
      real, dimension(km) :: dz

      logical, allocatable, dimension(:,:) :: SMASK

!---------------------------------------------------------------------
!---------------------------------------------------------------------

      allocate (WORK(imt,jmt), SMASK(imt,jmt), AYZ(km,ny)
     &        , psi0(km), WORK3D(imt,jmt,km))

!---------------------------------------------------------------------
!     ...find j index of southernmost ocean point in basin
!---------------------------------------------------------------------

      jj = 0
      do j=1,jmt
        if (ANY(KMTB(:,j).ne.0)) then
          jj = j
          go to 888
        endif
      enddo
888   print *, " "
      print *, " southernmost j = ",jj
      southern_lat = 0.5*(TLAT(1,jj) + TLAT(1,jj-1)) ! ydeg(jj-1)
      print *, " southernmost lat = ",southern_lat

!---------------------------------------------------------------------
!     ...zonal and vertical integration to find streamfunction
!        across southernmost grid circle in basin
!---------------------------------------------------------------------

      SMASK = KMTB == 0 .and. cshift(KMTB,shift=+1,dim=2) > 0

      WORK = -DZU(:,:,km)*V(:,:,km)*DXU
      psi0(km) = sum(WORK(:,jj-1),mask=SMASK(:,jj-1))
      do k = km,2,-1
        WORK = -DZU(:,:,k-1)*V(:,:,k-1)*DXU
        psi0(k-1) = psi0(k)
     &         + sum(WORK(:,jj-1),mask=SMASK(:,jj-1))
      enddo

      print *, " psi0:"
      do k = 1,km
        print *, k,psi0(k)*1.e-12
        WORK3D(:,:,k) = W(:,:,k)*TAREA
      enddo

      do j = 1,ny
        if(mod(j,5) == 0) write(*,*)' merid psi, j = ',j

        do k = 1,km
!         WORK = W(:,:,k)*TAREA
          SMASK = (TLAT < lats(j) .and. k <= KMTB)

!         PSI(k,j) = sum(WORK,mask=SMASK)
          PSI(k,j) = sum(WORK3D(:,:,k),mask=SMASK)
!         AYZ(k,j) = sum(TAREA,mask=SMASK .and.
!    &      (   (k.le.cshift(KMTB,shift=+1,dim=1)
!    &           .and.(.not.cshift(SMASK,shift=+1,dim=1)))
!    &      .or.(k.le.cshift(KMTB,shift=-1,dim=1)
!    &           .and.(.not.cshift(SMASK,shift=-1,dim=1)))
!    &      .or.(k.le.cshift(KMTB,shift=+1,dim=2)
!    &           .and.(.not.cshift(SMASK,shift=+1,dim=2)))
!    &      .or.(k.le.cshift(KMTB,shift=-1,dim=2)
!    &           .and.(.not.cshift(SMASK,shift=-1,dim=2)))))
        enddo

        if(lats(j) >= southern_lat) PSI(:,j) = PSI(:,j) + psi0

      enddo

!---------------------------------------------------------------------
!     ... smooth grid-point noise over y
!---------------------------------------------------------------------

      do j = 2,ny-1
        PSI(:,j) = 0.25*(PSI(:,j-1)+PSI(:,j+1)) + 0.5*PSI(:,j)
      enddo

      PSI = PSI*1.0e-12  ! normalize to Sv

      deallocate (WORK, SMASK, AYZ, psi0, WORK3D)

      end
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------

      subroutine nml_defaults (
     &   imt, jmt, km,ysouth_mht, ynorth_mht, dy_mht
     &  ,kmt_file, kmt_atl_file, kmt_indopac_file, in_depths, grid_file
     &  ,u_file, v_file , u_first_record, v_first_record
     &  ,uet_file, vnt_file, uet_first_record, vnt_first_record
     &  , out_file, verbose, do_msf, do_mht, pbcfile)

      implicit none

      integer :: imt, jmt, km, u_first_record, v_first_record
     &  , uet_first_record, vnt_first_record

      real :: ysouth_mht, ynorth_mht, dy_mht

      character (len=100) :: kmt_file, kmt_atl_file, kmt_indopac_file
     &, in_depths, grid_file, u_file, v_file, out_file
     &, uet_file, vnt_file, pbcfile
     
      logical :: verbose, do_msf, do_mht

!---------------------------------------------------------------------
!  set default values for namelist variables
!---------------------------------------------------------------------

      imt = 100
      jmt = 116
      km = 25
      ysouth_mht = -75.
      ynorth_mht = 75.
      dy_mht = 1.
      kmt_file = 'unknown'
      pbcfile = 'unknown'
      kmt_atl_file = 'unknown'
      kmt_indopac_file = 'unknown'
      in_depths = 'unknown'
      grid_file = 'unknown'
      u_file = 'unknown'
      v_file = 'unknown'
      u_first_record  = 1
      v_first_record = 1
      uet_file = 'unknown'
      vnt_file = 'unknown'
      uet_first_record  = 1
      vnt_first_record = 1
      out_file = 'unknown'
      verbose = .true.
      do_msf = .false.
      do_mht = .false.

      end
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------

      subroutine meridional_heat(UET,VNT,KMTB,TLAT,TAREA,dz
     &  ,imt,jmt,km,lats,ny,mht)

!---------------------------------------------------------------------
!     calculate mean meridional heat transport
!---------------------------------------------------------------------

      implicit none

      integer :: imt,jmt,km,n,i,j,jj,ny,k

      integer, dimension(imt,jmt) :: KMTB

      real :: y0,dy,y,mht0,southern_lat

      real, dimension(km) :: dz
      real, dimension(ny) :: mht, lats
      real, dimension(imt,jmt,km) :: UET,VNT
      real, dimension(imt,jmt) :: TLAT, TAREA
      real, allocatable, dimension(:,:) :: WORK1,WORK2,WORKX,WORKY

      logical, allocatable, dimension(:,:) :: SMASK

!---------------------------------------------------------------------
!---------------------------------------------------------------------

      allocate (WORK1(imt,jmt),WORK2(imt,jmt),WORKX(imt,jmt)
     &         ,WORKY(imt,jmt),SMASK(imt,jmt))

!---------------------------------------------------------------------
!     ...find vertically-integrated heat transport
!---------------------------------------------------------------------

      WORKX = 0.
      WORKY = 0.
      do k=1,km
        WORKX = WORKX + UET(:,:,k)*TAREA*dz(k)*4.186e-15 ! pw
        WORKY = WORKY + VNT(:,:,k)*TAREA*dz(k)*4.186e-15 ! pw
      enddo

!     ...find divergence of vertically-integrated heat transport

!     call w_shift(WORK1,WORKX)
!     call s_shift(WORK2,WORKY)
      WORK1 = cshift(WORKX, shift=-1, dim=1)
      WORK2 = cshift(WORKY, shift=-1, dim=2)
      WORK1 = WORKX - WORK1 + WORKY - WORK2

!     ...find j index of southernmost ocean point in basin

      jj = 0
      do j=1,jmt
        if (ANY(KMTB(:,j) /= 0)) then
          jj = j
          go to 888
        endif
      enddo
888   print *, " "
      print *, " southernmost j = ",jj
      southern_lat = 0.5*(TLAT(1,jj) + TLAT(1,jj-1)) ! ydeg(jj-1)
      print *, " southernmost lat = ",southern_lat

!     ...zonal and  integration to find  heat transport
!        across southernmost grid circle in basin

      SMASK = ( KMTB == 0 .and. cshift(KMTB,shift=+1,dim=2) > 0 )

      mht0 = sum(WORKY(:,jj-1),mask=SMASK(:,jj-1))

      print *, " mht0:", mht0

!     ...loop over latitudes

      do j = 1,ny
        SMASK = (TLAT < lats(j) .and. KMTB > 0)
        mht(j) = sum(WORK1,mask=SMASK)
        if(lats(j) > southern_lat)  mht(j) = mht(j) + mht0
      enddo

      deallocate (WORK1,WORK2,WORKX,WORKY,SMASK)

      end  ! meridional_heat

