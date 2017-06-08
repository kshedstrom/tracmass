SUBROUTINE readfields

  USE netcdf
  USE mod_param
  USE mod_vel
!  USE mod_coord   !FC
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile
  USE mod_seed, only: nff! LD ADDED, for nff
#ifdef tempsalt
  USE mod_tempsalt
  USE mod_dens
#endif

  IMPLICIT none
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  ! = Variables for filename generation
  CHARACTER                                  :: dates(62)*17
  CHARACTER (len=200)                        :: dataprefix, dstamp
  INTEGER                                    :: intpart1 ,intpart2
  INTEGER                                    :: ndates
  INTEGER                                    :: yr1 ,mn1 ,dy1,hr
  INTEGER                                    :: yr2 ,mn2 ,dy2

  ! = Loop variables
  INTEGER                                    :: t ,i ,j ,k ,kk ,tpos

  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile
  CHARACTER (len=50)                         :: varName

  ! = Variables for converting from S to Z
  REAL*8,       ALLOCATABLE, DIMENSION(:)    :: sc_r,Cs_r
  REAL*8,       ALLOCATABLE, DIMENSION(:)    :: sc_w,Cs_w
  INTEGER                                    :: hc
#ifdef arctic_cod
  REAL*8,       ALLOCATABLE, DIMENSION(:,:,:):: u_east,v_north
#endif

  ! = Input fields from GCM
  REAL*8,       ALLOCATABLE, DIMENSION(:,:)    :: ssh,dzt0
  ! ===   ===   ===


  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( ssh(imt,jmt), dzt0(imt,jmt) )
     allocate ( sc_r(km), Cs_r(km) )
     allocate ( sc_w(km), Cs_w(km) )
  end if alloCondUVW
  alloCondDZ: if(.not. allocated (dzu)) then
     allocate ( dzu(imt,jmt,km), dzv(imt,jmt,km) )
  end if alloCondDZ
#ifdef arctic_cod
  if(.not. allocated (u_east)) then
     allocate ( u_east(imt,jmt,km), v_north(imt,jmt,km) )
  end if
#endif
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  sc_r = 0
  Cs_r = 0
  sc_w = 0
  Cs_w = 0


  call datasetswap
  call updateClock

  ! === update the time counting ===
  intpart1    = mod(ints,24)
  intpart2    = floor((ints)/24.)
  dstamp      = 'arctic2_avg_XXXX-XX-XXT00:00:00.nc'

  write (dstamp(13:16),'(i4.2)') currYear
  write (dstamp(18:19),'(i2.2)') currMon
  write (dstamp(21:22),'(i2.2)') currDay

  dataprefix  = trim(inDataDir) // dstamp
  tpos        = intpart1+1
  print *, "curr JD:", currJDtot, "read file:", dataprefix

#ifdef arctic_cod
  u_east    = get3DfieldNC(trim(dataprefix) ,  'u_eastward')
  v_north   = get3DfieldNC(trim(dataprefix) , 'v_northward')
! should be on cell centers here.
  do k=1,km
    uvel(:imt,:,k) = u_east(:,:,k)*cos(ang(:imt,:)) + &
                     v_north(:,:,k)*sin(ang(:imt,:))
    vvel(:imt,:,k) = -u_east(:,:,k)*sin(ang(:imt,:)) +&
                     v_north(:,:,k)*cos(ang(:imt,:))
  enddo
! reusing arrays
  u_east(1:imt-1,:,:)  = 0.5*(uvel(1:imt-1,:,:) + uvel(2:imt,:,:))
  v_north(:imt,1:jmt-1,:) = 0.5*(vvel(:imt,1:jmt-1,:) + vvel(:imt,2:jmt,:))
  uvel(:imt,:,:) = u_east
  vvel(:imt,:,:) = v_north
#else
  uvel      = get3DfieldNC(trim(dataprefix) ,   'u')
  vvel      = get3DfieldNC(trim(dataprefix) ,   'v')
#endif
  ssh       = get2DfieldNC(trim(dataprefix) ,'zeta')
#ifdef explicit_w
  wvel      = get3DfieldNC(trim(dataprefix) ,'omega')
#endif
  where (uvel > 1000)
     uvel = 0
  end where
  where (vvel > 1000)
     vvel = 0
  end where
  where (ssh > 1000)
     ssh = 0
  end where

  hs(:imt,:jmt,2) = ssh(:imt,:jmt)

#ifdef explicit_w
  wflux(:,:,:,2) = 0.
  do j=1,jmt
    do i=1,imt
      wflux(i,j,0:km-1,2) = wvel(i,j,1:km)*dxdy(i,j)
    end do
  end do
#endif

  Cs_w = get1DfieldNC (trim(dataprefix), 'Cs_w')
  sc_w = get1DfieldNC (trim(dataprefix), 's_w')
  Cs_r = get1DfieldNC (trim(dataprefix), 'Cs_r')
  sc_r = get1DfieldNC (trim(dataprefix), 's_rho')
  hc   = getScalarNC (trim(dataprefix), 'hc')

  dzt(:,:,:,1)=dzt(:,:,:,2)

  z_w(:,:,0,2) = depth(:imt,:)

  do k=1,km
    dzt0 = (hc*sc_r(k) + depth*Cs_r(k)) / (hc + depth)
    z_r(:,:,k,2) = ssh(:imt,:) + (ssh(:imt,:) + depth(:imt,:)) * dzt0(:imt,:)
    dzt(:,:,k,2) = ssh(:imt,:) + (ssh(:imt,:) + depth(:imt,:)) * dzt0(:imt,:)
  enddo
  dzt0(:imt,:) = dzt(:,:,km,2)
  dzt(:,:,1:km-1,2)=dzt(:,:,2:km,2)-dzt(:,:,1:km-1,2)
  dzt(:,:,km,2) = ssh(:imt,:) - dzt0(:imt,:)

  dzu(1:imt-1,:,:) = dzt(1:imt-1,:,:,2)*0.5 + dzt(2:imt,:,:,2)*0.5
  dzv(:,1:jmt-1,:) = dzt(:,1:jmt-1,:,2)*0.5 + dzt(:,2:jmt,:,2)*0.5

  do k=1,km
     uflux(:,:,k,2)    = uvel(:imt,:,k) * dzu(:,:,k) * dyu(:imt,:)
     vflux(:,1:jmt,k,2)   = vvel(:imt,:,k) * dzv(:,:,k) * dxv(:imt,:)
  end do

  if (nff .le. 0) then
     uflux = -uflux
     vflux = -vflux
  end if

#ifdef tempsalt
  tem(:,:,:,2)      = get3DfieldNC(trim(dataprefix) ,   'temp')
  sal(:,:,:,2)      = get3DfieldNC(trim(dataprefix) ,   'salt')
  !rho(:,:,:,2)      = get3DfieldNC(trim(dataprefix) ,   'rho')
#ifdef larval_fish
  srflux(:,:,2)     = get2DfieldNC(trim(dataprefix) ,   'swrad')
! Note: this works as long as surface AKt is zero.
  ak2(:,:,:)        = get3DfieldNC(trim(dataprefix) ,   'AKt')
  akt(:,:,0:km-1,2) = ak2(:,:,:)
#endif
#endif

  return

end subroutine readfields
