subroutine interp_ll(ib,jb,kb,x1,y1,z1,ns)

!     computes latitude, longitude, depth at position of trajectory
!     by interpolating data from the center of eight nearest boxes

USE mod_grid
USE mod_vel
IMPLICIT none

REAL*8  :: x1,y1,z1,ax,ay,az

REAL    :: tpp,tpm,tmp,tmm
REAL    :: spp,spm,smp,smm
REAL    :: rppp,rppm,rpmp,rpmm,rmpp,rmpm,rmmp,rmmm

INTEGER :: ib,jb,kb,ip,im,jp,jm,kp,kn,ns
! determining nearest centers of boxes
      if(x1.le.dble(ib)-dble(.5)) then
       ip=ib
       im=ib-1
       if(im.eq.0) im=imt
      else
       ip=ib+1
       im=ib
       if(ip.gt.imt) ip=1
      endif
      if(y1.le.dble(jb)-dble(.5)) then
       jp=jb
       jm=jb-1
       if(jm.eq.0) jm=1
      else
       jp=jb+1
       jm=jb
       if(jp.gt.jmt) jp=jmt
      endif

      if(z1.le.dble(kb)-dble(.5)) then
       kp=kb
       kn=kb-1
       if(kn.le.0) kn=1
      else
       kp=kb+1
       if(kp.gt.km) kp=km
       kn=kb
      endif

      ax=dble(ip)-x1
      if(ax.gt.100.d0) then
!       print *,ax,ip,im,x1,ib
       ax=ax-dble(imt)
      elseif(ax.lt.-100.d0) then
       ax=ax+dble(imt)
      endif
      ay=(dble(jp)-y1)
      az=(dble(kp)-z1)

      tpp=lon_rho(ip,jp)
      spp=lat_rho(ip,jp)
      rppp=z_r(ip,jp,kp,ns)
      if(tpp==0. .and. spp==0.) then
       tpp=lon_rho(ip,jp)
       spp=lat_rho(ip,jp)
       rppp=z_r(ip,jp,kn,ns)
      endif

      rppm=z_r(ip,jp,kn,ns)
      if(rppm==0.) then
       rppm=z_r(ip,jp,kn,ns)
      endif

      tpm=lon_rho(ip,jm)
      spm=lat_rho(ip,jm)
      rpmp=z_r(ip,jm,kp,ns)
      if(tpm==0. .and. spm==0.) then
       tpm=lon_rho(ip,jp)
       spm=lat_rho(ip,jp)
       rpmp=z_r(ip,jp,kn,ns)
      endif

      rpmm=z_r(ip,jm,kn,ns)
      if(rpmm==0.) then
       rpmm=z_r(ip,jp,kn,ns)
      endif

      tmp=lon_rho(im,jp)
      smp=lat_rho(im,jp)
      rmpp=z_r(im,jp,kp,ns)
      if(tmp==0. .and. smp==0.) then
       tmp=lon_rho(ip,jp)
       smp=lat_rho(ip,jp)
       rmpp=z_r(ip,jp,kn,ns)
      endif

      rmpm=z_r(im,jp,kn,ns)
      if(rmpm==0.) then
       rmpm=z_r(ip,jp,kn,ns)
      endif

      tmm=lon_rho(im,jm)
      smm=lat_rho(im,jm)
      rmmp=z_r(im,jm,kp,ns)
      if(tmm==0. .and. smm==0.) then
       tmm=lon_rho(ip,jp)
       smm=lat_rho(ip,jp)
       rmmp=z_r(ip,jp,kn,ns)
      endif

      rmmm=z_r(im,jm,kn,ns)
      if(rmmm==0.) then
       rmmm=z_r(ip,jp,kn,ns)
      endif



      lon=tpp*(1.-ax)*(1.-ay) &
        + tmp*    ax *(1.-ay) &
        + tpm*(1.-ax)*    ay  &
        + tmm*    ax *    ay

      lat=spp*(1.-ax)*(1.-ay) &
        + smp*    ax *(1.-ay) &
        + spm*(1.-ax)*    ay  &
        + smm*    ax *    ay

      zed=rppp*(1.-ax)*(1.-ay)*(1.-az) &
        + rmpp*    ax *(1.-ay)*(1.-az) &
        + rpmp*(1.-ax)*    ay *(1.-az) &
        + rmmp*    ax *    ay *(1.-az) &
        + rppm*(1.-ax)*(1.-ay)*    az  &
        + rmpm*    ax *(1.-ay)*    az  &
        + rpmm*(1.-ax)*    ay *    az  &
        + rmmm*    ax *    ay *    az

return
end subroutine interp_ll
