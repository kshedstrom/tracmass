#ifdef tempsalt

subroutine interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,ns)

!     computes temperature, salinity, density at position of trajectory
!     by interpolating data from the center of eight nearest boxes
!
!     This subroutine should be improved in order to include time interpolation

USE mod_grid
USE mod_dens
USE mod_vel
USE mod_tempsalt
USE mod_loopvars
IMPLICIT none

REAL*8  :: x1,y1,z1,ax,ay,az,az2

REAL    :: tppp,tppm,tpmp,tpmm,tmpp,tmpm,tmmp,tmmm
REAL    :: sppp,sppm,spmp,spmm,smpp,smpm,smmp,smmm
REAL    :: rppp,rppm,rpmp,rpmm,rmpp,rmpm,rmmp,rmmm
REAL    :: zppp,zppm,zpmp,zpmm,zmpp,zmpm,zmmp,zmmm
REAL    :: temp,salt,dens
REAL    :: srfpp,srfpm,srfmp,srfmm
REAL    :: hsbpp,hsbpm,hsbmp,hsbmm

INTEGER :: ib,jb,kb,ip,im,jp,jm,kp,kn,kn2,ns,kp2
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

      kp2=kb
      kn2=kb-1

      ax=dble(ip)-x1
      if(ax.gt.100.d0) then
!       print *,ax,ip,im,x1,ib
       ax=ax-dble(imt)
!       stop 49678
      elseif(ax.lt.-100.d0) then
       ax=ax+dble(imt)
!       stop 49679
      endif
      ay=(dble(jp)-y1)
      az=(dble(kp)-z1)
      az2=(dble(kp2)-z1)

! temperature, salinity, density calculation
      tppp=tem(ip,jp,kp,ns)*mask(ip,jp)
      sppp=sal(ip,jp,kp,ns)*mask(ip,jp)
      rppp=rho(ip,jp,kp,ns)*mask(ip,jp)
      if(tppp==0. .and. sppp==0.) then
       tppp=tem(ip,jp,kn,ns)
       sppp=sal(ip,jp,kn,ns)
       rppp=rho(ip,jp,kn,ns)
      endif

      tppm=tem(ip,jp,kn,ns)*mask(ip,jp)
      sppm=sal(ip,jp,kn,ns)*mask(ip,jp)
      rppm=rho(ip,jp,kn,ns)*mask(ip,jp)
      if(tppm==0. .and. sppm==0.) then
       tppm=tem(ip,jp,kn,ns)
       sppm=sal(ip,jp,kn,ns)
       rppm=rho(ip,jp,kn,ns)
      endif

      tpmp=tem(ip,jm,kp,ns)*mask(ip,jm)
      spmp=sal(ip,jm,kp,ns)*mask(ip,jm)
      rpmp=rho(ip,jm,kp,ns)*mask(ip,jm)
      if(tpmp==0. .and. spmp==0.) then
       tpmp=tem(ip,jp,kn,ns)
       spmp=sal(ip,jp,kn,ns)
       rpmp=rho(ip,jp,kn,ns)
      endif

      tpmm=tem(ip,jm,kn,ns)*mask(ip,jm)
      spmm=sal(ip,jm,kn,ns)*mask(ip,jm)
      rpmm=rho(ip,jm,kn,ns)*mask(ip,jm)
      if(tpmm==0. .and. spmm==0.) then
       tpmm=tem(ip,jp,kn,ns)
       spmm=sal(ip,jp,kn,ns)
       rpmm=rho(ip,jp,kn,ns)
      endif

      tmpp=tem(im,jp,kp,ns)*mask(im,jp)
      smpp=sal(im,jp,kp,ns)*mask(im,jp)
      rmpp=rho(im,jp,kp,ns)*mask(im,jp)
      if(tmpp==0. .and. smpp==0.) then
       tmpp=tem(ip,jp,kn,ns)
       smpp=sal(ip,jp,kn,ns)
       rmpp=rho(ip,jp,kn,ns)
      endif

      tmpm=tem(im,jp,kn,ns)*mask(im,jp)
      smpm=sal(im,jp,kn,ns)*mask(im,jp)
      rmpm=rho(im,jp,kn,ns)*mask(im,jp)
      if(tmpm==0. .and. smpm==0.) then
       tmpm=tem(ip,jp,kn,ns)
       smpm=sal(ip,jp,kn,ns)
       rmpm=rho(ip,jp,kn,ns)
      endif

      tmmp=tem(im,jm,kp,ns)*mask(im,jm)
      smmp=sal(im,jm,kp,ns)*mask(im,jm)
      rmmp=rho(im,jm,kp,ns)*mask(im,jm)
      if(tmmp==0. .and. smmp==0.) then
       tmmp=tem(ip,jp,kn,ns)
       smmp=sal(ip,jp,kn,ns)
       rmmp=rho(ip,jp,kn,ns)
      endif

      tmmm=tem(im,jm,kn,ns)*mask(im,jm)
      smmm=sal(im,jm,kn,ns)*mask(im,jm)
      rmmm=rho(im,jm,kn,ns)*mask(im,jm)
      if(tmmm==0. .and. smmm ==0.) then
       tmmm=tem(ip,jp,kn,ns)
       smmm=sal(ip,jp,kn,ns)
       rmmm=rho(ip,jp,kn,ns)
      endif

      temp=tppp*(1.-ax)*(1.-ay)*(1.-az) &
        + tmpp*    ax *(1.-ay)*(1.-az) &
        + tpmp*(1.-ax)*    ay *(1.-az) &
        + tmmp*    ax *    ay *(1.-az) &
        + tppm*(1.-ax)*(1.-ay)*    az  &
        + tmpm*    ax *(1.-ay)*    az  &
        + tpmm*(1.-ax)*    ay *    az  &
        + tmmm*    ax *    ay *    az

      salt=sppp*(1.-ax)*(1.-ay)*(1.-az) &
        + smpp*    ax *(1.-ay)*(1.-az) &
        + spmp*(1.-ax)*    ay *(1.-az) &
        + smmp*    ax *    ay *(1.-az) &
        + sppm*(1.-ax)*(1.-ay)*    az  &
        + smpm*    ax *(1.-ay)*    az  &
        + spmm*(1.-ax)*    ay *    az  &
        + smmm*    ax *    ay *    az

      dens=rppp*(1.-ax)*(1.-ay)*(1.-az) &
        + rmpp*    ax *(1.-ay)*(1.-az) &
        + rpmp*(1.-ax)*    ay *(1.-az) &
        + rmmp*    ax *    ay *(1.-az) &
        + rppm*(1.-ax)*(1.-ay)*    az  &
        + rmpm*    ax *(1.-ay)*    az  &
        + rpmm*(1.-ax)*    ay *    az  &
        + rmmm*    ax *    ay *    az

return
end subroutine interp
#endif
