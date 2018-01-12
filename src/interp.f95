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
REAL    :: mpp, mpm, mmp, mmm, maskden

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
      mpp = mask(ip,jp)
      tppp=tem(ip,jp,kp,ns)*mpp
      sppp=sal(ip,jp,kp,ns)*mpp
      rppp=rho(ip,jp,kp,ns)*mpp

      tppm=tem(ip,jp,kn,ns)*mpp
      sppm=sal(ip,jp,kn,ns)*mpp
      rppm=rho(ip,jp,kn,ns)*mpp

      mpm = mask(ip,jm)
      tpmp=tem(ip,jm,kp,ns)*mpm
      spmp=sal(ip,jm,kp,ns)*mpm
      rpmp=rho(ip,jm,kp,ns)*mpm

      tpmm=tem(ip,jm,kn,ns)*mpm
      spmm=sal(ip,jm,kn,ns)*mpm
      rpmm=rho(ip,jm,kn,ns)*mpm

      mmp = mask(im,jp)
      tmpp=tem(im,jp,kp,ns)*mmp
      smpp=sal(im,jp,kp,ns)*mmp
      rmpp=rho(im,jp,kp,ns)*mmp

      tmpm=tem(im,jp,kn,ns)*mmp
      smpm=sal(im,jp,kn,ns)*mmp
      rmpm=rho(im,jp,kn,ns)*mmp

      mmm = mask(im,jm)
      tmmp=tem(im,jm,kp,ns)*mmm
      smmp=sal(im,jm,kp,ns)*mmm
      rmmp=rho(im,jm,kp,ns)*mmm

      tmmm=tem(im,jm,kn,ns)*mmm
      smmm=sal(im,jm,kn,ns)*mmm
      rmmm=rho(im,jm,kn,ns)*mmm

      maskden = mpp*(1.-ax)*(1.-ay) &
              + mmp*    ax *(1.-ay) &
              + mpm*(1.-ax)*    ay  &
              + mmm*    ax *    ay

      temp=(tppp*(1.-ax)*(1.-ay)*(1.-az) &
          + tmpp*    ax *(1.-ay)*(1.-az) &
          + tpmp*(1.-ax)*    ay *(1.-az) &
          + tmmp*    ax *    ay *(1.-az) &
          + tppm*(1.-ax)*(1.-ay)*    az  &
          + tmpm*    ax *(1.-ay)*    az  &
          + tpmm*(1.-ax)*    ay *    az  &
          + tmmm*    ax *    ay *    az) &
          / maskden

      salt=(sppp*(1.-ax)*(1.-ay)*(1.-az) &
          + smpp*    ax *(1.-ay)*(1.-az) &
          + spmp*(1.-ax)*    ay *(1.-az) &
          + smmp*    ax *    ay *(1.-az) &
          + sppm*(1.-ax)*(1.-ay)*    az  &
          + smpm*    ax *(1.-ay)*    az  &
          + spmm*(1.-ax)*    ay *    az  &
          + smmm*    ax *    ay *    az) &
          / maskden

      dens=(rppp*(1.-ax)*(1.-ay)*(1.-az) &
          + rmpp*    ax *(1.-ay)*(1.-az) &
          + rpmp*(1.-ax)*    ay *(1.-az) &
          + rmmp*    ax *    ay *(1.-az) &
          + rppm*(1.-ax)*(1.-ay)*    az  &
          + rmpm*    ax *(1.-ay)*    az  &
          + rpmm*(1.-ax)*    ay *    az  &
          + rmmm*    ax *    ay *    az) &
          / maskden

return
end subroutine interp
#endif
