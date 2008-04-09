!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x

subroutine coordinat

!USE param_variables
USE mod_param
USE mod_coord
USE mod_grid
USE mod_name
IMPLICIT none

INTEGER i,j,k,kk,mois(12)
REAL*8 pi,radian,radius,rlatt,rlatu,rmax,smax,tmax
REAL :: a

#if defined occam25 || occ66
REAL*8 b1,b2,b3,b4,z
#endif

#if defined occam25
data b1,b2,b3,b4/ 13635.703538701,12023.336474772,15.503284687122,7.2929323809417 /
#endif
#if defined occ66
data b1,b2,b3,b4/10725.0d0,  10275.0d0, 36.0d0,     13.0d0/
#endif

data mois/31,28,31,30,31,30,31,31,30,31,30,31/

!_________________________________ time lengths __________________________________________
tyear=365.25d0 * 24.d0 * 3600.d0
tday=24.d0 * 3600.d0

! month lengths
do i=1000,3000
do k=1,12
idmax(k,i)=mois(k)
enddo
if((mod(i,4).eq.0 .and. mod(i,100).ne.0) .or.  mod(i,400).eq.0 ) idmax(2,i)=29
!print *,i,(idmax(k,i),k=1,12)
enddo

! handling quantities concerning coordinates
radius = 6371229.d0 ! earth radius
pi = 2.d0 * dasin(1.d0)
radian = pi/180.d0
deg=radius*radian

! x,y gridsize in degrees

rmin=19.d0
rmax=28.5d0

#if defined occam25 || occ66
dx = 0.25d0
dy = 0.25d0
rmin=19.d0
rmax=28.5d0
#ifdef mod1
stlon1 = dx
stlat1 = -78.d0
#endif
#ifdef mod2
stlon2 = dx
stlat2 = -54.50d0
#endif
#endif

#ifdef rco
dx = 1./15.d0
dy = 1./30.d0
!stlon1= 9.0d0+0.25d0*dx
!stlat1=53.8d0+1.25d0*dy
!print *,'coord stlon,stlat=',stlon1,stlat1
rmin=-2.d0
rmax=10.d0
#endif

#if defined for || sim
dx = 185.2d0/deg
dy = dx
#ifdef sim
stlon1=16.d0+31.d0/60.d0
stlat1=57.d0+20.d0/60.d0
#endif
#ifdef for
stlon1=17.d0+59.d0/60.d0
stlat1=60.d0+ 8.d0/60.d0
#endif
rmin=-2.d0
rmax=10.d0
#endif

#ifdef tes
dx = 0.5d0
dy = 0.5d0
stlon1=0.
stlat1=0.
rmin=20.d0
rmax=30.d0
#endif

#ifdef tun
dx = 1./32.d0
dy = 1./32.d0
stlon1=9.72d0-0.5d0*dx 
stlat1=32.5d0-0.5d0*dy
rmin=20.d0
rmax=30.d0
#endif

dr=(rmax-rmin)/dble(MR-1)

#ifdef streamts 
tmin=-2.d0
! Values for the Baltic 
#if defined rco || for || sim  
tmax=25.d0 
smin= 0.d0
smax=15.d0
! Values for the world ocean but bad for brakish water
#else               
tmax=30.d0 
smin=20.d0
smax=40.d0
#endif
dtemp=(tmax-tmin)/dble(MR-1)
dsalt=(smax-smin)/dble(MR-1)
#endif

! cosines relating to corners of grid box: ...u
! cosines relating to center of grid box:   ...t
do j=1,JMT
 rlatu=stlat1+dy*(j-1)+dy
 rlatt=rlatu-0.5d0*dy
 csu(j)=dcos(rlatu*radian)
 cst(j)=dcos(rlatt*radian)
#if defined for || sim 
 csu(j)=1.d0
 cst(j)=1.d0
#endif
enddo

! z coordinates (z>0 going up) for layers in meters 
! bottom layer: k=0; surface layer: k=KM and zw=0
! dz = layer thickness  

#ifdef rco
dz(41)=   3.d0
dz(40)=   3.d0
dz(39)=   3.d0
dz(38)=   3.d0
dz(37)=   3.d0
dz(36)=   3.d0
dz(35)=   3.d0
dz(34)=   3.d0
dz(33)=   3.d0
dz(32)=   3.d0
dz(31)=   3.d0
dz(30)=   3.d0
dz(29)=   3.d0
dz(28)=   3.007080d0
dz(27)=   3.063581d0
dz(26)=   3.175872d0
dz(25)=   3.342542d0
dz(24)=   3.561495d0
dz(23)=   3.829976d0
dz(22)=   4.144610d0
dz(21)=   4.501440d0
dz(20)=   4.895979d0
dz(19)=   5.323265d0
dz(18)=   5.777925d0
dz(17)=   6.254241d0
dz(16)=   6.746222d0
dz(15)=   7.247683d0
dz(14)=   7.752317d0
dz(13)=   8.253778d0
dz(12)=   8.745760d0
dz(11)=   9.222075d0
dz(10)=   9.676735d0
dz( 9)=   10.10402d0
dz( 8)=   10.49856d0
dz( 7)=   10.85539d0
dz( 6)=   11.17002d0
dz( 5)=   11.43851d0
dz( 4)=   11.65746d0
dz( 3)=   11.82413d0
dz( 2)=   11.93642d0
dz( 1)=   11.99292d0
#endif

#if defined for || sim 
 dz(39)=2.5d0
 do j=38,9,-1
  dz(j)=5.d0
 enddo
 do j=8,2,-1
  dz(j)=10.d0
 enddo
 dz(1)=20.d0
#endif

#if defined occam25 
do k=0,KM
 zw(k)=b1*dble(KM-k)+b2*b4*dlog(dcosh((dble(KM-k)-b3)/b4)/dcosh(-b3/b4))
 zw(k)=1.d-2*zw(k) 
enddo
zw(KM)=0.d0
do k= 1, KM
 dz(k)  = zw(k-1) - zw(k)
enddo
#endif

#if defined occ66
! depth coord in meters from surface and down (OCCAM def)
do k=0,KM
 z = k
 zw(k) = b1*z + b2*b4*dlog(dcosh((z-b3)/b4)/dcosh(-b3/b4))
 zw(k)=1.d-2*zw(k) 
enddo

! layer depths in meters from bottom and up (trajecory def)
do k= 1, KM
 kk=KM+1-k
 dz(kk)  = zw(k) - zw(k-1) 
enddo
!print *,zw
!print *,dz
#endif


#if defined occam25 || occ66 || rco || tes || for || sim 
do i=1,IMT
 do j=1,JMT
  dxdy(i,j)=dx*cst(j)*dy*deg**2
 enddo
enddo
#endif

#if defined ifs
! This is in hPA but should perhaps be in Pa i.e. to multiply with 100!!!!!!!!!!!!!
zw( 0)=0.
zw( 1)=5.
zw( 2)=30.
zw( 3)=75.
zw( 4)=150.
zw( 5)=250.
zw( 6)=350.
zw( 7)=450.
zw( 8)=600.
zw( 9)=775.
zw(10)=887.5
zw(11)=962.5 ! but this level is overriten to take account of the surface pressure

do k=1,KM
 dz(k)=zw(k)-zw(k-1)
! print *,'k=',k,' dz(k)=',dz(k)
enddo

dx = 1.125
dy = 1000000. ! irregular and should be dy(j)
stlon1=0. 
stlat1=-90000.

! geopotential height in km
rmin=0.d0
rmax=200.d0
dr=(rmax-rmin)/dble(MR-1)
! temperature
tmin=-100.d0
tmax=50.d0 
dtemp=(tmax-tmin)/dble(MR-1)
! specific humidity
smin= 0.d0
smax=25.d0
dsalt=(smax-smin)/dble(MR-1)

! read in the guassian latitude coordinates
print *,'topo=',directory
!print *,'topo=',directory//'/topo/ifsn80lat.txt'
open(12,file=directory//'/topo/ifsn80lat.txt')
99 format(50x,f9.5)
do j=1,JMT-1
read(12,99) a
phi(j)=-a
enddo
phi(0  )=-90. ! add the north pole as a fictive grid point
phi(JMT)= 90. ! add the south pole as a fictive grid point

do j=1,JMT
dyt(j)=(phi(j)-phi(j-1))*deg
!print *,j,phi(j),phi(j-1),dyt(j)
enddo

do j=1,JMT
 rlatt=0.5*(phi(j)+phi(j-1))
 rlatu=phi(j)
 csu(j)=dcos(rlatu*radian)
 cst(j)=dcos(rlatt*radian)
enddo

do i=1,IMT
 do j=1,JMT
  dxdy(i,j)=dx*deg*cst(j)*dyt(j)
 enddo
enddo

kmt=KM  ! all grid cells are active in an AGCM


#endif

return
end subroutine coordinat
