
PROGRAM main
USE mod_param
USE mod_name
USE mod_time 
USE mod_domain
USE mod_buoyancy

IMPLICIT none

INTEGER i,j,n

!____________________________________________________________________________
print *,'TRACMASS trajectory code starts at'
call system('date')
!___________________________ print preprocess options _______________________
#if defined orca
print *,'ORCA GCM fields'
#elif defined rco
print *,'RCO GCM fields'
#elif defined tes
print *,'ACADEMIC TEST fields'
#elif defined simp
print *,'Simpevarp GCM fields'
#elif defined fors
print *,'Forsmark GCM fields'
#elif defined occam66
print *,'OCCAM 1/4 deg GCM fields'
#elif defined occam083
print *,'OCCAM 1/12 deg GCM fields'
#elif defined sigma
print *,'GCM with sigma coordinates fields'
#endif

#if defined tempsalt
print *,'with temperature and salinity fields'
#endif

#if defined turb
print *,'with sub-grid parameterisation'
#endif

#if defined rerun
print *,'Rerun in order to store the Lagrangian stream functions in the different basins'
#endif

#if defined streamxy
print *,'Lagrangian horizontal stream function stored'
#endif

#if defined streamv
print *,'Lagrangian vertical depth stream function stored'
#endif

#if defined streamr
#if defined streamts
print *,'Lagrangian density, temperature and salinity stream function stored'
#else
print *,'Lagrangian density stream function stored'
#endif
#endif

#if defined streamxy
print *,'Lagrangian horizontal stream function stored'
#endif

#if defined tracer
print *,'Lagrangian trajectory particle tracer stored'
#endif

!________________________________________________________________________________________
call init_params
call coordinat

partQuant=1


tseas=1.d0 * real(ngcm)*3600.d0 ! time step between data sets

if(intstep.gt.0) then ! forward 
 intstart=intmin          
 intend  =intmax
elseif(intstep.lt.0) then ! backward
 intstart=intmax
 intend  =intmin
 intspin =-intspin
 intrun  =-intrun    
endif

 if(nqua.eq.1) then ! number of trajectories (per time resolution)
! num=NTRACMAX
 num=partQuant
elseif(nqua.eq.2) then 
 voltr=partQuant 
elseif(nqua.eq.3) then 
 voltr=partQuant
endif

!mask=-1.  ! define start section with ist1,ist2,jst1,jst2
open(21,file=directory//'kmt/kmt',form='unformatted')
!open(21,file=directory//'kmt/mask',form='unformatted')
!open(21,file=directory//'kmt/maskust',form='unformatted')
!open(21,file=directory//'kmt/maskferrysyd',form='unformatted')
!open(21,file=directory//'kmt/maskferrynord',form='unformatted')
read(21)mask
close(21)

do i=1,IMT
 do j=1,JMT
  if(mask(i,j).ne.0 .and. mask(i,j).le.4 .and. j.lt.215) mask(i,j)=-1  ! entire shallow Baltic south of 61N
 enddo
enddo

!mask=-1  ! Finska viken

if(kriva.ne.0) open(56,file=directory//'orm/traj.'//name) ! trajectory path
open(57,file=directory//'orm/traj.ut.'//name)         ! exit position
open(58,file=directory//'orm/traj.in.'//name)         ! entrence position


!??????????????????????????????????? END ???????????????????????????????

! trajectory loops 
call loop

stop
ENDPROGRAM main

!______________ END OF MAIN PROGRAM _______________________________