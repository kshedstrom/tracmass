#ifdef larval_fish
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
!
! Subroutine of a sedimentation model for the Baltic Sea by Hanna Corell.
! The routine calculates the settling velocity of the particle size modelled.
! The formulas are taken from
! F. Shepard, Submarine Geology, Harper InternationalEdition, 1948
! p. Nielsen, Coastal bottom boundary layers ..., World Scientific, 1992
!
!      subroutine fishvel(rhof,D,visc,wfish) ! If called from main.F
!      by Kate Hedstrom
! CMP 11/13/12
! Modified to add egg density as a function of age group
! Age group is based on stage
! Stage is a function of age and temperature
! CMP 1/29/13
! Size at hatch is a random number from normal distr with mean and std from data
! CMP 5/10/13
! Added a function for temp-dep growth
! CMP 5/14/13
! Egg, yolksac, and preflexion swim towards middle of the mixed layer

! CDV 9/15/16
! Modified code from pollock to Arctic cod
! hatch time based on Ben Laurel's laboratory experiments; converted from hours to days

!CDV 9/15/16
! Updated random hatch length: Arctic cod size at hatch is a random number from normal distribution with mean=5.70 mm and SD=0.48 (Ben Laurel, unpublished data)

!CDV 9/15/16
! commented out Colleen's age-dependent routine: assumed only temperature-dependent growth based on B. Laurel's experiments; added stage-specific temperature-dependent growth models (hatch-10 mm: yolksac/preflexion; 10-25 mm: postflexion; 26-45 mm: age-0 early juvenile)

!CDV 9/15/16
! changed swimming speed based on NE Arctic cod (Vikebo et al. 2007), Arcto-Norwegian cod (Sundby and Fossum 1990), and modeled Arctic cod (Thanasekkos and Fortier 2012)

!CDV 10/06/16
! updated file with Colleen's comments (hatch_hrs to hatch_days, deleted age-dependent growth routine)

 SUBROUTINE fishvel2(temp,salt,rho)     ! If called from loop.F
 USE mod_param
 USE mod_fish
 USE mod_loopvars
 !USE mod_particle
 USE mod_time
 USE mod_traj
 USE mod_vel
 IMPLICIT none


 REAL*8   :: hatchJD
 REAL*8   :: age, egg_hatch, yolk_len, length
 REAL*8   :: hatch_days, hatch_len, g
 REAL*8   :: r1, r2, r, theta
 REAL*8   :: wmax, zday, znit
 REAL*8   :: temp, salt, rho
 INTEGER  :: clock, seed
 !DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793238462

 length = fish(ntrac,i_length)
 CALL SYSTEM_CLOCK(COUNT=clock)
 seed = clock
 CALL RANDOM_SEED(seed)
 !call init_random_seed()

!--------------DEVELOPMENT / GROWTH-----------------------
!EGGS
 If (stage(ntrac) .eq. f_egg)then
   !Time to hatch
   egg_hatch = fish(ntrac,i_hatchtime)
   hatch_days = 28.2817+20.7194*exp(-0.5590*temp) !days to hatch at current temp
   egg_hatch = egg_hatch + (1.0 / hatch_days) * dtmin/tday !accum time to hatch in days
   !egg_hatch=1  ! for simulations starting from hatching locations, comment out formula for egg hatch, use 1 instead
   fish(ntrac,i_hatchtime) = egg_hatch

   If (egg_hatch .ge. 1)then
      stage(ntrac) = f_yolk
      !Length at hatch
      !draw yolk length from a normal distribution with mean=5.70 and stdev=0.48
      CALL RANDOM_NUMBER(r1)
      CALL RANDOM_NUMBER(r2)
      r =(-2.0d0*log(r1))**0.5
      theta = 2.0d0*PI*r2
      yolk_len = 5.70 + (0.48**2)*r*sin(theta)
      fish(ntrac,i_hatchlength) = yolk_len
      fish(ntrac,i_length) = yolk_len
      fish(ntrac,i_jd) = currJDtot         !Use JD to calc dph
   Endif
 Endif !Egg if

! commented out age-dependent routine, updated for 3 growth models
!YOLKSAC/PREFLEXION: from hatch to 10 mm
If (stage(ntrac) .eq. f_yolk)then
   hatchJD = fish(ntrac,i_jd)
   length = fish(ntrac,i_length)
   g = (0.0735+(0.0149*temp)+(-0.0013*(temp**2))) * dtmin / tday	!growth in mm from B. Laurel (hatch to 10mm, unpublished data)
   If (g .ge. 0)then
      length = length + g
   Endif
   fish(ntrac,i_length) = length
   If (length .ge. 10.0)then
      stage(ntrac) = f_post
   Endif
 Endif !Yolk if

!POSTFLEXION: growth from 10-25 mm
 If (stage(ntrac) .eq. f_post)then
   hatchJD = fish(ntrac,i_jd)
   length = fish(ntrac,i_length)
   g = (0.0369+(0.0583*temp)+(-0.0044*(temp**2))) * dtmin / tday 	!growth in mm from B. Laurel (10-15 mm,unpublished data)
   If (g .ge. 0)then
      length = length + g
   Endif
   fish(ntrac,i_length) = length
   If (length .ge. 25.0)then
      stage(ntrac) = f_ejuv
   Endif
 Endif !Pre if


!TRANSFORMATION-EARLY JUVENILES: growth from 26-70 mm (45 mm to demersal age-0/end of simulation)
 If (stage(ntrac) .eq. f_ejuv)then
   hatchJD = fish(ntrac,i_jd)
   length = fish(ntrac,i_length)
   g = (0.1377+(0.0311*temp)+(0.0041*(temp**2))+(-0.0004*(temp**3))) * dtmin / tday 	!growth in mm from B. Laurel (40-60 mm,unpublished data)
   If (g .ge. 0)then
      length = length + g
   Endif
   fish(ntrac,i_length) = length
   If (length .ge. 45.0)then
      stage(ntrac) = f_ljuv
   Endif
 Endif !Post if



!----------------VERTICAL BEHAVIOR----------------
!Non-DVM larvae are directed to the middle of surface layer (-5 m): initial simulation with all stages at surface

 !length = fish(ntrac,i_length)
 !wfish=0.0

!POSTFLEXION
 !If (stage(ntrac) .eq. f_post)then
 ! zday = -5.
 ! znit = -5.
 !Endif !Post if

!TRANSFORMATION
 !If (stage(ntrac) .eq. f_ejuv)then
!   zday = -5.
!   znit = -5.
! Endif !Trans if

! wmax = length*0.3E-3 ! changed to slower swimming speed for Arctic cod
 !If(light > 0.0)then
!   wfish = wmax *( - tanh(0.2*(depth1 - zday)))
 !Else
!   wfish = wmax *( - tanh(0.2*(depth1 - znit)))
 !End if

 return
END SUBROUTINE fishvel2
#endif
