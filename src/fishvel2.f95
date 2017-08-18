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

! CDV 08/15/17
! make sure to call fishvel2 exactly once per day or the time units will be off: need to convert from mm/day to 1 hr timestep in growth equation (add * dtmin/3600 or leave blank)

 SUBROUTINE fishvel2
 USE mod_param
 USE mod_fish
 USE mod_loopvars
! USE mod_time
 USE mod_traj
! USE mod_vel
 USE mod_grid, ONLY: zed
 IMPLICIT none


 REAL*8   :: length
 REAL*8   :: wmax
 !DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793238462

 length = fish(ntrac,i_length)
 wmax = length*0.3E-3 ! changed to slower swimming speed for Arctic cod
 If(light > 0.0)then
   wfish = wmax *( - tanh(0.2*(zed - z_day)))
 Else
   wfish = wmax *( - tanh(0.2*(zed - z_nit)))
 End if

 return
END SUBROUTINE fishvel2
#endif
