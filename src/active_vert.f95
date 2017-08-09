#ifdef fishvel
 !=======================================================================
 ! Subroutine for active vertical behavior of particles
 !=======================================================================

 SUBROUTINE active_vert(z1)

 USE mod_param
 USE mod_loopvars
 USE mod_fish
 IMPLICIT none

 REAL*8   :: length, z_day, z_nit, wmax, wfish
 REAL*8   :: light

 !----------Calc light from lon and lat-----------------

 !----------------Vertical Migration -------------------
 length = fish(ntrac,i_length)
 wfish=0.0

 If (stage(ntrac) .eq. f_pre)then
    z_day = -10.
    z_nit = -10.
 ElseIf (stage(ntrac) .eq. f_post)then
    z_day = -20.
    z_nit = -5.
 Else
    If (stage(ntrac) .eq. f_egg)then
        length = 4.5
    End If
    z_day = -10. !Hsbl/2.
    z_nit = -10.
 End If

 wmax = -length*0.5E-3
 !If(light > 0.0)then
    wfish = wmax *( - tanh(0.2*(z1(ib,jb,kb) - z_day)))
 !Else
    !wfish = wmax *( - tanh(0.2*(z1(ib,jb,kb) - z_nit)))
 !Endif
 print*, length, wmax, z_day, wfish

!-----------------------------------------------------

 return
END SUBROUTINE active_vert
#endif
