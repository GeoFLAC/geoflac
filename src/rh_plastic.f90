!------ Elasto-Plastic

subroutine plastic(bulkm,rmu,coh,phi,psi,depls,ipls,diss,hardn,s11,s22,s33,s12,de11,de22,de33,de12,&
     ten_off,ndim)
!$ACC routine seq
implicit none

integer, intent(in) :: ndim
real*8, intent(in) :: bulkm, rmu, coh, phi, psi, hardn, de11, de22, de33, de12, ten_off
real*8, intent(inout) :: s11, s22, s33, s12
real*8, intent(out) :: depls, diss
integer, intent(out) :: ipls
real*8, parameter :: pi = 3.14159265358979323846d0
real*8, parameter :: degrad = pi/180.d0
real*8, parameter :: c4d3 = 4.d0/3.d0
real*8, parameter :: c2d3 = 2.d0/3.d0
! press_add formaely was passed by a parameter. in my case it is always zero.
real*8, parameter :: press_add = 0.d0

real*8 sphi, spsi, anphi, anpsi, amc, e1, e2, x1, ten_max, &
     s11i, s22i, s12i, s33i, sdif, s0, rad, si, sii, s1, s2, s3, psdif, &
     fs, alams, dep1, dep3, depm, cs2, si2, dc2, dss
integer icase


! ------------------------------
! Initialization section
! ------------------------------
depls = 0.d0
diss = 0.d0
ipls = 0

sphi  = dsin(phi * degrad)
spsi  = dsin(psi * degrad)
anphi = (1.d0+ sphi) / (1.d0- sphi)
anpsi = (1.d0+ spsi) / (1.d0- spsi)
amc   = 2.0d0 * coh * sqrt (anphi)
e1    = bulkm + c4d3 * rmu
e2    = bulkm - c2d3 * rmu
x1    = (e1 - e2*anpsi + e1*anphi*anpsi - e2*anphi)

if (phi.eq. 0.d0) then
    ten_max=ten_off
else
    ten_max=min(ten_off,coh/(tan(phi*degrad)))
end if

! ---------------
! Running section
! ---------------

!---- get new trial stresses from old, assuming elastic increment
!---- add press (which is positive press = - (sxx+syy)*0.5, 
!---- which has 2 components: add pressure due to application of forces from the top 
!---- and subtract pressure of the fluid
s11i = s11 + (de22 + de33) *e2  + de11 *e1 - press_add
s22i = s22 + (de11 + de33) *e2  + de22 *e1 - press_add
s12i = s12 + de12 * 2.0d0 * rmu
s33i = s33 + (de11 + de22) *e2  + de33 *e1 - press_add
sdif = s11i - s22i
s0   = 0.5d0 * (s11i + s22i)
rad  = 0.5d0 * sqrt(sdif*sdif + 4.0d0 *s12i*s12i)
! principal stresses
si  = s0 - rad
sii = s0 + rad
psdif = si - sii

!--------------------------------------------------------- 
!                         3D version 
!--------------------------------------------------------- 
if (ndim.eq.3) then
    !-- determine case ---
    if (s33i .gt. sii) then
        !- s33 is minor p.s. --
        icase = 3
        s1 = si
        s2 = sii
        s3 = s33i
    elseif (s33i .lt. si) then
        !- s33 is major p.s. --
        icase = 2
        s1 = s33i
        s2 = si
        s3 = sii
    else
        !- s33 is intermediate --
        icase = 1
        s1 = si
        s2 = s33i
        s3 = sii
    endif
endif

!------------------------------------------------------- 
!         2D version 
!------------------------------------------------------- 
if (ndim.eq.2) then 
    icase = 1
    s1 = si 
    s2 = s33i
    s3 = sii
endif 

!--------------------------------------------------------
! Check for tensional failure before the shear failure
!-------------------------------------------------------
 
!----- general tension failure


if (s1 .ge. ten_max) then
    ipls = -5
    s11 = ten_max + press_add
    s22 = ten_max + press_add
    s12 = 0.0d0
    s33 = ten_max + press_add
    return
endif

!- uniaxial tension ... intermediate p.s. ---
if (s2 .ge. ten_max .and. ndim .eq.3) then
    ipls = -6
    s2 = ten_max
    s3 = ten_max
endif

!- partial failure (only if s3 is greater than ten_max) 
if (s3 .ge. ten_max) then
    s3 = ten_max 
    ipls = -7
endif

!- check for shear yield (if fs<0 -> plastic flow)
fs = s1 - s3 * anphi + amc
if (fs .lt. 0.0d0) then
    !-- yielding in shear ----
    if (icase .eq. 1) ipls = -2
    if (icase .eq. 2) ipls = -3
    if (icase .eq. 3) ipls = -4
    alams = fs/(x1+hardn)
    s1 = s1 - alams * (e1 - e2 * anpsi )
    s2 = s2 - alams * e2 * (1.0d0 - anpsi )
    s3 = s3 - alams * (e2 - e1 * anpsi )

    ! Increment of the plastic strain (2nd Invariant)
    dep1 = alams
    dep3 = -alams*anpsi

    ! FOR 2D caculations
    depm = 0.5d0*(dep1+dep3)
    depls = 0.5d0*abs(dep1-dep3)

    ! Dissipation rate
    diss = s1*dep1+s3*dep3
else
    !-- no failure at all (elastic behaviour)
    s11 = s11i + press_add
    s22 = s22i + press_add
    s33 = s33i + press_add
    s12 = s12i
    return
endif

!- general tension failure?
if (s1 .ge. ten_max) then
    ipls = -5
    s11 = ten_max + press_add
    s22 = ten_max + press_add
    s12 = 0.0d0
    s33 = ten_max + press_add
    return
endif

!- uniaxial tension ... intermediate p.s. ---
if (s2 .ge. ten_max .and. ndim .eq. 3) then
    ipls = -6
    s2 = ten_max
    s3 = ten_max
else if (s3 .ge. ten_max) then
    !- uniaxial tension ... minor p.s. ---
    ipls = -7
    s3 = ten_max 
endif

!- direction cosines
if ( psdif .eq. 0.d0 ) then
    cs2 = 1.d0
    si2 = 0.d0
else
    cs2 = sdif / psdif
    si2 = 2.0d0 * s12i / psdif
endif

!- resolve back to global axes
select case (icase)
case (1)
    dc2 = (s1-s3) * cs2
    dss = s1 + s3
    s11 = 0.5d0 * (dss + dc2)
    s22 = 0.5d0 * (dss - dc2)
    s12 = 0.5d0 * (s1 - s3) * si2
    s33 = s2
case (2)
    dc2 = (s2-s3) * cs2
    dss = s2 + s3
    s11 = 0.5d0 * (dss + dc2)
    s22 = 0.5d0 * (dss - dc2)
    s12 = 0.5d0 * (s2 - s3) * si2
    s33 = s1
case (3)
    dc2 = (s1-s2) * cs2
    dss = s1 + s2
    s11 = 0.5d0 * (dss + dc2)
    s22 = 0.5d0 * (dss - dc2)
    s12 = 0.5d0 * (s1 - s2) * si2
    s33 = s3
end select

s11 = s11 + press_add
s22 = s22 + press_add
s33 = s33 + press_add

return
end



!==================================================================
! Prepare plastic properties depending on softening, weighted by phase ratio

subroutine pre_plast (i,j,coh,phi,psi,hardn)
!$ACC routine seq
use arrays
use params
implicit none
integer, intent(in) :: i, j
double precision, intent(out) :: coh, phi, psi, hardn
integer :: iph
double precision :: pls_curr, f, c, d, h, dpl, fma
pls_curr = aps(j,i)

phi = 0
coh = 0
psi = 0

! Strain-Hardening
hardn = 0

do iph = 1, nphase
    if(phase_ratio(iph,j,i) .lt. 0.01d0) cycle

    if(pls_curr < plstrain1(iph)) then
        ! no weakening yet
        f = fric1(iph)
        c = cohesion1(iph)
        d = dilat1(iph)
        h = 0
    else if (pls_curr < plstrain2(iph)) then
        ! Find current properties from linear interpolation
        dpl = (pls_curr - plstrain1(iph)) / (plstrain2(iph) - plstrain1(iph))
        f =  fric1(iph) + (fric2(iph) - fric1(iph)) * dpl
        d = dilat1(iph) + (dilat2(iph) - dilat1(iph)) * dpl
        c = cohesion1(iph) + (cohesion2(iph) - cohesion1(iph)) * dpl
        h = (cohesion2(iph) - cohesion1(iph)) / (plstrain2(iph) - plstrain1(iph))
    else
        ! saturated weakening
        f = fric2(iph)
        c = cohesion2(iph)
        d = dilat2(iph)
        h = 0
    endif

    ! using harmonic mean on friction and cohesion
    ! using arithmatic mean on dilation and hardening
    phi = phi + phase_ratio(iph,j,i) / f
    coh = coh + phase_ratio(iph,j,i) / c
    psi = psi + phase_ratio(iph,j,i) * d
    hardn = hardn + phase_ratio(iph,j,i) * h

enddo

phi = 1 / phi
coh = 1 / coh
fma = min(fmagma(j,i)+fmagma2(j,i),fmagma_max)
if (itype_melting == 1) then
    phi = phi * (1 - (1 - weaken_ratio_plastic) * fma / fmagma_max)
    coh = coh * (1 - (1 - weaken_ratio_plastic) * fma / fmagma_max)
endif

return
end subroutine pre_plast
