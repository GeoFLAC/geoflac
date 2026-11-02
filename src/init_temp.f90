! Initiate temperature profile

subroutine init_temp
use arrays
use params
implicit none
integer :: i, j, k, n, i1, i2, ixc, iwidth, kk, u, ios
character(len=256) :: iomsg
double precision :: ratio, age_1n, tp1n, tp2n, y, amp, pert, pert2

!  Read distribution of temperatures from the dat file
if (irtemp .gt. 0) then
    open(newunit=u, file=tempfile, status='old', iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
        call SysMsg('INIT_TEMP: Cannot open file with temperatures initial distrib: ' // trim(iomsg))
        stop 21
    endif
    do i = 1, nx
        do j = 1, nz
            read(u, *, iostat=ios) temp(j,i)
!     if(temp(j,i).ge.1000.d0) temp(j,i) = 1000.d0
            if (ios /= 0) then
                call SysMsg('INIT_TEMP: Error reading file with temperatures initial distrib!')
                stop 21
            endif
        enddo
    enddo
    close(u)
else
    !!  geotherm of a given age accross the box with variable age
    do n = 1, nzone_age
        if (n /= 1) then
            if (iph_col_trans(n-1) == 1) cycle
        endif

        do i = ixtb1(n), ixtb2(n)
            if (iph_col_trans(n) == 1) then
                i1 = ixtb1(n)
                i2 = ixtb2(n)
                ratio = (cord(1,i,1) - cord(1,i1,1)) / (cord(1,i2,1) - cord(1,i1,1))
                age_1n = age_1(n) + (age_1(n+1) - age_1(n)) * ratio
                tp1n = tp1(n) + (tp1(n+1) - tp1(n)) * ratio
                tp2n = tp2(n) + (tp2(n+1) - tp2(n)) * ratio
            else
                age_1n = age_1(n)
                tp1n = tp1(n)
                tp2n = tp2(n)
            endif
            call init_geotherm_profile(n, i, age_1n, tp1n, tp2n)
        enddo
    enddo
endif

! DISTRIBUTE SOURCES in elements
do j = 1,nz-1
    y = -( cord(j+1,1,2)+cord(j,1,2) )/2 / 1000
    source(j,1:nx-1) = hs*exp(-y/hr)
end do

do i = 1, inhom
    ! Initial gaussian temperature perturbation

    ! vertical gaussian
    ! x between (ix1, ix2), y between (iy1, iy2)
    ! gaussian dist. in x direction
    ! linear dist in y direction
    ! xinitaps: amplitude of gaussian
    ! inphase: not used
    if (igeom(i).eq.11) then
        ixc  = (ix1(i)+ix2(i))/2
        iwidth = (ix2(i)-ix1(i))
        amp = xinitaps(i)
        do j = ix1(i),ix2(i)
            pert = amp*exp(-(float(j-ixc)/(0.25d0*float(iwidth)))**2)
            do k = iy1(i),iy2(i)
                pert2 = 1.0d0*(k-iy1(i)) / (iy2(i) - iy1(i))
                temp(k,j) = min(t_bot, temp(k,j)+pert*pert2)
            enddo
        enddo
    endif

    ! slant gaussian
    ! x between (ix1, ix2) at top, shift 1-grid to right for every depth grid
    ! z between (iy1, iy2)
    ! xinitaps: amplitude of gaussian
    ! inphase: not used
    if (igeom(i).eq.13) then
        ixc  = (ix1(i)+ix2(i))/2
        iwidth = (ix2(i)-ix1(i))
        amp = xinitaps(i)
        do k = iy1(i),iy2(i)
            kk = k - iy1(i)
            do j = ix1(i),ix2(i)
                pert = amp*exp(-(float(j-ixc)/(0.25d0*float(iwidth)))**2)
                temp(k,j+kk) = max(t_top, min(t_bot, temp(k,j+kk)+pert))
                !print *, k, j, pert
            enddo
        enddo
    endif

    ! slant gaussian
    ! x between (ix1, ix2) at top, shift 1-grid to left for every depth grid
    ! z between (iy1, iy2)
    ! xinitaps: amplitude of gaussian
    ! inphase: not used
    if (igeom(i).eq.14) then
        ixc  = (ix1(i)+ix2(i))/2
        iwidth = (ix2(i)-ix1(i))
        amp = xinitaps(i)
        do k = iy1(i),iy2(i)
            kk = k - iy1(i)
            do j = ix1(i),ix2(i)
                pert = amp*exp(-(float(j-ixc)/(0.25d0*float(iwidth)))**2)
                temp(k,j-kk) = max(t_top, min(t_bot, temp(k,j-kk)+pert))
                !print *, k, j, pert
            enddo
        enddo
    endif
enddo

!call RedefineTemp

return
end subroutine init_temp


subroutine init_geotherm_profile(n, i, age_1n, tp1n, tp2n)
    !$ACC routine vector
    use arrays
    use params
    implicit none
    integer, intent(in) :: n, i
    double precision, intent(in) :: age_1n, tp1n, tp2n
    integer :: j, k, kl
    double precision :: diffusivity, age, f, yL0, age_init, tau_d, y, tss, tt, ymoho, cond, dens_c, rr, tm, qm, bot_dep, temp1, ylayer1, ylayer2
    double precision :: cond_c, cond_m, auc, alc, am, xdz, Fpart
    double precision :: Q_loc(nz)
    double precision, parameter :: pi = 3.14159265358979323846d0


    if (ictherm(n)==1) then
        !! Oceanic geotherm (half space cooling model, T&S 3rd ed. Eq(4.113))
        diffusivity = 1.d-6
        age = age_1n
        !$ACC loop
        do j = 1,nz
            f = (cord(1,i,2)-cord(j,i,2)) / sqrt(4 * diffusivity * age * 1.d6 * sec_year)
            temp(j,i) = t_top + (t_bot - t_top) * erf(f)
            !print *, j, age, -cord(j,i,2), temp(j,i)
        enddo
    elseif (ictherm(n)==2) then
        !! Oceanic geotherm (plate cooling cooling model, T&S 3rd ed. Eq(4.130))
        diffusivity = 1.d-6
        age = age_1n
        yL0 = tp1n    ! plate thickness in km
        age_init = age*sec_year*1d6
        tau_d = yL0*yL0*1d6 / (pi*pi*diffusivity)
        !$ACC loop
        do j = 1,nz
            ! depth in km
            y = (cord(1,i,2) - cord(j,i,2)) * 1d-3
            tss = t_top + (t_bot - t_top) * y / yL0
            tt = 0.d0
            do k = 1,100
                tt = tt + 1.d0/k * exp(-k*k*age_init/tau_d) * sin(pi*k*y/yL0)
            enddo
            temp(j,i) = tss + 2.d0/pi*(t_bot-t_top)*tt
            if(temp(j,i)>t_bot .or. y>yL0) temp(j,i) = t_bot
        enddo
    elseif (ictherm(n)==12) then
        !! Continental geotherm (plate cooling model with radiogenic heating)
        !
        ! Starting from the steady state (ss) solution as in T&S 3rd ed. Eq(4.30)
        ! Let the ss moho temperature be tm and ss heatflux be qm.
        ! qm = cond * (t_bot - tm) / (yL0 - ymoho)
        ! Substitute qm to 4.30 to solve for tm
        age = age_1n
        yL0 = tp1n    ! plate thickness in km
        ymoho = tp2n  ! crust thickness in km
        cond = 3.3d0
        dens_c = 2700.d0
        diffusivity = 1.d-6
        !   write(*,*) rzbo, hs, hr
        age_init = age*sec_year*1d6
        tau_d = yL0*yL0*1d6 / (pi*pi*diffusivity)
        rr = ymoho / (yL0 - ymoho)
        tm = (t_top + rr*t_bot + dens_c*hs*hr*hr*1d6/cond*(1d0-exp(-ymoho/hr))) / (1 + rr)
        qm = cond * (t_bot - tm) / (yL0 - ymoho)
        !$ACC loop
        do j = 1,nz
            ! depth in km
            y = (cord(1,i,2) - cord(j,i,2)) * 1d-3

            ! ss part with radiogenic heat
            if (y <= ymoho) then
                tss = t_top + qm/cond*y + (dens_c*hs*hr*hr*1.d+6/cond)*(1d0-exp(-y/hr))
            elseif (y <= yL0) then ! below moho, inside lithosphere
                tss = tm + qm/cond*(y-ymoho)
            else
                tss = t_bot
            endif
            ! time-dependent part
            ! see T&S 3rd ed. Eq(4.130)
            tt = 0.d0
            do k = 1,100
                tt = tt + 1.d0/k * exp(-k*k*age_init/tau_d) * sin(pi*k*y/yL0)
            enddo
            temp(j,i) = tss + 2.d0/pi*(t_bot-t_top)*tt
            if(temp(j,i)>t_bot .or. y>yL0) temp(j,i) = t_bot
            !       write(*,*) j,y,tss,yL0,tt
        enddo
    elseif (ictherm(n)==21) then
        !! Constant geotherm gradient at top layer, then T=t_bot all the way to the bottom
        bot_dep = age_1n   ! bottom of the top layer
        do j = 1,nz
            y = (cord(1,i,2)-cord(j,i,2))*1.d-3
            temp(j,i) = y * (t_bot-t_top) / bot_dep
            if(temp(j,i).gt.t_bot) temp(j,i) = t_bot
        enddo
    elseif (ictherm(n)==22) then
        !! Constant geotherm gradient at top two layers, then T=t_bot all the way to the bottom
        temp1 = age_1n
        ylayer1 = tp1n  ! depth in km
        ylayer2 = tp2n  ! depth in km
        !$ACC loop
        do j = 1,nz
            y = (cord(1,i,2)-cord(j,i,2))*1.d-3
            if (y <= ylayer1) then
                temp(j,i) = t_top + (temp1 - t_top) * y / ylayer1
            else
                temp(j,i) = temp1 + (t_bot - temp1) * (y - ylayer1) / (ylayer2 - ylayer1)
            endif
            if(temp(j,i).gt.t_bot) temp(j,i) = t_bot
        enddo
    elseif (ictherm(n)==13) then
        !! Continental geotherm (Hasterok & Chapman 2011) with local zone parameters.
        !! Unlike the other ictherm branches (closed-form formulas), this one
        !! integrates the steady-state 1-D heat equation directly, layer by
        !! layer down the column: each layer has its own conductivity and
        !! radiogenic heat generation rate, and Q_loc tracks the conductive
        !! heat flux at the top of the current element, decreasing with depth
        !! as each layer's own heat generation is subtracted off.
        !! age_1n -> Surface heat flow (mW/m^2)
        !! tp1n   -> Upper crust thickness (km)
        !! tp2n   -> Lower crust thickness (km)
        cond_c = 2.3d0
        cond_m = 3.3d0
        Fpart = 0.74d0 ! Partition coefficient: fraction of surface heat flow
                       ! attributed to upper-crustal radiogenic sources
        alc = 0.4d-6 ! Lower crust heat generation
        am = 0.02d-6 ! Mantle heat generation

        ! Upper-crust heat generation rate that reproduces the given surface
        ! heat flow through the given crust thickness
        auc = (1.0d0 - Fpart) * age_1n*1.0d-3 / (tp1n*1000.0d0)

        Q_loc(1) = age_1n*1.0d-3
        temp(1,i) = t_top
        kl = 0 ! becomes 1 once the conductive profile reaches t_bot,
               ! switching the rest of the column to isothermal mantle

        do j = 1, nz-1
            ! Depth in km
            y = (cord(1,i,2)-cord(j,i,2))*1.0d-3
            ! Distance between elements
            xdz = abs(cord(j+1,i,2)-cord(j,i,2))

            ! Upper crust: conductive temperature step plus the local
            ! radiogenic contribution, then reduce the flux carried
            ! downward by the heat generated in this element
            if (y.lt.tp1n) then
                temp(j+1,i) = temp(j,i) + (Q_loc(j)/cond_c)*xdz - (auc)/(2.0d0*cond_c)*(xdz**2)
                Q_loc(j+1) = Q_loc(j) - auc*xdz
                ! source(:,nx) doesn't exist (element array, one fewer column
                ! than nodes), hence the i<nx guard on every source(j,i) write
                ! below
                if (i .lt. nx) source(j,i) = auc*1.0d-3
            endif

            ! Lower crust: same as above with the lower-crust conductivity
            ! and heat generation
            if (y.ge.tp1n.and.y.lt.(tp2n+tp1n)) then
                temp(j+1,i) = temp(j,i) + (Q_loc(j)/cond_c)*xdz - (alc)/(2.0d0*cond_c)*(xdz**2)
                Q_loc(j+1) = Q_loc(j) - alc*xdz
                if (i .lt. nx) source(j,i) = alc*1.0d-3
            endif

            ! Mantle lithosphere: same integration with mantle conductivity
            ! and heat generation, until the conductive temperature reaches
            ! the mantle reference temperature (t_bot)
            if (y.ge.(tp2n+tp1n).and.kl.eq.0) then
                temp(j+1,i) = temp(j,i) + (Q_loc(j)/cond_m)*xdz - (am)/(2.0d0*cond_m)*(xdz**2)
                Q_loc(j+1) = Q_loc(j) - am*xdz
                if (i .lt. nx) source(j,i) = am*1.0d-3
                if (temp(j+1,i).ge.t_bot) kl = 1
            endif

            ! Asthenosphere: below the lithosphere the mantle is convecting,
            ! not conducting, so clamp to the (isothermal) reference
            ! temperature instead of continuing the conductive integration
            if (kl.eq.1) then
                temp(j+1,i) = t_bot
                if (i .lt. nx) source(j,i) = 0.0d0
            endif
        enddo
    else
        !call sysmsg('init_temp: ictherm not supported!')
        stop 1
    endif
end subroutine init_geotherm_profile


subroutine sidewalltemp(i1, i2)
  !$ACC routine(init_geotherm_profile) vector

  use arrays, only : temp, cord, source
  use params
  implicit none
  ! This subroutine is intended for remeshing.

  integer, intent(in) :: i1, i2
  double precision :: age_1n, tp1n, tp2n
  integer :: n, i

  if(i1 == 1) then
      ! left sidewall
      n = 1
  else
      ! right sidewall
      n = nzone_age
  endif
  age_1n = age_1(n)
  tp1n = tp1(n)
  tp2n = tp2(n)

  !$ACC parallel loop async(1)
  do i = i1, i2
    call init_geotherm_profile(n, i, age_1n, tp1n, tp2n)
  enddo

  if(i1 == 1) then
      !$ACC parallel loop async(1)
      do i = i1, i2-1
          source(1:nz-1,i) = source(1:nz-1,i2)
      enddo
  else
      !$ACC parallel loop async(1)
      do i = i1, i2-1
          source(1:nz-1,i) = source(1:nz-1,i1-1)
      enddo
  endif
  return
end subroutine sidewalltemp
