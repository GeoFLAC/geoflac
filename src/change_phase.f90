subroutine change_phase
!$ACC routine(newphase2marker) worker
!$ACC routine(count_phase_ratio) seq

USE marker_data
use arrays
use params
use phases

implicit none

integer :: jj, j, i, iph, &
           jbelow, k, kinc, kk, n, iph1, iph2, jocmoho
double precision :: yy, dep2, dep3, depth, press, quad_area, &
                    tmpr, trtmpr, trpres, trpres2, &
                    solidus, pmelt, total_phase_ratio, dt_melt

! max. depth (m) of eclogite phase transition, no serpentinization below it
real*8, parameter :: max_basalt_depth = 150.d3
! min. temperature (C) of eclogite phase transition
real*8, parameter :: min_eclogite_temp = 400.d0
real*8, parameter :: min_eclogite_depth = 20d3
real*8, parameter :: mantle_density = 3000.d0
real*8, parameter :: max_melting_depth = 200.d3

! temperature (C) of serpentine phase transition
real*8, parameter :: serpentine_temp = 550.d0

! max. depth (m) to check MOR basalt
real*8, parameter :: check_oc_depth = 30.d3
real*8, parameter :: new_oc_thickness = 7.d3
real*8, parameter :: mor_temp = 1215.d0
real*8, parameter :: max_mant_melting_depth = 150.d3

!$ACC kernels async(1)
itmp = 0  ! indicates which element has phase-changed markers
!$ACC end kernels

! MORB (old version)
!$ACC serial async(1)
do j = 1, nz-1
    dep2 = 0.25d0*(cord(j,1,2)+cord(j+1,1,2)+cord(j,2,2)+cord(j+1,2,2))
    if (cord(1,1,2) - dep2 >= check_oc_depth) exit
end do
jj = min(max(2, j), nz-1)
!$ACC end serial
!$ACC serial async(1)
do j = 1, nz-1
    dep2 = 0.25d0*(cord(j,1,2)+cord(j+1,1,2)+cord(j,2,2)+cord(j+1,2,2))
    if (cord(1,1,2) - dep2 >= new_oc_thickness) exit
end do
jocmoho = min(max(2, j), nz-1)
!$ACC end serial

!!$ACC parallel loop async(1)
!do i = 1, nx-1
!  tmpr = 0.25d0*(temp(jj,i)+temp(jj+1,i)+temp(jj,i+1)+temp(jj+1,i+1))
!  if ( tmpr > mor_temp ) then
!    call newphase2marker(1,jocmoho,i,i,kocean1)
!  end if
!end do



!open (1,file='change_ckeck.0')
!$OMP parallel private(kk,i,j,k,n,tmpr,depth,iph,press,jbelow,trpres,trpres2,kinc,quad_area,yy,solidus,trtmpr)
!$OMP do schedule(guided)
!$ACC parallel loop async(1)
do kk = 1 , nmarkers
    if (mark_dead(kk).eq.0) cycle
    
    ! from ntriag, get element number
    n = mark_ntriag(kk)
    k = mod(n - 1, 2) + 1
    j = mod((n - k) / 2, nz-1) + 1
    i = (n - k) / 2 / (nz - 1) + 1

    ! interpolate y-coordinate and temperature on markers
    if (k .eq. 1) then
       yy = cord(j,i,2)*mark_a1(kk) + cord(j+1,i,2)*mark_a2(kk) + cord(j,i+1,2)*(1-mark_a1(kk)-mark_a2(kk))
       tmpr = temp(j,i)*mark_a1(kk) + temp(j+1,i)*mark_a2(kk) + temp(j,i+1)*(1-mark_a1(kk)-mark_a2(kk))
    else
       yy = cord(j,i+1,2)*mark_a1(kk) + cord(j+1,i,2)*mark_a2(kk) + cord(j+1,i+1,2)*(1-mark_a1(kk)-mark_a2(kk))
       tmpr = temp(j,i+1)*mark_a1(kk) + temp(j+1,i)*mark_a2(kk) + temp(j+1,i+1)*(1-mark_a1(kk)-mark_a2(kk))
    endif

    ! depth below the surface in m
    depth = 0.5d0*(cord(1,i,2)+cord(1,i+1,2)) - yy

    ! # of markers inside quad
    kinc = nmark_elem(j,i)

    ! If location of this element is too deep, this marker is already
    ! too deep in the mantle, where there is no significant phase change.
    if (depth > max_melting_depth) cycle

    iph = mark_phase(kk)
    !write(1,'(i3)') iph
    ! Rules of phase changes
    select case(iph)
    case (kcont1, kcont2)
        ! subduction below continent, continent becomes weaker to
        ! facilitate further subduction
        !do jbelow = min(j+1,nz-1), min(j+3,nz-1)
        !    if(phase_ratio(kocean1,jbelow,i) > 0.8d0 .or. &
        !         phase_ratio(kocean2,jbelow,i) > 0.8d0 .or. &
        !         phase_ratio(karc1,jbelow,i) > 0.8d0 .or. &
        !         phase_ratio(ksed2,jbelow,i) > 0.8d0) then
        !        !$ACC atomic write
        !        !$OMP atomic write
        !        itmp(j,i) = 1
        !        mark_phase(kk) = kweak
        !        exit
        !    endif
        !enddo

        ! XXX: middle crust with high dissipation becomes weaker,
        ! this helps with localization
        !if(tmpr > 300.d0 .and. tmpr < 400.d0 &
        !     .and. stressII(j,i)*strainII(j,i) > 4.d6) then
        !    !$ACC atomic write
        !    !$OMP atomic write
        !    !itmp(j,i) = 1
        !    mark_phase(kk) = kweakmc
        !endif
        cycle
    case (kmant1, kmant2)
        ! subuducted oceanic crust below mantle, mantle is serpentinized
        if(depth > max_basalt_depth) cycle

        ! Phase diagram taken from Ulmer and Trommsdorff, Science, 1995
        ! Fixed points (730 C, 2.1 GPa) (500 C, 7.5 GPa)
        trpres = 2.1d9 + (7.5d9 - 2.1d9) * (tmpr - 730.d0) / (500.d0 - 730.d0)
        ! Fixed points (730 C, 2.1 GPa) (670 C, 0.6 GPa)
        trpres2 = 2.1d9 + (0.6d9 - 2.1d9) * (tmpr - 730.d0) / (670.d0 - 730.d0)
        press = mantle_density * g * depth
        solidus = -5.02d0 * abs(press*1d-9)**2 + 1.32d2 * abs(press*1d-9) + 1.09d3 + 35.d0
        if (press < trpres .and. press > trpres2) then
            do jbelow = min(j+1,nz-1), min(j+nelem_serp,nz-1)
                if(phase_ratio(kocean1,jbelow,i) > 0.8d0 .or. &
                    phase_ratio(kocean2,jbelow,i) > 0.8d0 .or. &
                    phase_ratio(ksed1,jbelow,i) > 0.8d0 .or. &
                    phase_ratio(ksed2,jbelow,i) > 0.8d0) then
                    !$ACC atomic write
                    !$OMP atomic write
                    itmp(j,i) = 1
                    mark_phase(kk) = kserp
                    exit
                endif
            enddo
        !--------------------------!
        ! hydrous mantle -> dehydrous mantle
        else if (tmpr >= solidus + 35d0) then
            !$ACC atomic write
            !$OMP atomic write
            itmp(j,i) = 1
            mark_phase(kk) = kdrymant
        endif
    case (koceandry, kocean1, kocean2)
        ! basalt -> eclogite
        ! phase change pressure
        ! Phase Diagram taken from Hacker, JGR, 2003 (Figure 8 or Figure 1)
        ! Fixed points (400 C, 5.1 GPa) (536 C, 0 GPa)
        trpres = -0.3d9 + 2.2d6*tmpr
        trpres2 = 20.1d9 - 0.0375d9*tmpr
        press = mantle_density * g * depth
        if (tmpr < min_eclogite_temp .or. depth < min_eclogite_depth .or. press < trpres .or. press < trpres2) cycle
        !$ACC atomic write
        !$OMP atomic write
        itmp(j,i) = 1
        mark_phase(kk) = keclg
    case (kserp)
        ! dehydration, serpentinite -> hydrated mantle
        ! Phase diagram taken from Ulmer and Trommsdorff, Science, 1995
        ! Fixed points (730 C, 2.1 GPa) (500 C, 7.5 GPa)
        trpres = 2.1d9 + (7.5d9 - 2.1d9) * (tmpr - 730.d0) / (500.d0 - 730.d0)
        ! Fixed points (730 C, 2.1 GPa) (650 C, 0.2 GPa)
        trpres2 = 2.1d9 + (0.2d9 - 2.1d9) * (tmpr - 730.d0) / (650.d0 - 730.d0)
        press = mantle_density * g * depth
        if (tmpr < serpentine_temp .or. (press < trpres .and. press > trpres2)) cycle
        !$ACC atomic write
        !$OMP atomic write
        itmp(j,i) = 1
        mark_phase(kk) = kmant1
    case (ksed1)
        ! compaction, uncosolidated sediment -> sedimentary rock
        if (tmpr < 80d0 .and. depth < 2d3) cycle
        !$ACC atomic write
        !$OMP atomic write
        itmp(j,i) = 1
        mark_phase(kk) = ksed2
    case (ksed2)
        ! sedimentary rock -> metamorphic sedimentary rock
        if (tmpr < 250d0 .and. depth < 7d3) cycle
        !$ACC atomic write
        !$OMP atomic write
        itmp(j,i) = 1
        mark_phase(kk) = kmetased
    case (kmetased)
        ! dehydration, sedimentary rock -> schist
        ! from metapelite KFMASH petrogenic grid,
        ! invariant points i5 (710.8 C, 0.880 GPa) i3 (662.3 C, 0.704 GPa)
        ! Fig 1, Wei, Powell, Clarke, J. Metamorph. Geol., 2004.
        trpres = 0.88d9 - (0.88d9 - 0.704d9) * (710.8d0 - tmpr) / (710.8d0 - 662.3d0)
        press = mantle_density * g * depth
        if (press < trpres ) cycle
        !$ACC atomic write
        !$OMP atomic write
        itmp(j,i) = 1
        mark_phase(kk) = kschist
    case (kschist)
        ! schist -> amphibolite
        trpres = 0.0d9 - (0.0d9 - 1.2d9) * (500.0d0 - tmpr) / (500.0d0 - 600.0d0)
        press = mantle_density * g * depth
        if (press > trpres ) cycle
        !$ACC atomic write
        !$OMP atomic write
        itmp(j,i) = 1
        mark_phase(kk) = kamphi
    case (khydmant)
        ! dehydration of chlorite
        ! Phase diagram from Grove et al. Nature, 2009
        trtmpr = 880 - 35d-9 * (depth - 62d3)**2
        if (tmpr < trtmpr) cycle
        !$ACC atomic write
        !$OMP atomic write
        itmp(j,i) = 1
        mark_phase(kk) = kmant1

    end select

enddo
!$OMP end do
!$OMP end parallel
!close(1)

!$ACC kernels async(2)
! storing plastic strain in temporary array
dummye(1:nz-1,1:nx-1) = aps(1:nz-1,1:nx-1)
!$ACC end kernels
!$ACC wait(2)

!$OMP parallel do private(iph, kinc)
!$ACC parallel loop collapse(2) async(1)
! recompute phase ratio of those changed elements
do i = 1, nx-1
    do j = 1, nz-1

        ! skip unchanged element
        if (itmp(j,i) == 0) cycle

        ! the phase of this element is the most abundant marker phase
        call count_phase_ratio(j,i)

        ! When phase change occurs, the mineral would recrystalize and lost
        ! all plastic strain associated with this marker.
        kinc = nmark_elem(j,i)
        aps(j,i) = max(aps(j,i) - dummye(j,i) / float(kinc), 0d0)
    enddo
enddo
!$OMP end parallel do

if (itype_melting == 1) then
    !$OMP parallel do private(tmpr, yy, depth, solidus, pmelt, total_phase_ratio, press, dt_melt)
    !$ACC parallel loop collapse(2) async(1)
    do i = 1, nx-1
        do j = 1, nz-1
            fmelt2(j,i) = 0
            ! hydrous mantle melting
            ! solidus from Katz, 2003 Geochemistrym Geophysics, Geosystems
            total_phase_ratio = phase_ratio(kmant1,j,i) + phase_ratio(kmant2,j,i)
            if (total_phase_ratio > 0.6d0 .and. cord(j,i,2) > -max_mant_melting_depth) then
                tmpr = 0.25d0 * (temp(j,i)+temp(j,i+1)+temp(j+1,i)+temp(j+1,i+1))
                yy = 0.25d0 * (cord(j,i,2)+cord(j,i+1,2)+cord(j+1,i,2)+cord(j+1,i+1,2))
                depth = 0.5d0*(cord(1,i,2)+cord(1,i+1,2)) - yy
                press = abs(mantle_density * g * depth * 1d-9) ! GPa

                solidus = -5.02d0*press**2+1.32d2*press+1.09d3
                if (tmpr > solidus) then
                    ! fraction of partial melting (enthalpy formulation)
                    ! F = cp * (T - T_solidus) / L_fusion, capped at 10%
                    if (latent_heat_magma > 0.d0) then
                        pmelt = min(cp(kmant1) * (tmpr - solidus) / latent_heat_magma, 0.1d0)
                    else
                        pmelt = min((tmpr - solidus) / 35.d0 * 0.1d0, 0.1d0)
                    endif
                    fmelt2(j,i) = pmelt * total_phase_ratio

                    ! Pull temperature back to solidus (latent heat buffering)
                    ! Subtract overshoot (T_apparent - T_solidus) from element's 4 corner nodes
                    dt_melt = tmpr - solidus
                    !$ACC atomic update
                    !$OMP atomic update
                    temp(j,i) = temp(j,i) - dt_melt
                    !$ACC atomic update
                    !$OMP atomic update
                    temp(j,i+1) = temp(j,i+1) - dt_melt
                    !$ACC atomic update
                    !$OMP atomic update
                    temp(j+1,i) = temp(j+1,i) - dt_melt
                    !$ACC atomic update
                    !$OMP atomic update
                    temp(j+1,i+1) = temp(j+1,i+1) - dt_melt
                endif
            endif
         enddo
    enddo
    !$OMP end parallel do

    !$OMP parallel do private(tmpr, yy, depth, solidus, pmelt, total_phase_ratio, dt_melt)
    !$ACC parallel loop collapse(2) async(1)
    do i = 1, nx-1
        do j = 1, nz-1
            fmelt(j,i) = 0

            ! sedimentary rock melting
            ! solidus from Nichols, 1994 Nature
            total_phase_ratio = phase_ratio(ksed1,j,i) + phase_ratio(ksed2,j,i) + &
                              phase_ratio(kmetased,j,i) + phase_ratio(kschist,j,i)
            if (total_phase_ratio > 0.6d0 .and. cord(j,i,2) > -max_melting_depth) then
                tmpr = 0.25d0 * (temp(j,i)+temp(j,i+1)+temp(j+1,i)+temp(j+1,i+1))

                ! depth below the surface in m
                yy = 0.25d0 * (cord(j,i,2)+cord(j,i+1,2)+cord(j+1,i,2)+cord(j+1,i+1,2))
                depth = 0.5d0*(cord(1,i,2)+cord(1,i+1,2)) - yy

                solidus = max(680+0.6d-3*(depth-140d3), 930-313*(1-exp(-depth/7d3)))
                if (tmpr > solidus) then
                    ! fraction of partial melting (enthalpy formulation)
                    ! F = cp * (T - T_solidus) / L_fusion, capped at 10%
                    if (latent_heat_magma > 0.d0) then
                        pmelt = min(cp(ksed2) * (tmpr - solidus) / latent_heat_magma, 0.1d0)
                    else
                        pmelt = min((tmpr - solidus) / 50.d0 * 0.1d0, 0.1d0)
                    endif
                    fmelt(j,i) = pmelt * total_phase_ratio
                    !print *, j, i, tmpr, pmelt

                    ! Pull temperature back to solidus (latent heat buffering)
                    ! Subtract overshoot (T_apparent - T_solidus) from element's 4 corner nodes
                    dt_melt = tmpr - solidus
                    !$ACC atomic update
                    !$OMP atomic update
                    temp(j,i) = temp(j,i) - dt_melt
                    !$ACC atomic update
                    !$OMP atomic update
                    temp(j,i+1) = temp(j,i+1) - dt_melt
                    !$ACC atomic update
                    !$OMP atomic update
                    temp(j+1,i) = temp(j+1,i) - dt_melt
                    !$ACC atomic update
                    !$OMP atomic update
                    temp(j+1,i+1) = temp(j+1,i+1) - dt_melt
                endif
            endif
        enddo
    enddo
    !$OMP end parallel do

    !$OMP parallel do private(tmpr, yy, depth, solidus, pmelt, total_phase_ratio, press, dt_melt)
    !$ACC parallel loop collapse(2) async(1)
    do i = 1, nx-1
        do j = 1, nz-1

            ! basalt and eclogite rock melting
            ! solidus from Gutscher, 2000 Geology
            total_phase_ratio = phase_ratio(kocean1,j,i) + phase_ratio(kocean2,j,i) &
                                + phase_ratio(koceandry,j,i) + phase_ratio(keclg,j,i)
            if (total_phase_ratio > 0.6d0 .and. cord(j,i,2) > -max_melting_depth) then
                tmpr = 0.25d0 * (temp(j,i)+temp(j,i+1)+temp(j+1,i)+temp(j+1,i+1))

                ! depth below the surface in m
                yy = 0.25d0 * (cord(j,i,2)+cord(j,i+1,2)+cord(j+1,i,2)+cord(j+1,i+1,2))
                depth = 0.5d0*(cord(1,i,2)+cord(1,i+1,2)) - yy
                press = 3000*10*depth/1d9

                !shaded area from Gutscher, 2000 Geology
                if (press < 1d0) then
                    solidus = 1050d0 - 420d0*(1d0 - exp(-press*3.3d0))
                elseif (press > 2.7d0) then
                    solidus = (press + 14d0)*43d0
                else
                    solidus=630d0 + 13d0*press**2
                endif

                if (tmpr > solidus) then
                    ! fraction of partial melting (enthalpy formulation)
                    ! F = cp * (T - T_solidus) / L_fusion, capped at 10%
                    if (latent_heat_magma > 0.d0) then
                        pmelt = min(cp(kocean1) * (tmpr - solidus) / latent_heat_magma, 0.1d0)
                    else
                        pmelt = min((tmpr - solidus) / 20.d0 * 0.1d0, 0.1d0)
                    endif
                    !$ACC atomic update
                    !$OMP atomic update
                    fmelt(j,i) = fmelt(j,i) + pmelt * total_phase_ratio
                    
                    ! Pull temperature back to solidus (latent heat buffering)
                    ! Subtract overshoot (T_apparent - T_solidus) from element's 4 corner nodes
                    dt_melt = tmpr - solidus
                    !$ACC atomic update
                    !$OMP atomic update
                    temp(j,i) = temp(j,i) - dt_melt
                    !$ACC atomic update
                    !$OMP atomic update
                    temp(j,i+1) = temp(j,i+1) - dt_melt
                    !$ACC atomic update
                    !$OMP atomic update
                    temp(j+1,i) = temp(j+1,i) - dt_melt
                    !$ACC atomic update
                    !$OMP atomic update
                    temp(j+1,i+1) = temp(j+1,i+1) - dt_melt
                endif
            endif
        enddo
    enddo
    !$OMP end parallel do

    !$OMP parallel do private(tmpr, yy, depth, jj, solidus, pmelt, dt_melt)
    !$ACC parallel loop async(1)
    do i = 1, nx-1
        do j = nz-1, 1, -1
            ! flux melting in the mantel wedge occurs above serpertine or chlorite
            if (phase_ratio(kserp,j,i) + phase_ratio(khydmant,j,i) > 0.6d0 .and. &
                cord(j,i,2) > -max_melting_depth) then

                ! search the mantle above for regions above solidus
                do jj = j, 1, -1
                    tmpr = 0.25d0 * (temp(jj,i)+temp(jj,i+1)+temp(jj+1,i)+temp(jj+1,i+1))

                    ! depth below the surface in m
                    yy = 0.25d0 * (cord(jj,i,2)+cord(jj,i+1,2)+cord(jj+1,i,2)+cord(jj+1,i+1,2))
                    depth = 0.5d0*(cord(1,i,2)+cord(1,i+1,2)) - yy

                    ! Water-saturated solidus from Grove et al., Nature, 2009
                    if (depth > 80.d3) then
                        solidus = 800
                    else
                        solidus = 800 + 6.2e-8 * (depth - 80.d3)**2
                    endif
                    !write(333, *) 'chl', depth, solidus
                    if (tmpr > solidus) then
                        ! fraction of partial melting (enthalpy formulation)
                        ! F = cp * (T - T_solidus) / L_fusion, capped at 10%
                        if (latent_heat_magma > 0.d0) then
                            pmelt = min(cp(kmant1) * (tmpr - solidus) / latent_heat_magma, 0.1d0)
                        else
                            pmelt = min((tmpr - solidus) / 50.d0 * 0.1d0, 0.1d0)
                        endif
                        !$ACC atomic update
                        !$OMP atomic update
                        fmelt(jj,i) = fmelt(jj,i) + pmelt * (phase_ratio(kmant1, jj, i)  &
                                                             + phase_ratio(kmant2, jj, i) &
                                                             + phase_ratio(kserp, jj, i))
                        !print *, jj, i, tmpr, pmelt
                        
                        ! Pull temperature back to solidus (latent heat buffering)
                        ! Subtract overshoot (T_apparent - T_solidus) from element's 4 corner nodes
                        dt_melt = tmpr - solidus
                        !$ACC atomic update
                        !$OMP atomic update
                        temp(jj,i) = temp(jj,i) - dt_melt
                        !$ACC atomic update
                        !$OMP atomic update
                        temp(jj,i+1) = temp(jj,i+1) - dt_melt
                        !$ACC atomic update
                        !$OMP atomic update
                        temp(jj+1,i) = temp(jj+1,i) - dt_melt
                        !$ACC atomic update
                        !$OMP atomic update
                        temp(jj+1,i+1) = temp(jj+1,i+1) - dt_melt
                    endif
                enddo
                ! no need to look up further
                exit
            endif
        enddo
    enddo
    !$OMP end parallel do
endif

return
end subroutine change_phase
