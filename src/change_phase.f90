subroutine change_phase
!$ACC routine(newphase2marker) worker
!$ACC routine(count_phase_ratio) seq
USE marker_data
use arrays
use params
use phases

implicit none

double precision, external :: stressI,stressII, strainII, srateII, xLith_dp
double precision, external :: CalcF, CalcT, DeltaTwat, calcsolid

double precision :: yy, depth, press, &
                    tmpr, trtmpr, trpres, trpres2, ystime,&
                    solidus, pmelt, total_phase_ratio,rogh,dP

double precision :: xpress, dumFF, dumP, DTsol, dumsol, dumF, dumT, arbar, xpTpP, hP, xpTpF, AvDF, AvDT
integer :: jj, j, i, iph, jbelow, k, kinc, kk, n, m, iblk, jblk

! max. depth (m) of eclogite phase transition, no serpentinization below it
real*8, parameter :: max_basalt_depth = 150.d3
! min. temperature (C) of eclogite phase transition
real*8, parameter :: min_eclogite_temp = 400.d0
real*8, parameter :: min_eclogite_depth = 20d3
real*8, parameter :: mantle_density = 3000.d0
real*8, parameter :: max_melting_depth = 200.d3

! temperature (C) of serpentine phase transition
real*8, parameter :: serpentine_temp = 550.d0

ystime = time/365./24./3600./1.e6 

!$ACC kernels async(1)
itmp = 0  ! indicates which element has phase-changed markers
!$ACC end kernels

! Lithostatic pressure
if (itype_melting.eq.2) then
    !$OMP parallel do private(i,j,dP,rogh)
    do i = 1, nx-1
        dummye(:,i) = 0.
        rogh = 0.
        do j = 1, nz-1
            dP = xLith_dp(j,i)
            dummye(j,i) = rogh + 0.5*dP ! pressure at the middle of the element
            rogh = rogh + dP
        enddo
    enddo
endif

!$OMP parallel private(kk,i,j,k,n,tmpr,depth,iph,press, &
!$OMP                    jbelow,trpres,trpres2,kinc,yy)
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
    
    if (itype_melting.eq.2) then
        ! Lithostatic pressure
        press = dummye(j,i) ! pre-calculated pressure
        ! press = 0.
        ! rogh  = 0.
        ! do jj = 1,j
        !     densT= den(k) * (1.-alfa(k)*tmpr)
        !     dh  = 0.5*(cord(jj,i,2)-cord(jj+1,i,2) + cord(jj,i+1,2)-cord(jj+1,i+1,2))
        !     dPT = den(k) * (1.-alfa(k)*tmpr)*g*dh
        !     dP = dPT*(1.-beta(k)*rogh)/(1.+beta(k)/2.*dPT)
        !     press = rogh + 0.5*dP
        !     rogh = rogh + dP
        ! enddo
    else
        press = mantle_density * g * depth
    endif

    ! # of markers inside quad
    kinc = nmark_elem(j,i)

    ! If location of this element is too deep, this marker is already
    ! too deep in the mantle, where there is no significant phase change.
    if (depth > 200.d3) cycle

    iph = mark_phase(kk)

    ! Rules of phase changes
    select case(iph)
    case (kcont1, kcont2)
        ! XXX: middle crust with high dissipation becomes weaker,
        ! this helps with localization
        if (itype_melting.eq.2) then
            if(i.ge.10.and.i.le.490) then
            if(depth.le.20.e3) then
            if(tmpr >= 250. .and. tmpr <= 300.) then
            if(stressII(j,i)*strainII(j,i) > 1.e7) then   ! 1.e7
                !$ACC atomic write
                !$OMP atomic write
                itmp(j,i) = 1
                mark_phase(kk) = kweakmc
            endif
            endif
            endif
            endif
        endif            
    case (kmant1,kmant2)
        ! subuducted oceanic crust below mantle, mantle is serpentinized
        ! if(depth > max_basalt_depth) cycle

        ! Phase diagram taken from Ulmer and Trommsdorff, Science, 1995
        ! Fixed points (730 C, 2.1 GPa) (500 C, 7.5 GPa)
        trpres = 2.1d9 + (7.5d9 - 2.1d9) * (tmpr - 730.d0) / (500.d0 - 730.d0)
        ! Fixed points (730 C, 2.1 GPa) (670 C, 0.6 GPa)
        trpres2 = 2.1d9 + (0.6d9 - 2.1d9) * (tmpr - 730.d0) / (670.d0 - 730.d0)
        ! press = mantle_density * g * depth
        if (.not. (press < trpres .and. press > trpres2)) cycle
        do jbelow = min(j+1,nz-1), min(j+nelem_serp,nz-1)
            if(phase_ratio(kocean1,jbelow,i) > 0.8d0 .or. &
                phase_ratio(kocean2,jbelow,i) > 0.8d0 .or. &
                phase_ratio(ksed1,jbelow,i) > 0.8d0) then
                !$ACC atomic write
                !$OMP atomic write
                itmp(j,i) = 1
                mark_phase(kk) = kserp
                exit
            endif
        enddo

        if (itype_melting.eq.2) then
            ! after 1 Myr, form mantle detachments
            if(ystime.gt.1.) then
            if(tmpr.gt.800. .and. tmpr.lt.1000. .and. &
                stressII(j,i)*srateII(j,i).ge.200.e-8) then
                !$ACC atomic write
                !$OMP atomic write
                itmp(j,i) = 1
                mark_phase(kk) = khtmsz 
            endif
            endif
        endif
    case (kocean0, kocean1, kocean2, ksills)
        ! basalt -> eclogite
        ! phase change pressure
        ! Phase Diagram taken from Hacker, JGR, 2003 (Figure 8 or Figure 1)
        ! Fixed points (400 C, 5.1 GPa) (536 C, 0 GPa)
        trpres = -0.3d9 + 2.2d6*tmpr
        trpres2 = 20.1d9 - 0.0375d9*tmpr
        ! press = mantle_density * g * depth
        if (tmpr < min_eclogite_temp .or. depth < min_eclogite_depth .or. &
            press < trpres .or. press < trpres2) cycle
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
        ! press = mantle_density * g * depth
        if (tmpr < serpentine_temp .or. (press < trpres .and. press > trpres2)) cycle
        !$ACC atomic write
        !$OMP atomic write
        itmp(j,i) = 1
        mark_phase(kk) = kmant1
    case (ksed2)
        ! compaction, uncosolidated sediment -> sedimentary rock
        if (tmpr < 80d0 .or. depth < 2d3) cycle
        !$ACC atomic write
        !$OMP atomic write
        itmp(j,i) = 1
        mark_phase(kk) = ksed1
    case (ksed1)
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
        ! press = mantle_density * g * depth
        if (press < trpres ) cycle
        !$ACC atomic write
        !$OMP atomic write
        itmp(j,i) = 1
        mark_phase(kk) = kschist
    case (khydmant)
        ! dehydration of chlorite
        ! Phase diagram from Grove et al. Nature, 2009
        trtmpr = 880 - 35d-9 * (depth - 62d3)**2
        if (tmpr < trtmpr) cycle
        !$ACC atomic write
        !$OMP atomic write
        itmp(j,i) = 1
        mark_phase(kk) = kmant1        
    case (khtmsz)
        if (itype_melting.eq.2) then
            if (tmpr.le.serpentine_temp) then
            if (aps(j,i).ge.0.1) then
                !$ACC atomic write
                !$OMP atomic write
                itmp(j,i) = 1
                mark_phase(kk) = kserp
            endif
            endif
        endif
    end select

enddo
!$OMP end do
!$OMP end parallel

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

if (itype_melting .eq. 1) then
    !$OMP parallel do private(tmpr, yy, depth, solidus, pmelt, total_phase_ratio)
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
                    ! fraction of partial melting
                    ! 10% of melting at solidus + 50 C
                    ! Hirschmann, 2000 G3.
                    pmelt = min((tmpr - solidus) / 50 * 0.1d0, 0.1d0)
                    fmelt(j,i) = pmelt * total_phase_ratio
                    !print *, j, i, tmpr, pmelt
                endif
            endif
        enddo
    enddo
    !$OMP end parallel do

    !$OMP parallel do private(tmpr, yy, depth, solidus, pmelt, total_phase_ratio, press)
    !$ACC parallel loop collapse(2) async(1)
    do i = 1, nx-1
        do j = 1, nz-1

            ! basalt and eclogite rock melting
            ! solidus from Gutscher, 2000 Geology
            total_phase_ratio = phase_ratio(kocean1,j,i) + phase_ratio(kocean2,j,i) &
                                + phase_ratio(kocean0,j,i) + phase_ratio(keclg,j,i)
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
                    ! fraction of partial melting
                    ! XXX: assuming 10% of melting at solidus + 20 C
                    pmelt = min((tmpr - solidus) / 20 * 0.1d0, 0.1d0)
                    !$ACC atomic update
                    !$OMP atomic update
                    fmelt(j,i) = fmelt(j,i) + pmelt * total_phase_ratio
                endif
            endif
        enddo
    enddo
    !$OMP end parallel do

    !$OMP parallel do private(tmpr, yy, depth, jj, solidus, pmelt)
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
                    if (tmpr > solidus) then
                        ! fraction of partial melting
                        ! 10% of melting at solidus + 50 C
                        ! Hirschmann, 2000 G3.
                        pmelt = min((tmpr - solidus) / 50 * 0.1d0, 0.1d0)
                        !$ACC atomic update
                        !$OMP atomic update
                        fmelt(jj,i) = fmelt(jj,i) + pmelt * (phase_ratio(kmant1, jj, i)  &
                                                             + phase_ratio(kmant2, jj, i) &
                                                             + phase_ratio(kserp, jj, i))
                        !print *, jj, i, tmpr, pmelt
                    endif
                enddo
                ! no need to look up further
                exit
            endif
        enddo
    enddo
    !$OMP end parallel do

elseif (itype_melting .eq. 2) then

    !Update Solidus and Liquidus
    call calc_solidus

    !===========================================
    !Nodes to Elements, T and P calc, and Initialization
    !===========================================
    !$OMP parallel do private(i,j,xpress)
    do j = 1,nz-1-10
        do i = 1+5,nx-1-5
            avT(j,i) = 0.25*(temp(j,i) + temp(j,i+1) + temp(j+1,i) + temp(j+1,i+1))
            !  Lithostatic pressure
            xpress = -stressI(j,i)/1.e9 ! Pressure in GPa   
            
            !===========================
            ! Initialize
            !===========================
            if(nloop.eq.1) then
                Eff_melt(j,i) = 0.

                if (phase_ratio(kmant2,j,i).gt.0.95) then  
                    ! Calculate initial melt
                    Eff_melt(j,i) = CalcF(xpress,avT(j,i),watercont,xmodalcpx)
                    fmelt(j,i) = Eff_melt(j,i) ! Only for the first time step. 
                end if
            end if
            
            !=========================================
            ! Initialize upper and lower crustal melt
            !=========================================
            ! Upper and lower crust solidus (Gerya p372)
            ! if (icrust_melt.eq.1) then
            !     if (phase_ratio(kcont1,j,i).gt.0.95) then
            !         Eff_melt_lc(j,i) = 0.0
            !     endif
            ! endif
        end do
    end do
    !$OMP end parallel do


    !=================================
    ! Calculate changes in melt and temps
    !=================================
    do iblk = 1,2
    do jblk = 1,2
    !$OMP parallel do private(i,j,m,xpress,dumFF,dumP,DTsol,dumsol,hP,dumF,dumT, &
    !$OMP                     arbar,xpTpP,xDF,xDT,AvDF,AvDT)
    do j = jblk,nz-1-34,2
        do i = iblk+10,nx-1-10,2
            xpress = -stressI(j,i)/1.e9 ! Pressure in GPa   
        
            if (phase_ratio(kmant2,j,i).gt.0.95.and.nloop.ne.1.and.fmagma(j,i).eq.0.) then         
                dumFF = 0.
                dumP = xpress
                DTsol = deltaTwat(watercont,dumP,dumFF)
                dumsol = tsol(j,i) 

                if (avT(j,i).gt.(dumsol-DTsol)) then
                        
                    hP = (zpressm(j,i)/1.e9) - xpress !Amount of decompression in GPa
                    if (hP.gt.1e-12) then
                    
                        do m = 1,4 !4th order solver
                            dumP = (zpressm(j,i)/1.e9) - hP*dumC(m)
            
                            if (m.eq.1) then
                                dumF = Eff_melt(j,i) 
                                dumT = avT(j,i)
                            else
                                dumF = Eff_melt(j,i) - dumC(m)*hP*xDF(m-1) !Increasing F
                                dumT = avT(j,i) - dumC(m)*hP*xDT(m-1) ! Decreasing T
                            end if

                            arbar = dumF*(68./2.9) + (1.-dumF)*(40./3.3)

                            ! Eqn 22 Katz 2003
                            xpTpP = ((calcT((dumP+(1.e-14)),watercont,xmodalcpx,dumF)) &
                                -(calcT((dumP-(1.e-14)),watercont,xmodalcpx,dumF))) &
                                /(2.e-14)                               

                            dumFF = 0.
                            DTsol = deltaTwat(watercont,dumP,dumFF) 
                            dumsol = calcsolid(dumP)
                            
                            if (dumT.gt.(dumsol-DTsol)) then
                                dumF = dmax1(dumF,1.e-9) !Work out from a reasonable value

                                ! Eqn 21 Katz 2003
                                xpTpF = ((calcT(dumP,watercont,xmodalcpx,(dumF+(1.e-14)))) &
                                    -(calcT(dumP,watercont,xmodalcpx,(dumF-(1.e-14))))) &
                                    /(2.e-14)

                                ! Melt Derivative - dFdP - Eqn 20 Katz2003
                                xDF(m) = (-1000.0*xpTpP/(dumT+273.0) + arbar)/(1000.0*xpTpF/ &
                                    (dumT+273.0)+ 300.0)
                            else
                                xDF(m) = 0.
                            end if

                            ! Temp derivative - dTdP - Eqn 23 Katz 2003
                            xDT(m) = (dumT+273.0)*(arbar - 300.0*xDF(m))/1000.0 
                        end do
                
                        AvDF = (xDF(1)/6. + xDF(2)/3. + xDF(3)/3. + xDF(4)/6.)
                        AvDT = (xDT(1)/6. + xDT(2)/3. + xDT(3)/3. + xDT(4)/6.)

                        if (AvDF.le.0.) then
                                
                            Eff_melt(j,i) = Eff_melt(j,i) + abs(AvDF*hP) ! Keep 5% melt in asthenosphere
                            fmelt(j,i) =  abs(AvDF*hP) ! 95% is migrating?
                            if (Eff_melt(j,i).ge.1.) then
                                Eff_melt(j,i) = 0.99 
                                fmelt(j,i) = 0.
                                AvDT = 0.0
                            end if

                            ! Proportionally distribute the element change in temperature to the nodes
                            if (AvDF.ne.0.) then
                                deltaTLH(j,i)= AvDT*hP
                                avT(j,i) = avT(j,i) - deltaTLH(j,i)
                                temp(j,i) = temp(j,i) - deltaTLH(j,i)*(temp(j,i)/avT(j,i)/4.)
                                temp(j,i+1) = temp(j,i+1) - deltaTLH(j,i)*(temp(j,i+1)/avT(j,i)/4.)
                                temp(j+1,i) = temp(j+1,i) - deltaTLH(j,i)*(temp(j+1,i)/avT(j,i)/4.)
                                temp(j+1,i+1) = temp(j+1,i+1) - deltaTLH(j,i)*(temp(j+1,i+1)/avT(j,i)/4.)
                            end if
                        end if

                    end if
                end if
            end if
        end do
    end do
    !$OMP end parallel do
    end do
    end do
endif

return
end subroutine change_phase


!==========================================================
! Solidus Subroutine
!==========================================================
! Four different solidus may be chosed from for melt production models
subroutine calc_solidus
use arrays
use params
implicit none
integer :: i,j
double precision :: c1,c2,c3,c4,Tolm
double precision :: stressone,xPs,dPdTm,tint,Tsm,deltam
double precision, external :: stressI

!$OMP Parallel private(i,j,stressone,xPs,dPdTm,tint,Tsm,deltam)

!==========================================================
! Calculate liquidus
!==========================================================

! McKenzie and Bickle liquidus
!

if (isolidus.eq.4) then
    !$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            stressone = -stressI(j,i)/1.e9
            tliq(j,i) =  1736.2 + 4.343*stressone + 180.*atan(stressone/2.2169)
        end do
    end do
    !$OMP end do

! Katz 2003 DEFAULT
else

    !$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            stressone = -stressI(j,i)/1.e9
            tliq(j,i) =  1780. + 45.*stressone - 2.*(stressone**2)
        end do
    end do 
    !$OMP end do
end if

if (icrust_melt.eq.1) then
! Liquidus for upper and lower crust
    !$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            stressone = -stressI(j,i)/1.e9
            tliquc(j,i) =  (1262. + 0.09*stressone*1000.) - 273.
            tliqlc(j,i) =  (1423 +  0.105*stressone*1000.) - 273.
        end do
    end do 
    !$OMP end do
endif

!==========================================================
! Calculate Solidus
!==========================================================

! Upper and lower crust solidus (Gerya p372)
    if (icrust_melt.eq.1) then
    !$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            ! Change pressure to GPa
            stressone = -stressI(j,i)/1.e6
            if (stressone.le.1.2) then
                tsoluc(j,i) = (889.+17900./(stressone+54.)+20200./(stressone+54.)**2.)-273.
            else
                tsoluc(j,i) = (831.+0.06*stressone)-273.
            end if
            if (stressone.le.1.6) then
                tsollc(j,i) = (973.+70400./(stressone+354.)+77800000./(stressone+354.)**2.)-273. 
            else
                tsollc(j,i) = (935.+0.0035*stressone+0.0000062*stressone*stressone)-273.
            end if
        end do
    end do
    !$OMP end do
endif

! Katz et al (2003) DEFAULT
if (isolidus.eq.1) then
    !$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            stressone = -stressI(j,i)/1.e9
            tsol(j,i) = (-5.104)*(stressone**2) + 132.899*(stressone) + 1085.7
        end do
    end do
    !$OMP end do

! Hirschmann et al (2000) & Hirschmann et al (2009) 
else if (isolidus.eq.2) then
    !$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            ! Change pressure to GPa
            stressone = -stressI(j,i)/1.e9
            if (stressone.le.10.) then
                tsol(j,i) = -5.140*stressone*stressone + 132.899*stressone + 1120.661 
            elseif (stressone.le.23.5) then
                tsol(j,i) = -1.092*((stressone-10.)**2) + 32.390*(stressone-10.) + 1935.
            else
                tsol(j,i) = 26.35*(stressone-23.5) + 2175
            end if
        end do
    end do
!$OMP end do
! Herzberg et al (2000) 
else if (isolidus.eq.3) then

!$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            stressone = -stressI(j,i)/1.e9
            if (stressone.le.2.7) then
                tsol(j,i) = 132*(stressone) + 1102
            else
                tsol(j,i) = 390*(log(stressone)) - 5.7*(stressone) + 1086
            end if
        end do
    end do
!$OMP end do

! McKenzie & Bickle (1988) Solidus 
else if (isolidus.eq.4) then
    c1 = 1100.
    c2 = 136.
    c3 = 0.0004968
    c4 = 0.012
    Tolm = 0.01

    !$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            stressone = -stressI(j,i)/1.e9
            tsol(j,i) = 1000. 
            xPs = ((tsol(j,i)-c1)/c2)+(c3*exp(c4*(tsol(j,i)-c1)))
            dPdTm = (c3*c4*exp(c4*(tsol(j,i)-c1)) + (1/c2))
            tint = xPs - (dPdTm*tsol(j,i))
            Tsm = (stressone-tint)/dPdTm
            deltam = abs(Tsm-tsol(j,i))
            tsol(j,i) = Tsm
            
            do while (deltam.ge.Tolm)
                xPs = ((tsol(j,i)-c1)/c2)+(c3*exp(c4*(tsol(j,i)-c1)))
                dPdTm= (c3*c4*exp(c4*(tsol(j,i)-c1)))+(1/c2)
                tint = xPs - (dPdTm*tsol(j,i))
                Tsm = (stressone-tint)/dPdTm
                deltam = abs(Tsm-tsol(j,i))
                tsol(j,i) = Tsm
            end do
        end do
    end do
    !$OMP end do
endif
!$OMP end parallel 

return
end subroutine calc_solidus
