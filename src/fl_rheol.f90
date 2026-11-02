!  Rheology (Update stresses depending on rheology)
!  Calculate total finite strain and plastic strain  
    
subroutine fl_rheol
!$ACC routine(pre_plast) seq
!$ACC routine(elastic) seq
!$ACC routine(maxwell) seq
!$ACC routine(plastic) seq
!$ACC routine(Eff_visc) seq
use arrays
use params
implicit none

double precision :: depl(4)
double precision :: s11p(4),s22p(4),s12p(4),s33p(4),s11v(4),s22v(4),s12v(4),s33v(4)
double precision :: bulkm,rmu,coh,phi,psi, &
                    stherm,hardn,vis, &
                    de11,de22,de12,de33,dv, &
                    diss, poiss, &
                    quad_area, s0, s0a,s0b, &
                    sII_plas, sII_visc, young
double precision :: vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
                    x1,y1,x2,y2,x3,y3,x4,y4, &
                    em,eda,edb, &
                    s11,s22,s12,srII,srI,srs2,stII, &
                    ar1,ar2,ar3,ar4
double precision :: shpdx_loc(3, 4), shpdz_loc(3, 4), sr_loc(3, 4)
double precision :: s11_start(4), s22_start(4), s12_start(4), s33_start(4)
double precision :: s11_final, s22_final, s12_final, s33_final
double precision :: tr_s_start, tr_s_final
double precision :: ee11_start, ee22_start, ee33_start, ee12_start
double precision :: ee11_final, ee22_final, ee33_final, ee12_final
double precision :: de11_el, de22_el, de33_el, de12_el
double precision :: de11_ne, de22_ne, de33_ne, de12_ne
double precision :: s11_avg, s22_avg, s12_avg, s33_avg
double precision :: q_shear_k, q_shear_elem
double precision :: lambda, mu, denom_elastic
double precision :: Eff_visc
integer :: i, j, k, iph, irh, &
           ipls

! Increment average time
dtavg = dtavg + dt
!$ACC update device(dtavg) async(2)

!$OMP Parallel Private(i,j,k,iph,irh,bulkm,rmu,coh,phi,psi, &
!$OMP                  stherm,hardn,vis, &
!$OMP                  de11,de22,de12,de33,dv, &
!$OMP                  s11p,s22p,s12p,s33p, &
!$OMP                  s11v,s22v,s12v,s33v, &
!$OMP                  depl,ipls,diss, &
!$OMP                  sII_plas,sII_visc, &
!$OMP                  quad_area,s0a,s0b,s0, &
!$OMP                  vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
!$OMP                  x1,y1,x2,y2,x3,y3,x4,y4, &
!$OMP                  em,eda,edb,s11,s22,s12,srII,srI,srs2,stII, &
!$OMP                  ar1,ar2,ar3,ar4,shpdx_loc,shpdz_loc,sr_loc, &
!$OMP                  s11_start,s22_start,s12_start,s33_start, &
!$OMP                  s11_final,s22_final,s12_final,s33_final, &
!$OMP                  tr_s_start,tr_s_final, &
!$OMP                  ee11_start,ee22_start,ee33_start,ee12_start, &
!$OMP                  ee11_final,ee22_final,ee33_final,ee12_final, &
!$OMP                  de11_el,de22_el,de33_el,de12_el, &
!$OMP                  de11_ne,de22_ne,de33_ne,de12_ne, &
!$OMP                  s11_avg,s22_avg,s12_avg,s33_avg, &
!$OMP                  q_shear_k,q_shear_elem, &
!$OMP                  lambda,mu,denom_elastic)
!$OMP do schedule(guided)
!$ACC parallel loop gang vector collapse(2) private(depl,s11p,s22p,s12p,s33p,s11v,s22v,s12v,s33v,shpdx_loc,shpdz_loc,sr_loc,s11_start,s22_start,s12_start,s33_start) async(1)
do i = 1,nx-1
    do j = 1,nz-1

        ! Node velocities:
        vx1 = vel(j  ,i  ,1)
        vy1 = vel(j  ,i  ,2)
        vx2 = vel(j+1,i  ,1)
        vy2 = vel(j+1,i  ,2)
        vx3 = vel(j  ,i+1,1)
        vy3 = vel(j  ,i+1,2)
        vx4 = vel(j+1,i+1,1)
        vy4 = vel(j+1,i+1,2)

        ! Node coordinates:
        x1 = cord(j  ,i  ,1)
        y1 = cord(j  ,i  ,2)
        x2 = cord(j+1,i  ,1)
        y2 = cord(j+1,i  ,2)
        x3 = cord(j  ,i+1,1)
        y3 = cord(j  ,i+1,2)
        x4 = cord(j+1,i+1,1)
        y4 = cord(j+1,i+1,2)

        ar1 = area(j, i, 1)
        ar2 = area(j, i, 2)
        ar3 = area(j, i, 3)
        ar4 = area(j, i, 4)

        ! Triangle A (k=1)
        shpdx_loc(1, 1) = (y2 - y3) * ar1
        shpdx_loc(2, 1) = (y3 - y1) * ar1
        shpdx_loc(3, 1) = (y1 - y2) * ar1

        shpdz_loc(1, 1) = (x3 - x2) * ar1
        shpdz_loc(2, 1) = (x1 - x3) * ar1
        shpdz_loc(3, 1) = (x2 - x1) * ar1

        ! Triangle B (k=2)
        shpdx_loc(1, 2) = (y2 - y4) * ar2
        shpdx_loc(2, 2) = (y4 - y3) * ar2
        shpdx_loc(3, 2) = (y3 - y2) * ar2

        shpdz_loc(1, 2) = (x4 - x2) * ar2
        shpdz_loc(2, 2) = (x3 - x4) * ar2
        shpdz_loc(3, 2) = (x2 - x3) * ar2

        ! Triangle C (k=3)
        shpdx_loc(1, 3) = (y2 - y4) * ar3
        shpdx_loc(2, 3) = (y4 - y1) * ar3
        shpdx_loc(3, 3) = (y1 - y2) * ar3

        shpdz_loc(1, 3) = (x4 - x2) * ar3
        shpdz_loc(2, 3) = (x1 - x4) * ar3
        shpdz_loc(3, 3) = (x2 - x1) * ar3

        ! Triangle D (k=4)
        shpdx_loc(1, 4) = (y4 - y3) * ar4
        shpdx_loc(2, 4) = (y3 - y1) * ar4
        shpdx_loc(3, 4) = (y1 - y4) * ar4

        shpdz_loc(1, 4) = (x3 - x4) * ar4
        shpdz_loc(2, 4) = (x1 - x3) * ar4
        shpdz_loc(3, 4) = (x4 - x1) * ar4

        ! (1) A element:
        sr_loc(1,1) = vx1 * shpdx_loc(1, 1) + vx2 * shpdx_loc(2, 1) + vx3 * shpdx_loc(3, 1)
        sr_loc(2,1) = vy1 * shpdz_loc(1, 1) + vy2 * shpdz_loc(2, 1) + vy3 * shpdz_loc(3, 1)
        sr_loc(3,1) = 0.5d0*(vx1 * shpdz_loc(1, 1) + vx2 * shpdz_loc(2, 1) + vx3 * shpdz_loc(3, 1) + &
                             vy1 * shpdx_loc(1, 1) + vy2 * shpdx_loc(2, 1) + vy3 * shpdx_loc(3, 1))

        !(2) B element:
        sr_loc(1,2) = vx3 * shpdx_loc(1, 2) + vx2 * shpdx_loc(2, 2) + vx4 * shpdx_loc(3, 2)
        sr_loc(2,2) = vy3 * shpdz_loc(1, 2) + vy2 * shpdz_loc(2, 2) + vy4 * shpdz_loc(3, 2)
        sr_loc(3,2) = 0.5d0*(vx3 * shpdz_loc(1, 2) + vx2 * shpdz_loc(2, 2) + vx4 * shpdz_loc(3, 2) + &
                             vy3 * shpdx_loc(1, 2) + vy2 * shpdx_loc(2, 2) + vy4 * shpdx_loc(3, 2))

        ! (3) C element:
        sr_loc(1,3) = vx1 * shpdx_loc(1, 3) + vx2 * shpdx_loc(2, 3) + vx4 * shpdx_loc(3, 3)
        sr_loc(2,3) = vy1 * shpdz_loc(1, 3) + vy2 * shpdz_loc(2, 3) + vy4 * shpdz_loc(3, 3)
        sr_loc(3,3) = 0.5d0*(vx1 * shpdz_loc(1, 3) + vx2 * shpdz_loc(2, 3) + vx4 * shpdz_loc(3, 3) + &
                             vy1 * shpdx_loc(1, 3) + vy2 * shpdx_loc(2, 3) + vy4 * shpdx_loc(3, 3))

        ! (4) D element:
        sr_loc(1,4) = vx1 * shpdx_loc(1, 4) + vx4 * shpdx_loc(2, 4) + vx3 * shpdx_loc(3, 4)
        sr_loc(2,4) = vy1 * shpdz_loc(1, 4) + vy4 * shpdz_loc(2, 4) + vy3 * shpdz_loc(3, 4)
        sr_loc(3,4) = 0.5d0*(vx1 * shpdz_loc(1, 4) + vx4 * shpdz_loc(2, 4) + vx3 * shpdz_loc(3, 4) + &
                             vy1 * shpdx_loc(1, 4) + vy4 * shpdx_loc(2, 4) + vy3 * shpdx_loc(3, 4))

        if ( mix_strain .eq. 1 ) then
            do k = 1, 3, 2
                em = 0.5d0*(sr_loc(1,k)+sr_loc(2,k)+sr_loc(1,k+1)+sr_loc(2,k+1))
                eda = sr_loc(1,k)-sr_loc(2,k)
                edb = sr_loc(1,k+1)-sr_loc(2,k+1)
                sr_loc(1,k) = 0.5d0*(em+eda)
                sr_loc(2,k) = 0.5d0*(em-eda)
                sr_loc(1,k+1) = 0.5d0*(em+edb)
                sr_loc(2,k+1) = 0.5d0*(em-edb)
            enddo
        endif

        ! Integration for averaging of strain rate and dissipation function
        s11 = 0.25d0 * sum(sr_loc(1,:))
        s22 = 0.25d0 * sum(sr_loc(2,:))
        s12 = 0.25d0 * sum(sr_loc(3,:))
        srII = 0.5d0 * sqrt((s11-s22)**2 + 4.d0*s12*s12)
        srI = 0.5d0 * (s11+s22)
        srs2 = (s11-srI)*(s11-srI) + (s22-srI)*(s22-srI) + 2.d0*s12*s12
        se2sr(j,i,1) = se2sr(j,i,1) + s11*dt
        se2sr(j,i,2) = se2sr(j,i,2) + s22*dt
        se2sr(j,i,3) = se2sr(j,i,3) + s12*dt

        s11 = 0.25d0 * sum(stress0(j,i,1,:))
        s22 = 0.25d0 * sum(stress0(j,i,2,:))
        s12 = 0.25d0 * sum(stress0(j,i,3,:))
        stII = 0.5d0 * sqrt((s11-s22)**2 + 4.d0*s12*s12)
        if( srII .ne. 0.d0 ) sshrheat(j,i) = sshrheat(j,i) + stII/srII*srs2*dt

        ! iphase (j,i) is number of a phase NOT a rheology
        iph = iphase(j,i)
        irh = irheol(iph)

        ! Elastic modules & viscosity & plastic properties
        bulkm = rl(iph) + 2.d0*rm(iph)/3.d0
        rmu   = rm(iph)

        ! Thermal stresses
        stherm = 0.d0
        if (istress_therm .gt. 0) stherm = -alfa(iph)*bulkm*(temp(j,i)-temp0(j,i))

        ! Preparation of plastic properties
        if (irh .eq. 6 .or. irh .ge. 11) call pre_plast(i,j,coh,phi,psi,hardn)
              
        ! Re-evaluate viscosity
        if (irh .eq. 3 .or. irh .ge. 11) then 
            if( mod(nloop,ifreq_visc) .eq. 0) visn(j,i) = Eff_visc(j,i)
        endif
        vis = visn(j,i)

        do k = 1,4
            s11_start(k) = stress0(j,i,1,k)
            s22_start(k) = stress0(j,i,2,k)
            s12_start(k) = stress0(j,i,3,k)
            s33_start(k) = stress0(j,i,4,k)
        enddo

        ! Cycle by triangles
        do k = 1,4
            de11 = sr_loc(1,k)*dt
            de22 = sr_loc(2,k)*dt
            de12 = sr_loc(3,k)*dt
            de33 = 0.d0
            ! Dike-intrusion eigenstrain (arc melting pathway, see
            ! fl_therm.f90): freshly-released dike volume dilates this
            ! element beyond what its node geometry (dvol) alone shows,
            ! the same trick already used for thermal-expansion stress
            ! (stherm) above -- the constitutive law is told the element
            ! dilated by this much extra, relaxing/lowering stress, and
            ! the actual node displacement then emerges over subsequent
            ! timesteps via the normal force-balance solve.
            dv = dvol(j,i,k) - dike_released(j,i)
            s11p(k) = stress0(j,i,1,k) + stherm 
            s22p(k) = stress0(j,i,2,k) + stherm 
            s12p(k) = stress0(j,i,3,k) 
            s33p(k) = stress0(j,i,4,k) + stherm
            s11v(k) = s11p(k)
            s22v(k) = s22p(k)
            s12v(k) = s12p(k)
            s33v(k) = s33p(k)

            if (irh .eq. 1) then
                call elastic(bulkm,rmu,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de12)
                stress0(j,i,1,k) = s11p(k)
                stress0(j,i,2,k) = s22p(k)
                stress0(j,i,3,k) = s12p(k)
                stress0(j,i,4,k) = s33p(k)

            elseif (irh .eq. 3) then
                call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),de11,de22,de33,de12,dv,ndim,dt)
                stress0(j,i,1,k) = s11v(k)
                stress0(j,i,2,k) = s22v(k)
                stress0(j,i,3,k) = s12v(k)
                stress0(j,i,4,k) = s33v(k)

            elseif (irh .eq. 6) then
                call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,ten_off,ndim)
                stress0(j,i,1,k) = s11p(k)
                stress0(j,i,2,k) = s22p(k)
                stress0(j,i,3,k) = s12p(k)
                stress0(j,i,4,k) = s33p(k)

            elseif (irh .ge. 11) then 
                call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,ten_off,ndim)
                call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),de11,de22,de33,de12,dv,ndim,dt)
            endif
        enddo

        if( irh .ge. 11 ) then
            sII_plas = (s11p(1)+s11p(2)+s11p(3)+s11p(4)-s22p(1)-s22p(2)-s22p(3)-s22p(4))**2 &
                     + 4.d0*(s12p(1)+s12p(2)+s12p(3)+s12p(4))**2

            sII_visc = (s11v(1)+s11v(2)+s11v(3)+s11v(4)-s22v(1)-s22v(2)-s22v(3)-s22v(4))**2 &
                     + 4.d0*(s12v(1)+s12v(2)+s12v(3)+s12v(4))**2

            if (sII_plas .lt. sII_visc) then
                do k = 1, 4
                    stress0(j,i,1,k) = s11p(k)
                    stress0(j,i,2,k) = s22p(k)
                    stress0(j,i,3,k) = s12p(k)
                    stress0(j,i,4,k) = s33p(k)
                end do
            else 
                do k = 1, 4
                    stress0(j,i,1,k) = s11v(k)
                    stress0(j,i,2,k) = s22v(k)
                    stress0(j,i,3,k) = s12v(k)
                    stress0(j,i,4,k) = s33v(k)
                end do
            endif
        endif

        ! Averaging of isotropic stresses for pair of elements
        if (mix_stress .eq. 1 ) then
            quad_area = 1.d0/ar1 + 1.d0/ar2
            s0a=0.5d0*(stress0(j,i,1,1)+stress0(j,i,2,1))
            s0b=0.5d0*(stress0(j,i,1,2)+stress0(j,i,2,2))
            s0=(s0a/ar2+s0b/ar1)/quad_area
            stress0(j,i,1,1) = stress0(j,i,1,1) - s0a + s0
            stress0(j,i,2,1) = stress0(j,i,2,1) - s0a + s0
            stress0(j,i,1,2) = stress0(j,i,1,2) - s0b + s0
            stress0(j,i,2,2) = stress0(j,i,2,2) - s0b + s0

            quad_area = 1.d0/ar3 + 1.d0/ar4
            s0a=0.5d0*(stress0(j,i,1,3)+stress0(j,i,2,3))
            s0b=0.5d0*(stress0(j,i,1,4)+stress0(j,i,2,4))
            s0=(s0a/ar4+s0b/ar3)/quad_area
            stress0(j,i,1,3) = stress0(j,i,1,3) - s0a + s0
            stress0(j,i,2,3) = stress0(j,i,2,3) - s0a + s0
            stress0(j,i,1,4) = stress0(j,i,1,4) - s0b + s0
            stress0(j,i,2,4) = stress0(j,i,2,4) - s0b + s0
        endif

        ! Remove thermal stress
        do k = 1,4
            stress0(j,i,1,k) = stress0(j,i,1,k) - stherm 
            stress0(j,i,2,k) = stress0(j,i,2,k) - stherm 
            stress0(j,i,4,k) = stress0(j,i,4,k) - stherm 
        enddo

        ! ACCUMULATED PLASTIC STRAIN (SUMMATION OVER TRIANGLES)
        if (irh .eq. 6 .or. irh .ge. 11) then
            quad_area = ar1 + ar2 + ar3 + ar4

            aps(j,i) = aps(j,i) + (depl(1)*ar1 + depl(2)*ar2 + depl(3)*ar3 + depl(4)*ar4) / quad_area 

            if( aps(j,i) .lt. 0.d0 ) aps(j,i) = 0.d0

            if (tau_heal .ne. 0.d0) &
                 aps (j,i) = aps (j,i)/(1.d0+dt/tau_heal)
        end if

        ! TOTAL FINITE STRAIN (Reuses s11, s22, s12)
        strain(j,i,1) = strain(j,i,1) + s11*dt
        strain(j,i,2) = strain(j,i,2) + s22*dt
        strain(j,i,3) = strain(j,i,3) + s12*dt

        ! Non-elastic shear heating calculation
        if (ishearh .eq. 1) then
            q_shear_elem = 0.d0
            do k = 1, 4
                s11_final = stress0(j,i,1,k)
                s22_final = stress0(j,i,2,k)
                s12_final = stress0(j,i,3,k)
                s33_final = stress0(j,i,4,k)

                lambda = bulkm - 2.d0/3.d0 * rmu
                mu = rmu

                denom_elastic = 3.d0 * lambda + 2.d0 * mu
                if (abs(denom_elastic) .lt. 1.d-30) denom_elastic = sign(1.d-30, denom_elastic)

                ! Elastic strain at start of step
                tr_s_start = s11_start(k) + s22_start(k) + s33_start(k)
                ee11_start = (s11_start(k) - (lambda / denom_elastic) * tr_s_start) / (2.d0 * mu)
                ee22_start = (s22_start(k) - (lambda / denom_elastic) * tr_s_start) / (2.d0 * mu)
                ee33_start = (s33_start(k) - (lambda / denom_elastic) * tr_s_start) / (2.d0 * mu)
                ee12_start = s12_start(k) / (2.d0 * mu)

                ! Elastic strain at end of step
                tr_s_final = s11_final + s22_final + s33_final
                ee11_final = (s11_final - (lambda / denom_elastic) * tr_s_final) / (2.d0 * mu)
                ee22_final = (s22_final - (lambda / denom_elastic) * tr_s_final) / (2.d0 * mu)
                ee33_final = (s33_final - (lambda / denom_elastic) * tr_s_final) / (2.d0 * mu)
                ee12_final = s12_final / (2.d0 * mu)

                ! Incremental elastic strain of this step
                de11_el = ee11_final - ee11_start
                de22_el = ee22_final - ee22_start
                de33_el = ee33_final - ee33_start
                de12_el = ee12_final - ee12_start

                ! Non-elastic strain increment
                de11_ne = sr_loc(1,k)*dt - de11_el
                de22_ne = sr_loc(2,k)*dt - de22_el
                de12_ne = sr_loc(3,k)*dt - de12_el
                de33_ne = 0.d0 - de33_el

                s11_avg = 0.5d0 * (s11_start(k) + s11_final)
                s22_avg = 0.5d0 * (s22_start(k) + s22_final)
                s12_avg = 0.5d0 * (s12_start(k) + s12_final)
                s33_avg = 0.5d0 * (s33_start(k) + s33_final)

                q_shear_k = s11_avg * de11_ne + s22_avg * de22_ne + 2.d0 * s12_avg * de12_ne + s33_avg * de33_ne
                q_shear_elem = q_shear_elem + 0.25d0 * q_shear_k
            enddo

            sshrheat(j,i) = sshrheat(j,i) + q_shear_elem
        endif

    enddo
enddo
!$OMP end do
!$OMP end parallel

! Periodic time-averaging update
!$ACC wait(2)
if( nsrate .eq. ifreq_avgsr ) then
    !$OMP parallel do
    !$ACC parallel loop collapse(2) async(1)
    do i = 1,nx-1
        do j = 1, nz-1
            e2sr(j,i) = (0.5d0 * sqrt((se2sr(j,i,1)-se2sr(j,i,2))**2 + 4.d0*se2sr(j,i,3)*se2sr(j,i,3))) / dtavg
            shrheat(j,i) = sshrheat(j,i) / dtavg
            se2sr(j,i,:) = 0.d0
            sshrheat(j,i) = 0.d0
        end do
    end do
    !$OMP end parallel do
    dtavg = 0.d0
    nsrate = 0
elseif( nsrate .eq. -1 ) then
    !$OMP parallel do
    !$ACC parallel loop collapse(2) async(1)
    do i = 1,nx-1
        do j = 1, nz-1
            e2sr(j,i) = (0.5d0 * sqrt((se2sr(j,i,1)-se2sr(j,i,2))**2 + 4.d0*se2sr(j,i,3)*se2sr(j,i,3))) / dtavg
            shrheat(j,i) = sshrheat(j,i) / dtavg
        end do
    end do
    !$OMP end parallel do
endif

nsrate = nsrate + 1

return
end
