
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
                    de11,de22,de12,de33,dv,dv1p,diss, &
                    quad_area, s0, s0a,s0b,sII_plas,sII_visc
double precision :: Eff_visc
integer :: i, j, k, iph, irh, &
           ipls


!$OMP Parallel Private(i,j,k,iph,irh,bulkm,rmu,coh,phi,psi, &
!$OMP                  stherm,hardn,vis, &
!$OMP                  de11,de22,de12,de33,dv,dv1p, &
!$OMP                  s11p,s22p,s12p,s33p, &
!$OMP                  s11v,s22v,s12v,s33v, &
!$OMP                  depl,ipls,diss, &
!$OMP                  sII_plas,sII_visc, &
!$OMP                  quad_area,s0a,s0b,s0)
!$OMP do schedule(guided)
!$ACC parallel loop gang vector collapse(2) private(depl,s11p,s22p,s12p,s33p,s11v,s22v,s12v,s33v) async(1)
do i = 1,nx-1
    do j = 1,nz-1
        ! iphase (j,i) is number of a phase NOT a rheology
        iph = iphase(j,i)
        irh = irheol(iph)

        ! Elastic modules & viscosity & plastic properties
        bulkm = rl(iph) + 2.d0*rm(iph)/3.d0
        rmu   = rm(iph)

        ! Thermal stresses (alfa_v = 3.e-5 1/K)
        stherm = 0.d0
        if (istress_therm.gt.0) stherm = -alfa(iph)*bulkm*(temp(j,i)-temp0(j,i))


        ! Preparation of plastic properties
        if (irh.eq.6 .or. irh.ge.11) call pre_plast(i,j,coh,phi,psi,hardn)
              
        ! Re-evaluate viscosity
        if (irh.eq.3 .or. irh.eq.12) then 
            if( mod(nloop,ifreq_visc).eq.0) visn(j,i) = Eff_visc(j,i)
        endif
        vis = visn(j,i)

        dv1p = 0.d0
        if (itype_melting.eq.2) dv1p = -new_intrusion(j,i)

        ! Cycle by triangles
        do k = 1,4

            ! Incremental strains
            de11 = strainr(1,k,j,i)*dt
            de22 = strainr(2,k,j,i)*dt
            de12 = strainr(3,k,j,i)*dt
            de33 = 0.d0
            dv = dvol(j,i,k)

            if (itype_melting.eq.2) dv = dv + dv1p

            s11p(k) = stress0(j,i,1,k) + stherm 
            s22p(k) = stress0(j,i,2,k) + stherm 
            s12p(k) = stress0(j,i,3,k) 
            s33p(k) = stress0(j,i,4,k) + stherm
            s11v(k) = s11p(k)
            s22v(k) = s22p(k)
            s12v(k) = s12p(k)
            s33v(k) = s33p(k)

            if (irh.eq.1) then
                ! elastic
                call elastic(bulkm,rmu,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de12)
                stress0(j,i,1,k) = s11p(k)
                stress0(j,i,2,k) = s22p(k)
                stress0(j,i,3,k) = s12p(k)
                stress0(j,i,4,k) = s33p(k)

            elseif (irh.eq.3) then
                ! viscous
                call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),de11,de22,de33,de12,dv,&
                     ndim,dt)
                stress0(j,i,1,k) = s11v(k)
                stress0(j,i,2,k) = s22v(k)
                stress0(j,i,3,k) = s12v(k)
                stress0(j,i,4,k) = s33v(k)

            elseif (irh.eq.6) then
                ! plastic
                call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,&
                     dv1p,ten_off,ndim)
                stress0(j,i,1,k) = s11p(k)
                stress0(j,i,2,k) = s22p(k)
                stress0(j,i,3,k) = s12p(k)
                stress0(j,i,4,k) = s33p(k)

            elseif (irh.ge.11) then 
                ! Mixed rheology (Maxwell or plastic)
                call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                    s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,dv1p,&
                    ten_off,ndim)
                call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),&
                    de11,de22,de33,de12,dv,&
                    ndim,dt)
            endif
        enddo

        if( irh.ge.11 ) then
            ! deside - elasto-plastic or viscous deformation
            sII_plas = (s11p(1)+s11p(2)+s11p(3)+s11p(4)-s22p(1)-s22p(2)-s22p(3)-s22p(4))**2 &
                     + 4*(s12p(1)+s12p(2)+s12p(3)+s12p(4))**2

            sII_visc = (s11v(1)+s11v(2)+s11v(3)+s11v(4)-s22v(1)-s22v(2)-s22v(3)-s22v(4))**2 &
                     + 4*(s12v(1)+s12v(2)+s12v(3)+s12v(4))**2

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
        
            ! For A and B couple:
            ! area(n,it) is INVERSE of "real" DOUBLE area (=1./det)
            quad_area = 1.d0/area(j,i,1) + 1.d0/area(j,i,2)
            s0a=0.5d0*(stress0(j,i,1,1)+stress0(j,i,2,1))
            s0b=0.5d0*(stress0(j,i,1,2)+stress0(j,i,2,2))
            s0=(s0a/area(j,i,2)+s0b/area(j,i,1))/quad_area
            stress0(j,i,1,1) = stress0(j,i,1,1) - s0a + s0
            stress0(j,i,2,1) = stress0(j,i,2,1) - s0a + s0
            stress0(j,i,1,2) = stress0(j,i,1,2) - s0b + s0
            stress0(j,i,2,2) = stress0(j,i,2,2) - s0b + s0

            ! For C and D couple:
            quad_area = 1.d0/area(j,i,3) + 1.d0/area(j,i,4)
            s0a=0.5d0*(stress0(j,i,1,3)+stress0(j,i,2,3))
            s0b=0.5d0*(stress0(j,i,1,4)+stress0(j,i,2,4))
            s0=(s0a/area(j,i,4)+s0b/area(j,i,3))/quad_area
            stress0(j,i,1,3) = stress0(j,i,1,3) - s0a + s0
            stress0(j,i,2,3) = stress0(j,i,2,3) - s0a + s0
            stress0(j,i,1,4) = stress0(j,i,1,4) - s0b + s0
            stress0(j,i,2,4) = stress0(j,i,2,4) - s0b + s0
        endif

        if (irh.eq.6 .or. irh.ge.11) then
            !  ACCUMULATED PLASTIC STRAIN
            ! Average the strain for pair of the triangles
            ! Note that area (n,it) is inverse of double area !!!!!
            aps(j,i) = aps(j,i) &
                + 0.5d0*( depl(1)/area(j,i,2)+depl(2)/area(j,i,1) ) / (1.d0/area(j,i,1)+1.d0/area(j,i,2)) &
                + 0.5d0*( depl(3)/area(j,i,4)+depl(4)/area(j,i,3) ) / (1.d0/area(j,i,3)+1.d0/area(j,i,4))
            if( aps(j,i) .lt. 0.d0 ) aps(j,i) = 0.d0
            if (av_intrusion(j,i).ne.0.) aps(j,i) = 0.

            !	write(*,*) depl(1),depl(2),depl(3),depl(4),area(j,i,1),area(j,i,2),area(j,i,3),area(j,i,4)

            ! LINEAR HEALING OF THE PLASTIC STRAIN
            if (tau_heal .ne. 0.d0) &
                 aps (j,i) = aps (j,i)/(1.d0+dt/tau_heal)
        end if

        ! TOTAL FINITE STRAIN
        strain(j,i,1) = strain(j,i,1) + 0.25d0*dt*(strainr(1,1,j,i)+strainr(1,2,j,i)+strainr(1,3,j,i)+strainr(1,4,j,i))
        strain(j,i,2) = strain(j,i,2) + 0.25d0*dt*(strainr(2,1,j,i)+strainr(2,2,j,i)+strainr(2,3,j,i)+strainr(2,4,j,i))
        strain(j,i,3) = strain(j,i,3) + 0.25d0*dt*(strainr(3,1,j,i)+strainr(3,2,j,i)+strainr(3,3,j,i)+strainr(3,4,j,i))

    enddo
enddo
!$OMP end do
!$OMP end parallel

return
end
