!---------------------------------------------------------------
!      Output
!---------------------------------------------------------------

subroutine outflac
use arrays
use params
implicit none
integer, parameter :: kindr=4, kindi=4
integer :: nrec, nwords, i, j, ii, iph
double precision :: sxx, szz, sxz, syy, dtmpr, dl, strainI, strainII, srateII, stressI, stressII, Eff_dens, Eff_conduct
real(kindr), allocatable :: D1d(:), De(:,:), Dn2(:,:,:)
real(kindr) :: rtime

! define record number and write it to contents
if( lastout .eq. 1 ) then
    nrec = 1
    open (1,file='_contents.0')
else
    open (1,file='_contents.0',status='old',err=5)
    do while (.TRUE.)
        read( 1, *, end=10 ) nrec
    end do
5   continue
    open (1,file='_contents.0',position='append')
    nrec = 0
10  continue
    nrec = nrec + 1
    backspace(1)
endif
write( 1, '(i6,1x,i8,1x,f10.6)' ) nrec, nloop, time/sec_year/1.d6
close(1)

! Time
open (1,file='time.0',access='direct',recl=kindr)
rtime = real(time)
write (1,rec=nrec) rtime
close (1) 


! Coordinates in [km]
allocate( Dn2(nz,nx,2) )

nwords = nz*nx*2
Dn2(1:nz,1:nx,1:2) = real(cord(1:nz,1:nx,1:2) / 1000)
open (1,file='mesh.0',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) Dn2
close (1)

! Velocities in [cm/year]
if( io_vel.eq.1 ) then
    Dn2(1:nz,1:nx,1:2) = real(vel(1:nz,1:nx,1:2) * sec_year * 100)
    open (1,file='vel.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) Dn2
    close (1)
endif

! Temperature in [Celsius]
if( io_temp.eq.1 ) then
    nwords = nz*nx
    Dn2(1:nz,1:nx,1) = real(temp(1:nz,1:nx))
    open (1,file='temperature.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) Dn2(1:nz,1:nx,1)
    close (1)
endif

! Cooling rate (dT/dt) in [C/Myr]
nwords = nz*nx
Dn2(1:nz,1:nx,1) = real((temp(1:nz,1:nx) - temp0(1:nz,1:nx)) / (dt / (sec_year * 1.d6)))
open (1,file='coolingrate.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) Dn2(1:nz,1:nx,1)
close (1)


nwords = nz*nx
Dn2(1:nz,1:nx,1) = real(xoriginal(1:nz,1:nx))
open (1,file='xoriginal.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) Dn2(1:nz,1:nx,1)
close (1)
Dn2(1:nz,1:nx,1) = real(zoriginal(1:nz,1:nx))
open (1,file='zoriginal.0',access='direct',recl=nwords*kindr)
write (1,rec=nrec) Dn2(1:nz,1:nx,1)
close (1)


deallocate( Dn2 )


! 2-D (nx-1)*(nz-1) arrays - elements defined
allocate( De(nz-1,nx-1) )

nwords = (nz-1)*(nx-1)

! Strain rate II
if( io_srII.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            if( e2sr(j,i).ne.0.d0 ) then
                De(j,i) = real(dlog10( e2sr(j,i) ))
            else
                De(j,i) = 0
            endif
        enddo
    enddo
    open (1,file='srII.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)

    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(se2sr(j,i,1)/dtavg)
        end do
    end do
    open (1,file='srxx.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)

    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(se2sr(j,i,2)/dtavg)
        end do
    end do
    open (1,file='srzz.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)

    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(se2sr(j,i,3)/dtavg)
        end do
    end do
    open (1,file='srxz.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)
endif


! Strain
if( io_eII.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(strainII(j,i))
        end do
    end do
    open (1,file='eII.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)

    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(strainI(j,i))
        end do
    end do
    open (1,file='eI.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)

    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(strain(j,i,1))
        end do
    end do
    open (1,file='exx.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)

    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(strain(j,i,2))
        end do
    end do
    open (1,file='ezz.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)

    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(strain(j,i,3))
        end do
    end do
    open (1,file='exz.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)

endif

! Density
if( io_density.eq.1 ) then
   do i = 1, nx-1
   do j = 1, nz-1
          De(j,i) = real(Eff_dens(j,i))
   enddo
   enddo
    open (1,file='density.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)

   do i = 1, nx-1
   do j = 1, nz-1
          De(j,i) = real(0.5d0/area(j,i,1) + 0.5d0/area(j,i,2))
   enddo
   enddo
    open (1,file='area.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)
endif


! APS
if( io_aps.eq.1 ) then
    De(1:nz-1,1:nx-1) = real(aps(1:nz-1,1:nx-1))
    open (1,file='aps.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Stress II in [Pa]
if( io_sII.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(stressII(j,i))
        end do
    end do
    open (1,file='sII.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Sxx in [Pa]
if( io_sxx.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            sxx = 0.25d0 * (stress0(j,i,1,1)+stress0(j,i,1,2)+stress0(j,i,1,3)+stress0(j,i,1,4) )
            De(j,i) = real( sxx-stressI(j,i) )
        end do
    end do
    open (1,file='sxx.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Szz in [Pa]
if( io_szz.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            szz = 0.25d0 * (stress0(j,i,2,1)+stress0(j,i,2,2)+stress0(j,i,2,3)+stress0(j,i,2,4) )
            De(j,i) = real( szz-stressI(j,i) )
        end do
    end do
    open (1,file='szz.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Sxz and Syy in [Pa]
if( io_sxz.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            sxz = 0.25d0 * (stress0(j,i,3,1)+stress0(j,i,3,2)+stress0(j,i,3,3)+stress0(j,i,3,4))
            De(j,i) = real(sxz)
        end do
    end do
    open (1,file='sxz.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)

    do i = 1, nx-1
        do j = 1, nz-1
            syy = 0.25d0 * (stress0(j,i,4,1)+stress0(j,i,4,2)+stress0(j,i,4,3)+stress0(j,i,4,4))
            De(j,i) = real( syy-stressI(j,i) )
        end do
    end do
    open (1,file='syy.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)
endif


! Pressure in [Pa]
if( io_pres.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(stressI(j,i))
        end do
    end do
    open (1,file='pres.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif


! Temperature
if( io_temp.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            De(j,i) = real(0.25d0*( temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1) ))
        end do
    end do
    open (1,file='temp.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)
endif



! Phase
if( io_phase.eq.1 ) then
    open (1,file='phase.0',access='direct',recl=nwords*kindi)
    write (1,rec=nrec) iphase(1:nz-1,1:nx-1)
    close (1)
endif


! Viscosities (log)
if( io_visc.eq.1 ) then
    De(1:nz-1,1:nx-1) = real(dlog10( visn(1:nz-1,1:nx-1) ))
    open (1,file='visc.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif



! Heat sources
if( io_src.eq.1 ) then
    De(1:nz-1,1:nx-1) = real(source(1:nz-1,1:nx-1))
    open (1,file='src.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif



! Melt fraction in the mantle wedge
if( io_src.eq.1 ) then
    De(1:nz-1,1:nx-1) = real(fmelt(1:nz-1,1:nx-1))
    open (1,file='fmelt.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)

    De(1:nz-1,1:nx-1) = real(fmelt2(1:nz-1,1:nx-1))
    open (1,file='fmelt2.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)
    
    De(1:nz-1,1:nx-1) = real(fmagma(1:nz-1,1:nx-1))
    open (1,file='fmagma.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)
    
    De(1:nz-1,1:nx-1) = real(fmagma2(1:nz-1,1:nx-1))
    open (1,file='fmagma2.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) De
    close (1)
endif



! Energy dissipation
if( io_diss.eq.1 ) then
    do i = 1, nx-1
        do j = 1, nz-1
            if(ishearh.ne.0) then
               iph = iphase(j,i)
               De(j,i) = real(shrheat(j,i)/den(iph)/hs)
            else
               De(j,i) = 0
            endif
        enddo
    enddo
    open (1,file='diss.0',access='direct',recl=nwords*kindr) 
    write (1,rec=nrec) De
    close (1)
endif

deallocate( De )


! 1-D nx array - nodes defined
allocate( D1d(nx) )
nwords = nx

! Surface heat flow
if( io_hfl.eq.1 ) then
    do i = 1,nx
        ii = min(i,nx-1)
        dtmpr = temp(2,i) - temp(1,i)
        dl = -(cord(2,i,2)-cord(1,i,2))/1000
        D1d(i) = real(Eff_conduct(1,ii) * dtmpr/dl)
    end do
    open (1,file='hfl.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) D1d
    close (1)
endif


! Topo
if( io_topo.eq.1 ) then
    do i = 1,nx
        D1d(i) = real(cord(1,i,2)/1000)
    end do
    open (1,file='topo.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) D1d
    close (1)

    do i = 1,nx
        D1d(i) = real(dtopo(i)/dt*1000*sec_year)  ! mm/yr
    end do
    open (1,file='dtopo.0',access='direct',recl=nwords*kindr)
    write (1,rec=nrec) D1d
    close (1)

    ! extrusion
    do i = 1,nx-1
        D1d(i) = real(extrusion(i))
    end do
    open (1,file='extrusion.0',access='direct',recl=(nwords-1)*kindr)
    write (1,rec=nrec) D1d(1:nx-1)
    close (1)
endif


deallocate( D1d )


return 
end
