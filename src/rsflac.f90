
! Re-starting FLAC

subroutine rsflac
use arrays
USE marker_data
use params
implicit none

integer, parameter :: kindr=8, kindi=4

real(kindr), allocatable :: dum1(:),dum2(:,:)
integer(kindi), allocatable :: dum11(:), idum2(:,:)
real*8 rtime, rdt, time_my
character*200 msg
integer :: i, iph, j, k, n, nrec, nwords

! TODO: include tracer information for restart
if (iint_tracer.eq.1) then
    stop 'Must disable tracers in restart'
endif

open( 1, file='_contents.rs', status='old' )
read( 1, * ) nrec, nloop, time_my, nmarkers, nmtracers
close(1)


! Read time and dt
open (1,file='time.rs',access='direct',recl=2*8) 
read (1,rec=nrec) rtime, rdt
close (1)
time = rtime
dt = rdt
!$ACC update device(nrec,nloop,nmarkers,nmtracers)
!$ACC update device(time,dt)

dvol = 0

! Coordinates and velocities
nwords = nz*nx*2

open (1,file='cord.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) cord
close (1)

open (1,file='dhacc.rs',access='direct',recl=(nx-1)*kindr)
read (1,rec=nrec) dhacc(1:nx-1)
close (1)

open (1,file='extr_acc.rs',access='direct',recl=(nx-1)*kindr)
read (1,rec=nrec) extr_acc(1:nx-1)
close (1)

open (1,file='vel.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) vel
close (1)

! Strain
nwords = 3*(nz-1)*(nx-1)

open (1,file='strain.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) strain
close (1)


! Stress
nwords = 4*4*(nx-1)*(nz-1)

open (1,file='stress.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) stress0
close (1)


! 2-D (nx*nz) arrays - nodes defined
nwords = nz*nx

! Temperature
open (1,file='temp.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) temp
close (1)


! 2-D (nx-1)*(nz-1) arrays - elements defined
allocate( dum2(nz-1,nx-1) )

nwords = (nz-1)*(nx-1)

! Phases
allocate( idum2(nz-1,nx-1) )
open (1,file='phase.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) idum2
close (1)
iphase(1:nz-1,1:nx-1) = idum2(1:nz-1,1:nx-1)
deallocate( idum2 )

! Check if viscous rheology present
ivis_present = 0
do i = 1,nx-1
    do j = 1, nz-1
        iph = iphase(j,i)
        if( irheol(iph).eq.3 .or. irheol(iph).ge.11 ) ivis_present = 1
    end do
end do
!$ACC update device(ivis_present)

! Plastic strain
open (1,file='aps.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum2
close (1)
aps(1:nz-1,1:nx-1) = dum2(1:nz-1,1:nx-1)

! Heat sources
open (1,file='source.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) dum2
close (1)
source(1:nz-1,1:nx-1) = dum2(1:nz-1,1:nx-1)

deallocate( dum2 )

!$ACC update device(cord,dhacc,extr_acc,vel,strain,stress0,temp, &
!$ACC               iphase,aps,source)

if (iint_marker.eq.1) then


! Markers
nwords= nmarkers
allocate (dum1(nmarkers))
! Markers
open (1,file='xmarker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark_x(i) = dum1(i)
enddo


open (1,file='ymarker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark_y(i) = dum1(i)
enddo


open (1,file='xa1marker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark_a1(i) = dum1(i)
enddo


open (1,file='xa2marker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark_a2(i) = dum1(i)
enddo


open (1,file='xagemarker.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dum1
close (1)
do i = 1,nmarkers
mark_age(i) = dum1(i)
enddo


allocate(dum11(nmarkers))

open (1,file='xIDmarker.rs',access='direct',recl=nwords*kindi)
read (1,rec=nrec) dum11
close (1)
do i = 1,nmarkers
mark_ID(i) = dum11(i)
enddo


open (1,file='xntriagmarker.rs',access='direct',recl=nwords*kindi)
read (1,rec=nrec) dum11
close (1)
do i = 1,nmarkers
mark_ntriag(i) = dum11(i)
enddo

open (1,file='xphasemarker.rs',access='direct',recl=nwords*kindi)
read (1,rec=nrec) dum11
close (1)
do i = 1,nmarkers
mark_phase(i) = dum11(i)
enddo

open (1,file='xdeadmarker.rs',access='direct',recl=nwords*kindi)
read (1,rec=nrec) dum11
close (1)
do i = 1,nmarkers
mark_dead(i) = dum11(i)
enddo

deallocate(dum11)

! recount marker phase
mark_id_elem(:,:,:) = 0
nmark_elem(:,:) = 0
print *, nmarkers
do n = 1, nmarkers
    if(mark_dead(n) .eq. 0) cycle

     if(mark_ntriag(n).lt.1 .or. mark_ntriag(n).gt.2*(nx-1)*(nz-1)) then
         print *, 'Wrong marker ntriag', mark_ID(n), mark_ntriag(n)
         stop 999
     endif

    ! from ntriag, get element number
    k = mod(mark_ntriag(n) - 1, 2) + 1
    j = mod((mark_ntriag(n) - k) / 2, nz-1) + 1
    i = (mark_ntriag(n) - k) / 2 / (nz - 1) + 1

    !if(mark_ntriag(n) .ne. 2 * ( (nz-1)*(i-1)+j-1) + k) write(*,*), mark_ntriag(n), i,j,k

    if(nmark_elem(j,i) == max_markers_per_elem) then
        write(msg,*) 'Too many markers at element:', i, j, nmark_elem(j,i)
        call SysMsg(msg)
        cycle
    endif

    ! recording the id of markers belonging to the element
    nmark_elem(j, i) = nmark_elem(j, i) + 1
    mark_id_elem(nmark_elem(j, i), j, i) = n
enddo

!$ACC update device(mark_a1, mark_a2, mark_x, mark_y, mark_age, &
!$ACC               mark_dead, mark_ntriag, mark_phase, mark_ID, &
!$ACC               nmark_elem, mark_id_elem)

call marker2elem

endif

! Pressure at the bottom: pisos 
if( nyhydro .eq. 2 ) then
    open(1,file='pisos.rs')
    read(1,*) pisos
    close (1)
endif

! Calculate AREAS (Important: iphase is needed to calculate area!)
call init_areas
!$ACC update device(area)

! Distribution of REAL masses to nodes
call rmasses

! Boundary conditions
call init_bc

! Inertial masses and time steps (elastic, maxwell and max_thermal)
call dt_mass

if( ivis_present.eq.1 ) call init_visc

! Initiate parameters for stress averaging
dtavg=0
nsrate=-1

!Initialization
!$ACC update device(cord, temp, vel, stress0, force, balance, amass, rmass, &
!$ACC               area, dvol, strain, bc, ncod, junk2, xmpt, tkappa, &
!$ACC               iphase, mark_id_elem, nmark_elem, &
!$ACC               nopbou, ncodbou, idtracer, phase_ratio, dtopo, dhacc, extrusion, &
!$ACC               andesitic_melt_vol, extr_acc, strainr, aps, visn, e2sr, &
!$ACC               temp0, source, shrheat, bcstress, &
!$ACC               pt, barcord, cold, cnew, numtr, &
!$ACC               se2sr, sshrheat, dtavg, nsrate)

return
end
