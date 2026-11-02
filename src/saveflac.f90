!---------------------------------------------------------------
!      Saving state
!---------------------------------------------------------------

subroutine saveflac
use arrays
use params
USE marker_data
implicit none

integer, parameter :: kindr=8, kindi=4
integer nrec, nwords, u
real*8 rtime, rdt
integer(kindi), allocatable :: ibuf(:)   ! widening buffer for byte-sized marker fields

! define record number and write it to contents

open (1,file='_contents.save')
nrec = 1
write( 1, '(i4,1x,i8,1x,f6.2,1x,i9,1x,i9)' ) nrec, nloop, time/sec_year/1.d6, &
     nmarkers,0
close(1)

! Time and dt
open (1,file='time.rs',access='direct',recl=2*8) 
rtime = time
rdt = dt
write (1,rec=nrec) rtime, rdt
close (1) 


! Coordinates and velocities
nwords = nz*nx*2

open (1,file='cord.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) cord
close (1)

open (1,file='dhacc.rs',access='direct',recl=(nx-1)*kindr)
write (1,rec=nrec) dhacc(1:nx-1)
close (1)

open (1,file='dhacc_corr.rs',access='direct',recl=nx*kindr)
write (1,rec=nrec) dhacc_correction(1:nx)
close (1)

open (1,file='extr_acc.rs',access='direct',recl=(nx-1)*kindr)
write (1,rec=nrec) extr_acc(1:nx-1)
close (1)

open (1,file='q_init.rs',access='direct',recl=nx*kindr)
write (1,rec=nrec) q_init
close (1)

open (1,file='vel_flex_old.rs',access='direct',recl=nx*kindr)
write (1,rec=nrec) vel_flex_old
close (1)

open (1,file='vel.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) vel
close (1)


! Strain
nwords = 3*(nz-1)*(nx-1)
open (1,file='strain.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) strain
close (1)


! Stress
nwords = 4*4*(nx-1)*(nz-1)
open (1,file='stress.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) stress0
close (1)


! 2-D (nx*nz) arrays - nodes defined
nwords = nz*nx

! Temperature
open (1,file='temp.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) temp
close (1)

! Original location
open (1,file='xoriginal.rs',access='direct',recl=nwords*kindr)
write (1,rec=nrec) xoriginal
close (1)
open (1,file='zoriginal.rs',access='direct',recl=nwords*kindr)
write (1,rec=nrec) zoriginal
close (1)

! 2-D (nx-1)*(nz-1) arrays - elements defined
nwords = (nz-1)*(nx-1)

! Phases
open (1,file='phase.rs',access='direct',recl=nwords*kindr)
write (1,rec=nrec) iphase
close (1)

! Plastic strain
open (1,file='aps.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) aps
close (1)

! Magma
open (1,file='fmagma.rs',access='direct',recl=nwords*kindr)
write (1,rec=nrec) fmagma
close (1)

! Magma2
open (1,file='fmagma2.rs',access='direct',recl=nwords*kindr)
write (1,rec=nrec) fmagma2
close (1)

! Heat sources
open (1,file='source.rs',access='direct',recl=nwords*kindr)
write (1,rec=nrec) source
close (1)

! Dike-intrusion eigenstrain/marker tracking (persistent parts only;
! dike_added/dike_released are same-timestep scratch, not saved)
open (1,file='dike_backlog.rs',access='direct',recl=nwords*kindr)
write (1,rec=nrec) dike_backlog
close (1)
open (1,file='dike_marker_vol.rs',access='direct',recl=nwords*kindr)
write (1,rec=nrec) dike_marker_vol
close (1)
open (1,file='mor_marker_vol.rs',access='direct',recl=nwords*kindr)
write (1,rec=nrec) mor_marker_vol
close (1)

! Minimum element area of the initial grid (see setflac.f90); constant for
! the whole run, so a plain formatted file (not indexed by nrec) suffices.
open(newunit=u, file='Amin.rs')
write(u,*) Amin
close(u)

! Markers
! Euler (x,y) coordinates are not saved: they're not part of the
! persistent marker state and are regenerated via bar2euler() from
! mark_a1/mark_a2/mark_ntriag + cord whenever needed after a restart.
nwords = nmarkers
nrec = 1
open (1,file='marker1.rs',access='direct',recl=nwords*kindr)
write (1,rec=nrec) mark_a1(1:nmarkers)
nrec = nrec + 1
write (1,rec=nrec) mark_a2(1:nmarkers)
nrec = nrec + 1
write (1,rec=nrec) mark_age(1:nmarkers)
nrec = nrec + 1
close (1)

! marker2.rs stores each field as a full-width (kindi=4) record regardless
! of the in-memory type, so mark_dead/mark_phase (1 byte each) are widened
! into ibuf before writing; mark_ntriag is already full width.
allocate(ibuf(nmarkers))
nrec = 1
open (1,file='marker2.rs',access='direct',recl=nwords*kindi)
ibuf = mark_dead(1:nmarkers)
write (1,rec=nrec) ibuf
nrec = nrec + 1
write (1,rec=nrec) mark_ntriag(1:nmarkers)
nrec = nrec + 1
ibuf = mark_phase(1:nmarkers)
write (1,rec=nrec) ibuf
nrec = nrec + 1
close (1)
deallocate(ibuf)
return
end
