!---------------------------------------------------------------
!      Saving state
!---------------------------------------------------------------

subroutine saveflac
use arrays
use params
USE marker_data
implicit none

integer, parameter :: kindr=8, kindi=4
integer nrec, nwords
real*8 rtime, rdt

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

open (1,file='extr_acc.rs',access='direct',recl=(nx-1)*kindr)
write (1,rec=nrec) extr_acc(1:nx-1)
close (1)

open (1,file='xlab.rs',access='direct',recl=nx*kindr)
write (1,rec=nrec) xlab(1:nx)
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

!melt
open (1,file='melt.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) Eff_melt
close (1)

!stored_intrusion
open (1,file='stored_intru.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) stored_intrusion
close (1)

!dv_intr
open (1,file='dv_intr.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) dv_intr
close (1)

! Heat sources
open (1,file='source.rs',access='direct',recl=nwords*kindr) 
write (1,rec=nrec) source
close (1)

! Markers
call bar2euler
!$ACC wait

nwords = nmarkers
nrec = 1
open (1,file='marker1.rs',access='direct',recl=nwords*kindr)
write (1,rec=nrec) mark_a1(1:nmarkers)
nrec = nrec + 1
write (1,rec=nrec) mark_a2(1:nmarkers)
nrec = nrec + 1
write (1,rec=nrec) mark_x(1:nmarkers)
nrec = nrec + 1
write (1,rec=nrec) mark_y(1:nmarkers)
nrec = nrec + 1
write (1,rec=nrec) mark_age(1:nmarkers)
nrec = nrec + 1
write (1,rec=nrec) mark_zpressn(1:nmarkers)
nrec = nrec + 1
write (1,rec=nrec) mark_Emeltm(1:nmarkers)
nrec = nrec + 1
close (1)

nrec = 1
open (1,file='marker2.rs',access='direct',recl=nwords*kindi)
write (1,rec=nrec) mark_dead(1:nmarkers)
nrec = nrec + 1
write (1,rec=nrec) mark_ntriag(1:nmarkers)
nrec = nrec + 1
write (1,rec=nrec) mark_phase(1:nmarkers)
nrec = nrec + 1
write (1,rec=nrec) mark_ID(1:nmarkers)
nrec = nrec + 1
close (1)
return
end
