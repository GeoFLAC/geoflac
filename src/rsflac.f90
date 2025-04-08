
! Re-starting FLAC

subroutine rsflac
use arrays
use params
USE marker_data
implicit none

integer, parameter :: kindr=8, kindi=4
integer :: nrec, nwords, i, j, k, iph, n
real*8 rtime, rdt, time_my
character*200 msg

open( 1, file='_contents.rs', status='old' )
read( 1, * ) nrec, nloop, time_my, nmarkers, i
close(1)


! Read time and dt
open (1,file='time.rs',access='direct',recl=2*8) 
read (1,rec=nrec) rtime, rdt
close (1)
time = rtime
dt = rdt

dvol = 0

! Coordinates and velocities
nwords = nz*nx*2

open (1,file='cord.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) cord
close (1)

! min element width and thickness
dxmin = minval(cord(2:nx,1,1) - cord(1:nx-1,1,1))
dzmin = minval(cord(1,1:nz-1,2) - cord(1,2:nz,2))

open (1,file='dhacc.rs',access='direct',recl=(nx-1)*kindr)
read (1,rec=nrec) dhacc(1:nx-1)
close (1)

open (1,file='extr_acc.rs',access='direct',recl=(nx-1)*kindr)
read (1,rec=nrec) extr_acc(1:nx-1)
close (1)

open (1,file='xlab.rs',access='direct',recl=nx*kindr)
read (1,rec=nrec) xlab(1:nx)
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
nwords = (nz-1)*(nx-1)

! Phases
open (1,file='phase.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) iphase
close (1)


! Check if viscous rheology present
call check_visc_rheol

! Plastic strain
open (1,file='aps.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) aps
close (1)

!melt
open (1,file='melt.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) Eff_melt
close (1)

!stored_intrusion
open (1,file='stored_intru.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) stored_intrusion
close (1)

!dv_intr
open (1,file='dv_intr.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dv_intr
close (1)

! Magma
open (1,file='fmagma.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) fmagma
close (1)

! Heat sources
open (1,file='source.rs',access='direct',recl=nwords*kindr) 
read (1,rec=nrec) source
close (1)

! Markers
nwords = nmarkers
nrec = 1
open (1,file='marker1.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) mark_a1(1:nmarkers)
nrec = nrec + 1
read (1,rec=nrec) mark_a2(1:nmarkers)
nrec = nrec + 1
read (1,rec=nrec) mark_x(1:nmarkers)
nrec = nrec + 1
read (1,rec=nrec) mark_y(1:nmarkers)
nrec = nrec + 1
read (1,rec=nrec) mark_age(1:nmarkers)
nrec = nrec + 1
read (1,rec=nrec) mark_zpressn(1:nmarkers)
nrec = nrec + 1
read (1,rec=nrec) mark_Emeltm(1:nmarkers)
nrec = nrec + 1
close (1)

nrec = 1
open (1,file='marker2.rs',access='direct',recl=nwords*kindi)
read (1,rec=nrec) mark_dead(1:nmarkers)
nrec = nrec + 1
read (1,rec=nrec) mark_ntriag(1:nmarkers)
nrec = nrec + 1
read (1,rec=nrec) mark_phase(1:nmarkers)
nrec = nrec + 1
read (1,rec=nrec) mark_ID(1:nmarkers)
nrec = nrec + 1
close (1)

! recount marker phase
mark_id_elem(:,:,:) = 0
nmark_elem(:,:) = 0
print *, '# of markers:', nmarkers
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

call marker2elem

! Pressure at the bottom: pisos 
if( nyhydro .eq. 2 ) then
    open(1,file='pisos.rs')
    read(1,*) pisos
    close (1)
endif
!$ACC update device(pisos) async(1)

! Calculate AREAS (Important: iphase is needed to calculate area!)
call init_areas

! Boundary conditions
call init_bc

temp0 = temp
shrheat = 0
sshrheat = 0
dtopo = 0
extrusion = 0
fmelt = 0
se2sr = 1d-16
e2sr = 1d-16

! MOR magma migration related
zpressm = 0.
new_intrusion = 0.
av_intrusion = 0.
icrust_melt = 1
watercont = 0.000000 !wt fraction
xmodalcpx = 0.03 !wt fraction
xDF(1) = 0.0 !Initiatilize DF
xDT(1) = 0.0 !Initialize DT
dumC(1) = 0.0 !Runge Kutta
dumC(2) = 0.5 !Runge Kutta
dumC(3) = 0.5 !Runge Kutta
dumC(4) = 1.0 !Runge Kutta
istart_profile = 50

call update_acc

! Distribution of REAL masses to nodes
call rmasses

if( ivis_present.eq.1 ) call init_visc

! Inertial masses and time steps (elastic, maxwell and max_thermal)
call dt_mass

! Initiate parameters for stress averaging
dtavg=0
nsrate=-1

return
end
