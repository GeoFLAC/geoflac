
! Re-starting FLAC

subroutine rsflac
use arrays
use params
USE marker_data
implicit none

integer, parameter :: kindr=8, kindi=4
integer :: nrec, nwords, i, j, k, iph, n, u, ios
real*8 rtime, rdt, time_my
character*200 msg
integer(kindi), allocatable :: ibuf(:)   ! widening buffer for byte-sized marker fields

open(newunit=u, file='_contents.rs', status='old' )
read(u, * ) nrec, nloop, time_my, nmarkers, i
close(u)


! Read time and dt
open(newunit=u, file='time.rs', access='direct', recl=2*8) 
read(u, rec=nrec) rtime, rdt
close(u)
time = rtime
dt = rdt

dvol = 0

! Coordinates and velocities
nwords = nz*nx*2

open(newunit=u, file='cord.rs', access='direct', recl=nwords*kindr) 
read(u, rec=nrec) cord
close(u)

! min element width and thickness
dxmin = minval(cord(1,2:nx,1) - cord(1,1:nx-1,1))
dzmin = minval(cord(1:nz-1,1,2) - cord(2:nz,1,2))

open(newunit=u, file='dhacc.rs', access='direct', recl=(nx-1)*kindr)
read(u, rec=nrec) dhacc(1:nx-1)
close(u)

open(newunit=u, file='dhacc_corr.rs', access='direct', recl=nx*kindr)
read(u, rec=nrec) dhacc_correction(1:nx)
close(u)
!$ACC update device(dhacc, dhacc_correction) async(1)

open(newunit=u, file='extr_acc.rs', access='direct', recl=(nx-1)*kindr)
read(u, rec=nrec) extr_acc(1:nx-1)
close(u)

open(newunit=u, file='q_init.rs', access='direct', recl=nx*kindr, status='old', iostat=ios)
if (ios == 0) then
    read(u, rec=nrec) q_init
    close(u)
else
    q_init = 0.d0
endif

open(newunit=u, file='vel_flex_old.rs', access='direct', recl=nx*kindr, status='old', iostat=ios)
if (ios == 0) then
    read(u, rec=nrec) vel_flex_old
    close(u)
else
    vel_flex_old = 0.d0
endif

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

! Original location
open (1,file='xoriginal.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) xoriginal
close (1)
open (1,file='zoriginal.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) xoriginal
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

! Magma
open (1,file='fmagma.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) fmagma
close (1)

! Magma2
open (1,file='fmagma2.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) fmagma2
close (1)

! Heat sources
open (1,file='source.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) source
close (1)

! Dike-intrusion eigenstrain/marker tracking (persistent parts only)
open (1,file='dike_backlog.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dike_backlog
close (1)
open (1,file='dike_marker_vol.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) dike_marker_vol
close (1)
open (1,file='mor_marker_vol.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) mor_marker_vol
close (1)

! Minimum element area of the initial grid (see setflac.f90)
open(newunit=u, file='Amin.rs')
read(u,*) Amin
close(u)
!$ACC update device(Amin) async(1)

! Markers
! Euler (x,y) coordinates are not part of the restart file; see the
! matching comment in saveflac.f90.
nwords = nmarkers
nrec = 1
open (1,file='marker1.rs',access='direct',recl=nwords*kindr)
read (1,rec=nrec) mark_a1(1:nmarkers)
nrec = nrec + 1
read (1,rec=nrec) mark_a2(1:nmarkers)
nrec = nrec + 1
read (1,rec=nrec) mark_age(1:nmarkers)
nrec = nrec + 1
close (1)

! marker2.rs stores each field as a full-width (kindi=4) record regardless
! of the in-memory type, so mark_dead/mark_phase (1 byte each) are read
! into the 4-byte ibuf and narrowed explicitly; mark_ntriag is already
! full width and can be read directly.
allocate(ibuf(nmarkers))
nrec = 1
open (1,file='marker2.rs',access='direct',recl=nwords*kindi)
read (1,rec=nrec) ibuf
mark_dead(1:nmarkers) = int(ibuf, 1)
nrec = nrec + 1
read (1,rec=nrec) mark_ntriag(1:nmarkers)
nrec = nrec + 1
read (1,rec=nrec) ibuf
mark_phase(1:nmarkers) = int(ibuf, 1)
nrec = nrec + 1
close (1)
deallocate(ibuf)

! recount marker phase
mark_id_elem(:,:,:) = 0
nmark_elem(:,:) = 0
print *, '# of markers:', nmarkers
do n = 1, nmarkers
    if(mark_dead(n) .eq. 0) cycle

     if(mark_ntriag(n).lt.1 .or. mark_ntriag(n).gt.2*(nx-1)*(nz-1)) then
         print *, 'Wrong marker ntriag', n, mark_ntriag(n)
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
!$ACC kernels async(1)
shrheat = 0
sshrheat = 0
dtopo = 0
extrusion = 0
fmelt = 0
fmelt2 = 0
se2sr = 1d-16
e2sr = 1d-16
!$ACC end kernels

call update_acc

! Distribution of REAL masses to nodes
call rmasses

if( ivis_present.eq.1 ) call init_visc

! Inertial masses and time steps (elastic, maxwell and max_thermal)
call dt_mass

! Initiate parameters for stress averaging
dtavg=0
nsrate=-1
istart_profile = 50
! dtavg/nsrate are reset after update_acc's bulk push above (which must
! run before rmasses/init_visc/dt_mass, since those need nx, nz, and the
! other params scalars/small arrays already on device); push these two
! explicitly so the device copy doesn't retain a stale pre-reset value.
!$ACC update device(dtavg, nsrate) async(1)

return
end
