subroutine outmarker
USE marker_data
use arrays
use params
implicit none
integer, parameter :: kindr=4, kindi=4
integer :: i, l, nrec, nwords
real(kindr) :: D1d(nmarkers)
integer(kindi) :: D1i(nmarkers)
double precision :: xtmp(nmarkers), ytmp(nmarkers)

character*100 fn

! mark_a1/mark_a2 are still needed below (as the true barycentric
! coordinates), so get Euler coordinates into a local buffer instead of
! letting bar2euler overwrite them in place.
call bar2euler_xy(xtmp, ytmp)
!$ACC wait

nrec = 0
! define record number and write it to contents
if( lastout .eq. 1 ) then
    nrec = 1
    open (1,file='_markers.0')
else
    open (1,file='_markers.0',status='old',err=5)

    do while (.TRUE.)
        read( 1, *, end=10 ) nrec
    end do
    5 continue
    open (1,file='_markers.0',position='append')
    nrec = 0
    10 continue
    nrec = nrec + 1
    backspace(1) ! Neede for further writing since EOF has been reached.
endif
write( 1, '(i6,1x,i8,1x,i8,1x,f10.6)' ) nrec, nloop,nmarkers, time/sec_year/1.d6
close(1)

!! Since the number of markers changes with time, the marker data cannot be
!! output as a single unformatted file. Output different files for each record.

! Coordinates  [km]
nwords = nmarkers
write(fn,'(A,I6.6,A)') 'marker1.', nrec, '.0'
open (1,file=fn,access='direct',recl=nwords*kindr)

do i = 1, nmarkers
    D1d(i)= real(xtmp(i) * 1d-3)
enddo
write (1,rec=1) D1d

do i = 1,nmarkers
    D1d(i)= real(ytmp(i) * 1d-3)
enddo
write (1,rec=2) D1d

! Age [Myrs]
do i = 1,nmarkers
    D1d(i)= real(mark_age(i) / sec_year / 1.d6)
enddo
write (1,rec=3) D1d

! Barycentric coordinates
do i = 1,nmarkers
    D1d(i)= real(mark_a1(i))
enddo
write (1,rec=4) D1d

do i = 1,nmarkers
    D1d(i)= real(mark_a2(i))
enddo
write (1,rec=5) D1d

close (1)


write(fn,'(A,I6.6,A)') 'marker2.', nrec, '.0'
open (1,file=fn,access='direct',recl=nwords*kindi)

do l = 1,nmarkers
    D1i(l)= mark_dead(l)
enddo
write (1,rec=1) D1i

do l = 1,nmarkers
    D1i(l)= mark_phase(l)
enddo
write (1,rec=2) D1i

do l = 1,nmarkers
    D1i(l)= mark_ntriag(l)
enddo
write (1,rec=3) D1i

close (1)

return
end subroutine outmarker
