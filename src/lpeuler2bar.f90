subroutine lpeuler2bar
!$ACC routine(euler2bar) seq
USE marker_data
use arrays
use params

include 'precision.inc'

!$ACC kernels async(1)
!$OMP sections
!$OMP section
mark_id_elem(:,:,:) = 0
!$OMP section
nmark_elem(:,:) = 0
!$OMP section
Emeltcounter(:,:) = 0.
!$OMP section
zpresscounter(:,:) = 0.
!$OMP end sections
!$ACC end kernels

!$OMP parallel do private(n,k,j,i,xx,yy,bar1,bar2,ntr,inc)
!$ACC parallel loop async(1)
do n = 1 , nmarkers
    if (mark_dead(n).eq.0) cycle

    ! from ntriag, get element number
    k = mod(mark_ntriag(n) - 1, 2) + 1
    j = mod((mark_ntriag(n) - k) / 2, nz-1) + 1
    i = (mark_ntriag(n) - k) / 2 / (nz - 1) + 1

    xx = mark_x(n)
    yy = mark_y(n)
    call euler2bar(xx,yy,bar1,bar2,ntr,i,j,inc)
    if (inc.eq.0) then
        !write(*,*) i,j,ntr,xx,yy
        mark_dead(n) = 0
        cycle
    endif

    mark_a1(n) = bar1
    mark_a2(n) = bar2
    mark_ntriag(n) = ntr

    if(nmark_elem(j, i) == max_markers_per_elem) then
        !write(msg,*) 'Too many markers at elements in lpeuler2bar:', i, j, n
        !call SysMsg(msg)
        mark_dead(n) = 0
        cycle
    endif
    !$OMP atomic capture
    !$ACC atomic capture
    ! recording the id of markers belonging to surface elements
    nmark_elem(j, i) = nmark_elem(j, i) + 1
    kk = nmark_elem(j, i)
    !$ACC end atomic
    !$OMP end atomic
    mark_id_elem(kk, j, i) = n
    !$ACC atomic update
    !$OMP atomic update
    zpresscounter(j,i) = zpresscounter(j,i) + mark_zpressn(n)
    !$ACC atomic update
    !$OMP atomic update
    Emeltcounter(j,i) = Emeltcounter(j,i) + mark_Emeltm(n)
enddo
!$OMP end parallel do

return
end subroutine lpeuler2bar
