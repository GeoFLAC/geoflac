subroutine lpeuler2bar
!$ACC routine(euler2bar) seq
USE marker_data
use arrays
use params
implicit none
character*200 msg
integer :: n, i, j, k, kk, ntr, inc
double precision :: xx, yy, bar1, bar2

!$ACC kernels async(1)
mark_id_elem(:,:,:) = 0
nmark_elem(:,:) = 0
!$ACC end kernels

!$OMP parallel do private(n,k,j,i,xx,yy,bar1,bar2,ntr,inc)
!$ACC parallel loop async(1)
do n = 1 , nmarkers
    if (mark_dead(n).eq.0) cycle

    ! from ntriag, get element number
    k = mod(mark_ntriag(n) - 1, 2) + 1
    j = mod((mark_ntriag(n) - k) / 2, nz-1) + 1
    i = (mark_ntriag(n) - k) / 2 / (nz - 1) + 1

    ! bar2euler() (called just before this, with no output args) wrote the
    ! marker's Euler coordinates into mark_a1/mark_a2 as scratch space,
    ! since they're about to be overwritten below with the real new
    ! barycentric coordinates anyway.
    xx = mark_a1(n)
    yy = mark_a2(n)
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
enddo
!$OMP end parallel do

return

end subroutine lpeuler2bar
