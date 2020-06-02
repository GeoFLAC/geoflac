MODULE marker_data
  integer, parameter :: max_markers = 10000000

  !!! maximum number of ELEMENTS !!!
  integer, parameter :: max_markers_per_elem=32

  SAVE
  double precision, allocatable :: mark_a1(:), mark_a2(:) ! baricentric coordinates
  double precision, allocatable :: mark_x(:), mark_y(:)   ! Euler coordinates
  double precision, allocatable :: mark_age(:)            ! creation time
  integer, allocatable :: mark_dead(:)
  integer, allocatable :: mark_ntriag(:)      ! number of FE-triangle
  integer, allocatable :: mark_phase(:)
  integer, allocatable :: mark_ID(:)          ! unique ID-number

  integer, allocatable :: nphase_counter(:,:,:), ntopmarker(:), itopmarker(:,:)


  !$ACC declare create(mark_a1, mark_a2, mark_x, mark_y, mark_age, &
  !$ACC                mark_dead, mark_ntriag, mark_phase, mark_ID, &
  !$ACC                nphase_counter, ntopmarker, itopmarker)
  contains

  subroutine allocate_markers(nz, nx)
    implicit none

    integer, intent(in) :: nz, nx
    integer :: max_markers
    max_markers = nz * nx * max_markers_per_elem

    allocate(mark_a1(max_markers), mark_a2(max_markers), &
             mark_x(max_markers), mark_y(max_markers), &
             mark_age(max_markers), &
             mark_dead(max_markers), &
             mark_ntriag(max_markers), &
             mark_phase(max_markers), &
             mark_ID(max_markers))

    allocate(nphase_counter(20, nz-1, nx-1))
    allocate(ntopmarker(nx))
    allocate(itopmarker(max_markers_per_elem, nx-1))

    !$ACC update device(mark_a1, mark_a2, mark_x, mark_y, mark_age, &
    !$ACC               mark_dead, mark_ntriag, mark_phase, mark_ID, &
    !$ACC               nphase_counter, ntopmarker, itopmarker)

  end subroutine

  subroutine add_marker(x, y, iph, age, kk, j, i, inc)
    ! Add a marker at physical coordinate (x, y), with phase iph and age, to
    ! element (j, i). The current (before adding thsi marker) marker size
    ! is kk. If (x, y) is not within the element, inc is set to 0 and
    ! marker not added. Otherwise, marker is added to "mark" array and kk
    ! incremented by 1.
    !

    !$ACC routine seq
    use arrays
    use params
    use euler2bar
    implicit none
    integer, intent(in) :: iph, j, i
    integer, intent(out) :: inc
    integer, intent(inout) :: kk
    double precision, intent(in) :: x, y, age
    integer :: ntr
    double precision :: bar1, bar2
  
    !character*200 msg
  
    call check_inside(x , y, bar1, bar2, ntr, i, j, inc)
    if(inc.eq.0) return
  
    if(j == 1) then
        if(ntopmarker(i) == max_markers_per_elem) then
            !write(msg,*) 'Too many markers at surface elements:', i, ntopmarker(i)
            !call SysMsg(msg)
            !call SysMsg('Marker skipped, not added!')
            inc = 0
            return
        endif
        ! recording the id of markers belonging to surface elements
        !$ACC atomic update
        ntopmarker(i) = ntopmarker(i) + 1
        !$ACC atomic write
        itopmarker(ntopmarker(i), i) = kk + 1
    end if

    !$ACC atomic update
    kk = kk + 1

    if(kk > max_markers) then
      !call SysMsg('ADD_MARKER: # of markers exceeds max. value. Please increase mark array size.')
      stop 15
    endif

    mark_x(kk) = x
    mark_y(kk) = y
    mark_dead(kk) = 1
    mark_ID(kk) = kk
    mark_a1(kk) = bar1
    mark_a2(kk) = bar2
    mark_age(kk) = age
    mark_ntriag(kk) = ntr
    mark_phase(kk) = iph
  
    !$ACC atomic update
    nphase_counter(iph,j,i) = nphase_counter(iph,j,i) + 1
  
  end subroutine add_marker
  

  subroutine newphase2marker (j1, j2, i1, i2, iph)
    !$ACC routine seq
    use arrays
    use params
    implicit none
    
    integer :: j1, j2, i1, i2, iph, &
               kk, n, k, j, i
    
    ! reset the markers within elements in the rectangular region
    
    do kk = 1 , nmarkers
        if (mark_dead(kk).eq.0) cycle
        n = mark_ntriag(kk)
        k = mod(n - 1, 2) + 1
        j = mod((n - k) / 2, nz-1) + 1
        i = (n - k) / 2 / (nz - 1) + 1
    
        if(j>=j1 .and. j<=j2 .and. i>=i1 .and. i<=i2) then
            !$ACC atomic update
            nphase_counter(mark_phase(kk),j,i) = nphase_counter(mark_phase(kk),j,i) - 1
            !$ACC atomic write
            mark_phase(kk) = iph
            !$ACC atomic update
            nphase_counter(iph,j,i) = nphase_counter(iph,j,i) + 1
        endif
    enddo
    
    iphase(j1:j2,i1:i2) = iph
    phase_ratio(:,j1:j2,i1:i2) = 0.d0
    phase_ratio(iph,j1:j2,i1:i2) = 1.d0
    
    return
    end subroutine newphase2marker
    
END MODULE marker_data
