MODULE marker_data
  integer :: max_markers

  !!! maximum number of ELEMENTS !!!
  integer, parameter :: max_markers_per_elem=32

  SAVE
  double precision, allocatable :: mark_a1(:), mark_a2(:) ! baricentric coordinates
  double precision, allocatable :: mark_x(:), mark_y(:)   ! Euler coordinates
  double precision, allocatable :: mark_age(:)            ! creation time
  double precision, allocatable :: mark_zpressn(:)        ! normal pressure
  double precision, allocatable :: mark_Emeltm(:)         ! melt
  integer, allocatable :: mark_dead(:)
  integer, allocatable :: mark_ntriag(:)      ! number of FE-triangle
  integer, allocatable :: mark_phase(:)
  integer, allocatable :: mark_ID(:)          ! unique ID-number

  integer, allocatable :: mark_id_elem(:,:,:), nmark_elem(:,:)
  !$ACC declare create(max_markers)

  contains

  subroutine allocate_markers(nz, nx)
    implicit none

    integer, intent(in) :: nz, nx

    max_markers = nz * nx * max_markers_per_elem
    !$ACC update device(max_markers) async(1)

    allocate(mark_a1(max_markers), mark_a2(max_markers), &
             mark_x(max_markers), mark_y(max_markers), &
             mark_age(max_markers), &
             mark_dead(max_markers), &
             mark_ntriag(max_markers), &
             mark_phase(max_markers), &
             mark_ID(max_markers), &
             mark_zpressn(max_markers), &
             mark_Emeltm(max_markers))

    allocate(mark_id_elem(max_markers_per_elem, nz-1, nx-1))
    allocate(nmark_elem(nz-1, nx-1))

    return
  end subroutine

  
  subroutine add_marker(x, y, iph, xxpres, xxmelt, age, j, i, inc)
    !$ACC routine seq
    !$ACC routine(check_inside) seq
    ! Add a marker at physical coordinate (x, y), with phase iph and age, to
    ! element (j, i). The current (before adding thsi marker) marker size
    ! is kk. If (x, y) is not within the element, inc is set to 0 and
    ! marker not added. Otherwise, marker is added to "mark" array and kk
    ! incremented by 1.

    use arrays
    use params
    implicit none

    integer :: iph, j, i, inc
    double precision :: x, y, age, xxpres, xxmelt
    integer :: ntr, kk, nm
    double precision :: bar1, bar2

    call check_inside(x , y, bar1, bar2, ntr, i, j, inc)
    if(inc.eq.0) return

    !$OMP atomic capture
    !$ACC atomic capture
    nmarkers = nmarkers + 1
    kk = nmarkers
    !$ACC end atomic
    !$OMP end atomic

    !$OMP atomic capture
    !$ACC atomic capture
    nmark_elem(j,i) = nmark_elem(j,i) + 1
    nm = nmark_elem(j,i)
    !$ACC end atomic
    !$OMP end atomic

    if(nm > max_markers_per_elem .or. kk > max_markers) then
        !write(msg*) 'Too many markers at element:', i, j, nm
        !call SysMsg(msg)
        call SysMsg('Marker skipped, not added!')
          
        !$OMP atomic update
        !$ACC atomic update
        nmarkers = nmarkers - 1

        !$OMP atomic update
        !$ACC atomic update
        nmark_elem(j,i) = nmark_elem(j,i) - 1
    
        inc = -1
        return
    endif

    ! recording the id of markers belonging to the element
    mark_id_elem(nm,j,i) = kk

    mark_x(kk) = x
    mark_y(kk) = y
    mark_dead(kk) = 1
    mark_ID(kk) = kk
    mark_a1(kk) = bar1
    mark_a2(kk) = bar2
    mark_age(kk) = age
    mark_ntriag(kk) = ntr
    mark_phase(kk) = iph
    mark_zpressn(kk) = xxpres
    mark_Emeltm(kk) = xxmelt

    return
  end subroutine add_marker


  subroutine newphase2marker (j1, j2, i1, i2, iph)
    !$ACC routine worker
    use arrays
    use params
    implicit none
    
    integer :: j1, j2, i1, i2, iph, &
               kk, n, j, i

    !$ACC loop collapse(2)
    !$OMP parallel do private(i,j,n,kk)
    ! reset the markers within elements in the rectangular region
    do i = i1, i2
      do j = j1, j2

        do n = 1 , nmark_elem(j,i)
          kk = mark_id_elem(n,j,i)
          mark_phase(kk) = iph
        enddo

        iphase(j,i) = iph
        phase_ratio(:,j,i) = 0.d0
        phase_ratio(iph,j,i) = 1.d0

      enddo
    enddo
    !$OMP end parallel do

    return
  end subroutine newphase2marker


  subroutine count_phase_ratio(j, i)
    !$ACC routine seq
    use params
    use arrays

    integer, intent(in) :: j, i
    integer :: n, kk, ncounters(maxph), iph, nm

    if (nmark_elem(j,i) == 0) stop 'No markers in element'

    ncounters = 0
    do n = 1, nmark_elem(j,i)
      kk = mark_id_elem(n,j,i)
      iph = mark_phase(kk)
      ncounters(iph) = ncounters(iph) + 1
    enddo

    do kk = 1, nphase
      phase_ratio(kk,j,i) = ncounters(kk) / (nmark_elem(j,i) * 1.0d0)
    enddo

    ! the phase of this element is the most abundant marker phase
    iph = 1
    nm = 0
    do kk = 1, nphase
       if (ncounters(kk) > nm) then
           nm = ncounters(kk)
           iph = kk
       end if
    end do
    iphase(j,i) = iph

    return
  end subroutine count_phase_ratio

END MODULE marker_data
