! -*- F90 -*-

module arrays
  implicit none

  save

  ! fortran array allocatable
  double precision, allocatable:: cord(:,:,:), temp(:,:), vel(:,:,:), stress0(:,:,:,:), &
       force(:,:,:), amass(:,:), rmass(:,:), &
       area(:,:,:), dvol(:,:,:), strain(:,:,:), bc(:,:,:), &
       xoriginal(:,:), zoriginal(:,:)

  integer, allocatable :: ncod(:,:,:), iphase(:,:), &
      nopbou(:,:), ncodbou(:,:), jmoho(:)

  double precision, allocatable :: phase_ratio(:,:,:), &
      dtopo(:), dhacc_correction(:), dhacc(:), extrusion(:), &
      fmelt(:,:), fmelt2(:,:), fmagma(:,:), fmagma2(:,:), extr_acc(:), &
      flux(:,:,:,:), &
      aps(:,:),visn(:,:),e2sr(:,:), &
      temp0(:,:),source(:,:),shrheat(:,:), &
      bcstress(:,:), q_init(:), vel_flex_old(:)

  double precision, allocatable :: se2sr(:,:,:), sshrheat(:,:)

  ! Dike-intrusion eigenstrain tracking, arc pathway (fmelt/fmagma).
  ! dike_added: this-timestep's raw accumulation from migration (scratch,
  !   zeroed and consumed within a single fl_therm.f90 call).
  ! dike_backlog: persistent, banked volumetric strain not yet released
  !   (rate-capped at 10% of Amin per timestep); survives remeshing.
  ! dike_released: this-timestep's actual released eigenstrain, computed by
  !   fl_therm.f90 and consumed by fl_rheol.f90 the same timestep (scratch).
  ! dike_marker_vol: persistent accumulator of released volume, until enough
  !   has arrived to add a new marker; survives remeshing.
  double precision, allocatable :: dike_added(:,:), dike_backlog(:,:), &
      dike_released(:,:), dike_marker_vol(:,:)

  ! Same idea for the MOR pathway (fmelt2/fmagma2), but no rate cap/backlog:
  ! just accumulates fmagma2's own additions until enough for a new marker.
  double precision, allocatable :: mor_marker_vol(:,:)

  ! remeshing arrays
  double precision, allocatable :: pt(:,:,:), barcord(:,:,:), &
            cold(:,:,:), cnew(:,:,:), &
            cordo(:,:,:), dhnew(:), extnew(:)
  integer, allocatable :: numtr(:,:)

  ! temporary array
  double precision, allocatable :: dummyn(:,:), dummye(:,:), stmpn(:)
  integer, allocatable :: itmp(:,:)

contains

  subroutine allocate_arrays(nz, nx, nphase)
    implicit none

    integer, intent(in) :: nz, nx, nphase
    integer, parameter :: ntri = 4    ! triangles per element

    allocate(cord(nz, nx, 2))
    allocate(temp(nz, nx))
    allocate(vel(nz, nx, 2))
    allocate(stress0(nz-1, nx-1, 4, ntri))
    allocate(force(nz, nx, 2))
    allocate(amass(nz, nx))
    allocate(rmass(nz, nx))
    allocate(area(nz-1, nx-1, ntri))
    allocate(dvol(nz-1, nx-1, ntri))
    allocate(strain(nz-1, nx-1, 3))
    allocate(bc(nz, nx, 2))
    allocate(xoriginal(nz, nx))
    allocate(zoriginal(nz, nx))

    allocate(ncod(nz, nx, 2))
    allocate(iphase(nz-1, nx-1))
    allocate(nopbou(((nz-1)+(nx-1))*2, 4))
    allocate(ncodbou(((nz-1)+(nx-1))*2, 3))
    allocate(jmoho(nx-1))

    allocate(phase_ratio(nphase, nz-1, nx-1))
    allocate(dtopo(nx))
    allocate(dhacc(nx-1))
    allocate(dhacc_correction(nx))
    allocate(extrusion(nx-1))
    allocate(fmelt(nz-1, nx-1))
    allocate(fmelt2(nz-1, nx-1))
    allocate(fmagma(nz-1, nx-1))
    allocate(fmagma2(nz-1, nx-1))
    allocate(extr_acc(nx-1))
    allocate(flux(2, 2, nz-1, nx-1))
    allocate(aps(nz-1, nx-1))
    allocate(visn(nz-1, nx-1))
    allocate(e2sr(nz-1, nx-1))
    allocate(temp0(nz, nx))
    allocate(source(nz-1, nx-1))
    allocate(shrheat(nz-1, nx-1))
    allocate(bcstress(((nz-1)+(nx-1))*2, 3))

    allocate(se2sr(nz-1, nx-1, 3))
    allocate(sshrheat(nz-1,nx-1))

    allocate(pt((nz-1)*(nx-1)*2, 2, 3), barcord(nz, nx, 3), &
             cold(nz, nx, 2), cnew(nz, nx, 2))
    allocate(numtr(nz, nx))
    allocate(cordo(nz,nx,2))
    allocate(dhnew(nx-1), extnew(nx-1))

    ! tmp arrays used in subroutines
    allocate(dummyn(nz, nx))
    allocate(dummye(nz-1, nx-1))
    allocate(stmpn(max(nx,nz)))
    allocate(itmp(nz, nx))
    allocate(q_init(nx))
    allocate(vel_flex_old(nx))
    vel_flex_old = 0.d0

    ! dike-intrusion eigenstrain/marker tracking
    allocate(dike_added(nz-1, nx-1))
    allocate(dike_backlog(nz-1, nx-1))
    allocate(dike_released(nz-1, nx-1))
    allocate(dike_marker_vol(nz-1, nx-1))
    allocate(mor_marker_vol(nz-1, nx-1))
    dike_added = 0.d0
    dike_backlog = 0.d0
    dike_released = 0.d0
    dike_marker_vol = 0.d0
    mor_marker_vol = 0.d0

  end subroutine allocate_arrays

end module arrays
