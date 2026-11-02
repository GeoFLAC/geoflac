
subroutine init_visc
  !$ACC routine(Eff_visc) seq
  use arrays
  use params
  implicit none
  integer :: i, j
  double precision :: Eff_visc
  !$OMP parallel do
  !$ACC parallel loop collapse(2) async(1)
  do i = 1,nx-1
      do j = 1,nz-1
          visn(j,i) = Eff_visc(j,i)
      end do
  end do
  !$OMP end parallel do
  return
end subroutine init_visc
