!================================
! First invariant of strain
function strainI(iz,ix)
use arrays
implicit none
double precision :: strainI
integer :: iz, ix

strainI = 0.5 * ( strain(iz,ix,1) + strain(iz,ix,2) )

return
end function strainI


!================================
! Second invariant of strain
function strainII(iz,ix)
use arrays
implicit none
double precision :: strainII
integer :: iz, ix

strainII = 0.5 * sqrt((strain(iz,ix,1)-strain(iz,ix,2))**2 + 4*strain(iz,ix,3)**2)

return
end function strainII


!================================
! Second invariant of strain rate
function srateII(iz,ix)
use arrays
implicit none
double precision :: srateII
integer :: iz, ix
double precision :: s11, s22, s12

s11 = 0.25 * (strainr(1,1,iz,ix)+strainr(1,2,iz,ix)+strainr(1,3,iz,ix)+strainr(1,4,iz,ix))
s22 = 0.25 * (strainr(2,1,iz,ix)+strainr(2,2,iz,ix)+strainr(2,3,iz,ix)+strainr(2,4,iz,ix))
s12 = 0.25 * (strainr(3,1,iz,ix)+strainr(3,2,iz,ix)+strainr(3,3,iz,ix)+strainr(3,4,iz,ix))
srateII = 0.5 * sqrt((s11-s22)**2 + 4*s12*s12)

return
end function srateII


!================================
! First invariant of stress (pressure)
function stressI(iz,ix)
use arrays
implicit none
double precision :: stressI
integer :: iz, ix
double precision :: s11, s22, s33

s11 = 0.25 * (stress0(iz,ix,1,1)+stress0(iz,ix,1,2)+stress0(iz,ix,1,3)+stress0(iz,ix,1,4))
s22 = 0.25 * (stress0(iz,ix,2,1)+stress0(iz,ix,2,2)+stress0(iz,ix,2,3)+stress0(iz,ix,2,4))
s33 = 0.25 * (stress0(iz,ix,4,1)+stress0(iz,ix,4,2)+stress0(iz,ix,4,3)+stress0(iz,ix,4,4))
stressI = (s11+s22+s33)/3

return
end function stressI


!================================
! Second invariant of stress
function stressII(iz,ix)
use arrays
implicit none
double precision :: stressII
integer :: iz, ix
double precision :: s11, s22, s12, s33

s11 = 0.25 * (stress0(iz,ix,1,1)+stress0(iz,ix,1,2)+stress0(iz,ix,1,3)+stress0(iz,ix,1,4))
s22 = 0.25 * (stress0(iz,ix,2,1)+stress0(iz,ix,2,2)+stress0(iz,ix,2,3)+stress0(iz,ix,2,4))
s12 = 0.25 * (stress0(iz,ix,3,1)+stress0(iz,ix,3,2)+stress0(iz,ix,3,3)+stress0(iz,ix,3,4))
s33 = 0.25 * (stress0(iz,ix,4,1)+stress0(iz,ix,4,2)+stress0(iz,ix,4,3)+stress0(iz,ix,4,4))
stressII = 0.5 * sqrt((s11-s22)**2 + 4*s12*s12)

return
end function stressII


!==================================================
! Get the largest eigenvalue and its eigenvector (with its x-component
! fixed to 1) of the deviatoric of a symmetric 2x2 matrix
subroutine eigen2x2(a11, a22, a12, eigval1, eigvecy1)
implicit none
double precision :: a11, a22, a12, eigval1, eigvecy1
double precision :: adif

adif = 0.5 * (a11 - a22)
eigval1 = sqrt(adif**2 + a12**2)
eigvecy1 = (eigval1 - adif) / a12

end subroutine eigen2x2


  !==================================================
subroutine SysMsg( message )
use params
implicit none
character* (*) message

!$OMP critical (sysmsg1)
open( 13, file='sys.msg', position="append", action="write" )
write(13, * ) "Loops:", nloop, "Time[Ma]:", time/sec_year/1.e+6
write(13, * ) message
close(13)
!$OMP end critical (sysmsg1)

return
end


!=======================================
FUNCTION ERFCC(X)
implicit none
double precision :: ERFCC
double precision :: X, Z, T

Z=ABS(X)      
T=1./(1.+0.5*Z)
ERFCC=T*EXP(-Z*Z-1.26551223+T*(1.00002368+T*(.37409196+ &
      T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+ &
      T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
IF (X.LT.0.) ERFCC=2.-ERFCC

RETURN
END FUNCTION ERFCC


!===========================================
! Protected dlog10
function dlog10ab( val )
implicit none
double precision :: dlog10ab, val

if ( val .le. 0. ) then
    dlog10ab = 1.e-30
else
    dlog10ab = dlog10( val )
endif

return
end function dlog10ab


!===============================
! Gaussian elemination
!===============================
subroutine Gauss(A,N,B)
implicit none
integer :: N
double precision :: A(N,N),B(N)
integer :: IPIV(N),INDXR(N),INDXC(N)
double precision :: BIG, dum, pivinv
integer :: i, j, k, icol, irow, l, ll
ipiv = 0      

do i=1,N
    BIG=0.
    do j=1,N
        if( IPIV(j).ne.1 ) then
            do k=1,N
                if( IPIV(k).eq.0 ) then
                    if( abs(A(j,k)).gt.BIG ) then
                        BIG = abs( A(j,k) )
                        IROW=j
                        ICOL=k
                    endif
                else if( IPIV(k).gt.1 ) then
                    write(*,*) 'Singular matrix'
                    stop
                endif
            end do
        endif
    end do

    IPIV(ICOL) = IPIV(ICOL) + 1
    if( IROW.ne.ICOL ) then
        do l=1,N
            dum = A(IROW,l)
            A(IROW,l) = A(ICOL,l)
            A(ICOL,l) = dum
        end do
        dum = B(IROW)
        B(IROW)=B(ICOL)
        B(ICOL)=dum
    end if

    INDXR(i) = IROW
    INDXC(i) = ICOL
    if( A(ICOL,ICOL).eq.0. ) then
        write(*,*) 'Singular matrix.'
        stop
    endif
    pivinv = 1./A(ICOL,ICOL)
    A(ICOL,ICOL) = 1.
    do l=1,N
        A(ICOL,l) = A(ICOL,l) * pivinv
    end do
    B(ICOL) = B(ICOL) * pivinv
    do ll=1,N
        if( ll.ne.ICOL ) then
            dum = A(ll,ICOL)
            A(ll,ICOL) = 0.
            do l=1,N
                A(ll,l) = A(ll,l) - A(ICOL,l)*dum
            end do
            B(ll) = B(ll) - B(ICOL)*dum
        endif
    end do
end do

do l=N,1,-1
    if( INDXR(l).ne.INDXC(l) ) then
        do k=1,N
            dum = A(k,INDXR(l))
            A(k,INDXR(l)) = A(k,INDXC(l))
            A(k,INDXC(l)) = dum
        end do
    end if
end do


return
end subroutine Gauss
