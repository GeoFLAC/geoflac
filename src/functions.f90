!================================
! First invariant of strain
function strainI(iz,ix)
use arrays
include 'precision.inc'

strainI = 0.5d0 * ( strain(iz,ix,1) + strain(iz,ix,2) )

return
end function strainI


!================================
! Second invariant of strain
function strainII(iz,ix)
use arrays
include 'precision.inc'

strainII = 0.5d0 * sqrt((strain(iz,ix,1)-strain(iz,ix,2))**2 + 4*strain(iz,ix,3)**2)

return
end function strainII


!================================
! Second invariant of strain rate
function srateII(iz,ix)
use arrays
include 'precision.inc'

s11 = 0.25d0 * (strainr(1,1,iz,ix)+strainr(1,2,iz,ix)+strainr(1,3,iz,ix)+strainr(1,4,iz,ix))
s22 = 0.25d0 * (strainr(2,1,iz,ix)+strainr(2,2,iz,ix)+strainr(2,3,iz,ix)+strainr(2,4,iz,ix))
s12 = 0.25d0 * (strainr(3,1,iz,ix)+strainr(3,2,iz,ix)+strainr(3,3,iz,ix)+strainr(3,4,iz,ix))
srateII = 0.5d0 * sqrt((s11-s22)**2 + 4*s12*s12)

return
end function srateII


!================================
! First invariant of stress (pressure)
function stressI(iz,ix)
use arrays
include 'precision.inc'

s11 = 0.25d0 * (stress0(iz,ix,1,1)+stress0(iz,ix,1,2)+stress0(iz,ix,1,3)+stress0(iz,ix,1,4))
s22 = 0.25d0 * (stress0(iz,ix,2,1)+stress0(iz,ix,2,2)+stress0(iz,ix,2,3)+stress0(iz,ix,2,4))
s33 = 0.25d0 * (stress0(iz,ix,4,1)+stress0(iz,ix,4,2)+stress0(iz,ix,4,3)+stress0(iz,ix,4,4))
stressI = (s11+s22+s33)/3d0

return
end function stressI


!================================
! Second invariant of stress
function stressII(iz,ix)
use arrays
include 'precision.inc'

s11 = 0.25d0 * (stress0(iz,ix,1,1)+stress0(iz,ix,1,2)+stress0(iz,ix,1,3)+stress0(iz,ix,1,4))
s22 = 0.25d0 * (stress0(iz,ix,2,1)+stress0(iz,ix,2,2)+stress0(iz,ix,2,3)+stress0(iz,ix,2,4))
s12 = 0.25d0 * (stress0(iz,ix,3,1)+stress0(iz,ix,3,2)+stress0(iz,ix,3,3)+stress0(iz,ix,3,4))
s33 = 0.25d0 * (stress0(iz,ix,4,1)+stress0(iz,ix,4,2)+stress0(iz,ix,4,3)+stress0(iz,ix,4,4))
stressII = 0.5d0 * sqrt((s11-s22)**2 + 4*s12*s12)

return
end function stressII


!==================================================
! Get the largest eigenvalue and its eigenvector (with its x-component
! fixed to 1) of the deviatoric of a symmetric 2x2 matrix
subroutine eigen2x2(a11, a22, a12, eigval1, eigvecy1)
include 'precision.inc'

adif = 0.5d0 * (a11 - a22)
eigval1 = sqrt(adif**2 + a12**2)
eigvecy1 = (eigval1 - adif) / a12

end subroutine eigen2x2


  !==================================================
subroutine SysMsg( message )
use params
include 'precision.inc'
character* (*) message

!$OMP critical (sysmsg1)
open( 13, file='sys.msg', position="append", action="write" )
write(13, * ) "Loops:", nloop, "Time[Ma]:", time/sec_year/1.d+6
write(13, * ) message
close(13)
!$OMP end critical (sysmsg1)

return
end

!===========================================
! calcsolid
!===========================================
function calcsolid(zP)
use params
implicit none
double precision :: calcsolid, zP
double precision :: dumtm, zPs, dPdTm, tint, Tsm, deltam, dutm

! Hirschmann et al (2000) & Hirschmann et al (2009) 
if (isolidus.eq.2) then
    if (zP.le.10.) then
            calcsolid = -5.140*(zP**2) + 132.899*zP + 1120.661
    elseif (zP.le.23.5) then
            calcsolid = -1.092*((zP-10.)**2) + 32.390*(zP-10.) + 1935.
    else
            calcsolid = 26.35*(zP-23.5) + 2175
    end if

! Herzberg et al (2000) 
else if (isolidus.eq.3) then
    if (zP.le.2.7) then
            calcsolid = 132.*(zP) + 1102.
    else
            calcsolid = 390.*(log(zP)) - 5.7*(zP) + 1086.
    end if

! McKenzie & Bickle (1988) Solidus
! else if (isolidus.eq.4) then
!     dumtm = 1000.
!     zPs = ((dumtm-1100.)/136.)+(0.0004968*exp(0.012*(dumtm-1100.)))
!     dPdTm = (0.0004968*0.012*exp(0.012*(dumtm-1100.)) + (1/136.))
!     tint = zPs - (dPdTm*dumtm)
!     Tsm = (zP-tint)/dPdTm
!     deltam = abs(Tsm-dumtm)
!     dumtm = Tsm
!     do while (deltam.ge.0.01)
!         zPs = ((dumtm-1100.)/136.)+(0.0004968*exp(0.012*(dumtm-1100.)))
!         dPdTm = (0.0004968*0.012*exp(0.012*(dumtm-1100.)) + (1/136.))
!         tint = zPs - (dPdTm*dutm)
!         Tsm = (zP-tint)/dPdTm
!         deltam = abs(Tsm-dumtm)
!         dumtm = Tsm
!     end do

!     calcsolid = dumtm

! Katz et al (2003) DEFAULT
else
    calcsolid = (-5.104)*(zP**2) + 132.899*(zP) + 1085.7
    !calcsolid = 1000.
end if

return
end function calcsolid

!===========================================
! calcliq
!===========================================
function calcliq(zzP)
use params
implicit none
double precision :: calcliq, zzP

! McKenzie and Bickle liquidus
if (isolidus.eq.4) then
    calcliq =  1736.2 + 4.343*zzP + 180.*atan(zzP/2.2169)

! Katz 2003 DEFAULT
else
    calcliq =  1780. + 45.*zzP - 2.*(zzP**2)
end if

return
end function calcliq
    
!===========================================
! DeltaTwat
!===========================================
function deltaTwat(X,xP,F)
implicit none
double precision :: deltaTwat, X, xP, F
double precision :: xm, xsat

xsat = 0.12*(xP**0.6)+0.01*xP !Eqn 17 Katz2003
xm = X /(0.01 + F*(1.0-0.01)) !Eqn 18 Katz2003
if (xm.ge.xsat) then !Eqn 15 Katz2003
    xm = xsat
end if

deltaTwat = 43.0*((xm*100.0)**(0.75)) !Eqn 16 Katz 2003

return
end function deltaTwat
    
!===========================================
! Calc F
!===========================================
function CalcF(xxP,T,X,CpxM)
use arrays
implicit none
double precision :: CalcF, xxP, T, X, CpxM
double precision, external :: deltaTwat, calcsolid, calcliq
double precision :: Rcpx, Fcpxout, Tsolf, Tlhliq, Tliqf, Tcpxout, xDT, guessF, altguessF, deltaF
double precision :: Twm, Twc, Twl, dumFF, Tws
integer :: n


Tsolf = calcsolid(xxP) 
Tliqf = calcliq(xxP)
Rcpx = 0.5 + 0.08*xxP !Eqn 7 Katz03
Fcpxout = CpxM/Rcpx !Eqn 6 Katz03
Tlhliq = 1475. + 80.*xxP - 3.2*(xxP**2) !Eqn5 Katz03
Tcpxout = (Fcpxout**(2./3.))*(Tlhliq-Tsolf) + Tsolf !Eqn9 Katz03
dumFF = 0.0


if (X.gt.0) then
    xDT = deltaTwat(X,xxP,dumFF)
    Tws = Tsolf - xDT !Eqn11 Katz03
    xDT = deltaTwat(X,xxP,Fcpxout)
    Twm = Tlhliq - xDT !Eqn12 Katz03
    Twc = Tcpxout - xDT
    dumFF = 1.0
    xDT = deltaTwat(X,xxP,dumFF)
    Twl = Tliqf-xDT !Eqn13 Katz03
else
    Tws = Tsolf
    Twm = Tlhliq
    Twc = Tcpxout
    Twl = Tliqf
end if 

if (T.lt.Tws) then
    CalcF = 0.0
end if

if (T.ge.Tws .AND. T.lt.Twc) then
    if (X.eq.0) then
        CalcF = ((T-Tsolf)/(Tlhliq-Tsolf))**(1.5) !Eqn2&3 Katz03
    else
        guessF = ((T-Tws)/(Twm-Tws))**(1.5)
        xDT = deltaTwat(X,xxP,guessF)

        do while ((Tsolf-xDT).ge.T)                
            guessF = guessF*(0.9)
            xDT = deltaTwat(X,xxP,guessF)
        end do

        altguessF = ((T-(Tsolf-xDT))/(Tlhliq-Tsolf))**(1.5)
        deltaF = abs(guessF - altguessF)
        
        n = 0
        !Bisection method to solve for F
        do while (deltaF.gt.1e-7 .AND. n.le.40)
            n = n+1

            if (guessF.lt.altguessF) then
                guessF = guessF+(deltaF*0.25)
                xDT = deltaTwat(X,xxP,guessF)
                
                do while ((Tsolf-xDT).ge.T)
                    guessF = guessF*(0.95)
                    xDT = deltaTwat(X,xxP,guessF)
                end do
                
                altguessF = ((T-(Tsolf-xDT))/(Tlhliq-Tsolf))**(1.5)
                deltaF = abs(guessF-altguessF)
            else
                guessF = guessF-(deltaF*0.25)     
                xDT = deltaTwat(X,xxP,guessF)                           
                do while ((Tsolf-xDT).ge.T)
                    guessF = guessF*(0.95)
                    xDT = deltaTwat(X,xxP,guessF)
                end do

                altguessF = ((T-(Tsolf-xDT))/(Tlhliq-Tsolf))**(1.5)
                deltaF = abs(guessF-altguessF)
            end if
        end do
        
        CalcF = guessF
    end if
end if

if (T.ge.Twc .AND. T.le.Twl) then
    if (X.eq.0) then
        CalcF = Fcpxout + (1.-Fcpxout)*(((T-Tcpxout)/(Tliqf-Tcpxout))**(1.5)) !Eqn8 Katz03
    else
        guessF = Fcpxout + (1.-Fcpxout)*(((T-Twc)/(Twl-Twc))**(1.5))
        xDT = deltaTwat(X,xxP,guessF)

        do while ((Tcpxout-xDT).ge.T)
            guessF = guessF*(0.95)
            xDT = deltaTwat(X,xxP,guessF)
        end do

        altguessF = Fcpxout + (1.-Fcpxout)*(((T-(Tcpxout-xDT))/(Tliqf-Tcpxout))**(1.5))
        deltaF = abs(guessF - altguessF)

        n = 0
        do while (deltaF.gt.1e-7 .AND. n.le.40)
            n = n+1

            if (guessF.lt.altguessF) then
                guessF = guessF+(deltaF*0.25)
                xDT = deltaTwat(X,xxP,guessF)

                do while ((Tcpxout-xDT).ge.T)
                    guessF = guessF*(0.95)
                    xDT = deltaTwat(X,xxP,guessF)
                end do

                altguessF = Fcpxout + (1.-Fcpxout)*(((T-(Tcpxout-xDT))/(Tliqf-Tcpxout))**(1.5))
                deltaF = abs(guessF-altguessF)
            else
                guessF = guessF-(deltaF*0.25)
                xDT = deltaTwat(X,xxP,guessF)

                do while ((Tcpxout-xDT).ge.T)
                    guessF = guessF*(0.95)
                    xDT = deltaTwat(X,xxP,guessF)
                end do

                altguessF = Fcpxout + (1.-Fcpxout)*(((T-(Tcpxout-xDT))/(Tliqf-Tcpxout))**(1.5))
                deltaF = abs(guessF-altguessF)
            end if
        end do
        CalcF = guessF
    end if
end if

if (T.gt.Twl) then
    CalcF = 1.0
end if

return
end function CalcF
    
!===========================================
!CalcT
!===========================================
function CalcT(xPx,X,CpxM,F)
use arrays
implicit none
double precision :: Rcpx, Fcpxout, Tlhliq, Tcpxout, xDT, xT, xPx, X, CpxM, F, dumFF
double precision :: Tsolf, Tliqf, CalcT
double precision, external :: deltaTwat, calcsolid, calcliq

Tsolf = calcsolid(xPx)
Tliqf = calcliq(xPx)
Rcpx = 0.5 + 0.08*xPx !Eqn 7 Katz03
Fcpxout = CpxM/Rcpx !Eqn 6 Katz03
Tlhliq = 1475.0 + 80.0*xPx - 3.2*(xPx**2.) !Eqn5 Katz03
Tcpxout = (Fcpxout**(2.0/3.0))*(Tlhliq-Tsolf) + Tsolf !Eqn9 Katz03

if (F.le.0.0) then
    xT = Tsolf
    if (X.gt.0) then
        dumFF = 0.0
        xDT = deltaTwat(X,xPx,dumFF)
        CalcT = xT-xDT
    else
        CalcT = xT
    end if
end if

if (F.gt.0.0 .AND. F.le.Fcpxout) then
    xT = (F**(2.0/3.0))*(Tlhliq-Tsolf)+Tsolf
    if (X.gt.0.) then
            xDT = deltaTwat(X,xPx,F)
            CalcT = xT-xDT
    else
            CalcT = xT
    end if
end if

if (F.gt.Fcpxout .AND. F.le.1.0) then
    xT = (((F-Fcpxout)/(1-Fcpxout))**(2.0/3.0))*(Tliqf-Tcpxout)+Tcpxout
    if (X.gt.0.) then
        xDT = deltaTwat(X,xPx,F)
        CalcT = xT-xDT
    else
        CalcT = xT
    end if
end if

if (F.gt.1.0) then
    xT = Tliqf
    if (X.gt.0) then
        dumFF = 1.0
        xDT = deltaTwat(X,xPx,dumFF)
        CalcT = xT-xDT
    else
        CalcT = xT
    end if
end if

return
end function CalcT 
    
!===============================
! Gaussian elemination
!===============================
subroutine Gauss(A,N,B)
include 'precision.inc'
dimension A(N,N),B(N),IPIV(N),INDXR(N),INDXC(N)

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
