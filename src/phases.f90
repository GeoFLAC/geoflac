! -*- F90 -*-
module phases
!! Phases the got special treatment in the code

    integer, parameter :: kocean0 = 1
    integer, parameter :: kcont1 = 2
    integer, parameter :: kocean1 = 3
    integer, parameter :: kmant1 = 4
    integer, parameter :: kschist = 5
    integer, parameter :: kcont2 = 6
    integer, parameter :: kocean2 = 7
    integer, parameter :: kmant2 = 8
    integer, parameter :: kserp = 9
    integer, parameter :: ksed1 = 10
    integer, parameter :: ksed2 = 11
    integer, parameter :: kweak = 12
    integer, parameter :: keclg = 13
    integer, parameter :: karc1 = 14
    integer, parameter :: kweakmc = 15
    integer, parameter :: khydmant = 16
    integer, parameter :: kmetased = 17
    integer, parameter :: ksills = 18
    integer, parameter :: khtmsz = 19
    integer, parameter :: kmeltlc = 20
    integer, parameter, dimension(6) :: mantle_phases = (/kmant1, kmant2, kserp, khtmsz, ksills, kocean1/)
    integer, parameter, dimension(7) :: surface_phases = (/kcont1, kcont2, kweakmc, ksed1, ksed2, kocean2, kocean0/)
    
end module phases
