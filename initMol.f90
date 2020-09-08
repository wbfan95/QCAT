subroutine initMolecule
use para
implicit none

    write(*,*); write(*,*) 'Initializing molecules...'
    
    allocate( coordA(3,NA), coordB(3,NB) )
    allocate( eleMassA(NA), eleMassB(NB) )
    allocate( eleNameA(NA), eleNameB(NB) )
    allocate( pA(3,NA), pB(3,NB) )
    
    call splitCoord
    call ABremoveCOM ! shift center of mass to origin
    call ABremoveRot ! remove rotation by inertia tensor
    call saveCoord
    call initProp
    call calcR0 ! calculate the r0

end subroutine


subroutine saveCoord
use para
implicit none

    allocate(coordAInit(3,NA), coordBInit(3,NB))
    coordAInit(1:3,1:NA) = coordA(1:3,1:NA)
    coordBInit(1:3,1:NB) = coordB(1:3,1:NB)

end subroutine


subroutine ABremoveRot
use para
implicit none

    write(*,*) 'Rotate to principal axis...'
    call removeRotate(NA, coordA, eleMassA, ilinearA)
    if (ilinearA) write(*,*) ' Molecule A is linear! '
    call removeRotate(NB, coordB, eleMassB, ilinearB)
    if (ilinearB) write(*,*) ' Molecule B is linear! '
    
    if (ibug) then
        write(*,*) ' Molecule A : '
        call logXYZ(NA, coordA, eleNameA)
        write(*,*) ' Molecule B : '
        call logXYZ(NB, coordB, eleNameB)
    end if

end subroutine

! subroutine : calculate R0
! ref - ANT 2019 manual P71
! Method - 2.2 * (diameter A + diameter B)
subroutine calcR0
use para
implicit none
real(8) :: maxA, maxB

! find the maximum of each components (x,y,z)
! then find the MAX of these miaximums
! the final MAX is the diameter
! 2020-08-25 12:19:32 Wenbin, FAN @ SHU
    call calcMaxR(NA, coordA, maxA)
    call calcMaxR(NB, coordB, maxB)
    r0 = (maxA + maxB) * 2.2d0
    write(*,"(' R0 (Ang): ', F10.6, ' , dA : ', F10.6, ' , dB: ', F10.6)") r0*b2a, maxA*b2a, maxB*b2a

end subroutine


subroutine calcMaxR(Natoms, coord, r)
implicit none
integer, intent(in) :: Natoms
real(8), intent(in) :: coord(3,Natoms)
real(8), intent(out) :: r
integer i
real(8) comp(3) ! max in each component

    do i = 1, 3 ! x,y,z
        comp(i) = maxval( coord(i,1:Natoms) )
    end do
    r = maxval( comp(1:3) )

end subroutine 


subroutine removeRotate(Natoms, coord, eleMass, linear)
implicit none
integer, intent(in) :: Natoms
real(8), intent(inout) :: coord(3, Natoms)
real(8), intent(in) :: eleMass(Natoms)
logical, intent(out) :: linear
real(8) :: IT(3,3) ! xyz, xyz

    call computeInertiaTensor(Natoms, coord, eleMass, IT)
    call rotateToPrincipalAxis(Natoms, coord, IT, linear)

end subroutine


subroutine rotateToPrincipalAxis(Natoms, coord, IT, linear)
implicit none
integer, intent(in) :: Natoms
real(8), intent(inout) :: coord(3, Natoms)
real(8), intent(in) :: IT(3,3) ! xyz, xyz
logical, intent(out) :: linear
real(8) :: ITev(3,3), ITevalue(3)
integer i

    call diagMat(3, IT, ITevalue, ITev)
    do i = 1, Natoms
        coord(1:3,i) = matmul( coord(1:3,i), ITev(1:3, 1:3) )
    end do
    !write(*,*) ITevalue
    if ( minval( abs(ITevalue) ) .lt. 1D-5 ) then
        linear = .true.
    else
        linear = .false.
    end if

end subroutine


subroutine computeInertiaTensor(Natoms, coord, eleMass, IT)
implicit none
integer, intent(in) :: Natoms
real(8), intent(in) :: coord(3, Natoms)
real(8), intent(in) :: eleMass(Natoms)
real(8), intent(out) :: IT(3,3) ! xyz, xyz
integer i
real(8) :: x,y,z,m

    it(1:3,1:3) = 0d0
    do i = 1, Natoms
        m = eleMass(i)
        x = coord(1,i); y = coord(2,i); z = coord(3,i)
        IT(1,1) = IT(1,1) + m * (y*y + z*z)
        IT(2,1) = IT(2,1) - m * x*y
        IT(3,1) = IT(3,1) - m * x*z
        IT(1,2) = IT(1,2) - m * y*x
        IT(2,2) = IT(2,2) + m * (x*x + z*z)
        IT(3,2) = IT(3,2) - m * y*z
        IT(1,3) = IT(1,3) - m * z*x
        IT(2,3) = IT(2,3) - m * z*y
        IT(3,3) = IT(3,3) + m * (x*x + y*y)
    end do

end subroutine


! subroutine : shift A, B's COM to origin seperately
subroutine ABremoveCOM
use para
implicit none
integer :: i

    write(*,*) 'Shift center of mass to origin...'
    call removeCOM(NA, coordA, eleMassA)
    call removeCOM(NB, coordB, eleMassB)
    
    if (ibug) then
        write(*,*) ' Molecule A : '
        call logXYZ(NA, coordA, eleNameA)
        write(*,*) ' Molecule B : '
        call logXYZ(NB, coordB, eleNameB)
    end if

end subroutine

! subroutine : shift COM of molecule to origin
subroutine removeCOM(Natoms, coord, eleMass)
implicit none
integer, intent(in) :: Natoms
real(8), intent(inout) :: coord(3, Natoms)
real(8), intent(in) :: eleMass(Natoms)
real(8) :: COM(3) ! xyz
integer :: i

    call getCOM(Natoms, coord, eleMass, COM)
    do i = 1, 3
        coord(i,1:Natoms) = coord(i,1:Natoms) - COM(i)
    end do

end subroutine


subroutine getCOM(Natoms, coord, eleMass, COM)
implicit none
integer, intent(in) :: Natoms
real(8), intent(in) :: coord(3, Natoms)
real(8), intent(in) :: eleMass(Natoms)
real(8), intent(out) :: COM(3)
real(8) totalMass
integer i

    totalMass = sum( eleMass )
    do i = 1, 3
        COM(i) = sum( coord(i,1:Natoms) * eleMass(1:Natoms) ) / totalMass
    end do

end subroutine

subroutine splitCoord
use para
implicit none
integer :: i, j

    do i = 1, NA
        coordA(1:3,i) = coord(1:3, order(i))
        pA(1:3,i) = p(1:3, order(i))
        eleMassA(i) = eleMass( order(i) )
        eleNameA(i) = eleName( order(i) )
    end do
    
    do j = 1, NB ! order in B molecule
        i = j + NA ! real order in PES rank
        coordB(1:3,j) = coord(1:3, order(i))
        pB(1:3,i) = p(1:3, order(i))
        eleMassB(j) = eleMass( order(i) )
        eleNameB(j) = eleName( order(i) )
    end do

end subroutine


subroutine combineCoord
use para
implicit none
integer i

    do i = 1, NA
        coord(1:3, order(i) ) = coordA(1:3, i)
        p(1:3, order(i) ) = pA(1:3, i)
    end do
    do i = 1, NB
        coord(1:3, order(i+NA) ) = coordB(1:3, i)
        p(1:3, order(i+NA) ) = pB(1:3, i)
    end do

end subroutine


subroutine initProp
use para
implicit none


    massA = sum( eleMassA(1:NA) )
    massB = sum( eleMassB(1:NB) )
    muAB = massA * massB / (massA + massB)
    
    ! rotational inertia
    call calcRI(NA, coordA, eleMassA, RIA)
    call calcRI(NB, coordB, eleMassB, RIB)
    write(*,"(' Rotational inertia of A,B (au·Ang^2): ' 2F10.2)") RIA*b2a*b2a, RIB*b2a*b2a
    
    ! degree of freedom
    if (NA .ge. 2) then
        DOFA = 3 * NA - 3
    else
        DOFA = 1
    end if
    
    if (NB .ge. 2) then
        DOFB = 3 * NB - 3
    else
        DOFB = 1
    end if
    DOF = 3 * Natoms - 3
    write(*,"(' Degree of freedom (DOF) of A,B : ', 2I5)") DOFA, DOFB 

end subroutine


subroutine calcRI(Natoms, coord, eleMass, RI)
implicit none
integer, intent(in) :: Natoms
real(8), intent(in) :: coord(3, Natoms)
real(8), intent(in) :: eleMass(Natoms)
real(8), intent(out) :: RI
integer i
real(8) x,y,z, m, r

    RI = 0d0
    do i = 1, Natoms
        x = coord(1,i); y = coord(2,i); z = coord(3,i)
        m = eleMass(i)
        RI = RI + (x*x + y*y + z*z) * m
    end do

end subroutine