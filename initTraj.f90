subroutine initTraj
use para
implicit none

    call loadInitMol
    call ABrandRot ! random rotation
    call BrandCoord ! random shift molecule B
    
    call ABrandThermal ! initial momentum
    call ABrotState ! rotational states
    call ABvibState ! vibrational states
    
    call ABremoveTranslation
    call ABremoveAngularMomentum
    
    call unloadInitMol

end subroutine


subroutine ABremoveAngularMomentum
use para
implicit none
real(inum) :: initTemp

    write(*,*) 'Removing overall angular momentum...'
    call splitCoord
    if (NA .gt. 1) call removeAngularMomentum(NA, coordA, pA, eleMassA)
    if (NB .gt. 1) call removeAngularMomentum(NB, coordB, pB, eleMassB)
    call combineCoord
    call getTemp(Natoms, eleMass, p, DOF, initTemp)
    write(*,'(A,F7.2)') ' Temperature (K) : ', initTemp * au2K

end subroutine


subroutine removeAngularMomentum(Natoms, coord, p, eleMass)
use constant
implicit none
integer :: Natoms
real(inum), intent(in) :: coord(3, Natoms)
real(inum), intent(inout) :: p(3, Natoms)
real(inum), intent(in) :: eleMass(Natoms)
integer i
real(inum) q(3,Natoms), COM(3)
real(inum) AM(3), J(3) ! angular momentum, sum AM
real(inum) IT(3,3), ITinv(3,3)
real(inum) omega(3), omegaDotq(3)

    ! remove center of mass
    call getCOM(Natoms, coord, eleMass, COM)
    do i = 1, 3
        q(i,1:Natoms) = coord(i,1:Natoms) - COM(i)
    end do
    
    ! calculate total angular momentum
    J(1:3) = 0d0
    do i = 1, Natoms
        call cross_product(q(1:3,i), p(1:3,i), AM(1:3))
        J(1:3) = J(1:3) + AM(1:3)
    end do
    !!!!! wrong code : calculate angular momentum
    !!!AM = matmul( q, p ) ! L = p\dot r, L_x = p_x x
    !!!do i = 1, 3
    !!!    J(i) = sum( AM(i,1:Natoms) )
    !!!end do
    
    call computeInertiaTensor(Natoms, q, eleMass, IT)
    call pinv(3, IT, ITinv)
    omega = matmul(ITinv, J)
    !#write(*,*) 'ITinv'
    !#write(*,'(3F12.6)') ITinv
    !#write(*,*) 'coord'
    !#write(*,'(3F12.6)') q
    do i = 1, Natoms
        call cross_product(omega(1:3), q(1:3,i), omegaDotq(1:3))
        p(1:3,i) = p(1:3,i) - eleMass(i) * omegaDotq(1:3)
    end do
    !#write(*,*) 'omegaDotq'
    !#write(*,'(3F12.6)') omegaDotq

end subroutine


subroutine ABremoveTranslation
use para
implicit none
real(inum) initTemp

    if(ibug) write(*,*) 'Momentum before removing translation (au) :'
    if(ibug) call logXYZ(Natoms, p, eleName)
    call splitCoord
    call removeTranslation(NA, pA)
    call removeTranslation(NB, pB)
    call combineCoord
    if(ibug) write(*,*) 'Momentum after removing translation (au) :'
    if(ibug) call logXYZ(Natoms, p, eleName)
    
    call getTemp(Natoms, eleMass, p, DOF, initTemp)
    write(*,'(A,F7.2)') ' Temperature (K) : ', initTemp * au2K

end subroutine


subroutine removeTranslation(Natoms, p)
use constant
implicit none
integer, intent(in) :: Natoms
real(inum), intent(inout) :: p(3, Natoms)
real(inum) COP(3) ! center of momentum
integer i

    do i = 1, 3
        COP(i) = sum( p(i,1:Natoms) )
        p(i,1:Natoms) = p(i,1:Natoms) - COP(i)
    end do

end subroutine


subroutine ABvibState
implicit none

    call combineCoord
    call vibInit
    call vibGen

end subroutine

! subroutine : generate initial vibrational states
! freq - Hartree, freq * au2cm1 = cm^{-1}
subroutine vibGen
use para
implicit none
real(inum) :: Evib(nMode), xi(nMode), N(nMode)
real(inum) :: EHO(nMode), rturn(nMode), sigmaX(nMode), sigmaP(nMode)
! displacement of coordinate and momentum, `ksi` is random number
real(inum) :: qvdis(3,Natoms), pvdis(3,Natoms), ksi(2,nMode)
real(inum) ppm ! momentum positive or minus, random number
integer i, j

    if (nMode .eq. 0) call logError('The number of vibrational mode is zero! ')
    call random_number(xi)
    do i = 1, nMode
        Evib(i) = - kB * temp * dlog( 1 - xi(i) ) / au2J
        N(i) = Evib(i)*au2J / (freq(i)*au2Hz) / hbar
        EHO(i) = hbar * freq(i) * au2Hz * (N(i) + 0.5d0) / au2J
        rturn(i) = dsqrt(2d0 * EHO(i)*au2J ) / (freq(i)*au2Hz) / dsqrt(au2kg) / au2m
        sigmaX(i) = dsqrt(hbar / (2d0 * au2kg * freq(i)*au2Hz)) / au2m
        sigmaP(i) = hbar / (2d0 * sigmaX(i)*au2m) / (au2kg*au2m/au2s)
    end do
    if (ibug) write(*,*) 'Vibrational energy (kcal/mol) : '
    if (ibug) write(*,'(9F7.2)') Evib * au2kcm
    if (ibug) write(*,*) 'Vibrational quantum number : '
    if (ibug) write(*,'(9F7.2)') N
    if (ibug) write(*,*) 'Harmonic energy (kcal/mol) : '
    if (ibug) write(*,'(9F7.2)') EHO * au2kcm
    if (ibug) write(*,*) 'Turning points (Ang) : '
    if (ibug) write(*,'(9F7.2)') rturn * b2a
    if (ibug) write(*,*) 'Variance of harmonic oscillator (Ang) :'
    if (ibug) write(*,'(9F7.2)') sigmaX * b2a
    if (ibug) write(*,*) 'Variance of momentum (au) *100:'
    if (ibug) write(*,'(9F7.2)') sigmaP * 100d0
    
    call random_number(ksi(1:2,1:nMode))
    do i = 1, Natoms
    do j = 1, 3
        qvdis(j,i) = sum( &
            fVec((i-1)*3+j,1:nMode) * sigmaX(1:nMode) * &
            dsqrt(-2d0 * dlog(ksi(1,1:nMode)) ) * &
            dcos(2d0 * pi * ksi(2,1:nMode)) &
        )
        pvdis(j,i) = sum( &
            fVec((i-1)*3+j,1:nMode) * & !1d0/(2d0*sigmaX(1:nMode)) * &
            sigmaP(1:nMode) * & 
            dsqrt(-2d0 * dlog(ksi(1,1:nMode)) ) * &
            dsin(2d0 * pi * ksi(2,1:nMode)) &
        )
    end do
    end do
    !!!! verify that all eigenvectors are normalized. 
    !!!do i = 1, nMode
    !!!    write(*,'(9999F7.2)') sum(fVec(:,i)*fVec(:,i))
    !!!end do
    !!!! print all eigenvectors
    !!!do i = 1, nMode
    !!!    write(*,'(9999F7.2)') fVec(:,i)
    !!!end do
    
    ! unscale the mass-scaled coordinate and momentum
    do i = 1, Natoms
        qvdis(1:3, i) = qvdis(1:3, i) / dsqrt(eleMass(i))
        pvdis(1:3, i) = pvdis(1:3, i) * dsqrt(eleMass(i))
    end do
    
    if (ibug) write(*,*) 'Vibrational displacement of coordinates (Ang) : '
    if (ibug) call logXYZ(Natoms, qvdis * b2a, eleName)
    if (ibug) write(*,*) 'Vibrational momentum (a.u.) : '
    if (ibug) call logXYZ(Natoms, pvdis, eleName)
    !!!#write(*,*) Natoms; write(*,*); 
    !!!#call logXYZ(Natoms, (coord + qvdis) * b2a, eleName)
    !!!#write(*,*) Natoms; write(*,*); 
    !!!#call logXYZ(Natoms, (coord) * b2a, eleName)
    
    do i = 1, Natoms
        do j = 1, 3
            call random_number(ppm)
            if ( ppm .gt. 0.5d0 ) then
                p(j,i) = p(j,i) + pvdis(j,i)
            else 
                p(j,i) = p(j,i) - pvdis(j,i)
            end if
        end do
    end do
    coord = coord + qvdis

end subroutine

! subroutine : calculate valid frequency and vibrational vectors
subroutine vibInit
use para
implicit none
real(inum) :: energy, grad(3,Natoms)
real(inum) :: evec(3*Natoms, 3*Natoms), frequency(3*Natoms)
real(inum) :: tmp, tmpVec(3*Natoms)
integer :: i,j
logical swapped
    
    write(*,*); write(*,*) 'Generate vibrational states...'
    if (ibug) call logXYZ(Natoms, coord*b2a, eleName)
    !#call pot(Natoms, coord, energy)
    call freqAna(Natoms, coord, eleMass, frequency, evec)
    !#write(*,*) 'test!'
    !#write(*,*) 'Wavenumber (cm^{-1}) : '
    !#write(*,'(6F9.2)') frequency * au2cm1
    
    ! sort the frequencies by absolute value
    ! to remove six smallest freq
    do i = 1, Natoms*3
        swapped = .false.
        do j = i+1, Natoms*3
            if (abs(frequency(i)) .gt. abs(frequency(j))) then
                tmp = frequency(i)
                frequency(i) = frequency(j)
                frequency(j) = tmp
                tmpVec(:) = evec(:,i)
                evec(:,i) = evec(:,j)
                evec(:,j) = tmpVec(:)
                swapped = .true. 
            end if
        end do
        if (.not. swapped) exit
    end do
    if (ibug) write(*,*) 'Small wavenumbers (cm^{-1}) : '
    if (ibug) write(*,'(6F9.2)') frequency(1:6) * au2cm1
    frequency(1:6) = -1d16 ! a huge negative value represents invalid
    
    ! sort frequencies from small to big
    ! to remove negative freq
    do i = 1, Natoms * 3
        swapped = .false.
        do j = i+1, Natoms*3
            if (frequency(i) .gt. frequency(j)) then
                tmp = frequency(i)
                frequency(i) = frequency(j)
                frequency(j) = tmp
                tmpVec(:) = evec(:,i)
                evec(:,i) = evec(:,j)
                evec(:,j) = tmpVec(:)
                swapped = .true. 
            end if
        end do
    end do
    write(*,*) 'Wavenumber (cm^{-1}, 3N-6) : '
    write(*,'(6F9.2)') frequency(7:) * au2cm1
    ! count the valid frequencies
    do i = 1, Natoms*3
        if ( frequency(i) .gt. 0d0 ) then
            nMode = Natoms * 3 - i + 1
            exit
        end if
    end do
    write(*,'(A,I4)') ' Effective number of modes : ', nMode
    write(*,'(A,F9.2)') ' Recommended time step (fs) : ', &
        1d0 / maxval(frequency) / au2cm1 / cLight * 1D12
    ! au -> cm-1 -> period
        
    !#write(*,'(6F9.2)') freq*au2cm1
    allocate( freq(nMode), fVec(3*Natoms,nMode) )
    freq(1:nMode) = frequency(3*Natoms-nMode+1:3*Natoms)
    fVec(1:Natoms*3,1:nMode) = evec(:,3*Natoms-nMode+1:3*Natoms)

end subroutine


subroutine ABrotState
use para
implicit none
real(8) :: rotpA(3,NA), rotpB(3,NB)
real(8) :: rotTemp

    write(*,*); write(*,*) 'Generate rotational states...'
    !write(*,*) 'Molecule A'
    if (NA .gt. 1) then
    if (ilinearA) then
        call linearRot(NA, coordA, eleMassA, tempAu, RIA, rotpA)
    else
        call nonLinearRot(NA, coordA, eleMassA, tempAu, RIA, rotpA)
    end if
    call getTemp(NA, eleMassA, rotpA, DOFA, rotTemp)
    write(*,'(A,F8.2)') ' A Rotational temperature (K) :', rotTemp * au2K
    end if
    !write(*,*) 'Molecule B'
    if (NB .gt. 1) then
    if (ilinearB) then
        call linearRot(NB, coordB, eleMassB, tempAu, RIB, rotpB)
    else
        call nonLinearRot(NB, coordB, eleMassB, tempAu, RIB, rotpB)
    end if
    call getTemp(NB, eleMassB, rotpB, DOFB, rotTemp)
    write(*,'(A,F8.2)') ' B Rotational temperature (K) :', rotTemp * au2K
    end if
    
    pA = pA + rotpA
    pB = pB + rotpB

end subroutine

! subroutine : generate rotational mometum for non-linear molecule
! temp - temperature in au
subroutine nonLinearRot(Natoms, coord, eleMass, temp, RI, p)
use constant
implicit none
integer, intent(in) :: Natoms
real(8), intent(in) :: coord(3,Natoms), eleMass(Natoms), temp, RI
real(8), intent(out) :: p(3,Natoms)
real(8) xi(3), Erot(3), omega(3)
! remove rotation and get principal axis
! q0 - COM removed, q - rotated to principal axis
real(8) q0(3,Natoms), q(3,Natoms), COM(3), IT(3,3) ,ITev(3,3), ITevalue(3), ITevInv(3,3)
integer i, j
real(8), external :: deg
real(8) theta ! angle in principal axis
!(3) ! random direction in three axis
real(8) pv ! vertical momentum
real(8) p1(3), p0(3) ! rotated momentum, unrotated momentum
integer anotherTwoAxis(2), ax1, ax2, ax3

    call random_number(xi(1:3))
    do i = 1, 3 ! xyz
        Erot(i) = - temp * dlog( 1 - xi(i) ) * 0.5d0 ! kB = 1
        omega(i) = dsqrt( 2d0 * Erot(i) / RI )
    end do
    write(*,'(A, 3F12.6)') ' Erot (kcal/mol) in x,y,z : ', Erot(1:3) * au2kcm
    if (ibug) write(*,'(A, 3F12.6)') ' Omega (rad/fs) in x,y,z : ', omega(1:3) / au2fs
    
    !!! remove COM
    q0(1:3, 1:Natoms) = coord(1:3, 1:Natoms)
    call getCOM(Natoms, q0, eleMass, COM)
    do i = 1, 3
        q0(i,1:Natoms) = q0(i,1:Natoms) - COM(i)
    end do
    !#if(ibug) write(*,*) 'Move to origin: '
    !#if(ibug) call logXYZ(Natoms, q0, 0)
    
    ! find the rotational axis
    call computeInertiaTensor(Natoms, q0, eleMass, IT)
    call diagMat(3, IT, ITevalue, ITev)
    call pinv(3, ITev, ITevInv)
    !#write(*,*); write(*,'(3F12.6)') ITev 
    !#write(*,*); write(*,'(3F12.6)') ITevInv
    !#write(*,*); write(*,'(3F12.6)') matmul( ITev, ITevInv )
    
    !#q(1:3,1:Natoms) = matmul( q(1:3,1:Natoms), ITev(1:3,1:3) )
    do i = 1, Natoms
        q(1:3,i) = matmul( q0(1:3,i), ITev(1:3, 1:3) )
    end do
    !#call logXYZ(Natoms, q, 0)
    
    !!!! random direction of momentum
    !!!call random_number(theta(1:3))
    !!!theta(1:3) = theta(1:3) * 2d0 * pi ! convert to 0~2pi
    !!!write(*,'(A,3F8.2)') ' Random angle (deg) in x,y,z : ', theta(1:3) / pi * 180d0
    
    ! generate momentum in each principal axis
    p = 0d0
    do j = 1, 3 ! three pincipal axis
        call getAnotherTwoAxis(j, anotherTwoAxis(1:2))
        !#write(*,*) anotherTwoAxis
        ax1 = j; ax2 = anotherTwoAxis(1); ax3 = anotherTwoAxis(2)
        do i = 1, Natoms
            p1(1:3) = 0d0; p0(1:3) = 0d0
            pv = eleMass(i) * dsqrt(q(ax2,i)*q(ax2,i) + q(ax3,i)*q(ax3,i)) * omega(ax1)
            theta = atan2(q(ax3,i), q(ax2,i)) + pi/2d0 ! plus 90 degree, which means omega direction
            !#write(*,*) pv, deg(theta)
            p1(ax2) = pv * dcos(theta)
            p1(ax3) = pv * dsin(theta)
            call unRotMomentum(ITevInv, p1, p0, q(1:3,i), q0(1:3,i) )
            p(1:3,i) = p(1:3,i) + p0(1:3)
            !# write(*,'(3F12.6)') p1(1:3) 
        end do
    end do
    
end subroutine


!!! subroutine : rotate momentum to original orientation according to inverted matrix of rotation
!!! p1(3) - rotated momentum
!!! matInv(3,3) - inverted matrix of rotation
!!! q0(3) - unrotated coordinates, check the correctness of inverted rotation
!!! p0(3) - unrotated momentum
subroutine unRotMomentum(matInv, p1, p0, q1, q0)
implicit none
real(8), intent(in) :: matInv(3,3), p1(3), q1(3), q0(3)
real(8), intent(out) :: p0(3)
integer i
real(8) pBond(3) ! the end of momentum bond, the start is q(1)
real(8) p1l, p0l ! length of p0 and p1
real(8) qInv(3) ! inverted coordinate, check the inversion matrix

    pBond(1:3) = q1(1:3) + p1(1:3)
    call vecLen(p1(1:3), p1l)
    
    qInv(1:3) = matmul( q1(1:3), matInv(1:3,1:3) )
    !# write(*,'(3F12.6)') q1, q0, qInv
    if ( maxval( qInv(1:3) - q0(1:3) ) .gt. 1d-5 ) call logError('Rotation failed! ')
    pBond(1:3) = matmul( pBond(1:3), matInv(1:3,1:3) )
    p0(1:3) = pBond(1:3) - q0(1:3)
    call vecLen(p0(1:3), p0l)
    p0(1:3) = p0(1:3) / p0l * p1l
    !# write(*,'(3F12.6)') p1, p0
    !# write(*,'(3F12.6)') p0l, p1l, p0l / p1l

end subroutine


! subroutine : generate rotational momentum of linear molecule
! temp - temperature in au
! RI - rotational inertia
subroutine linearRot(Natoms, coord, eleMass, temp, RI, p)
use constant
implicit none
integer, intent(in) :: Natoms
real(8), intent(in) :: coord(3,Natoms), eleMass(Natoms), temp, RI
real(8), intent(out) :: p(3,Natoms)
real(8) xi, Erot, J, omega
! move to origin
real(8) :: IT(3,3), ITevalue(3), ITev(3,3), ITevInv(3,3)
real(8) COM(3), q(3,Natoms) ! q is alternative of coord
integer i, rotAxis, anotherRotAxis(2)
real(8), external :: randab
real(8) theta ! random direction of p(vertical), in rad
real(8) pVertical ! momentum vertical, temporary for Natoms loop
real(8) pBond(3,2) ! shift momentum to the atom, rotate the momentum bond ! xyz, start->end
real(8) pTransLen ! scale transformed momentum to pVertical
real(8) linRotTemp

    call random_number(xi)
    Erot = - temp * dlog( 1 - xi ) ! kb = 1d0
    J = (-1d0 + dsqrt(1 + 8d0 * Erot * RI) ) * 0.5d0
    omega = dsqrt( 2d0 * Erot / RI )
    write(*,"(' Erot (kcal/mol) =', F12.6, ', J =', F12.6)") Erot*au2kcm, J
    if(ibug) write(*,"(' I (au·Ang^2) =', E9.2, ', omega (rad/fs) =', E9.2)") RI*b2a*b2a, omega/au2fs
    
    ! move to origin
    call getCOM(Natoms, coord, eleMass, COM)
    do i = 1, 3
        q(i,1:Natoms) = coord(i,1:Natoms) - COM(i)
    end do
    !#if(ibug) write(*,*) 'Move to origin: '
    !#if(ibug) call logXYZ(Natoms, q, 0)

    ! find the rotational axis
    call computeInertiaTensor(Natoms, coord, eleMass, IT)
    call diagMat(3, IT, ITevalue, ITev)
    call pinv(3, ITev, ITevInv)
    !#write(*,*); write(*,'(3F12.6)') ITev 
    !#write(*,*); write(*,'(3F12.6)') ITevInv
    !#write(*,*); write(*,'(3F12.6)') matmul( ITev, ITevInv )

    ! rotate to principal axis
    !#call logXYZ(Natoms, coord, 0)
    do i = 1, Natoms
        q(1:3,i) = matmul( q(1:3,i), ITev(1:3, 1:3) )
    end do
    !q = matmul(ITev, q)
    if(ibug) write(*,*) 'Rotate to principal axis: '
    if(ibug) call logXYZ(Natoms, q, 0)
    !!! Note : the bond lengths (in `q` and `coord`) will not change after transformation. 
    
    ! check linear
    if ( minval(ITevalue) .gt. 1D-5) then
        call logXYZ(Natoms, coord, 0)
        call logError('This molecule is not linear. Check the coordinates above! ')
    end if
    ! find the axis
    do i = 1, 3
        if( ITevalue(i) .lt. 1D-5 ) then
            if(ibug) write(*,'(A,I2)') ' Orientation (xyz->123): ', i
            rotAxis = i
            exit
        end if
    end do
    ! 
    call getAnotherTwoAxis(rotAxis, anotherRotAxis)
    if (ibug) write(*,'(A,2I2)') ' Another two axis : ', anotherRotAxis
    
    ! generate a random direction
    theta = randab(0, 2d0*pi)
    if (ibug) write(*,"(' Random direction of momentum (degree) : ', F8.2)"), theta / pi * 180d0
    p = 0d0
    do i = 1, Natoms
        pVertical = eleMass(i) * coord(rotAxis,i) * omega
        if(ibug) write(*,'(A,F12.6)') ' Momentum vertical : ', pVertical
        p(anotherRotAxis(1),i) = pvertical * dcos(theta)
        p(anotherRotAxis(2),i) = pvertical * dsin(theta)
    end do
    !call logXYZ(Natoms, p, 0)
    !! We can see that the momentums in a diatom molucule are same but with opposite direction. 
    
    !#write(*,'(3F12.6)') p
    !#call getTemp(Natoms, eleMass, p, 3, linRotTemp)
    !#write(*,*) linRotTemp * au2K
    
    !!! rotate the momentum to before (p*ITevInv)
    !do i = 1, Natoms
    !    p(1:3,i) = matmul( p(1:3,i), ITevInv(1:3,1:3) )
    !end do

    !!! bug : the temperature doesn't match. ! 2020-08-27 16:38:38 Wenbin, FAN @ SHU
    !!! bug fix note : the inversion matrix cound preserve the 2-norm (length). Then we scale the new momentum to before. ! 2020-08-28 15:25:33 Wenbin, FAN @ SHU
    !!! appendix MMA code : 
    !!!!! X = RandomReal[{-10, 10}, {8, 3}];
    !!!!! T = Eigenvectors[RandomReal[{0, 1}, {3, 3}]] // Chop;
    !!!!! TT = Inverse[T];
    !!!!! X.T.TT - X // Chop
    !!!!! A = RandomReal[{0, 1}, {2, 3}];
    !!!!! A // Norm
    !!!!! A.TT // Norm
    
    !!! transform momentum in atom
    !!call normEigenVec(3, ITevInv)
    !#call normEigenVec(3, reshape((/1d0,2d0,3d0,4d0,5d0,6d0,7d0,8d0,9d0/),(/3,3/)))

    !!#write(*,*) '4'
    !!#write(*,*)'coord rotated'; 
    !!#do i = 1, Natoms
    !!#write(*,'(A, 3F12.6)') 'O', q(:,i)
    !!#end do
    !!#!write(*,*)'p original'; 
    !!#do i = 1, Natoms
    !!#write(*,'(A, 3F12.6)') 'H', p(1:3,i)+q(1:3,i)
    !!#end do
    !!#write(*,*) '4'
    !!#write(*,*)'coord original';
    !!#do i = 1, Natoms
    !!#write(*,'(A, 3F12.6)') 'O', coord(:,i)
    !!#end do
    do i = 1, Natoms
        pBond(1:3,1) = q(1:3,i)
        pBond(1:3,2) = q(1:3,i) + p(1:3,i)
        pBond(1:3,1) = matmul( pBond(1:3,1), ITevInv(1:3,1:3) )
        pBond(1:3,2) = matmul( pBond(1:3,2), ITevInv(1:3,1:3) )
        p(1:3,i) = pBond(1:3,2) - pBond(1:3,1)
        !#write(*,'(A, 3F12.6)') 'p(1:3,i)', p(1:3,i)
        call vecLen(p(1:3,i), pTransLen)
        p(1:3,i) = p(1:3,i) / pTransLen * abs(pVertical)
        !#write(*,*) 'pTransLen', pTransLen, pVertical
        !!#write(*,'(A, 6F12.6)') 'H', pBond(1:3,2), p(1:3,i) + coord(1:3,i)
    end do
    !#write(*,'(3F12.6)') p
    !#call getTemp(Natoms, eleMass, p, 3, linRotTemp)
    !#write(*,*) linRotTemp * au2K

end subroutine


subroutine ABrandThermal
use para
implicit none
real(inum) :: rotTempA, rotTempB

    write(*,*); write(*,*) 'Generate initial momentum...'
    call randThermal(NA, eleMassA, tempAU, pA)
    call randThermal(NB, eleMassB, tempAU, pB)

    call getTemp(NA, eleMassA, pA, DOFA, rotTempA)
    call getTemp(NB, eleMassB, pB, DOFB, rotTempB)
    write(*,'(A,2F8.2)') ' Temperature of A, B (K): ', rotTempA * au2K, rotTempB * au2K

end subroutine


subroutine randThermal(Natoms, eleMass, temp, p)
implicit none
integer, intent(in) :: Natoms
real(8), intent(in) :: eleMass(Natoms), temp
real(8), intent(out) :: p(3, Natoms)
integer i, j
real(8), external :: randNorm
real(8) sqrtTM

    do i = 1, Natoms
        sqrtTM = dsqrt( temp / eleMass(i) ) ! kB = 1d0
    do j = 1, 3 ! xyz
        p(j,i) = randNorm(0, sqrtTM)
    end do
    end do

end subroutine


subroutine ABrandRot
use para
implicit none
real(8) mat(3,3)
integer i

    write(*,*) 'Randomly rotate molecules...'
    call rotMat(mat)
    do i = 1, NA
        coordA(1:3,i) = matmul( coordA(1:3,i), mat(1:3,1:3) )
    end do
    call rotMat(mat)
    do i = 1, NB
        coordB(1:3,i) = matmul( coordB(1:3,i), mat(1:3,1:3) )
    end do
    
    if (ibug) then
        write(*,*) ' Molecule A : '
        call logXYZ(NA, coordA, eleNameA)
        write(*,*) ' Molecule B : '
        call logXYZ(NB, coordB, eleNameB)
    end if

end subroutine


subroutine BrandCoord
use para
implicit none
real(inum) :: xi ! random b0, impact factor
real(inum) :: b1, b2

    if (b0_opt) then ! if b0 provided
        b1 = b0
        b2 = dsqrt(r0*r0 - b1*b1)
    else
        call random_number(xi)
        b1 = r0 * dsqrt(xi)
        b2 = r0 * dsqrt(1 - xi)
    end if
    write(*,'(A,2F12.6)') ' Impact factor (Ang) b1, b2: ', b1*b2a, b2*b2a
    !#call logXYZ(NB, coordB, eleNameB)
    call transCoord(NB, coordB, (/0d0, b1, -b2/))
    !#call logXYZ(NB, coordB, eleNameB)

end subroutine


! subroutine : generate a random rotation matrix using quaternion method
! ref : F. J. Vesely, J. Comput. Phys. 47, 291 (1982). 
!       A. R. Leach, Molecular Modelling: Principles and Applications (Pearson Education, 2001). 420-422. 
subroutine rotMat(mat)
use constant
implicit none
real(8), intent(out) :: mat(3,3)
real(8), external :: randab
real(8) :: x1,x2,x3,x4, q0,q1,q2,q3, s1,s2, ratio
real(8) evalue(3), ev(3,3)
real(8) theta, phi, varphi ! three Euler angles
real(8), external :: deg ! rad to degree

    s1 = 1d5; s2 = 1d5 ! a big value
    do while (s1 .gt. 1d0)
        x1 = randab(-1d0,1d0)
        x2 = randab(-1d0,1d0)
        s1 = x1*x1 + x2*x2
    end do
    
    do while (s2 .gt. 1d0)
        x3 = randab(-1d0,1d0)
        x4 = randab(-1d0,1d0)
        s2 = x3*x3 + x4*x4
    end do
    
    q0 = x1; q1 = x2
    ratio = dsqrt( (1-s1)/s2 )
    q2 = x3 * ratio
    q3 = x4 * ratio
    
    mat(1,1) = q0*q0 + q1*q1 - q2*q2 - q3*q3
    mat(2,1) = 2d0 * (q1*q2 + q0*q3)
    mat(3,1) = 2d0 * (q1*q3 - q0*q2)
    mat(1,2) = 2d0 * (q1*q2 - q0*q3)
    mat(2,2) = q0*q0 - q1*q1 + q2*q2 - q3*q3
    mat(3,2) = 2d0 * (q2*q3 + q0*q1)
    mat(1,3) = 2d0 * (q1*q3 + q0*q2)
    mat(2,3) = 2d0 * (q2*q3 - q0*q1)
    mat(3,3) = q0*q0 - q1*q1 - q2*q2 + q3*q3
    !# write(*,'(3F12.6)') mat
    
    ! compute Euler angles ! computationally expensive because of no use of repetitive variable
    if (ibug) then
        theta = 2d0*atan(dsqrt( (q1*q1+q2*q2)/(q0*q0+q3*q3) ))
        phi = atan2(q2,q1) + atan2(q3,q0)
        varphi = -atan2(q2,q1) + atan2(q3,q0)
        write(*,"(' Euler angles (degree): ', 3F9.2)") deg(theta), deg(phi), deg(varphi)
    end if

end subroutine


subroutine loadInitMol
use para
implicit none

    coordA(1:3,1:NA) = coordAInit(1:3,1:NA)
    coordB(1:3,1:NB) = coordBInit(1:3,1:NB)
    pA(1:3,1:NA) = 0d0
    pB(1:3,1:NB) = 0d0

end subroutine


subroutine unloadInitMol
use para
implicit none

    deallocate(freq, fVec)

end subroutine