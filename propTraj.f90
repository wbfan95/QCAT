module prop
use para, only: temp, Natoms, eleMass, eleName, &
                traj, &
                DOF, coord, p, tstep, &
                thermostat, Nstep, Nprint, Nxyz
use constant
implicit none

    ! Nose-Hoover
    real(inum) :: NHM(2), NHv(2), NHx(2)
    
    ! properties of trajectory
    real(inum) :: Ek ! kinetic energy
    
    ! prop
    integer :: step
    integer :: istop
    real(inum) :: tstep2
    real(inum) :: tstep4
    real(inum) :: tstep8

end module


subroutine propTraj
use prop
implicit none

    write(*,*); write(*,*) 'Start propagating trajectory...'
    write(*,*) ' Step          Epot          Ek        Etot       Temp.'
    call propTrajInit
    do while (istop .gt. 0)
        if (thermostat) call chain
        call velocityVerlet
        if (thermostat) call chain
        if (mod(step,Nprint) .eq. 0) call propPrint
        if (mod(step,Nxyz) .eq. 0) call xyzPrint
        call stopCheck
    end do
    
    close(999)

end subroutine


subroutine xyzPrint
use prop
implicit none
integer i

    write(999, '(I4)') Natoms
    write(999, '(I6)') step
    do i = 1, Natoms
        write(999, '(A2,3F12.6)') eleName(i), coord(1:3,i)*b2a
    end do

end subroutine


subroutine propPrint
use prop
implicit none
real(inum) ene, tempNow

    call pot(Natoms, coord, ene)
    call getTemp(Natoms, eleMass, p, DOF, tempNow)
    
    write(*,'(I8,4F12.2)') step, &
        ene*au2kcm, Ek*au2kcm, (ene+Ek)*au2kcm, &
        tempNow*au2K

end subroutine


subroutine stopCheck
use prop
implicit none
integer :: i
real(inum) :: maxDist(3)

    step = step + 1

    ! max distance in Cartesian coodinates
    do i = 1, 3 ! xyz
        maxDist(i) = maxval(coord(i,:))
    end do
    if ( maxval(maxDist(1:3)) .gt. 10d0 ) then
        istop = -1
        write(*,*) 'Prop end! A, B are far away from each other! '
        return
    end if
    
    ! max step
    if (step .gt. Nstep) then
        istop = -1
        write(*,*) 'Prop end! Max step reached! '
        return
    end if

end subroutine


subroutine velocityVerlet
use prop
implicit none
integer i
real(inum) grad(3, Natoms), ene

    Ek = 0d0
    
    do i = 1, Natoms
        coord(1:3,i) = coord(1:3,i) + p(1:3,i) / eleMass(i) * tstep2
    end do
    
    call calcGrad(Natoms, coord, grad)
    !do i = 1, Natoms
    !    write(*,'(A,3F12.6)') eleName(i), grad(:,i)*b2a
    !end do
    !call pot(Natoms, coord, ene)
    !write(*,*) 'Energy: ', ene*au2eV
    if (maxval(grad) - minval(grad) .lt. 1d-12) then
        istop = -100
        write(*,*) 'Prop end! Zero gradient was found! '
    end if
    
    do i = 1, Natoms
        p(1:3,i) = p(1:3,i) - grad(1:3,i)*tstep
        coord(1:3,i) = coord(1:3,i) + p(1:3,i) / eleMass(i) * tstep2
        Ek = Ek + dot_product( p(1:3,i), p(1:3,i) ) / 2d0 / eleMass(i)
    end do

end subroutine


! subroutine : Nose-Hoover chain
! source : Frenkel & Smit, Understanding Molecular Simulation, Page 536 - 544. 
!!! Note : This algorithm amazed me but idk how does it work. :-(
subroutine chain
use prop
implicit none
integer :: i
real(inum) G1, G2, s!, tempNow

    G2 = (NHM(1) * NHv(1)*NHv(1) - temp) / NHM(2)
    NHv(2) = NHv(2) + G2*tstep4
    NHv(1) = NHv(1) * dexp(-NHV(2)*tstep8)
    G1 = (2d0*Ek - DOF*temp) / NHM(1)
    
    NHv(1) = NHv(1) + G1*tstep4
    NHv(1) = NHv(1) * dexp(-NHv(2)*tstep8)
    NHx(1) = NHx(1) + NHv(1)*tstep2
    NHx(2) = NHx(2) + NHv(2)*tstep2
    
    s = dexp(-NHv(1)*tstep2)
    !#write(*,*) 's, ', s
    p(1:3,1:Natoms) = p(1:3,1:Natoms) * s
    Ek = Ek*s*s
    
    NHv(1) = NHv(1) * dexp(-NHv(2)*tstep8)
    G1 = (2d0*Ek - DOF*temp) / NHM(1)
    NHv(1) = NHv(1) + G1*tstep4
    NHv(1) = NHv(1) * dexp(-NHv(2)*tstep8)
    G2 = (NHM(1) * NHv(1)*NHv(1) - temp) / NHM(2)
    NHv(2) = NHv(2) + G2*tstep4
    
    !call getTemp(Natoms, eleMass, p, DOF, tempNow)
    !write(*,*) 'Temp: ', tempNow*au2K

end subroutine


subroutine propTrajInit
use prop
implicit none
integer i
real(inum) sysfreq
character(len=4) ntrajChar

    ! propagation
    tstep2 = tstep / 2d0
    tstep4 = tstep / 4d0
    tstep8 = tstep / 8d0
    
    ! Nose-Hoover initial condition
    sysfreq = 1d0
    NHM(1:2) = 2d0*DOF*temp*sysfreq*sysfreq
    NHv(1:2) = 0d0
    NHx(1:2) = 0d0
    
    ! get Ek
    Ek = 0d0
    do i = 1, Natoms
        Ek = Ek + dot_product( p(1:3,i), p(1:3,i) ) / 2d0 / eleMass(i)
    end do
    !#write(9876,*) Ek*au2kcm
    
    step = 1
    istop = 100
    
    write(ntrajChar, '(I0.4)') traj
    open(999, file='traj//traj_'//ntrajChar//'.xyz', status='replace', action='write')

end subroutine

