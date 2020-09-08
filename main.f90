program main
implicit none

call initProgram
call initMolecule
call batchTraj

!#call test

end program

!subroutine test
!use para
!implicit none
!real(8) :: q(3,7), hessian(21,21), evec(21,21), freq(21)
!character(4) :: eName(7)
!integer i
!
!    open(999, file='test.xyz', status='old', action='read')
!    read(999,*); read(999,*)
!    do i = 1,7
!        read(999,*) eName(i), q(1:3,i)
!        !write(*,*) q(1:3,i)
!    end do
!    q = q * a2b
!    call freqAna(7, q, eleMass, freq, evec)
!    write(*,'(6F9.2)') freq * au2cm1
!
!end subroutine


subroutine initProgram
!use para
implicit none

call logHeader
call readInput

call random_seed()
call initPES

end subroutine


subroutine batchTraj
use para
implicit none
integer :: traj

    write(*,*); write(*,*) 'Batching trajectories...'
    do traj = 1, ntraj
        write(*,'(A,I8)') ' Trajectory ', traj
        call initTraj
        
        call loghline
    end do

end subroutine






