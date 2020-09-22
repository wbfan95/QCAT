program main
implicit none

call initProgram
call initMolecule
call batchTraj

end program


subroutine initProgram
implicit none

call logHeader
call readInput

call random_seed()
call initPES

end subroutine


subroutine batchTraj
use para
implicit none

    write(*,*); write(*,*) 'Batching trajectories...'
    do traj = 1, ntraj
        write(*,'(A,I8)') ' Trajectory ', traj
        call initTraj ! initialize trajectory
        call propTraj ! propagate trajectory
        
        call loghline
    end do

end subroutine






