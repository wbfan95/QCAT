module constant
implicit none

    ! precision
    integer, parameter :: inum = 8 ! double precision
    integer, parameter :: icha = 1024 ! default length of character
    
    ! constant
    real(inum), parameter :: pi = dacos(-1d0)
    
    ! physical constant
    real(inum), parameter :: kb = 1.380649d-23 ! au/K, kB = 1.380649E-23 J/K
    real(inum), parameter :: clight = 2.99792458d8 ! c, light speed, unit : m/s
    real(inum), parameter :: hbar = 1.0545718176461565e-34 ! Plank constant / 2 / Pi, J·s
    
    ! physical unit convert
    real(inum), parameter :: au2eV = 27.211386245988d0 ! eg : 1 a.u. = au2eV eV = 27.21~ eV
    real(inum), parameter :: au2kcm = 627.509474063056d0 ! kcal/mol
    real(inum), parameter :: au2K = 315775.02480407d0
    real(inum), parameter :: au2fs = 0.024188843265857d0 ! time not period
    real(inum), parameter :: au2s = 2.4188843265857d-17
    real(inum), parameter :: au2J = 4.3597447222071d-18
    real(inum), parameter :: au2cm1 = 219474.6313632d0 ! cm^{-1}
    real(inum), parameter :: au2kg = 1.6605390666e-27
    real(inum), parameter :: au2m = 5.29177210903e-11 ! Bohr to meter
    real(inum), parameter :: au2Hz = 6.579683920502e+15 ! Hertz = [second]^{-1}
    real(inum), parameter :: a2b = 1.88972612462577d0 ! Ang to Bohr
    real(inum), parameter :: b2a = 0.529177210903d0 ! Bohr to Ang
    
    ! system precision
    real(inum), parameter :: linearTh = 1d-5
    real(inum), parameter :: delta = 1d-5 ! finite displacement, Bohr
    logical ibug ! debug option

end module

module fileIO
use constant, only : icha
implicit none

    !! log file
    !integer, parameter :: logFileID = 1
    !integer, parameter :: flog = 1 ! a abbreviation of `logFileID`, means Fortran log (flog)
    !character(icha), parameter :: logFile = 'stdout.log' 

end module


module para
use constant
use fileIO
implicit none

    ! structure
    integer :: Natoms
    integer :: NA, NB ! number of atom A and B
    integer, allocatable :: order(:) ! record order to PES order
    real(inum), allocatable :: coord(:,:)
    character(2), allocatable :: eleName(:)
    integer, allocatable :: eleNum(:)
    real(inum), allocatable :: eleMass(:)
    real(inum) :: massA, massB, muAB ! muAB is reduced mass of A and B $\mu = \frac{m_Am_B}{m_A+m_B}$
    ! A, B seperately
    real(inum), allocatable :: coordA(:,:), coordB(:,:)
    real(inum), allocatable :: eleMassA(:), eleMassB(:)
    character(2), allocatable :: eleNameA(:), eleNameB(:)
    real(inum) :: ITA(3,3), ITB(3,3) ! inertia tensor
    real(inum) :: RIA, RIB ! rotational inertia
    logical ilinearA, ilinearB ! 1: linear, 0: non-linear
    integer DOF, DOFA, DOFB ! degree of freedom. single atom: 1, diatom: 2, multiple atom: 3N-3
    real(inum), allocatable :: coordAInit(:,:), coordBInit(:,:) ! save initialized molecule here
    
    ! trajectory
    real(inum), allocatable :: p(:,:) ! momentum
    real(inum), allocatable :: pA(:,:), pB(:,:) ! momentum
    integer :: nMode ! number of valid vibrational modes
    real(inum), allocatable :: freq(:), fVec(:,:)
    
    ! system
    character(3) :: ensemble
    real(inum) :: temp, tempAU
    real(inum) :: tstep, tstepAU
    logical :: b0_opt
    real(inum) :: b0, r0 ! r0 = 2.2*(diameter A + diameter B) ! ANT 2019 manual P71
    
    ! output
    integer :: nprint, ntraj


end module


! calculate eigenvalue and eigenvector
! N : size of matrix
! matIn : matrix(N,N)
! evReal : real eigenvalue(N)
! ev : eigenvector(N,N)
subroutine diagMat(N, matIn, evReal, ev)
use lapack95, only: geev
integer, intent(in) :: N
real(8), intent(in) :: matIn(N,N)
real(8) :: evReal(N), evImag(N)
real(8) :: evl(N,N), evr(N,N)
real(8) :: work(N*100)
real(8), intent(out) :: ev(N,N)
integer :: info, i

    call dgeev('V','V',N,matIn,N,evReal,evImag,evl,N,evr,N,work,N*100,info)
    ev(1:N,1:N) = evl(1:N,1:N)
    if (info /= 0) then
        write(*,*) matIn
        write(*,*) 'Diagonalization failed! Check the matrix above. '
        stop
    end if
    !#write(*,*) 'eigenvalue imaginary: '
    !#write(*,*) evImag

end subroutine


subroutine pinv(N, A, pinvA)!, info)
use lapack95, only: gesvd
implicit none
integer, intent(in) :: N
real(8), intent(in) :: A(N,N)
real(8), intent(out) :: pinvA(N,N)
real(8) U(N,N), S(N), VT(N,N)!, WW(3)
real(8) S1(N,N), V(N,N), VS1(N,N)!, SM(N,N) ! S matrix -> SM
real(8) work(N*100)!, lwork
!integer, intent(out) :: info
real(8) saveA(N,N) ! preserve the value of A, because A will be overwritten ( WTF? )
integer info
integer i

U = 0d0
VT = 0d0
saveA(1:N,1:N) = A(1:N,1:N)
call dgesvd('A', 'A', N, N, saveA, N, S, U, N, VT, N, work, N*100, info)
if (info /= 0) then
    write(*,*) '[ERROR] Matrix inversion error! '
    stop
end if
! pseudo inversion for S
S1 = 0d0
!SM = 0d0
do i = 1, N
if (S(i) .gt. 1d-10) S1(i,i) = 1/S(i)
!SM(i,i) = S(i)
end do

! A-1 = V * S-1 * U^T
VS1 = matmul(transpose(VT), S1)
pinvA = matmul(VS1, transpose(U))

!# write(*,*) 'A: '
!# write(*,'(3F12.6)') A
!# write(*,*) 'Inv A: '
!# write(*,'(3F12.6)') pinvA
!# write(*,*) 'U: '
!# write(*,'(3F12.6)') U
!# write(*,*) 'VT: '
!# write(*,'(3F12.6)') VT
!# write(*,*) 'SM: '
!# write(*,'(3F12.6)') SM
!# write(*,*) 'U Sigma VT: '
!# write(*,'(3F12.6)') matmul( matmul( transpose(VT), SM ), transpose(U) )
!# write(*,*)

end subroutine

! subroutine : cross product
subroutine cross_product(a, b, cross)
real(8), intent(in) :: a(3), b(3)
real(8), intent(out) :: cross(3)

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
  
end subroutine

!!! subroutine : normalize all column vector
!!! it should be checked in the situation that summation is near zero. 
subroutine normEigenVec(N, mat)
implicit none
integer, intent(in) :: N
real(8), intent(inout) :: mat(N,N)
integer i, j
real(8) colSum

    do i = 1, N
        colSum = 0d0
        do j = 1, N
            colSum = colSum + mat(i, j) * mat(i, j)
        end do
        if (colSum .ne. 0d0) then 
        mat(i, 1:N) = mat(i, 1:N) / dsqrt(colSum)
        end if
    end do
    
    !!do i = 1, N
    !!    colSum = 0d0
    !!    do j = 1, N
    !!        colSum = colSum + mat(j, i) * mat(j, i)
    !!    end do
    !!    if (colSum .ne. 0d0) then 
    !!    mat(1:N, i) = mat(1:N, i) / dsqrt(colSum)
    !!    end if
    !!end do

end subroutine

! random from lower bound to upper bound
function randab(lowBound, upBound)
implicit none
real(8) randab
real(8) :: lowBound, upBound

    ! swap the upper and lower bound if the order is wrong. 
    if (lowBound .gt. upBound) then ! if l=100, u=1
    lowBound = lowBound + upBound   ! l = l + u = 101
    upBound  = lowBound - upBound   ! u = l - u = 101 - 1 = 100
    lowBound = lowBound - upBound   ! l = l - u = 101 - 100 = 1 ! Bingo!
    end if

    call random_number(randab)
    randab = randab * (upBound - lowBound) + lowBound

end function


function randNorm(opt_mean, opt_std)
implicit none
real(8) randNorm
real(8), intent(in), optional :: opt_mean
real(8), intent(in), optional :: opt_std
real(8) zeta(2), mean, std
real(8), parameter :: pi = dacos(-1d0)


    if (present(opt_mean) ) then
        mean = opt_mean
    else
        mean = 0d0
    end if
    
    if (present(opt_std) ) then
    std = opt_std
    else
    !!std = mean * 0.125d0 ! \sigma = \bar x / 4 ! expect that rand(x) satisfy 0.75 \bar x < rand(x) < 1.25 \bar x
    std = 1d0
    end if

    call random_number(zeta)
    randNorm = dsqrt(- 2d0 * log(zeta(1)) ) * dcos( 2d0*pi*zeta(2) )
    randNorm = mean + std * randNorm

end function

subroutine getTemp(Natoms, eleMass, p, DOF, temp)
implicit none
integer, intent(in) :: Natoms, DOF
real(8), intent(in) :: eleMass(Natoms)
real(8), intent(in) :: p(3, Natoms)
real(8), intent(out) :: temp
integer i,j

    temp = 0d0
    do i = 1, Natoms
    do j = 1, 3
        temp = temp + p(j,i) * p(j,i) / eleMass(i)
    end do
    end do
    temp = temp / DOF
    !#write(*,*) 'Temperature : ', temp
    
    if (temp .lt. 0d0) then
        write(*,*) 'momentum: '
        write(*,'(3F12.6)') p
        write(*,*) 'mass: '
        write(*,'(6F12.6)') eleMass
        write(*,*) 'DOF: '
        write(*,'(I)') DOF
        call logError('Temperature is less than zero! ')
    end if

end subroutine


! subroutine : get anthor two axis, input x (1) and output yz (2,3)...
subroutine getAnotherTwoAxis(i, o)
implicit none
integer, intent(in) :: i
integer, intent(out) :: o(2)
integer j, c ! count

    c = 1
    do j = 1, 3
        if (j .eq. i) then
            continue
        else
            o(c) = j
            c = c + 1
            if (c .eq. 3) exit
        end if
    end do

end subroutine

!!! subroutine : calculate the length of vector
subroutine vecLen(vec, length)
implicit none
real(8), intent(in) :: vec(3)
real(8), intent(out) :: length
integer :: i

    length = 0d0
    do i = 1, 3
        length = length + vec(i) * vec(i)
    end do
    length = dsqrt(length)

end subroutine

! subroutine : shift coordinate
subroutine transCoord(Natoms, coord, vec)
implicit none
integer, intent(in) :: Natoms
real(8), intent(inout) :: coord(3, Natoms)
real(8), intent(in) :: vec(3)
integer :: i

    do i = 1, 3 ! x,y,z
        coord(i,1:Natoms) = coord(i,1:Natoms) + vec(i)
    end do

end subroutine


!!! function : convert rad to degree
function deg(rad)
use constant
implicit none
real(8), intent(in) :: rad
real(8) deg

    deg = rad / pi * 180d0

end function

!!! subroutine : calculate Euler angles from rotational matrix
subroutine getEulerAngle(mat)
use constant
implicit none
real(8), intent(in) :: mat(3,3)
real(8) :: alpha, beta, gamma
real(8), external :: deg

    beta = asin(-mat(3,1))
    alpha = asin( mat(2,1) / cos(beta) )
    gamma = asin( mat(3,2) / cos(beta) )
    write(*,'(A,3F9.2)') ' Rotational Euler angle (degree) : ', deg(alpha), deg(beta), deg(gamma)

end subroutine
