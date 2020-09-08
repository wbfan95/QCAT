! initialize the PES
subroutine initPES
call pes_init
end subroutine

! subroutine : get energy from potential energy surface
! coord - coordinate in Bohr
! energy - Hartree
subroutine pot(Natoms, coord, energy)
use constant
implicit none
integer, intent(in) :: Natoms
real(inum), intent(in) :: coord(3, Natoms)
real(inum), intent(out) :: energy
real(inum) :: v1,v2,v3

    call ohch4pipNN(coord*b2a, energy, v1,v2,v3)
    energy = energy / au2eV

end subroutine

subroutine pot1d(Natoms, coord1D, energy)
use constant
implicit none
integer, intent(in) :: Natoms
real(inum), intent(in) :: coord1D(3*Natoms)
real(inum), intent(out) :: energy

    call pot(Natoms, reshape(coord1D, (/3,Natoms/)), energy)

end subroutine


! subroutine : calculate gradient
! grad - Hartree/Bohr = a.u.
subroutine calcGrad(Natoms, coord, grad)
use constant
implicit none
integer, intent(in) :: Natoms
real(inum), intent(in) :: coord(3, Natoms)
real(inum), intent(out) :: grad(3, Natoms)
real(inum) :: q(3,Natoms), ef, eb ! energy forward, energy backward
integer :: i, j

    q(1:3,1:Natoms) = coord(1:3,1:Natoms)
    do i = 1, Natoms
        do j = 1, 3
            q(j,i) = q(j,i) + delta
            call pot(Natoms, q, ef)
            q(j,i) = q(j,i) - delta - delta
            call pot(Natoms, q, eb)
            q(j,i) = q(j,i) + delta
            grad(j,i) = (ef - eb) / delta / 2d0
        end do
    end do 

end subroutine


subroutine freqAna(Natoms, coord, eleMass, freq, evec)
use constant
implicit none
integer, intent(in) :: Natoms
real(inum), intent(in) :: coord(3,Natoms), eleMass(Natoms)
real(inum), intent(out) :: evec(3*Natoms,3*Natoms), freq(3*Natoms)
real(inum) :: hessian(3*Natoms,3*Natoms), hessianMW(3*Natoms,3*Natoms)
real(inum) :: eval(3*Natoms)
integer :: i

    call calcHessian(Natoms, coord, hessian)
    call calcMWHessian(Natoms, eleMass, hessian, hessianMW)
    call diagMat(Natoms*3, hessianMW, eval, evec)
    !#call syev(hessian, eval)
    
    eval(1:3*Natoms) = eval(1:3*Natoms) * au2J/au2m**2/au2kg
    do i = 1, Natoms*3
        if ( eval(i) .lt. 0d0 ) then
            freq(i) = - dsqrt( -eval(i) ) / pi / 2d0
        else
            freq(i) = dsqrt(eval(i)) / pi / 2d0
        end if
    end do
    freq = freq / clight / 1d2 / au2cm1 ! to cm-1, to a.u.
    !#write(*, '(8F9.2)') freq * au2cm1

end subroutine

subroutine calcMWHessian(Natoms, eleMass, hessian, hessianMW)
implicit none
integer, intent(in) :: Natoms
real(8), intent(in) :: eleMass(Natoms), hessian(3*Natoms,3*Natoms)
real(8), intent(out) :: hessianMW(3*Natoms,3*Natoms)
integer :: i!,j,m,n

    hessianMW = 0d0
    !do i = 1, Natoms
    !    do j = 1, Natoms
    !        hessianMW( (i-1)*3+1:i*3, (j-1)*3+1:j*3 ) = &
    !        hessian( (i-1)*3+1:i*3, (j-1)*3+1:j*3 ) / &
    !        dsqrt( eleMass(i) * eleMass(j) )
    !    end do
    !end do
    
    do i = 1, Natoms
        hessianMW(1:Natoms*3,(i-1)*3+1:i*3) = hessian(1:Natoms*3,(i-1)*3+1:i*3) / dsqrt(eleMass(i))
    end do
    do i = 1, Natoms
        hessianMW((i-1)*3+1:i*3,1:Natoms*3) = hessianMW((i-1)*3+1:i*3,1:Natoms*3) / dsqrt(eleMass(i))
    end do
    
    !#write(*,*) 'New method MW: (*100) '
    !#write(*,'(21F6.1)') hessianMW*100

end subroutine


subroutine calcHessian(Natoms, coord, hessian)
use constant
implicit none
integer, intent(in) :: Natoms
real(inum), intent(in) :: coord(3, Natoms)
real(inum), intent(out) :: hessian(3*Natoms, 3*Natoms)
integer :: i,j,k, N3
real(inum) :: q(3*Natoms), secDer

    N3 = Natoms * 3
    q(1:N3) = reshape(coord(1:3, 1:Natoms), (/N3/))
    hessian(1:N3,1:N3) = 0d0
    
    ! upper part include diagonal elements
    do i = 1, N3
        do j = i, N3
            call calcSecDer(Natoms, q, i, j, secDer)
            hessian(i,j) = secDer
            !#write(*,*) i,j
        end do
    end do
    
    ! lower part
    do i = 2, N3
        do j = 1, i - 1
            hessian(i,j) = hessian(j,i)
        end do
    end do


end subroutine

! subroutine : calculate second derivate of potential energy
subroutine calcSecDer(Natoms, coord1D, i, j, secDer)
use constant
implicit none
integer, intent(in) :: Natoms, i, j
real(inum), intent(in) :: coord1D(Natoms*3)
real(inum), intent(out) :: secDer
integer :: N3
real(inum) q(Natoms*3), epp, epm, emp, emm ! e_plus_plus, e_plus_minus, ...

    !#write(*,'(3F12.6)') coord1D
    N3 = Natoms * 3
    q(1:N3) = coord1D(1:N3)
    ! +x, +y
    q(i) = q(i) + delta
    q(j) = q(j) + delta
    call pot1D(Natoms, q, epp)
    q(i) = coord1D(i)
    q(j) = coord1D(j)
    
    ! +x, -y
    q(i) = q(i) + delta
    q(j) = q(j) - delta
    call pot1D(Natoms, q, epm)
    q(i) = coord1D(i)
    q(j) = coord1D(j)
    
    ! -x, +y
    q(i) = q(i) - delta
    q(j) = q(j) + delta
    call pot1D(Natoms, q, emp)
    q(i) = coord1D(i)
    q(j) = coord1D(j)
    
    ! -x, -y
    q(i) = q(i) - delta
    q(j) = q(j) - delta
    call pot1D(Natoms, q, emm)
    q(i) = coord1D(i)
    q(j) = coord1D(j)
    
    secDer = (epp - epm - emp + emm) / (4d0 * delta * delta)
    !#write(*,'(3F12.6)') q
    !#write(*,'(5E12.4)') epp, epm, emp, emm, secDer

end subroutine