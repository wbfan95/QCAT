! input file ID : 478 ! input -> ipt -> 478
! scratch file ID : 479 ! 478+1

subroutine readInput
use para
implicit none

    call readInputFile
    call readAtomCoord
    call readParameter

end subroutine



subroutine readInputFile
use constant, only : icha
implicit none
character(icha) inputFileName
logical inputFileStatus
integer inputFileReadStatus, inputFileLineCount
character(icha) uppercase, lineContent ! temporary content of a line

    call getarg(1, inputFileName)
    
    ! input file was not given
    if (inputFileName .eq. '') then
        call logError('No input file! Show me the input file as `./QCAT inputfile`! ')
    end if
    
    ! check existence of input file
    inquire(file=trim(inputFileName), exist=inputFileStatus)
    if (inputFileStatus) then
        open(478, file=trim(inputFileName), status='old', action='read') ! input -> ipt -> 478
        write(*,*) 'Input file : ', trim(inputFileName)
    else
        call logError("Your input file `" // trim(inputFileName) // "` doesn't exist! ")
    end if
    
    ! count the line, transfer contents to scratch file
    open(479, status='scratch', action='readwrite')
    inputFileLineCount = 0
    do while (.true.)
        read(478, '(A)', iostat=inputFileReadStatus) lineContent
        if (inputFileReadStatus/=0) exit
        inputFileLineCount = inputFileLineCount + 1
        lineContent = uppercase(lineContent) ! convert to uppercase
        write(479, *) trim(lineContent)
        !#write(*,*) trim(lineContent)
    end do
    ! check whether the input file is blank
    if (inputFileLineCount == 0) then
        call logError('Your input file is blank! ')
    else
        write(*,"(' Lines in input file : ', I4)") inputFileLineCount
        rewind(478) ! go back to the beginning of input file
    end if


end subroutine


subroutine readAtomCoord
use para
implicit none
integer i,j,k
character(2) atomType ! 'A' or 'B'
character(icha) orderLine ! the content of line recording the order


    write(*,*); write(*,*) 'Reading input file...'
    ! read Natoms
    read(478,*) NA  ! read number of atom A
    do i = 1, NA    ! 
        read(478,*) ! 
    end do          ! 
    read(478,*) NB  ! read number of atom B
    do i = 1, NB    ! 
        read(478,*) ! 
    end do          ! 
    Natoms = NA + NB
    
    ! read order of atoms
    allocate( order(Natoms) )
    read(478,'(A)') orderLine
    if (orderLine .ne. '') then
        read(orderLine, *) order(1:Natoms)
    else ! if the order line is blank, normal order will be generated
        order(1:Natoms) = (/ (j, j=1,Natoms) /)
    end if
    !# write(*,*) order(1:Natoms)
    rewind(478)
    
    ! read coordinates of atoms
    allocate( coord(3,Natoms), p(3,Natoms), eleName(Natoms) )
    read(478,*) ! NA
    do i = 1, NA
        read(478, *) eleName( order(i) ), ( coord(j,order(i)), j=1,3 )
    end do
    read(478,*) ! NB
    do i = NA + 1, Natoms
        read(478, *) eleName( order(i) ), ( coord(j,order(i)), j=1,3 )
    end do
    rewind(478)
    
    ! write atom coordinates
    coord(1:3,1:Natoms) = coord(1:3,1:Natoms) * a2b ! atomic unit used
    allocate( eleNum(Natoms), eleMass(Natoms) )
    call eleNameParse ! in the `initMol.f90`
    write(*,*) 'Atom coordinates :'
    write(*,*) '  Type  Ele.    x/Ang       y/Ang       z/Ang      mass/au'
    write(*,*) '  ----  ----  ---------   ---------   ---------   ---------'
    do i = 1, Natoms
        ! judge A or B
        do j = 1, NA
            if ( i .eq. order(j) ) then
                atomType = 'A'
                exit
            else
                atomType = 'B'
            end if
        end do
        ! log all coordinates: atom type, element name, coordinates
        write(*,'(4X,A2,4X,A2,3F12.6,F12.4)') &
            trim(atomType), trim(eleName(i)), coord(1:3,i)*b2a, eleMass(i)
    end do
    
    close(478) ! close unit 478, and all parameters will be read from unit 479 (upper case)

end subroutine


subroutine readParameter
use para
implicit none
character(icha) :: locOut
integer :: lresult

    ! debug
    call locPara(479, 'DEBUG', locOut, lresult)
    if (lresult .gt. 0) then; ibug = .true. ; write(*,*) '*** Debug mode! ***'
    else; ibug = .false. ; end if
    
    ! b0, impact factor
    call locPara(479, 'B0', locOut, lresult)
    if (lresult .gt. 0) then; read(locOut, *) b0; b0_opt = .true.
    else
        b0_opt = .false.
        write(*,*) '`B0` will be generated randomly! '
    end if
    
    ! ensemble
    call locPara(479, 'ENS', locOut, lresult)
    if (lresult .gt. 0) then; read(locOut, *) ensemble
    else
        ensemble = 'NVE'
        write(*,*) '`ENSEMBLE` was set to default NVE! '
    end if
    write(*,"(' Ensemble    : ', A3)") ensemble
    if ( ensemble .eq. 'NVE' ) then
        thermostat = .false.
    else
        thermostat = .true.
    end if
    
    ! print gap
    call locPara(479, 'NP', locOut, lresult)
    if (lresult .gt. 0) then; read(locOut, *) Nprint
    else
        Nprint = 100
        write(*,*) '`NPRINT` was set to default 100! '
    end if
    write(*,"(' N_print     : ', I8)") Nprint
    
    ! number of steps
    call locPara(479, 'NS', locOut, lresult)
    if (lresult .gt. 0) then; read(locOut, *) Nstep
    else
        Nstep = 1000
        write(*,*) '`NSTEP` was set to default 1000! '
    end if
    write(*,"(' N_step      : ', I8)") Nstep
    
    ! frequency of printing .xyz
    call locPara(479, 'NX', locOut, lresult)
    if (lresult .gt. 0) then; read(locOut, *) Nxyz
    else
        Nxyz = 10
        write(*,*) '`NXYZ` was set to default 10! '
    end if
    write(*,"(' N_xyz       : ', I8)") Nxyz
    
    ! number of trajectories
    call locPara(479, 'NT', locOut, lresult)
    if (lresult .gt. 0) then; read(locOut, *) Ntraj
    else
        Ntraj = 1
        write(*,*) '`NTRAJ` was set to default 1! '
    end if
    write(*,"(' N_traj      : ', I8)") Ntraj
    
    ! temperature
    call locPara(479, 'TE', locOut, lresult)
    if (lresult .gt. 0) then; read(locOut, *) tempSI
    else
        tempSI = 298.15d0
        write(*,*) '`TEMPERATURE` was set to default 298.15K! '
    end if
    temp = tempSI / au2K
    write(*,"(' Temperature : ', F8.2 '  K,  ', F12.6, ' au')") tempSI, temp
    
    ! timestep
    call locPara(479, 'TI', locOut, lresult)
    if (lresult .gt. 0) then; read(locOut, *) tstepSI
    else
        tstepSI = 0.01d0
        write(*,*) '`TIMESTEP` was set to default 0.1fs! '
    end if
    tstep = tstepSI / au2fs
    write(*,"(' Timestep    : ', F8.2 ' fs,  ', F12.6, ' au')") tstepSI, tstep
    
    close(479)

end subroutine

! subroutine : locate the value of parameter
! line - content of line
! label - parameter name
! output - value
! lresult - result, `< 0` fail, `>= 0` succeed! 
subroutine locPara(fileID, label, output, lresult)
implicit none
integer, intent(in) :: fileID
character(len=*), intent(in) :: label
character(1024), intent(out) :: output
integer :: lresult
integer :: fileStatus
character(1024) :: lineContent, lineJudge
integer :: lenIn

    rewind(fileID)
    lenIn = len(trim(label)) ! lenght of `label`
    do while (.true.)
        read(fileID, '(A)', iostat=fileStatus) lineContent
        if (fileStatus/=0) exit
        lineJudge = adjustl(lineContent)
        if ( trim(lineJudge(1:lenIn)) .eq. trim(label) ) then
            call parseLine(lineContent, '=', output, lresult)
            lresult = 100
            return
        end if
    end do
    lresult = -1

end subroutine

subroutine parseLine(line, sep, output, lresult)
implicit none
character(len=*), intent(in) :: line, sep
character(1024), intent(out) :: output
integer, intent(out) :: lresult
integer :: ppos, length

    length = len( trim(line) )
    ppos = scan(trim(line), trim(sep), back=.true.)
    if ( ppos .gt. 0 .and. ppos .le. length ) then
        !item = line(1:ppos-1)
        output = line(ppos+1:length)
        lresult = 100
        return
    else
        lresult = -1
    end if

end subroutine


! obtained from https://stackoverflow.com/a/10759613 ! 2020-08-24 17:11:51 Wenbin, FAN @ SHU
function uppercase(strIn) result(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page

     implicit none

     character(len=*), intent(in) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i,j

     do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

end function uppercase