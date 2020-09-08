subroutine logHeader
implicit none

    call logHline
    call logging("             QCAT : A simple code for quasi-classical trajectory")
    call logging()
    call logging("    Author : Wenbin, FAN (fanwenbin@shu.edu.cn)")
    call logging("    Date : 2020-Aug-24 ")
    call logHline

end subroutine


subroutine logging(content)
use fileIO
implicit none
character(len=*), intent(in) :: content

    write(*,*) content
    !write(logFileID, *) content

end subroutine

! log a horizontal line, width is 80 characters
subroutine logHline
use fileIO
implicit none

    call logging("--------- --------- --------- --------- --------- --------- --------- ---------")

end subroutine

! log error and exit whole program
subroutine logError(content)
implicit none
character(len=*), intent(in) :: content

    call logging('[ERROR] '//trim(content))
    stop

end subroutine

subroutine logWarn(content)
implicit none
character(len=*), intent(in) :: content

    call logging('[WARNING] '//trim(content))

end subroutine

subroutine logXYZ(Natoms, coord, eleName)
integer, intent(in) :: Natoms
real(8), intent(inout) :: coord(3, Natoms)
character(2), intent(in) :: eleName(Natoms)
integer :: i

    do i = 1, Natoms
        write(*,'(2X,A2,3F12.6)') trim(eleName(i)), coord(1:3,i)
    end do

end subroutine