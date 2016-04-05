subroutine print2d(inarray,filename)
  implicit none
  include 'runhydro.h'

!************************************************************
!*
!*  Global Variables

  real, dimension(numr,numz,numphi) :: pot, rho
  common /poisson/ pot, rho

  real, dimension(numr,numz,numphi) :: psi, enth

!*
!************************************************************ 
!*
!*   Local variables
  
  integer :: i,j,k
  character(len=*) :: filename  
  
  open(unit=10,file='res.txt')
    do j=1,numz
       do i=1,numr  
          write(10,*) i,j,rho(i,j,1) 
       enddo
       write(10,*)
    enddo
  close(10)
  
end subroutine print2d
