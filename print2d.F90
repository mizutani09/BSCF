subroutine print2d(inarray,filename)
  implicit none
  include 'runhydro.h'


  real, dimension(numr,numz,numphi) :: inarray
  integer :: i,j,k
  character(len=*) :: filename  
  
  open(unit=10,file=filename)
    do j=1,numz
       do i=1,numr  
          write(10,*) i,j,inarray(i,j,1) 
       enddo
       write(10,*)
    enddo
  close(10)
  
  print*,filename
  
end subroutine print2d
