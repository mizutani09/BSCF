subroutine print1d(inarray,axis,inum,filename)
  implicit none
  include 'runhydro.h'


  real, dimension(numr,numz,numphi) :: inarray
  integer :: i,j,k,inum
  character(len=*) :: filename,axis  
  
  
  if (axis=="x") then
    open(unit=10,file=filename)
      do i=1,numz  
        write(10,*) inarray(inum,i,1) 
      enddo
    close(10) 
    
  elseif (axis=="y") then
    open(unit=10,file=filename)
      do i=1,numr  
        write(10,*) inarray(i,inum,1) 
      enddo
    close(10) 
  
  else
    print*,"Sum Ting Wong in print1d!"
  
  endif
  
  print*,filename
  
end subroutine print1d
