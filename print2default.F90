subroutine print2default(inarray)
  implicit none
  include 'runhydro.h'

  real, dimension(numr,numz,numphi) :: inarray
  integer :: i,j,k
  character(len=100) :: filename
  character*20 char_np1, char_ax, char_by, char_numr, char_bx,char_ix, &
  char_np2, char_mu1, char_mu2, char_rcore, char_numz   
  
!!Convert numbers to strings
  if (ax.lt.100) then
    write (char_ax, "(I2)") ax
  else
    write (char_ax, "(I3)") ax
  endif

  if ((by.lt.100).and.(by.gt.9)) then
    write (char_by, "(I2)") by
  elseif (by.lt.10) then  
    write (char_by, "(I1)") by
  else
    write (char_by, "(I3)") by
  endif

  if ((bx.lt.100).and.(bx.gt.10)) then
    write (char_bx, "(I2)") bx
  elseif (bx.lt.10) then  
    write (char_bx, "(I1)") bx
  else
    write (char_bx, "(I3)") bx
  endif 

  write (char_np1, "(F4.1)") np1
  write (char_np2, "(F4.1)") np2
  write (char_numr, "(I3)") numr
  write (char_numz, "(I3)") numz
  write (char_ix, "(I2)") ix
  write (char_mu1,"(F4.1)") mu1
  write (char_mu2,"(F4.1)") mu2
  write (char_rcore,"(F5.2)") ix*1.0/ax  

!!Make filename	
  if (bx==2) then	
    filename='Bi_'//trim(char_ix)//"_"//trim(char_np1)//'w'//trim(char_mu1)&
    //'_'//trim(char_np2)//'w'//trim(char_mu2)//'_'//trim(char_ax)//"x"//&
    trim(char_by)//"_"//trim(char_numr)//".2"   
  else
    filename='Bi_-'//trim(char_ix)//"_"//trim(char_np1)//'w'//trim(char_mu1)&
    //'_'//trim(char_np2)//'w'//trim(char_mu2)//'_'//trim(char_ax)//"x"//&
    trim(char_by)//"_"//trim(char_numr)//".2" 
  endif
  
!!Write array
  open(unit=10,file=filename)
    do j=1,numz
       do i=1,numr
          write(10,*) i,j,inarray(i,j,1)
       enddo
       write(10,*)
    enddo
  close(10)

  print*, trim(filename)

end subroutine print2default
