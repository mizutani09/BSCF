subroutine rhoavg(rav)
!Returns average value of density (rav) from common array rho 
  implicit none
  include 'runhydro.h'
 
!*
!*  Global Variables

   real, dimension(numr,numz,numphi) :: pot, rho
   common /poisson/ pot, rho
!*
!*  Local Variables  
   real :: m, dr, rav,pi,vol,r, counter
   integer :: i,j,count
    
    
  print*, ">>> rhoavg"
  Pi=3.14159265359
  count=0
   
   
   dr=1.0/(ax-1.5)
   vol=0.0
   
   counter=2.0
   do i=2,ax
     do j=2,by
        r=(i-1.5)*dr
        if (rho(i,j,1).gt.0.0) then
          vol=vol+2*pi*r*dr**2
          count=count+1
        endif
     enddo
     counter=counter+1.0
   enddo 
   
   vol=vol*2
   
   call findmass(m)
   
   
   rav=m/vol
   
   print*, "rhoavg vol=",vol,"volcount", count
   print*, "rhoavg density=",rav
   
   
end subroutine rhoavg
   
