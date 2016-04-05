subroutine findvol(vol)
!Returns average value of density (rav) and total volume (vol) 
!from common array rho 
  implicit none
  include 'runhydro.h'
 
!*
!*  Global Variables

   real, dimension(numr,numz,numphi) :: pot, rho
   common /poisson/ pot, rho
!*
!*  Local Variables  
   real :: m, dr, vol,r,Re, counter
   integer :: i,j,count
    
    
!  print*, ">>> rhoavg"
  count=0
   
   
   dr=1.0/(ax-1.5)
   Re=1.0!(ax-1.5)/(numr-3.0)
   vol=0.0
   
   counter=2.0
   do i=2,ax
     do j=2,numz  !was by
        r=(counter-1.5)*dr
        if (rho(i,j,1).gt.0.0) then
          vol=vol+2*pi*r*dr**2
          !count=count+1
        endif
     enddo
     counter=counter+1.0
   enddo 
   
   vol=vol*2/Re**3
   
!   print*, "rhoavg vol=",vol,"volcount", count
!   print*, "rhoavg density=",rav
   
   
end subroutine findvol
   
