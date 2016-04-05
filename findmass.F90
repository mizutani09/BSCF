subroutine findmass(rho_2i,m_core,m)
!Returns total mass (m) from common array rho
  implicit none
  include 'runhydro.h'
 
!*  Global Variables
   real, dimension(numr,numz,numphi) :: pot, rho
   common /poisson/ pot, rho
!*  
   real :: m, dr,r,Re,rho_2i,m_core, counter
   integer :: i,j, count
!*   
   
   dr=1.0/(ax-1.5)
   Re=1!(ax-1.5)/(numr-3.0)
   m=0.0
   m_core=0.0
   count=0
   
   counter=2.0
   do i=2,ax
     do j=2,numz   !was by
        r=(counter-1.5)*dr        
        m=m+rho(i,j,1)*2*pi*r*dr**2
        if (rho(i,j,1).gt.rho_2i) then         
          m_core=m_core+rho(i,j,1)*2*pi*r*dr**2
        endif
     enddo
     counter=counter+1.0
   enddo
   m_core=m_core*2/Re**3
   m=m*2/Re**3
   
!   print*,"m_core", m_core
!   print*, "m", m
   
   
!   print*,"findmass mass=", m, "masscount", count 
end subroutine findmass
   
