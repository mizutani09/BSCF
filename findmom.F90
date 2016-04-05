subroutine findmom(rho_2i,mom,kt,kc,ke)
!Returns moment of inertia (rav) NOT angular momentum 
!from common array rho
  implicit none
  include 'runhydro.h'
 
!*
!*  Global Variables

   real, dimension(numr,numz,numphi) :: pot, rho
   common /poisson/ pot, rho
!*
!*  Local Variables  
   real, dimension(numr,numz,numphi) :: psi
   real :: m, dr, mom, r, Re, count
   integer :: i,j
   real, intent(in) :: rho_2i
   real, intent(out) :: kt, kc, ke
   
   
   dr=1.0/(ax-1.5)
   Re=1.0!(ax-1.5)/(numr-3.0)
   m=0.0
   mom=0.0
   kt=0.0
   kc=0.0
   ke=0.0
   
   count=2.0
   do i=2,ax
     do j=2,numz  !was by
        r=(count-1.5)*dr
        m=rho(i,j,1)*2*pi*r*dr**2
        mom=mom+m*r**2
          kt = kt + m*r**2
        if (rho(i,j,1).gt.rho_2i) then
          kc = kc + m*r**2
        else
          ke = ke + m*r**2
        endif
     enddo
     count=count+1.0
   enddo
   
   mom=mom*2/Re**5

   kt = 2*kt
   kc = 2*kc
   ke = 2*ke
!   print*,"findmom mom=",mom
   
end subroutine findmom
