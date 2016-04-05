subroutine virial(T,W,P,omega,rho_2i)
  implicit none
  include 'runhydro.h'

  real, dimension(numr,numz,numphi) :: pot, rho
  common /poisson/ pot, rho  

  real, dimension(numr,numz,numphi) :: enth
  common /vir/enth  
  
  real :: T, W, P, omega, r, m, dr, press, Re, rho_2i, count
  integer :: i,j
  
  
   dr=1.0/(ax-1.5)
   Re=1.0!(ax-1.5)/(numr-3.0)
!  print*, dr
  !Find rotational energy T
  T=0.0
  count=2.0
  do i=2,ax
     do j=2,numz  !was by
        r=(count-1.5)*dr
        m=rho(i,j,1)*2*pi*r*dr**2
        T=T+0.5*m*r**2*omega**2
     enddo
     count=count+1.0
   enddo  
   T=T*2/Re**5

  !Find Potential energy W
  W=0.0
  count=2.0
  do i=2,ax
     do j=2,numz  !was by
        r=(count-1.5)*dr
        W=W-0.5*pot(i,j,1)*rho(i,j,1)*2*pi*r*dr**2
     enddo
     count=count+1.0
  enddo   	  
  W=W*2/Re**3  

  !Find pressure energy P
  P=0.0
  count=2.0
  do i=2,ax
     do j=2,numz  !was by
        r=(count-1.5)*dr
        if (rho(i,j,1).gt.rho_2i) then
          press=rho(i,j,1)*enth(i,j,1)/(1.0+np1)
        else
          press=rho(i,j,1)*enth(i,j,1)/(1.0+np2) 
        endif
        P=P+press*2*pi*r*dr**2
     enddo
     count=count+1.0
   enddo  
  P=P*2/Re**3 
  return
end subroutine virial
