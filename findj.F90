subroutine findj(m, J2, J4, J6)
!Returns zonal harmonics coefficients of gravitational potential
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
   real :: m, dr, J2, J4, J6, r, Re, c1, c2, dtheta, theta, Rcyl, dz, &
          cos_theta1, cos_theta2, x1, x2, jj2, jj4, jj6
   integer :: i,j
   
   
   dr=1.0/(ax-1.5)
   dz = dr
!   dtheta=pi/
   Re=1.0
   
   J2=0.0
   J4=0.0
   J6=0.0
   jj2=0
   jj4=0
   jj6=0
   
   c1=2.0

   do i = 2, numr!ax
     c2=2.0
     do j = 2, numz  !was by
        
        r = ( (c1-1.5)**2 +(c2-1.5)**2 )**(0.5)  * dr
        Rcyl = (c1-1.5) /(ax-1.5)
!        theta = atan((c2-1.5)/(c1-1.5))
        theta = atan((c1-1.5)/(c2-1.5))
        x1 = cos(theta)! Rcyl/r
        x2 = cos(pi-theta)
         
        JJ2 = JJ2 + rho(i,j,1)*(1+3*cos(2*theta))/4 * r**2 * Rcyl 

        j2 = j2+rho(i,j,1) * ( 3* x1**2 -1 )/2.0 * &
                  r**2 * Rcyl +                      & 
                  rho(i,j,1) * ( 3* x2**2 -1 )/2.0 * &
                  r**2 * Rcyl 
 
        JJ4 = JJ4 + rho(i,j,1)*(9+20*cos(2*theta)+35*cos(4*theta))/64 * r**4 * Rcyl

        j4 = j4+rho(i,j,1) * ( 35*x1**4 - 30*      &
                  x1**2 + 3 )/8.0 * r**4 * Rcyl +    &
                  rho(i,j,1) * ( 35*x2**4 - 30*      &
                  x2**2 + 3 )/8.0 * r**4 * Rcyl 

        JJ6 = JJ6 + rho(i,j,1)*(50+105*cos(2*theta)+126*cos(4*theta)+231*cos(6*theta))/512 &
                 * r**6 * Rcyl 

        j6 = j6+rho(i,j,1) * ( 231*x1**6 - 315*x1**4 +&
                  105*x1**2 -5 )/16.0 * r**6 * Rcyl +   &
                  rho(i,j,1) * ( 231*x2**6 - 315*x2**4 +&
                  105*x2**2 -5 )/ 16.0 * r**6 * Rcyl 

        c2=c2+1.0
     enddo
     c1=c1+1.0
   enddo


   J2 = J2 * (-2 * pi)/ (m * (Re)**2) * dr * dz
   J4 = J4 * (-2 * pi)/ (m * (Re)**4) * dr * dz
   J6 = J6 * (-2 * pi)/ (m * (Re)**6) * dr * dz

   JJ2 = JJ2 * (-4 * pi)/ (m * (Re)**2) * dr * dz
   JJ4 = JJ4 * (-4 * pi)/ (m * (Re)**4) * dr * dz
   JJ6 = JJ6 * (-4 * pi)/ (m * (Re)**6) * dr * dz

   
end subroutine findj
