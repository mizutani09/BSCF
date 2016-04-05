program rho_generator
  implicit none

!Initialize  
  integer:: i,j,k,numr,numphi,numz
  real, allocatable :: rho(:,:,:),phi(:,:,:)
  double precision:: radius, Pi, angle, M, G,r, slope, c
  real:: den
  numr=100
  numz=100
  numphi=1
  radius=90.0
  Pi=3.14159265359
  den=1.0
  G=1
  
  
!Allocate arrays  
  allocate(rho(numr,numz,numphi))
  allocate(phi(numr,numz,numphi))

  
  
!Create rho array  
  slope=-80.0/90	
  c=80.0	
  do i=1,numr
    do j=1,numz
      do k=1,numphi
!        angle=k*2*Pi/numphi
!        r=((i*cos(angle))**2+(i*sin(angle))**2+(j+10)**2)**0.5
!        if (r**2.lt.radius**2) then
!          rho(i,j,k)=den   
        if (j .lt. slope * i + c) then
          rho(i,j,k)=den        
        else
          rho(i,j,k)=0.0
        endif
!	if (j==numz/2) then
!	  print*, angle
!	endif
      enddo
    enddo
  enddo
 
!Find mass
  M=den*4/3.0*Pi*radius**3	
	
  print*, M	
	
!Create phi array which is solution to the Poisson solver
  do i=1,numr
    do j=1,numz
      do k=1,numphi
        angle=k*2*Pi/numphi
        r=((i*cos(angle))**2+(i*sin(angle))**2+(j)**2)**0.5
        if (r**2.lt.radius**2) then
          phi(i,j,k)= -G*M*(3*radius**2-r**2)/2/radius**3          	  
        else
          phi(i,j,k)= -G*M/r
        endif
        
!	if (j==numz/2) then
!	  print*, angle
!	endif
      enddo
    enddo
  enddo	
	
	
	
	
	
  
!Write rho in a binary file
  open(unit=19,file= '/Users/kundan/Desktop/Patrick_2D_Poisson/serial/den',form='unformatted',status='unknown')
  write(19) rho
  close(19)

  
  
!  print*,"postwrite"
!  open(unit=10,file='/Users/kundan/Desktop/Patrick_2D_Poisson/serial/den',form='unformatted',status='unknown')
!  read(10) readrho
!  close(10)      
!  print*,"readrho"

	
!Write equatorial cross section file for gnuplot
  open(unit=10,file='zz.txt')   
        do i=1,numr          
          do k=1,numphi
            angle=k*2*Pi/numphi            
              write(10,*) i*cos(angle),i*sin(angle),rho(i,numz/2,k) !Eq cross section              	    
          enddo
          write(10,*)
        enddo   
  close(10)
  print*,"Equatorial density file zz.txt printed"
  
!Write vertical cross section file for gnuplot  
  open(unit=10,file='rho.txt')
      do j=1,numz
        do i=1,numr
              angle=1*2*Pi/numphi
              write(10,*) i*cos(angle),j,rho(i,j,1) !Vertical cross section
        enddo
        write(10,*)
      enddo
  close(10)  
  print*,"Vertical density file rho.txt printed"
  
!Write vertical cross section file for gnuplot  
  open(unit=10,file='pot.txt')
      do j=1,numz
        do i=1,numr
              angle=1*2*Pi/numphi
              write(10,*) i*cos(angle),j,phi(i,j,1) !Vertical cross section
        enddo
        write(10,*)
      enddo
  close(10)  
  print*,"Vertical potential file pot.txt printed"  
  
end program rho_generator
