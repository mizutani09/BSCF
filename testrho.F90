subroutine testrho(id,gh,den)
  implicit none
  include 'runhydro.h'

  real, dimension(numr,numz,numphi) :: pot, rho
  common /poisson/ pot, rho
  
!Initialize  
  integer:: id, gh 
  real:: den, M, vol, mom

!  id=1
!  gh=1
!  den=1.0

!id=1 => Constant density sphere with radius ax
!id=2 => Constant density cylinder with radius ax and half-height by
!id=3 => Constant density half-cone with half-base ax and half-height by	
!gh=1 => Adds a ghost structure to the density
!den = The constant value of density
	
  print*, ">>>Testrho"
  
  rho=0.0
  
  if (id==1) then
    call sph(den,m,vol,mom)
  elseif (id==2) then
    call cyl(den,m,vol,mom)
  elseif (id==3) then
    call lin(den,m,vol,mom)
  endif
	
  if (gh==1) then
    call ghost(den) 
  endif
  
  
  call print2d(rho,"test2d")
  call print1d(rho,"y",2,"tes1d")
  print*, "Total mass=", m, "Total volume=",vol,"MI=",mom
  
end subroutine testrho





!!!! Functions!!!!
!!Sphere
subroutine sph(den,m,vol,mom)
  implicit none
  include 'runhydro.h'
  real, dimension(numr,numz,numphi) :: pot, rho
  common /poisson/ pot, rho
  
  integer:: i,j,k,count
  real :: radius, m, vol, mom,angle,den,r,c
  
  radius=(ax-1.5)	
  c=0
  count=0
  
  rho=0.0
  
  do i=1,numr
    do j=1,numz
        r=((i-1.5)**2+(j+c-1.5)**2)**0.5
        if ((r**2.le.radius**2).and.(j.ge.2).and.(i.ge.2)) then
          rho(i,j,1)=den      
          count=count+1
          
          
          
        else
          rho(i,j,1)=0.0
        endif
        
                  
      !    if (j==82) then
      !      print*,rho(i,j,1), i,j
      !    endif
        
        
        
    enddo
  enddo
  
  radius=(ax-1.5-1)/(numr-1)
  vol=4.0/3*pi*radius**3
  m=den*vol
  mom=2.0/5*m*radius**2
       print*, "anacount", count
end subroutine sph


!!Cylinder
subroutine cyl(den,m,vol,mom)
  implicit none
  include 'runhydro.h'
  real, dimension(numr,numz,numphi) :: pot, rho
  common /poisson/ pot, rho
  
  integer:: i,j,k, count
  real :: radius, height, m, vol, mom,den

  count=0
  
  do i=2,numr
    do j=2,numz
      if ((i.le.ax).and.(j.le.by)) then
        rho(i,j,1)=den
        
        count=count+1
        
      else
        rho(i,j,1)=0.0
      endif
    enddo
  enddo
  
  radius=(ax-1.5-1)/(numr-1)
  height=(by-1.5-1)/(numr-1)
  
  vol=pi*radius**2*height*2
  m=den*vol
  mom=(m*radius**2)/2.0  
  print*, "anacount", count
end subroutine cyl


!!Cone
subroutine lin(den,m,vol,mom)
  implicit none
  include 'runhydro.h'
  real, dimension(numr,numz,numphi) :: pot, rho
  common /poisson/ pot, rho

  integer:: i,j,k
  real :: slope, c, count
  real :: radius, height, m, vol, mom, den
  
  slope=-by*1.0/ax	
  c=by*1.0  
  
  
  count=0
  radius=(ax-1.5)
  
  do i=2,numr
    do j=2,numz
      do k=1,numphi
        if (j .le. slope * i + c) then
          rho(i,j,k)=den     
          count=count+1
        else
          rho(i,j,k)=0.0
        endif
      enddo
    enddo
  enddo
  
  radius=(ax-1.5-1)/(numr-1)
  height=(by-1.5-1)/(numr-1)
  
  vol=pi/3.0*radius**2*height*2
  m=den*vol
  mom=0.3*(m*radius**2)  !Not *2 because vol takes care of it  
  print*, "anacount", count
  
end subroutine lin


!!Ghost
subroutine ghost(den)
  implicit none
  include 'runhydro.h'
  real, dimension(numr,numz,numphi) :: pot, rho
  common /poisson/ pot, rho
  
  integer:: i,j,k
  real :: slope, c, den
  
  slope=1	  
  
  do i=1,numr
    do j=1,numz
      if (j .lt. slope *(i-ax-1)+1) then
          rho(i,j,1)=den      
      endif
    enddo
  enddo
  
end subroutine ghost
