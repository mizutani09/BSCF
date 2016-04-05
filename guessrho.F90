subroutine guessrho
  implicit none
  include 'runhydro.h'

  real, dimension(numr,numz,numphi) :: pot, rho
  common /poisson/ pot, rho
  
  real, dimension(numr,numz,numphi) :: phi, phi1
!Initialize  
  integer:: i,j,k
  double precision:: M, G,r, slope, c, rb,imx,imy,m1,m2
  double precision:: c1,c2
  real:: den

  den=1.0
  G=1.0
  rb=by*1.0/ax

  
!Create rho array  
  if (bx.eq.2) then 
  !Object is not a toroid
    if ((rb.lt.0.4).and.(rb.ge.0.15)) then
    !Object is flat enough that linear guess fails to reach the boundary
      do i=2,numr
        do j=2,numz
          if ((i.le.ax).and.(j.le.by)) then
            rho(i,j,1)=den
          else
            rho(i,j,1)=0.0
          endif
        enddo
      enddo

    elseif (rb.lt.0.15) then
    !Object is so flat that even square/cylindrical guess may not work
      imx=(ax-bx)/2.0
      imy=ax*0.15

      m1=(imy-by)/(imx-bx)
      m2=(imy-ay)/(imx-ax)
      c1=by*1.0-m1*bx
      c2=ay*1.0-m1*ax

      do i=1,numr
        do j=1,numz
          if ((j.le.m1*i+c1).and.(i.le.ax).and.(j.le.0.15*ax)) then

            rho(i,j,1)=den
          else
            rho(i,j,1)=0.0
          endif
        enddo
      enddo
 
    else
    !Object is not highly distorted, linear guess
      slope=-by*1.0/ax	
      c=by*1.0  
      do i=1,numr
        do j=1,numz
          do k=1,numphi
            if (j .lt. slope * i + c) then
              rho(i,j,k)=den       
            else
              rho(i,j,k)=0.0
            endif
          enddo
        enddo
      enddo
    endif

  else
  !Guess for a toroid
    do i=1,numr
      do j=1,numz
        if ((j.lt.0.15*ax).and.(i.lt.ax).and.(i.gt.bx)) then
          rho(i,j,1)=den
        else
          rho(i,j,1)=0.0
        endif
      enddo
    enddo

  endif



  return
end subroutine guessrho
