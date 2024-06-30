!************************************************************
!*
!*  MAIN
!*
!************************************************************

program main
      implicit none
      include 'runhydro.h'

!************************************************************
!*
!*  Global Variables

      real, dimension(numr,numz,numphi) :: pot, rho
      common /poisson/ pot, rho

      real, dimension(numr,numz,numphi) :: psi

      real, dimension(numr,numz,numphi) :: enth
      common /vir/enth

      ! real, allocatable :: rho3d(:,:,:)
      ! real, allocatable :: pres3d(:,:,:)
!*
!************************************************************
!*
!*   Local variables
      real :: w_rot, phi_a, phi_b, h_a, h_b, psi_a,psi_b,phi_c
      real :: rho_c, rho_norm, h_max
      integer :: i,j,k,count
      real :: cpu1,cpu2, p_max,cput, kappa1, kappa2
      real :: phi_i, psi_i, rho_2i, gamma1, gamma2, h_2i
      real :: c1,c2,omega_sq,d_c1,d_c2,d_omega_sq,c1_old,c2_old,&
              omega_sq_old, VC, omega
      real :: T, W, P
      real :: k1,k2, re, rho_1i, h_1i, h_norm, rho_2i_norm
      character*20 char_it
      real, dimension(numr,numz,numphi) :: pres
!*
!************************************************************


      call cpu_time(cpu1)
      print*, "BSCF Started!!"

      ! allocate(rho3d(numr,numz,hydrophi))
      ! allocate(pres3d(numr,numz,hydrophi))


      rho_c = 1.0
      gamma1=1+1.0/np1
      gamma2=1+1.0/np2

!Guess the initial density
      call guessrho
!      call print2d(rho,"rho.2")
!      call print1d(rho,"y",2,"rho")

!Find rotational potential
      do i=1,numr
        do j=1,numz
          do k=1,numphi
            w_rot=(i-1.5)/(ax-1.5)
            psi(i,j,k)=-w_rot**2/2.0
          enddo
        enddo
      enddo

!      call print1d(psi,"y",2,"psi")

!Normalization
      Re=1.0!(ax-1.5)/(numr-3.0)



!!!!!!!Iterate till Convergence!!!!!!!
      d_c1=1
      d_c2=1
      d_omega_sq=1
      VC=1
      c1=0
      c2=0
      c1_old=0
      c2_old=0

      omega_sq=0
      count=0

      do while ((d_c1 .gt. eps).or.(d_c2.gt.eps).or.(d_omega_sq.gt.eps))!&
                !.or.(VC.gt.eps))
        count=count+1


!Poisson solve for density
        call poisson_solve
        pot=pot/Re**2
!        call print1d(pot,"y",2,"pot")

!Find the constants c1, c2 and omega_sq
        phi_a=pot(ax,ay,1)
        phi_b=pot(bx,by,1)
        psi_a=psi(ax,ay,1)
        psi_b=psi(bx,by,1)

        phi_i=pot(ix,2,1)
        psi_i=psi(ix,2,1)
        rho_2i=rho(ix,2,1)

        rho_1i=rho_2i*mu1/mu2

!Edited for torous        c2=phi_b

	c2=(phi_a*psi_b-phi_b*psi_a)/(psi_b-psi_a)
        omega_sq=(c2-phi_a)/psi_a

        h_2i=c2-phi_i-omega_sq*psi_i

        h_1i=h_2i*(np1+1)/(np2+1)*rho_2i/rho_1i

        c1=h_1i+phi_i+omega_sq*psi_i

!Get enthalpy
        do i=1,numr
          do j=1,numz
            if (rho(i,j,1).gt.rho_2i) then
              enth(i,j,1)=  c1 - pot(i,j,1) - omega_sq* psi(i,j,1)
              !enth(i,j,1)=enth(i,j,1)/h_norm
            else
              enth(i,j,1)=  c2 - pot(i,j,1) - omega_sq* psi(i,j,1)
            endif
          enddo
        enddo

        h_max=maxval(enth)
        rho_2i_norm=mu2/mu1*(h_1i/h_max)**np1

!Find the new normalized density
        do i=1,numr
          do j=1,numz
	      if (enth(i,j,1).gt.0) then
	         if (rho(i,j,1).gt.rho_2i) then
                   rho(i,j,1)=(enth(i,j,1)/h_max)**np1
                 else
              	   rho(i,j,1)=rho_2i_norm*(enth(i,j,1)/h_2i)**np2
                 endif
	      else
                rho(i,j,1)=0.0
              endif
          enddo
        enddo

do i=1, numz
   rho(1,i,1)=rho(2,i,1)    !Freaking oscillations!
enddo

        rho_norm=maxval(rho)

        rho=rho/rho_norm

        omega=sqrt(abs(omega_sq))
        call virial(T,W,P,omega,rho_2i)

        VC=abs(2*T-W+3*P)/abs(W)

        d_c1=abs((c1_old-c1)/c1)
        d_c2=abs((c2_old-c2)/c2)
        d_omega_sq=abs((omega_sq_old-omega_sq)/omega_sq)


        c1_old=c1
        c2_old=c2
        omega_sq_old=omega_sq

        write (char_it, "(I3)") count
        print*, "Iteration number = ",trim(char_it), ", VC = ", VC
!        print*,"T = ",T, " W = ",W," P = ",P,"omega_sq = ",omega_sq
        !print*,"c1 = ",c1, "c2 = ",c2, "omega_sq = ", omega_sq
        !print*,"d_c1 = ",d_c1, "dc_2 = ", d_c2, "d_omega_sq = ", d_omega_sq

             print*, d_c1, d_c2, d_omega_sq, VC

     enddo


        print*, d_c1, d_c2, d_omega_sq, VC
!     print*,"old rho_2i",rho_2i, "old h_max", h_max

     rho_2i=rho(ix,2,1)
     h_max=maxval(enth)

     call cpu_time(cpu2)
     cput=(cpu2-cpu1)/60.0


!Calculate and print pressure for interpolation
    kappa1 = rho_c*h_max/(np1+1.0)/rho_c**(gamma1)
    kappa2 = kappa1*rho_1i**gamma1/rho_2i**gamma2

        do i=1,numr
          do j=1,numz
             if (rho(i,j,1).gt.rho_2i) then
               pres(i,j,1) = kappa1 * rho(i,j,1)**(gamma1)
             else
               pres(i,j,1) = kappa2 * rho(i,j,1)**(gamma2)
             endif
          enddo
        enddo


    !  do i=1,numr
    !     do j=1,numz
    !        do k=1,256
    !           pres3d(i,j,k)  = pres(i,j,1)
    !        enddo
    !     enddo
    !  enddo


    ! open(unit=8,file='pressure.bin',                                   &
    !     form='unformatted',convert='BIG_ENDIAN',status='unknown')
    !    write(8) pres3d
    ! close(8)



     call getinfo(omega_sq,h_max,rho_2i,count,cput)


!     call print1d(enth,"y",2,"enth")
!     call print1d(pot,"y",2,"pot")
    !  do i=1,numr
    !     do j=1,numz
    !        do k=1,hydrophi
    !           rho3d(i,j,k)  = rho(i,j,1)
    !        enddo
    !     enddo
    !  enddo

    !  do i=1,numr
    !     do j=1,numz
    !        do k=1,hydrophi
    !           if (rho3d(i,j,k).lt. 1d-10) then
    !              rho3d(i,j,k) = 1d-10
    !           endif
    !        enddo
    !     enddo
    !  enddo

    !  call print2default(rho)
    !  call print1default(rho,"x",2)
    !  call print1default(rho,"y",2)

  !  open(unit=12,file="star1_confirm")
  !        do j=1,numz
  !          do i=1,numr
  !            write(12,*) i,j,rho(i,j,1)
  !          enddo
  !          write(12,*)
  !        enddo
  ! close(12)
  ! print*,"File star1_confirm printed"

! !Write binary output file for code initial_conditions_fc.F90
!     open(unit=8,file='density.bin',                                   &
!         form='unformatted',convert='BIG_ENDIAN',status='unknown')
!        write(8) rho3d
!     close(8)


  open(unit=12,file="star1")
         do j=1,numz
           do i=1,numr
             write(12,*) i,j,rho(i,j,1)
           enddo
           write(12,*)
         enddo
  close(12)
  print*,"File star1 printed"


  ! open(unit=12,file="star2")
  !        do j=1,numz
  !          do i=1,numr
  !            write(12,*) i,j,rho(i,j,hydrophi/2)
  !          enddo
  !          write(12,*)
  !        enddo
  ! close(12)
  ! print*,"File star2 printed"

  open(unit=12,file="pres1")
         do j=1,numz
           do i=1,numr
             write(12,*) i,j,pres(i,j,1)
           enddo
           write(12,*)
         enddo
  close(12)
  print*,"File pres1 printed"


  ! open(unit=12,file="pres2")
  !        do j=1,numz
  !          do i=1,numr
  !            write(12,*) i,j,pres(i,j,hydrophi/2)
  !          enddo
  !          write(12,*)
  !        enddo
  ! close(12)
  ! print*,"File pres2 printed"


    !  print*,"Binary file density.bin printed"
     print*,"==========================================================================="



      stop
end program main

