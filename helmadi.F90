       subroutine helmadi(nsteps)
       implicit none
       include 'runhydro.h'
       include 'pot.h' 
!********************************************************************
!*
!  helmadi solves Poisson's equation for the self-gravity potential
!  First, Fourier transform the density and initial potential to
!  decouple the azimuthal direction.  Then use ADI method to 
!  solve for the potential and inverse Fourier transform the
!  potential to get back to physical space.
!*
!********************************************************************
!*
!*  Subroutine Arguemtns

       integer :: nsteps

!*
!********************************************************************
!*
!*  Global Variables

       real, dimension(numr,numz,numphi) :: potp, rhop
       common /potarrays/ potp, rhop

       real, dimension(numr) :: ar, cr, alphar
       real, dimension(numr,numphi) :: brb
       common /ADI_R_sweep/ ar, cr, alphar, brb

       real :: az, cz
       real, dimension(numr) :: alphaz, betaz
       real, dimension(numz) :: bzb
       real, dimension(numr,numphi) :: elambdazb
       common /ADI_Z_sweep/ az, cz, alphaz, betaz, bzb, elambdazb

       real :: gamma, piinv, four_pi
       common /pot_constants/ gamma, piinv, four_pi

       real, dimension(numr) :: rhf, r, rhfinv, rinv
       real, dimension(numz) :: zhf
       real, dimension(numphi) :: phi
       common /grid/ rhf, r, rhfinv, rinv, zhf, phi

       real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
       real, dimension(numz) :: zhf_g
       common /global_grid/ rhf_g, r_g, rhfinv_g, rinv_g, zhf_g

!       integer :: isym
!       integer, dimension(3) :: boundary_condition
!       common /boundary_conditions/ boundary_condition

!*
!********************************************************************
!*
!*  Local Variables
       
  
       real, dimension(numr,numz,numphi) :: ffrho, ffpot

       real, dimension(numr,numz,numphi) :: knownr, rhor, potr

       real, dimension(numz) :: bz

       real, dimension(numr,numphi) :: br

       real, dimension(numr,numphi) :: elambdaz

       integer :: i, j, k, l, m, index

       real :: alph, dtt, dtt_mtwogamma, angle

       real, dimension(nsteps) :: dt

!* 
!********************************************************************
! initialize the local variables
       ffrho =  0.0
       ffpot = 0.0
       knownr = 0.0
       rhor = 0.0
       potr = 0.0
       bz = 0.0
       br = 0.0
       elambdaz= 0.0
       i = 0
       j = 0
       k = 0
       l = 0
       m = 0
       index = 0
       alph = 0.0
       dtt = 0.0
       dtt_mtwogamma = 0.0 
       dt = 0.0

!       print*,"Running helmadi"
!  setup the pseudo timesteps for the ADI iterations
       dt(nsteps) = 4.0*(rinv_g(numr)*rinv_g(numr))
       alph = (r_g(numr)*rinv_g(3))**(2.0/(nsteps-1.0))
       do i = 2, nsteps
          index = nsteps + 1 - i
          dt(index) = alph*dt(index+1)
       enddo 
!       write(*,*) 'HELMADI: dt ',dt
!       write(*,*) 'HELMADI: gamma ',gamma
!       write(*,*) 'HELMADI: four_pi ',four_pi


	
	
	
!  fourier transform the density and re-order the output
!       call realft(rhop,numr,numz,numphi,+1)
       ffrho(:,:,1) = rhop(:,:,1)
!       ffrho(:,:,2:numphi_by_two) = rhop(:,:,3:numphi:2)
!       ffrho(:,:,numphi_by_two+1) = rhop(:,:,2)
!       ffrho(:,:,numphi_by_two+2:numphi) = - rhop(:,:,4:numphi:2)

!  fourier transform the potential and re-order
!       call realft(potp,numr,numz,numphi,+1)
       ffpot(:,:,1) = potp(:,:,1)
!       ffpot(:,:,2:numphi_by_two) = potp(:,:,3:numphi:2)
!       ffpot(:,:,numphi_by_two+1) = potp(:,:,2)
!       ffpot(:,:,numphi_by_two+2:numphi) = - potp(:,:,4:numphi:2)

!  swap arrays around to get data aligned for ADI iteration
       potr = ffpot
       rhor = ffrho 

!       open(unit=10,file='/scratch1/phmotl/ffpot',form='unformatted',status='unknown')
!       write(10) ffpot
!       close(10)
!       open(unit=11,file='/scratch1/phmotl/ffrho',form='unformatted',status='unknown')
!       write(11) ffrho
!       close(10)

!  ADI iteration cycle
       do i = 1, nsteps
!       do i = 1, 1

          dtt = dt(i)

          dtt_mtwogamma = dtt - 2.0*gamma
         
          !write(*,*) 'HELMADI: dtt_mtwogamma ',dtt_mtwogamma
 
          ! Radial ADI sweep

          br = brb + dtt

          if( isym /= 1 ) then
          !  the conditional looks a little stange, here is the deal...
          !  in this data decomposition, the vertical index is block
          !  distributed across numz_procs, for the global k index of
          !  two need to make a special case if isym is 2 or 3 which
          !  is that you don't include the potential at k = 1 in
          !  knownr with k = 2 and the weighting of potr with k=2
          !  is different
             knownr = 0.0
             do l = philwb, phiupb
                do k = zlwb+1, zupb
                   do j = rlwb, rupb
                      knownr(j,k,l) = -four_pi * rhor(j,k,l) +dtt_mtwogamma*potr(j,k,l) +gamma*(potr(j,k+1,l) + potr(j,k-1,l))
                   enddo
                enddo
             enddo
             ! the special case
             k = 2
             do l = philwb, phiupb
                do j = rlwb, rupb
                   knownr(j,k,l) = - four_pi * rhor(j,k,l) +(dtt-gamma)*potr(j,k,l) +gamma*potr(j,k+1,l)
                enddo
             enddo
          else
             do l = philwb, phiupb
                do k = zlwb, zupb
                   do j = rlwb, rupb
                      knownr(j,k,l) = - four_pi * rhor(j,k,l) +dtt_mtwogamma*potr(j,k,l) +gamma*(potr(j,k+1,l) + potr(j,k-1,l))
                   enddo
                enddo
             enddo
          endif

          !  add in the boundary potential on side to knownr
          knownr(numr-1,zlwb:zupb,philwb:phiupb) = knownr(numr-1,zlwb:zupb,philwb:phiupb) &
          -alphar(numr-1)*potr(numr,zlwb:zupb,philwb:phiupb)

!          open(unit=12,file='/scratch1/phmotl/knownr1',form='unformatted',status='unknown')
!          write(12) knownr
!          close(12)
!          open(unit=13,file='/scratch1/phmotl/potr1',form='unformatted',status='unknown')
!          write(13) potr
!          close(13)

          !  solve the system of equations
          call tridagr(ar,br,cr,knownr,potr)

!          open(unit=14,file='/scratch1/phmotl/knownr2',form='unformatted',status='unknown')
!          write(14) knownr
!          close(14)
!          open(unit=15,file='/scratch1/phmotl/potr2',form='unformatted',status='unknown')
!          write(15) potr
!          close(15)

          !  Vertical ADi sweep

          bz = bzb + dtt

          elambdaz = elambdazb + dtt

          ! now the special case is independent of isym and exists for
          ! the global radial index j = 2.  We have swapped having all
          ! j values in local memory to having all k values in local 
          ! memory with j block distributed across numz_procs so pe's
          ! on the bottom of the pe grid hold the special case radial
          ! index.  The special case is that knownz with j of 2 doesn't
          ! include the influence of potz with j = 1
          knownr = 0.0
          do l = philwb, phiupb
          
             do k = zlwb, zupb
                do j = rlwb+1, rupb
                   knownr(j,k,l) = - four_pi * rhor(j,k,l) +elambdaz(j,l)*potr(j,k,l) &
                   -alphaz(j)*potr(j+1,k,l) -betaz(j)*potr(j-1,k,l)
                enddo
             enddo
          enddo
          ! the sepcial case
          j = 2
          do l = philwb, phiupb
             do k = zlwb, zupb
                knownr(j,k,l) = - four_pi * rhor(j,k,l) +elambdaz(j,l)*potr(j,k,l) -alphaz(j)*potr(j+1,k,l)
             enddo
          enddo

          ! add boundary potential at top and (if isym = 1)
          ! bottom of the grid to knownz
          knownr(rlwb:rupb,zupb,philwb:phiupb) = knownr(rlwb:rupb,zupb,philwb:phiupb) +gamma*potr(rlwb:rupb,zupb+1,philwb:phiupb)

          if( isym == 1 ) then
             knownr(rlwb:rupb,zlwb,philwb:phiupb) =knownr(rlwb:rupb,zlwb,philwb:phiupb) +gamma*potr(rlwb:rupb,zlwb-1,philwb:phiupb)
          endif

!          open(unit=16,file='/scratch1/phmotl/knownr3',form='unformatted',status='unknown')
!          write(16) knownr
!          close(16)
!          open(unit=17,file='/scratch1/phmotl/potr3',form='unformatted',status='unknown')
!          write(17) potr
!          close(17)

          ! solve the system of equations
          call tridagz(az,bz,cz,knownr,potr)
             ! add the logical flag iam_on_top because
             ! the data in j is ditributed across numz_procs
             ! and will in general have to be padded to be
             ! evenly divisible.  In tridagz, use iam_on_top
             ! to set upper limits on j in do loops so we
             ! don't operate on padded array entries that
             ! are zero and unphysical

!          open(unit=17,file='/scratch1/phmotl/knownr4',form='unformatted',status='unknown')
!          write(17) knownr
!          close(17)
!          open(unit=18,file='/scratch1/phmotl/potr4',form='unformatted',status='unknown')
!          write(18) potr
!          close(18)

       enddo	! END of ADI cycle

       ffpot = potr

       ! inverse fourier transform the potential
       potp(:,:,1) = ffpot(:,:,1)
!       potp(:,:,2) = ffpot(:,:,numphi_by_two+1)
!       potp(:,:,3:numphi:2) = ffpot(:,:,2:numphi_by_two)
!       potp(:,:,4:numphi:2) = - ffpot(:,:,numphi_by_two+2:numphi)
!       call realft(potp,numr,numz,numphi,-1)
!       potp = 2.0*numphiinv*potp	! f-t normalization

!       open(unit=19,file='/Users/kundan/Desktop/Patrick_2D_Poisson/serial/potential',form='unformatted',status='unknown')
!	open(unit=19,file='potential',form='unformatted',status='unknown')
!       write(19) potp
!       close(19)

!       print*, "The symmetry used is ", isym
       
!       open(unit=10,file='ss.txt')
!         do j=1,numz
!           do i=1,numr
!              angle=1*2*Pi/numphi
!              write(10,*) i*cos(angle),j,potp(i,j,1) !Vertical cross section
!           enddo
!           write(10,*)
!         enddo
!       close(10)
  
  
!       print*,"Resultant potential file ss.txt printed"
       
       
       return
       end subroutine helmadi
