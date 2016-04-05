       subroutine potsetup
       implicit none
       include 'runhydro.h'
       include 'pot.h'
!********************************************************************
!*
!  potsetup does intitialization for the poisson solver
!  package.  In order of occurance this is what potsetup
!  does:
!
!     -> initialize trigonometric functions used by bessel
!
!     -> call tm and sm to fill in the Green function tmr and smz
!        respectively
!
!     -> setup the geometric arrays for the ADI solver
!
!*
!********************************************************************
!*
!*  Global Variables

       real, dimension(numr,numz,numr,mmax) :: tmr
       real, dimension(numr,numz,numz,mmax) :: smz
       common /green_functions/ tmr, smz

       real, dimension(numphi,mmax) :: bes_cos, bes_sin
       common /bessel_trig/ bes_cos, bes_sin

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

      real ::  grav
      common /constants/ grav

       real :: dr, dz, dphi, drinv, dzinv, dphiinv
       common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

       real, dimension(numr) :: rhf, r, rhfinv, rinv
       real, dimension(numz) :: zhf
       real, dimension(numphi) :: phi
       common /grid/ rhf, r, rhfinv, rinv, zhf, phi

       real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
       real, dimension(numz) :: zhf_g
       common /global_grid/ rhf_g, r_g, rhfinv_g, rinv_g, zhf_g

       real, dimension(numr,numz,numphi) :: pot, rho
       common /poisson/ pot, rho

!       integer, dimension(3) :: boundary_condition
!       common /boundary_conditions/  boundary_condition

!*
!********************************************************************
!*
!*  Local Variables
  
       real :: mindex, lindex, drinv2, dphiinv2

       real, dimension(numphi) :: elm, m1mode

       real, dimension(numr) :: betar

       integer :: j, k, l, m, lstop, index, mode

!*
!********************************************************************
!  initialize the local variables

       mindex = 0.0
       lindex = 0.0
       drinv2 = 0.0
       dphiinv2 = 0.0
       elm = 0.0
       m1mode = 0.0
       betar = 0.0
       j = 0
       k = 0
       l = 0
       m = 0
       lstop = 0
       index = 0
       mode = 0 

       piinv = 1.0 / pi
       four_pi = 4.0 * pi
!       write(*,*) 'POTSETUP: pi ',pi
!       write(*,*) 'POTSETUP: piinv ',piinv
!       write(*,*) 'POTSETUP: four_pi ', four_pi

!  setup the trig arrays used by bessel, use lindex and mindex to
!  avoid casts which are relatively expensive
       mindex = 1.0
       lindex = 1.0
       if( isym == 3 ) then
          do m = 1, mmax
             do l = philwb, phiupb
                bes_cos(l,m) = cos(dphi*(mindex-1.0)*(2.0*lindex-1.0))
                bes_sin(l,m) = sin(dphi*(mindex-1.0)*(2.0*lindex-1.0))
                lindex = lindex + 1.0
             enddo
             lindex = 1.0
             mindex = mindex + 1.0
          enddo
       else
          do m = 1, mmax
             do l = philwb, phiupb
                bes_cos(l,m) = cos(0.5*dphi*(mindex-1.0)*(2.0*lindex-1.0))
                bes_sin(l,m) = sin(0.5*dphi*(mindex-1.0)*(2.0*lindex-1.0))
                lindex = lindex + 1.0
             enddo
             lindex = 1.0
             mindex = mindex + 1.0
          enddo
       endif
!       open(unit=10,file='bes_sin',form='unformatted',status='unknown')
!       write(10) bes_sin
!       close(10)
!       open(unit=11,file='bes_cos',form='unformatted',status='unknown')
!       write(11) bes_cos
!       close(11)

!  call tm and sm to fill in the Green functions for the top (and bottom)
!  potential - tmr and the Green function for the side potential - smz
      call tm(tmr)
      call sm(smz)
!      open(unit=10,file='/scratch1/phmotl/tmr',form='unformatted',status='unknown')
!      read(10) tmr
!      close(10)
!      open(unit=11,file='/scratch1/phmotl/smz',form='unformatted',status='unknown')
!      read(11) smz
!      close(11)

!  the rest of the code in potsetup sets up arrays that hold geometric
!  information for the discretization of Poisson's equation.  The
!  arrays are all used by the ADI solver
       gamma = dzinv * dzinv
       drinv2 = drinv * drinv 
       dphiinv2 = dphiinv * dphiinv
!       write(*,*) 'POTSETUP: gamma ',gamma
!       write(*,*) 'POSTETUP: drinv2 ',drinv2
!       write(*,*) 'POTSETUP: dphiinv2 ',dphiinv2

       lstop = numphi_by_two + 1
       do l = philwb, phiupb
          if( isym == 3 ) then
             if( l <= lstop ) then
                mode = (l-1)*(l-1)
             else
                mode = (l-lstop)*(l-lstop)     
             endif
          else
             if( l <= lstop ) then
                mode = l - 1
             else
                mode = l - lstop
             endif
          endif
          m1mode(l) = (-1.0)**mode
          elm(l) = cos(mode*dphi)
       enddo

!       open(unit=10,file='m1mode',form='unformatted',status='unknown')
!       write(10) m1mode
!       close(10)
!       open(unit=11,file='elm',form='unformatted',status='unknown')
!       write(11) elm
!       close(11)
      
       do j = rlwb, rupb

          alphar(j) = -r_g(j+1)*rhfinv_g(j)*drinv2

          betar(j) = -r_g(j)*rhfinv_g(j)*drinv2

       enddo

!       open(unit=12,file='alphar',form='unformatted',status='unknown')
!       write(12) alphar
!       close(12)
!       open(unit=13,file='betar',form='unformatted',status='unknown')
!       write(13) betar
!       close(13)

       ! have to used indexed global radius array to initialize 
       ! alphaz and betaz because they will be used when the radial
       ! data is block distributed across numz_procs and numr_dd
       ! is not in general equal to numr_dd_pad
       do j = rlwb, rupb

          alphaz(j) = -r_g(j+1)*rhfinv_g(j)*drinv2

          betaz(j) = -r_g(j)*rhfinv_g(j)*drinv2

       enddo

!       open(unit=14,file='alphaz',form='unformatted',status='unknown')      
!       write(14) alphaz
!       close(14)
!       open(unit=15,file='betaz',form='unformatted',status='unknown')
!       write(15) betaz
!       close(15)
 
       ar = betar

       cr = alphar

       az = - gamma
 
       cz = - gamma

       do l = philwb, phiupb
          do j = rlwb + 1, rupb
             brb(j,l) = 2.0*drinv2 - 2.0*(elm(l)-1.0)*dphiinv2*rhfinv_g(j)*rhfinv_g(j)
          enddo
       enddo
       if( isym == 3 ) then
          do l = philwb, phiupb
             brb(2,l) = -alphar(2) - 2.0*betar(2) - 2.0*(elm(l)-1.0)*dphiinv2*rhfinv_g(2)*rhfinv_g(2)
          enddo
       else
          do l = philwb, phiupb
             brb(2,l) = -alphar(2) + (m1mode(l)-1.0)*betar(2) -2.0*(elm(l)-1.0)*dphiinv2*rhfinv_g(2)*rhfinv_g(2)
          enddo
       endif
!       open(unit=17,file='brb',status='unknown',form='unformatted')
!       write(17) brb
!       close(17) 
       do k = zlwb, zupb
          bzb(k) = 2.0*gamma
       enddo
       if( isym /= 1 ) bzb(2) = gamma
!       open(unit=18,file='bzb',form='unformatted',status='unknown')
!       write(18) bzb
!       close(18)
       do l = philwb, phiupb
          do j = rlwb, rupb
             elambdazb(j,l) = -2.0*drinv2 + 2.0*(elm(l)-1.0)*dphiinv2*rhfinv_g(j)*rhfinv_g(j)
          enddo
       enddo
       if( isym == 3 ) then
          do l = philwb, phiupb
             elambdazb(2,l) = alphaz(2) + 2.0*betaz(2) + 2.0*(elm(l)-1.0)*dphiinv2*rhfinv(2)*rhfinv(2)
          enddo
       endif
       if( isym /= 3 ) then
          do l = philwb, phiupb
             elambdazb(2,l) = alphaz(2) - (m1mode(l)-1.0)*betaz(2) +2.0*(elm(l)-1.0)*dphiinv2*rhfinv(2)*rhfinv(2)
          enddo
       endif
!       open(unit=19,file='elambdazb',form='unformatted',status='unknown')
!       write(19) elambdazb
!       close(19) 
!       print*,"mmax potsetup",mmax
       return
       end subroutine potsetup
