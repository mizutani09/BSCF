       subroutine bessel
       implicit none
       include 'runhydro.h'
       include 'pot.h'
!****************************************************************
!*
!  bessel calculates the boundary values of the potential
!  using Howie's expansion in terms of bessel functions.
!  In essence you multiply the density field by the appropriate
!  cylindrical Green's function and sum it up over the grid.
!*
!****************************************************************
!*
!*  Global Variables

       real, dimension(numr,numz,numphi) :: pot, rho
       common /poisson/ pot, rho

       real, dimension(numr,numz,numphi) :: potp, rhop
       common /potarrays/ potp, rhop

       real, dimension(numr,numz,numr,mmax) :: tmr
       real, dimension(numr,numz,numz,mmax) :: smz
       common /green_functions/ tmr, smz

       real, dimension(numphi,mmax) :: bes_cos, bes_sin
       common /bessel_trig/ bes_cos, bes_sin

       real, dimension(numr) :: rhf, r, rhfinv, rinv
       real, dimension(numz) :: zhf
       real, dimension(numphi) :: phi
       common /grid/ rhf, r, rhfinv, rinv, zhf, phi

       real :: dr, dz, dphi, drinv, dzinv, dphiinv
       common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv
    
!       integer :: isym
!       integer, dimension(3) :: boundary_condition
!       common /boundary_conditions/ boundary_condition

!*
!****************************************************************
!*
!*  Local Variables

       real, dimension(numr,numz) :: TMPC, TMPS

       real, dimension(numr,numphi) :: phitTMP, pott, phitTMPC, phitTMPS

       real, dimension(numr,numphi) :: phibTMP, potb, phibTMPC, phibTMPS

       real, dimension(numz,numphi) :: phisTMP, pots, phisTMPC, phisTMPS

       real, dimension(numr,mmax) :: StC, StS

       real, dimension(numr,mmax) :: SbC, SbS

       real, dimension(numz,mmax) :: SsC, SsS

       real, dimension(numr,mmax) :: sum_top_C, sum_top_S

       real, dimension(numr,mmax) :: sum_bot_C, sum_bot_S

       real, dimension(numz,mmax) :: sum_sid_C, sum_sid_S

       real :: factor, angle

       integer :: i, j, k, l, m, lwrb, uprb

       
!*
!****************************************************************
! initialize the local variables
       TMPC = 0.0
       TMPS = 0.0
       phitTMP = 0.0
       pott = 0.0
       phitTMPC = 0.0
       phitTMPS = 0.0
       phibTMP = 0.0
       potb = 0.0
       phibTMPC = 0.0
       phibTMPS = 0.0
       phisTMP = 0.0
       pots = 0.0
       phisTMPC = 0.0
       phisTMPS = 0.0
       StC = 0.0
       StS = 0.0
       SbC = 0.0
       SbS = 0.0
       SsC = 0.0
       SsS = 0.0
       sum_top_C = 0.0
       sum_top_S = 0.0
       sum_bot_C = 0.0
       sum_bot_S = 0.0
       sum_sid_C = 0.0
       sum_sid_S = 0.0
       factor = 0.0
       i = 0
       j = 0
       k = 0
       l = 0 
       m = 0
       lwrb = 0
       uprb = 0
       
!       print*,"Running Bessel"
!  factor is the common multiplier for converting summation of
!  Green function times the density to a potential
       if( isym == 3 ) then
          factor = - 2.0 * dr * dz * dphi
       else
          factor = - dr * dz * dphi
       endif

!  evaluate the m=0 contribution to top, side and bottom slices
!  of the potential
       do k = zlwb, zupb
          do j = rlwb, rupb
             do l = philwb, phiupb
                TMPC(j,k) = TMPC(j,k) + rho(j,k,l)
             enddo
          enddo
       enddo
       if( isym == 1 ) then
          do j = 2, rupb
             StC(j,1) = sum(tmr(rlwb:rupb,zlwb:zupb,j,1)* TMPC(rlwb:rupb,zlwb:zupb))
             SbC(j,1) = sum(tmr(rlwb:rupb,zupb:zlwb:-1,j,1)*TMPC(rlwb:rupb,zlwb:zupb))
          enddo
          do k = 2, zupb 
             SsC(k,1) = sum(smz(rlwb:rupb,zlwb:zupb,k,1)*TMPC(rlwb:rupb,zlwb:zupb))
          enddo
       else
          do j = 2, rupb
             StC(j,1) = sum(tmr(rlwb:rupb,zlwb:zupb,j,1)* TMPC(rlwb:rupb,zlwb:zupb))
          enddo 
          do k = 2, zupb
             SsC(k,1) = sum(smz(rlwb:rupb,zlwb:zupb,k,1)*TMPC(rlwb:rupb,zlwb:zupb))
          enddo
       endif

!  now compute the contributions to the boundary potential for
!  modes with m > 0
       if( isym == 1 ) then
          do m = 2, mmax
             TMPC = 0.0
             TMPS = 0.0
             do k = zlwb, zupb
                do j = rlwb, rupb
                   do l = philwb, phiupb
                      TMPC(j,k) = TMPC(j,k) + rho(j,k,l)*bes_cos(l,m)
                      TMPS(j,k) = TMPS(j,k) + rho(j,k,l)*bes_sin(l,m)
                   enddo
                enddo
             enddo
             do j = 2, rupb
                StC(j,m) = sum(tmr(rlwb:rupb,zlwb:zupb,j,m)* TMPC(rlwb:rupb,zlwb:zupb))
                StS(j,m) = sum(tmr(rlwb:rupb,zlwb:zupb,j,m)* TMPS(rlwb:rupb,zlwb:zupb))
                SbC(j,m) = sum(tmr(rlwb:rupb,zupb:zlwb:-1,j,m)* TMPC(rlwb:rupb,zlwb:zupb))
                SbS(j,m) = sum(tmr(rlwb:rupb,zupb:zlwb:-1,j,m)* TMPS(rlwb:rupb,zlwb:zupb))
             enddo
             do k = 2, zupb
                SsC(k,m) = sum(smz(rlwb:rupb,zlwb:zupb,k,m)* TMPC(rlwb:rupb,zlwb:zupb))
                SsS(k,m) = sum(smz(rlwb:rupb,zlwb:zupb,k,m)* TMPS(rlwb:rupb,zlwb:zupb))
             enddo
          enddo
       else
          do m = 2, mmax
             TMPC = 0.0
             TMPS = 0.0
             do k = zlwb, zupb
                do j = rlwb, rupb
                   do l = philwb, phiupb
                      TMPC(j,k) = TMPC(j,k) + rho(j,k,l)*bes_cos(l,m)
                      TMPS(j,k) = TMPS(j,k) + rho(j,k,l)*bes_sin(l,m)
                   enddo
                enddo
             enddo
             do j = 2, rupb
                StC(j,m) = sum(tmr(rlwb:rupb,zlwb:zupb,j,m)* TMPC(rlwb:rupb,zlwb:zupb))
                StS(j,m) = sum(tmr(rlwb:rupb,zlwb:zupb,j,m)* TMPS(rlwb:rupb,zlwb:zupb))
             enddo
             do k = 2, zupb
                SsC(k,m) = sum(smz(rlwb:rupb,zlwb:zupb,k,m)* TMPC(rlwb:rupb,zlwb:zupb))
                SsS(k,m) = sum(smz(rlwb:rupb,zlwb:zupb,k,m)* TMPS(rlwb:rupb,zlwb:zupb))
             enddo
          enddo
       endif

       if( isym == 1 ) then
          sum_top_C = StC
          sum_top_S = StS
          sum_bot_C = SbC
          sum_bot_S = SbS
          sum_sid_C = SsC
          sum_sid_S = SsS
       else
          sum_top_C = StC
          sum_top_S = StS
          sum_sid_C = SsC
          sum_sid_S = SsS
       endif

       do l = philwb, phiupb
          do j = rlwb, rupb
             phitTMP(j,l) = sum_top_C(j,1)
          enddo
       enddo
       do m = 2, mmax
          do l = philwb, phiupb
             do j = rlwb, rupb
                phitTMPC(j,l) = sum_top_C(j,m)*bes_cos(l,m)
                phitTMPS(j,l) = sum_top_S(j,m)*bes_sin(l,m)
             enddo
          enddo
          do j = rlwb, rupb
             phitTMP(j,:) = phitTMP(j,:) + 2.0*phitTMPC(j,:) + 2.0*phitTMPS(j,:)
          enddo
       enddo
       pott(rlwb:rupb,:) = factor * phitTMP(rlwb:rupb,:)
       if( isym == 3 ) then
          pott(1,:) = pott(2,:)
       else
          do l = 1, numphi_by_two
             pott(1,l) = pott(2,l+numphi_by_two)
             pott(1,l+numphi_by_two) = pott(2,l)
          enddo
       endif
       potp(:,numz,:) = pott
       pot(:,numz,:) = pott

       do l = philwb, phiupb
          do k = zlwb, zupb
             phisTMP(k,l) = sum_sid_C(k,1)
          enddo
       enddo
       do m = 2, mmax
          do l = philwb, phiupb
             do k = zlwb, zupb
                phisTMPC(k,l) = sum_sid_C(k,m)*bes_cos(l,m)
                phisTMPS(k,l) = sum_sid_S(k,m)*bes_sin(l,m)
             enddo
          enddo
          do k = zlwb, zupb
             phisTMP(k,:) = phisTMP(k,:) + 2.0*phisTMPC(k,:) + 2.0*phisTMPS(k,:)
         enddo
       enddo
       pots(zlwb:zupb,:) = factor*phisTMP(zlwb:zupb,:)
       if( isym /= 1 ) then
          pots(1,:) = pots(2,:)
       endif
       potp(numr,:,:) = pots
       pot(numr,:,:) = pots

       if( isym == 1 ) then
          do l = philwb, phiupb
             do j = rlwb, rupb
                phibTMP(j,l) = sum_bot_C(j,1)
             enddo
          enddo
          do m = 2, mmax
             do l = philwb, phiupb
                do j = rlwb, rupb
                   phibTMPC(j,l) = sum_bot_C(j,m)*bes_cos(l,m)
                   phibTMPS(j,l) = sum_bot_S(j,m)*bes_sin(l,m)
                enddo
             enddo
             do j = rlwb, rupb
                phibTMP(j,:) = phibTMP(j,:) + 2.0*phibTMPC(j,:) + 2.0*phibTMPS(j,:)
             enddo
          enddo
          potb(rlwb:rupb,:) = factor * phibTMP(rlwb:rupb,:)
          do l = 1, numphi_by_two
             potb(1,l) = potb(2,l+numphi_by_two)
             potb(1,l+numphi_by_two) = potb(2,l)
          enddo
          potp(:,1,:) = potb
          pot(:,1,:) = potb
       endif      ! done calculating bottom boundary potential 
       
       	       
!       open(unit=10,file='ii.txt')
!!         do j=1,numz
!           do i=1,numr
!              angle=1*2*Pi/numphi
!              write(10,*) i*cos(angle),j,potp(i,j,1) !Vertical cross section
!           enddo
!           write(10,*)
!         enddo
!       close(10)
  
  
!       print*,"Intermediate potential file ii.txt printed"

!       open(unit=10,file='pott',form='unformatted',status='unknown')
!       open(unit=11,file='pots',form='unformatted',status='unknown')
!       write(10) pott
!       write(11) pots
!       close(10)
!       close(11)
!       print*,"mmax bessel",mmax
!       print*,"numphi_by_two",numphi_by_two
       return
       end subroutine bessel
