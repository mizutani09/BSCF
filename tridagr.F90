      subroutine tridagr(ar,br,cr,knownr,potr)
      implicit none
      include 'runhydro.h'
      include 'pot.h'
!***************************************************************************
!*
!  tridagr solves for potr from the linear tridiagonal system of equations
!  with br being the diagonal elements of the matrix, ar and cr are the
!  off-diagonal elements and knownr is the right hand side.  The code
!  comes from section 2.4, page 43 of Numerical Recipes in Fortran, 2nd ed.
!*
!***************************************************************************
!*
!*  Subroutine Arguments

       real, dimension(numr) :: ar, cr
 
       real, dimension(numr,numphi) :: br

       real, dimension(numr,numz,numphi) :: knownr, potr

!*
!***************************************************************************
!*
!*  Local Variables

       real, dimension(numr,numphi) :: bet, gam  

       integer :: j, k, l

!*
!***************************************************************************
!  initialize the local variables
       gam = 0.0
       bet = 0.0
       j = 0
       k = 0
       l = 0

      ! setup
      do l = 1, numphi
         bet(2,l) = br(2,l)
      enddo
      do l = 1, numphi
         do k = zlwb, zupb
            potr(2,k,l) = knownr(2,k,l) / bet(2,l)
         enddo
      enddo

      !  decomposition and forward substitution
      do l = 1, numphi
         do j = 3, numr-1
            gam(j,l) = cr(j-1) / bet(j-1,l)
            bet(j,l) = br(j,l) - ar(j)*gam(j,l)
            do k = zlwb, zupb
               potr(j,k,l) = (knownr(j,k,l)-ar(j)*potr(j-1,k,l))/ bet(j,l)
            enddo
         enddo
      enddo

      ! back subsitution
      do l = 1, numphi
         do k = zlwb, zupb
            do j = numr-2, 2, -1
               potr(j,k,l) = potr(j,k,l) - gam(j+1,l)*potr(j+1,k,l)
            enddo
         enddo
      enddo

      return
      end subroutine tridagr


