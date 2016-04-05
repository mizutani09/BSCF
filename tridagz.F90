      subroutine tridagz(az,bz,cz,knownz,potz)
      implicit none
      include 'runhydro.h'
      include 'pot.h'
!***************************************************************************
!*
!  tridagz solves for potz from the linear tridiagonal system of equations
!  with bz being the diagonal elements of the matrix, az and cz are the
!  off-diagonal elements and knownz is the right hand side.  The code
!  comes from section 2.4, page 43 of Numerical Recipes in Fortran, 2nd ed.
!*
!***************************************************************************
!*
!*  Subroutine Arguments

       real :: az, cz
 
       real, dimension(numz) :: bz

       real, dimension(numr,numz,numphi) :: knownz, potz

       logical :: iam_on_top

!*
!***************************************************************************
!*
!*  Local Variables

       real, dimension(numz) :: bet, gam  

       integer :: j, k, l, rad_lwr_bnd, rad_upr_bnd

!*
!***************************************************************************
!  initialize the local vraiables
       bet = 0.0
       gam = 0.0
       j = 0
       k = 0
       l = 0
      rad_lwr_bnd = 2
      rad_upr_bnd = rupb
!      if( iam_on_top ) then
!         rad_upr_bnd = numr_dd_pad - 1 - pad
!      else
!         rad_upr_bnd = numr_dd_pad - 1
!      endif

      ! setup
      bet(2) = bz(2)
      do l = 1, numphi
         do j = rad_lwr_bnd, rad_upr_bnd
            potz(j,2,l) = knownz(j,2,l) / bet(2)
         enddo
      enddo

      !  decomposition and forward substitution
      do k = 3, numz-1
         gam(k) = cz / bet(k-1)
         bet(k) = bz(k) - az*gam(k)
         do l = 1, numphi
            do j = rad_lwr_bnd, rad_upr_bnd
               potz(j,k,l) = (knownz(j,k,l) - az*potz(j,k-1,l))/  bet(k)
            enddo
         enddo
      enddo

      ! back subsitution
      do l = 1, numphi
         do k = numz-2, 2, -1
            do j = rad_lwr_bnd, rad_upr_bnd
               potz(j,k,l) = potz(j,k,l) - gam(k+1)*potz(j,k+1,l)
            enddo
         enddo
      enddo

      return
      end subroutine tridagz


