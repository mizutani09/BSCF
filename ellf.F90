      FUNCTION ellf(ak)
      implicit none
!  From section 6.11, page 206 of Numerical Recipes
!  in Fortran 2nd ed.  ellf returns the Legendre 
!  elliptic integral of the first kind, F(phi,kappa)
!  using Carlson's function Rf.  See source code for rf.f.
!  For the hydrocode implementaion only:
!  All calls to ellf have phi = pi / 2, therefor don't pass
!  phi and simplify the calculations in ellf to take advantage
!  of this fact
      REAL ellf,ak
!CU    USES rf
      REAL rf
      ellf=rf(0.0,(1.-ak)*(1.+ak),1.)
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software .
