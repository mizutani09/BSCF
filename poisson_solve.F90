subroutine poisson_solve
      implicit none
      include 'runhydro.h'
!************************************************************
!*
!*  Global Variables

      real, dimension(numr,numz,numphi) :: pot, rho
      common /poisson/ pot, rho

!*
!************************************************************      
	
      call setup

      call potsetup

      call potential_solver
 
      return
end subroutine poisson_solve

