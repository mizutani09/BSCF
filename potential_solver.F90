       subroutine potential_solver
       implicit none
       include 'runhydro.h'
       include 'pot.h'
!*************************************************************
!*
!  potential_solver is the dirver porgram to find the
!  gravitational potential due to the fluid
!* 
!*************************************************************
!*
!*   Global Variables

       real, dimension(numr,numz,numphi) :: pot, rho
       common /poisson/ pot, rho

       real, dimension(numr,numz,numphi) :: potp, rhop
       common /potarrays/ potp, rhop

!      integer :: isym
!      integer, dimension(3) :: boundary_condition
!      common /boundary_conditions/ boundary_condition

!*
!*************************************************************
!*
!*  Local Variables

       integer :: nsteps, l, temp

!*
!*************************************************************
!  initialize the local variable
 !      print*,"Running potential_solver"
       nsteps = 0
       l = 0
       temp = 0

       ! copy the density array into a workspace copy
       rhop = rho
       ! need to zero out sections of rhop to avoid problems
       ! with material piling up in the boundary zones for
       ! the dirichlet boundary conditions
       rhop(:,zupb:zupb+1,:) = 0.0
       rhop(rupb:rupb+1,:,:) = 0.0
       if( isym == 1 ) then
          rhop(:,1,:) = 0.0
       endif

       nsteps = 20
       potp = 0.0
       
       ! solve for the boundary values of the potential
       ! bessel fills in both potp and pot with the
       ! boundary values
       call bessel 

       call helmadi(nsteps)
  
       ! fill in the potential with the solution
       ! if you wanted an external potential in 
       ! addition to self gravity add it in here
       pot(rlwb:rupb,zlwb:zupb,philwb:phiupb) =potp(rlwb:rupb,zlwb:zupb,philwb:phiupb)

       if( isym /= 1 ) then
          pot(:,1,:) = pot(:,2,:)
       endif

       if( isym == 3 ) then
          pot(1,:,:) = pot(2,:,:)
       else
          do l = 1, numphi_by_two
             pot(1,:,l) = pot(2,:,l+numphi_by_two) 
             pot(1,:,l+numphi_by_two) = pot(2,:,l)
          enddo
       endif

       return
       end subroutine potential_solver
