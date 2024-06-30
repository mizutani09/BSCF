
!Gridsize
       integer, parameter :: numr = 120
       integer, parameter :: numz = 120
       integer, parameter :: numphi = 1
       integer, parameter :: hydrophi = 256
       real, parameter :: eps = 1e-3


!Polytropic index
       real, parameter :: np1 = 3.0
       real, parameter :: np2 = 1.7

       real, parameter :: mu1 = 1.3
       real, parameter :: mu2 = 0.5

!Specify Boundary points A and B for a 2d configuration
        integer, parameter :: ix = 12
        integer, parameter :: ax = 120
        integer, parameter :: ay = 2

        integer, parameter :: bx = 2
        integer, parameter :: by = 120

!Not sure what these are
       integer, parameter :: rlwb = 2, rupb = numr - 1

       integer, parameter :: zlwb = 2, zupb = numz - 1

       integer, parameter :: philwb = 1, phiupb = numphi

       integer, parameter :: numphi_by_two = numphi / 2

       real, parameter :: numphiinv = 1.0 / numphi

       real, parameter :: pi = acos(-1.0)
!Symmetry type
       integer, parameter :: isym = 2
