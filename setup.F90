      subroutine setup
      implicit none
      include 'runhydro.h'
!********************************************************************************
!*
!   setup does a whole bunch of things, highlights are
!   listed below:
!
!   --> read input file, fort.7, and broadcast the values to all pe's
!
!   --> initialize coordinate differentials and arrays for the volumes,
!       face areas and radii for the cell centered grid and vertex
!       centered grid
!
!   --> if an initial model call scfin3d (model_type = 1) which reads in
!       arrays for rho, s, and a.  Then set values for densmin and taumin
!       that will be used throughout the evolution and also initialize
!       the frame number counter to 1
!
!   --> if a continuation model (model_type = 1) read  
!*
!********************************************************************************
!*
!*   Global variables

      real, dimension(numr,numz,numphi) :: pot, rho
      common /poisson/ pot, rho
      
      
      real :: grav
      common /constants/  grav

      real :: dr, dz, dphi, drinv, dzinv, dphiinv
      common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

      real, dimension(numr) :: rhf, r, rhfinv, rinv
      real, dimension(numz) :: zhf
      real, dimension(numphi) :: phi
      common /grid/ rhf, r, rhfinv, rinv, zhf, phi

      real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
      real, dimension(numz) :: zhf_g
      common /global_grid/ rhf_g,r_g,rhfinv_g,rinv_g,zhf_g

!      integer :: isym
!      integer, dimension(3) :: boundary_condition
!      common /boundary_conditions/ isym, boundary_condition
!      common /boundary_conditions/ boundary_condition
      double precision:: radius
!*
!*********************************************************************************
!*
!*  Local variables
	
      integer :: i, j, k, l

      real :: x, xinv, angle

!*
!*********************************************************************************
!  initialize the local variables

       i = 0
       j = 0
       k = 0
       l = 0
       x = 0.0
       xinv = 0.0
       radius=10.0
!  Have root read the fort.7 file into arrays and then broadcast
!  those arrays to all other pe's
!      open(unit=20,file='/Users/kundan/Desktop/Patrick_2D_Poisson/serial/fort.7',form='formatted',status='old')
!      read(20,*) isym
!      read(20,*) boundary_condition(1:3) 
!      close(20)


!  initialize the run parameters
      grav = 1.0
!      write(*,*) 'SETUP: pi ', pi

!  if running with pi or equaorial symmetry then the boundary
!  condition at the bottom of the grid has to be a wall
!  condition
!      if( isym == 2 .or. isym == 3 ) then
!         boundary_condition(1) = 1
!      endif

!  setup coordinates, all pe's use their own portion of the domain
!  decomposition except for the global ( _g ) arrays which span the
!  entire grid

      !  set the coordinate differentials, numr, numz and numphi come
      !  from the runhydro.h header file
!      dr = 1.0 / (numr - 3.0)
       dr = 1.0/ (ax-1.5)    !Freaking oscillations	
      !dr = 1.0 / 127.0
      dz = dr
      dphi = 2.0 * pi * numphiinv
      drinv = 1.0 / dr
      dzinv = 1.0 / dz
      dphiinv = 1.0 / dphi

      !  define r array on every processor, use temp here to avoid 
      !  coercions from integer to floating types
      x = 1.0
      do j = rlwb-1,rupb+1
         r(j) = (x - 2.0)*dr
         x = x + 1.0
      enddo

      !  now define rhf from r
      x = 0.5 * dr 
      do j = rlwb-1,rupb+1
         rhf(j) = r(j) + x
      enddo

      !  and define the inverses of both arrays
      where(r /= 0.0) 
        rinv = 1.0/r
      elsewhere
        rinv = 0.0
      endwhere

      rhfinv = 1.0/rhf 

      ! setup the local zhf array
      x = 1.0
      if( isym /= 1 ) then
         do k = zlwb-1, zupb+1
            zhf(k) = (x-1.5)*dz
            x = x + 1.0
         enddo
      else
         do k = zlwb-1,zupb+1
            zhf(k) = (x - 0.5)*dz
            x = x + 1.0
         enddo 
      endif

      ! set up the azimuthal angle
      x = 0.0
      do l = philwb, phiupb
         phi(l) = x * dphi
         x = x + 1.0
      enddo

      ! global radius array
      x = 1.0
      do j = 1, numr
         r_g(j) = (x-2.0)*dr
         x = x + 1.0
      enddo

      ! global rhf array
      x = 0.5*dr
      do j = 1, numr
         rhf_g(j) = r_g(j) + x
      enddo

      ! define the inverse arrays for r and rhf
      where( r_g /= 0.0 ) 
         rinv_g = 1.0/r_g
      else where
         rinv_g = 0.0
      end where

      rhfinv_g = 1.0/rhf_g

      ! setup the global zhf array
      x = 1.0
      if( isym /= 1 ) then
         do k = 1, numz
            zhf_g(k) = (x-1.5)*dz
            x = x + 1.0
         enddo
      else
         do k = 1, numz
            zhf_g(k) = (x-0.5)*dz
            x = x + 1.0
         enddo
      endif
    
!      open(unit=10,file='r',form='unformatted',status='unknown')
!      write(10) r
!      close(10)
!      open(unit=11,file='rinv',form='unformatted',status='unknown')
!      write(11) rinv
!      close(11)
!      open(unit=12,file='rhf',form='unformatted',status='unknown')
!      write(12) rhf
!      close(12)
!      open(unit=13,file='rhfinv',form='unformatted',status='unknown')
!      write(13) rhfinv
!      close(13)
!      open(unit=14,file='r_g',form='unformatted',status='unknown')
!      write(14) r_g
!      close(14)
!      open(unit=15,file='rinv_g',form='unformatted',status='unknown')
!      write(15) rinv_g
!      close(15)
!      open(unit=16,file='rhf_g',form='unformatted',status='unknown')
!      write(16) rhf_g
!      close(16)
!      open(unit=17,file='rhfinv_g',form='unformatted',status='unknown')
!      write(17) rhfinv_g
!      close(17)
	
	
!      rho = 0.0

!      open(unit=10,file='/Users/kundan/Desktop/Patrick_2D_Poisson/serial/den',form='unformatted',status='unknown')
!      read(10) rho
!      close(10)         
      
!      open(unit=14,file='density_print')

!  do i=1,numr
!    do j=1,numz
!      do k=1,numphi
!        if (rho(i,j,k).gt.0) then
!          print*, rho(i,j,k), i,j,k
!        endif
!      enddo
!    enddo
!  enddo  

!      close(14)
!      print*,'density_print printed'
      
        


  
  
      return
      end subroutine setup
