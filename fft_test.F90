      program  fft_test
      implicit none
      include 'runhydro.h'

      real, dimension(numr,numz,numphi) :: rho


! read in the density array
      rho = 0.0
      open(unit=10,file='density_128',form='unformatted',status='unknown')
      read(10) rho
      close(10)

      call realft(rho,numr,numz,numphi,+1)

      open(unit=11,file='density_ftd',form='unformatted',status='unknown')
      write(11) rho
      close(11)

      call realft(rho,numr,numz,numphi,-1)

      open(unit=12,file='density_inv_ftd',form='unformatted', status='unknown')
      write(12) rho
      close(12)

      stop
      end program fft_test
