subroutine getinfo(h_0,h_max,rho_2i,count,cput)
!Prints summary of the run on the screen and in a file
  implicit none
  include 'runhydro.h'
  real::rav, mom, m, vol, h_0,am, mac_x, mac_y,  rb, h_max, p_max, cput,  &
        omega, T, W, P, rho_2i, m_core, frac_core, r_core, kappa1, kappa2,   &
        rho_1i, VC, stab, amsq, j2, j4, j6
  character(len=100) :: filename
  character*20 char_np1, char_ax, char_by, char_numr, char_numz, char_m
  character*20 char_vol, char_rav, char_mom, char_h_0, char_am, char_rb,     &
               char_p_max,char_mac_x,char_mac_y,char_cput,char_count,        &
               char_bx, char_T, char_W, char_3P, char_ix, char_np2, char_mu1,&
               char_mu2, char_rcore, char_m_core, char_frac_core,            &
               char_kappa1, char_kappa2, char_rho1i, char_rho2i, char_vc,    &
               char_omega, char_stab
  integer :: count
  double precision :: eggx, eggm, eggr, kt, kc, ke

   if (by==2) then
      rb=-(bx-1.5)/(ax-1.5)
   else 
      rb=(by-1.5)/(ax-1.5)
   endif

  call findmass(rho_2i,m_core,m)
  call findvol(vol)
  call findmom(rho_2i,mom,kt,kc,ke)
  kc = kc/m
  ke = ke/m
  kt = kt/m  


  frac_core = m_core/m

  rav=m/vol
  
  omega=(abs(h_0))**(0.5)
  am=mom*omega
  
  mac_y=h_0/4.0/pi/rav
  mac_x=am**2/(4.0*pi*m**(10.0/3)*rav**(-1/3.0))
  
  p_max=h_max/(1.0+np1)
  r_core=(ix-1.5)*1.0/(ax-1.5)

  call virial(T,W,P,omega,rho_2i)
  VC=abs(2*T-W+3*P)/abs(W)

  rho_1i=rho_2i*mu1/mu2
  
  kappa1 = h_max/(np1+1.0)*1.0**(1+1.0/np1)
  kappa2 = kappa1*(rho_1i)**(1+1.0/np1)/rho_2i**(1+1.0/np2)    
  
  stab = T/abs(W) 

  call findj(m,j2,j4,j6)

 
  !!Convert numbers to strings
  	  
  if (ax.lt.100) then
    write (char_ax, "(I2)") ax
  else
    write (char_ax, "(I3)") ax
  endif

  if ((by.lt.100).and.(by.gt.9)) then
    write (char_by, "(I2)") by
  elseif (by.lt.10) then
    write (char_by, "(I1)") by
  else
    write (char_by, "(I3)") by
  endif

  
  if ((bx.lt.100).and.(bx.gt.10)) then
    write (char_bx, "(I2)") bx
  elseif (bx.lt.10) then
    write (char_bx, "(I1)") bx
  else
    write (char_bx, "(I3)") bx
  endif
  
  write (char_np1, "(F5.1)") np1
  write (char_np2, "(F5.1)") np2
  write (char_numr, "(I3)") numr
  write (char_numz, "(I3)") numz
  write (char_ix, "(I2)") ix
  write (char_mu1,"(F5.1)") mu1
  write (char_mu2,"(F5.1)") mu2
  write (char_rcore,"(F7.4)") r_core
  write (char_m, "(F7.4)") m
  write (char_m_core, "(F7.4)") m_core
  write (char_frac_core, "(F7.4)") frac_core
  write (char_vol, "(F7.4)") vol
  write (char_rav, "(F7.4)") rav
  write (char_mom, "(F7.4)") mom
  write (char_h_0, "(F7.4)") h_0
  write (char_am, "(F7.4)") am
  write (char_rb, "(F7.4)") rb
  write (char_T, "(F7.4)") T
  write (char_W, "(F7.4)") W
  write (char_3P, "(F7.4)") 3*P
  write (char_p_max, "(F7.4)") p_max
  write (char_mac_x, "(F7.4)") mac_x
  write (char_mac_y, "(F7.4)") mac_y
  write (char_count, "(I2)") count
  write (char_cput, "(F8.4)") cput
  write (char_vc, "(F7.4)") VC
  write (char_kappa1, "(F7.4)") kappa1
  write (char_kappa2, "(F7.4)") kappa2
  write (char_rho1i, "(F7.4)") rho_1i  
  write (char_rho2i, "(F7.4)") rho_2i
  write (char_omega, "(F7.4)") omega   
  write (char_stab, "(F7.4)") stab 
 

  if (bx==2) then
!!Make filename
    filename='Bi_'//trim(char_ix)//"_"//trim(char_np1)//'w'//trim(char_mu1)&
    //'_'//trim(char_np2)//'w'//trim(char_mu2)//'_'//trim(char_ax)//"x"//  &
    trim(char_by)//"_"//trim(char_numr)//".info"
  else
    filename='Bi_-'//trim(char_ix)//"_"//trim(char_np1)//'w'//trim(char_mu1)&
    //'_'//trim(char_np2)//'w'//trim(char_mu2)//'_'//trim(char_ax)//"x"//   &
    trim(char_by)//"_"//trim(char_numr)//".info"
  endif

  print*,"================================SUMMARY===================================="
  print*,"n_core  = ", trim(char_np1), "  n_env = ", trim(char_np2)
  print*,"mu_core = ", trim(char_mu1), "  mu_env = ", trim(char_mu2)
  print*,"Equatorial r_core = ", trim(char_rcore)
  print*,"Core mass = ", trim(char_m_core)
  print*,"Core mass fraction = ", trim(char_frac_core)
  print*,"Resolution = ", trim(char_numr),"x", trim(char_numz)
  print*,"kappa_core = ",kappa1, "  kappa_env = ", kappa2
  print*, "Omega", omega
  print*, "Omega_sq", h_0
  print*,"Interface density: ", "Core = ", trim(char_rho1i), "  Env = ", trim(char_rho2i)
  print*,"b/a = ", trim(char_by), "/", trim(char_ax)
  print*,"  rb  ","  Omega_sq  ","  M     ", "  V   ","    J   ","    T   ", "   -W   "&
  ,"  3PI  ","   P_max  "
  print*,trim(char_rb),"   ",trim(char_h_0),"   ",trim(char_m),"  ",trim(char_vol),    &
  "  ",trim(char_am),"  ",trim(char_T),"  ",trim(char_W),"  ",trim(char_3P)," ",       &
  trim(char_p_max)

  print*, "Normalized moment of inertia = ", mom/m
  print*, "J2 = ", J2, "J4 = ", J4
  print*, "J6 = ", J6
  print*,"VC = ", VC
  print*,"cpu time =", trim(char_cput) , " min"

  print*,"==============================OUTPUT FILES================================="  

!!Data is in following format
!!np1 np2 mu1 mu2 ix by ax numr numz r_core rb m_core frac_core Omega_sq M V J T -W 3PI P_max 
 
!np1 np2 mu1 mu2 ix by ax numr numz rcore rb mcore fraccore h0 m vol am T W 3P pmax 
!kappa1 kappa2 rho1i rho2i VC count cput omega 

 
  ! open(unit=10,file=filename)
  ! write(10,*) trim(char_np1)," ",trim(char_np2)," ",trim(char_mu1)," ",trim(char_mu2),     &
  ! " ",trim(char_ix)," ",trim(char_by)," ",trim(char_ax), " ",trim(char_numr)," ",          &
  ! trim(char_numz)," ",trim(char_rcore)," ",trim(char_rb)," ",trim(char_m_core)," ",        &
  ! trim(char_frac_core)," ",trim(char_h_0)," ",trim(char_m)," ",trim(char_vol)," ",         &
  ! trim(char_am), " ",trim(char_T)," ", trim(char_W)," ",trim(char_3P)," ",trim(char_p_max),&
  ! " ",trim(char_kappa1)," ",trim(char_kappa1)," ",trim(char_rho1i)," ",trim(char_rho2i),   &
  ! " ",VC," ", trim(char_count), " ",trim(char_cput), " ", trim(char_omega), " ",           &
  ! trim(char_mac_x), " ", trim(char_mac_y), " ", trim(char_stab)
  ! close(10)

  open(unit=13,file="autoread.dat")
  write(13,*) trim(char_np1)," ",trim(char_np2)," ",trim(char_mu1)," ",trim(char_mu2),     &
  " ",trim(char_ix)," ",trim(char_by)," ",trim(char_ax), " ",trim(char_numr)," ",          &
  trim(char_numz)," ",trim(char_rcore)," ",trim(char_rb)," ",trim(char_m_core)," ",        &
  trim(char_frac_core)," ",trim(char_h_0)," ",trim(char_m)," ",trim(char_vol)," ",         &
  trim(char_am), " ",trim(char_T)," ", trim(char_W)," ",trim(char_3P)," ",trim(char_p_max),&
  " ",trim(char_kappa1)," ",trim(char_kappa2)," ",trim(char_rho1i)," ",trim(char_rho2i),   &
  " ",VC," ", trim(char_count), " ",trim(char_cput), " ", trim(char_omega)," ",            &
  trim(char_mac_x), " ", trim(char_mac_y), " ", trim(char_stab)
  close(13)


  open(unit=14,file="model_details",status="unknown")
  write(14,*) "================================SUMMARY===================================="
  write(14,*) "n_core  = ", trim(char_np1), "  n_env = ", trim(char_np2)
  write(14,*) "mu_core = ", trim(char_mu1), "  mu_env = ", trim(char_mu2)
  write(14,*) "Equatorial r_core = ", trim(char_rcore)
  write(14,*) "Core mass = ", trim(char_m_core)
  write(14,*) "Core mass fraction = ", trim(char_frac_core)
  write(14,*) "Resolution = ", trim(char_numr),"x", trim(char_numz)
  write(14,*) "kappa_core = ",kappa1, "  kappa_env = ", kappa2
  write(14,*) "Omega", omega
  write(14,*) "Omega_sq", h_0
  write(14,*) "Interface density: ", "Core = ", trim(char_rho1i), "  Env = ", trim(char_rho2i)
  write(14,*) "b/a = ", trim(char_by), "/", trim(char_ax), "ix = ", trim(char_ix)
  write(14,*) "rb  ","  Omega_sq  ","  M     ", "  V   ","    J   ","    T   ", "  -W   ",   &
  "  3PI  ","   P_max  "
  write(14,*) trim(char_rb)," ",trim(char_h_0),"   ",trim(char_m),"  ",               &
  trim(char_vol), "  ",trim(char_am),"  ",trim(char_T),"  ",trim(char_W),"  ",trim(char_3P), &
  " ", trim(char_p_max)

  write(14,*) "VC = ", VC
  write(14,*) "Normalized moment of inertia = ", mom/m
  write(14,*) "J2 = ", J2, abs(1.469643E-2 - J2)/1.469643E-2*100
  write(14,*) "J4 = ", J4, abs(-5.8714E-4 - J4)/5.8714E-4*100
  write(14,*) "J6 = ", J6, abs(3.425E-5 - J6)/3.425E-5*100
  write(14,*) "cpu time =", trim(char_cput) , " min"

  write(14,*) "============================================================================"

  open(unit=13,file="hachisu.dat")
  write(13,*) trim(char_rb)," ",trim(char_h_0)," ",trim(char_m),"",trim(char_vol)," ",    &
  trim(char_am), " ",trim(char_T)," ", trim(char_W)," ",trim(char_3P)," ",trim(char_p_max)
  close(13)


  ! open(unit=13,file="maclaurin.dat")
  ! write(13,*) mac_x, " ", mac_y, " ", stab
  ! close(13)

  !   eggx = log10(1/rho_1i)
  !   eggm = frac_core
  !   eggr = log10(1/r_core)

  !   open(unit=12,file='egg.dat',access='APPEND')
  !   write(12,*) eggx, eggm, eggr
  !   close(12)

  ! open(unit=13,file="sc.dat")
  ! write(13,*) r_core, frac_core
  ! close(13)

  ! open(unit=13,file="ru.dat")
  ! write(13,*) r_core, frac_core, m, h_0, kt, kc, ke
  ! close(13)
  
  print*, trim(filename)
  
  end subroutine getinfo
