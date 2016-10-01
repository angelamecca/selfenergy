! This program computes the single-particle properties of
!  PURE NEUTRON MATTER (1)
!
! using the CBF-effective interacion approach
!
!
program veff_PNM
  !use precision
  use typeveffdata
  use algebra
  use selfenergy_first
  use selfenergy_second
  use selfenergy
 ! use vegas_module
  use timing_module
  use representation
  use utilities  


 implicit none 

  interface
     subroutine deriv5p_new(x,Nx,dx,f,df,n0)
       implicit none
       integer Nx, n0
       double precision x(Nx), f(Nx), df
       double precision dx
     end subroutine deriv5p_new
     
     subroutine deriv5p_all(x, Nx,dx,f,df)
       implicit none
       integer Nx
       double precision x(Nx), f(Nx), df(Nx)
       double precision dx
     end subroutine deriv5p_all
  end interface


 
  integer, parameter :: r8 = selected_real_kind(15,9)   
  type(VEFFdata), dimension(:), allocatable  :: veffmat_in, veffmat
  
  integer :: irho, irhomax, ir,irmax, ipol
  character(len=64):: home, filetable
  
  
  real(kind=r8), dimension(:), allocatable :: kval, ep, eh, k_short, ek_hf, e0, der_ekhf, mstar_hf
  !real(kind=r8), dimension(:), allocatable :: ct
  real(kind=r8) :: e1, etot, tF, e1_check, etot_check, tF_check
  !real(kind=r8) :: dct
  !integer :: nct, ict
  integer :: nk, ik,  nEh, nEp, iE, nk_short
  real(kind=r8) :: kF, rho,  kmin, kmax, dk, kmin_short, kmax_short, dk_short
  real(kind=r8) :: pi, htc, ht2_m,  m
  real(kind=r8) :: ekF, emin, emax, dep, deh
  integer :: nu
  real(kind=r8), dimension(:,:), allocatable :: veff
  real(kind=r8), dimension(:,:), allocatable :: slater
  real(kind=r8), dimension(:),   allocatable :: kF_i, rho_i, x_i, r, y
  integer ::  ikF, iekFh, iekFp
  
  logical :: TSrep
  integer :: polchoice
  real(kind=r8) :: time1, time2
  character(len=64) :: numrho, filename, namepol, polnum
  !integer, parameter :: stencil = 5
  real(kind = r8), dimension(:,:),allocatable :: se_2p1h, se_2h1p
  real(kind = r8), parameter :: small= 1.e-8
  real(kind = r8), dimension(:), allocatable :: lambda
  real(kind = r8), dimension(:),allocatable :: spe, mstar, reh, rep, rehf, imh, imp


    
  TSrep =.false.
  polchoice = 1 ! 1 = PNM, 2 = SYM NUCLEAR, 3 = polarised
  
  if(polchoice.eq.1) then
     nu = 2
     write(*,*) '=====PURE NEUTRON MATTER ======'
     namepol = "pure_ek/"
  else if(polchoice.eq.2) then
     nu = 4
     write(*,*) ' ====== SYM NUCLEAR MATTER ======'
     namepol = "sym_ek/"
  else
     nu = 2
     write(*,*) '======= POLARISED NEUTRON MATTER====='
     allocate(lambda(1:npol))
     open(1,file= 'pol.in')
     do ipol = 1, npol
        read(1,*) lambda(ipol)
        write(*,*) lambda(ipol)
     end do
     write(polnum,'(I3.3)') int(lambda(3)*100.d0 )
     namepol = "pol_ek/"//trim(polnum)//"/"
  end if
  

    




  pi = acos(-1.d0)
  htc = 197.3 !
  m = 939.656 ! neutron mass in  MEV
  ht2_m =  htc*htc/m
  !nct = 100
  !dct = 2.d0/dble(nct)

  !files

  open(2,file=trim(namepol)//'etot.dat')
  write(2,*) '# rho, kF, TF, e1, etot'
  open(10,file=trim(namepol)//'spe_hf.dat')
  write(10,*) '# rho, kF, TF, e1, etot'
  open(16,file=trim(namepol)//'mstarHF.dat')
  write(16,*)  '# rho, kF, k(ikF),  e0(ikF)+ehf(ikF), mstar_hf(ikF)'  
  open(19,file=trim(namepol)//'mstar.in')
  
  home= 'Data/'
  filetable = trim(home)//"veffmat.data"

  call read_table_veff(filetable, veffmat_in, irhomax, nop, nu)
  ! read the effective interacion in 'Data/veffmat.data'
  
  allocate(veffmat(irhomax))
  do irho = 1, irhomax
     if(TSrep) then
        call  change_toV6(veffmat_in(irho),veffmat(irho))
     else
        veffmat(irho) = veffmat_in(irho)
     end if
  end do
  
  ! define matrices A & B -----> N.F. PhD Thesis  pg.
  call set_mat(npol,nop)

  ! set  
  allocate(x_i(npol), kF_i(npol),rho_i(npol))
  !allocate(ct(nct) )
  if(polchoice.eq.1) then
     !x_i(:) = 0.5d0
     x_i(1) = 0.00d0
     x_i(2) = 0.00d0
     x_i(3) = 0.50d0
     x_i(4) = 0.50d0
     write(*,*) ' PURE NEUTRON MATTER x_i : ', x_i(1),x_i(2), x_i(3), x_i(4)
  else if(polchoice.eq.2) then
     x_i(:) = 0.25d0
     write(*,*) ' SYMMETRICAL NUCLEAR MATTER x_i : ', x_i(1), x_i(2), x_i(3), x_i(4)
  else   
     
     do ipol = 1, npol
        x_i(ipol) = lambda(ipol)
     end do
     write(*,*) ' POLARISED NEUTRON MATTER x_i', x_i(1), x_i(2), x_i(3), x_i(4)
     deallocate(lambda)
  end if
        

  !do ict = 1, nct
  !ct(ict) = -1.d0 + dble(ict) * dct
  !end do

  

  do irho =  1 ,2 !irhomax
     
     
     call get_rho(veffmat(irho), rho)

     write(numrho,'(I3.3)') int( rho/0.16d0 * 100.d0)
     filename= trim(namepol)//'spe_hf_'//trim(numrho)//'.dat'
     open(9,file=filename)
     open(15,file=trim(namepol)//'mstar_hf_'//trim(numrho)//'.dat')

     rho_i(:) = rho * x_i(:)              
     kF = (6.d0 * pi**2.d0 * rho / dble(nu))**(1.d0/3.0d0) ! fm^{-1}
     kF_i(:)  = (6.0 * pi**2.d0 * rho_i(:) )**(1.d0/3.0d0) ! fm^{-1} !attention!!!!

     ! kgrid
     kmin =  0.d0 * kF
     kmax =  10.d0 * kF
     dk   = 0.001d0 * kF
     nk   = int ( (kmax-kmin)/ dk  + 0.5d0 )
     ikF  = int( (kF-kmin)/dk + 0.5)
     
     
     kmin_short =   0.d0 * kF
     kmax_short =  2.0d0 * kF
     dk_short   =  0.5d0 * kF
     nk_short  = int ( (kmax_short-kmin_short)/ dk_short  + 0.5d0 ) 
     nk_short = 1
     write(*,*) '*******ATTENZIONE nk_SHORT for second order = 1******'
    
     
     ekF  = ht2_m * kF**2.d0/(2.d0)   ! MeV
     emin = -2.d0 * ekF
     emax =  2.d0 * ekF
     deh  =  0.5d0 * ekF
     dep  =  0.5d0 * ekF
     nEh    = int ( (emax-emin)/deh + 0.5d0)
     nEp    = int ( (emax-emin)/dep + 0.5d0)
     iekFh  = int( (ekF - emin)/dEh + 0.5d0)
     iekFp  = int( (ekF - emin)/dEp + 0.5d0)

     write(*,*) 
     write(*,*) 'rho, kF, EkF, nk,nk_short,  neh nep ', rho, kF, EkF,   nk, nk_short,  nEh, nEp
     write(*,*) 
     
    
     
     
     call get_nx(veffmat(irho), irmax)  
     
     allocate(kval(1:nk),   ek_hf(1:nk), e0(1:nk), der_ekhf(1:nk), mstar_hf(1:nk) )
     allocate(k_short(1:nk_short), spe(1:nk_short), mstar(1:nk_short))
     allocate(ep(nEp),eh(nEh))
     allocate(slater(irmax,npol), r(irmax),y(irmax), veff(irmax,nop))
     !allocate(se_2p1h(nk,nEp), se_2h1p(nk,nEh))
     allocate(rehf(1:nk_short))
     allocate(rep(1:nk_short),reh(1:nk_short),imh(1:nk_short),imp(1:nk_short))

     
     do ik = 1, nk
        kval(ik) = kmin + dble(ik) * dk
        !write(99,*) kval(ik), kval(ik)/kF, kF
     end do
     

    
     do ik = 1 , nk_short
        k_short(ik) = kF !kmin_short +  dble(ik) * dk_short
        !write(*,*) k_short(ik), ik,  nk_short, dk_short
     end do

     do ie = 1, nEp
        ep(ie) = emin + dble(ie) * dep
        !write(*,*) ep(ie), ep(ie)/EkF
     end do
     do ie = 1, nEh
        eh(ie) = emin + dble(ie) * deh
        write(*,*) eh(ie), eh(ie)/EkF
     end do


     write(*,*)
     write(*,*) 'check ikf  : kF  = ',kF, 'k(ikF)   = ', kval(ikF)
     write(*,*) 'check iekF : eKF = ',EkF,  'Ep(iekF) = ', ep(iekFp),  'Eh(iekf) = ', eh(iekFh)
     
     
     
     do ipol = 1, npol
        call get_x(veffmat(irho), y)
        y(:) = kF_i(ipol) * y(:)  
        
        where(y(:).lt.small)
           slater(:,ipol) = 1.d0
        end where
        
        where(y(:) .gt.small)
           slater(:,ipol) = 3.d0 / y(:)**3.d0 * ( sin( y(:) ) - y(:)*cos(y(:) ) ) 
        end where
     end do
     
     

     e0(:) = ht2_m* kval(:)**2.d0/2.d0

     
     call selfe_hf(veffmat(irho), x_i,  slater,  kval, ek_hf)
     call total_ene(veffmat(irho), x_i, rho_i, slater, kval, TF, e1, etot,   ht2_m)
     call check_ene(rho,nu,  nk, kval, ek_hf,tf_check, e1_check, etot_check, ht2_m)

     
     
     
     write(10,*) '# rho, k, e0(k), ehf(k), e0(k)+ehf(k)'
     write(9, *) '# rho, k, e0(k), ehf(k), e0(k)+ehf(k)'
     write(15,*) '# rho, kF,  k, e0(k)+ehf(k), mstar_hf(k)'
    

     write(*,*)
     write(*,*) '======== Total energy ========'
     write(*,*) 'rho, kF, TF, V, etot' 
     write(*,'(1p,6e14.4)')  rho, kF, TF, e1      , etot 
     write(*,'(1p,6e14.4)')  rho, kF, TF, e1_check, etot_check
     write(*,*)'============================'
     write(*,*)
     
     write(2,'(1p,5e14.4)')  rho, kF, TF,       e1      , etot
     write(12,'(1p,5e14.4)') rho, kF, TF_check, e1_check, etot_check
    
     do ik = 1, nk
        write(9,'(1p,5e14.4)' ) rho, kval(ik), e0(ik), ek_hf(ik), e0(ik) + ek_hf(ik)
        write(10,'(1p,5e14.4)') rho, kval(ik), e0(ik), ek_hf(ik), e0(ik) + ek_hf(ik)
     end do

     
    
     close(9)
     
     call timing(time1) 
    ! call selfe_2nd_off(veffmat(irho), x_i, npol, k_short, eh, ep,  nk_short, nEh, nEp, m, htc, namepol)
     call timing(time2)


     call compute_spe(veffmat(irho), x_i, npol, slater,  k_short, spe, mstar, rep, reh, rehf, imh, imp,  nk_short,  m, htc, namepol)

     call compute_nk(veffmat(irho), x_i, npol, k_short, nk_short, ikF,  m, htc, namepol, 1, nk_short)
     write(*,*)
     write(*,*) ' check time : rho = ',rho, 'TIME = ', time2-time1
     write(*,*)


     !effective mass
    ! do ik = 1, nk
        call deriv5p_all(kval,nk, dk,ek_hF,der_ekhf)
      !  Write(*,*) d
     !end do

     mstar_hf(:) =     der_ekhf(:)/kval(:)/ht2_m
     mstar_hf(:) = 1.d0 / ( 1 + mstar_hf(:))

     do ik = 1, nk
        write(15,'(1p,5e14.4)') rho,kF,  kval(ik), e0(ik) + ek_hf(ik), mstar_hf(ik)
     end do


     write(16,'(1p,5e14.4)' ) rho, kF, kval(ikF), e0(ikF) + ek_hf(ikF), mstar_hf(ikF)
     write(19,'(1p,e10.4, e14.4)') rho, mstar_hf(ikF)
     close(15)
     deallocate(kval, k_short,spe, mstar,rep,reh,imh, imp, ek_hf,e0, ep, eh, der_ekhf, mstar_hf)
     deallocate(slater,r, y, veff)
     
     deallocate(rehf)

  end do

  
  deallocate(veffmat_in, veffmat)
  close(2)
  close(10)
  close(16)
  close(19)
end program veff_PNM
