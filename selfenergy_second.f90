
module selfenergy_second
  use selfenergy
  use selfenergy_first
  use precision
  use vegas_module
  
  
contains


  !====================================================================================

  subroutine selfe_2nd_off(veff_r, xi, npol, k, EEh, EEp,  nk, nEh, nEp, mass, htc, name)
    implicit none
    integer :: nmax, npol
    real(kind=my_kind) :: htc, mass
    type(VEFFdata) :: veff_r
    real(kind=my_kind), dimension(1:npol) :: xi
    integer :: nk, nEh, nEp
    real(kind=my_kind), dimension(1:nk) ::  k
    real(kind=my_kind), dimension(1:nEh) :: EEh 
    real(kind=my_kind), dimension(1:nEp) :: EEp 
    character(len=64) :: dove, filerho, name
    real(kind=my_kind) :: qmin, qmax, kF, EkF, nyq, qmax_tmp,  hr, pi, rho
    integer :: nu,i
    
    pi =  acos(-1.d0)
    
    
    call get_nu(veff_r,nu)
    call get_rho(veff_r,rho)
    call get_dx( veff_r,hr)
   
    write(filerho,'(I3.3)') int( rho/0.16d0 * 100.d0)
    dove = trim(name)//'selfe/rho_'//trim(filerho)//'/'
    write(*,*) 
    write(*,*) 'Results in  : ', dove
    write(*,*)

    
    kF =  ( 6.d0 * pi**2.d0 * rho / dble(nu) ) ** (1.d0/3.d0)
    EkF =  (htc* htc  / mass) * kF**2.d0 / 2.d0
    nyq =  1.d0/(2.d0*hr)
    qmin = 0.d0
    qmax_tmp = kF * sqrt(EEp(nEp)/EkF)
    write(*,*) '------> qmax ? ', qmax_tmp, 15.d0 *  kF, nyq
    qmax = min(nyq,alpha_cut*kF)

    allocate(qq(1:nq))

    call get_nop(veff_r, nmax)
    
    
    call set_grid(qq, nq, qmin, qmax, dq)
    allocate( ff(1:nq, 1:nmax), gg(1:nq,1:nmax), hh(1:nq, 1:nmax))
    call compute_mel_veff(veff_r, xi, npol, htc, dove,.true.)

    call compute_selfenergy(kF, EkF, nmax, k, EEh, EEp,  nk, nEh, nEp, mass, htc, dove)


    deallocate(qq,ff,gg,hh)
    return
    
  end subroutine selfe_2nd_off

! =========================================================================

  subroutine compute_mel_veff(veff_r, xi, npol, htc, name, iwrite)
      
  ! in this subroutine E and k must be in fm^{-1}
  ! M_par, E_par in fm^{-1}
  ! mass , Eval in MeV

    implicit none
    integer :: nmax, npol
    type(VEFFdata) :: veff_r
    real(kind=my_kind), dimension(1:npol) :: xi
    real(kind=my_kind), dimension(:,:), allocatable :: pol_d1, pol_d2, pol_e1, pol_e2
    real(kind=my_kind) :: pi, htc
    real(kind=my_kind), dimension(:,:,:), allocatable :: C
    integer :: n , m, iq,  nr
    real(kind=my_kind), dimension(:), allocatable :: rr
    real(kind=my_kind), dimension(:,:), allocatable :: vr, vj0, vsl
    real(kind=my_kind) :: sum_d1, sum_d2, sum_e1, sum_e2
    character(len=64) ::   name
    logical :: iwrite
    
    pi =  acos(-1.d0)

    call get_nop(veff_r, nmax)
    call get_nx(  veff_r,nr)

    allocate(rr(1:nr))
    allocate(vr(1:nr,1:nmax))
    
    call get_x(veff_r, rr)
    
    
    allocate(pol_d1(1:nmax,1:nmax),pol_d2(1:nmax,1:nmax))
    allocate(pol_e1(1:nmax,1:nmax),pol_e2(1:nmax,1:nmax))
    allocate(C(1:nmax,1:nmax,1:nmax))
    allocate(vj0(1:nq,1:nmax),vsl(1:nq,1:nmax))
    

     
    call define_algebra(C,nmax)    
    call sumpol_2nd_d(xi,npol, C, nmax, Amat, pol_d1, pol_d2)
    call sumpol_2nd_e(xi,npol, C, nmax, Amat, pol_e1, pol_e2)
    
    write(*,*) 'nq,nmax', nq, nmax
    vj0(:,:) = 0.d0
    vsl(:,:) = 0.d0

    call get_mat(veff_r,vr)
    do n = 1, nmax
       call integration_j0(rr,vr(:,n), qq, vj0(:,n))
       call integration_slater(rr,vr(:,n), qq, vsl(:,n))
    end do

    do m  = 1, nmax
       do iq = 1, nq
          sum_d1 = 0.d0
          sum_d2 = 0.d0
          sum_e1 = 0.d0
          sum_e2 = 0.d0
          do n = 1, nmax
             sum_d1 = sum_d1 + pol_d1(n,m) *  vj0(iq,n)
             sum_d2 = sum_d2 + pol_d2(n,m) * (vj0(iq,n) - 2.d0 * vsl(iq,n))
             sum_e1 = sum_e1 + pol_e1(n,m) *  vj0(iq,n)
             sum_e2 = sum_e2 + pol_e2(n,m) * (vj0(iq,n) - 2.d0 * vsl(iq,n))
          end do
          ! defined in fm^{-1}
          gg(iq,m) = 2.d0 * (sum_d1 + sum_d2)/ htc
          hh(iq,m) = 2.d0 * (sum_e1 + sum_e2)/ htc
          ff(iq,m) = 2.d0 * vj0(iq,m)/htc
       end do
    end do

    if(iwrite) then
       open(19,file=trim(name)//'veff_q1.dat')
       open(20,file=trim(name)//'veff_q2.dat')
       open(21,file=trim(name)//'veff_q3.dat')
       
       write(19,'(1p, 7e20.7)') (qq(iq), gg(iq,1),gg(iq,2), gg(iq,3), gg(iq,4), gg(iq,5),gg(iq,6),  iq = 1,nq)
       write(20,'(1p, 7e20.7)') (qq(iq), hh(iq,1),hh(iq,2), hh(iq,3), hh(iq,4), hh(iq,5),hh(iq,6),  iq = 1,nq)
       write(21,'(1p, 7e20.7)') (qq(iq), ff(iq,1),ff(iq,2), ff(iq,3), ff(iq,4), ff(iq,5),ff(iq,6),  iq = 1,nq)
       
       close(19)
       close(20)
       close(21)
    end if
    
    return
  end subroutine compute_mel_veff

  !==============================

  subroutine compute_selfenergy(kF, EkF, nmax, k, EEh, EEp,  nk, nEh, nEp, mass, htc, name)
    implicit none
    integer :: nmax, nk, nEh, nEp
    real(kind=my_kind) :: kF, EkF, mass, htc
    real(kind=my_kind), dimension(1:nk) :: k
    real(kind=my_kind), dimension(1:nEh) :: EEh 
    real(kind=my_kind), dimension(1:nEp) :: EEp
    integer, dimension(nEp,nk) :: failed_im_p, failed_re_p
    integer, dimension(nEh,nk) :: failed_im_h, failed_re_h
    character(len=64) :: name, filenum, filename
    integer :: ik, iE
    
    allocate(ReSEh(nEh,nk), ReSEh_er(nEh, nk))
    allocate(ReSEp(nEp,nk), ReSEp_er(nEp, nk))
    allocate(ImSEh(nEh,nk), ImSEh_er(nEh, nk))
    allocate(ImSEp(nEp,nk), ImSEp_er(nEp, nk))

    failed_im_p(:,:) = 0
    failed_im_h(:,:) = 0
    failed_re_p(:,:) = 0
    failed_re_h(:,:) = 0

    call set_number_operators(nmax)
    
    part = .false.
    !call imaginary_SE(ImSEh,ImSEh_er,EEh,nEh,k,nk,name,failed_im_h,mass,kF,htc,1,nk,1,n)
    !call real_SE(ReSEh,ReSEh_er,EEh,nEh,k,nk,name,failed_re_h,mass,kF,htc,1,nk,1,nEh)      
    part = .true.
    !call imaginary_SE(ImSEp,ImSEp_er,EEp,nEp,k,nk,name,failed_im_p,mass,kF,htc,1,nk,1,nEp)
    call real_SE(ReSEp,ReSEp_er,EEp,nEp,k,nk,name,failed_re_h,mass,kF,htc,1,nk,1,nEp)
    

    do ik = 1 , nk
       write(filenum,'(I3.3)')  int(k(ik)/kF * 100)
       open(1,file=trim(name)//'imp'//trim(filenum)//'.dat')
       open(2,file=trim(name)//'imh'//trim(filenum)//'.dat')
       open(3,file=trim(name)//'rep'//trim(filenum)//'.dat')
       open(4,file=trim(name)//'reh'//trim(filenum)//'.dat')

       write(1,'(1p,2e18.6e3,1p,i5.1,1p,2e18.6e3)')  &
             (ImSEp(iE,ik),ImSEp_er(iE,ik),failed_im_p(iE,ik),k(ik)/kF,EEp(iE)/EkF, iE=1, nEp)
       write(2,'(1p,2e18.6e3,1p,i5.1,1p,2e18.6e3)')  &
            (ImSEh(iE,ik),ImSEh_er(iE,ik),failed_im_h(iE,ik),k(ik)/kF,EEh(iE)/EkF, iE=1,nEh)
       write(3,'(1p,2e18.6e3,1p,i5.1,1p,2e18.6e3)') &
            (ReSEp(iE,ik),ReSEp_er(iE,ik),failed_re_p(iE,ik),k(ik)/kF,EEp(iE)/EkF, iE=1, nEp)
       write(4,'(1p,2e18.6e3,1p,i5.1,1p,2e18.6e3)')  &
            (ReSEh(iE,ik),ReSEh_er(iE,ik),failed_re_h(iE,ik),k(ik)/kF,EEh(iE)/EkF, iE=1, nEh)

       close(1)
       close(2)
       close(3)
       close(4)
    end do
    
    deallocate(ReSEh, ReSEp, ReSEh_er, ReSEp_er)
    deallocate(ImSEh, ImSEp, ImSEh_er, ImSEp_er)
      
    return
  end subroutine compute_selfenergy

!
! ==============================================================================================================
!
   subroutine compute_spe(veff_r, xi, npol, slater,  k, spe, mstar,  rep, reh, rehf, imh, imp, nk,  mass, htc, dove)
    implicit none
    type(VEFFdata) :: veff_r
    integer :: nk, npol, nmax
    real(kind=my_kind) :: htc, mass
    real(kind=my_kind), dimension(:) ::  k,  xi, spe, mstar
    real(kind=my_kind), dimension(:),allocatable ::  e0, ek
    real(kind=my_kind), dimension(:,:)   :: slater
    real(kind=my_kind), dimension(:) :: rehf, rep, reh, imh, imp
    !real(kind=my_kind), dimension(:), allocatable :: rep_er, reh_er, imh_er, imp_er
    real(kind=my_kind), dimension(:,:),allocatable:: tmh, tmh_er, tmp, tmp_er
    integer, dimension(:,:),allocatable :: fail_imp, fail_rep
    integer, dimension(:,:),allocatable :: fail_imh, fail_reh
    character(len=64) :: dove, name, name1, filerho, filenum, filename
    real(kind=my_kind) :: qmin, qmax, kF, EkF, nyq, qmax_tmp,  hr, pi, rho
    real(kind=my_kind) :: diff, ein
    real(kind=my_kind), parameter :: small = 1.d-6
    integer :: nu, i
    integer :: ik, it
    integer, parameter :: nt = 100
    real(kind=my_kind), dimension(:), allocatable :: kmass, emass, mHF
    
    pi =  acos(-1.d0)
 
    call get_nu(veff_r,nu)
    call get_rho(veff_r,rho)
    call get_dx( veff_r,hr)
   
    write(filerho,'(I3.3)') int( rho/0.16d0 * 100.d0)
    name = trim(dove)//'selfe/rho_'//trim(filerho)//'/'
        
    write(*,*) 
    write(*,*) 'SPE results in  : ', dove
    write(*,*)

    nk = size(k(:))
    allocate(e0(1:nk), ek(1:nk))
   
    
    kF =  ( 6.d0 * pi**2.d0 * rho / dble(nu) ) ** (1.d0/3.d0)
    EkF =  (htc* htc  / mass) * kF**2.d0 / 2.d0
    nyq =  1.d0/(2.d0*hr)
    qmin = 0.d0
    qmax_tmp = k(nk)
    write(*,*) '------', qmax_tmp, 15.d0 *  kF, nyq
    qmax = min(nyq,15.d0*kF)
    
    allocate(qq(1:nq))
    
    
    call set_grid(qq, nq, qmin, qmax, dq)   
    
    call get_nop(veff_r, nmax)
    
    allocate( ff(1:nq, 1:nmax), gg(1:nq,1:nmax), hh(1:nq, 1:nmax))
    call compute_mel_veff(veff_r, xi, npol, htc, dove,.false.)

        
    
    !allocate(rep_er(1:nk),reh_er(1:nk),imh_er(1:nk),imp_er(1:nk))
    allocate(tmh(1:nk,1:nk), tmh_er(1:nk, 1:nk))
    allocate(tmp(1:nk,1:nk), tmp_er(1:nk, 1:nk))
    allocate(fail_reh(1:nk,1:nk),fail_rep(1:nk,1:nk))
    allocate(fail_imh(1:nk,1:nk),fail_imp(1:nk,1:nk))
    allocate(mHF(1:nk), emass(1:nk), kmass(1:nk) )

    tmh(:,:) = 0.d0
    tmp(:,:) = 0.d0
    
    
    call set_number_operators(nmax)    
    e0(:) = (htc* k(:) )**2.d0/(2.d0*mass)

    !write(filerho,'(I3.3)') int( rho/0.16d0 * 100.d0)
    !dove = trim(name)//'selfe/rho_'//trim(filerho)//'/'

    filename = trim(name)//"conv.dat"

    write(*,*) 'conver in ', name, filename
    
    open(19,file=filename)

     
    write(*,*) name, htc, mass
    call selfe_hf(veff_r, xi, slater,  k, rehf)
    do ik = 1, nk
       ek(ik) = e0(ik)
       part = .false.
       !call real_SE(tmh,tmh_er,ek,nk,k,nk,name,fail_reh,mass,kF,htc,ik,ik,ik,ik)
       reh(ik) = tmh(ik,ik)
       part = .true.
       !call real_SE(tmp,tmp_er,ek,nk,k,nk,name,fail_rep,mass,kF,htc,ik,ik,ik,ik)
       rep(ik) = tmp(ik,ik)
       
       ek(ik) = e0(ik) + rehf(ik) + reh(ik) + rep(ik)
       write(*,*) 're', e0(ik), rehf(ik),  reh(ik), rep(ik), ek(ik), ek(ik)/ekF
       do it = 1, nt
          ein = ek(ik)
          part = .false.
          !call real_SE(tmh,tmh_er,ek,nk,k,nk,name,fail_reh,mass,kF,htc,ik,ik,ik,ik)
          reh(ik) = tmh(ik,ik)
          part = .true.
          !call real_SE(tmp,tmp_er,ek,nk,k,nk,name,fail_rep,mass,kF,htc,ik,ik,ik,ik)
          rep(ik) = tmp(ik,ik)
          ek(ik) = e0(ik) + rehf(ik) + reh(ik) + rep(ik)
          diff = abs(ek(ik)-ein)
          if( diff.lt.small) then
             write(*,*)
             write(*,*) 'ho trovato!', ein, ek(ik), e0(ik), rehf(ik), reh(ik), rep(ik)
             write(*,*)
             exit
          else
             write(*,*)
             write(19,*) ' step  ', it, diff, reh(ik), rep(ik)
             write(*,*)
          end if
       end do

       spe(ik) = ek(ik)
    end do

    close(19)

    do ik= 1, nk
       part=.false.
       !call imaginary_SE(tmh,tmh_er,ek,nk,k,nk,name,fail_imh,mass,kF,htc,ik,ik,ik,ik)
       imh(ik)    = tmh(ik,ik)
       part =  .true.
       !call imaginary_SE(tmp,tmp_er,ek,nk,k,nk,name,fail_imh,mass,kF,htc,ik,ik,ik,ik)
       imp(ik)    = tmp(ik, ik)
    end do

    filename = trim(name)//"spe.dat"
    open(7,file=filename)
    write(7,*) '# k, spe, imh, imp, rehf, reh, rep'
    write(7,'(1p,7e20.7)')( k(ik), spe(ik), imh(ik), imp(ik), rehf(ik),reh(ik), rep(ik), ik = 1,nk)
    close(7)



   call compute_mstar(veff_r, xi, npol, slater, rho, kF, EkF, k, spe, nk, kmass, emass, mstar, mHF,  mass, htc, dove, 1, nk)

   
   filename = trim(name)//"mstar_"//trim(filerho)//".dat"
   open(32, file = filename)
   write(32,*) '# k, spe, mstar, mHF'
   write(32,'(1p,5e20.7)') ( k(ik), spe(ik), mstar(ik), mHF(ik), k(ik)/kF,  ik = 1, nk)
   close(32)   
   
   deallocate(e0, ek, mHF,emass, kmass)
   deallocate(tmp, tmh, tmp_er, tmh_er)
   deallocate(fail_reh,fail_rep,fail_imh,fail_imp)
   deallocate(qq,ff,gg,hh)

    
    return
  end subroutine compute_spe

!
!=================================================================
  !
  
  subroutine compute_mstar(veff_r, xi, npol, slater, rho, kF, EkF, k, E, nk, kmass, emass, mstar, mhF,  mass, htc, name, ik_in, ik_fin)
    implicit none

    interface
       subroutine deriv5p_new(x,Nx,dx,f,df,n0)
         implicit none
         integer Nx, n0
         double precision x(Nx), f(Nx), df
         double precision dx
       end subroutine deriv5p_new
    end interface
    
    integer :: nk, nE,  npol, nmax
    integer :: ik_in, ik_fin
    type(VEFFdata) :: veff_r
    real(kind=my_kind), dimension(1:npol) ::   xi
    real(kind=my_kind), dimension(:,:)   :: slater
    real(kind=my_kind) :: htc, mass, kF, EkF, rho
    real(kind=my_kind), dimension(1:nk) :: k, E
    real(kind=my_kind), dimension(:), allocatable :: dk_rehf
    real(kind=my_kind), dimension(:), allocatable :: dk_rep, dk_reh
    real(kind=my_kind), dimension(:), allocatable :: de_rep, de_reh
    real(kind=my_kind), dimension(:) :: kmass, emass, mstar, mHF
    real(kind=my_kind), dimension(:), allocatable :: kmassp, kmassh, kmasshf, emassp, emassh, mstarh, mstarp
    integer, parameter :: imax = 5
    real(kind=my_kind), dimension(1:imax,1:imax) :: tmp, tmp_er, tmh, tmh_er
    real(kind=my_kind), dimension(1:imax) :: eval, kval
    real(kind=my_kind), dimension(1:imax) :: krh, krp, erh, erp, krhf
    !real(kind=my_kind), dimension(1:imax) :: dkrh, dkrp, dkrhf, derh, derp
    real(kind=my_kind) :: dee, dkk
    integer :: i, i0, ik, iE
    character(len=64):: dove, name, filerho    
    integer, dimension(1:imax,1:imax) :: failp,failh
  
    
    dkk = 0.125d0*kF
    dee = 0.250d0*EkF
    !i0 = imax/2 + 1

    write(filerho,'(I3.3)') int( rho/0.16d0 * 100.d0)
    dove = trim(name)//'selfe/rho_'//trim(filerho)//'/'
    
    allocate(dk_rehf(1:nk), dk_rep(1:nk), dk_reh(1:nk))
    allocate(de_rep(1:nk), de_reh(1:nk))
    allocate(kmassp(1:nk), kmassh(1:nk),kmasshf(1:nk), emassp(1:nk), emassh(1:nk))
    allocate(mstarh(1:nk), mstarp(1:nk) ) 
    
    do ik = ik_in, ik_fin
       if (k(ik).lt.dkk) then
          i0 = 1
       else if(k(ik).lt.2.d0*dkk) then
          i0 = 2
       else
          i0 = imax/2 + 1
       end if
       
       do i=1, imax
          kval(i)  = k(ik) + dble(i-i0) * dkk
          eval(i)  = E(ik) + dble(i-i0) * dee
       end do
       
       call selfe_hf(veff_r, xi, slater,  kval, krhf)
       call deriv5p_new(kval,imax,dkk,krhf, dk_rehf(ik),i0)
       
       tmh(:,:) = 0.d0
       tmp(:,:) = 0.d0
       
       part = .false.
       !call real_SE(tmh,tmh_er,eval,imax,kval,imax,dove,failh,mass,kF,htc,i0,i0,1,imax)
       krh(:) = tmh(i0,1:imax)
       !call real_SE(tmh,tmh_er,eval,imax,kval,imax,dove,failh,mass,kF,htc,1,imax,i0,i0)
       erh(:) = tmh(1:imax,i0)
       part = .true.
       !call real_SE(tmp,tmp_er,eval,imax,kval,imax,dove,failp,mass,kF,htc,i0,i0,1,imax)
       krp(:) = tmp(i0,1:imax)
       !call real_SE(tmp,tmp_er,eval,imax,kval,imax,dove,failp,mass,kF,htc,1,imax,i0,i0)
       erp(:) = tmp(1:imax,i0)
       
       call deriv5p_new(eval,imax,dee,erh, de_reh(ik),i0)
       call deriv5p_new(eval,imax,dee,erp, de_rep(ik),i0)
       call deriv5p_new(kval,imax,dkk,krh, dk_reh(ik),i0)
       call deriv5p_new(kval,imax,dkk,krp, dk_rep(ik),i0)
       
    end do
 
       
    kmasshf(ik_in:ik_fin) = 1.0d0/(1.d0 + dk_rehf(ik_in:ik_fin)* mass/k(ik_in:ik_fin) )
    
    emassh(ik_in:ik_fin) = 1.d0 - de_reh(ik_in:ik_fin)
    emassp(ik_in:ik_fin) = 1.d0 - de_rep(ik_in:ik_fin)
    emass (ik_in:ik_fin) = emassh(ik_in:ik_fin) +  emassp(ik_in:ik_fin) - 1.d0
    
    kmassh(ik_in:ik_fin) = 1.0d0/(1.d0 + dk_reh(ik_in:ik_fin)* mass/k(ik_in:ik_fin) )
    kmassp(ik_in:ik_fin) = 1.0d0/(1.d0 + dk_rep(ik_in:ik_fin)* mass/k(ik_in:ik_fin) )
    kmass (ik_in:ik_fin) =  1.d0 / ( (dk_reh(ik_in:ik_fin) + dk_rep(ik_in:ik_fin) + dk_rehf(ik_in:ik_fin) ) * mass/k(ik_in:ik_fin) )
    

    

    mstarh(ik_in:ik_fin) = emassh(ik_in:ik_fin) * kmassh(ik_in:ik_fin)
    mstarp(ik_in:ik_fin) = emassp(ik_in:ik_fin) * kmassp(ik_in:ik_fin)
    mstar(ik_in:ik_fin) =   emass(ik_in:ik_fin) * kmass(ik_in:ik_fin)
    mHF(ik_in:ik_fin)   = kmasshf(ik_in:ik_fin)


    open(9,file=trim(dove)//"der_e.dat")
    open(10,file=trim(dove)//"der_k.dat")
    open(11,file=trim(dove)//"emass.dat")
    open(12,file=trim(dove)//"kmass.dat")
    open(13,file=trim(dove)//"totmass.dat")
    write(9,*)  'k, ek, de_reh, de_rep'
    write(10,*) 'k, ek, dk_rehf, dk_reh, dk_rep' 
    write(11,*) 'k, em_h, em_p, em_tot k/kF'
    write(12,*) 'k, km_hf, km_h, km_p, km_tot k/kF'
    write(13,*) 'k, mstar_hf, mstar_h, mstar_p, mstar, k/kF'
    
    write(9,'(1p,4e20.7)') ( k(ik), E(ik), de_reh(ik), de_rep(ik), ik = ik_in, ik_fin)
    write(10,'(1p,5e20.7)') ( k(ik),E(ik), dk_rehf(ik),dk_reh(ik), de_rep(ik)  ,  ik = ik_in, ik_fin)
    write(11,'(1p,5e20.7)') ( k(ik), emassh(ik) , emassp(ik), emass(ik),  k(ik)/kF           ,  ik = ik_in, ik_fin)
    write(12,'(1p,6e20.7)') ( k(ik), kmasshf(ik), kmassh(ik), kmassp(ik), kmass(ik), k(ik)/kF,  ik = ik_in, ik_fin)
    write(13,'(1p,6e20.7)') ( k(ik), mHF(ik),     mstarh(ik), mstarp(ik), mstar(ik), k(ik)/kF,  ik = ik_in, ik_fin)
    close(9)
    close(10)
    close(11)
    close(12)
    close(13)
    
    deallocate(dk_rehf, dk_rep, dk_reh)
    deallocate(de_rep, de_reh)
    deallocate(kmassp, kmassh, emassp, emassh,mstarh,mstarp)
    

    return
  end subroutine compute_mstar
!
! ==============================================================================================================
!

  subroutine compute_nk(veff_r, xi, npol, k, nk, ikF,  mass, htc, name, ik_in, ik_fin)
    implicit none
    
    interface
       subroutine deriv5p_new(x,Nx,dx,f,df,n0)
         implicit none
         integer Nx, n0
         double precision x(Nx), f(Nx), df
         double precision dx
       end subroutine deriv5p_new
    end interface


    integer :: nmax, npol
    type(VEFFdata) :: veff_r
    real(kind=my_kind), dimension(1:npol) :: xi
    real(kind=my_kind) :: qmin, qmax, nyq, qmax_tmp,  hr, pi, rho
    integer :: nu
    integer :: nk
    integer :: ik_in, ik_fin
    real(kind=my_kind) :: htc, mass, kF, EkF
    real(kind=my_kind), dimension(1:nk) :: k
    real(kind=my_kind), dimension(:), allocatable :: e0, de_rep, de_reh, nmin, nmag
    integer, parameter :: imax = 5
    real(kind=my_kind), dimension(1:imax,1:imax) :: tmp, tmp_er, tmh, tmh_er
    real(kind=my_kind), dimension(1:imax) :: eval, kval
    real(kind=my_kind), dimension(1:imax) :: erh, erp
    real(kind=my_kind) ::  dee, dkk
    integer :: i, i0, ik, ikF
    character(len=64)::  name, dove, filename, filerho  
    integer, dimension(1:imax,1:imax) :: failp,failh
    real(kind = my_kind) :: norm, sum  

    pi = acos(-1.d0)



    call get_nu(veff_r,nu)
    call get_rho(veff_r,rho)
    call get_dx( veff_r,hr)
   
    write(filerho,'(I3.3)') int( rho/0.16d0 * 100.d0)
    dove = trim(name)//'selfe/rho_'//trim(filerho)//'/'
    write(*,*) 
    write(*,*) 'Results in  : ', dove
    write(*,*)


    allocate(e0(1:nk))
    allocate(de_rep(1:nk), de_reh(1:nk))
    allocate(nmin(1:nk), nmag(1:nk))
    
    e0(:) = (htc*k(:))**2.d0/ 2.d0/mass
    

    kF =  ( 6.d0 * pi**2.d0 * rho / dble(nu) ) ** (1.d0/3.d0)
    EkF =  (htc* htc  / mass) * kF**2.d0 / 2.d0
    nyq =  1.d0/(2.d0*hr)
    qmin = 0.d0
    qmax_tmp = kF * sqrt(e0(nk)/EkF)
    write(*,*) '------', qmax_tmp, 15.d0 *  kF, nyq
    qmax = min(nyq,15.d0*kF)

    allocate(qq(1:nq))

    call get_nop(veff_r, nmax)
    
    
    call set_grid(qq, nq, qmin, qmax, dq)
    allocate( ff(1:nq, 1:nmax), gg(1:nq,1:nmax), hh(1:nq, 1:nmax))
    call compute_mel_veff(veff_r, xi, npol, htc, dove,.false.)
    call set_number_operators(nmax)

    
    dee = 0.250d0*EkF
    dkk = 0.125d0 *kF

   
    
    do ik = ik_in, ik_fin
       
       if (k(ik).lt.dkk) then
          i0 = 1
       else if(k(ik).lt.2.d0*dkk) then
          i0 = 2
       else
          i0 = imax/2 + 1
       end if

       
       do i=1, imax
          kval(i)  = k(ik) + dble(i-i0) * dkk
          eval(i)  = e0(ik) + dble(i-i0) * dee
       end do
       
        
       tmh(:,:) = 0.d0
       tmp(:,:) = 0.d0
       
       part = .false.
       call real_SE(tmh,tmh_er,eval,imax,kval,imax,dove,failh,mass,kF,htc,i0,i0,1,imax)
       erh(:) = tmh(i0,1:imax)
       part = .true.
       call real_SE(tmp,tmp_er,eval,imax,kval,imax,dove,failp,mass,kF,htc,i0,i0,1,imax)
       erp(:) = tmp(i0,1:imax)
       
       call deriv5p_new(eval,imax,dee,erh, de_reh(ik),i0)
       call deriv5p_new(eval,imax,dee,erp, de_rep(ik),i0)
       
       
    end do

    where(k(:).le.kF)
       nmin(:) = 1.d0 + de_rep(:)
       nmag(:) = 0.d0
    end where
    where(k(:).ge.kF)
       nmag(:) = - de_reh(:)
       nmin(:) = 0.d0
    end where
    

    sum = 0.d0
    if(ik_fin.gt.ikF) then
       do ik = ik_in, ikF-1
          sum = sum +(  k(ik)**2.d0 * nmin(ik) + k(ik+1)**2.d0 * nmin(ik+1) ) /2.d0
       end do
       do ik = ikF, ik_fin-1
          sum = sum +(  k(ik)**2.d0 * nmag(ik) + k(ik+1)**2.d0 * nmag(ik+1) ) /2.d0
       end do
    end if
    
    norm = 4.d0 * pi * kF**3.d0/3.d0
    sum = 4.d0 * pi * sum/norm
    
    filename= trim(dove)//"nk_"//trim(filerho)//".dat"

    write(*,*)
    write(*,*) ' n(k) in ', filename
    write(*,*)
    
    open(15,file=filename)
    write(15,*) '# k n(k) k/kF'
    write(15,*) '# norm =', sum, ' kmax/kF', k(ik_fin)/kF
    do ik = ik_in, ik_fin
       if(k(ik).le.kF) then
          write(15,'(1p,3e20.7)') k(ik), nmin(ik), k(ik)/kF
       end if
       if(k(ik).ge.kF) then
          write(15,'(1p,3e20.7)') k(ik), nmag(ik), k(ik)/kF
       end if
    end do
    close(15)

    deallocate(qq,ff,gg,hh)
    deallocate(e0, de_rep, de_reh, nmin, nmag)
    return
  end subroutine compute_nk

  

  
 
end module selfenergy_second


!========================================================================================
!========================================================================================
!========================================================================================

