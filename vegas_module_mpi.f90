! modified 10/06/2016
! INSTRUCTIONS:
! if veff and mass are given in fm^-1 ----> pass htc = 1
! if veff and mass are given in MeV ------> pass htc = 197.3
! call the subroutine set_number_operators(nop) to set the number  of operators  in veff
! (this can be changed....)
! kF must  be in fm^-1
! in the internal functions everything is in fm and fm^-1,  attention!


module vegas_module
  use precision
  use Boundary

  real(kind=my_kind), dimension(:,:), allocatable:: ImSEh, ImSEh_er, ImSEp, ImSEp_er
  real(kind=my_kind), dimension(:,:), allocatable:: ReSEh, ReSEh_er, ReSEp, ReSEp_er
  real(kind=my_kind), dimension(:,:), allocatable :: ff, gg, hh
  real(kind=my_kind), dimension(:), allocatable :: qq
  real(kind=my_kind) :: dq
  integer, parameter :: nq = 100   ! poi vedi 
  integer :: nmax_par
  
contains

  subroutine set_number_operators(nop)
    implicit none
    integer, intent(in)  :: nop   
    nmax_par = nop
    return
  end subroutine set_number_operators

  
  subroutine myvegas(ndim, ncomp, integrand, userdata, nvec,epsrel, epsabs, verbose, seed,&
       mineval, maxeval, nstart, nincrease, nbatch, gridno, statefile,&
       neval, fail, integral, error, prob)
    implicit none
    integer ndim, ncomp, nvec, last, seed, mineval, maxeval
    double precision epsrel, epsabs
    double precision, dimension(1:4) :: userdata
    integer verbose, nregions, fail
    integer nstart, nincrease, nbatch, gridno
    character*(*) statefile
    double precision integral(ncomp), error(ncomp), prob(ncomp)
    integer neval
    integer spin
    parameter (spin = -1)
        
    interface
       integer function integrand(ndim, xx, ncomp, ff,userdata)
         implicit none
         integer ndim, ncomp
         double precision xx(ndim), ff(ncomp)
         double precision, dimension(1:4)::  userdata
       end function integrand
    end interface

    ! cuba 4.1 -----> spin
    call vegas(ndim, ncomp, integrand, userdata, nvec,epsrel, epsabs, verbose, seed,&
         mineval, maxeval, nstart, nincrease, nbatch, gridno, statefile, spin, &
         neval, fail, integral, error, prob)   
    return
  end subroutine myvegas
  

end module vegas_module


!========================================================================================






!========================================================================================
!         REAL PART
!========================================================================================

subroutine real_SE(ReSE,ReSE_er,EE,nE,kk,nk,dove,failed,mass,kF,htc,i_in, i_fin)
  use Boundary
  use vegas_module
  !use precision
  
  implicit none
  integer, parameter :: r8 = selected_real_kind(15,9)
  !.....SET VEGAS param
  integer  :: ndim, ncomp, nvec, last, seed, mineval, maxeval
  double precision :: epsrel, epsabs
  double precision, dimension(1:5) :: userdata
  parameter (ndim = 5)
  parameter (ncomp = 1)
  !parameter (userdata = 0.d0)
  parameter (nvec = 1)
  parameter (epsrel = 1D-3)
  parameter (epsabs = 1D-9)
  parameter (last = 4)
  parameter (seed = 0)
  parameter (mineval = 0)
  parameter (maxeval =  1000000000)  
  integer nstart, nincrease, nbatch, gridno
  character*(*) statefile
  parameter (nstart = 50000)
  parameter (nincrease = 500)
  parameter (nbatch = 1000)
  parameter (gridno = 0)
  parameter (statefile = "")
  
  interface
     integer function integrandREAL(ndim, xx, ncomp, ff,userdata)
       implicit none
       integer ndim, ncomp
       double precision::  xx(ndim), ff(ncomp)
       double precision, dimension(1:4) :: userdata
     end function integrandREAL

     subroutine itot_to_ix_iy(itot, nx, ix, iy )
       implicit none
       integer :: itot, ix, iy, nx
     end subroutine itot_to_ix_iy
     
  end interface
 
  double precision :: integral1(ncomp), error1(ncomp), prob1(ncomp)
  double precision :: integral2(ncomp), error2(ncomp), prob2(ncomp)
  double precision :: integral3(ncomp), error3(ncomp), prob3(ncomp)
  integer :: verbose, nregions, fail1, fail2, fail3
  integer :: neval, neval1, neval2, neval3
  character(len=16) ::  env
  integer :: c
  !....END SET VEGAS

 
  real(kind=r8) :: const
  integer       :: in, ik, nE, nk,i_in, i_fin, ntot, ii
  !integer       :: ik_in, ik_fin, in_in, in_fin
  real(kind=r8) :: kcut, kpole
  !real(kind=r8) :: ReSE(nE,nk), ReSE_er(nE,nk),  EE(nE), kk(nk)
  real(kind=r8) :: ReSE(:), ReSE_er(:),  EE(nE), kk(nk) 
  integer       :: failed(:)
  logical       :: im
  real(kind=r8) :: pi, mass, htc, kF, EkF, Ene
  character(len=64)  :: dove, filenum, nameV




  ntot = size(ReSE)
  
  im      = .false. 
  verbose = 0
  pi = acos(-1.d0)
  const   =  2.d0*pi  * mass / (2.d0*pi)**4
  ! there is no (2\pi) factor in definition of F,G,H (2pi)^6 ----> (2pi)^4

  EkF = (htc *kF)**2.d0/2.d0/mass
  
  kcut    = 12.5d0*kF
  
  
  userdata(3) = mass/htc
  userdata(4) = kF

  
  do ii = i_in, i_fin
     call itot_to_ix_iy(ii, nE, in, ik )
     
     userdata(1) = kk(ik)
     write(filenum,'(I3.3)')  int(kk(ik)/kF * 100)    
     
     if (part) then
        !write(6,*) '---------> RE Particle', kk(ik)/kF, ik
        nameV ='repV'//trim(filenum)
        open(17,file=trim(adjustl(dove))//trim(nameV)//'.dat')
        open(18,file=trim(adjustl(dove))//trim(nameV)//'_1.dat')
        open(19,file=trim(adjustl(dove))//trim(nameV)//'_2.dat')
        open(20,file=trim(adjustl(dove))//trim(nameV)//'_3.dat')
     else
        !write(6,*) '---------> RE Hole    ', kk(ik)/kF, ik
        nameV ='rehV'//trim(filenum)//'.dat'
        open(17,file=trim(adjustl(dove))//nameV)            
     end if
     
     
     Ene = EE(in)
     userdata(2) = EE(in)/htc
     
     kpole = max(kF,kF+sqrt(max (0.d0,kF**2 + 2.d0*mass/htc*Ene/htc) ) )
     
     !----------------> PARTICLES  
     if(part) then
        write(*,*)
        write(*,*) 'REAL particle', kk(ik),  ene, ik, in
        write(*,*)
        lowerIM(1) =  kF
        upperIM(1) =  kpole
        lowerIM(2) =  kF
        upperIM(2) =  kpole
        lowerIM(3) = -1.d0
        upperIM(3) =  1.d0
        lowerIM(4) = -1.d0
        upperIM(4) =  1.d0
        lowerIM(5) =  0.d0
        upperIM(5) =  2.d0*pi
        
        !print *, "---------- Vegas RE pole pole-----------"
        !write(6,*) 'pole pole', upperIM(1),  upperIM(2)
        call myvegas(ndim, ncomp, integrandREAL, userdata, nvec,epsrel, epsabs, verbose, seed,&
             mineval, maxeval, nstart, nincrease, nbatch, gridno, statefile,&
             neval1, fail1, integral1, error1, prob1)
        
        !print *, EE(n), n , EE(n)/EkF
        !print *, "neval    =", neval1
        !print *, "fail     =", fail1
        !print '(F25.12," +- ",F25.12,"   p = ",F8.3)',(integral1(c), error1(c), prob1(c), c = 1, ncomp)                 
        !print *, " "
        
        write(18,155) const*integral1(1),const*error1(1),fail1,kk(ik)/kF,EE(in)/EkF,neval1,kpole
        flush(18)
        
        
        lowerIM(1) =  kF
        upperIM(1) =  kpole
        lowerIM(2) =  kpole
        upperIM(2) =  kcut
        lowerIM(3) = -1.d0
        upperIM(3) =  1.d0
        lowerIM(4) = -1.d0
        upperIM(4) =  1.d0
        lowerIM(5) =  0.d0
        upperIM(5) =  2.d0*pi
        
        !print *, "---------- Vegas RE pole cut-----------"
        !write(6,*) 'pole cut',  upperIM(1),  upperIM(2)
        
        call myvegas(ndim, ncomp, integrandREAL, userdata, nvec,epsrel, epsabs, verbose, seed,&
             mineval, maxeval, nstart, nincrease, nbatch,gridno, statefile, &
             neval2, fail2, integral2, error2, prob2)
        
        !print *, EE(n), n , EE(n)/EkF
        !print *, "neval    =", neval2
        !print *, "fail     =", fail2
        !print '(F25.12," +- ",F25.12,"   p = ",F8.3),(integral2(c), error2(c), prob2(c), c = 1, ncomp)
        
        write(19,155) const*integral2(1),const*error2(1),fail2,kk(ik)/kF,EE(in)/EkF,neval2,kpole
        flush(19)
        
        !print *, " "
        
        lowerIM(1) =  kpole
        upperIM(1) =  kcut
        lowerIM(2) =  kpole
        upperIM(2) =  kcut
        lowerIM(3) = -1.d0
        upperIM(3) =  1.d0
        lowerIM(4) = -1.d0
        upperIM(4) =  1.d0
        lowerIM(5) =  0.d0
        upperIM(5) =  2.d0*pi
        
        !print *, "---------- Vegas RE cut cut -----------"
        !write(6,*) 'cut cut',  upperIM(1),  upperIM(2)
        
        call myvegas(ndim, ncomp, integrandREAL, userdata, nvec,epsrel, epsabs, verbose, seed,&
             mineval, maxeval, nstart, nincrease, nbatch,gridno, statefile, &
             neval3, fail3, integral3, error3, prob3)
        
        !print *, EE(n), n , EE(n)/EkF
        !print *, "neval    =", neval3
        !print *, "fail     =", fail3
        !print '(F25.12," +- ",F25.12,"   p = ",F8.3)',(integral3(c), error3(c), prob3(c), c = 1, ncomp)
        !print *, " "
        write(20,155) const*integral3(1),const*error3(1),fail3,kk(ik)/kF,EE(in)/EkF,neval3,kpole
        flush(20)
        !                                                     
        ReSE(ii) = const*(integral1(1)+2.d0*integral2(1)+integral3(1))
        ReSE_er(ii) = const*sqrt(error1(1)**2+error2(1)**2+error3(1)**2)
        Failed(ii) = fail1 + fail2 + fail3
        neval     = neval1+ neval2 + neval3
        
        ! write(6,156) ReSE(in, ik),ReSE_er(in,ik),Failed(in,ik),kk(ik)/kF,Ene/EkF,neval           
        write(17,156) ReSE(ii),ReSE_er(ii),Failed(ii),kk(ik)/kF,EE(in)/EkF,neval
        flush(17)     
        
        !----------------> HOLES  
     else
        write(*,*)
        write(*,*) 'REAL holes', kk(ik), EE(in),  ik,  in
        write(*,*)
        lowerIM(1) =  0.d0
        upperIM(1) =  kF
        lowerIM(2) =  0.d0
        upperIM(2) =  kF
        lowerIM(3) = -1.d0
        upperIM(3) =  1.d0
        lowerIM(4) = -1.d0
        upperIM(4) =  1.d0
        lowerIM(5) =  0.d0
        upperIM(5) =  2.d0*pi
        
        !print *, "---------- Vegas RE hole-----------"
        call myvegas(ndim, ncomp, integrandREAL, userdata, nvec,epsrel, epsabs, verbose, seed,&
             mineval, maxeval, nstart, nincrease, nbatch,gridno, statefile,&
             neval, fail1, integral1, error1, prob1)
        
        !print *, EE(n), n , EE(n)/EkF
        !print *, "neval    =", neval
        !print *, "fail     =", fail1
        !print '(F25.12," +- ",F25.12,"   p = ",F8.3)',  (integral1(c), error1(c), prob1(c), c = 1, ncomp)
        !print *, " "
        
        ReSE(ii) = const*integral1(1)
        ReSE_er(ii)= const*error1(1)
        Failed(ii)= fail1
        neval = neval1
        write(17,156) ReSE(ii),ReSE_er(ii),Failed(ii),kk(ik)/kF,EE(in)/EkF,neval
        flush(17)
     end if
     
     !end do
     close(17)
     close(18)
     close(19) 
     close(20)
     !     
     
     !call printMC(dove,ReSe,EE,nE,qq,nq,k,kF,EkF,part,im,ReSE_er,failed)
  end do
  
  return
  
155 format(2e18.8e3,i12.1,2e18.8e3,i12.1,e18.8e3) 
156 format(2e18.8e3,i12.1,2e18.8e3,i12.1)     
  
end subroutine real_SE
!===========================================================================




!========================================================================================
! IMAGINARY PART
!========================================================================================

subroutine imaginary_SE(ImSE,ImSE_er,EE,nE,kk,nk,dove,failed,mass,kF,htc,i_in, i_fin)      
  use Boundary
  use vegas_module
  !use precision
  
  implicit none
  integer, parameter :: r8 = selected_real_kind(15,9)
  !.....SET VEGAS param
  integer ndim, ncomp, nvec, last, seed, mineval, maxeval
  double precision epsrel, epsabs
  double precision, dimension(1:4) :: userdata
  parameter (ndim = 4)
  parameter (ncomp = 1)
  ! parameter (userdata = 0.d0)
  parameter (nvec = 1)
  parameter (epsrel = 1D-3)
  parameter (epsabs = 1D-9)
  parameter (last = 4)
  parameter (seed = 0)
  parameter (mineval = 0)
  parameter (maxeval = 1000000000)
  
  integer nstart, nincrease, nbatch, gridno, gridno1, gridno2
  character*(*) statefile
  parameter (nstart = 50000)
  parameter (nincrease = 500)
  parameter (nbatch = 1000)
  parameter (gridno1 = 10)
  parameter (gridno2 = 0)
  parameter (statefile = "")
  
  interface
     integer function integrandIM(ndim, xx, ncomp, ff,userdata)
       implicit none
       integer ndim, ncomp
       double precision xx(ndim), ff(ncomp)
       double precision, dimension(1:4) :: userdata
     end function integrandIM

     subroutine itot_to_ix_iy(itot, nx, ix, iy )
       implicit none
       integer :: itot, ix, iy, nx
     end subroutine itot_to_ix_iy
     
    
  end interface
  
  double precision integral(ncomp), error(ncomp), prob(ncomp)
  integer verbose, nregions, neval, fail
  character(len=16) env
  integer c
  !....END SET VEGAS
  
  real(kind=r8) :: pp_csieta, pmin, pmax
  integer       :: in, ik, nE, nk, i_in, i_fin, ii, ntot
  !integer       :: in_in, in_fin, ik_in, ik_fin
  real(kind=r8) :: ImSE(:),ImSE_er(:), EE(nE), kk(nk)! , the_err(nE)
  real(kind=r8) :: const, mass, pi, kF, htc, EkF, Ene
  integer       :: failed(:)
  logical       ::  im
  character(len=64) ::   dove, name, filenum


  ntot = size(ImSE)
  im = .true.
  verbose = 0
  pi = acos(-1.d0)
  const =  2.d0*pi * pi * mass / (2.d0*pi)**4 ! there is no (2 \pi) factor in definition of FF, GG, HH
  
  EkF = (htc *kF)**2.d0/2.d0/mass
  
  userdata(3) = mass/htc
  userdata(4) = kF
  
  
  write(*,*) 'imaginary part results in ', dove
  
  
  open(7,file='fail.dat')            
  
  do ii = i_in, i_fin
     
     call itot_to_ix_iy(ii, nE, in, ik )     
     userdata(1) = kk(ik)
     
     write(filenum,'(I3.3)')  int(kk(ik)/kF * 100)    
     
     if (part) then
        ! write(6,*) '---------> Im Particle', kk(ik)/kF, ik
        name ='impV'//trim(filenum)
        open(17,file=trim(adjustl(dove))//trim(name)//'.dat')
        !write(*,*) trim(adjustl(dove))//trim(name)
     else
        !write(6,*) '---------> Im Hole    ', kk(ik)/kF, ik
        name ='imhV'//trim(filenum)
        open(17,file=trim(adjustl(dove))//trim(name)//'.dat')
        !write(*,*) trim(adjustl(dove))//trim(name)
     end if

     Ene = EE(in)
     userdata(2) = EE(in)/htc
     pp_csieta   = max(0.d0, 2.d0*mass/htc*Ene/htc)
     pp_csieta   = sqrt(pp_csieta)    
     
     
     
     if(part) then
        write(*,*)
        write(6,*) '---------> IM  Particle', kk(ik), EE(in), ik, in
        write(*,*)
        pmin = kF
        pmax = max(kF,pp_csieta)
        gridno = gridno1
     else
        write(*,*)
        write(6,*) '---------> Im Hole', kk(ik), EE(in), ik, in
        write(*,*)
        pmin = 0.d0
        pmax = kF
        gridno = gridno2
     end if
     ! write(*,*) 'ene/ekf', Ene/EkF, in, pp_csieta, 'pmin, pmax', pmin/kF, pmax/kF
     !     
     !     
     lowerIM(1) =  pmin
     upperIM(1) =  pmax
     lowerIM(2) =  pmin
     upperIM(2) =  pmax
     lowerIM(3) = -1.d0
     upperIM(3) =  1.d0
     lowerIM(4) = -1.d0
     upperIM(4) =  1.d0
     !
     call myvegas(ndim, ncomp, integrandIM, userdata, nvec,  &
          epsrel, epsabs, verbose, seed, mineval, maxeval, &
          nstart, nincrease, nbatch,gridno, statefile,     &
          neval, fail, integral, error, prob)
     !     
     !     print *, EE(in), in 
     !     print *, "neval    =", neval
     !     print *, "fail     =", fail
     !     print '(F25.12," +- ",F25.12,"   p = ",F8.3)',
     !     &           (integral(c), error(c), prob(c), c = 1, ncomp)
     
     !     print *, " "
     IMSE(ii)    = const*integral(1)
     IMSE_er(ii) = const*error(1)
     Failed(ii)   = fail
     write(17,156) ImSE(ii),ImSE_er(ii),Failed(ii),kk(ik)/kF,EE(in)/EkF,neval
     flush(17)     
  end do
  close(17)
  
  !call printMC(dove,ImSe,EE,nE,qq,nq,k,kF,EkF,part,im,ImSE_er,failed)
  
  !  
  
156 format(2e18.8e3,i12.1,2e18.8e3,i12.1)     
  
  return
end subroutine imaginary_SE
!     
!===============================================================================
!



!===========================================================================
!===========================================================================
!===========================================================================
!===========================================================================
!===========================================================================
 
  



integer function integrandIM(ndim, xx, ncomp, ff, userdata)
  use Boundary
  
  implicit none
  integer, parameter :: r8 = selected_real_kind(15,9) 
  integer i,ndim, ncomp
  real(kind=r8) :: xx(*), ff(*)
  real(kind=r8) ::  x(ndim),range(10)
  real(kind=r8) ::  jacobian, ffint_imh, ffint_imp 
  double precision, dimension(1:4) :: userdata
  
  
  jacobian = 1.d0
  do i = 1,ndim
     range(i) = upperIM(i) - lowerIM(i)
     jacobian = jacobian*range(i) 
     x(i) = lowerIM(i) + range(i)*xx(i)
  end do
     
  
  if(part) then
     ff(1) = jacobian*ffint_imp(x(1),x(2),x(3),x(4),userdata)        
  else
     ff(1) = jacobian*ffint_imh(x(1),x(2),x(3),x(4),userdata)        
  end if
  
  integrandIM = 0
  
  return
end function integrandIM



!===============================================================================
!
!     Imaginary part Sigma 2h1p
!
!===============================================================================

real(kind=selected_real_kind(15,9)) function ffint_imh(q1,q2,cq1,cq2,userdata)
  use Boundary
  use vegas_module
  
  implicit none
  integer, parameter :: r8 = selected_real_kind(15,9)
  real(kind=r8) :: q1, q2, cq1, cq2
  real(kind=r8) :: k1, k2, sq1, sq2, fint
  real(kind=r8) :: num, den, sphi0, cphi0
  real(kind=r8) :: ff1, ff3,  Gt, Ft, Gu, Fu, Hu, Ht,  t, u, Mtu
  real(kind=r8) :: sum
  integer :: in, i
  real(kind = r8) :: delta, energy, mass, kF
  double precision, dimension(1:4) :: userdata
  
  ffint_imh = 0.d0
  
  k1     = userdata(1)
  energy = userdata(2)
  mass   = userdata(3)
  kF     = userdata(4)
  
  sq1 = sqrt(1.d0 - cq1**2.d0)
  sq2 = sqrt(1.d0 - cq2**2.d0)
  
  num = -2.d0 * mass * energy - k1**2.d0 +  2.d0*k1*( q1*cq1 + q2*cq2 )
  den =  2.d0*q1*sq1*q2*sq2
  
  
  if(sq1.eq.0.d0 .or. sq2.eq.0.d0 ) then
     if(sq1.eq.0.d0) then
        sphi0 = sq2
     else
        sphi0=  sq1
     end if
     cphi0 = sqrt(1.d0-sphi0**2.d0)
     delta = - num + 2.d0 * q1*q2 * cphi0
     if(delta.ne.0.d0) return
     ff1 = 1.d0
  else
     sphi0 = num/den - (cq1*cq2)/(sq1 *sq2 )
     if (abs(sphi0).gt.1.d0) return
     if (abs(sphi0).ge.1.d0) then
        ff1 = 1.d0
     else
        cphi0 = sqrt( 1.d0 - sphi0**2 )
        ff1 =  2.d0/cphi0/den
     end if
  end if
  
  k2 = sqrt( q1**2.d0 + q2**2.d0 + k1**2.d0 - 2.d0*k1*( q1 * cq1 + q2 * cq2 ) &
       + 2.d0*q1*q2*( sq1*sq2*sphi0 + cq1*cq2 ) )     
  if(k2.lt.kF) return
 

  ff3 = q1**2 * q2**2
   
  t = sqrt( q1**2 + k1**2 - 2.d0*k1*q1*cq1)
  u = sqrt( q2**2 + k1**2 - 2.d0*k1*q2*cq2)

  sum = 0.d0
  
  do in = 1, nmax_par
     ft = fint(qq,ff(:,in),nq,dq,t,1)
     fu = fint(qq,ff(:,in),nq,dq,u,1)
     gt = fint(qq,gg(:,in),nq,dq,t,1)
     gu = fint(qq,gg(:,in),nq,dq,u,1)
     ht = fint(qq,hh(:,in),nq,dq,t,1)
     hu = fint(qq,hh(:,in),nq,dq,u,1)
     sum = sum +  (gt*ft + gu *fu) - ( ht*fu + hu*ft)
  end do

  Mtu = sum

  
  ffint_imh =  ff1*ff3*Mtu
  
  return
end  function ffint_imh
!     
!===============================================================================
! 

!===============================================================================
!
!     Imaginary part Sigma 2p1h
!
!===============================================================================


real(kind=selected_real_kind(15,9)) function ffint_imp(q1,q2,cq1,cq2,userdata)
  use vegas_module
  
  implicit none
  integer, parameter :: r8 = selected_real_kind(15,9)
  real(kind=r8) ::  q1, q2, cq1, cq2
  real(kind=r8) ::  k1, k2, sq1, sq2, fint
  real(kind=r8) ::  num, den, sphi0, cphi0
  real(kind=r8) :: ff1, ff3,  Gt, Ft, Gu, Fu, Hu, Ht,  t, u, Mtu
  real(kind=r8) :: sum
  integer:: in
  real(kind=r8) :: delta, energy, mass, kF
  double precision, dimension(1:4) :: userdata
  
  ffint_imp = 0.d0
  
  k1     = userdata(1)      
  energy = userdata(2)
  mass   = userdata(3)
  kF     = userdata(4)
  
  sq1 = sqrt(1.d0 - cq1**2)
  sq2 = sqrt(1.d0 - cq2**2)
  
  num = -2.d0* mass *energy- k1**2 +  2.d0*k1*( q1*cq1 + q2*cq2 )
  den =  2.d0*q1*sq1*q2*sq2
  
  if(sq1.eq.0.d0 .or. sq2.eq.0.d0 ) then

     if(sq1.eq.0.d0) then
        sphi0 = sq2
     else
        sphi0 =  sq1
     end if
     
     cphi0 = sqrt(1.d0-sphi0**2.d0)
     delta = num - 2.d0 * q1*q2 * cphi0
     if(delta.ne.0.d0) return
     ff1 = 1.d0
  else
     sphi0 = num/den - (cq1*cq2)/(sq1 *sq2 )
     if (abs(sphi0).gt.1.d0) return
     if (abs(sphi0).ge.1.d0)then
        ff1 = 1.d0
     else
        cphi0 = sqrt( 1.d0 - sphi0**2 )
        ff1 =  2.d0/cphi0/den
     end if
  end if
  
  k2 = sqrt(q1**2 + q2**2 + k1**2 - 2.d0*k1*( q1 * cq1 + q2 * cq2 ) &
       + 2.d0*q1*q2*( sq1*sq2*sphi0 + cq1*cq2 ))
  
  if(k2.gt.kF) return
  
  ff3 = q1**2 * q2**2
  
  t = sqrt( q1**2 + k1**2 - 2.d0*k1*q1*cq1)
  u = sqrt( q2**2 + k1**2 - 2.d0*k1*q2*cq2)
  
  sum = 0.d0
  do in = 1, nmax_par
     ft = fint(qq,ff(:,in),nq,dq,t,1)
     fu = fint(qq,ff(:,in),nq,dq,u,1)
     gt = fint(qq,gg(:,in),nq,dq,t,1)
     gu = fint(qq,gg(:,in),nq,dq,u,1)
     ht = fint(qq,hh(:,in),nq,dq,t,1)
     hu = fint(qq,hh(:,in),nq,dq,u,1)

     sum = sum +  (gt*ft + gu *fu) - ( ht*fu + hu*ft)
  end do

  Mtu = sum

  !      const =  2.d0*pi * pi * M / (2.d0*pi)**6
  
  ffint_imp =  ff1*ff3*Mtu
  return
end  function ffint_imp
!     
!===============================================================================
!     


integer function integrandREAL(ndim, xx, ncomp, ff,userdata)
  use Boundary  
 
  implicit none
  integer, parameter :: r8 = selected_real_kind(15,9)
  integer ::  i, ndim, ncomp
  real(kind=r8) ::  xx(ndim), ff(ncomp)
  real(kind=r8) ::  x(ndim),range(ndim)
  double precision, dimension(1:4) :: userdata
  real(kind=r8) ::  jacobian, ffint_realp, ffint_realh
  
  
  
  jacobian = 1.d0
  do i = 1,ndim
     range(i) = upperIM(i) - lowerIM(i)
     jacobian = jacobian*range(i) 
     x(i) = lowerIM(i) + range(i)*xx(i)
  end do
  
  if(part) then
     ff(1) = jacobian*ffint_realp(x(1),x(2),x(3),x(4),x(5),userdata)   
  else
     ff(1) = jacobian*ffint_realh(x(1),x(2),x(3),x(4),x(5),userdata)   
  end if
  
  
  integrandREAL = 0
  
  return
end function integrandREAL


!===============================================================================

!===============================================================================
!  
!     Real part Sigma 2h1p
!
!===============================================================================

real(kind=selected_real_kind(15,9)) function ffint_realh(q1,q2,cq1,cq2,phi,userdata)
  use vegas_module
  
  implicit none
  integer, parameter :: r8 = selected_real_kind(15,9)
  real(kind=r8) ::  q1, q2, cq1, cq2, phi
  real(kind=r8) :: small, fint
  real(kind=r8) :: ff1, ff3, t, u, ft, fu, gt, gu, ht, hu, Mtu
  real(kind=r8) :: sq1, sq2, k1, k2
  real(kind=r8) :: den
  real(kind=r8) :: sum, energy, mass, kF
  integer :: in
  double precision, dimension(1:4) :: userdata
  
  ffint_realh = 0.d0
  small = 1.e-3
  
  k1     = userdata(1)
  energy = userdata(2)
  mass   = userdata(3)
  kF     = userdata(4)
  
  sq1 = sqrt(1.d0 - cq1**2)
  sq2 = sqrt(1.d0 - cq2**2)   
  
  k2 =   sqrt(q1**2 + q2**2 + k1**2 - 2.d0*k1*(q1*cq1 + q2*cq2 ) &
       + 2.d0*q1*q2*(sq1*sq2*sin(phi) + cq1*cq2) )
  
  if (k2.lt.kF) return
  den = 2.d0 * mass * energy -  q1**2 - q2**2 + k2**2 
  ff1 = den /  (den**2 + small**2)
  !write(6,*), 'realh', Eval, q1, q2, k2, den**2, small**2
  !write(6,*), 'realh', den**2, small**2
  
  ff3 = q1**2 * q2**2
  
  t = sqrt( q1**2 + k1**2 - 2.d0*k1*q1*cq1)
  u = sqrt( q2**2 + k1**2 - 2.d0*k1*q2*cq2)
  
  sum = 0.d0
  do in = 1, nmax_par
     ft = fint(qq,ff(:,in),nq,dq,t,1)
     fu = fint(qq,ff(:,in),nq,dq,u,1)
     gt = fint(qq,gg(:,in),nq,dq,t,1)
     gu = fint(qq,gg(:,in),nq,dq,u,1)
     ht = fint(qq,hh(:,in),nq,dq,t,1)
     hu = fint(qq,hh(:,in),nq,dq,u,1)
     
     sum = sum +  (gt*ft + gu *fu) - ( ht*fu + hu*ft)
  end do
  
  Mtu = sum
  
  
  ffint_realh = ff1*ff3*Mtu
  return
end  function ffint_realh

!===============================================================================
!     
!     Real part Sigma 2p1h
!     
!===============================================================================

real(kind=selected_real_kind(15,9)) function ffint_realp(q1,q2,cq1,cq2,phi,userdata)
  use vegas_module
  
  implicit none
  integer, parameter :: r8 = selected_real_kind(15,9)
  real(kind=r8) ::  q1, q2, cq1, cq2, phi
  real(kind=r8) ::  small, fint
  real(kind=r8) ::  ff1, ff3, t, u, ft, fu, gt, gu, ht, hu, Mtu
  real(kind=r8) ::  sq1, sq2, k1, k2
  real(kind=r8) ::  den
  real(kind=r8) ::  sum, energy, kF, mass
  double precision, dimension(1:4) :: userdata
  integer :: in
  
  ffint_realp = 0.d0
  small = 1.e-3

  k1     = userdata(1)
  energy = userdata(2)
  mass   = userdata(3)
  kF     = userdata(4)
  
  !      write(6,*) ' kk nella funzione', k_par, k1
  !     
  sq1 = sqrt(1.d0 - cq1**2)
  sq2 = sqrt(1.d0 - cq2**2)   
  !     
  k2 =   sqrt(q1**2 + q2**2 + k1**2  - 2.d0*k1*(q1*cq1 + q2*cq2) &
       + 2.d0*q1*q2*(sq1*sq2*sin(phi) + cq1*cq2) )
  !     
  if (k2.gt.kF) return
  den =  -2.d0*mass*energy +  q1**2 + q2**2 - k2**2 
  ff1 = den /  (den**2 + small**2)
  !write(6,*), 'realp', Eval, q1, q2, k2, den**2, small**2
  !write(6,*), 'realp',  den**2, small**2
  
  
  ff3 = q1**2 * q2**2
  
  t = sqrt( q1**2 + k1**2 - 2.d0*k1*q1*cq1)
  u = sqrt( q2**2 + k1**2 - 2.d0*k1*q2*cq2)
  
  sum = 0.d0
  do in = 1, nmax_par
     ft = fint(qq,ff(:,in),nq,dq,t,1)
     fu = fint(qq,ff(:,in),nq,dq,u,1)
     gt = fint(qq,gg(:,in),nq,dq,t,1)
     gu = fint(qq,gg(:,in),nq,dq,u,1)
     ht = fint(qq,hh(:,in),nq,dq,t,1)
     hu = fint(qq,hh(:,in),nq,dq,u,1)
     
     sum = sum +  (gt*ft + gu *fu) - ( ht*fu + hu*ft)
  end do
  
  Mtu = sum
  
  
  ffint_realp = - ff1*ff3*Mtu
  return
end  function ffint_realp

!===============================================================================






!==================================================================  
! subroutines to map matrix into array and viceversa
! Mat(ix,ix) <----> vec(itot)
! ix = 1, nx, iy = 1, ny
! itot = (iy-1) * nx + ix
!==================================================================      
  
 
  subroutine matrix_to_vector(mat, nx, ny, vect, ntot )
    implicit none
    integer, parameter :: r8 = selected_real_kind(15,9)
    integer :: nx, ny, ntot
    real(kind=r8), dimension (1:nx,1:ny),intent(in)  :: mat
    real(kind=r8), dimension(1:ntot)    ,intent(out) :: vect
    integer :: ix, iy, itot
    
    if(ntot.eq. nx*ny) then
       write(*,*)  'vactor has wrong dimensions:'
       write(*,*) 'ntot must be = nx*ny=',nx*ny, 'instead of', ntot
       return
    end if
    
    do iy = 1, ny
       do ix = 1, nx
          itot = (iy-1) * nx + ix
          vect(itot) = mat(ix,iy)
       end do
    end do
    
    return
  end subroutine matrix_to_vector
  
  
  subroutine vector_to_matrix(vect, ntot, nx, ny, mat)
    implicit none
    integer, parameter :: r8 = selected_real_kind(15,9)
    integer :: nx, ny, ntot
    real(kind=r8), dimension (1:ntot), intent(in) :: vect
    real(kind=r8), dimension(1:nx,1:ny), intent(out) :: mat
    integer :: ix, iy, itot
    
    do iy = 1, ny
       do ix = 1, nx
          itot = (iy-1) * nx + ix
          mat(ix,iy) = vect(itot)
       end do
    end do
    
    return
    
  end subroutine vector_to_matrix
  
  
!================================================  
  
  subroutine itot_to_ix_iy(itot, nx, ix, iy )
    implicit none
    integer, parameter :: r8 = selected_real_kind(15,9)
    integer :: itot, ix, iy, nx
    real(kind=r8) :: r1
    integer       :: i1, controllo
    
    !r1 = fraction(dble(itot)/dble(nx))
    
    i1 = itot/nx
    r1 = dble(itot)/dble(nx) - dble(i1)
    if(r1.eq.0.d0) then
       ix = nx
       iy = i1
    else  
       ix = nint( dble(nx) * r1)
       iy = i1 + 1
    end if
    controllo  = (iy-1) * nx + ix
    if(controllo.ne.itot) then
       write(*,*) ' something wrong in mapping vector to matrix'
       write(*,*) 'controllo', itot, controllo, i1, r1, nx, ix
       stop
    end if
    
    return
  end subroutine itot_to_ix_iy
  
  
  subroutine ix_iy_to_itot(itot, nx, ix, iy)
    implicit none
    integer, parameter :: r8 = selected_real_kind(15,9)
    integer :: itot, ix, iy, nx
    
    itot = (iy-1) * nx + ix
    return
  end subroutine ix_iy_to_itot


!================================
