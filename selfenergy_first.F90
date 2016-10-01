module selfenergy_first
  use selfenergy
  use precision
 ! use algebra
 ! use typeveffdata
 ! use polarisation
 ! use integration
  
contains

  ! self energy HARTREE FOCK
  subroutine  selfe_hf(veff, x_i, slater,  kval, ehf)
    implicit none
    
    type(VEFFdata)::  veff
    integer :: nmax  , npol
    real(kind=my_kind), dimension(:)     ::  x_i
    real(kind=my_kind), dimension(:,:)   :: slater
    real(kind=my_kind), dimension(:)     :: kval, ehf
    real(kind=my_kind), dimension(:)    , allocatable :: Atot
    real(kind=my_kind), dimension(:,:)  , allocatable :: Btot
    real(kind=my_kind), dimension(:), allocatable :: vn, integrand, r, vr
    integer :: ir, nr, nn,  kk, nk 
    real(kind=my_kind) ::  pi, dr, rho
    real(kind=my_kind), dimension(:), allocatable :: Icos_3, Icos_4
    pi = acos(-1.d0)
    
    !nmax  veff%nop
    call get_dx(veff ,dr )   !dr   = veff%dr 
    call get_rho(veff,rho)   !rho  = veff%rho
    call get_nx(veff ,nr )   !nr = veff%nr
    call get_nop(veff,nmax) ! NMAX = VEFF%nop
    
    npol = size(x_i(:))
    !nct  = size(ct)
    nk   = size(kval)

    allocate(Atot(nmax),Btot(nr,nmax))

    call sum_xAx(Amat,npol,nmax, x_i,  Atot)
    call sum_xBxell(Bmat,npol, nmax, x_i, slater, nr, Btot)

    allocate(r(nr),integrand(nr), Icos_3(nr), Icos_4(nr))
    allocate(vn(nmax))

    call get_x(veff,r) !r(:) = veff%r
    
    do kk = 1, nk
       do ir = 1, nr
          Icos_3(ir) =  2.d0 * bessel(kval(kk)*r(ir))
          Icos_4(ir) = Integrate_cos( kval(kk),r(ir) )
       end do
       
       do nn = 1, nmax
          integrand(:) = 0.d0
          allocate(vr(nr))
          call get_mat_n(veff,nn,vr)  !vr(:) = veff%mat(:,n)

          if(nn.le.4) then
             integrand(:) = vr(:) * r(:)**2.d0 *  ( Atot(nn) * Icos_1 - Btot(:,nn) * Icos_3(:) )
          else 
             integrand(:) = vr(:) * r(:) **2.d0 * ( Atot(nn) * Icos_2 - Btot(:,nn) * Icos_4(:) )
          end if
          
          vn(nn) = dr*sum(integrand(:))
          
          deallocate(vr)
       end do
              
       ehf(kk) = 2.d0 *pi * rho * sum(vn(:))
       write(3,'(1p,8e14.5)') rho, kval(kk), vn(1),vn(2), vn(3), vn(4), vn(5), vn(6)
    end do
  end subroutine selfe_hf

  !====================================================================

   subroutine  total_ene(veff, x_i, rho_i, slater,  kval, TF, e1, etot, ht2_m)
    implicit none
    
    type(VEFFdata)::  veff
    integer :: nmax  , npol
    real(kind=my_kind), dimension(:)     ::  rho_i, x_i
    real(kind=my_kind), dimension(:), allocatable :: kF_i, TF_i
    real(kind=my_kind), dimension(:,:)   :: slater
    real(kind=my_kind), dimension(:)     :: kval
    real(kind=my_kind) :: tF, e1, etot, ht2_m
    real(kind=my_kind), dimension(:)    , allocatable :: Atot
    real(kind=my_kind), dimension(:,:)  , allocatable :: Btot
    real(kind=my_kind), dimension(:), allocatable :: vn, integrand, r, vr
    integer ::  nr, nn, nk
    real(kind=my_kind) ::  pi, dr, rho
    !real(kind= my_kind) :: kF
    pi = acos(-1.d0)
    
    !nmax  veff%nop
    call get_dx(veff ,dr )   !dr   = veff%dr 
    call get_rho(veff,rho)   !rho  = veff%rho
    call get_nx(veff ,nr )   !nr = veff%nr
    call get_nop(veff,nmax) ! NMAX = VEFF%nop


    npol = size(x_i(:))
    !nct  = size(ct)
    nk   = size(kval)

    allocate (kF_i(1:npol), TF_i(1:npol))
    kF_i(:) = (6.0d0*pi**2.d0 *  rho_i(:))**(1.d0/3.0d0)
    TF_i(:) = ht2_m * 3.d0 / 10.d0 * kF_i(:)**2.d0

    allocate(Atot(nmax),Btot(nr,nmax))

    call sum_xAx(Amat,npol,nmax, x_i,  Atot)
    call sum_xellBxell(Bmat,npol, nmax, x_i, slater, nr, Btot)

    allocate(r(nr),integrand(nr)) 
    allocate(vn(nmax))

    call get_x(veff,r) !r(:) = veff%r
    
    do nn = 1, nmax
       integrand(:) = 0.d0
       allocate(vr(nr))
       call get_mat_n(veff,nn,vr)  !vr(:) = veff%mat(:,n)
       
       if(nn.le.4) then
          integrand(:) = vr(:) * r(:)**2.d0 * ( Atot(nn) - Btot(:,nn)) * Icos_1 
       else 
          integrand(:) = vr(:) * r(:)**2.d0 * ( Atot(nn) - Btot(:,nn)) * Icos_2 
       end if
       vn(nn) = 0.d0    
       vn(nn) = dr*sum(integrand(:))
          
       deallocate(vr)
    end do

    TF  = 0.0d0
    
    TF = sum(x_i(:) *  TF_i(:))
    !TF  =  3./5. * (3.*pi*pi*rho)**(1./3.)
    !write(*,*) 2.d0 * pi * rho/2.d0 * vn(3)
    
    e1 = 2.d0 * pi * rho/2.d0 * sum(vn(:))
    etot = e1 + TF  ! controllare come mettere in modo consistente i fattori 1/2
   
  end subroutine total_ene
       
   !====================================================================

  subroutine check_ene(rho,nu,  nk, kval, ehf,tf, e1, etot, ht2_m)
    implicit none
    integer :: nk, nu
    real(kind=my_kind) :: rho, ht2_m, tf, e1, etot
    real(kind=my_kind), dimension(1:nk):: kval, ehf
    real(kind=my_kind), dimension(1:nk) :: integrand
    real(kind=my_kind) :: kF,dk, pi

    pi = acos(-1.d0)

    dk = kval(2)-kval(1)
    kF = (6.0*pi**2.d0*rho/dble(nu))**(1.0d0/3.0d0)
    integrand(:) = 0.0d0
    
    where(kval(:).le.kF)
       integrand(:) = kval(:)**2.d0 * ehf(:)
    end where


    e1= sum(integrand(:))
    e1 = e1 * dk *dble(nu) / ( 4.d0 * pi**2.d0 * rho)

    tF =  ht2_m * kF**2.d0 * 3.d0/10.d0 

    etot = tF + e1
    
   ! write(*,*)'check',rho, kF,  tF, e1, etot 
    
  end subroutine check_ene

end module selfenergy_first
