module selfenergy
  use precision
  use algebra
  use typeveffdata
  use polarisation
  use integration
  
  real(kind=my_kind), parameter  :: Icos_1 = 2.d0
  real(kind=my_kind), parameter  :: Icos_2 = 0.d0
  real(kind=my_kind), dimension(:,:,:), allocatable :: Amat, Bmat, Bbis

contains

    subroutine set_mat(npol,nmax)
    implicit none
    integer :: npol, nmax
    integer :: ipol, jpol, n
    
    allocate(Amat(npol,npol,nmax), Bmat(npol,npol,nmax))
    allocate(Bbis(npol,npol,nmax))
    
    call create_Amat(Amat,npol,nmax)
    call create_Bmat(Bmat,npol,nmax)
    call create_Bmatbis(Amat, Bbis,npol, nmax)
    open(11,file='mat_A.dat')
    open(12,file='mat_B.dat')
    open(13,file='mat_Bbis.dat')
    open(14,file='mat_Bdiff.dat')
    
    do n = 1, nmax
       write(11,*) 'n = ', n, ' npol = ',  npol
       write(12,*) 'n = ', n, ' npol = ',  npol
       write(13,*) 'n = ', n, ' npol = ',  npol
       write(14,*) 'n = ', n, ' npol = ',  npol
             
       do ipol = 1, npol
          write(11,'(1p,6e12.4)') (Amat(ipol,jpol,n), jpol= 1, npol)
          write(12,'(1p,6e12.4)') (Bmat(ipol,jpol,n), jpol= 1, npol)
          write(13,'(1p,6e12.4)') (Bbis(ipol,jpol,n), jpol= 1, npol)
          write(14,'(1p,6e12.4)') (Bmat(ipol,jpol,n) - Bbis(ipol,jpol,n), jpol= 1, npol)
       end do
    end do

    close(11)
    close(12)
    close(13)
    close(14)
    return
  end subroutine set_mat



   !====================================================================


  subroutine create_Bmatbis(Amat, Bbis,npol, nmax)
    implicit none
    integer :: npol, nmax
    real(kind=my_kind), dimension(1:npol,1:npol,1:nmax) :: Amat, Bbis
    real(kind=my_kind), dimension(:,:,:), allocatable :: C
    real(kind=my_kind), dimension(:,:), allocatable :: coef
    integer :: n, i, p, ipol, jpol
    
    allocate(C(1:nmax,1:nmax,1:nmax))
    allocate(coef(1:nmax, 1:nmax))
    call define_algebra(C,nmax)
    
    Bbis(:,:,:) = 0.d0
    coef(:,:) = 0.d0
    do p = 1,nmax
       do i = 1, 4
          do n = 1, nmax
             coef(n,p) = coef(n,p) +  C(n,i,p)
          end do
       end do
    end do            
    coef(:,:) = 0.25 * coef(:,:)

    

    do n = 1, nmax
       do jpol = 1, npol
          do ipol = 1, npol
             Bbis(ipol,jpol, n) = sum(coef(n,:)* Amat(ipol,jpol,:))
          end do
       end do
    end do
    
    return
  end subroutine create_Bmatbis



 !=======================================================
  
  
  subroutine set_grid(x, nx, xmin, xmax, dx)
    implicit none
    integer :: nx
    real(kind=my_kind), dimension(1:nx) :: x
    real(kind=my_kind) :: xmin, xmax, dx
    integer ix
    
    dx = (xmax - xmin) / dble(nx)
    do ix = 1, nx
       x(ix) = xmin + (dble(ix)- 0.5d0) * dx 
    end do
    return
  end subroutine set_grid
  !==================================================================================



  
!==================================================
!=====================FUNCTIONS====================
!==================================================
  


real(kind=my_kind) function integrate_cos(p,r)
  implicit none
  real(kind=my_kind) :: p, r, t
  real(kind=my_kind), parameter :: small= 1.e-8
  t = p*r
  if(t.lt.small) then
     integrate_cos = 0.d0
  else
     integrate_cos = 4.d0 * (3.d0 * t * cos(t) + (t**2.d0 - 3.d0) * sin(t) ) / t**3.d0
  end if
  return
end function integrate_cos


!============================================

real(kind=my_kind) function bessel(x)
  implicit none
  real(kind=my_kind) :: x
  real(kind=my_kind), parameter :: small= 1.e-8
  
  if(x.lt.small) then
       bessel = 1.d0
    else
       bessel = sin(x)/x
    end if
    
    return
  end function bessel
  !================================================

  
end module selfenergy
