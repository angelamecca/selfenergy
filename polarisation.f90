module polarisation
  use precision
  implicit none
  
contains

  !===============================================================
  subroutine sum_xAx(Amat,npol, nmax,  xi,  Atot)
    implicit none
    integer ::  npol, nmax
    real(kind=my_kind),dimension(1:npol) :: xi
    real(kind=my_kind),dimension(1:npol,1:npol,1:nmax):: Amat
    real(kind=my_kind),dimension(1:nmax) :: Atot
    integer :: nn, ii, jj
    
    if(size(xi).ne.npol) then
       write(*,*) ' set size(xi)  equal to npol = ', npol
       stop
    end if
    
    if(size(Atot).ne.nmax) then
       write(*,*) ' set size(Atot)', size(Atot), ' equal to nmax = ', nmax
       stop
    end if

    Atot(:) = 0.d0

    do nn = 1, nmax 
       do jj = 1, npol
          do ii = 1, npol
            Atot(nn)   = Atot(nn) + xi(ii) * xi(jj) * Amat(ii,jj,nn)   
          end do
       end do
    end do


    !do nn = 1, nmax
       !write(*,*) ' Atot', Atot(nn), nn
    !end do

    
    return
  end subroutine sum_xAx
!===============================================================
  
  subroutine sum_xA(Amat,npol, nmax,  xi,  Apart)
    implicit none
    integer ::  npol, nmax
    real(kind=my_kind),dimension(1:npol) :: xi
    real(kind=my_kind),dimension(1:npol,1:npol,1:nmax):: Amat
    real(kind=my_kind),dimension(1:nmax) :: Apart
    integer :: nn, ii, jj

    if(size(xi).ne.npol) then
       write(*,*) ' set size(xi)  equal to npol = ', npol
       stop
    end if
    
    if(size(Apart).ne.nmax) then
       write(*,*) ' set size(Apart)', size(Apart), ' equal to nmax = ', nmax
       stop
    end if

    Apart(:) = 0.d0

    do nn = 1, nmax 
       do jj = 1, npol
          do ii = 1, npol
             Apart(nn)   = Apart(nn) + xi(ii)  * Amat(ii,jj,nn)
             !write(*,*)'controllo',ii,xi(ii),jj,   nn, Amat(ii,jj,nn), Apart(nn)
          end do
       end do
    end do    


    
    return
  end subroutine sum_xA

  !=========================================================================
  !=========================================================================
     
  subroutine sum_xBxell(Bmat,npol, nmax, xi, sl, nr,  Btot)
    implicit none
    integer :: nmax, npol
    integer :: nr
    real(kind=my_kind), dimension(1:npol) :: xi
    real(kind=my_kind), dimension(1:nr, 1:npol) :: sl
    real(kind=my_kind), dimension(1:npol,1:npol,1:nmax) :: Bmat
    real(kind=my_kind), dimension(1:nr, 1:nmax) :: Btot
    integer :: nn, ii, jj
  
    if(size(xi).ne.npol) then
       write(*,*) ' set size(xi)  equal to npol = ', npol
       stop
    end if
    
    if(size(Btot(1,:)).ne.nmax) then
       write(*,*) ' set size(Btot(ir,:))',size(Btot(1,:)),'  equal to nmax = ', nmax
       stop
    end if

    Btot(:,:) = 0.d0
    
    do nn = 1, nmax 
       do jj = 1, npol
          do ii = 1, npol
             !do ir = 1, nr
                Btot(:,nn)  = Btot(:,nn) + xi(ii) * xi(jj) * Bmat(ii,jj,nn)* sl(:,jj)   
             !end do
          end do
       end do
    end do
      
    return
  end subroutine sum_xBxell
  
!============
  subroutine sum_xellBxell(Bmat,npol, nmax, xi, sl, nr,  Btot)
    implicit none
    integer :: nmax, npol
    integer :: nr
    real(kind=my_kind), dimension(1:npol) :: xi
    real(kind=my_kind), dimension(1:nr, 1:npol) :: sl
    real(kind=my_kind), dimension(1:npol,1:npol,1:nmax) :: Bmat
    real(kind=my_kind), dimension(1:nr, 1:nmax) :: Btot
    integer :: nn, ii, jj
    
    if(size(xi).ne.npol) then
       write(*,*) ' set size(xi)  equal to npol = ', npol
       stop
    end if
    
    if(size(Btot(1,:)).ne.nmax) then
       write(*,*) ' set size(Btot(ir,:))',size(Btot(1,:)),'  equal to nmax = ', nmax
       stop
    end if

    Btot(:,:) = 0.d0
    
    do nn = 1, nmax 
       do jj = 1, npol
          do ii = 1, npol
             !do ir = 1, nr
                Btot(:,nn)  = Btot(:,nn) + xi(ii) * xi(jj) * Bmat(ii,jj,nn)* sl(:,jj) * sl(:,ii)  
             !end do
          end do
       end do
    end do
      
    return
  end subroutine sum_xellBxell
  !================================================


  subroutine sum_xBell(Bmat,npol, nmax, xi, sl, nr,  Bpart)
    implicit none
    integer :: nmax, npol
    integer :: nr
    real(kind=my_kind), dimension(1:npol) :: xi
    real(kind=my_kind), dimension(1:nr, 1:npol) :: sl
    real(kind=my_kind), dimension(1:npol,1:npol,1:nmax) :: Bmat
    real(kind=my_kind), dimension(1:nr, 1:nmax) :: Bpart
    integer :: nn, ii, jj
    
   
    if(size(xi).ne.npol) then
       write(*,*) ' set size(xi)  equal to npol = ', npol
       stop
    end if  
    
    if(size(Bpart(1,:)).ne.nmax) then
       write(*,*) ' set size(Btot(ir,:))',size(Bpart(1,:)),'  equal to nmax = ', nmax
       stop
    end if

    Bpart(:,:) = 0.d0
    
    do nn = 1, nmax 
       do jj = 1, npol
          do ii = 1, npol
             !do ir = 1, nr
                Bpart(:,nn)  = Bpart(:,nn) + xi(ii)  * Bmat(ii,jj,nn)* sl(:,jj)   
             !end do
          end do
       end do
    end do  
    
    return
  end subroutine sum_xBell
  
  
  !=========================================================================
  !=========================================================================  

  subroutine sumpol_2nd_d(xi,npol, C, nmax, Amat, pol_d1, pol_d2)
    implicit none
    integer :: npol, nmax
    real(kind=my_kind), dimension(1:npol) :: xi
    real(kind=my_kind), dimension(1:nmax,1:nmax,1:nmax) :: C
    real(kind=my_kind), dimension(1:npol,1,npol,1:nmax) :: Amat
    real(kind=my_kind), dimension(1:nmax,1:nmax):: pol_d1, pol_d2
    real(kind=my_kind), dimension(:), allocatable :: Apart    
    integer :: n, m, p
    
    allocate(Apart(1:nmax)) 
 
    pol_d1(:,:) = 0.d0
    pol_d2(:,:) = 0.d0
    
    call sum_xA(Amat,npol, nmax,  xi,  Apart)

    !direct term d1
    do n = 1, nmax
       do m = 1, nmax 
          pol_d1(n,m) = sum(C(n,m,:) * Apart(:))
       end do
    end do
    
    !derect term d2
    do m = 1, nmax
       do n = 1, nmax
          do p = 5, 6
             pol_d2(n,m) = pol_d2(n,m) +  C(n,m,p) * Apart(p)
          end do
       end do
    end do
    
    return
  end subroutine sumpol_2nd_d

  !=======================================================================================

  subroutine sumpol_2nd_e(xi,npol, C, nmax, Amat, pol_e1, pol_e2)
    implicit none
    integer :: npol, nmax
    real(kind=my_kind), dimension(1:npol) :: xi
    real(kind=my_kind), dimension(1:nmax,1:nmax,1:nmax) :: C
    real(kind=my_kind), dimension(1:npol,1,npol,1:nmax) :: Amat
    real(kind=my_kind), dimension(1:nmax,1:nmax):: pol_e1, pol_e2
    real(kind=my_kind), dimension(:,:), allocatable :: coef
    real(kind=my_kind), dimension(:), allocatable :: Apart    
    integer :: n, m, p, q, i
    
    allocate(Apart(1:nmax)) 
    allocate(coef(1:nmax,1:nmax))
 
    pol_e1(:,:) = 0.d0
    pol_e2(:,:) = 0.d0
    coef(:,:)= 0.d0
    call sum_xA(Amat,npol, nmax,  xi,  Apart)
 
    !exchange term e1   
    do p = 1,nmax
       do i = 1, 4
          do n = 1, nmax
             coef(n,p) = coef(n,p) +  C(n,i,p)
          end do
       end do
    end do
    coef(:,:) = 0.25 * coef(:,:)
    
    do m = 1, nmax
       do n = 1,nmax
          do  p = 1, nmax
             do q = 1, nmax
                pol_e1(n,m) = pol_e1(n,m) + coef(n,p)*  C(m,p,q) * Apart(q)
             end do
             do q = 5, 6
                pol_e2(n,m) = pol_e2(n,m) + coef(n,p) * C(m,p,q) * Apart(q)
             end do
          end do
       end do
    end do
    return
  end subroutine sumpol_2nd_e



  
end module polarisation
