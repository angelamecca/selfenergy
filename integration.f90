!=======================================================================
!
!     
!=======================================================================
!     

module integration
  implicit none
  integer, parameter, private :: r8 = selected_real_kind(15,9)
  
contains
  !=======================================================================
  
  !\int r^2* v(r) * j0(qr)
  subroutine integration_j0( rr,vr, qq,  vj0)
    implicit none   
    real(kind=r8), dimension(:) :: rr, vr, qq,  vj0
    real(kind=r8), dimension(:), allocatable :: fr
    integer :: nr, nq, iq
    real(kind=r8) :: q, pi
    
    !opt= .true.
    pi = acos(-1.d0)
    nr = size(rr)
    nq = size(qq)
    
    allocate(fr(nr))
    
    fr(:) = rr(:) * vr(:)
    do iq = 1, nq
       q = qq(iq)
       call filon_formula(rr, fr, q, .true., vj0(iq))
       vj0(iq) =  vj0(iq) / q
    end do

   
    deallocate(fr)
    return
  end subroutine integration_j0
  
  !=======================================================================
  ! \int v(r) r^2 \ell(qr)
  subroutine integration_slater( rr,vr, qq,  vsl)
    implicit none
    real(kind=r8), dimension(:) :: rr, vr, qq,  vsl
    real(kind=r8), dimension(:), allocatable :: fr_s, fr_c, vsl_s, vsl_c
    integer :: nr, nq, iq
    real(kind=r8) :: h, q, pi
    
    !opt= .true.
    pi = acos(-1.d0)
    nr = size(rr)
    nq = size(qq)
    h  = rr(2) - rr(1)
    
    allocate(fr_s(nr), fr_c(nr))
    allocate(vsl_s(nq),vsl_c(nq))
    
    fr_s(:) =  vr(:) / rr(:)
    fr_c(:) =  vr(:)
    
    do iq = 1, nq
       q = qq(iq)
       call filon_formula(rr, fr_s, q, .true., vsl_s(iq))
       call filon_formula(rr, fr_c, q, .false., vsl_c(iq))
       vsl(iq) = 3.0 * (  vsl_s(iq)/ q**3.d0  - vsl_c(iq)/q**2.d0)
       !write(19,*) qq(iq),  3.d0* vsl_s(iq)/q**3.d0, 3.d0 *vsl_c(iq)/q**2.d0, vsl(iq)
    end do

    deallocate(fr_s,fr_c,vsl_s,vsl_c)
    return    
  end subroutine integration_slater



  
  !=======================================================================
  ! FILON'S Integration formula
  ! \int f(x) sin(kx) dx       if  opt = true
  ! \int f(x) cos(kx) dx       if  opt = false
  !=======================================================================
  subroutine filon_formula(rr,fr, q, opt, Iq)
    implicit none
    logical opt
    real(kind(1.d0)), dimension(:) :: rr, fr !, qq, fq
    real(kind(1.d0)) :: Iq, q
    integer(kind=4) nr
    real(kind(1.d0)) :: h, theta, even, odd
    real(kind(1.d0)) :: r_1, r_n
    real(kind(1.d0)) :: alpha, beta, gamma, snt, cst, sn2t
    real(kind(1.d0)) :: one, two, three, pi
    integer ::  j
    
    pi = acos(-1.d0)
    nr = size(fr)
    h = abs(rr(2) -rr(1))
    
    r_1 = rr(1)
    r_n = rr(nr)
    
    theta = h * q
    snt = sin(theta)
    sn2t = sin(2.d0*theta)
    cst = cos(theta)
    alpha = 1.d0 / theta + sn2t/ 2.d0/theta**2 - 2.d0 * snt*snt / theta**3
    beta  = 2.d0 *  ( (1.d0 + cst * cst)/theta**2  - sn2t/theta**3 )
    gamma = 4.d0 * ( snt / theta**3  - cst /theta**2)
    
    
    even = 0.d0
    odd  = 0.d0
    
    if(opt) then ! sin
       do j = 1, nr/2
          even = even + fr(2*j)   * sin(q*rr(2*j))
          odd  = odd  + fr(2*j-1) * sin(q*rr(2*j-1))
       end do
       even = even - 1.d0/2.d0 * ( fr(nr)* sin(q*r_n) + fr(1) * sin(q*r_1) )
       
       one   = alpha * (fr(1) * cos(q *r_1) - fr(nr) * cos(q * r_n) ) 
       two   = beta * even 
       three = gamma * odd
       
    else
       do j = 1, nr/2
          even = even + fr(2*j)   * cos(q*rr(2*j))
          odd  = odd  + fr(2*j-1) * cos(q*rr(2*j-1))
       end do
       even = even - 1.d0/2.d0 * ( fr(nr)* cos(q*r_n) + fr(1) * cos(q*r_1) )
       
       one   = alpha * (  fr(nr) * sin(q * r_n) - fr(1) * sin(q *r_1) ) 
       two   = beta * even
       three = gamma * odd
    end if
    
    Iq =   h * (one + two + three ) 
    return
  end subroutine filon_formula


  !=======================================================================
  subroutine FilonFourier(rr,fr,qq,fq,opt)
    implicit none
    logical opt
    !  real(kind(1.d0)), dimension(:) :: rr(1:nr), fr(1:nr), qq(1:nq), fq(1:nq) 
    real(kind(1.d0)), dimension(:) :: rr, fr, qq, fq 
    integer(kind=4) nr, nq
    real(kind(1.d0)), dimension(:), allocatable :: gr
    real(kind(1.d0)) :: h, theta, even, odd, q, r_1, r_n
    real(kind(1.d0)) :: alpha, beta, gamma, snt, cst, sn2t
    real(kind(1.d0)) :: one, two, three, pi
    integer ::  i, j
    
    pi = acos(-1.d0)
    nr = size(fr)
    nq = size(fq)
    h = abs(rr(2) -rr(1))
    
    allocate(gr(nr))

    gr = rr(:) *  fr(:)
    r_1 = rr(1)
    r_n = rr(nr)
    
    do i = 1,nq
  
       q = qq(i)
       theta = h * q
       snt = sin(theta)
       sn2t = sin(2.d0*theta)
       cst = cos(theta)
       alpha = 1.d0 / theta + sn2t/ 2.d0/theta**2 - 2.d0 * snt*snt / theta**3
       beta  = 2.d0 *  ( (1.d0 + cst * cst)/theta**2  - sn2t/theta/theta**2 )
       gamma = 4.d0 * ( snt / theta**3  - cst /theta**2)
       even = 0.d0
       odd  = 0.d0
       do j = 1, nr/2
          even = even + gr(2*j)   * sin(q*rr(2*j))
          odd  = odd  + gr(2*j-1) * sin(q*rr(2*j-1))
       end do
       even = even - 1.d0/2.d0 * ( gr(nr)* sin(q*r_n) + gr(1) * sin(q*r_1) )
       
       one   = alpha * (gr(1) * cos(q *r_1) - gr(nr) * cos(q * r_n) ) 
       two   = beta * even 
       three = gamma * odd
       fq(i) =   (one + two + three ) / q
    end do
    
    if(opt) then
       fq = 4.d0 * pi * h * fq
    else 
       fq = 4.d0 * pi * h * fq/(2.d0*pi)**3
    end if
    
    deallocate(gr)
    return
    
  end subroutine FilonFourier
  
end module integration
