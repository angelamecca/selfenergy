module algebra
  use precision
  implicit none
  integer, parameter :: nop = 6
  integer, parameter :: npol = 4
contains
  
  subroutine create_Amat(Amat,npol,nop)
    implicit none
    integer :: npol, nop
    real(kind=my_kind), dimension(1:npol,1:npol,1:nop) :: Amat 
    
    !A1
    Amat(:,:,1) = +1.0d0
    
    !A2
    Amat(1,1,2)  = +1.0d0
    Amat(1,2,2)  = +1.0d0
    Amat(1,3,2)  = -1.0d0
    Amat(1,4,2)  = -1.0d0
    
    Amat(2,1,2)  = +1.0d0
    Amat(2,2,2)  = +1.0d0
    Amat(2,3,2)  = -1.0d0
    Amat(2,4,2)  = -1.0d0

    Amat(3,1,2)  = -1.0d0
    Amat(3,2,2)  = -1.0d0
    Amat(3,3,2)  = +1.0d0
    Amat(3,4,2)  = +1.0d0

    Amat(4,1,2)  = -1.0d0
    Amat(4,2,2)  = -1.0d0
    Amat(4,3,2)  = +1.0d0
    Amat(4,4,2)  = +1.0d0

    !A3
    Amat(1,1,3)  = +1.0d0
    Amat(1,2,3)  = -1.0d0
    Amat(1,3,3)  = +1.0d0
    Amat(1,4,3)  = -1.0d0
    
    Amat(2,1,3)  = -1.0d0
    Amat(2,2,3)  = +1.0d0
    Amat(2,3,3)  = -1.0d0
    Amat(2,4,3)  = +1.0d0
    
    Amat(3,1,3)  = +1.0d0
    Amat(3,2,3)  = -1.0d0
    Amat(3,3,3)  = +1.0d0
    Amat(3,4,3)  = -1.0d0

    Amat(4,1,3)  = -1.0d0
    Amat(4,2,3)  = +1.0d0
    Amat(4,3,3)  = -1.0d0
    Amat(4,4,3)  = +1.0d0

    !A4
    Amat(1,1,4)  = +1.0d0
    Amat(1,2,4)  = -1.0d0
    Amat(1,3,4)  = -1.0d0
    Amat(1,4,4)  = +1.0d0
    
    Amat(2,1,4)  = -1.0d0
    Amat(2,2,4)  = +1.0d0
    Amat(2,3,4)  = +1.0d0
    Amat(2,4,4)  = -1.0d0

    Amat(3,1,4)  = -1.0d0
    Amat(3,2,4)  = +1.0d0
    Amat(3,3,4)  = +1.0d0
    Amat(3,4,4)  = -1.0d0

    Amat(4,1,4)  = +1.0d0
    Amat(4,2,4)  = -1.0d0
    Amat(4,3,4)  = -1.0d0
    Amat(4,4,4)  = +1.0d0

    ! A5 solo matrice, la dipendenza da cos(theta) nell'integrazione
    Amat(:,:,5) = Amat(:,:,3)

    !A6 solo matrice, la dipendenza da cos(theta) nell'integrazione
    Amat(:,:,6) = Amat(:,:,4)

    return
  end subroutine create_Amat
  
  subroutine create_Bmat(Bmat,npol,nop)
    implicit none
    integer :: nop , npol 
    real(kind=my_kind), dimension(1:npol,1:npol,1:nop) :: Bmat
    integer :: i 
        
    !B1
    Bmat(:,:,:) = 0.0d0
    do i = 1,npol
       Bmat(i,i,1) = 1.0d0
    end do

    !B2
    Bmat(1,1,2)  = +1.0d0
    Bmat(1,2,2)  =  0.0d0
    Bmat(1,3,2)  = +2.0d0
    Bmat(1,4,2)  =  0.0d0
    
    Bmat(2,1,2)  =  0.0d0
    Bmat(2,2,2)  = +1.0d0
    Bmat(2,3,2)  =  0.0d0
    Bmat(2,4,2)  = +2.0d0

    Bmat(3,1,2)  = +2.0d0
    Bmat(3,2,2)  =  0.0d0
    Bmat(3,3,2)  = +1.0d0
    Bmat(3,4,2)  =  0.0d0

   
    Bmat(4,1,2)  =  0.0d0
    Bmat(4,2,2)  = +2.0d0
    Bmat(4,3,2)  =  0.0d0
    Bmat(4,4,2)  = +1.0d0

    !B3
    Bmat(1,1,3)  = +1.0d0
    Bmat(1,2,3)  = +2.0d0
    Bmat(1,3,3)  =  0.0d0
    Bmat(1,4,3)  =  0.0d0
    
    Bmat(2,1,3)  = +2.0d0
    Bmat(2,2,3)  = +1.0d0
    Bmat(2,3,3)  =  0.0d0
    Bmat(2,4,3)  =  0.0d0
    
    Bmat(3,1,3)  =  0.0d0
    Bmat(3,2,3)  =  0.0d0
    Bmat(3,3,3)  = +1.0d0
    Bmat(3,4,3)  = +2.0d0

    Bmat(4,1,3)  =  0.0d0
    Bmat(4,2,3)  =  0.0d0
    Bmat(4,3,3)  = +2.0d0
    Bmat(4,4,3)  = +1.0d0

    !B4
    Bmat(1,1,4)  = +1.0d0
    Bmat(1,2,4)  = +2.0d0
    Bmat(1,3,4)  = +2.0d0
    Bmat(1,4,4)  = +4.0d0
   
    Bmat(2,1,4)  = +2.0d0
    Bmat(2,2,4)  = +1.0d0
    Bmat(2,3,4)  = +4.0d0
    Bmat(2,4,4)  = +2.0d0

    Bmat(3,1,4)  = +2.0d0
    Bmat(3,2,4)  = +4.0d0
    Bmat(3,3,4)  = +1.0d0
    Bmat(3,4,4)  = +2.0d0

    Bmat(4,1,4)  = +4.0d0
    Bmat(4,2,4)  = +2.0d0
    Bmat(4,3,4)  = +2.0d0
    Bmat(4,4,4)  = +1.0d0

    ! B5 solo matrice, la dipendenza da cos(theta) nell'integrazione
    Bmat(:,:,5) = 0.d0
    do i = 1, npol
       Bmat(i,i,5) = 1.d0
    end do
    
    Bmat(1,2,5) = -1.0d0
    Bmat(2,1,5) = -1.0d0
    Bmat(3,4,5) = -1.0d0
    Bmat(4,3,5) = -1.0d0
    
    !B6 solo matrice, la dipendenza da cos(theta) nell'integrazione
    Bmat(1,1,6)  = +1.0d0
    Bmat(1,2,6)  = -1.0d0
    Bmat(1,3,6)  = +2.0d0
    Bmat(1,4,6)  = -2.0d0
   
    Bmat(2,1,6)  = -1.0d0
    Bmat(2,2,6)  = +1.0d0
    Bmat(2,3,6)  = -2.0d0
    Bmat(2,4,6)  = +2.0d0

    Bmat(3,1,6)  = +2.0d0
    Bmat(3,2,6)  = -2.0d0
    Bmat(3,3,6)  = +1.0d0
    Bmat(3,4,6)  = -1.0d0

    Bmat(4,1,6)  = -2.0d0
    Bmat(4,2,6)  = +2.0d0
    Bmat(4,3,6)  = -1.0d0
    Bmat(4,4,6)  = +1.0d0

    return
  end subroutine create_Bmat

  
  
  
  
  subroutine define_algebra(C,nop)
    implicit none
    integer :: nop
    real(kind=my_kind), dimension(1:nop,1:nop,1:nop) :: C  
    integer :: n, m, p    
    
    C(:,:,:) = 0.d0
    ! Cn1_p
    !do n = 1, nop       
    !   C(n,1,n) = 1.d0
    !end do

    C(1,1,1) =  1.d0
    
    ! C^{n2}_p
    C(1,2,2) =  1.d0
    C(2,2,1) =  3.0d0
    C(2,2,2) = -2.0d0
    
    ! C^{n3}_p
    C(1,3,3) =  1.d0
    C(2,3,4) =  1.d0
    C(3,3,1) =  3.d0
    C(3,3,3) = -2.d0

    ! C^{n4}_p
    C(1,4,4) =  1.d0    
    C(2,4,3) =  3.d0
    C(2,4,4) = -2.d0
    C(3,4,2) =  3.d0
    C(3,4,4) = -2.d0
    C(4,4,1) =  9.d0
    C(4,4,2) = -6.d0
    C(4,4,3) = -6.d0
    C(4,4,4) = +4.d0
    
    !C^{n5}_p
    C(1,5,5) =  1.d0
    C(2,5,6) =  1.d0
    C(3,5,5) =  1.d0
    C(4,5,6) =  1.d0
    C(5,5,1) =  6.d0
    C(5,5,3) =  2.d0
    C(5,5,5) = -2.d0
    
    ! C^{n6}_p
    C(1,6,6) =   1.d0
    C(2,6,5) =   3.d0
    C(2,6,6) =  -2.d0
    C(3,6,6) =   1.d0
    C(4,6,5) =   3.d0
    C(4,6,6) =  -2.d0
    C(5,6,2) =   6.d0
    C(5,6,4) =   2.d0
    C(5,6,6) =  -2.d0
    C(6,6,1) =  18.d0
    C(6,6,2) = -12.d0
    C(6,6,3) =   6.d0
    C(6,6,4) =  -4.d0
    C(6,6,5) =  -6.d0
    C(6,6,6) =  +4.d0
    
    
    
    
    do m = 1, nop-1
       do n = m+1, nop
          do p = 1, nop
             C(n,m,p) = C(m,n,p)
             !write(*,*) n,m,p, C(n,m,p), C(m,n,p)
          end do
          !write(*,*) 'def ',n,m,' da ', m, n
          !
       end do
       !
    end do
      
    return
  end subroutine define_algebra
  
  
  
  
end module algebra
