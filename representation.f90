module representation
  use typeveffdata
  implicit none
  integer, parameter, private :: r8 = selected_real_kind(15,9)
  
contains
  
  subroutine change_toTS(veff_v6,veff_TS)
    implicit none
    type(VEFFdata), intent(in) :: veff_v6
    type(VEFFdata), intent(out) :: veff_TS
    real(kind=r8) :: rho, dr
    real(kind=r8), dimension(:,:), allocatable :: vTS, v6
    real(kind=r8), dimension(:), allocatable :: r
    integer :: nr, nop, nu
    integer :: ir, i, j
    real(kind=r8), dimension(1:4,1:4) :: M
    real(kind=r8), dimension(1:2,1:2) :: Mten   

    call get_rho(veff_v6,rho)
    call get_nu(veff_v6, nu)
    call get_dx(veff_v6,dr)
    call get_nx(veff_v6, nr)
    allocate(r(nr))
    call get_x(veff_v6, r)
    
    call get_nop(veff_v6,nop)
    allocate(v6(nr,nop),vTS(nr,nop))
    v6(:,:)  = 0.d0
    vTS(:,:) = 0.d0
    call get_mat(veff_v6, v6)


    M(1,1) = +1.0d0
    M(1,2) = -3.0d0
    M(1,3) = -3.0d0
    M(1,4) = +9.0d0

    M(2,1) = +1.0d0
    M(2,2) = +1.0d0
    M(2,3) = -3.0d0
    M(2,4) = -3.0d0

    M(3,1) = +1.0d0
    M(3,2) = -3.0d0
    M(3,3) = +1.0d0
    M(3,4) = -3.0d0

    M(4,1) = +1.0d0
    M(4,2) = +1.0d0
    M(4,3) = +1.0d0
    M(4,4) = +1.0d0


    Mten(1,1) = +1.0d0
    Mten(1,2) = -3.0d0
    Mten(2,1) = +1.0d0
    Mten(2,2) = +1.0d0

    
    do i = 1, 4
       do j = 1, 4
          do ir = 1, nr
             vTS(ir,i) = vTS(ir,i) + M(i,j) * v6(ir,j)
          end do
       end do
    end do
    
    do i = 5, 6
       do j = 5, 6
          do ir = 1, nr
             vTS(ir,i) =vTS(ir,i) +  Mten(nop-i,nop-j) * v6(ir,j)
          end do
       end do
    end do
    
    
    call define_VEFFdata(rho, nu, dr,nr,r,nop,vTS,veff_TS)
    
  end subroutine change_toTS
  



  subroutine change_toV6(veff_TS,veff_v6)
    implicit none
    type(VEFFdata), intent(in) :: veff_TS
    type(VEFFdata), intent(out) :: veff_v6
    real(kind=r8) :: rho, dr
    real(kind=r8), dimension(:,:), allocatable :: vTS, v6
    real(kind=r8), dimension(:), allocatable :: r
    integer :: nr, nop, nu
    integer :: ir, i, j
    real(kind=r8), dimension(1:4,1:4) :: M
    real(kind=r8), dimension(1:2,1:2) :: Mten

    call get_rho(veff_TS,rho)
    call get_nu(veff_TS,nu)
    call get_dx(veff_TS,dr)

    call get_nx(veff_TS, nr)
    allocate(r(nr))
    call get_x(veff_TS, r)
    
    call get_nop(veff_TS,nop)
    allocate(v6(nr,nop),vTS(nr,nop))
    v6(:,:)  = 0.d0
    vTS(:,:) = 0.d0
    call get_mat(veff_TS, vTS)


    M(1,1) = +1.0d0
    M(1,2) = +3.0d0
    M(1,3) = +3.0d0
    M(1,4) = +9.0d0

    M(2,1) = -1.0d0
    M(2,2) = +1.0d0
    M(2,3) = -3.0d0
    M(2,4) = +3.0d0

    M(3,1) = -1.0d0
    M(3,2) = -3.0d0
    M(3,3) = +1.0d0
    M(3,4) = +3.0d0

    M(4,1) = +1.0d0
    M(4,2) = -1.0d0
    M(4,3) = -1.0d0
    M(4,4) = +1.0d0


    Mten(1,1) = +1.0d0
    Mten(1,2) = +3.0d0
    Mten(2,1) = -1.0d0
    Mten(2,2) = +1.0d0


    M(:,:) = 1.d0/16.d0 * M(:,:)
 
    do i = 1, 4
       do j = 1, 4
          do ir= 1, nr
             v6(ir,i)  = v6(ir,i)  + M(i,j) * vTS(ir,j)   !v6(ir,i) = M(i,j) * vTS(ir,j)
          end do
       end do
    end do
    
    do i = 5, 6
       do j = 5, 6
          do ir = 1, nr
             v6(ir,i) = v6(ir,i) + Mten(nop-i,nop-j) * vts(ir,j)
          end do
       end do
    end do
    
    call define_VEFFdata(rho,nu, dr,nr,r,nop,v6,veff_v6)
    
  end subroutine change_toV6
  
end module representation
