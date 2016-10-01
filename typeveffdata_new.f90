module typeveffdata
  use precision
  implicit none
  
  !integer, parameter, private :: r8 = selected_real_kind(15,9)
  
  type VEFFdata                                                    ! my type 
     real(kind=my_kind) :: rho                                     ! density
     integer ::  nu                                                ! degeneracy
     real(kind=my_kind) :: dx                                      ! step rgrid
     integer :: nx                                                 ! # points r grid
     real(kind=my_kind), dimension(:),  pointer :: px => null()    ! rgrid
     integer :: nop                                                ! # of operators
     real(kind=my_kind), dimension(:,:), pointer :: pmat => null() ! Effective interaction  veff(ir, n)
  end type VEFFdata
  
  
  interface assignment (=)
     module procedure copy_veffdata
  end interface assignment (=)
  
contains

  subroutine define_VEFFdata(rho, nu, dx,nx,x,nop,mat,obj)
    implicit none
    real(kind=my_kind) :: rho, dx
    integer :: nx, nop, nu
    real(kind=my_kind), dimension (1:nx)    ,  target  :: x
    real(kind=my_kind), dimension (1:nx,1:nop), target:: mat
    type(VEFFdata) :: obj
    logical, save :: first_entry = .true.
      
    if(first_entry) then
       nullify(obj%px)
       nullify(obj%pmat)
       first_entry = .false.
    end if
  
    
    obj%rho   = rho
    obj%nu    = nu
    obj%dx    = dx
    obj%nx    = nx
    if(.not.associated(obj%px)) then
       allocate( obj%px(1:nx))
    end if
    obj%px     = x
    obj%nop   = nop
    if(.not.associated(obj%pmat)) then
       allocate( obj%pmat(1:nx,1:nop))
    end if
    obj%pmat    = mat



    return
    
  end subroutine define_VEFFdata



  subroutine get_rho(obj,rho)
    implicit none
    type(VEFFdata), intent(in) ::  obj
    real(kind=my_kind) :: rho
    rho = obj%rho
    return
  end subroutine get_rho

  subroutine get_nu(obj,nu)
    implicit none
    type(VEFFdata), intent(in) ::  obj
    integer :: nu
    nu = obj%nu
    return
  end subroutine get_nu
  
  subroutine get_dx(obj,dx)
    implicit none
    type(VEFFdata), intent(in) ::  obj
    real(kind=my_kind) :: dx
    dx = obj%dx
    return
  end subroutine get_dx

  subroutine get_nx(obj,nx)
    implicit none
    type(VEFFdata) ::  obj
    integer:: nx
    nx = obj%nx
    return
  end subroutine get_nx
  
  subroutine get_nop(obj, nop) 
    implicit none
    type(VEFFdata), intent(in) ::  obj
    integer :: nop
    nop = obj%nop
    return
  end subroutine get_nop

  subroutine get_x(obj, x) 
    implicit none
    type(VEFFdata), intent(in) ::  obj
    real(kind=my_kind), dimension(:) :: x
    
    x = obj%px
    
    return
  end subroutine  get_x

  
  subroutine get_mat(obj, mat) 
    implicit none
    type(VEFFdata), intent(in) ::  obj
    real(kind=my_kind), dimension(:,:):: mat

       mat(:,:) = obj%pmat

    return
  end subroutine  get_mat

  
  subroutine get_mat_n(obj, n, mat_n) 
    implicit none
    type(VEFFdata), intent(in) ::  obj
    integer, intent(in) :: n
    real(kind=my_kind), dimension(:):: mat_n
    integer :: i, imax
    
  
       imax = size(mat_n)
       do i = 1 , imax
          mat_n(i) = obj%pmat(i,n)
       end do
   
    return
  end subroutine  get_mat_n
  
  
  
  subroutine copy_veffdata(vout, vin)
    implicit none
    type(VEFFdata), intent(in)    :: vin
    type(VEFFdata), intent(out)   :: vout
    real(kind=my_kind) :: rho_in, dx_in
    integer :: nu_in,  nx_in, nop_in
    real(kind=my_kind), dimension(:),   allocatable :: x_in
    real(kind=my_kind), dimension(:,:), allocatable :: mat_in
    
    call get_rho(vin,rho_in)
    call get_nu( vin, nu_in)
    call get_dx( vin, dx_in)
    call get_nx( vin, nx_in)
    allocate(x_in(1:nx_in))
    call get_x( vin, x_in)
    call get_nop(vin, nop_in)
    allocate(mat_in(1:nx_in,1:nop_in))
    call get_mat(vin, mat_in)

    call define_VEFFdata(rho_in, nu_in, dx_in,nx_in,x_in,nop_in,mat_in,vout)

    return    
  end subroutine copy_veffdata
  
  
  
end module typeveffdata
