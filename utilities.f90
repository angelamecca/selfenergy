module utilities
  use typeveffdata
  integer, parameter, private :: r8 = selected_real_kind(15,9)   


contains
  
  subroutine read_table_veff(filename, veffmat_in, irhomax, nop, nu)
    implicit none
    
    character(len=64) :: filename
    integer :: irhomax, irho, nop, nu
    integer :: io, i, ir
    real(kind=r8) ::  tmp1, tmp2
    integer  itmp3
    type(VEFFdata), dimension(:), allocatable  :: veffmat_in
    real(kind=r8), dimension(:,:), allocatable :: v_val
    real(kind=r8), dimension(:), allocatable   :: r_val

    open(1,file=trim(filename))     
    irhomax = 0
    do
       read(1,*,iostat=io) tmp1, tmp2, itmp3
       if(io.eq.0) then
          irhomax = irhomax + 1
          do i = 1, itmp3
             read(1,*)
          end do
       else if(io.gt.0) then
          write(*,*) ' Problem in reading veffmat.dat'
          stop
       else
          exit
       end if
    end do
    
    rewind(1)
    
    write(*,*)
    write(*,*) ' # of densities read : ',  irhomax
    write(*,*)
    
    
    
    allocate(veffmat_in(1:irhomax))
    
    do irho = 1, irhomax
       read(1,*) tmp1, tmp2, itmp3
       allocate(r_val(itmp3), v_val(itmp3,nop))
       do ir = 1, itmp3
          read(1,*)  r_val(ir), v_val(ir,1), v_val(ir,2), v_val(ir,3), v_val(ir,4), v_val(ir,5), v_val(ir,6)
       end do
       call define_VEFFdata(tmp1, nu, tmp2, itmp3, r_val, nop, v_val, veffmat_in(irho))        
       deallocate(r_val, v_val)
    end do
    
    close(1)
    return
  end subroutine read_table_veff
  
end module utilities
