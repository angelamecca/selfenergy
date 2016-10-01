module set_mpi
  use mpi
  implicit none
  integer, parameter, private :: r8 = selected_real_kind(15,9)

  integer, parameter :: master = 0
  integer :: ierr, islave
  integer :: myid, nprocs, nsteps
  integer , save :: nslaves
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, parameter :: tag_send = 2000, tag_rec=3000
  !
  integer :: ixy, nxy
  integer :: ixy_start, ixy_end,  nxy_send
  integer :: ixy_start_rec, ixy_end_rec,  nxy_rec
  real(kind=r8), dimension(:), allocatable :: im_vect, im_vect_err, im_vect_fail
  real(kind=r8), dimension(:), allocatable :: tmp_im_vect, tmp_im_vect_err, tmp_im_vect_fail

contains
  subroutine initialize()
    implicit none
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)

    nslaves = nprocs -1
    
    return
  end subroutine initialize

  subroutine finalize()
    implicit none
    call MPI_FINALIZE(ierr)
    return
  end subroutine finalize

  ! MASTER : SEND & REC
  subroutine master_send_int(n, id)
    implicit none
    integer :: n, id
    call MPI_SEND(n, 1, MPI_INT, id, tag_send, MPI_COMM_WORLD, ierr)
    return
  end subroutine master_send_int

  subroutine master_send_real_vect(x, n1, nx, id)
    implicit none
    real(kind=r8), dimension(:) :: x
    integer :: n1, nx, id
    call MPI_SEND(x(n1), nx, MPI_DOUBLE_PRECISION, id, tag_send,MPI_COMM_WORLD, ierr)
    return
  end subroutine master_send_real_vect

  subroutine master_receive_int(n,id)
    implicit none
    integer :: n, id
    call MPI_RECV(n,1, MPI_INT, id, tag_rec, MPI_COMM_WORLD, status, ierr)
    return
  end subroutine master_receive_int

  subroutine master_receive_real_vect(x,n1,nx,id)
    implicit none
    real(kind=r8), dimension(:) :: x
    integer :: n1, nx, id
    call MPI_RECV(x(n1),nx, MPI_INT, id, tag_rec, MPI_COMM_WORLD, status, ierr)
    return
  end subroutine master_receive_real_vect


  ! SLAVE SEND & REC
  subroutine slave_receive_int(n)
    implicit none
    integer :: n
    call MPI_RECV(n,1, MPI_INT, master, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)    
    return
  end subroutine slave_receive_int

  subroutine slave_receive_real_vect(x,n1,nx)
    implicit none
    integer :: n1, nx
    real(kind=r8), dimension(:) :: x
    call MPI_RECV(x(n1),nx, MPI_INT, master, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)    
    return
  end subroutine slave_receive_real_vect

  subroutine slave_send_int(n)
    implicit none
    integer :: n
    call MPI_SEND(n,1, MPI_INT, master, tag_rec, MPI_COMM_WORLD, ierr)
    return
  end subroutine slave_send_int

  subroutine slave_send_real_vect(x,n1,nx)
    implicit none
    integer :: n1, nx
    real(kind=r8), dimension(:) :: x
    call MPI_SEND(x(n1),nx, MPI_DOUBLE_PRECISION, master, tag_rec,MPI_COMM_WORLD, ierr)
    return
  end subroutine slave_send_real_vect


! BCAST
  subroutine bcast_int(n)
    implicit none
    integer :: n
    call MPI_BCAST (n,1,MPI_INT,master,MPI_COMM_WORLD,ierr)
    return
  end subroutine bcast_int
  
  subroutine bcast_real(x)
    implicit none
    real(kind=r8) :: x
    call MPI_BCAST (x,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    return
  end subroutine bcast_real

  subroutine bcast_log(val)
    implicit none
    logical :: val
    call MPI_BCAST (val,1,MPI_LAND,master,MPI_COMM_WORLD,ierr)
   end subroutine bcast_log

   subroutine bcast_char(name)
     implicit none
     character :: name
     call MPI_BCAST (name,len(name),MPI_INT,master,MPI_COMM_WORLD,ierr)
   end subroutine bcast_char

   
   subroutine bcast_int_vect(ivec)
     implicit none
     integer, dimension(:) :: ivec
     call MPI_BCAST (ivec,size(ivec),MPI_INT,master,MPI_COMM_WORLD,ierr)
     !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     return
   end subroutine bcast_int_vect
  
  subroutine bcast_real_vect1d(rvec)
    implicit none
    integer :: i
    real(kind=r8), dimension(:) :: rvec
    call MPI_BCAST (rvec,size(rvec),MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
   ! call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    return
  end subroutine bcast_real_vect1d

  subroutine bcast_real_vect2d(rvec)
    implicit none
    real(kind=r8), dimension(:,:) :: rvec
    call MPI_BCAST (rvec,size(rvec),MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    ! call MPI_BARRIER(MPI_COMM_WORLD,ierr)       
    return
  end subroutine bcast_real_vect2d

  subroutine bcast_real_vect3d(rvec)
    implicit none
    real(kind=r8), dimension(:,:,:) :: rvec
    call MPI_BCAST (rvec,size(rvec),MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
   ! call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    return
  end subroutine bcast_real_vect3d






  !split
  
end module set_mpi
