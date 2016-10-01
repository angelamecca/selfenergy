  module precision
    implicit none
    integer, parameter :: sp_kind = kind(1.0)
    integer, parameter :: dp_kind = selected_real_kind(15,9)
# ifdef SINGLEPRECISION
    integer, parameter :: my_kind = sp_kind
#else
    integer, parameter :: my_kind = dp_kind
#endif
  end module precision
