!=======================================================================
!
module Boundary
  use precision
  implicit none
  real(kind=my_kind) :: lowerIM(10),upperIM(10)
  real(kind=my_kind) :: lowerRE(10),upperRE(10)
  logical  :: part
end module Boundary

!
!=======================================================================
!
