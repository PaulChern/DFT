!---------------------------------------------------------------------------------------------------
!
!	filename = consts.f90
!	authors = Edison Salazar
!		  Pedro Guarderas
!	date = 18/10/2013
!
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Module of constants
module consts
  implicit none
  ! Constants definition
  double precision :: pi_ = 4.0D0 * atan( 1.0D0 )
  double precision :: pi2_ = 8.0D0 * atan( 1.0D0 )
  double precision :: pi4_ = 16.0D0 * atan( 1.0D0 )

end module consts
