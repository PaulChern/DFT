!---------------------------------------------------------------------------------------------------
!
!	filename = grid.f90
!	authors = Edison Salazar
!		  Pedro Guarderas
!	date = 19/10/2013
!
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Module for grid generations
module grid
  implicit none

contains

!---------------------------------------------------------------------------------------------------
! Adaptative grid
subroutine grid_adap_exp( r, N, Rf, alpha )
  integer, intent( in ) :: N
  double precision, intent( in ) :: Rf, alpha
  double precision, intent( out ) :: r( N )
  double precision :: a, b, h
  integer :: i

  h = 1.0D0 / ( dble( N ) - 1.0D0 )
  a = Rf / ( exp( alpha ) - 1.0D0 )
  b = Rf * exp( -alpha ) / ( exp(-alpha) - 1.0D0 )

  !$omp parallel do private(i)
  do i = 1,N
     r(i) = a * exp( alpha * h * dble( i - 1.0D0 ) ) + b
  end do

end subroutine grid_adap_exp

!---------------------------------------------------------------------------------------------------
! Uniform grid
subroutine grid_unif( r, N, a, b )
  integer, intent( in ) :: N
  double precision, intent( in ) :: a, b
  double precision, intent( out ) :: r(N)
  double precision :: h
  integer :: i

  h = ( b - a ) / ( dble( N ) - 1.0D0 )

  !$omp parallel do private(i)
  do i = 1, N
     r(i) = dble( i - 1.0D0 ) * h
  end do

end subroutine grid_unif

end module grid
