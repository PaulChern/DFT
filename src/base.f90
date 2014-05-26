!---------------------------------------------------------------------------------------------------
!
!	filename = base.f90
!	authors = Edison Salazar
!		  Pedro Guarderas
!	date = 14/09/2013
!
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Module for base functions
module base
  use consts
  implicit none

contains

!---------------------------------------------------------------------------------------------------
! Legendre polynomials
subroutine Legendre( x, P, N )
  implicit none
  integer i
  integer, intent( in ) :: N
  double precision, intent( in ) :: x
  double precision, intent( out ), dimension(0:N) :: P

  P(0) = 1.0
  P(1) = x
  do i = 1,N-1
     P(i+1) = ( ( 2.0D0 * dble(i) + 1.0D0 ) * x * P(i) - dble( i ) * P(i-1) ) / dble(i+1)
  end do
end subroutine Legendre

!---------------------------------------------------------------------------------------------------
! Spherical harmonics Y_l type
subroutine SphHrm( Y, N, th, Nth )
  implicit none
  integer, intent( in ) :: N, Nth
  double precision, intent( in ) :: th(Nth)
  double precision, intent( out ), dimension(0:N,Nth) :: Y
  integer :: l, i
  double precision :: x

  !!$omp parallel do private(i)
  do i = 1,Nth
     x = cos(th(i))
     Y(0,i) = sqrt( 1.0D0 / pi4_ )
     Y(1,i) = sqrt( 3.0D0 / pi4_ ) * x
     do l = 1,N-1
        Y(l+1,i) = sqrt( ( 2.0D0 * l + 3.0D0 ) / pi4_ ) * ( &
             ( 2.0D0 * dble(l) + 1.0D0 ) * x * Y(l,i) - dble( l ) * Y(l-1,i) ) / dble(l+1)
     end do
  end do
end subroutine SphHrm

!---------------------------------------------------------------------------------------------------
! Radial Hartree-Fock base function
function HFRb( r, n, e ) result( rb )
  implicit none
  double precision, intent( in ) :: r, n, e
  double precision :: rb
  
  rb = ( 1.0D0 / sqrt( dgamma( 2.0D0 * n + 1.0D0 ) * ( 2.0D0 * e )**( -2.0D0 * n - 1.0D0 ) ) ) * &
       ( r**( n - 1.0D0 ) ) * exp( -e * r )

end function HFRb

!---------------------------------------------------------------------------------------------------
function DrHFRb( r, n, e ) result( grb )
  implicit none
  double precision, intent( in ) :: r, n, e
  double precision :: grb
  
  grb = ( 1.0D0 / sqrt( dgamma( 2.0D0 * n + 1.0D0 ) * ( 2.0D0 * e )**( -2.0D0 * n - 1.0D0 ) ) ) * &
       ( ( n - 1.0D0 ) / r - e ) * ( r**( n - 1.0D0 ) ) * exp( -e * r )

end function DrHFRb

function SLTnorm( n, e ) result( sltn )
  implicit none
  double precision, intent( in ) :: n, e
  double precision :: sltn

!   sltn = ( ( 2.0D0 * e )**n ) * sqrt( 2.0D0 * e / dgamma( 2.0D0 * n + 1 ) )
  sltn = gamma( n + 1.0D0 ) * e**( -n - 1.0D0 )

end function SLTnorm

end module base
