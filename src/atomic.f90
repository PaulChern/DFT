!---------------------------------------------------------------------------------------------------
!
!	filename = atomic.f90
!	authors = Edison Salazar
!		  Pedro Guarderas
!	date = 01/08/2013
!
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Module for atomic calculus
! Based in notes of Eduardo Ludeña
module atomic
  use base
  use consts

  implicit none

contains

!---------------------------------------------------------------------------------------------------
subroutine radialDFT( p, r, Nr, c, nc, e, ne, d, n, m )
  implicit none
  integer, intent( in ) :: n, m, Nr, d(n+1)
  double precision, intent( in ) :: c(m), ne(n), nc(m), e(m), r(Nr)
  double precision, intent( out ) :: p(Nr)
  double precision :: sum, nrmj
  integer i, j, k, Mi, Ni

  do k = 1,Nr
    Mi = 1
    Ni = 0
    do i = 1,n
      Mi = Mi + d(i)
      Ni = d( i + 1 ) + Mi - 1
      sum = 0.0D0
      do j = Mi,Ni
        nrmj = 1.0D0 / sqrt( STOnorm( 2.0D0 * nc(j), 2.0D0 * e(j) ) )
        sum = sum + c(j) * nrmj * ( r(k)**( nc(j) - 1.0D0 ) ) * exp( -e(j) * r(k) )
      end do
      p(k) = p(k) + ne(i) * ( r(k) * sum )**2.0D0
    end do
  end do
  
end subroutine radialDFT

end module atomic

