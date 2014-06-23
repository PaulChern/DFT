!---------------------------------------------------------------------------------------------------
!
!	filename = quad.f90
!	authors = Edison Salazar
!		  Pedro Guarderas
!	date = 18/10/2013
!
!---------------------------------------------------------------------------------------------------

module quad
  use consts
  implicit none
contains

function integrateRTH( rho, r, th, Nr, Nth ) result( intgrt )
  integer, intent( in ) :: Nr, Nth
  double precision, intent( in ) :: r(Nr), th(Nth), rho(Nr,Nth)
  double precision :: intgrt
  integer :: i, j
  
  intgrt = 0.0D0
  !#omp parallel do shared(Nr,Nth,intgrt,rho,r,th) collapse(2)
  do i = 1,(Nr-1)
     do j = 1,(Nth-1)
        intgrt = intgrt + pi2_ * rho( i, j ) * ( ( r(i)**2.0D0 )* sin( th(j) ) ) * &
             ( r(i+1) - r(i) ) * ( th(j+1) - th(j) )
     end do
  end do
  
end function integrateRTH

end module quad
