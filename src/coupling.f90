!---------------------------------------------------------------------------------------------------
!
! filename = coupling.f90
! authors = Edison Salazar
!           Pedro Guarderas
! date = 11/01/2014
!
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Module for constructing Khon-Sham equations
module coupling
  use base
  implicit none

contains 

!---------------------------------------------------------------------------------------------------
! Building Khon-Shan system of equations
subroutine buildKS( H, T, nc, e, ne, l, d, n, m )
  implicit none
  integer, intent( in ) :: n, m, d(n+1)
  double precision, intent( in ) :: l(n), nc(m), e(m), ne(n)
  double precision, intent( out ) ::  H(m,m), T(m,m)
  integer :: i, j, k, Mi, Ni
  double precision :: nrmj, nrmk

  !$omp parallel do private(i,j) shared(H,T) collapse(2)
  do i = 1,m
     do j = 1,m
        T(i,j) = 0.0D0
        H(i,j) = 0.0D0
     end do
  end do

  Mi = 1
  Ni = 0

  do i = 1,n

      Mi = Mi + d( i )
      Ni = d( i + 1 ) + Mi - 1

      do j = Mi,Ni

        nrmj = 1.0D0 / sqrt( SLTnorm( 2.0D0 * nc(j), 2.0D0 * e(j) ) )

        do k = Mi,Ni

           nrmk = 1.0D0 / sqrt( SLTnorm( 2.0D0 * nc(k), 2.0D0 * e(k) ) )

!            T( j, k ) = 0.5D0 * ne(i) * ( ( l(i) * ( l(i) + 1.0D0 ) - nc(k) * ( nc(k) - 1.0D0 ) ) * &
!                 SLTnorm( nc(j) + nc(k) - 1.0D0, e(j) + e(k) ) +  &
!                 2.0D0 * e(k) * nc(k) * SLTnorm( nc(j) + nc(k), e(j) + e(k) ) - &
!                 e(k) * e(k) * SLTnorm( nc(j) + nc(k) + 1.0D0, e(j) + e(k) ) ) / &
!                 ( nrmj * nrmk )
                
          T( j, k ) = 0.5D0 * ne(i) * nrmj * nrmk * ( &
            ( nc(i) - 1.0D0 ) * ( nc(j) - 1.0D0 ) * &
            SLTnorm( nc(i) + nc(j) - 2.0D0, e(i) + e(j) ) - &
            
            ( ( nc(i) - 1.0D0 ) * e(j) + ( nc(j) - 1.0D0 ) * e(i) ) * &
            SLTnorm( nc(i) + nc(j) - 1.0D0, e(i) + e(j) ) + &
            
            e(i) * e(j) * &
            SLTnorm( nc(i) + nc(j), e(i) + e(j) ) + &
            
            l(i) * ( l(i) + 1.0D0 ) * &
            SLTnorm( nc(i) + nc(j) - 2.0D0, e(i) + e(j) ) )

           H( j, k ) = dgamma( nc(j) + nc(k) + 1.0D0 ) / ( e(j) + e(k) )**( nc(j) + nc(k) + 1.0D0 )

        end do
     end do
  end do

end subroutine buildKS

!---------------------------------------------------------------------------------------------------
! Compute the potential
subroutine effectivePotential( Vef, V, p, r, th, Nr, Nth )
  implicit none
  integer, intent( in ) :: Nr, Nth
  double precision, intent( in ) :: r(Nr), th(Nth), p(Nr,Nth), V( Nr, Nth, Nr, Nth )
  double precision, intent( out ) :: Vef(Nr,Nth)
  integer :: i, j, k, l
  double precision :: dx, dy
  
  !$omp parallel do private(i,j,k,l,dx,dy) shared(Vef,V,p,th) collapse(4)
  do i = 1,Nr
    do j = 1,Nth
      do k = 2,Nr
        do l = 2,Nth
          if ( i /= k .and. j /= l ) then
            dx = r( k ) - r( k - 1 );
            dy = th( l ) - th( l - 1 );
            Vef( i, j ) = Vef( i, j ) + p( k, l ) * V( i, j, k, l ) * &
              sin( th(l) ) * dx * dy;
          end if
        end do
      end do
    end do
  end do

end subroutine effectivePotential

end module coupling
