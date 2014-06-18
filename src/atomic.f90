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

!---------------------------------------------------------------------------------------------------
! Create generalized density matrix
!
! P = density matrix constructed by blocks [ ne(i) C(i) C(i)^t ]
! c = coefficients vector
! ne = number of electrons by orbital
! d = base dimensions by orbital
! n = number of orbitals
! m = total matrix size
!
! subroutine GDM( P, c, ne, d, n, m )
!   implicit none
!   integer, intent( in ) :: n, m, d(n+1)
!   double precision, intent( in ) :: c(m), ne(n)
!   double precision, intent( out ) :: P(m,m)
!   integer :: i, Mi, Ni, j, k
!   
!   !$omp parallel do collapse(2)
!   do j = 1,m
!      do k = 1,m
!         P( j, k ) = 0.0D0
!      end do
!   end do
! 
!   ! Construction by blocks
!   Mi = 1
!   Ni = 0
!   do i = 1,n
!      Mi = Mi + d( i )
!      Ni = d( i + 1 ) + Mi - 1
!      !$omp parallel do private(j,k) shared(Mi,Ni,ne,c,P) collapse(2)
!      do j = Mi,Ni
!         do k = Mi,Ni
!            P( j, k ) = ne(i) * c(j) * c(k)
!         end do
!      end do
!   end do
! 
! end subroutine GDM
! 
! !---------------------------------------------------------------------------------------------------
! ! Create generalized base matrix in a given grid
! !
! ! RB = base vectors
! ! r = radial grid
! ! theta = angular grid
! ! l = quantic number
! ! nc = electronic configuration
! ! e = energy values
! ! d = base dimensions by orbital
! ! Nr = radial grid size
! ! Nth = angular grid size
! ! n = number of orbitals
! ! m = total matrix size
! !
! subroutine RBS( Rb, r, th, l, nc, e, d, Nr, Nth, n, m )
!   implicit none
!   integer, intent( in ) :: Nr, Nth, n, m, d(n+1)
!   double precision, intent( in ) :: r(Nr), th(Nth), nc(m), e(m), l(n)
!   double precision, intent( out ) :: Rb( m, Nr, Nth )
!   double precision :: Y( 0:int(l(n)), Nth )
!   integer :: i, j, i_r, i_th, Mi, Ni, li
! 
!   !$omp parallel do collapse(3)
!   do i_r = 1,Nr
!      do i_th = 1,Nth
!         do i = 1,m
!            Rb( i, i_r, i_th ) = 0.0D0
!         end do
!      end do
!   end do
! 
!   ! spherical harmonics matrix calculation
!   call SphHrm( Y, int(l(n)), th, Nth )
! 
!   !$omp parallel do private(i,j,Mi,Ni,li) shared(n,d,l,Rb,r,nc,e,Y) collapse(2)
!   do i_r = 1,Nr
!      do i_th = 1,Nth
!         Mi = 1
!         Ni = 0
!         do i = 1,n
!            Mi = Mi + d( i )
!            Ni = d( i + 1 ) + Mi - 1
!            li = int( l(i) )
!            do j = Mi,Ni
!               Rb( j, i_r, i_th ) = STOradial( r( i_r ), nc(j), e(j) ) * Y( li, i_th )
!            end do
!         end do
!      end do
!   end do
!   
! end subroutine RBS
! 
! !---------------------------------------------------------------------------------------------------
! ! Density function calculated with the generalized density matrix
! ! 
! ! rho = density 
! ! RP = left product base and density matrix R^t * P
! ! R = base
! ! m = total matrix size
! ! Nr = radial grid size
! ! Nth = angular grid size
! !
! subroutine DFT( rho, RbP, P, Rb, d, n, m, r, th, Nr, Nth )
!   implicit none
!   integer, intent( in ) :: n, m, Nr, Nth, d(n+1)
!   double precision, intent( out ) :: rho(Nr,Nth), RbP(Nr,Nth,m)
!   double precision, intent( in ) :: Rb(m,Nr,Nth), P(m,m), r(Nr), th(Nth)
!   double precision :: sum
!   integer :: i, j, k, Mi, Ni, i_r, i_th
! 
!   !$omp parallel do collapse(3)
!   do i_r = 1,Nr
!      do i_th = 1,Nth
!         do i = 1,m
!            RbP( i_r, i_th, i ) = 0.0D0
!         end do
!      end do
!   end do
! 
!   ! R^t * P product construction by blocks
!   !$omp parallel do private(i,j,k,Mi,Ni,sum) shared(n,d,Rb,P,RbP) collapse(2)
!   do i_r = 1,Nr
!      do i_th = 1,Nth
!         Mi = 1
!         Ni = 0
!         do i = 1,n
!            Mi = Mi + d( i )
!            Ni = d( i + 1 ) + Mi - 1
!            do j = Mi,Ni
!               sum = 0.0D0
!               do k = Mi,Ni
!                  sum = sum + Rb( k, i_r, i_th ) * P( k, j )
!               end do
!               RbP( i_r, i_th, j ) = sum
!            end do
!         end do
!      end do
!   end do
!   
!   !$omp parallel do private(j,sum) shared(m,r,th,RbP,Rb,rho) collapse(2)
!   do i_r = 1,Nr
!      do i_th = 1,Nth
!         sum = 0.0D0
!         do j = 1,m
!            sum = sum + RbP( i_r, i_th, j ) * Rb( j, i_r, i_th )
!         end do
!         rho(i_r,i_th) = sum * sin( th( i_th ) ) * r( i_r )**2
!      end do
!   end do
!   
! end subroutine DFT
! 
! !---------------------------------------------------------------------------------------------------
! ! Enhancement factor An
! !subroutine AnFactor( An, RP, P, Rb, rho, r, th, l, nc, e, d, Nr, Nth, n, m )
! subroutine AnFactor( An, RP, Rb, rho, r, th, l, nc, e, d, Nr, Nth, n, m )
!   implicit none
!   integer, intent( in ) :: Nr, Nth, m, n, d(n+1)
!   double precision, intent( in ) :: RP(Nr,Nth,m), Rb(m,Nr,Nth), rho(Nr,Nth)
! !  double precision, intent( in ) :: P(m,m)
!   double precision, intent( in ) :: r(Nr), th(Nth), l(n), nc(m), e(m)
!   double precision, intent( out ) :: An(Nr,Nth)
!   double precision :: DrR(m,Nr,Nth), DthR(m,Nr,Nth), Y( 0:int(l(n)), Nth ), LR(m,Nr,Nth)
! !   double precision :: DrRP(Nr,Nth,m), DthRP(Nr,Nth,m)
!   double precision :: Gr, Gth, Wzck, T
!   integer :: i_r, i_th, i, j, Mi, Ni, li
! !   double precision ::  sum_r, sum_th
! !   integer :: k
! 
!   !$omp parallel do shared(LR,DrR,DthR) collapse(3)
!   do i_r = 1,Nr
!      do i_th = 1,Nth
!         do i = 1,m
!            LR( i, i_r, i_th ) = 0.0D0
!            DrR( i, i_r, i_th ) = 0.0D0
!            DthR( i, i_r, i_th ) = 0.0D0
! !            DrRP( i_r, i_th, i ) = 0.0D0
! !            DthRP( i_r, i_th, i ) = 0.0D0
!         end do
!      end do
!   end do
! 
!   ! spherical harmonics matrix calculation
!   call SphHrm( Y, int(l(n)), th, Nth )
! 
!   !$omp parallel do private(i,j,Mi,Ni,li) shared(n,d,l,Rb,LR,DrR,DthR,th,nc,e,Y) collapse(2)
!   do i_r = 1,Nr
!      do i_th = 1,Nth
!         Mi = 1
!         Ni = 0
!         do i = 1,n
!            Mi = Mi + d( i )
!            Ni = d( i + 1 ) + Mi - 1
!            li = int( l(i) )
!            do j = Mi,Ni
!               
!               LR( j, i_r, i_th ) = ( -( ( nc(j) - 1.0D0 ) / r(i_r) - e(j) )**2.0D0 + &
!                    l(i) * ( l(i) + 1.0D0 ) / ( r(i_r)**2.0D0 ) ) * Rb( j, i_r, i_th )
! 
! !               LR( j, i_r, i_th ) = ( ( -e(j) * r(i_r) + ( nc(j) - 1.0D0 ) )**2.0D0 + &
! !                    l(i) * ( l(i) + 1.0D0 ) ) * Rb( j, i_r, i_th ) / ( r(i_r)**2.0D0 )
! 
! !               LR( j, i_r, i_th ) = ( -nc(j) + ( nc(j) - e(j) * r(i_r) )**2.0D0 - &
! !                    l(i) * ( l(i) + 1.0D0 ) ) * Rb( j, i_r, i_th ) / ( r(i_r)**2.0D0 )
! !              
!               DrR( j, i_r, i_th ) = ( ( nc(j) - 1.0D0 ) / r(i_r) - e(j) ) * Rb( j, i_r, i_th )
! 
!               if ( li > 0 .and. sin( th(i_th) ) > 0.0D0 ) then
!                  DthR( j, i_r, i_th ) = l(i) * ( cos( th(i_th) ) * Rb( j, i_r, i_th ) - &
!                       STOradial( r( i_r ), nc(j), e(j) ) * Y( li - 1, i_th ) ) / sin( th(i_th) )
!               end if
! 
!            end do
!         end do
!      end do
!   end do
! 
! !   !$omp parallel do private(i,j,k,Mi,Ni,sum_r,sum_th) shared(n,d,P,DrR,DthR,DrRP,DthRP) collapse(2)
! !   do i_r = 1,Nr
! !      do i_th = 1,Nth
! !         Mi = 1
! !         Ni = 0
! !         do i = 1,n
! !            Mi = Mi + d( i )
! !            Ni = d( i + 1 ) + Mi - 1
! !            do j = Mi,Ni
! !               sum_r = 0.0D0
! !               sum_th = 0.0D0
! !               do k = Mi,Ni
! !                  sum_r = sum_r + DrR( k, i_r, i_th ) * P( k, j )
! !                  sum_th = sum_th + DthR( k, i_r, i_th ) * P( k, j )
! !               end do
! !               DrRP( i_r, i_th, j ) = sum_r
! !               DthRP( i_r, i_th, j ) = sum_th
! !            end do
! !         end do
! !      end do
! !   end do
!   
! !   !$omp parallel do &
! !   !$omp& private(i,j,Mi,Ni,T,Wzck,Gr,Gth,sum_r,sum_th) &
! !   !$omp& shared(n,d,RP,DrR,DthR,DrRP,DthRP,rho,An) &
! !   !$omp& collapse(2)
! 
!   !$omp parallel do &
!   !$omp& private(i,j,Mi,Ni,T,Wzck,Gr,Gth) &
!   !$omp& shared(n,d,RP,LR,DrR,DthR,rho,An) &
!   !$omp& collapse(2)
!   do i_r = 1,Nr
!      do i_th = 1,Nth
!         T = 0.0D0
!         Gr = 0.0D0
!         Gth = 0.0D0
!         Mi = 1
!         Ni = 0
! !         sum_r = 0.0D0
! !         sum_th = 0.0D0
!         do i = 1,n
!            Mi = Mi + d( i )
!            Ni = d( i + 1 ) + Mi - 1
!            do j = Mi,Ni
!               T = T + RP( i_r, i_th, j ) * LR( j, i_r, i_th )
!               Gr = Gr + RP( i_r, i_th, j ) * DrR( j, i_r, i_th )
!               Gth = Gth + RP( i_r, i_th, j ) * DthR( j, i_r, i_th )
! !               sum_r = sum_r + DrRP( i_r, i_th, j ) * DrR( j, i_r, i_th )
! !               sum_th = sum_th + DthRP( i_r, i_th, j ) * DthR( j, i_r, i_th )
!            end do
!         end do
!         Wzck = ( Gr**2.0D0 + ( r(i_r)**(-2.0D0) ) * Gth**2.0D0 ) / rho( i_r, i_th )
! !        An( i_r, i_th ) = ( T - pi2_ * Wzck ) / ( ( pi2_ * rho( i_r, i_th ) )**( 5.0D0 / 3.0D0 ) )
! !         T = sum_r**2.0D0 + ( r(i_r)**(-2.0D0) ) * sum_th**2.0D0
! !        An( i_r, i_th ) =  0.5D0 * T
! 
!      end do
!   end do
! end subroutine AnFactor

end module atomic

