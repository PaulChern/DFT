!---------------------------------------------------------------------------------------------------
!
!	filename = coupling.f90
!	authors = Edison Salazar
!		  Pedro Guarderas
!	date = 11/01/2014
!
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Module for constructing Khon-Sham equations
module makeks
  implicit none

contains 

!---------------------------------------------------------------------------------------------------
! Building Khon-Shan system of equations
subroutine buildKS( nc, l, d, n, m )
  implicit none
  integer, intent( in ) :: n, m, d(n+1)
  double precision, intent( in ) :: l(n), nc(m)
  double precision, intent( out ) ::  E(m,m), T(m,m)
  integer :: i, j, k, Mi, Ni
  double precision :: nrmj, nrmk, li

  Mi = 1
  Ni = 0
  do i = 1,n
     Mi = Mi + d( i )
     Ni = d( i + 1 ) + Mi - 1
     !$omp parallel do private(j,k,nrmj,nrmk) shared(Mi,Ni,e,nc,l,T,E) collapse(2)
     do j = Mi,Ni
        nrmj = SLTnorm( n(j), e(j) )
        do k = Mi,Ni

           nrmk = SLTnorm( n(k), e(k) )

           T( j, k ) = 0.5D0 * ( ( l(i) * ( l(i) + 1.0D0 ) - nc(k) * ( nc(k) - 1.0D0 ) ) * &
                SLTnorm( nc(j) + nc(k) - 1.0D0, e(j) + e(k) ) +  &
                2.0D0 * e(k) * nc(k) * SLTnorm( nc(j) + nc(k), e(j) + e(k) ) - &
                e(k) * e(k) * SLTnorm( nc(j) + nc(k) + 1.0D0, e(j) + e(k) ) ) / &
                ( nrmj * nrmk )

           E( j, k ) = gamma( nc(j) + nc(k) + 1.0D0 ) / ( e(j) + e(k) )**( nc(j) + nc(k) + 1.0D0 )

        end do
     end do
  end do

end subroutine buildKS

end module makeks
