!---------------------------------------------------------------------------------------------------
!
!	filename = linear.f90
!	authors = Edison Salazar
!		  Pedro Guarderas
!	date = 20/09/2013
!
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Module for linear algebra
module linear
  implicit none

contains

!---------------------------------------------------------------------------------------------------
! Printing a matrix
subroutine WriteMatrix( unit, A, m, n )
  implicit none
  integer, intent( in ) :: unit, m, n
  double precision, intent( in ) :: A( m, n )
  integer :: i, j

  do i = 1,m
     do j = 1,n
        write( unit, "(F50.40,1X)", advance="no" ) A( i, j )
     end do
     write( unit, * ) 
  end do

end subroutine WriteMatrix

!---------------------------------------------------------------------------------------------------
! Printing a vector
subroutine WriteVector( unit, A, n )
  implicit none
  integer, intent( in ) :: unit, n
  double precision, intent( in ) :: A( n )
  integer :: i

  do i = 1,n
     write( unit, "(F50.40,1X)", advance="no" ) A( i )
     write( unit, * ) 
  end do

end subroutine WriteVector

!---------------------------------------------------------------------------------------------------
! Matrix multiplication
subroutine MatMul( C, A, B, m, n, p, q, TA, TB )
  implicit none
  double precision, allocatable, intent( out ) :: C(:,:)
  double precision, intent( in ) :: A(m,n), B(p,q)
  integer, intent( in ) :: m, n, p, q
  logical, intent( in ) :: TA, TB

  double precision :: sum
  integer :: i, j, k

  if ( TA .and. TB ) then
     if ( m == q ) then
        allocate( C(n,p) )
        !$omp parallel do private(sum,k) shared(n,p,m,A,B,C) collapse(2)
        do i = 1,n
           do j = 1,p
              sum = 0.0D0
              do k = 1,m
                 sum = sum + A(k,i) * B(j,k)
              end do
              C(i,j) = sum
           end do
        end do
     end if
  else if ( .not. TA .and. TB ) then
     if ( n == q ) then
        allocate( C(m,p) )
        !$omp parallel do private(sum,k) shared(n,p,m,A,B,C) collapse(2)
        do i = 1,m
           do j = 1,p
              sum = 0.0D0
              do k = 1,n
                 sum = sum + A(i,k) * B(j,k)
              end do
              C(i,j) = sum
           end do
        end do
     end if
  else if ( TA .and. .not. TB ) then
     if ( m == p ) then
        allocate( C(n,q) )
        !$omp parallel do private(sum,k) shared(n,p,m,A,B,C) collapse(2)
        do i = 1,n
           do j = 1,q
              sum = 0.0D0
              do k = 1,m
                 sum = sum + A(k,i) * B(k,j)
              end do
              C(i,j) = sum
           end do
        end do
     end if
  else if ( .not. TA .and. .not. TB ) then
     if ( n == p ) then
        allocate( C(m,q) )
        !$omp parallel do private(sum,k) shared(n,p,m,A,B,C) collapse(2)
        do i = 1,m
           do j = 1,q
              sum = 0.0D0
              do k = 1,n
                 sum = sum + A(i,k) * B(k,j)
              end do
              C(i,j) = sum
           end do
        end do
     end if
  end if

end subroutine MatMul

!---------------------------------------------------------------------------------------------------
! Vector norm
function norm( v, n ) result( nrm )
  integer, intent( in ) :: n
  double precision, intent( in ) :: v(n)
  double precision :: nrm
  integer :: i

  nrm = 0.0D0
  !#omp parallel do shared(n,nrm,v)
  do i = 1,n
     nrm = nrm + v(i) * v(i)
  end do
  
  nrm = sqrt( nrm )
  
end function norm

!---------------------------------------------------------------------------------------------------
! Scalar product
function scalar( v, w, n ) result( scl )
  integer, intent( in ) :: n
  double precision, intent( in ) :: v(n), w(n)
  double precision :: scl
  integer :: i

  scl = 0.0D0
  !#omp parallel do shared(n,scl,v,w)
  do i = 1,n
     scl = scl + v(i) * w(i)
  end do

end function scalar

!---------------------------------------------------------------------------------------------------
! Cholesky decomposition
subroutine CHBfac( L, A, n )
  integer, intent( in ) :: n
  double precision, intent( in ) :: A(n,n)
  double precision, intent( out ) :: L(n,n)
  integer :: i, j, k

  !$omp parallel do shared(L) collapse(2)
  do i = 1,n
    do j = 1,n
      if ( i >= j ) then
        L(i,j) = A(i,j)
      else
        L(i,j) = 0.0D0
      end if
    end do
  end do

  do j = 1,n
    L(j,j) = sqrt(L(j,j))
    do i = (j+1),n
      L(i,j) = L(i,j) / L(j,j)
    end do
    do k = (j+1),n
      do i = k,n
        L(i,k) = L(i,k) - L(i,j) * L(k,j)
      end do
    end do
  end do

end subroutine CHBfac

!---------------------------------------------------------------------------------------------------
! QR factorization no destructive version
subroutine QRfac( Q, R, A, n )
  integer, intent( in ) :: n
  double precision, intent( in ) :: A(n,n)
  double precision, intent( out ) :: Q(n,n), R(n,n)
  integer :: i, j, k

  !$omp parallel do shared(Q,R,A) collapse(2)
  do i = 1,n
    do j = 1,n
      Q(i,j) = A(i,j)
      R(i,j) = 0.D0
    end do
  end do

  do j = 1,n
     R(j,j) = norm( Q(1:n,1), n )
     if ( R(j,j) == 0 ) then
        write( *, * ) "Error matrix has linearly dependent columns"
     else
        !$omp parallel do shared(Q,R)
        do i = 1,n
           Q(i,j) = Q(i,j) / R(j,j)
        end do
     end if

     do k = j+1,n
        R(j,k) = scalar( A(1:n,j), A(1:n,k), n )
        do i = 1,n
           Q(i,k) = Q(i,k) - Q(i,j) * R(j,k)
        end do
     end do
  end do

end subroutine

end module linear
