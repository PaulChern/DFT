!---------------------------------------------------------------------------------------------------
!
! filename = plotting.f90
! authors = Edison Salazar
!     Pedro Guarderas
! date = 20/05/2014
!
!---------------------------------------------------------------------------------------------------

! Plotting routines with DISLIN
module plotting
!   use dislin
  implicit none
  
contains
subroutine contour( A, m, n, p, xmin, xmax, ymin, ymax, zmin, zmax )
  implicit none
  integer, intent( in ) :: m, n, p
  double precision, intent( in ) :: A(m,n), xmin, xmax, ymin, ymax, zmin, zmax
  double precision :: hx, hy, hz
  
  hx = abs( ( xmax - xmin ) / dble( m - 1 ) )
  hy = abs( ( ymax - ymin ) / dble( n - 1 ) )
  hz = abs( ( zmax - zmin ) / dble( p - 1 ) )

  call METAFL('CONS')
  call DISINI()
  call PAGERA()
  call HWFONT()

  call TITLIN('Function',1)

  call NAME('X-axis','r')
  call NAME('Y-axis','theta')
  call NAME('Z-axis','An')

  call INTAX()
  call AUTRES(m,n)
  call AXSPOS(500,900)
!   call AX3LEN(2200,1400,1400)

  call GRAF3( xmin, xmax, xmin, hx, ymin, ymax, ymin, hy, zmin, zmax, zmin, hz )
  call CRVMAT( A, m, n, 1, 1 )

  call HEIGHT(50)
  call TITLE()
  call MPAEPL(3)
  call DISFIN()

end subroutine contour

end module plotting

