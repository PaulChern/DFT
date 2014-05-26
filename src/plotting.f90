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

!---------------------------------------------------------------------------------------------------
subroutine contour( A, m, n, xmin, xmax, sx, ymin, ymax, sy, zmin, zmax, sz, inter )
  implicit none
  integer, intent( in ) :: m, n, sx, sy, sz, inter
  double precision, intent( in ) :: A(m,n)
  real, intent( in ) :: xmin, xmax, ymin, ymax, zmin, zmax
  real :: hx, hy, hz
  
  hx = sngl( xmax - xmin ) / real( sx - 1 )
  hy = sngl( ymax - ymin ) / real( sy - 1 )
  hz = sngl( zmax - zmin ) / real( sz - 1 )

  call METAFL('CONS')
  call DISINI()
  call PAGERA()
  call HWFONT()

  call TITLIN('Function',2)

  call NAME('X-axis','r')
  call NAME('Y-axis','theta')
  call NAME('Z-axis','An')

  call INTAX()
  call AUTRES(m,n)
!   call AXSPOS(300,1850)
  call AX3LEN(1500,1500,1500)
  
  call GRAF3( xmin, xmax, xmin, hx, ymin, ymax, ymin, hy, zmin, zmax, zmin, hz )
  call CRVMAT( sngl(A), m, n, inter, inter )

  call HEIGHT(50)
  call TITLE()
  call DISFIN()

end subroutine contour

!---------------------------------------------------------------------------------------------------
subroutine surface( A, x, y, m, n, &
                    xmin, xmax, sx, ymin, ymax, sy, zmin, zmax, sz, vx, vy, vz )
  implicit none
  integer, intent( in ) :: m, n, sx, sy, sz
  character ( len = 60 ) :: tit1, tit2
  double precision, intent( in ) :: A(m,n), x(m), y(n)
  real, intent( in ) :: xmin, xmax, ymin, ymax, zmin, zmax, vx, vy, vz
  real :: hx, hy, hz
  
  hx = sngl( xmax - xmin ) / real( sx - 1 )
  hy = sngl( ymax - ymin ) / real( sy - 1 )
  hz = sngl( zmax - zmin ) / real( sz - 1 )
 
  tit1 = 'Function'
  tit2 = 'Surface'

  call SETPAG('DA4P')
  call METAFL('CONS')
  call DISINI()
  call PAGERA()
  call COMPLX()
!   call AXSPOS(10,10)
!   call AXSLEN(2500,2000)

  call NAME('X-axis','X')
  call NAME('Y-axis','Y')
  call NAME('Z-axis','Z')
 
  call TITLIN( tit1, 2 )
  call TITLIN( tit2, 4 )

  call VIEW3D( vx, vy, vz,'ANGLE' )
  call GRAF3D( xmin, xmax, xmin, hx, ymin, ymax, ymin, hy, zmin, zmax, zmin, hz )
  call HEIGHT(50)
  call TITLE()

  call SHDMOD('SMOOTH','SURFACE')
  call SURSHD( sngl(x), m, sngl(y), n, sngl(A) )

  call DISFIN()

end subroutine surface


end module plotting

