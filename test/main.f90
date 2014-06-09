!---------------------------------------------------------------------------------------------------
!
!	filename = main.f90
!	authors = Edison Salazar
!		  Pedro Guarderas
!	date = 01/08/2013
!
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!
! Main program
!
!---------------------------------------------------------------------------------------------------
program principal

  use atomic
  use consts
  use coupling
  use grid
  use linear
  use plotting
  use quad

!---------------------------------------------------------------------------------------------------
!
! Atomic parameters
!
!  m   : Total number of bases
!  n   : Number of orbitals by type
!  l   : Second Quantic numbers for every type of orbital 
!  c   : Coefficients for the bases
!  nc  : Principal quantic number for all the base orbitals
!  e   : Coefficients for the base orbitals
!  ne  : Quantity of electrons for every orbital
!  d   : Quantity of coefficients for kind of orbital, starting in the second component
!
!---------------------------------------------------------------------------------------------------

  character( len = 100 ) :: infile
  integer :: n, m
  integer, allocatable :: d(:)
  double precision, allocatable :: l(:), c(:), nc(:), e(:), ne(:)
  double precision, allocatable :: P(:,:), T(:,:), U(:,:), y(:)

!---------------------------------------------------------------------------------------------------
! Reading parameters

  call getarg( 1, infile )

  open( unit = 10, file = infile )

  read( 10, * ) n
  read( 10, * ) m

  allocate( l(n), c(m), nc(m), e(m), ne(n), d(n+1), y(m) )
  
  read( 10, * ) l
  read( 10, * ) c
  read( 10, * ) nc
  read( 10, * ) e
  read( 10, * ) ne
  read( 10, * ) d
 
!---------------------------------------------------------------------------------------------------
  write(*,*) "DFT computations"
  write(*,*)

  allocate( T(m,m), P(m,m), U(m,m) )
	
  call buildKS( P, T, nc, e, ne, l, d, n, m )
  call MatMulVec( y, T, c, m, m )
  write(*,*) "Exact integral of the kinetic energy T: ", scalar( y, c, m )
  write(*,*)
  call MatMulVec( y, P, c, m, m )
  write(*,*) "Exact integral of density: ", scalar( y, c, m )
  write(*,*)
  call CHBfac( U, P, m )
  
  open( unit = 15, file = "kinetic.data" )
  call WriteMatrix( 15, T, m, m )
  close( unit = 15 )
  
  open( unit = 16, file = "density.data" )
  call WriteMatrix( 16, P, m, m )
  close( unit = 16 )
  
  open( unit = 17, file = "deschol.data" )
  call WriteMatrix( 17, U, m, m )
  close( unit = 17 )
  
  call contour( T, m, m, 0.0, real(m), 10, 0.0, real(m), 10, -446.843, 220.607, 10, 1 )
  call contour( P, m, m, 0.0, real(m), 10, 0.0, real(m), 10, 0.0, 3.0, 10, 1 )

!   call contour( U, m, m, 0.0, real(m), 10, 0.0, real(m), 10, -0.06284573, 12.16492, 10, 1 )
 
!---------------------------------------------------------------------------------------------------
  deallocate( l, c, nc, e, ne, d )
  deallocate( T, P, U, y )
 
  stop

end program principal
