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
  integer :: n, m, Nr
  integer, allocatable :: d(:)
  double precision, allocatable :: l(:), c(:), nc(:), e(:), ne(:)
  double precision, allocatable :: DFT(:,:), T(:,:), y(:), p(:), r(:)
  double precision, dimension(3) :: Ct, W, J
  double precision :: Tapp

!---------------------------------------------------------------------------------------------------
! Reading parameters
  write(*,*) "Reading file with parameters"
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

  Ct = (/ 3.26422D0, -0.02631D0, 0.000498D0 /)
  W = (/ 5.0D0/3.0D0, 4.0D0/3.0D0, 11.0D0/9.0D0 /)
  J = (/ 1.0D0, 2.0D0, 3.0D0 /)
  
!---------------------------------------------------------------------------------------------------
! Grid generations
  Nr = 1000
  write(*,*) "Building radial mesh with size: ", Nr
  allocate( r(Nr), p(Nr) )
  call grid_adap_exp( r, Nr, 5.0D0, 5.0D0 )
  !call grid_unif( r, Nr, 0.0D0, 5.0D0 )
 
!---------------------------------------------------------------------------------------------------
  write(*,*) "DFT computations"

  allocate( T(m,m), DFT(m,m) )
	
  write(*,*) "Building DFT matrices"
  call buildKS( DFT, T, nc, e, ne, l, d, n, m )
  
  write(*,*) "Computing exact kinetic energy"
  open( unit = 14, file = 'results.data' )
  call MatMulVec( y, T, c, m, m )
  write(14,*) "Exact integral of the kinetic energy T: ", scalar( y, c, m )
  
  write(*,*) "Computing exact density"
  call MatMulVec( y, DFT, c, m, m )
  write(14,*) "Exact integral of density: ", scalar( y, c, m )
  
  
  write(*,*) "Computing radial density"
  call radialDFT( p, r, Nr, c, nc, e, ne, d, n, m )
  
  write(*,*) "Computing Kinectic Energy approximation"
  call HomogeneousApprox( Tapp, Ct, W, J, p, r, Nr, 3 )
  write(14,*) "Kinectic Energy approximation: ",  Tapp
  close( unit = 14 )

!---------------------------------------------------------------------------------------------------
  write(*,*) "Writing results in *.data files"
  
  write(*,*) "Writing kinetic energy"
  open( unit = 15, file = "kinetic.data" )
  call WriteMatrix( 15, T, m, m )
  close( unit = 15 )
  
  write(*,*) "Writing density matrix"
  open( unit = 16, file = "density.data" )
  call WriteMatrix( 16, DFT, m, m )
  close( unit = 16 )
  
  write(*,*) "Writing radial grid"
  open( unit = 17, file = "radial_grid.data" )
  call WriteVector( 17, r, Nr )
  close( unit = 17 )
 
  write(*,*) "Writing radial density"
  open( unit = 18, file = "radial_density.data" )
  call WriteVector( 18, p, Nr )
  close( unit = 18 )
  
  
!---------------------------------------------------------------------------------------------------
  write(*,*) "Plotting"
  call contour( T, m, m, 0.0, real(m), 10, 0.0, real(m), 10, -446.843, 220.607, 10, 1 )
  call contour( DFT, m, m, 0.0, real(m), 10, 0.0, real(m), 10, 0.0, 3.0, 10, 1 )


!---------------------------------------------------------------------------------------------------
  deallocate( l, c, nc, e, ne, d )
  deallocate( T, DFT, y )
  deallocate( r, p )
 
  stop

end program principal
