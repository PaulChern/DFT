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
  use grid
  use linear
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

  integer :: Nr, Nth, i
  double precision :: Rf
  double precision, allocatable :: r(:), th(:)
  double precision, allocatable :: P(:,:), Rb(:,:,:), RP(:,:,:)
  double precision, allocatable :: rho(:,:), An(:,:)

!---------------------------------------------------------------------------------------------------
! Reading parameters

  call getarg( 1, infile )

  open( unit = 10, file = infile )

  read( 10, * ) n
  read( 10, * ) m

  allocate( l(n), c(m), nc(m), e(m), ne(n), d(n+1) )
  
  read( 10, * ) l
  read( 10, * ) c
  read( 10, * ) nc
  read( 10, * ) e
  read( 10, * ) ne
  read( 10, * ) d

!---------------------------------------------------------------------------------------------------
! Grid generation

  Rf = 5.0D0
  Nr = 100
  Nth = 100

  ! Grid generation
  allocate( r( Nr ), th( Nth ) )

  write(*,*) "--------------------------------------------------------------------------------"

  write(*,*) "Grid generation for radial coordinates"
  write(*,*)
  call grid_adap_exp( r, Nr, Rf, 1.5D0 )

  write(*,*) "Grid generation for angular coordinates"
  write(*,*)
  call grid_unif( th, Nth, 0.0D0, pi_ )

!---------------------------------------------------------------------------------------------------
! DFT calculation
  write(*,*) "--------------------------------------------------------------------------------"
  write(*,*) "DFT"
  write(*,*)
  allocate( P(m,m), Rb(m,Nr,Nth), rho(Nr,Nth), RP(Nr,Nth,m), An(Nr,Nth) )

  write(*,*) "Density matrix P construction"
  write(*,*)
  call GDM( P, c, ne, d, n, m )

  write(*,*) "Base matrix R construction"
  write(*,*)
  call RBS( Rb, r, th, l, nc, e, d, Nr, Nth, n, m )

  write(*,*) "DFT Density function calculation"
  write(*,*)
  call DFT( rho, RP, P, Rb, d, n, m, Nr, Nth )

  write(*,*) "Enhancement factor An calculation"
  write(*,*)
  call AnFactor( An, RP, Rb, rho, r, th, l, nc, e, d, Nr, Nth, n, m )
!  call AnFactor( An, RP, P, Rb, rho, r, th, l, nc, e, d, Nr, Nth, n, m )

  write(*,*) "Integral of the density: ", integrate( rho, r, th, Nr, Nth )
  write(*,*) 

  write(*,*) "Integral of the enhancement factor T: ", integrate( An, r, th, Nr, Nth )
  write(*,*) 

!---------------------------------------------------------------------------------------------------
! Writing results
  write(*,*) "--------------------------------------------------------------------------------"
  write(*,*) "Writing grids radial an angular"
  write(*,*)
  open( unit = 10, file = "radial_grid.data" )
  call WriteVector( 10, r, Nr )
  close( unit = 10 )

  open( unit = 11, file = "angular_grid.data" )
  call WriteVector( 11, th, Nth )
  close( unit = 11 )

  write(*,*) "Writing density matrix P"
  write(*,*)
  open( unit = 12, file = "density_matrix.data" )
  call WriteMatrix( 12, P, m, m )
  close( unit = 12 )

  write(*,*) "Writing rho"
  write(*,*)
  open( unit = 13, file = "rho.data" )
  call WriteMatrix( 13, rho, Nr, Nth )
  close( unit = 13 )

  write(*,*) "Writing An"
  write(*,*)
  open( unit = 14, file = "enhancement.data" )
  call WriteMatrix( 14, An, Nr, Nth )
  close( unit = 14 )

!---------------------------------------------------------------------------------------------------
  deallocate( P, Rb, rho, RP, An, r, th )
  deallocate( l, c, nc, e, ne, d )
  
  stop

end program principal
