
PROGRAM ccdt
  USE parallel
  USE ang_mom_functions
  USE hdf5
  USE constants
  USE single_particle_orbits
  USE chiral_constants
  USE contracts, ONLY: debug_contracts
  USE mem_tracker
  ! USE deltafull_parameters
  
  IMPLICIT NONE
  CHARACTER(len=100) :: tmp2,tmp3
  REAL(dp) :: startwtime, endwtime
  INTEGER, DIMENSION (10) :: Nocclist = (/ 0,2,14,38,54,66,114,162,186,246 /)
  CHARACTER(len=256) :: input(3)=[character(len=256) :: '&args','','/']
  CHARACTER(len=256) :: message
  INTEGER :: i, slash, ios, error, cc_approx0
  
  ! ! Define arguments and their default values
  ! !
  ! Nmax, Nocc, Pocc, rho
  ! cc_approx: 0=CCD, 1=CCD(T), 2=CCDT1, -1=MBPT3
  ! tnf_approx: 0=2nf, 1=2nf+3nfNO, 2=2nf+3nfNO+3nfT3
  ! pre_gs: 1=load gs_file
  ! cc_file: file 'tag' (use cc_file='"tag"' or cc_file="'tag'")
  ! add_n, add_p: EOM parameter
  ! eom_approx: 0=CCD, 1=CCDT1
  ! test: 1=int, 2=hbar, 3=gs, 4=eom
  ! !
  
  ! Setup Parallelization and HDF5
  CALL h5open_f(error)
  CALL mpi_init(ierror)
  CALL mpi_comm_rank(mpi_comm_world,iam,ierror)
  CALL mpi_comm_size(mpi_comm_world,num_procs,ierror)
  master = 0
  startwtime = MPI_WTIME()
  ! For OMP parallelization
  threads = omp_get_max_threads()
  
  ! read arguments from command line
  DO i = 1, command_argument_count()
     CALL get_command_argument(i,input(2))
     read(input,nml=args,iostat=ios,iomsg=message)
     IF (ios .ne. 0) then
        IF ( iam == 0 ) write(*,*)'ERROR:'//trim(message)
        stop
     end IF
  end DO

  pre_gs0 = .false.
  IF ( pre_gs == 1 ) pre_gs0 = .true.
  debug_contracts = ( debug_contracts_flag /= 0 )
  below_ef_n = Nocclist(Nocc+1)
  below_ef_p = Nocclist(Pocc+1)
  below_ef = below_ef_n + below_ef_p
  IF ( cc_approx < 2 .and. tnf_approx == 2 ) tnf_approx = 1
  IF ( t3_cut < 0 ) t3_cut = 0
  IF ( eom3_cut < 0 ) eom3_cut = 0

  ! input/output file
  ! get interaction name
  slash = index(cc_file, '/', .true.) ! true for last slash in string
  if (slash > 0) then
     interaction = trim(cc_file(slash+1:))
  else
     interaction = trim(cc_file)
  end if
  ! write CC file name
  cc_file0 = cc_file
  cc_approx0 = cc_approx
  IF ( cc_approx == 1 ) cc_approx0 = 0
  IF ( cc_approx >= 0 ) then
     write(tmp2,'(3(I2.2,A1),I3.3,2(A1,I1.1))') Nmax,'_',Nocc,'_',Pocc,'_',int(1000*rho),'_',cc_approx0,'_',tnf_approx
     IF ( cc_approx == 2 .and. t3_cut < 100 ) then
        write(tmp3,'(I2.2)') t3_cut
        tmp2 = trim(tmp2) // '_' // trim(tmp3)
     end IF
     cc_file = trim(cc_file0) // '_' // trim(tmp2) // '.h5'
  end IF
  
  IF ( iam == 0 ) then
     write(6,*)
     write(6,*) '...Starting Coupled-Cluster Infinite Matter Calculation!...'
     write(6,*)
     write(6,'(A16,5x,I4)') 'Neutron Shells:', Nocc
     write(6,'(A15,6x,I4)') 'Proton Shells:', Pocc
     write(6,'(A12,9x,I4)') '# Neutrons:', below_ef_n
     write(6,'(A11,10x,I4)') '# Protons:', below_ef_p
     write(6,'(A14,7x,I4)') 'Total Shells:', Nmax
     write(6,'(A9,13x,F5.3,4x,A5)') 'Density:', rho, 'fm^-3'
     write(6,'(A13,11x,A20)') 'Interaction:', interaction
     IF ( cc_approx == 0 ) write(6,'(A11,11x,A3)') 'CC-Approx:', 'CCD'
     IF ( cc_approx == 1 ) write(6,'(A11,11x,A6)') 'CC-Approx:', 'CCD(T)'
     IF ( cc_approx == 2 ) write(6,'(A11,11x,A6)') 'CC-Approx:', 'CCDT-1'
     IF ( cc_approx == -1 ) write(6,'(A11,11x,A3)') 'MBPT3'
     IF ( tnf_approx == 0 ) write(6,'(A12,10x,A3)') '3NF-Approx:', '2NF'
     IF ( tnf_approx == 1 ) write(6,'(A12,10x,A11)') '3NF-Approx:', '2NF+3NF(no)'
     IF ( tnf_approx == 2 ) write(6,'(A12,10x,A19)') '3NF-Approx:', '2NF+3NF(no)+3NF(T3)'
     IF ( pre_gs == 0 ) write(6,'(A14,9x,A2)') 'Load CC-File:', 'No'
     IF ( pre_gs == 1 ) write(6,'(A14,8x,A3)') 'Load CC-File:', 'Yes'
     IF ( abs(add_n) > 0 ) write(6,'(A15,6x,I4)') 'Neutron PR/PA:', add_n
     IF ( abs(add_p) > 0 ) write(6,'(A14,7x,I4)') 'Proton PR/PA:', add_p
     IF ( cc_approx >= 0 ) write(6,*) 'CC-File:             ', trim(cc_file)
     write(6,*)
  end IF

  ! Setup Memory Tracking
  CALL initialize_memory
  
  ! Setup Model Space and Two-Body Structures
  IF ( iam == 0 ) write(6,*) '...Setting up Model Space...'
  CALL allocate_sp_data
  CALL setup_sp_data
  CALL setup_structures
  CALL setup_cross_structures
  CALL mpi_barrier(mpi_comm_world, ierror)
  IF ( iam == 0 ) write(6,*) '...Setting up Model Space Done!'
  CALL mem_report('base structures')
  
  ! Setup Interaction
  IF ( iam == 0 ) write(6,*) '...Setting up Bare Interaction...'
  CALL commons_to_angmom
  CALL init_chp_constants
  CALL precalc_chp_functions
  CALL mpi_barrier(mpi_comm_world, ierror)

  
  CALL setup_2b_interaction
  IF ( test > 0 ) CALL setup_tests
  IF ( test == 4 ) then
     test = 0
     test0 = 4
  end IF
  IF ( iam == 0 ) write(6,*) '...Setting up Bare Interaction Done!'
  CALL mem_report('interactions')
  
  
  ! Perform Coupled Cluster Iterations
  IF ( cc_approx == -1 ) then
     CALL mbpt
  else
     CALL ccdt_iter
  end IF
  IF ( add_n == 100 .or. add_p == 100 ) then
     CALL compute_eom_states
  else IF ( add_n > 0 .or. add_p > 0 ) then
     CALL compute_paeom_states
  else IF ( add_n < 0 .or. add_p < 0 ) then
     CALL compute_preom_states
  end IF

  
  endwtime = MPI_WTIME()
  IF ( iam == 0 ) write(6,*)
  IF ( iam == 0 ) write(6,*) 'Total execution time for CCD(T) code', endwtime - startwtime
  IF ( iam == 0 ) write(6,*)

  
  ! Finalize Parallelization and HDF5
  CALL mpi_finalize(ierror)
  CALL h5close_f(error)
  
END PROGRAM ccdt
