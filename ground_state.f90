
SUBROUTINE ccdt_iter
  USE parallel
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE diis_mod
  USE hdf5_wrapper
  USE configurations
  USE build_status
  USE contracts
  USE mpi_check
  USE mem_tracker
  
  USE mpi_check
  IMPLICIT NONE
  INTEGER :: iterations, max_iterations, count
  REAL(dp) :: sigma, tolerance
  COMPLEX(dpc) :: e1, e2, e3
  REAL(dp) :: startwtime, endwtime
  LOGICAL :: selfconsistency

  ! Preconditions: basis, channels, mappings, interactions, and Fock must be built.
  CALL assert_built('basis',        'ccdt_iter')
  CALL assert_built('channels_2b',  'ccdt_iter')
  CALL assert_built('mappings',     'ccdt_iter')
  CALL assert_built('interactions', 'ccdt_iter')
  CALL assert_built('fock',         'ccdt_iter')
  CALL assert_holes_ordered('ccdt_iter')
  
  ! Setup t-amplitudes
  IF ( iam == 0 ) write(6,*) '...Setting up T amplitudes...'
  CALL setup_t_amplitudes
  IF ( cc_approx > 1 .and. .not. pre_gs0 ) then
     CALL setup_t3_amplitudes(.false.)
  end IF
  diis_subspace = 10
  diis_step = 10
  IF ( .not. pre_gs0 ) CALL diis_setup
  IF ( iam == 0 ) write(6,*) '...Setting up T amplitudes done!'
  CALL mem_report('T amplitudes')
  
  
  ! Setup Hbar structures
  IF ( .not. pre_gs0 ) then
     IF ( iam == 0 ) write(6,*) '...Setting up Hbar structures...'
     CALL allocate_hbar_t2_iter
     IF ( iam == 0 ) write(6,*) '...Setting up Hbar structures done!'
     CALL mem_report('interactions + Hbar')
  end IF
  
  
  ! Total Memory
  IF ( iam == 0 ) write(6,*) '...Solving T Amplitudes Ready!'
  CALL mem_report('total')
  
  
  ! Start Solving T-amplitudes
  IF ( pre_gs0 ) then
     IF ( iam == 0 ) write(6,*)
     IF ( iam == 0 ) write(6,*) 'Reading T amplitudes from file... ', trim(cc_file)
     CALL read_gs_hdf5
     CALL mpi_barrier(mpi_comm_world,ierror)

     CALL cc_energy(e2, .true.)
     IF ( iam == 0 ) write(6,*)
     
  ELSE
     ! Here is the main iteration loop
     max_iterations = 100
     tolerance = 1.e-7
     IF ( cc_approx > 1 ) tolerance = 1.e-5
     cc_level = 0.d0
     cc_scale = 0.25d0
     sigma = 1000.d0
     selfconsistency = .FALSE.

     IF ( iam == 0 ) THEN
        write(6,*)
        write(6,*) '=== CCD Iterations ==='
        write(6,*) '  Max iterations:        ', max_iterations
        write(6,*) '  Convergence tolerance: ', tolerance
        write(6,*)
     END IF
     
     CALL cc_energy(e1, .true.)
     
     count = 0
     nstep = 0
     iterations = 1
     t3_switch = .FALSE.
     e3 = 0.d0

     DO WHILE ( .not.selfconsistency .and. iterations <= max_iterations-1 )

        startwtime = MPI_WTIME()
        CALL build_hbar_t2_iter
        endwtime = MPI_WTIME()
        IF ( iam == 0 ) write(6,'(A29,F12.2,A5)') ' ...Computing Hbar...        ', endwtime - startwtime, ' sec.'

        startwtime = MPI_WTIME()
        CALL t2_eqn
        CALL t2_denom
        endwtime = MPI_WTIME()
        IF ( iam == 0 ) write(6,'(A29,F12.2,A5)') ' ...Computing T2...          ', endwtime - startwtime, ' sec.'
        
        count = count + 1
        nstep = nstep + 1
        ! CALL linear_t2
        CALL diis_t2(count)
        CALL t2_cross_recouple
        CALL cc_energy(e2, .false.)
        eccdt = e2
        sigma = abs(e1 - e2)
        e1 = e2

        IF ( iam == 0 ) write(6,'(A14,I4,A8,F16.10,A8,ES11.3)') &
             ' ...Iteration ', iterations, '   E/A =', REAL(e2)/below_ef, '   dE = ', sigma
        IF ( iam == 0 ) write(6,*)
        
        iterations = iterations + 1
        IF ( ( iterations == max_iterations .or. abs(sigma) < tolerance ) .and. &
             iterations > min(max_iterations-1,10) )  then
           selfconsistency = .TRUE.
        ELSE
           selfconsistency = .FALSE.
        end IF
     end DO

     IF ( cc_approx > 0 ) then
        IF ( iam == 0 ) write(6,*)
        IF ( iam == 0 ) write(6,'(A,F16.10,A)') '  CCD Energy   = ', REAL(eccdt), ' MeV'
        IF ( iam == 0 ) write(6,'(A,F16.10,A)') '  CCD Energy/A = ', REAL(eccdt)/below_ef, ' MeV'
        IF ( iam == 0 ) write(6,*)
     end IF
     
     IF ( cc_approx > 1 ) then
        max_iterations = 100
        tolerance = 1.e-6
        cc_level = 0.d0
        cc_scale = 0.25d0
        sigma = 1000.d0
        selfconsistency = .FALSE.
        
        diis_subspace = 10
        diis_step = 10
        CALL diis_take_down
        CALL diis_setup
        
        count = 0
        nstep = 0
        iterations = 1
        t3_switch = .TRUE.

        IF ( iam == 0 ) THEN
           write(6,*)
           write(6,*) '=== CCDT Iterations ==='
           write(6,*) '  Max iterations:        ', max_iterations
           write(6,*) '  Convergence tolerance: ', tolerance
           write(6,*)
        END IF

        DO WHILE ( .not.selfconsistency .and. iterations <= max_iterations-1 )

           startwtime = MPI_WTIME()
           CALL t3_eqn
           CALL t3_denom
           endwtime = MPI_WTIME()
           IF ( iam == 0 ) write(6,'(A29,F12.2,A5)') ' ...Computing T3...          ', endwtime - startwtime, ' sec.'
                      
           startwtime = MPI_WTIME()
           CALL build_hbar_t2_iter
           endwtime = MPI_WTIME()
           IF ( iam == 0 ) write(6,'(A29,F12.2,A5)') ' ...Computing Hbar...        ', endwtime - startwtime, ' sec.'
           
           startwtime = MPI_WTIME()
           CALL t2_eqn
           CALL t2_t3_eqn
           CALL t2_denom
           endwtime = MPI_WTIME()
           IF ( iam == 0 ) write(6,'(A29,F12.2,A5)') ' ...Computing T2...          ', endwtime - startwtime, ' sec.'
           
           count = count + 1
           nstep = nstep + 1
           ! CALL linear_t2
           CALL diis_t2(count)
           CALL t2_cross_recouple
           CALL cc_energy(e2, .false.)
           eccdt = e2
           sigma = abs(e1 - e2)
           e1 = e2
           ! IF ( tnf_switch ) CALL cc_energy_3b(e3) ! only for testing

           IF ( iam == 0 ) write(6,'(A14,I4,A8,F16.10,A8,ES11.3)') &
                ' ...Iteration ', iterations, '   E/A =', REAL(e2)/below_ef, '   dE = ', sigma
           IF ( iam == 0 ) write(6,*)

           iterations = iterations + 1
           IF ( ( iterations == max_iterations .or. abs(sigma) < tolerance ) .and. &
                iterations > min(max_iterations-1,10) )  then
              selfconsistency = .TRUE. 
           ELSE
              selfconsistency = .FALSE.
           end IF
           
        end DO
        IF ( cc_approx > 1 .and. tnf_approx > 1 ) then
           CALL cc_energy_3b(e3, .true.)
           eccdt = e3
        end IF
     end IF
     
     ! Write to file
     IF ( iam == 0 ) write(6,*)
     IF ( iam == 0 ) write(6,'(A34,31x,A)') 'Writing CC amplitudes to file...', trim(cc_file)
     IF ( iam == 0 ) write(6,*)
     CALL print_gs_hdf5
  end IF
  IF ( .not. pre_gs0 ) then
     CALL deallocate_hbar_t2_iter
     hbar_built = .FALSE.
     IF (ALLOCATED(t2_diis)) CALL diis_take_down
     IF ( cc_approx > 1 ) then
        CALL deallocate_t3_amplitudes ! comment out for testing
        t3_switch = .FALSE.
     end IF
  end IF
  IF ( cc_approx == 1 ) then
     t3_switch = .TRUE.
     CALL setup_t3_amplitudes(.false.)
     CALL cc_energy_lambda(e3, .true.)
     eccdt = e3
     IF ( iam == 0 ) write(6,*)
  end IF
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  IF ( iam == 0 ) THEN
     write(6,*)
     write(6,'(A,F16.10,A)') '  Final Energy   = ', REAL(eccdt), ' MeV'
     write(6,'(A,F16.10,A)') '  Final Energy/A = ', REAL(eccdt)/below_ef, ' MeV'
     write(6,*)
  END IF
  
END SUBROUTINE ccdt_iter



!
!   \     /   \     /     \  /   \  /     \     /  \     /   \     /  \     /     \     /  \     /   \     /  \     /   \     /  \     /
!  i\    /a  j\    /b    i\ /a  j\ /b    a\   i/   \    /    \    /   \b   /j    i\   a/   \b   /j  a\   i/   \j   /b  a\   i/   \b   /j
!   \   /     \   /       \/     \/       \   /   b\   /j   a\   /i   \   /       \   /    \   /     \   /    \   /     \   /    \   /  
!   \  /      \  /   <--  ========  +  xxxx  /     \  /   +  \  /     \  xxxx  +  \  =======  /   +  \  xxxxxxx  /   +  \  x:::::x  /   
!   \ /       \ /                        c\ /      \ /       \ /      \ /k        \ /c    d\ /       \ /k    l\ /       \ /k    c\ /    
!   \/        \/                          \/       \/        \/       \/          \/       \/        \/       \/        \/       \/     
!   ------------                          -----------        -----------          -----------        -----------        -----------     
!
! <ab|t|ij> <-- <ab|v|ij> + P(ab).<cb|t|ij>.<a|x|c> - P(ij).<ab|t|ik>.<k|x|j> + (1/2).<ab|v|cd>.<cd|t|ij> + (1/2).<ab|t|kl>.<kl|x|ij> - P(ab|ij).<ac|t|kj>.<kb|x*|ic>
!
SUBROUTINE t2_eqn
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE parallel
  USE configurations
  USE contracts
  
  USE mpi_check
  IMPLICIT NONE
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max
  INTEGER :: a,b,i,j,c,k, bra,ket,bra2,ket2
  INTEGER :: dim1,dim2,dim3, phase
  COMPLEX(dpc) :: v1b

  CALL assert_built('t2',   't2_eqn')
  CALL assert_built('hbar', 't2_eqn')
  
  IF ( iam == 0 ) WRITE(6,*) "...Computing T2 amplitudes..."

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle

     ! <ab|t|ij>
     t2_ccm_eqn(ch)%cval = 0.d0
  end DO

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)
     
     IF ( check_my_channel_v2b_pphh(ch) == 0 ) cycle
     bra_min = mapping_v2b_pphh(iam+1,ch,2)
     bra_max = mapping_v2b_pphh(iam+1,ch,3)

     ! <ab|t|ij> <-- <ab|v|ij>
     t2_ccm_eqn(ch)%cval(bra_min:bra_max,:) = t2_ccm_eqn(ch)%cval(bra_min:bra_max,:) &
          + v2b_pphh(ch)%cval(bra_min:bra_max,:)

     ! <ab|t|ij> <-- + P(ab).<cb|t|ij>.<a|x|c>
     DO bra = bra_min, bra_max
        a   = lookup_2b_configs(3,ch)%ival2(1,bra)
        b   = lookup_2b_configs(3,ch)%ival2(2,bra)
        DO c = below_ef+1, tot_orbs
           IF ( c == a .or. c == b ) cycle
           IF ( ch /= pp_channel_2b%ival2(c,b) ) cycle
           
           phase = 1
           bra2 = pp_config_2b%ival2(c,b)
           IF ( b < c ) phase = -phase
           v1b = hbar1b_I2(a,c)
           
           t2_ccm_eqn(ch)%cval(bra,:) = t2_ccm_eqn(ch)%cval(bra,:) &
                + 2.d0 * phase * t2_ccm(ch)%cval(bra2,:)*v1b
        end DO
     end DO
     
     ! <ab|t|ij> <-- - P(ij).<ab|t|ik>.<k|x|j>
     DO ket = 1, ket_confs
        i   = lookup_2b_configs(1,ch)%ival2(1,ket)
        j   = lookup_2b_configs(1,ch)%ival2(2,ket)
        DO k = 1, below_ef
           IF ( k == i .or. k == j ) cycle
           IF ( ch /= hh_channel_2b%ival2(i,k) ) cycle
           
           phase = 1
           ket2 = hh_config_2b%ival2(i,k)
           IF ( k < i ) phase = -phase
           v1b = hbar1b_I3(k,j)
           
           t2_ccm_eqn(ch)%cval(bra_min:bra_max,ket) = t2_ccm_eqn(ch)%cval(bra_min:bra_max,ket) &
                - 2.d0 * phase * t2_ccm(ch)%cval(bra_min:bra_max,ket2)*v1b
        end DO
     end DO
          
     ! <ab|t|ij> <-- + (1/2).<ab|t|kl>.<kl|x|ij>
     dim1 = bra_max - bra_min + 1
     dim2 = ket_confs
     dim3 = ket_confs
     CALL ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), t2_ccm(ch)%cval(bra_min:bra_max,:), dim1, &
          hbar2b_I4(ch)%cval, dim3, dcmplx(1.d0,0.d0), t2_ccm_eqn(ch)%cval(bra_min:bra_max,:), dim1 )     
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)
     
     IF ( check_my_channel_v2b_pppp(ch) == 0 ) cycle
     bra_min = mapping_v2b_pppp(iam+1,ch,2)
     bra_max = mapping_v2b_pppp(iam+1,ch,3)
  
     ! <ab|t|ij> <-- + (1/2).<ab|v|cd>.<cd|t|ij>
     dim1 = bra_max - bra_min + 1
     dim2 = ket_confs
     dim3 = bra_confs
     CALL ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), v2b_pppp(ch)%cval(bra_min:bra_max,:), dim1, &
          t2_ccm(ch)%cval, dim3, dcmplx(1.d0,0.d0), t2_ccm_eqn(ch)%cval(bra_min:bra_max,:), dim1 )
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(3,ch) == 0 ) cycle
     IF ( number_2bcross(2,ch) == 0 ) cycle
     bra_confs = number_2bcross(3,ch)
     ket_confs = number_2bcross(2,ch)

     ! <a-j|t|i-b>
     t2_ccm_eqn_cross(ch)%cval = 0.d0
     
     IF ( check_my_channel_v2b_phhp_cross(ch) == 0 ) cycle
     bra_min = mapping_v2b_phhp_cross(iam+1,ch,2)
     bra_max = mapping_v2b_phhp_cross(iam+1,ch,3)
     
     ! <a-j|t|i-b> <-- -P(ab|ij).<ac|t|kj>.<kb|x*|ic>
     !                 -P(ab|ij).<a-j|t|k-c>.<k-c|x*|i-b>
     dim1 = bra_max - bra_min + 1
     dim2 = ket_confs
     dim3 = ket_confs
     CALL ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(-1.d0,0.d0), t2_ccm_cross(ch)%cval(bra_min:bra_max,:), dim1, &
          hbar2b_I5e_cross(ch)%cval, dim3, dcmplx(1.d0,0.d0), t2_ccm_eqn_cross(ch)%cval(bra_min:bra_max,:), dim1 )
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)
  ! Allreduce t2_ccm_eqn in "t2_denom"
  
end SUBROUTINE t2_eqn

!
!   \     /   \     /        \     / \     /         \     / \     /
!  i\    /a  j\    /b        \    /  \    /b         \    /  \    /j
!   \   /     \   /         a\  i/  j\   =====      a\  i/  b\   =====
!   \  /      \  /    <--    \  /    \  /   /\   +   \  /    \  /   /\
!   \ /       \ /            \ /     \ /d c| |k      \ /     \ /k c| |l
!   \/        \/             \/      \/    \/        \/      \/    \/
!   ------------            -----------------       -----------------
!
!  <ab|t|ij> <-- + (1/2).P(ab).<cda|t|kji>.<kb|v|cd> - (1/2).P(ij).<bca|t|kli>.<kl|v|jc>
!
SUBROUTINE t2_t3_eqn
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE parallel
  USE contracts
  
  USE mpi_check
  IMPLICIT NONE
  INTEGER :: ch3, ch1,ch2, ket_confs
  INTEGER :: bra0,ket0, bra,ket, bra_min,bra_max
  INTEGER :: a,b,c, i,j,k, c1,k1, dim1,dim2,dim3
  INTEGER :: aind, kind1, cind1, iind
  REAL(dp) :: phase
  COMPLEX(dpc), allocatable :: temp_mtx(:,:), temp_v2b(:,:)

  CALL assert_built('t2', 't2_t3_eqn')
  CALL assert_built('t3', 't2_t3_eqn')
  
  IF ( iam == 0 ) WRITE(6,*) "...Computing T2 <- T3..."  

  !
  ! <ab|t|ij> <-- + (1/2).P(ab).<cda|t|kji>.<kb|v|cd>
  !               - (1/2).P(ab).<cda|t|ijk>.<kb|v|cd>
  DO ch3 = ch3_min, ch3_max
     DO aind = climits_t3(ch3,1), climits_t3(ch3,2)
        a       = clist_t3(ch3)%ival2(aind,1)
        ch1     = clist_t3(ch3)%ival2(aind,2)
        bra_min = mapping_t3(ch3)%ival2(aind,1)
        bra_max = mapping_t3(ch3)%ival2(aind,2)
        IF ( bra_min <= 0 ) CYCLE
        IF ( bra_max < bra_min ) CYCLE
        dim1 = number_2b(2,ch1)
        IF ( dim1 == 0 ) cycle

        ALLOCATE( temp_v2b(bra_min:bra_max, dim1) )
        temp_v2b = 0.d0
        DO bra = bra_min, bra_max
           bra0 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
           temp_v2b(bra,:) = conjg(v2b_pphp(ch1)%cval(bra0,:))
        end DO
        
        DO kind1 = 1, klimit_t3(ch3)
           IF ( .not. ASSOCIATED(t3_ccm(ch3)%val2(aind,kind1)%cval) ) cycle
           k         = klist_t3(ch3)%ival2(kind1,1)
           ch2       = klist_t3(ch3)%ival2(kind1,2)
           ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
           IF ( ket_confs <= 0 ) CYCLE
           
           ! <kb|v|cd>.<cda|t|ijk>
           dim2 = ket_confs
           dim3 = bra_max-bra_min+1
           
           ALLOCATE( temp_mtx(dim1,dim2) )
           temp_mtx = 0.d0
           CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), temp_v2b(bra_min:bra_max,:), dim3, &
                t3_ccm(ch3)%val2(aind,kind1)%cval(bra_min:bra_max,:), dim3, dcmplx(1.d0,0.d0), temp_mtx, dim1 )
           
           DO bra = 1, dim1
              k1  = lookup_2b_configs(2,ch1)%ival2(1,bra)
              b   = lookup_2b_configs(2,ch1)%ival2(2,bra)
              IF ( k1 /= k ) cycle
              IF ( ch2 /= pp_channel_2b%ival2(a,b) ) cycle
              phase = 1
              bra0 = pp_config_2b%ival2(a,b)
              IF ( b < a ) phase = -phase
              DO ket  = 1, dim2
                 ket0 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
                 t2_ccm_eqn(ch2)%cval(bra0,ket0) = t2_ccm_eqn(ch2)%cval(bra0,ket0) - phase * temp_mtx(bra,ket)
              end DO
           end DO
           DEALLOCATE( temp_mtx )
        end DO
        DEALLOCATE( temp_v2b )
     end DO
  end DO
  
  !
  ! <ab|t|ij> <-- - (1/2).P(ij).<bca|t|kli>.<kl|v|jc>
  !               - (1/2).P(ij).<abc|t|kli>.<kl|v|jc>
  DO ch3 = ch3_min, ch3_max
     DO iind = 1, klimit_t3(ch3)
        i         = klist_t3(ch3)%ival2(iind,1)
        ch2       = klist_t3(ch3)%ival2(iind,2)
        ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
        IF ( ket_confs <= 0 ) CYCLE
        dim2 = number_2b(2,ch2)
        dim3 = ket_confs
        IF ( dim2 == 0 ) cycle

        ALLOCATE( temp_v2b(dim2, dim3) )
        temp_v2b = 0.d0
        DO ket = 1, dim3
           ket0 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
           temp_v2b(:,ket) = conjg(v2b_hphh(ch2)%cval(:,ket0))
        end DO
        
        DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
           IF ( .not. ASSOCIATED(t3_ccm(ch3)%val2(cind1,iind)%cval) ) cycle
           c       = clist_t3(ch3)%ival2(cind1,1)
           ch1     = clist_t3(ch3)%ival2(cind1,2)
           bra_min = mapping_t3(ch3)%ival2(cind1,1)
           bra_max = mapping_t3(ch3)%ival2(cind1,2)
           IF ( bra_min <= 0 ) CYCLE
           IF ( bra_max < bra_min ) CYCLE
        
           ! <abc|t|kli>.<kl|v|jc>
           dim1 = bra_max-bra_min+1

           ALLOCATE ( temp_mtx(bra_min:bra_max,dim2) )
           temp_mtx = dcmplx(0.d0,0.d0)
           CALL ZGEMM ( 'n', 't', dim1, dim2, dim3, dcmplx(1.d0,0.d0), t3_ccm(ch3)%val2(cind1,iind)%cval(bra_min:bra_max,:), dim1, &
                temp_v2b, dim2, dcmplx(1.d0,0.d0), temp_mtx, dim1 )
           
           DO ket = 1, dim2
              j   = lookup_2b_configs(2,ch2)%ival2(1,ket)
              c1  = lookup_2b_configs(2,ch2)%ival2(2,ket)
              IF ( c1 /= c ) cycle
              IF ( ch1 /= hh_channel_2b%ival2(i,j) ) cycle
              phase = 1
              ket0 = hh_config_2b%ival2(i,j)
              IF ( j < i ) phase = -phase
              DO bra = bra_min, bra_max
                 bra0 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
                 t2_ccm_eqn(ch1)%cval(bra0,ket0) = t2_ccm_eqn(ch1)%cval(bra0,ket0) - phase * temp_mtx(bra,ket)
              end DO
           end DO
           DEALLOCATE( temp_mtx )           
        end DO
        DEALLOCATE( temp_v2b )
     end DO
  end DO
  
end SUBROUTINE t2_t3_eqn


SUBROUTINE t2_denom
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE ang_mom_functions
  USE contracts
  
  USE mpi_check
  IMPLICIT NONE
  INTEGER :: ch, bra, ket, bra_confs, ket_confs
  INTEGER :: a,b,i,j
  COMPLEX(dpc) :: denom

  CALL assert_built('t2',   't2_denom')
  CALL assert_built('hbar', 't2_denom')

  CALL t2_add_cross
  
  CALL mpi_barrier(mpi_comm_world,ierror)
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,t2_ccm_eqn(ch)%cval,size(t2_ccm_eqn(ch)%cval),mpi_complex16,mpi_sum,mpi_comm_world,ierror)
     CALL check_mpi(ierror, 't2_denom', 'allreduce')
  end DO

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)
     
     DO bra = 1, bra_confs
        a = lookup_2b_configs(3,ch)%ival2(1,bra)
        b = lookup_2b_configs(3,ch)%ival2(2,bra)
        DO ket = 1, ket_confs
           i = lookup_2b_configs(1,ch)%ival2(1,ket)
           j = lookup_2b_configs(1,ch)%ival2(2,ket)

           denom = hbar1b_I3(i,i) + hbar1b_I3(j,j) - hbar1b_I2(a,a) - hbar1b_I2(b,b) - 2.d0*cc_level
           t2_ccm_eqn(ch)%cval(bra,ket) = t2_ccm_eqn(ch)%cval(bra,ket)/denom
           ! IF ( abs(t2_ccm_eqn(ch)%cval(bra,ket)) > 1.e-6 ) write(6,*) 'T2: ', ch,bra,ket, t2_ccm_eqn(ch)%cval(bra,ket)
        end DO
     end DO
  end DO
  
  IF ( test == 3 ) CALL build_t2_eqn_test
    
END SUBROUTINE t2_denom


!
!     \    /   \    /   \    /      a\ /k  b\   j/  \c   /i   j\ /c  i\   a/  \k   /b
!     \   /    \   /    \   /        \/     \   /   \   /      \/     \   /   \   /
!    a\  /i   b\  /j   c\  /k   <--  ========  /    \  /   +   ========  /    \  /
!     \ /      \ /      \ /                d\ /     \ /              l\ /     \ /
!     \/       \/       \/                  \/      \/                \/      \/
!     --------------------                  ----------                ----------
!
! <abc|t|ijk> <-- - P(c/ab|k/ij).<cd|t|ij>.<ab|v|kd> + P(c/ab|k/ij).<ab|t|lk>.<lc|v|ij>
!
SUBROUTINE t3_eqn
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE parallel
  USE t3_diagrams
  USE omp_lib
  USE contracts

  USE mpi_check
  IMPLICIT NONE
  INTEGER :: ch1,ch2,ch3
  INTEGER :: bra, bra0, bra_min,bra_max, bra_confs,ket_confs
  INTEGER :: a,b,c, k, cind1,kind1
  INTEGER :: ch_ab, ch_cb, ch_ac, bra_ab, bra_cb, bra_ac
  INTEGER :: phase_ab, phase_cb, phase_ac
  REAL(dp) :: startwtime, endwtime
  TYPE (superblock_storage) :: pp_t3_temp

  CALL assert_built('t2', 't3_eqn')
  CALL assert_built('t3', 't3_eqn')

  IF ( iam == 0 ) WRITE(6,'(A29)',advance='no') " ...Computing T3...          "

  startwtime = MPI_WTIME()
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  ! <abc|t|ijk> <-- -P(c/ab|k/ij).<cd|t|ij>.<ab|v|kd>
  DO ch3 = ch3_min, ch3_max
     IF (climits_t3(ch3,2) < climits_t3(ch3,1)) CYCLE
     DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
        IF (ASSOCIATED(t3_ccm(ch3)%pack1(cind1)%buf)) THEN
           t3_ccm(ch3)%pack1(cind1)%buf = 0.d0
        END IF
     END DO
  END DO

  ! <abc|t|ijk> <-- -<cd|t|ij>.<ab|v|kd> + <ab|t|lk>.<lc|v|ij>
  ! <abc|t|ijk> <-- +<ad|t|ij>.<cb|v|kd> - <cb|t|lk>.<la|v|ij>
  ! <abc|t|ijk> <-- +<bd|t|ij>.<ac|v|kd> - <ac|t|lk>.<lb|v|ij>
  DO ch3 = ch3_min, ch3_max
     IF (climits_t3(ch3,2) < climits_t3(ch3,1)) CYCLE
     
     ! <abc|t|ijk> <-- -<cd|t|ij>.<ab|v|kd> + <ab|t|lk>.<lc|v|ij>
     DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
        IF ( .not. ASSOCIATED(t3_ccm(ch3)%pack1(cind1)%buf) ) cycle
        c       = clist_t3(ch3)%ival2(cind1,1)
        ch1     = clist_t3(ch3)%ival2(cind1,2)
        bra_min = mapping_t3(ch3)%ival2(cind1,1)
        bra_max = mapping_t3(ch3)%ival2(cind1,2)
        IF ( bra_min <= 0 ) CYCLE
        IF ( bra_max < bra_min ) CYCLE
        
        !$omp parallel default(shared) private(bra,bra0,a,b, phase_ab,bra_ab,ch_ab, kind1,k,ch2)
        !$omp do schedule(dynamic)
        DO bra = bra_min, bra_max
           DO kind1 = 1, klimit_t3(ch3)
              IF ( .not. ASSOCIATED(t3_ccm(ch3)%val2(cind1,kind1)%cval) ) cycle
              IF ( .not. allocated(pp_config_t3(ch3)%ival1(ch1)%ival1) ) CYCLE
              bra0 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
              a    = lookup_2b_configs(3,ch1)%ival2(1,bra0)
              b    = lookup_2b_configs(3,ch1)%ival2(2,bra0)
              IF ( a == c .or. b == c ) cycle
              phase_ab = 1
              bra_ab = bra0
              ch_ab = ch1
           
              k    = klist_t3(ch3)%ival2(kind1,1)
              ch2  = klist_t3(ch3)%ival2(kind1,2)
              ! <abc|t|ijk> <-- -<cd|t|ij>.<ab|v|kd>
              CALL t3_diag1(ch3, ch_ab,ch2, cind1,kind1, c,k, bra,bra_ab, phase_ab)
              ! <abc|t|ijk> <-- +<ab|t|lk>.<lc|v|ij>
              CALL t3_diag2(ch3, ch_ab,ch2, cind1,kind1, c,k, bra,bra_ab, phase_ab)
           end DO
        end DO        
        !$omp end do
        !$omp end parallel
     end DO
     
     ! <abc|t|ijk> <-- +<ad|t|ij>.<cb|v|kd> - <cb|t|lk>.<la|v|ij>
     DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
        IF ( .not. ASSOCIATED(t3_ccm(ch3)%pack1(cind1)%buf) ) cycle
        c       = clist_t3(ch3)%ival2(cind1,1)
        ch1     = clist_t3(ch3)%ival2(cind1,2)
        bra_min = mapping_t3(ch3)%ival2(cind1,1)
        bra_max = mapping_t3(ch3)%ival2(cind1,2)
        IF ( bra_min <= 0 ) CYCLE
        IF ( bra_max < bra_min ) CYCLE
        
        !$omp parallel default(shared) private(bra,bra0,a,b, phase_cb,bra_cb,ch_cb, kind1,k,ch2)
        !$omp do schedule(dynamic)
        DO bra = bra_min, bra_max
           DO kind1 = 1, klimit_t3(ch3)
              IF ( .not. ASSOCIATED(t3_ccm(ch3)%val2(cind1,kind1)%cval) ) cycle
              IF ( .not. allocated(pp_config_t3(ch3)%ival1(ch1)%ival1) ) CYCLE
              bra0 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
              a    = lookup_2b_configs(3,ch1)%ival2(1,bra0)
              b    = lookup_2b_configs(3,ch1)%ival2(2,bra0)
              IF ( a == c .or. b == c ) cycle
              ch_cb = pp_channel_2b%ival2(c,b)
              phase_cb = -1
              bra_cb = pp_config_2b%ival2(c,b)
              IF ( b < c ) phase_cb = -phase_cb
           
              k    = klist_t3(ch3)%ival2(kind1,1)
              ch2  = klist_t3(ch3)%ival2(kind1,2)
              ! <abc|t|ijk> <-- +<ad|t|ij>.<cb|v|kd>
              CALL t3_diag1(ch3, ch_cb,ch2, cind1,kind1, a,k, bra,bra_cb, phase_cb)
              ! <abc|t|ijk> <-- -<cb|t|lk>.<la|v|ij>
              CALL t3_diag2(ch3, ch_cb,ch2, cind1,kind1, a,k, bra,bra_cb, phase_cb)
           end DO
        end DO        
        !$omp end do
        !$omp end parallel
     end DO

     ! <abc|t|ijk> <-- +<bd|t|ij>.<ac|v|kd> - <ac|t|lk>.<lb|v|ij>
     DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
        IF ( .not. ASSOCIATED(t3_ccm(ch3)%pack1(cind1)%buf) ) cycle
        c       = clist_t3(ch3)%ival2(cind1,1)
        ch1     = clist_t3(ch3)%ival2(cind1,2)
        bra_min = mapping_t3(ch3)%ival2(cind1,1)
        bra_max = mapping_t3(ch3)%ival2(cind1,2)
        IF ( bra_min <= 0 ) CYCLE
        IF ( bra_max < bra_min ) CYCLE
        
        !$omp parallel default(shared) private(bra,bra0,a,b, phase_ac,bra_ac,ch_ac, kind1,k,ch2)
        !$omp do schedule(dynamic)
        DO bra = bra_min, bra_max
           DO kind1 = 1, klimit_t3(ch3)
              IF ( .not. ASSOCIATED(t3_ccm(ch3)%val2(cind1,kind1)%cval) ) cycle
              IF ( .not. allocated(pp_config_t3(ch3)%ival1(ch1)%ival1) ) CYCLE
              bra0 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
              a    = lookup_2b_configs(3,ch1)%ival2(1,bra0)
              b    = lookup_2b_configs(3,ch1)%ival2(2,bra0)
              IF ( a == c .or. b == c ) cycle
              ch_ac = pp_channel_2b%ival2(a,c)
              phase_ac = -1
              bra_ac = pp_config_2b%ival2(a,c)
              IF ( c < a ) phase_ac = -phase_ac
              
              k    = klist_t3(ch3)%ival2(kind1,1)
              ch2  = klist_t3(ch3)%ival2(kind1,2)
              ! <abc|t|ijk> <-- +<bd|t|ij>.<ac|v|kd>
              CALL t3_diag1(ch3, ch_ac,ch2, cind1,kind1, b,k, bra,bra_ac, phase_ac)
              ! <abc|t|ijk> <-- -<ac|t|lk>.<lb|v|ij>
              CALL t3_diag2(ch3, ch_ac,ch2, cind1,kind1, b,k, bra,bra_ac, phase_ac)
           end DO           
        end DO        
        !$omp end do
        !$omp end parallel
     end DO
     
  end DO

  CALL mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  IF ( iam == 0 ) WRITE(6,'(A29,A29,f12.6,A5)') bsp29, " ...Computing T3...          ", (endwtime - startwtime), ' sec.'

  IF ( cc_approx > 1 ) CALL antisymmetrize_t3

end SUBROUTINE t3_eqn


SUBROUTINE antisymmetrize_t3
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE parallel
  USE contracts

  USE mpi_check
  IMPLICIT NONE
  INTEGER :: ch3, ch2, ch_i,ch_j
  INTEGER :: ket,ket0,ket1, ket_confs, bra_min,bra_max
  INTEGER :: i,j,k, iind,jind,kind1, cind1
  REAL(dp) :: phase, startwtime, endwtime
  INTEGER, allocatable :: klist_inv(:)
  TYPE (superblock_storage) :: hh_t3_temp
  TYPE (superblock_storage) :: t3_temp

  CALL assert_built('t3', 'antisymmetrize_t3')

  IF ( iam == 0 ) WRITE(6,'(A29)',advance='no') " ...Antisymmetrizing T3...   "

  ALLOCATE( klist_inv(1:below_ef) )
  
  ! (k/ij)
  startwtime = MPI_WTIME()
  CALL mpi_barrier(mpi_comm_world,ierror)
  DO ch3 = ch3_min, ch3_max

     ! Build inverse k map
     klist_inv = 0
     DO kind1 = 1, klimit_t3(ch3)
        k = klist_t3(ch3)%ival2(kind1,1)
        klist_inv(k) = kind1
     end DO

     ! Build inverse hh->hh(t3_cut) map
     ALLOCATE( hh_t3_temp%ival1(channels_2b%number_confs) )
     DO ch2 = 1, channels_2b%number_confs
        ket_confs = number_2b(1,ch2)
        IF ( ket_confs == 0 ) cycle
        ALLOCATE( hh_t3_temp%ival1(ch2)%ival1(ket_confs) )
        hh_t3_temp%ival1(ch2)%ival1 = 0
        
        ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
        IF ( ket_confs <= 0 ) cycle
        DO ket  = 1, ket_confs
           ket0 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
           hh_t3_temp%ival1(ch2)%ival1(ket0) = ket
        end DO
     end DO
        
     DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
        bra_min = mapping_t3(ch3)%ival2(cind1,1)
        bra_max = mapping_t3(ch3)%ival2(cind1,2)
        IF ( bra_min <= 0 ) CYCLE
        IF ( bra_max < bra_min ) CYCLE
        ALLOCATE( t3_temp%val1(1:klimit_t3(ch3)) )
        
        ! Allocate and Fill t3_temp for all k_channels
        DO kind1 = 1, klimit_t3(ch3)
           IF ( .not. ASSOCIATED(t3_ccm(ch3)%val2(cind1,kind1)%cval) ) cycle
           ch2       = klist_t3(ch3)%ival2(kind1,2)
           ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
           IF ( ket_confs <= 0 ) cycle
           ALLOCATE( t3_temp%val1(kind1)%cval(bra_min:bra_max,ket_confs) )
           t3_temp%val1(kind1)%cval = 0.d0
           t3_temp%val1(kind1)%cval(bra_min:bra_max,:) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra_min:bra_max,:)
        end DO
        
        DO kind1 = 1, klimit_t3(ch3)
           IF ( .not. ASSOCIATED(t3_ccm(ch3)%val2(cind1,kind1)%cval) ) cycle
           k         = klist_t3(ch3)%ival2(kind1,1)
           ch2       = klist_t3(ch3)%ival2(kind1,2)
           ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
           IF ( ket_confs <= 0 ) cycle
           
           ! k <-> i
           DO ket  = 1, ket_confs
              ket0 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
              i    = lookup_2b_configs(1,ch2)%ival2(1,ket0)
              j    = lookup_2b_configs(1,ch2)%ival2(2,ket0)
              iind = klist_inv(i)
              IF ( iind == 0 ) cycle
              IF ( .not. allocated(t3_temp%val1(iind)%cval) ) cycle
              ch_i = klist_t3(ch3)%ival2(iind,2)
              IF ( ch_i == 0 ) cycle
              IF ( ch_i /= hh_channel_2b%ival2(k,j) ) cycle
              phase = 1
              ket1 = hh_config_2b%ival2(k,j)
              ket1 = hh_t3_temp%ival1(ch_i)%ival1(ket1)
              IF ( ket1 == 0 ) cycle
              IF ( j < k ) phase = -phase
              
              t3_ccm(ch3)%val2(cind1,kind1)%cval(bra_min:bra_max,ket) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra_min:bra_max,ket) &
                   - phase * t3_temp%val1(iind)%cval(bra_min:bra_max,ket1)
           end DO
           ! k <-> j
           DO ket  = 1, ket_confs
              ket0 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
              i    = lookup_2b_configs(1,ch2)%ival2(1,ket0)
              j    = lookup_2b_configs(1,ch2)%ival2(2,ket0)
              jind = klist_inv(j)
              IF ( jind == 0 ) cycle
              IF ( .not. allocated(t3_temp%val1(jind)%cval) ) cycle
              ch_j = klist_t3(ch3)%ival2(jind,2)
              IF ( ch_j == 0 ) cycle
              IF ( ch_j /= hh_channel_2b%ival2(i,k) ) cycle
              phase = 1
              ket1 = hh_config_2b%ival2(i,k)
              ket1 = hh_t3_temp%ival1(ch_j)%ival1(ket1)
              IF ( ket1 == 0 ) cycle
              IF ( k < i ) phase = -phase

              t3_ccm(ch3)%val2(cind1,kind1)%cval(bra_min:bra_max,ket) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra_min:bra_max,ket) &
                   - phase * t3_temp%val1(jind)%cval(bra_min:bra_max,ket1)
           end DO
        end DO
        
        DO kind1 = 1, klimit_t3(ch3)
           IF ( .not. allocated(t3_temp%val1(kind1)%cval) ) cycle
           DEALLOCATE( t3_temp%val1(kind1)%cval )
        end DO
        DEALLOCATE( t3_temp%val1 )
        
     end DO
     
     DO ch2 = 1, channels_2b%number_confs
        IF ( .not. allocated(hh_t3_temp%ival1(ch2)%ival1) ) cycle
        ket_confs = number_2b(1,ch2)
        IF ( ket_confs == 0 ) cycle
        DEALLOCATE( hh_t3_temp%ival1(ch2)%ival1 )
     end DO
     DEALLOCATE( hh_t3_temp%ival1 )
     
  end DO
  DEALLOCATE( klist_inv )

  CALL mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  IF ( iam == 0 ) WRITE(6,'(A29,A29,f12.6,A5)') bsp29, " ...Antisymmetrizing T3...   ", (endwtime - startwtime), ' sec.'
  
end SUBROUTINE antisymmetrize_t3


!
!     \    /   \    /   \    /        \    /   \    /  \    /       \    /   \    /  \    /
!     \   /    \   /    \   /        a\   /i  b\   /j  \k  /c      a\   /i  b\  j/   \c  /k
!    a\  /i   b\  /j   c\  /k   <--   \  /     \  /    \  xxxx  +   \  /     \  /    \  xxxx
!     \ /      \ /      \ /           \ /      \ /     \ /d         \ /      \ /     \ /l
!     \/       \/       \/            \/       \/      \/           \/       \/      \/
!     --------------------            -------------------           -------------------
!
! <abc|t|ijk> <-- + P(c/ab).<c|v|d>.<abd|t|ijk> - P(k/ij).<abc|t|ijl>.<l|v|k>
!
SUBROUTINE t3_denom
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE ang_mom_functions
  USE contracts
  
  USE mpi_check
  IMPLICIT NONE
  INTEGER :: ch3, ch1,ch2, bra,ket, bra0,ket0
  INTEGER :: ket_confs, bra_min,bra_max
  INTEGER :: a,b,c, i,j,k, cind1,kind1
  COMPLEX(dpc) :: denom

  CALL assert_built('t3',   't3_denom')
  CALL assert_built('hbar', 't3_denom')

  IF ( iam == 0 ) WRITE(6,*) "...Collecting T3..."

  DO ch3 = ch3_min, ch3_max
     DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
        c       = clist_t3(ch3)%ival2(cind1,1)
        ch1     = clist_t3(ch3)%ival2(cind1,2)
        bra_min = mapping_t3(ch3)%ival2(cind1,1)
        bra_max = mapping_t3(ch3)%ival2(cind1,2)
        IF ( bra_min <= 0 ) CYCLE
        IF ( bra_max < bra_min ) CYCLE
        DO kind1 = 1, klimit_t3(ch3)
           k         = klist_t3(ch3)%ival2(kind1,1)
           ch2       = klist_t3(ch3)%ival2(kind1,2)
           ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
           IF ( ket_confs <= 0 ) CYCLE
           IF ( .not. ASSOCIATED(t3_ccm(ch3)%val2(cind1,kind1)%cval) ) cycle
              
           ! <abc|t|ijk> <-- <abc|v|ijk>
           IF ( tnf_approx > 1 ) then
              IF ( ASSOCIATED(t3_ccm0(ch3)%val2(cind1,kind1)%cval) ) then
                 t3_ccm(ch3)%val2(cind1,kind1)%cval(bra_min:bra_max,:) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra_min:bra_max,:) &
                      + t3_ccm0(ch3)%val2(cind1,kind1)%cval(bra_min:bra_max,:)
              end IF
           end IF
           
           DO bra  = bra_min, bra_max
              bra0 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
              a    = lookup_2b_configs(3,ch1)%ival2(1,bra0)
              b    = lookup_2b_configs(3,ch1)%ival2(2,bra0)
              IF ( a == c .or. b == c ) cycle
              DO ket  = 1, ket_confs
                 ket0 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
                 i    = lookup_2b_configs(1,ch2)%ival2(1,ket0)
                 j    = lookup_2b_configs(1,ch2)%ival2(2,ket0)
                 IF ( i == k .or. j == k ) cycle
                 
                 IF ( tnf_approx > 1 ) then
                    denom = ( hbar1b_I3(i,i) + hbar1b_I3(j,j) + hbar1b_I3(k,k) &
                         - hbar1b_I2(a,a) - hbar1b_I2(b,b) - hbar1b_I2(c,c) ) - 3.d0*cc_level
                 ELSE
                    denom = ( fock_mtx(i,i) + fock_mtx(j,j) + fock_mtx(k,k) &
                         - fock_mtx(a,a) - fock_mtx(b,b) - fock_mtx(c,c) ) - 3.d0*cc_level
                 end IF
                 t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,ket) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,ket)/denom
              end DO
           end DO
        end DO
     end DO
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  IF ( test == 3 ) CALL build_t3_eqn_test
  
end SUBROUTINE t3_denom

