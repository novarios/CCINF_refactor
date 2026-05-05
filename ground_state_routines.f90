
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE setup_t_amplitudes
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE ang_mom_functions
  USE build_status
  USE contracts
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, bra_confs,ket_confs, bra,ket, ndim
  INTEGER :: a,b,i,j
  COMPLEX(dpc) :: matel, denom

  CALL assert_built('fock',         'setup_t_amplitudes')
  CALL assert_built('interactions', 'setup_t_amplitudes')
  CALL assert_holes_ordered('setup_t_amplitudes')
  
  ! setup t2-amplitude
  ALLOCATE( t2_ccm(channels_2b%number_confs) )
  ALLOCATE( t2_ccm_eqn(channels_2b%number_confs) )

  ndim = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)

     ALLOCATE( t2_ccm(ch)%cval(bra_confs,ket_confs) )
     ALLOCATE( t2_ccm_eqn(ch)%cval(bra_confs,ket_confs) )
     t2_ccm(ch)%cval = 0.d0
     t2_ccm_eqn(ch)%cval = 0.d0
     CALL mem_register('t2', REAL(2 * bra_confs*ket_confs * 16.d0, dp))
     
     DO bra = 1, bra_confs
        a = lookup_2b_configs(3,ch)%ival2(1,bra)
        b = lookup_2b_configs(3,ch)%ival2(2,bra)
        DO ket = 1, ket_confs
           i = lookup_2b_configs(1,ch)%ival2(1,ket)
           j = lookup_2b_configs(1,ch)%ival2(2,ket)
           
           matel = v2b_pphh(ch)%cval(bra,ket)
           denom = fock_mtx(i,i) + fock_mtx(j,j) - fock_mtx(a,a) - fock_mtx(b,b)
           t2_ccm(ch)%cval(bra,ket) = matel / denom

           ndim = ndim + 1
        end DO
     end DO
  end DO
  IF ( iam == 0 ) write(6,'(A14,31x,I14)') 'T2 dimension:', ndim

  ! setup t2-cross-amplitude
  ALLOCATE( t2_ccm_cross(channels_2bcross%number_confs) )
  ALLOCATE( t2_ccm_eqn_cross(channels_2bcross%number_confs) )

  ndim = 0
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(3,ch) == 0 ) cycle
     IF ( number_2bcross(2,ch) == 0 ) cycle
     bra_confs = number_2bcross(3,ch)
     ket_confs = number_2bcross(2,ch)

     ALLOCATE( t2_ccm_cross(ch)%cval(bra_confs,ket_confs) )
     ALLOCATE( t2_ccm_eqn_cross(ch)%cval(bra_confs,ket_confs) )
     t2_ccm_cross(ch)%cval = 0.d0
     t2_ccm_eqn_cross(ch)%cval = 0.d0
     CALL mem_register('t2', REAL(2 * bra_confs*ket_confs * 16.d0, dp))
     
     ndim = ndim + bra_confs*ket_confs
  end DO
  IF ( iam == 0 ) write(6,'(A20,25x,I14)') 'T2-cross dimension:', ndim
  ! [old memory print removed - replaced by mem_report]
  IF ( iam == 0 ) write(6,*)
  
  CALL t2_cross_recouple

  IF ( test > 1 ) CALL build_tamp_test

  t2_built = .TRUE.
  
END SUBROUTINE setup_t_amplitudes

SUBROUTINE setup_t3_amplitudes(fill)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE ang_mom_functions
  USE chiral_potentials
  USE build_status
  USE contracts
  use iso_fortran_env, only: error_unit
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: fill
  INTEGER :: number_channels, ch3, ch1,ch2
  INTEGER :: bra_min,bra_max, k1,k2,k3,k4
  INTEGER :: a,b,c,i,j,k, p,q,r, cind1, kind1
  INTEGER :: bra,ket, bra0,ket0, nrow
  INTEGER :: bra_confs,ket_confs
  INTEGER :: c_count, k_count
  INTEGER :: nx3,ny3,nz3,tz3
  INTEGER :: num_ch0, chmin, chmax, ch_pr_proc
  INTEGER(i8) :: ndim3, total, offset, nelem
  COMPLEX(dpc) :: v3b
  INTEGER :: istat

  CALL assert_built('t2',           'setup_t3_amplitudes')
  CALL assert_built('interactions', 'setup_t3_amplitudes')
  CALL assert_holes_ordered('setup_t3_amplitudes')
  
  cut_3b = t3_cut
  
  n3min = 1000
  n3max = -1000
  t3min = 1000
  t3max = -1000
  DO p = below_ef+1, tot_orbs-2
     DO q = p+1, tot_orbs-1
        DO r = q+1, tot_orbs
           IF ( all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r) > n3max ) &
                n3max = all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r)
           IF ( all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r) < n3min ) &
                n3min = all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r)
           IF ( all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r) > t3max ) &
                t3max = all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r)
           IF ( all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r) < t3min ) &
                t3min = all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r)
        end DO
     end DO
  end DO
  
  DO p = 1, below_ef-2
     DO q = p+1, below_ef-1
        DO r = q+1, below_ef
           IF ( all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r) > n3max ) &
                n3max = all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r)
           IF ( all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r) < n3min ) &
                n3min = all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r)
           IF ( all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r) > t3max ) &
                t3max = all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r)
           IF ( all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r) < t3min ) &
                t3min = all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r)
        end DO
     end DO
  end DO

  num_ch0 = ((t3max-t3min)/2 + 1)*(n3max-n3min + 1)**3
  ch_pr_proc = int((num_ch0 + (num_procs-1))/num_procs)
  ALLOCATE( num0_t3(1:2, n3min:n3max, n3min:n3max, n3min:n3max, t3min:t3max) )
  CALL mem_register('t3', REAL(num_ch0 * 4.d0, dp))
  num0_t3 = 0

  chmin = iam*ch_pr_proc + 1
  chmax = min(num_ch0, (iam+1)*ch_pr_proc)
  
  num_ch0 = 0
  number_channels = 0
  DO nx3 = n3min, n3max
     DO ny3 = n3min, n3max
        DO nz3 = n3min, n3max
           DO tz3 = t3min, t3max, 2
              
              num_ch0 = num_ch0 + 1
              IF ( num_ch0 < chmin .or. chmax < num_ch0 ) cycle
              
              CALL number_3b_confs(nx3, ny3, nz3, tz3, 2, ket_confs) ! 2 for hhh, cut here
              IF ( ket_confs <= 0 ) cycle
              CALL number_3b_confs(nx3, ny3, nz3, tz3, 1, bra_confs) ! 1 for ppp, cut here
              IF ( bra_confs <= 0 ) cycle
              num0_t3(1,nx3,ny3,nz3,tz3) = bra_confs
              num0_t3(2,nx3,ny3,nz3,tz3) = ket_confs
              number_channels = number_channels + 1
           end DO
        end DO
     end DO
  end DO
  CALL mpi_allreduce(mpi_in_place,number_channels,1,mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_t3_amplitudes', 'allreduce')
  CALL mpi_allreduce(mpi_in_place,num0_t3,size(num0_t3),mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_t3_amplitudes', 'allreduce')
  channels_t3%number_confs = number_channels
  
  ALLOCATE( channels_t3%config_NxNyNz_Tz(4*number_channels) )
  CALL mem_register('t3', REAL(4 * number_channels * 4.d0, dp))

  number_channels = 0
  ndim3 = 0
  DO nx3 = n3min, n3max
     DO ny3 = n3min, n3max
        DO nz3 = n3min, n3max
           DO tz3 = t3min, t3max, 2

              bra_confs = num0_t3(1,nx3,ny3,nz3,tz3)
              IF ( bra_confs <= 0 ) cycle
              ket_confs = num0_t3(2,nx3,ny3,nz3,tz3)
              IF ( ket_confs <= 0 ) cycle
              
              number_channels = number_channels + 1
              ch3 = number_channels

              ndim3 = ndim3 + int(bra_confs,8)*int(ket_confs,8)

              k1 = ch3*4 - 3
              k2 = ch3*4 - 2
              k3 = ch3*4 - 1
              k4 = ch3*4
              channels_t3%config_NxNyNz_Tz(k1) = nx3
              channels_t3%config_NxNyNz_Tz(k2) = ny3
              channels_t3%config_NxNyNz_Tz(k3) = nz3
              channels_t3%config_NxNyNz_Tz(k4) = tz3
           end DO
        end DO
     end DO
  end DO
  IF ( iam == 0 ) write(6,'(A14,31x,I17)') 'T3 dimension:', ndim3

  IF ( ndim3 == 0 ) then
     cc_approx = 0
     IF ( iam == 0 ) write(6,*) '...T3 cut is too small, reverting to CCD...'
     return
  end IF

  
  ! Setup T3 mapping
  CALL setup_proc_mappings0_t3

  
  ! Find c+ab combinations
  ALLOCATE( clist_t3(number_channels) )
  ALLOCATE( climit_t3(number_channels) ) ! max c
  CALL mem_register('t3', REAL(number_channels * 4.d0, dp))
  DO ch3 = 1, number_channels
     nx3 = channels_t3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_t3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_t3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_t3%config_NxNyNz_Tz(ch3*4)     
     
     c_count = 0
     DO c = below_ef+1, tot_orbs
        ch1 = get_ch_ch3(nx3,ny3,nz3,tz3,1,c)
        IF ( ch1 == 0 ) cycle
        c_count = c_count + 1
     end DO
     climit_t3(ch3) = c_count
     IF ( c_count == 0 ) then
        write(error_unit,*) "t3_c_count: ch3=", ch3, ", c=", c
        error stop
     end IF
     
     ALLOCATE( clist_t3(ch3)%ival2(c_count,2) ) ! c, ch1
     clist_t3(ch3)%ival2 = 0
     CALL mem_register('t3', REAL(2 * c_count * 4.d0, dp))
     
     c_count = 0
     DO c = below_ef+1, tot_orbs
        ch1 = get_ch_ch3(nx3,ny3,nz3,tz3,1,c)
        IF ( ch1 == 0 ) cycle
        c_count = c_count + 1
        clist_t3(ch3)%ival2(c_count, 1) = c
        clist_t3(ch3)%ival2(c_count, 2) = ch1
     end DO
  end DO

  ! Find k+ij combinations
  ALLOCATE( klist_t3(number_channels) )
  ALLOCATE( klimit_t3(number_channels) ) ! max k
  CALL mem_register('t3', REAL(number_channels * 4.d0, dp))
  DO ch3 = 1, number_channels
     nx3 = channels_t3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_t3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_t3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_t3%config_NxNyNz_Tz(ch3*4)     
     
     k_count = 0
     DO k = 1, below_ef
        ch2 = get_ch_ch3(nx3,ny3,nz3,tz3,2,k)
        IF ( ch2 == 0 ) cycle
        k_count = k_count + 1
     end DO
     klimit_t3(ch3) = k_count
     IF ( k_count == 0 ) then
        write(error_unit,*) "t3_k_count: ch3=", ch3, ", k=", k
        error stop
     end IF
     
     ALLOCATE( klist_t3(ch3)%ival2(k_count,2) ) ! k, ch2
     klist_t3(ch3)%ival2 = 0
     CALL mem_register('t3', REAL(2 * k_count * 4.d0, dp))
     
     k_count = 0
     DO k = 1, below_ef
        ch2 = get_ch_ch3(nx3,ny3,nz3,tz3,2,k)
        IF ( ch2 == 0 ) cycle
        k_count = k_count + 1
        klist_t3(ch3)%ival2(k_count, 1) = k
        klist_t3(ch3)%ival2(k_count, 2) = ch2
     end DO
  end DO

  ! Find 2b -> 2b(t3_cut) translations
  ALLOCATE( number_2b_t3(ch3_min:ch3_max) )
  ALLOCATE( pp_config_t3(ch3_min:ch3_max) )
  ALLOCATE( hh_config_t3(ch3_min:ch3_max) )  
  DO ch3 = ch3_min, ch3_max
     ALLOCATE( number_2b_t3(ch3)%ival2(1:2,channels_2b%number_confs) )
     ALLOCATE( pp_config_t3(ch3)%ival1(channels_2b%number_confs) )
     ALLOCATE( hh_config_t3(ch3)%ival1(channels_2b%number_confs) )
     number_2b_t3(ch3)%ival2= 0

     ! c+ab
     DO cind1 = 1, climit_t3(ch3)
        c    = clist_t3(ch3)%ival2(cind1,1)
        ch1  = clist_t3(ch3)%ival2(cind1,2)
        
        bra = 0 ! count c+ab(t3_cut) configs
        bra_confs = number_2b(3,ch1)
        DO bra0 = 1, bra_confs
           a = lookup_2b_configs(3,ch1)%ival2(1,bra0)
           b = lookup_2b_configs(3,ch1)%ival2(2,bra0)
           IF ( cut3b(a,b,c) ) cycle
           bra = bra + 1
        end DO
        number_2b_t3(ch3)%ival2(1,ch1) = bra
        IF ( .not. ALLOCATED( pp_config_t3(ch3)%ival1(ch1)%ival1 ) .and. bra > 0 ) then
           ALLOCATE( pp_config_t3(ch3)%ival1(ch1)%ival1(bra) )
           pp_config_t3(ch3)%ival1(ch1)%ival1 = 0
           bra = 0
           DO bra0 = 1, bra_confs
              a = lookup_2b_configs(3,ch1)%ival2(1,bra0)
              b = lookup_2b_configs(3,ch1)%ival2(2,bra0)
              IF ( cut3b(a,b,c) ) cycle
              bra = bra + 1
              pp_config_t3(ch3)%ival1(ch1)%ival1(bra) = bra0
           end DO
        end IF
     end DO
     
     ! k+ij
     DO kind1 = 1, klimit_t3(ch3)
        k    = klist_t3(ch3)%ival2(kind1,1)
        ch2  = klist_t3(ch3)%ival2(kind1,2)
        
        ket = 0 ! count k+ij(t3_cut) configs
        ket_confs = number_2b(1,ch2)
        DO ket0 = 1, ket_confs
           i = lookup_2b_configs(1,ch2)%ival2(1,ket0)
           j = lookup_2b_configs(1,ch2)%ival2(2,ket0)
           IF ( cut3b(i,j,k) ) cycle
           ket = ket + 1
        end DO
        number_2b_t3(ch3)%ival2(2,ch2) = ket
        IF ( .not. ALLOCATED( hh_config_t3(ch3)%ival1(ch2)%ival1 ) .and. ket > 0 ) then
           ALLOCATE( hh_config_t3(ch3)%ival1(ch2)%ival1(ket) )
           hh_config_t3(ch3)%ival1(ch2)%ival1 = 0
           ket = 0
           DO ket0 = 1, ket_confs
              i = lookup_2b_configs(1,ch2)%ival2(1,ket0)
              j = lookup_2b_configs(1,ch2)%ival2(2,ket0)
              IF ( cut3b(i,j,k) ) cycle
              ket = ket + 1
              hh_config_t3(ch3)%ival1(ch2)%ival1(ket) = ket0
           end DO
        end IF
     end DO
     
  end DO

  
  ! Setup T3 mapping
  CALL setup_proc_mappings1_t3
  CALL mem_report('T3 structures')



  ALLOCATE( t3_ccm(ch3_min:ch3_max) )
  IF ( tnf_approx > 1 ) ALLOCATE( t3_ccm0(ch3_min:ch3_max) )
  DO ch3 = ch3_min, ch3_max
     
     IF (climits_t3(ch3,2) < climits_t3(ch3,1)) CYCLE
     ALLOCATE( t3_ccm(ch3)%val2(climits_t3(ch3,1):climits_t3(ch3,2), klimit_t3(ch3)) )
     ALLOCATE( t3_ccm(ch3)%pack1(climits_t3(ch3,1):climits_t3(ch3,2)) )
     
     DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
        c       = clist_t3(ch3)%ival2(cind1,1)
        ch1     = clist_t3(ch3)%ival2(cind1,2)
        IF ( c <= 0 .or. ch1 <= 0 ) cycle
        bra_min = mapping_t3(ch3)%ival2(cind1,1)
        bra_max = mapping_t3(ch3)%ival2(cind1,2)
        IF ( bra_min <= 0 ) cycle
        IF ( bra_min > bra_max ) cycle
        nrow = bra_max - bra_min + 1

        total = 0
        DO kind1 = 1, klimit_t3(ch3)
           k   = klist_t3(ch3)%ival2(kind1,1)
           ch2 = klist_t3(ch3)%ival2(kind1,2)
           IF ( k <= 0 .or. ch2 <= 0 ) cycle
           ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
           IF ( ket_confs <= 0 ) cycle
           total = total + int(nrow,8) * int(ket_confs,8)
        END DO
        IF ( total <= 0 ) cycle

        ALLOCATE( t3_ccm(ch3)%pack1(cind1)%buf(total), stat=istat )
        IF ( istat /= 0 ) THEN
           WRITE(error_unit,'(A,I6,A,I8,A,I8,A,I14,A)') &
                'RANK ', iam, ': ALLOC FAILED for T3 buf, ch3=', ch3, &
                ' cind1=', cind1, ' requested=', total*16, ' bytes'
           CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierror)
        END IF
        t3_ccm(ch3)%pack1(cind1)%buf = 0.d0

        offset = 1
        DO kind1 = 1, klimit_t3(ch3)
           k   = klist_t3(ch3)%ival2(kind1,1)
           ch2 = klist_t3(ch3)%ival2(kind1,2)
           IF ( k <= 0 .or. ch2 <= 0 ) cycle
           ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
           IF ( ket_confs <= 0 ) cycle
           nelem = int(nrow,8) * int(ket_confs,8)
           
           t3_ccm(ch3)%val2(cind1,kind1)%cval(bra_min:bra_max,1:ket_confs) => &
                t3_ccm(ch3)%pack1(cind1)%buf(offset : offset + nelem - 1)

           CALL mem_register('t3', REAL(nrow*ket_confs * 16.d0, dp))
           
           offset = offset + nelem
        END DO

        if (offset /= total + 1_8) then
           write(error_unit,*) 'PACK OFFSET MISMATCH', ' ch3=',ch3,' cind1=',cind1, &
                ' offset=',offset,' total=',total
           error stop
        end if

     end DO
  end DO

  ! Allocate and fill t3_ccm0 (3-body matrix elements) for tnf_approx > 1
  IF ( tnf_approx > 1 ) THEN
     DO ch3 = ch3_min, ch3_max
        IF (climits_t3(ch3,2) < climits_t3(ch3,1)) CYCLE
        ALLOCATE( t3_ccm0(ch3)%val2(climits_t3(ch3,1):climits_t3(ch3,2), klimit_t3(ch3)) )
        ALLOCATE( t3_ccm0(ch3)%pack1(climits_t3(ch3,1):climits_t3(ch3,2)) )

        DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
           c       = clist_t3(ch3)%ival2(cind1,1)
           ch1     = clist_t3(ch3)%ival2(cind1,2)
           IF ( c <= 0 .or. ch1 <= 0 ) cycle
           bra_min = mapping_t3(ch3)%ival2(cind1,1)
           bra_max = mapping_t3(ch3)%ival2(cind1,2)
           IF ( bra_min <= 0 ) cycle
           IF ( bra_min > bra_max ) cycle
           nrow = bra_max - bra_min + 1

           total = 0
           DO kind1 = 1, klimit_t3(ch3)
              k   = klist_t3(ch3)%ival2(kind1,1)
              ch2 = klist_t3(ch3)%ival2(kind1,2)
              IF ( k <= 0 .or. ch2 <= 0 ) cycle
              ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
              IF ( ket_confs <= 0 ) cycle
              total = total + int(nrow,8) * int(ket_confs,8)
           END DO
           IF ( total <= 0 ) cycle

           ALLOCATE( t3_ccm0(ch3)%pack1(cind1)%buf(total), stat=istat )
           IF ( istat /= 0 ) THEN
              WRITE(error_unit,'(A,I6,A,I8,A,I8,A,I14,A)') &
                   'RANK ', iam, ': ALLOC FAILED for T3_CCM0 buf, ch3=', ch3, &
                   ' cind1=', cind1, ' requested=', total*16, ' bytes'
              CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierror)
           END IF
           t3_ccm0(ch3)%pack1(cind1)%buf = 0.d0

           offset = 1
           DO kind1 = 1, klimit_t3(ch3)
              k   = klist_t3(ch3)%ival2(kind1,1)
              ch2 = klist_t3(ch3)%ival2(kind1,2)
              IF ( k <= 0 .or. ch2 <= 0 ) cycle
              ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
              IF ( ket_confs <= 0 ) cycle
              nelem = int(nrow,8) * int(ket_confs,8)

              t3_ccm0(ch3)%val2(cind1,kind1)%cval(bra_min:bra_max,1:ket_confs) => &
                   t3_ccm0(ch3)%pack1(cind1)%buf(offset : offset + nelem - 1)

              CALL mem_register('t3_v3b', REAL(nrow*ket_confs * 16.d0, dp))

              ! Fill with 3-body matrix elements
              !$omp parallel default(shared) private(bra,ket,bra0,ket0, a,b,i,j, v3b)
              !$omp do schedule(dynamic)
              DO bra = bra_min, bra_max
                 DO ket = 1, ket_confs
                    bra0 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
                    a = lookup_2b_configs(3,ch1)%ival2(1,bra0)
                    b = lookup_2b_configs(3,ch1)%ival2(2,bra0)
                    IF ( a == c .or. b == c ) cycle
                    ket0 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
                    i   = lookup_2b_configs(1,ch2)%ival2(1,ket0)
                    j   = lookup_2b_configs(1,ch2)%ival2(2,ket0)
                    IF ( i == k .or. j == k ) cycle
                    v3b = v3int(a,b,c,i,j,k)
                    t3_ccm0(ch3)%val2(cind1,kind1)%cval(bra,ket) = v3b
                 end DO
              end DO
              !$omp end do
              !$omp end parallel

              offset = offset + nelem
           END DO
        end DO
     end DO
     CALL mem_report('T3 3-body matrix elements')
  END IF

  CALL mem_report('T3 amplitudes')
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  IF ( fill ) then
     CALL t3_eqn
     CALL t3_denom
  end IF
  IF ( test > 1 ) CALL build_tamp_test

  t3_built = .TRUE.
  
end SUBROUTINE setup_t3_amplitudes

SUBROUTINE deallocate_t3_amplitudes
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE ang_mom_functions
  USE chiral_potentials
  USE build_status
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch3, cind1,kind1
    
  DO ch3 = ch3_min, ch3_max
     DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
        IF ( ASSOCIATED(t3_ccm(ch3)%pack1(cind1)%buf) ) then
           DEALLOCATE( t3_ccm(ch3)%pack1(cind1)%buf )
        end IF
        DO kind1 = 1, klimit_t3(ch3)
           IF ( ASSOCIATED( t3_ccm(ch3)%val2(cind1,kind1)%cval) ) then
              NULLIFY( t3_ccm(ch3)%val2(cind1,kind1)%cval )
           end IF
        end DO
     end DO
  end DO
  DEALLOCATE( t3_ccm )
  IF ( tnf_approx > 1 ) THEN
     DO ch3 = ch3_min, ch3_max
        IF (climits_t3(ch3,2) < climits_t3(ch3,1)) CYCLE
        DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
           IF ( ASSOCIATED(t3_ccm0(ch3)%pack1(cind1)%buf) ) then
              DEALLOCATE( t3_ccm0(ch3)%pack1(cind1)%buf )
           end IF
           DO kind1 = 1, klimit_t3(ch3)
              IF ( ASSOCIATED( t3_ccm0(ch3)%val2(cind1,kind1)%cval) ) then
                 NULLIFY( t3_ccm0(ch3)%val2(cind1,kind1)%cval )
              end IF
           end DO
        end DO
     end DO
     DEALLOCATE( t3_ccm0 )
  END IF
  
  DO ch3 = ch3_min, ch3_max
     DEALLOCATE( clist_t3(ch3)%ival2 )
     DEALLOCATE( klist_t3(ch3)%ival2 )
  end DO
  DEALLOCATE( clist_t3, klist_t3 )
  
  DO ch3 = 1, channels_t3%number_confs
     IF ( ch3_min <= ch3 .and. ch3 <= ch3_max ) DEALLOCATE( mapping_t3(ch3)%ival2 )
  end DO
  DEALLOCATE( mapping_t3 )

  DEALLOCATE( climit_t3, klimit_t3 )
  DEALLOCATE( climits_t3 )
  DEALLOCATE( channels_t3%config_NxNyNz_Tz )

  t3_built = .FALSE.
  
end SUBROUTINE deallocate_t3_amplitudes



SUBROUTINE setup_t3_mbpt
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE ang_mom_functions
  USE chiral_potentials
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: number_channels, ch3, ch1,ch2
  INTEGER :: bra_min,bra_max, k1,k2,k3,k4
  INTEGER :: a,b,c,i,j,k, p,q,r, cind1, kind1
  INTEGER :: bra,ket, bra0,ket0
  INTEGER :: bra_confs,ket_confs
  INTEGER :: bra_confs2,ket_confs2,nconfs
  INTEGER :: c_count, k_count
  INTEGER :: nx3,ny3,nz3,tz3
  INTEGER :: nxp, nyp, nzp, tzp
  INTEGER :: Nx1, Ny1, Nz1, Tz1
  INTEGER(i8) :: ndim3
  COMPLEX(dpc) :: v3b

  cut_3b = t3_cut

  n3min = 1000
  n3max = -1000
  t3min = 1000
  t3max = -1000
  DO p = below_ef+1, tot_orbs-2
     DO q = p+1, tot_orbs-1
        DO r = q+1, tot_orbs
           IF ( all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r) > n3max ) &
                n3max = all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r)
           IF ( all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r) < n3min ) &
                n3min = all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r)
           IF ( all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r) > t3max ) &
                t3max = all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r)
           IF ( all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r) < t3min ) &
                t3min = all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r)
        end DO
     end DO
  end DO
  
  DO p = 1, below_ef-2
     DO q = p+1, below_ef-1
        DO r = q+1, below_ef
           IF ( all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r) > n3max ) &
                n3max = all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r)
           IF ( all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r) < n3min ) &
                n3min = all_orbit%nx(p) + all_orbit%nx(q) + all_orbit%nx(r)
           IF ( all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r) > t3max ) &
                t3max = all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r)
           IF ( all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r) < t3min ) &
                t3min = all_orbit%tz(p) + all_orbit%tz(q) + all_orbit%tz(r)
        end DO
     end DO
  end DO

  number_channels = 0
  DO nx3 = n3min, n3max
     DO ny3 = n3min, n3max
        DO nz3 = n3min, n3max
           DO tz3 = t3min, t3max, 2
              CALL number_3b_red_confs(nx3, ny3, nz3, tz3, 2, ket_confs) ! 2 for hhh
              IF ( ket_confs <= 0 ) cycle
              CALL number_3b_red_confs(nx3, ny3, nz3, tz3, 1, bra_confs) ! 1 for ppp
              IF ( bra_confs <= 0 ) cycle
              number_channels = number_channels + 1
           end DO
        end DO
     end DO
  end DO
  channels_t3%number_confs = number_channels

  ALLOCATE( channels_t3%config_NxNyNz_Tz(4*number_channels) )
  CALL mem_register('t3', REAL(4 * number_channels * 4.d0, dp))
  ALLOCATE( number_3b_t3(number_channels, 1:2) )
  CALL mem_register('t3', REAL(2 * number_channels * 4.d0, dp))
  
  number_channels = 0
  ndim3 = 0
  DO nx3 = n3min, n3max
     DO ny3 = n3min, n3max
        DO nz3 = n3min, n3max
           DO tz3 = t3min, t3max, 2
              CALL number_3b_red_confs(nx3, ny3, nz3, tz3, 2, ket_confs) ! 2 for hhh
              IF ( ket_confs <= 0 ) cycle
              CALL number_3b_red_confs(nx3, ny3, nz3, tz3, 1, bra_confs) ! 1 for ppp
              IF ( bra_confs <= 0 ) cycle

              number_channels = number_channels + 1
              ch3 = number_channels
              ! write(6,*) '$$$ ', ch3, bra_confs, ket_confs
              ndim3 = ndim3 + int(bra_confs,8)*int(ket_confs,8)

              k1 = ch3*4 - 3
              k2 = ch3*4 - 2
              k3 = ch3*4 - 1
              k4 = ch3*4
              channels_t3%config_NxNyNz_Tz(k1) = nx3
              channels_t3%config_NxNyNz_Tz(k2) = ny3
              channels_t3%config_NxNyNz_Tz(k3) = nz3
              channels_t3%config_NxNyNz_Tz(k4) = tz3
              number_3b_t3(ch3, 1) = bra_confs
              number_3b_t3(ch3, 2) = ket_confs
           end DO
        end DO
     end DO
  end DO
  IF ( iam == 0 ) write(6,'(A14,31x,I17)') 'T3 dimension:', ndim3

  IF ( ndim3 == 0 ) then
     cc_approx = 0
     IF ( iam == 0 ) write(6,*) '...T3 cut is too small, reverting to CCD...'
     return
  end IF
  
  ! Setup T3 mapping
  CALL setup_proc_mappings_t3_red

  ! Setup T3 configs
  ALLOCATE( lookup_t3_configs(1:2,ch3_min:ch3_max) )
  DO ch3 = ch3_min, ch3_max
     nx3 = channels_t3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_t3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_t3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_t3%config_NxNyNz_Tz(ch3*4)

     bra_confs = number_3b_t3(ch3, 1)
     ket_confs = number_3b_t3(ch3, 2)

     ALLOCATE( lookup_t3_configs(1,ch3)%ival2(3,bra_confs) )
     lookup_t3_configs(1,ch3)%ival2 = 0
     CALL mem_register('t3', REAL(3 * bra_confs * 4.d0, dp))
     ALLOCATE( lookup_t3_configs(2,ch3)%ival2(3,ket_confs) )
     lookup_t3_configs(2,ch3)%ival2 = 0
     CALL mem_register('t3', REAL(3 * ket_confs * 4.d0, dp))

     ! ppp
     nconfs = 0
     DO p = below_ef+1, tot_orbs
        nxp = all_orbit%nx(p)
        nyp = all_orbit%ny(p)
        nzp = all_orbit%nz(p)
        tzp = all_orbit%tz(p)
        DO ch1 = 1, channels_2b%number_confs
           IF ( number_2b(3,ch1) == 0 ) cycle
           Nx1 = channels_2b%config_NxNyNz_Tz(4*ch1-3)
           Ny1 = channels_2b%config_NxNyNz_Tz(4*ch1-2)
           Nz1 = channels_2b%config_NxNyNz_Tz(4*ch1-1)
           Tz1 = channels_2b%config_NxNyNz_Tz(4*ch1)
           IF ( Nx1 + nxp /= nx3 ) cycle
           IF ( Ny1 + nyp /= ny3 ) cycle
           IF ( Nz1 + nzp /= nz3 ) cycle
           IF ( 2*Tz1 + tzp /= tz3 ) cycle
           bra_confs2 = number_2b(3,ch1)
           DO bra = 1, bra_confs2
              q = lookup_2b_configs(3,ch1)%ival2(1,bra)
              r = lookup_2b_configs(3,ch1)%ival2(2,bra)
              IF ( cut3b(p,q,r) ) cycle
              IF ( p < q ) then
                 nconfs = nconfs + 1
                 lookup_t3_configs(1,ch3)%ival2(1,nconfs) = p
                 lookup_t3_configs(1,ch3)%ival2(2,nconfs) = q
                 lookup_t3_configs(1,ch3)%ival2(3,nconfs) = r
              end IF
           end DO
           exit
        end DO
     end DO

     ! hhh
     nconfs = 0
     DO p = 1, below_ef
        nxp = all_orbit%nx(p)
        nyp = all_orbit%ny(p)
        nzp = all_orbit%nz(p)
        tzp = all_orbit%tz(p)
        DO ch2 = 1, channels_2b%number_confs
           IF ( number_2b(1,ch2) == 0 ) cycle
           Nx1 = channels_2b%config_NxNyNz_Tz(4*ch2-3)
           Ny1 = channels_2b%config_NxNyNz_Tz(4*ch2-2)
           Nz1 = channels_2b%config_NxNyNz_Tz(4*ch2-1)
           Tz1 = channels_2b%config_NxNyNz_Tz(4*ch2)
           IF ( Nx1 + nxp /= nx3 ) cycle
           IF ( Ny1 + nyp /= ny3 ) cycle
           IF ( Nz1 + nzp /= nz3 ) cycle
           IF ( 2*Tz1 + tzp /= tz3 ) cycle
           ket_confs2 = number_2b(1,ch2)
           DO ket = 1, ket_confs2
              q = lookup_2b_configs(1,ch2)%ival2(1,ket)
              r = lookup_2b_configs(1,ch2)%ival2(2,ket)
              IF ( cut3b(p,q,r) ) cycle
              IF ( p < q ) then
                 nconfs = nconfs + 1
                 lookup_t3_configs(2,ch3)%ival2(1,nconfs) = p
                 lookup_t3_configs(2,ch3)%ival2(2,nconfs) = q
                 lookup_t3_configs(2,ch3)%ival2(3,nconfs) = r
              end IF
           end DO
           exit
        end DO
     end DO
  end DO
  
  CALL mem_report('T3 structures')
  
end SUBROUTINE setup_t3_mbpt



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE t2_cross_recouple
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE ang_mom_functions
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, ch_t2, bra_confs,ket_confs, bra_min,bra_max
  INTEGER :: bra,ket, t2bra,t2ket, phase
  INTEGER :: a,b,i,j
  
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(3,ch) == 0 ) cycle
     IF ( number_2bcross(2,ch) == 0 ) cycle
     bra_confs = number_2bcross(3,ch)
     ket_confs = number_2bcross(2,ch)

     t2_ccm_cross(ch)%cval = 0.d0
     
     IF ( check_my_channel_v2b_phhp_cross(ch) == 0 ) cycle
     bra_min = mapping_v2b_phhp_cross(iam+1,ch,2)
     bra_max = mapping_v2b_phhp_cross(iam+1,ch,3)
     
     DO bra = bra_min, bra_max
        a   = lookup_2bcross_configs(3,ch)%ival2(1,bra)
        j   = lookup_2bcross_configs(3,ch)%ival2(2,bra)
        DO ket = 1, ket_confs
           i   = lookup_2bcross_configs(2,ch)%ival2(1,ket)
           b   = lookup_2bcross_configs(2,ch)%ival2(2,ket)

           ch_t2 = pp_channel_2b%ival2(a,b)
           IF ( ch_t2 == 0 ) cycle
           IF ( ch_t2 /= hh_channel_2b%ival2(i,j) ) cycle
           
           phase = 1
           t2bra = pp_config_2b%ival2(a,b)
           t2ket = hh_config_2b%ival2(i,j)
           IF ( b < a ) phase = -phase
           IF ( j < i ) phase = -phase
           
           t2_ccm_cross(ch)%cval(bra,ket) = phase * t2_ccm(ch_t2)%cval(t2bra,t2ket)
        end DO
     end DO
  end DO

  CALL mpi_barrier(mpi_comm_world,ierror)
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(3,ch) == 0 ) cycle
     IF ( number_2bcross(2,ch) == 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,t2_ccm_cross(ch)%cval,size(t2_ccm_cross(ch)%cval),&
          mpi_complex16,mpi_sum,mpi_comm_world,ierror)
     CALL check_mpi(ierror, 't2_cross_recouple', 'allreduce')
  end DO
  
END SUBROUTINE t2_cross_recouple


SUBROUTINE t2_add_cross
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE ang_mom_functions
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, ch_t2
  INTEGER :: bra_confs,ket_confs, bra_min,bra_max
  INTEGER :: bra,ket, t2bra,t2ket
  INTEGER :: a,b,i,j
  COMPLEX(dpc) :: t2
  
  CALL mpi_barrier(mpi_comm_world,ierror)
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(3,ch) == 0 ) cycle
     IF ( number_2bcross(2,ch) == 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,t2_ccm_eqn_cross(ch)%cval,size(t2_ccm_eqn_cross(ch)%cval),&
          mpi_complex16,mpi_sum,mpi_comm_world,ierror)
     CALL check_mpi(ierror, 't2_add_cross', 'allreduce')
  end DO

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)
     
     IF ( check_my_channel_v2b_pphh(ch) == 0 ) cycle
     bra_min = mapping_v2b_pphh(iam+1,ch,2)
     bra_max = mapping_v2b_pphh(iam+1,ch,3)

     DO bra = bra_min, bra_max
        a   = lookup_2b_configs(3,ch)%ival2(1,bra)
        b   = lookup_2b_configs(3,ch)%ival2(2,bra)
        ! (ab)(ij)
        DO ket = 1, ket_confs
           i   = lookup_2b_configs(1,ch)%ival2(1,ket)
           j   = lookup_2b_configs(1,ch)%ival2(2,ket)
           
           ch_t2 = ph_channel_2bcross%ival2(a,j)
           IF ( ch_t2 == 0 ) cycle
           IF ( ch_t2 /= hp_channel_2bcross%ival2(i,b) ) cycle
           t2bra = ph_config_2bcross%ival2(a,j)
           t2ket = hp_config_2bcross%ival2(i,b)
           t2    = t2_ccm_eqn_cross(ch_t2)%cval(t2bra,t2ket)
           t2_ccm_eqn(ch)%cval(bra,ket) = t2_ccm_eqn(ch)%cval(bra,ket) + t2
        end DO
        ! P(ab)(ij)
        DO ket = 1, ket_confs
           i   = lookup_2b_configs(1,ch)%ival2(1,ket)
           j   = lookup_2b_configs(1,ch)%ival2(2,ket)
           
           ch_t2 = ph_channel_2bcross%ival2(b,j)
           IF ( ch_t2 == 0 ) cycle
           IF ( ch_t2 /= hp_channel_2bcross%ival2(i,a) ) cycle
           t2bra = ph_config_2bcross%ival2(b,j)
           t2ket = hp_config_2bcross%ival2(i,a)
           t2    = t2_ccm_eqn_cross(ch_t2)%cval(t2bra,t2ket)
           t2_ccm_eqn(ch)%cval(bra,ket) = t2_ccm_eqn(ch)%cval(bra,ket) - t2
        end DO
        ! (ab)P(ij)
        DO ket = 1, ket_confs
           i   = lookup_2b_configs(1,ch)%ival2(1,ket)
           j   = lookup_2b_configs(1,ch)%ival2(2,ket)
           
           ch_t2 = ph_channel_2bcross%ival2(a,i)
           IF ( ch_t2 == 0 ) cycle
           IF ( ch_t2 /= hp_channel_2bcross%ival2(j,b) ) cycle
           t2bra = ph_config_2bcross%ival2(a,i)
           t2ket = hp_config_2bcross%ival2(j,b)
           t2    = t2_ccm_eqn_cross(ch_t2)%cval(t2bra,t2ket)
           t2_ccm_eqn(ch)%cval(bra,ket) = t2_ccm_eqn(ch)%cval(bra,ket) - t2
        end DO
        ! P(ab)P(ij)
        DO ket = 1, ket_confs
           i   = lookup_2b_configs(1,ch)%ival2(1,ket)
           j   = lookup_2b_configs(1,ch)%ival2(2,ket)
           
           ch_t2 = ph_channel_2bcross%ival2(b,i)
           IF ( ch_t2 == 0 ) cycle
           IF ( ch_t2 /= hp_channel_2bcross%ival2(j,a) ) cycle
           t2bra = ph_config_2bcross%ival2(b,i)
           t2ket = hp_config_2bcross%ival2(j,a)
           t2    = t2_ccm_eqn_cross(ch_t2)%cval(t2bra,t2ket)
           t2_ccm_eqn(ch)%cval(bra,ket) = t2_ccm_eqn(ch)%cval(bra,ket) + t2
        end DO
     end DO
  end DO           

END SUBROUTINE t2_add_cross


SUBROUTINE print_gs_hdf5
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE hdf5_wrapper
  USE hdf5
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, error, bra_confs,ket_confs
  INTEGER(HID_T) :: file_id
  REAL(dp), ALLOCATABLE :: t2_amp(:,:)
  REAL(dp), DIMENSION(4) :: energies
  INTEGER, DIMENSION(7) :: params
  CHARACTER(len=13) :: dsetname0 = "CC_parameters"
  CHARACTER(len=16) :: dsetname1 = "CC_energies_real"
  CHARACTER(len=16) :: dsetname2 = "CC_energies_imag"
  CHARACTER(len=19) :: dsetname30 = "T2_Amplitudes_real_"
  CHARACTER(len=19) :: dsetname40 = "T2_Amplitudes_imag_"
  CHARACTER(len=100) :: dsetname3, dsetname4
  
  IF ( iam == 0 ) then
     CALL h5fcreate_f(cc_file,H5F_ACC_TRUNC_F,file_id,error)
     
     params = 0
     params(1) = Nmax
     params(2) = Nocc
     params(3) = Pocc
     params(4) = int(1000.d0*rho)
     params(5) = 0
     IF ( cc_approx == 2 ) params(5) = 2
     params(6) = tnf_approx
     params(7) = 1000
     IF ( cc_approx == 2 ) params(7) = t3_cut
     error = write_integer_vector(file_id, params, trim(dsetname0))
     
     energies = 0.d0
     energies(1) = real(eccdt)
     energies(2) = real(e0)
     energies(3) = real(ecorr2)
     IF ( cc_approx > 1 .and. tnf_approx > 1 ) then
        energies(4) = real(ecorr3)
     end if
     error = write_double_vector(file_id, energies, trim(dsetname1))
     
     energies = 0.d0
     energies(1) = aimag(eccdt)
     energies(2) = aimag(e0)
     energies(3) = aimag(ecorr2)
     IF ( cc_approx > 1 .and. tnf_approx > 1 ) then
        energies(4) = aimag(ecorr3)
     end if
     error = write_double_vector(file_id, energies, trim(dsetname2))

     DO ch = 1, channels_2b%number_confs
        IF ( number_2b(3,ch) == 0 ) cycle
        IF ( number_2b(1,ch) == 0 ) cycle
        bra_confs = number_2b(3,ch)
        ket_confs = number_2b(1,ch)
        write(dsetname3,'(a,i6.6)') dsetname30, ch
        write(dsetname4,'(a,i6.6)') dsetname40, ch
        
        ALLOCATE( t2_amp(bra_confs, ket_confs) )
        t2_amp = real(t2_ccm(ch)%cval)
        error = write_double_matrix(file_id, t2_amp, trim(dsetname3))
        t2_amp = aimag(t2_ccm(ch)%cval)
        error = write_double_matrix(file_id, t2_amp, trim(dsetname4))
        
        DEALLOCATE( t2_amp )
     end DO
     
     CALL h5fclose_f(file_id,error)
  end IF
  
END SUBROUTINE print_gs_hdf5

SUBROUTINE read_gs_hdf5
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE hdf5_wrapper
  USE hdf5
  USE mpi_check

  IMPLICIT NONE
  INTEGER :: cc_approx0, t3_cut0
  INTEGER :: ch, error, bra,ket, bra_confs,ket_confs
  INTEGER(HID_T) :: file_id
  REAL(dp), DIMENSION(4) :: energies_r, energies_i
  INTEGER, DIMENSION(7) :: params
  CHARACTER(len=13) :: dsetname0 = "CC_parameters"
  CHARACTER(len=16) :: dsetname1 = "CC_energies_real"
  CHARACTER(len=16) :: dsetname2 = "CC_energies_imag"
  CHARACTER(len=19) :: dsetname30 = "T2_Amplitudes_real_"
  CHARACTER(len=19) :: dsetname40 = "T2_Amplitudes_imag_"
  CHARACTER(len=100) :: dsetname3, dsetname4
  INTEGER, DIMENSION(:), POINTER :: dummyi1
  REAL(dp), DIMENSION(:), POINTER :: dummyr1_r, dummyr1_i
  REAL(dp), DIMENSION(:,:), POINTER :: dummyr2_r, dummyr2_i

  cc_approx0 = cc_approx
  t3_cut0 = t3_cut
  IF ( cc_approx < 2 ) then
     cc_approx0 = 0
     t3_cut0 = 1000
  end IF

  CALL mpi_barrier(mpi_comm_world,ierror)
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     t2_ccm(ch)%cval = 0.d0
  end DO
  
  IF ( iam == 0 ) then
     CALL h5fopen_f(cc_file,H5F_ACC_RDONLY_F,file_id,error)
     
     dummyi1 => NULL()
     error = read_integer_vector(file_id, dummyi1, trim(dsetname0))
     params = dummyi1(:)
     DEALLOCATE(dummyi1)
     
     write(6,*) params(1), params(2), params(3), params(4), params(5), params(6), params(7)
     IF ( params(1) /= Nmax .or. params(2) /= Nocc .or. params(3) /= Pocc .or. params(4) /= int(1000.d0*rho) &
          .or. params(5) /= cc_approx0 .or. params(6) /= tnf_approx .or. params(7) /= t3_cut0 ) then
        IF ( iam == 0 ) write(6,*)
        IF ( iam == 0 ) write(6,*) ' Mismatch between CC parameters from file and input'
        IF ( iam == 0 ) write(6,*)
        IF ( iam == 0 ) write(6,*) ' CC parameters from file: '
        IF ( iam == 0 ) write(6,'(A6,I4,3x,A6,I4,3x,A6,I4,3x,A5,f5.2,3x,A11,I4,3x,A12,I4,3x,A8,I5)') &
             'Nmax =',params(1),'Nocc =',params(2),'Pocc =',params(3),'rho =',real(params(4))/1000.d0,&
             'cc_approx =',params(5),'tnf_approx =',params(6),'t3_cut =',params(7)
        IF ( iam == 0 ) write(6,*) ' CC parameters from input: '
        IF ( iam == 0 ) write(6,'(A6,I4,3x,A6,I4,3x,A6,I4,3x,A5,f5.2,3x,A11,I4,3x,A12,I4,3x,A8,I5)') &
             'Nmax =',Nmax,'Nocc =',Nocc,'Pocc =',Pocc,'rho =',rho,&
             'cc_approx =',cc_approx,'tnf_approx =',tnf_approx,'t3_cut =',t3_cut
        IF ( iam == 0 ) write(6,*)
        call h5fclose_f(file_id, error)
        stop
     end IF
     
     dummyr1_r => NULL()
     error = read_double_vector(file_id, dummyr1_r, trim(dsetname1))
     energies_r = dummyr1_r(:)
     dummyr1_i => NULL()
     error = read_double_vector(file_id, dummyr1_i, trim(dsetname2))
     energies_i = dummyr1_i(:)
     DEALLOCATE(dummyr1_r, dummyr1_i)
     
     write(6,*)
     write(6,*) ' ECC = ', energies_r(1), energies_i(1)
     write(6,*) ' E0  = ', energies_r(2), energies_i(2)
     write(6,*) ' E2  = ', energies_r(3), energies_i(3)
     write(6,*) ' E3  = ', energies_r(4), energies_i(4)
     write(6,*)
     e0 = dcmplx(energies_r(2), energies_i(2))
     ecorr2 = dcmplx(energies_r(3), energies_i(3))
     ecorr3 = dcmplx(energies_r(4), energies_i(4))
     ecorr = ecorr2 + ecorr3
     eccdt = e0 + ecorr
     
     DO ch = 1, channels_2b%number_confs
        IF ( number_2b(3,ch) == 0 ) cycle
        IF ( number_2b(1,ch) == 0 ) cycle
        bra_confs = number_2b(3,ch)
        ket_confs = number_2b(1,ch)
        write(dsetname3,'(a,i6.6)') dsetname30, ch
        write(dsetname4,'(a,i6.6)') dsetname40, ch
        
        dummyr2_r => NULL()
        error = read_double_matrix(file_id, dummyr2_r, trim(dsetname3))
        dummyr2_i => NULL()
        error = read_double_matrix(file_id, dummyr2_i, trim(dsetname4))
        
        DO bra = 1, bra_confs
           DO ket = 1, ket_confs
              t2_ccm(ch)%cval(bra,ket) = dcmplx(dummyr2_r(bra,ket), dummyr2_i(bra,ket))
           end DO
        end DO
        DEALLOCATE(dummyr2_r, dummyr2_i)
     end DO
     
     CALL h5fclose_f(file_id,error)
  end IF

  CALL mpi_barrier(mpi_comm_world,ierror)
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     call mpi_bcast(t2_ccm(ch)%cval, size(t2_ccm(ch)%cval), mpi_complex16, master, mpi_comm_world, ierror)
     CALL check_mpi(ierror, 'read_gs_hdf5', 'bcast')
  end DO
  CALL t2_cross_recouple
  
END SUBROUTINE read_gs_hdf5
