
SUBROUTINE setup_preom_amplitudes(ndim)
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE ang_mom_functions
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER(i8), intent(out) :: ndim
  INTEGER :: ch,ch2, ket_confs
  INTEGER :: Nx,Ny,Nz,Tz
  INTEGER :: a, nxa,nya,nza,tza
  INTEGER :: i, nxi,nyi,nzi,tzi

  IF ( iam == 0 ) write(6,*) 'Setting up PR-EOM structures...'

  ! 1p + 1p2h dimension
  ndim = 0

  ! find PR-EOM indices
  ALLOCATE( r1_preom_ind(2) ) ! spin-up/spin-down
  r1_preom_ind = 0
  CALL mem_register('preom2', REAL(2 * 4.d0, dp))
  DO i = 1, below_ef
     nxi = all_orbit%nx(i)
     nyi = all_orbit%ny(i)
     nzi = all_orbit%nz(i)
     tzi = all_orbit%tz(i)
     IF ( nxi /= nx_preom ) cycle
     IF ( nyi /= ny_preom ) cycle
     IF ( nzi /= nz_preom ) cycle
     IF ( tzi /= tz_preom ) cycle
     ndim = ndim + 1
     r1_preom_ind(ndim) = i
  end DO
  
  ! setup 1h-amplitudes
  ALLOCATE( r1_preom(2) )
  ALLOCATE( r1_preom_eqn(2) )
  r1_preom = 0.d0
  r1_preom_eqn = 0.d0
  CALL mem_register('preom2', REAL(4 * 16.d0, dp))

  ! find PR-EOM indices
  ALLOCATE( r2_preom_ind(2*channels_2b%number_confs) )
  r2_preom_ind = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     Nx = channels_2b%config_NxNyNz_Tz(ch*4-3)
     Ny = channels_2b%config_NxNyNz_Tz(ch*4-2)
     Nz = channels_2b%config_NxNyNz_Tz(ch*4-1)
     Tz = channels_2b%config_NxNyNz_Tz(ch*4)
     DO a = below_ef+1, tot_orbs
        nxa = all_orbit%nx(a)
        nya = all_orbit%ny(a)
        nza = all_orbit%nz(a)
        tza = all_orbit%tz(a)
        IF ( Nx - nxa /= nx_preom ) cycle
        IF ( Ny - nya /= ny_preom ) cycle
        IF ( Nz - nza /= nz_preom ) cycle
        IF ( 2*Tz - tza /= tz_preom ) cycle
        ndim = ndim + number_2b(1,ch)
        IF ( all_orbit%sz(a) == -1 ) r2_preom_ind(2*ch-1) = a
        IF ( all_orbit%sz(a) == 1  ) r2_preom_ind(2*ch)   = a
     end DO
  end DO
  
  ! setup 2p1h-amplitudes
  ALLOCATE( r2_preom(channels_2b%number_confs) )
  ALLOCATE( r2_preom_eqn(channels_2b%number_confs) )
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     ket_confs = number_2b(1,ch)
     
     ALLOCATE( r2_preom(ch)%cval(2,ket_confs) )
     ALLOCATE( r2_preom_eqn(ch)%cval(2,ket_confs) )
     r2_preom(ch)%cval = 0.d0
     r2_preom_eqn(ch)%cval = 0.d0
     CALL mem_register('preom2', REAL(2 * 2*ket_confs * 16.d0, dp))
  end DO

  ! find PR-EOM indices
  ALLOCATE( r2_preom_ind_cross(2*channels_2bcross%number_confs) )
  r2_preom_ind_cross = 0
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     Nx = channels_2bcross%config_NxNyNz_Tz(ch2*4-3)
     Ny = channels_2bcross%config_NxNyNz_Tz(ch2*4-2)
     Nz = channels_2bcross%config_NxNyNz_Tz(ch2*4-1)
     Tz = channels_2bcross%config_NxNyNz_Tz(ch2*4)
     DO i = 1, below_ef
        nxi = all_orbit%nx(i)
        nyi = all_orbit%ny(i)
        nzi = all_orbit%nz(i)
        tzi = all_orbit%tz(i)
        IF ( nxi + Nx /= nx_preom ) cycle
        IF ( nyi + Ny /= ny_preom ) cycle
        IF ( nzi + Nz /= nz_preom ) cycle
        IF ( tzi + 2*Tz /= tz_preom ) cycle
        IF ( all_orbit%sz(i) == -1 ) r2_preom_ind_cross(2*ch2-1) = i
        IF ( all_orbit%sz(i) == 1  ) r2_preom_ind_cross(2*ch2)   = i
     end DO
  end DO
  
  ! setup 1p2h-cross-amplitudes
  ALLOCATE( r2_preom_cross(channels_2bcross%number_confs) )
  ALLOCATE( r2_preom_eqn_cross(channels_2bcross%number_confs) )
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_preom_ind_cross(2*ch2) == 0 ) cycle
     ket_confs = number_2bcross(2,ch2)
     
     ALLOCATE( r2_preom_cross(ch2)%cval(2,ket_confs) )
     ALLOCATE( r2_preom_eqn_cross(ch2)%cval(2,ket_confs) )
     r2_preom_cross(ch2)%cval = 0.d0
     r2_preom_eqn_cross(ch2)%cval = 0.d0
     CALL mem_register('preom2', REAL(2 * 2*ket_confs * 16.d0, dp))
  end DO

  CALL setup_proc_mappings_preom
  IF ( iam == 0 ) write(6,'(A19,26x,I14)') 'PR-EOM2 dimension:', ndim
  ! [old memory print removed - replaced by mem_report]
  IF ( iam == 0 ) write(6,*)
  
end SUBROUTINE setup_preom_amplitudes


SUBROUTINE setup_preom3_amplitudes(ndim3)
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE ang_mom_functions
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER(i8), INTENT(OUT) :: ndim3
  INTEGER :: number_channels, ch3,ch1,ch2
  INTEGER :: i,j,k, k1,k2,k3,k4, iind
  INTEGER :: nx3,ny3,nz3,tz3, Nx1,Ny1,Nz1,Tz1
  INTEGER :: ket, bra_confs,ket_confs
  INTEGER :: ket_min,ket_max
  LOGICAL :: ch_added

  cut_3b = eom3_cut
  
  n3min = 1000
  n3max = -1000
  t3min = 1000
  t3max = -1000
  DO i = 1, below_ef-2
     DO j = i+1, below_ef-1
        DO k = j+1, below_ef
           IF ( all_orbit%nx(i) + all_orbit%nx(j) + all_orbit%nx(k) > n3max ) &
                n3max = all_orbit%nx(i) + all_orbit%nx(j) + all_orbit%nx(k)
           IF ( all_orbit%nx(i) + all_orbit%nx(j) + all_orbit%nx(k) < n3min ) &
                n3min = all_orbit%nx(i) + all_orbit%nx(j) + all_orbit%nx(k)
           IF ( all_orbit%tz(i) + all_orbit%tz(j) + all_orbit%tz(k) > t3max ) &
                t3max = all_orbit%tz(i) + all_orbit%tz(j) + all_orbit%tz(k)
           IF ( all_orbit%tz(i) + all_orbit%tz(j) + all_orbit%tz(k) < t3min ) &
                t3min = all_orbit%tz(i) + all_orbit%tz(j) + all_orbit%tz(k)
        end DO
     end DO
  end DO
  
  number_channels = 0
  DO nx3 = n3min, n3max
     DO ny3 = n3min, n3max
        DO nz3 = n3min, n3max
           DO tz3 = t3min, t3max, 2
              CALL number_3b_confs(nx3, ny3, nz3, tz3, 2, ket_confs) ! 2 for hhh
              IF ( ket_confs <= 0 ) cycle
              
              DO ch1 = 1, channels_2b%number_confs
                 IF ( number_2b(3,ch1) == 0 ) cycle
                 Nx1 = channels_2b%config_NxNyNz_Tz(ch1*4-3)
                 Ny1 = channels_2b%config_NxNyNz_Tz(ch1*4-2)
                 Nz1 = channels_2b%config_NxNyNz_Tz(ch1*4-1)
                 Tz1 = channels_2b%config_NxNyNz_Tz(ch1*4)
                 IF ( nx3 - Nx1 /= nx_preom ) cycle
                 IF ( ny3 - Ny1 /= ny_preom ) cycle
                 IF ( nz3 - Nz1 /= nz_preom ) cycle
                 IF ( tz3 - 2*Tz1 /= tz_preom ) cycle
                 number_channels = number_channels + 1
                 exit
              end DO
              
           end DO
        end DO
     end DO
  end DO
  channels_preom3%number_confs = number_channels

  ALLOCATE( channels_preom3%config_NxNyNz_Tz(4*number_channels) )
  CALL mem_register('preom3', REAL(4 * number_channels * 4.d0, dp))
  ALLOCATE( ch1_preom3(number_channels) )
  ch1_preom3 = 0
  CALL mem_register('preom3', REAL(number_channels * 4.d0, dp))
  
  number_channels = 0
  ndim3 = 0
  DO nx3 = n3min, n3max
     DO ny3 = n3min, n3max
        DO nz3 = n3min, n3max
           DO tz3 = t3min, t3max, 2
              CALL number_3b_confs(nx3, ny3, nz3, tz3, 2, ket_confs) ! 1 for ppp
              IF ( ket_confs <= 0 ) cycle
              
              ch_added = .FALSE.
              DO ch1 = 1, channels_2b%number_confs
                 IF ( number_2b(3,ch1) == 0 ) cycle
                 bra_confs = number_2b(3,ch1)
                 Nx1 = channels_2b%config_NxNyNz_Tz(ch1*4-3)
                 Ny1 = channels_2b%config_NxNyNz_Tz(ch1*4-2)
                 Nz1 = channels_2b%config_NxNyNz_Tz(ch1*4-1)
                 Tz1 = channels_2b%config_NxNyNz_Tz(ch1*4)
                 IF ( nx3 - Nx1 /= nx_preom ) cycle
                 IF ( ny3 - Ny1 /= ny_preom ) cycle
                 IF ( nz3 - Nz1 /= nz_preom ) cycle
                 IF ( tz3 - 2*Tz1 /= tz_preom ) cycle
                 ndim3 = ndim3 + bra_confs*ket_confs
                 number_channels = number_channels + 1
                 ch1_preom3(number_channels) = ch1
                 ch_added = .TRUE.
                 exit
              end DO
              IF ( .not. ch_added ) cycle
              ch3 = number_channels
              
              k1 = ch3*4 - 3
              k2 = ch3*4 - 2
              k3 = ch3*4 - 1
              k4 = ch3*4
              channels_preom3%config_NxNyNz_Tz(k1) = nx3
              channels_preom3%config_NxNyNz_Tz(k2) = ny3
              channels_preom3%config_NxNyNz_Tz(k3) = nz3
              channels_preom3%config_NxNyNz_Tz(k4) = tz3
              
           end DO
        end DO
     end DO
  end DO
  IF ( iam == 0 ) write(6,'(A19,26x,I14)') 'PR-EOM3 dimension:', ndim3  
  
  
  ! Setup PA_EOM3 mapping
  CALL setup_proc_mappings_preom3
  CALL mem_report('PR-EOM3 structures')


  ndim3 = 0
  ALLOCATE( r3_preom(ch3_preom_min:ch3_preom_max) )
  DO ch3 = ch3_preom_min, ch3_preom_max
     ALLOCATE( r3_preom(ch3)%val1(klimits_preom3(ch3,1):klimits_preom3(ch3,2)) )    
     DO iind = klimits_preom3(ch3,1), klimits_preom3(ch3,2)
        i         = klist_preom3(ch3)%ival2(iind,1)
        ch2       = klist_preom3(ch3)%ival2(iind,2)
        ket_min   = mapping_preom3(ch3)%ival2(iind,1)
        ket_max   = mapping_preom3(ch3)%ival2(iind,2)
        ch1       = ch1_preom3(ch3)
        bra_confs = number_2b(3,ch1)
        ALLOCATE( r3_preom(ch3)%val1(iind)%cval(bra_confs,ket_min:ket_max) )
        r3_preom(ch3)%val1(iind)%cval = 0.d0
        CALL mem_register('preom3', REAL(bra_confs*(ket_max-ket_min+1) * 16.d0, dp))
        
        DO ket = ket_min, ket_max
           j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
           k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
           IF ( i >= j ) cycle
           ndim3 = ndim3 + bra_confs
        end DO
        
     end DO
  end DO
  CALL mpi_allreduce(mpi_in_place,ndim3,1,mpi_integer8,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_preom3_amplitudes', 'allreduce')
  IF ( iam == 0 ) write(6,'(A26,19x,I14)') 'PR-EOM3 vector dimension:', ndim3
  CALL mem_report('PR-EOM3 amplitudes')
  
end SUBROUTINE setup_preom3_amplitudes


SUBROUTINE populate_preom(ndim1, ndim2, vecin)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER(i8), INTENT(in) :: ndim1, ndim2
  COMPLEX(dpc), INTENT(in) :: vecin(ndim1:ndim2)
  INTEGER :: ch, ch1,ch2,ch3
  INTEGER :: bra,ket, bra_confs,ket_confs, ket_min,ket_max
  INTEGER :: i,j,k, b,c, iind
  INTEGER(i8) :: ii
  
  IF ( iam == 0 ) WRITE(6,*) "...Distributing right EOM vector..."

  r1_preom = 0.d0
  ii = 1
  IF ( eom%my_start <= ii .and. ii <= eom%my_stop ) then
     r1_preom(ii) = vecin(ii)
  end IF
  ii = ii + 1
  IF ( eom%my_start <= ii .and. ii <= eom%my_stop ) then
     r1_preom(ii) = vecin(ii)
  end IF
  CALL mpi_allreduce(mpi_in_place,r1_preom,size(r1_preom),mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'populate_preom', 'allreduce')
  
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     ket_confs = number_2b(1,ch)
     r2_preom(ch)%cval = 0.d0
     DO ket = 1, ket_confs
        DO bra = 1, 2
           ii = ii + 1
           IF ( eom%my_start <= ii .and. ii <= eom%my_stop ) then
              r2_preom(ch)%cval(bra,ket) = vecin(ii)
           end IF
        end DO
     end DO
     CALL mpi_allreduce(mpi_in_place,r2_preom(ch)%cval,size(r2_preom(ch)%cval),mpi_complex16,mpi_sum,mpi_comm_world,ierror)
     CALL check_mpi(ierror, 'populate_preom', 'allreduce')
  end DO
  
  CALL r2_preom_cross_recouple
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  IF ( eom_approx > 0 ) then

     ii = max(eom_ndim2, eom%my_start-1)
     DO ch3 = ch3_preom_min, ch3_preom_max
        ch1       = ch1_preom3(ch3)
        bra_confs = number_2b(3,ch1)
        DO iind = klimits_preom3(ch3,1), klimits_preom3(ch3,2)
           i       = klist_preom3(ch3)%ival2(iind,1)
           ch2     = klist_preom3(ch3)%ival2(iind,2)
           ket_min = mapping_preom3(ch3)%ival2(iind,1)
           ket_max = mapping_preom3(ch3)%ival2(iind,2)
           
           r3_preom(ch3)%val1(iind)%cval = 0.d0
           
           DO bra = 1, bra_confs
              b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
              c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
              DO ket = ket_min, ket_max
                 j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
                 k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
                 IF ( i >= j ) cycle
                 ii = ii + 1
                 r3_preom(ch3)%val1(iind)%cval(bra,ket) = vecin(ii)
              end DO
           end DO
        end DO
     end DO
     CALL mpi_barrier(mpi_comm_world,ierror)
     CALL expand_preom_2p3h
     
  end IF
  
  IF ( test == 4 ) CALL build_preom_test

end SUBROUTINE populate_preom


SUBROUTINE populate_preom_vec(ndim1, ndim2, vecout)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER(i8), INTENT(in) :: ndim1, ndim2
  COMPLEX(dpc), INTENT(out) :: vecout(ndim1:ndim2)
  INTEGER :: ch, ch3, ch1,ch2, i,j,k,b,c, iind
  INTEGER :: ket_min,ket_max, bra,ket
  INTEGER :: bra_confs,ket_confs
  INTEGER(i8) :: ii

  IF ( iam == 0 ) WRITE(6,*) "...Collecting right EOM vector..."
  
  CALL r2_preom_add_cross
  CALL mpi_barrier(mpi_comm_world,ierror)
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,r2_preom_eqn(ch)%cval,size(r2_preom_eqn(ch)%cval),mpi_complex16,mpi_sum,mpi_comm_world,ierror)
     CALL check_mpi(ierror, 'populate_preom_vec', 'allreduce')
  end DO
  
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF ( test == 4 ) then
     CALL build_r1_preom_eqn_test
     CALL build_r2_preom_eqn_test
     IF ( eom_approx > 0 ) CALL build_r3_preom_eqn_test
  end IF
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  vecout = 0.d0

  ii = 1
  IF ( eom%my_start <= ii .and. ii <= eom%my_stop ) then
     vecout(ii) = r1_preom_eqn(ii)
  end IF
  ii = ii + 1
  IF ( eom%my_start <= ii .and. ii <= eom%my_stop ) then
     vecout(ii) = r1_preom_eqn(ii)
  end IF

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     ket_confs = number_2b(1,ch)
     DO ket = 1, ket_confs
        DO bra = 1, 2
           ii = ii + 1
           IF ( eom%my_start <= ii .and. ii <= eom%my_stop ) then
              vecout(ii) = r2_preom_eqn(ch)%cval(bra,ket)
           end IF
        end DO
     end DO
  end DO

  
  IF ( eom_approx > 0 ) then

     ii = max(eom_ndim2, eom%my_start-1)
     DO ch3 = ch3_preom_min, ch3_preom_max
        ch1       = ch1_preom3(ch3)
        bra_confs = number_2b(3,ch1)
        DO iind = klimits_preom3(ch3,1), klimits_preom3(ch3,2)
           i       = klist_preom3(ch3)%ival2(iind,1)
           ch2     = klist_preom3(ch3)%ival2(iind,2)
           ket_min = mapping_preom3(ch3)%ival2(iind,1)
           ket_max = mapping_preom3(ch3)%ival2(iind,2)           
           DO bra = 1, bra_confs
              b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
              c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
              DO ket = ket_min, ket_max
                 j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
                 k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
                 IF ( i >= j ) cycle
                 ii = ii + 1
                 vecout(ii) = r3_preom(ch3)%val1(iind)%cval(bra,ket)
              end DO
           end DO
        end DO
     end DO
     
  end IF

end SUBROUTINE populate_preom_vec


SUBROUTINE expand_preom_2p3h
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE parallel
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch3,ch1,ch2, ket,ket1
  INTEGER :: bra_confs,ket_confs
  INTEGER :: ket_min,ket_max, i,j,k
  INTEGER :: iind,jind,kind1, ch_j
  INTEGER, allocatable :: klist_inv(:)
  TYPE (superblock_storage) :: r3_temp
  INTEGER :: group0
  
  IF ( iam == 0 ) WRITE(6,*) "...Expanding R3..."

  ALLOCATE( klist_inv(1:below_ef) )
  
  DO ch3 = ch3_preom_min, ch3_preom_max        
     ch1       = ch1_preom3(ch3)
     bra_confs = number_2b(3,ch1)
     
     ! Build inverse k map
     klist_inv = 0
     DO iind = 1, klimit_preom3(ch3)
        i = klist_preom3(ch3)%ival2(iind,1)
        klist_inv(i) = iind
     end DO
     
     ! Allocate r3_temp for all k_channels
     ALLOCATE( r3_temp%val1(klimit_preom3(ch3)) )
     DO iind = 1, klimit_preom3(ch3)
        ch2       = klist_preom3(ch3)%ival2(iind,2)
        ket_confs = number_2b(1,ch2)
        ALLOCATE( r3_temp%val1(iind)%cval(bra_confs,ket_confs) )
        r3_temp%val1(iind)%cval = 0.d0
     end DO
     
     ! Fill r3_temp for local c_channels
     DO iind = klimits_preom3(ch3,1), klimits_preom3(ch3,2)
        ket_min = mapping_preom3(ch3)%ival2(iind,1)
        ket_max = mapping_preom3(ch3)%ival2(iind,2)
        r3_temp%val1(iind)%cval(:,ket_min:ket_max) = r3_preom(ch3)%val1(iind)%cval(:,ket_min:ket_max)
     end DO
     
     ! Allreduce r3_temp for all c_channels
     group0 = group_preom3(iam+1, ch3)
     IF ( group0 > 0 ) then
        DO iind = 1, klimit_preom3(ch3)
           CALL mpi_allreduce(mpi_in_place,r3_temp%val1(iind)%cval,size(r3_temp%val1(iind)%cval),&
                mpi_complex16,mpi_sum,subcomm_preom3(group0),ierror)
           CALL check_mpi(ierror, 'expand_preom_2p3h', 'allreduce')
        end DO
     end IF

     DO iind = klimits_preom3(ch3,1), klimits_preom3(ch3,2)
        i         = klist_preom3(ch3)%ival2(iind,1)
        ch2       = klist_preom3(ch3)%ival2(iind,2)
        ket_min   = mapping_preom3(ch3)%ival2(iind,1)
        ket_max   = mapping_preom3(ch3)%ival2(iind,2)
        ket_confs = number_2b(1,ch2)
        DO ket = ket_min, ket_max
           j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
           k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
           IF ( i <= j ) cycle ! i(jk)

           jind = klist_inv(j)
           ch_j = klist_preom3(ch3)%ival2(jind,2)
           IF ( ch_j <= 0 ) cycle
           IF ( i == k ) cycle
           IF ( i > k ) then ! j(ki)
              ket1 = hh_config_2b%ival2(k,i)
              r3_preom(ch3)%val1(iind)%cval(:,ket) = r3_preom(ch3)%val1(iind)%cval(:,ket) &
                   + r3_temp%val1(jind)%cval(:,ket1)
           ELSE ! j(ik)
              ket1 = hh_config_2b%ival2(i,k)
              r3_preom(ch3)%val1(iind)%cval(:,ket) = r3_preom(ch3)%val1(iind)%cval(:,ket) &
                   - r3_temp%val1(jind)%cval(:,ket1)
           end IF

        end DO
     end DO
     
     DO kind1 = 1, klimit_preom3(ch3)
        DEALLOCATE( r3_temp%val1(kind1)%cval )
     end DO
     DEALLOCATE( r3_temp%val1 )     
  end DO
  DEALLOCATE( klist_inv )
  
end SUBROUTINE expand_preom_2p3h


SUBROUTINE r2_preom_cross_recouple
  use parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch,ch2, ket, bra2,ket2
  INTEGER :: ket_min,ket_max
  INTEGER :: b,i,j
  REAL(dp) :: phase

  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_preom_ind_cross(2*ch2) == 0 ) cycle
     
     r2_preom_cross(ch2)%cval = 0.d0
     
     IF ( check_my_channel_preom_hhp_cross(ch2) == 0 ) cycle
     ket_min = mapping_preom_hhp_cross(iam+1,ch2,2)
     ket_max = mapping_preom_hhp_cross(iam+1,ch2,3)

     DO bra2 = 1, 2
        i = r2_preom_ind_cross(2*(ch2-1)+bra2)
        DO ket2 = ket_min, ket_max
           j   = lookup_2bcross_configs(2,ch2)%ival2(1,ket2)
           b   = lookup_2bcross_configs(2,ch2)%ival2(2,ket2)

           ch = hh_channel_2b%ival2(i,j)
           IF ( ch == 0 ) cycle
           phase = 1
           ket = hh_config_2b%ival2(i,j)
           IF ( j < i ) phase = -phase
           IF ( all_orbit%sz(b) == -1 ) then
              r2_preom_cross(ch2)%cval(bra2,ket2) = phase*r2_preom(ch)%cval(1,ket)
           ELSE
              r2_preom_cross(ch2)%cval(bra2,ket2) = phase*r2_preom(ch)%cval(2,ket)
           end IF
        end DO
     end DO
  end DO
  
  CALL mpi_barrier(mpi_comm_world,ierror)
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_preom_ind_cross(2*ch2) == 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,r2_preom_cross(ch2)%cval,size(r2_preom_cross(ch2)%cval),&
          mpi_complex16,mpi_sum,mpi_comm_world,ierror)
     CALL check_mpi(ierror, 'r2_preom_cross_recouple', 'allreduce')
  end DO
  
END SUBROUTINE r2_preom_cross_recouple


SUBROUTINE r2_preom_add_cross
  use parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE ang_mom_functions
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch,ch2, bra,ket, ket2
  INTEGER :: ket_min,ket_max
  INTEGER :: b,i,j
  
  CALL mpi_barrier(mpi_comm_world,ierror)
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_preom_ind_cross(2*ch2) == 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,r2_preom_eqn_cross(ch2)%cval,size(r2_preom_eqn_cross(ch2)%cval),&
          mpi_complex16,mpi_sum,mpi_comm_world,ierror)
     CALL check_mpi(ierror, 'r2_preom_add_cross', 'allreduce')
  end DO

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     
     IF ( check_my_channel_preom_phh(ch) == 0 ) cycle
     ket_min = mapping_preom_phh(iam+1,ch,2)
     ket_max = mapping_preom_phh(iam+1,ch,3)

     DO ket = ket_min, ket_max
        i   = lookup_2b_configs(1,ch)%ival2(1,ket)
        j   = lookup_2b_configs(1,ch)%ival2(2,ket)
        DO bra = 1, 2
           b = r2_preom_ind(2*(ch-1)+bra)
           ch2 = hp_channel_2bcross%ival2(j,b)
           IF ( ch2 == 0 ) cycle
           ket2 = hp_config_2bcross%ival2(j,b)
           IF ( all_orbit%sz(i) == -1 ) then
              r2_preom_eqn(ch)%cval(bra,ket) = r2_preom_eqn(ch)%cval(bra,ket) + r2_preom_eqn_cross(ch2)%cval(1,ket2)
           ELSE
              r2_preom_eqn(ch)%cval(bra,ket) = r2_preom_eqn(ch)%cval(bra,ket) + r2_preom_eqn_cross(ch2)%cval(2,ket2)
           end IF
        end DO
        ! P(ij)
        DO bra = 1, 2
           b = r2_preom_ind(2*(ch-1)+bra)
           ch2 = hp_channel_2bcross%ival2(i,b)
           IF ( ch2 == 0 ) cycle
           ket2 = hp_config_2bcross%ival2(i,b)
           IF ( all_orbit%sz(j) == -1 ) then
              r2_preom_eqn(ch)%cval(bra,ket) = r2_preom_eqn(ch)%cval(bra,ket) - r2_preom_eqn_cross(ch2)%cval(1,ket2)
           ELSE
              r2_preom_eqn(ch)%cval(bra,ket) = r2_preom_eqn(ch)%cval(bra,ket) - r2_preom_eqn_cross(ch2)%cval(2,ket2)
           end IF
        end DO
     end DO
  end DO
  
END SUBROUTINE r2_preom_add_cross


SUBROUTINE deallocate_preom
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch1, ch2, ch3, kind1

  CALL deallocate_hbar1b_I2
  CALL deallocate_hbar1b_I3
  CALL deallocate_hbar2b_I4pr
  CALL deallocate_hbar2b_I5pr_cross
  CALL deallocate_hbar2b_I7pr
  CALL deallocate_hbar3b_preom
  
  DEALLOCATE( preom_eigs )
  
  DEALLOCATE( r1_preom_ind )
  DEALLOCATE( r1_preom )
  DEALLOCATE( r1_preom_eqn )
  
  DO ch1 = 1, channels_2b%number_confs
     IF ( number_2b(3,ch1) == 0 ) cycle
     IF ( r2_preom_ind(2*ch1) == 0 ) cycle
     DEALLOCATE( r2_preom(ch1)%cval )
     DEALLOCATE( r2_preom_eqn(ch1)%cval )
  end DO
  DEALLOCATE( r2_preom )
  DEALLOCATE( r2_preom_eqn )
  DEALLOCATE( r2_preom_ind ) !!

  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_preom_ind_cross(2*ch2) == 0 ) cycle
     DEALLOCATE( r2_preom_cross(ch2)%cval )
     DEALLOCATE( r2_preom_eqn_cross(ch2)%cval )
  end DO
  DEALLOCATE( r2_preom_cross )
  DEALLOCATE( r2_preom_eqn_cross )
  DEALLOCATE( r2_preom_ind_cross ) !!

  DEALLOCATE( mapping_preom_phh )
  DEALLOCATE( check_my_channel_preom_phh )
  DEALLOCATE( mapping_preom_hhp_cross )
  DEALLOCATE( check_my_channel_preom_hhp_cross )

  DEALLOCATE( eom%all_starts )
  DEALLOCATE( eom%all_stops )

  IF ( eom_approx > 0 ) then
     DEALLOCATE( channels_preom3%config_NxNyNz_Tz )
     DEALLOCATE( ch1_preom3 )
     DEALLOCATE( group_preom3 )
     DEALLOCATE( subcomm_preom3 )
  
     DO ch3 = ch3_preom_min, ch3_preom_max
        DO kind1 = klimits_preom3(ch3,1), klimits_preom3(ch3,2)
           DEALLOCATE( r3_preom(ch3)%val1(kind1)%cval )
        end DO
        DEALLOCATE( r3_preom(ch3)%val1 )
     end DO
     DEALLOCATE( r3_preom )

     DO ch3 = ch3_preom_min, ch3_preom_max
        DEALLOCATE( klist_preom3(ch3)%ival2 )     
     end DO
     DEALLOCATE( klist_preom3 )
     DEALLOCATE( klimit_preom3 )

     DO ch3 = ch3_preom_min, ch3_preom_max
        IF ( ALLOCATED(mapping_preom3(ch3)%ival2) ) then
           DEALLOCATE( mapping_preom3(ch3)%ival2 )
        end IF
     end DO
     DEALLOCATE( mapping_preom3 )
     DEALLOCATE( klimits_preom3 )
  end IF
  
end SUBROUTINE deallocate_preom
