
SUBROUTINE setup_paeom_amplitudes(ndim)
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE ang_mom_functions
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER(i8), intent(out) :: ndim
  INTEGER :: ch,ch2, bra_confs,ket_confs
  INTEGER :: Nx,Ny,Nz,Tz
  INTEGER :: a, nxa,nya,nza,tza
  INTEGER :: i, nxi,nyi,nzi,tzi

  IF ( iam == 0 ) write(6,*) 'Setting up PA-EOM structures...'

  ! 1p + 2p1h dimension
  ndim = 0

  ! find PA-EOM indices
  ALLOCATE( r1_paeom_ind(2) ) ! spin-up/spin-down
  r1_paeom_ind = 0
  CALL mem_register('paeom2', REAL(2 * 4.d0, dp))
  DO a = below_ef+1, tot_orbs
     nxa = all_orbit%nx(a)
     nya = all_orbit%ny(a)
     nza = all_orbit%nz(a)
     tza = all_orbit%tz(a)
     IF ( nxa /= nx_paeom ) cycle
     IF ( nya /= ny_paeom ) cycle
     IF ( nza /= nz_paeom ) cycle
     IF ( tza /= tz_paeom ) cycle
     ndim = ndim + 1
     r1_paeom_ind(ndim) = a
  end DO
  
  ! setup 1p-amplitudes
  ALLOCATE( r1_paeom(2) )
  ALLOCATE( r1_paeom_eqn(2) )
  r1_paeom = 0.d0
  r1_paeom_eqn = 0.d0
  CALL mem_register('paeom2', REAL(4 * 16.d0, dp))

  ! find PA-EOM indices
  ALLOCATE( r2_paeom_ind(2*channels_2b%number_confs) )
  r2_paeom_ind = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     Nx = channels_2b%config_NxNyNz_Tz(ch*4-3)
     Ny = channels_2b%config_NxNyNz_Tz(ch*4-2)
     Nz = channels_2b%config_NxNyNz_Tz(ch*4-1)
     Tz = channels_2b%config_NxNyNz_Tz(ch*4)
     DO i = 1, below_ef
        nxi = all_orbit%nx(i)
        nyi = all_orbit%ny(i)
        nzi = all_orbit%nz(i)
        tzi = all_orbit%tz(i)
        IF ( Nx - nxi /= nx_paeom ) cycle
        IF ( Ny - nyi /= ny_paeom ) cycle
        IF ( Nz - nzi /= nz_paeom ) cycle
        IF ( 2*Tz - tzi /= tz_paeom ) cycle
        ndim = ndim + number_2b(3,ch)
        IF ( all_orbit%sz(i) == -1 ) r2_paeom_ind(2*ch-1) = i
        IF ( all_orbit%sz(i) == 1  ) r2_paeom_ind(2*ch)   = i
     end DO
  end DO
  
  ! setup 2p1h-amplitudes
  ALLOCATE( r2_paeom(channels_2b%number_confs) )
  ALLOCATE( r2_paeom_eqn(channels_2b%number_confs) )
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     
     ALLOCATE( r2_paeom(ch)%cval(bra_confs,2) )
     ALLOCATE( r2_paeom_eqn(ch)%cval(bra_confs,2) )
     r2_paeom(ch)%cval = 0.d0
     r2_paeom_eqn(ch)%cval = 0.d0
     CALL mem_register('paeom2', REAL(2 * bra_confs*2 * 16.d0, dp))
  end DO

  ! find PA-EOM indices
  ALLOCATE( r2_paeom_ind_cross(2*channels_2bcross%number_confs) )
  r2_paeom_ind_cross = 0
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     Nx = channels_2bcross%config_NxNyNz_Tz(ch2*4-3)
     Ny = channels_2bcross%config_NxNyNz_Tz(ch2*4-2)
     Nz = channels_2bcross%config_NxNyNz_Tz(ch2*4-1)
     Tz = channels_2bcross%config_NxNyNz_Tz(ch2*4)
     DO a = below_ef+1, tot_orbs
        nxa = all_orbit%nx(a)
        nya = all_orbit%ny(a)
        nza = all_orbit%nz(a)
        tza = all_orbit%tz(a)
        IF ( nxa - Nx /= nx_paeom ) cycle
        IF ( nya - Ny /= ny_paeom ) cycle
        IF ( nza - Nz /= nz_paeom ) cycle
        IF ( tza - 2*Tz /= tz_paeom ) cycle
        IF ( all_orbit%sz(a) == -1 ) r2_paeom_ind_cross(2*ch2-1) = a
        IF ( all_orbit%sz(a) == 1  ) r2_paeom_ind_cross(2*ch2)   = a
     end DO
  end DO

  ! setup 2p1h-cross-amplitudes
  ALLOCATE( r2_paeom_cross(channels_2bcross%number_confs) )
  ALLOCATE( r2_paeom_eqn_cross(channels_2bcross%number_confs) )
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_paeom_ind_cross(2*ch2) == 0 ) cycle
     ket_confs = number_2bcross(2,ch2)
     
     ALLOCATE( r2_paeom_cross(ch2)%cval(2,ket_confs) )
     ALLOCATE( r2_paeom_eqn_cross(ch2)%cval(2,ket_confs) )
     r2_paeom_cross(ch2)%cval = 0.d0
     r2_paeom_eqn_cross(ch2)%cval = 0.d0
     CALL mem_register('paeom2', REAL(2 * 2*ket_confs * 16.d0, dp))
  end DO

  CALL setup_proc_mappings_paeom
  IF ( iam == 0 ) write(6,'(A19,26x,I14)') 'PA-EOM2 dimension:', ndim
  ! [old memory print removed - replaced by mem_report]
  IF ( iam == 0 ) write(6,*)
  
end SUBROUTINE setup_paeom_amplitudes


SUBROUTINE setup_paeom3_amplitudes(ndim3)
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
  INTEGER :: a,b,c, k1,k2,k3,k4, aind
  INTEGER :: nx3,ny3,nz3,tz3, Nx2,Ny2,Nz2,Tz2
  INTEGER :: bra, bra_confs,ket_confs
  INTEGER :: bra_min,bra_max
  LOGICAL :: ch_added

  cut_3b = eom3_cut

  n3min = 1000
  n3max = -1000
  t3min = 1000
  t3max = -1000
  DO a = below_ef+1, tot_orbs-2
     DO b = a+1, tot_orbs-1
        DO c = b+1, tot_orbs
           IF ( all_orbit%nx(a) + all_orbit%nx(b) + all_orbit%nx(c) > n3max ) &
                n3max = all_orbit%nx(a) + all_orbit%nx(b) + all_orbit%nx(c)
           IF ( all_orbit%nx(a) + all_orbit%nx(b) + all_orbit%nx(c) < n3min ) &
                n3min = all_orbit%nx(a) + all_orbit%nx(b) + all_orbit%nx(c)
           IF ( all_orbit%tz(a) + all_orbit%tz(b) + all_orbit%tz(c) > t3max ) &
                t3max = all_orbit%tz(a) + all_orbit%tz(b) + all_orbit%tz(c)
           IF ( all_orbit%tz(a) + all_orbit%tz(b) + all_orbit%tz(c) < t3min ) &
                t3min = all_orbit%tz(a) + all_orbit%tz(b) + all_orbit%tz(c)
        end DO
     end DO
  end DO
  
  number_channels = 0
  DO nx3 = n3min, n3max
     DO ny3 = n3min, n3max
        DO nz3 = n3min, n3max
           DO tz3 = t3min, t3max, 2
              CALL number_3b_confs(nx3, ny3, nz3, tz3, 1, bra_confs) ! 1 for ppp
              IF ( bra_confs <= 0 ) cycle
              
              DO ch2 = 1, channels_2b%number_confs
                 IF ( number_2b(1,ch2) == 0 ) cycle
                 Nx2 = channels_2b%config_NxNyNz_Tz(ch2*4-3)
                 Ny2 = channels_2b%config_NxNyNz_Tz(ch2*4-2)
                 Nz2 = channels_2b%config_NxNyNz_Tz(ch2*4-1)
                 Tz2 = channels_2b%config_NxNyNz_Tz(ch2*4)
                 IF ( nx3 - Nx2 /= nx_paeom ) cycle
                 IF ( ny3 - Ny2 /= ny_paeom ) cycle
                 IF ( nz3 - Nz2 /= nz_paeom ) cycle
                 IF ( tz3 - 2*Tz2 /= tz_paeom ) cycle
                 number_channels = number_channels + 1
                 exit
              end DO
              
           end DO
        end DO
     end DO
  end DO
  channels_paeom3%number_confs = number_channels
  
  ALLOCATE( channels_paeom3%config_NxNyNz_Tz(4*number_channels) )
  CALL mem_register('paeom3', REAL(4 * number_channels * 4.d0, dp))
  ALLOCATE( ch2_paeom3(number_channels) )
  ch2_paeom3 = 0
  CALL mem_register('paeom3', REAL(number_channels * 4.d0, dp))
  
  number_channels = 0
  ndim3 = 0
  DO nx3 = n3min, n3max
     DO ny3 = n3min, n3max
        DO nz3 = n3min, n3max
           DO tz3 = t3min, t3max, 2
              CALL number_3b_confs(nx3, ny3, nz3, tz3, 1, bra_confs) ! 1 for ppp
              IF ( bra_confs <= 0 ) cycle
              
              ch_added = .FALSE.
              DO ch2 = 1, channels_2b%number_confs
                 IF ( number_2b(1,ch2) == 0 ) cycle
                 ket_confs = number_2b(1,ch2)
                 Nx2 = channels_2b%config_NxNyNz_Tz(ch2*4-3)
                 Ny2 = channels_2b%config_NxNyNz_Tz(ch2*4-2)
                 Nz2 = channels_2b%config_NxNyNz_Tz(ch2*4-1)
                 Tz2 = channels_2b%config_NxNyNz_Tz(ch2*4)
                 IF ( nx3 - Nx2 /= nx_paeom ) cycle
                 IF ( ny3 - Ny2 /= ny_paeom ) cycle
                 IF ( nz3 - Nz2 /= nz_paeom ) cycle
                 IF ( tz3 - 2*Tz2 /= tz_paeom ) cycle
                 ndim3 = ndim3 + bra_confs*ket_confs
                 number_channels = number_channels + 1
                 ch2_paeom3(number_channels) = ch2
                 ch_added = .TRUE.
                 exit
              end DO
              IF ( .not. ch_added ) cycle
              ch3 = number_channels
              
              k1 = ch3*4 - 3
              k2 = ch3*4 - 2
              k3 = ch3*4 - 1
              k4 = ch3*4
              channels_paeom3%config_NxNyNz_Tz(k1) = nx3
              channels_paeom3%config_NxNyNz_Tz(k2) = ny3
              channels_paeom3%config_NxNyNz_Tz(k3) = nz3
              channels_paeom3%config_NxNyNz_Tz(k4) = tz3

           end DO
        end DO
     end DO
  end DO
  IF ( iam == 0 ) write(6,'(A19,26x,I14)') 'PA-EOM3 dimension:', ndim3  
  
  
  ! Setup PA_EOM3 mapping
  CALL setup_proc_mappings_paeom3
  CALL mem_report('PA-EOM3 structures')


  ndim3 = 0
  ALLOCATE( r3_paeom(ch3_paeom_min:ch3_paeom_max) )
  DO ch3 = ch3_paeom_min, ch3_paeom_max
     ALLOCATE( r3_paeom(ch3)%val1(climits_paeom3(ch3,1):climits_paeom3(ch3,2)) )    
     DO aind = climits_paeom3(ch3,1), climits_paeom3(ch3,2)
        a         = clist_paeom3(ch3)%ival2(aind,1)
        ch1       = clist_paeom3(ch3)%ival2(aind,2)
        bra_min   = mapping_paeom3(ch3)%ival2(aind,1)
        bra_max   = mapping_paeom3(ch3)%ival2(aind,2)
        ch2       = ch2_paeom3(ch3)
        ket_confs = number_2b(1,ch2)
        ALLOCATE( r3_paeom(ch3)%val1(aind)%cval(bra_min:bra_max,ket_confs) )
        r3_paeom(ch3)%val1(aind)%cval = 0.d0
        CALL mem_register('paeom3', REAL((bra_max-bra_min+1)*ket_confs * 16.d0, dp))
        
        DO bra = bra_min, bra_max
           b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
           c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
           IF ( a >= b ) cycle
           ndim3 = ndim3 + ket_confs
        end DO
        
     end DO
  end DO
  CALL mpi_allreduce(mpi_in_place,ndim3,1,mpi_integer8,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_paeom3_amplitudes', 'allreduce')
  IF ( iam == 0 ) write(6,'(A26,19x,I14)') 'PA-EOM3 vector dimension:', ndim3
  CALL mem_report('PA-EOM3 amplitudes')
  
end SUBROUTINE setup_paeom3_amplitudes


SUBROUTINE populate_paeom(ndim1, ndim2, vecin)
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
  INTEGER :: bra,ket, bra_confs,ket_confs, bra_min,bra_max
  INTEGER :: a,b,c, j,k, aind
  INTEGER(i8) :: ii
  
  IF ( iam == 0 ) WRITE(6,*) "...Distributing right EOM vector..."
  
  r1_paeom = 0.d0
  ii = 1
  IF ( eom%my_start <= ii .and. ii <= eom%my_stop ) then
     r1_paeom(ii) = vecin(ii)
  end IF
  ii = ii + 1
  IF ( eom%my_start <= ii .and. ii <= eom%my_stop ) then
     r1_paeom(ii) = vecin(ii)
  end IF
  CALL mpi_allreduce(mpi_in_place,r1_paeom,size(r1_paeom),mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'populate_paeom', 'allreduce')
  
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     r2_paeom(ch)%cval = 0.d0
     DO bra = 1, bra_confs
        DO ket = 1, 2
           ii = ii + 1
           IF ( eom%my_start <= ii .and. ii <= eom%my_stop ) then
              r2_paeom(ch)%cval(bra,ket) = vecin(ii)
           end IF
        end DO
     end DO
     CALL mpi_allreduce(mpi_in_place,r2_paeom(ch)%cval,size(r2_paeom(ch)%cval),mpi_complex16,mpi_sum,mpi_comm_world,ierror)
     CALL check_mpi(ierror, 'populate_paeom', 'allreduce')
  end DO
  
  CALL r2_paeom_cross_recouple
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  IF ( eom_approx > 0 ) then

     ii = max(eom_ndim2, eom%my_start-1)
     DO ch3 = ch3_paeom_min, ch3_paeom_max
        ch2       = ch2_paeom3(ch3)
        ket_confs = number_2b(1,ch2)
        DO aind = climits_paeom3(ch3,1), climits_paeom3(ch3,2)
           a       = clist_paeom3(ch3)%ival2(aind,1)
           ch1     = clist_paeom3(ch3)%ival2(aind,2)
           bra_min = mapping_paeom3(ch3)%ival2(aind,1)
           bra_max = mapping_paeom3(ch3)%ival2(aind,2)
           
           r3_paeom(ch3)%val1(aind)%cval = 0.d0
           
           DO bra = bra_min, bra_max
              b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
              c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
              IF ( a >= b ) cycle
              DO ket = 1, ket_confs
                 j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
                 k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
                 ii = ii + 1
                 r3_paeom(ch3)%val1(aind)%cval(bra,ket) = vecin(ii)
              end DO
           end DO
        end DO
     end DO
     CALL mpi_barrier(mpi_comm_world,ierror)
     CALL expand_paeom_3p2h
     
  end IF
  
  IF ( test == 4 ) CALL build_paeom_test

end SUBROUTINE populate_paeom


SUBROUTINE populate_paeom_vec(ndim1, ndim2, vecout)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER(i8), INTENT(in) :: ndim1, ndim2
  COMPLEX(dpc), INTENT(out) :: vecout(ndim1:ndim2)
  INTEGER :: ch, ch3, ch1,ch2, a,b,c,j,k, aind
  INTEGER :: bra_min,bra_max, bra,ket
  INTEGER :: bra_confs,ket_confs
  INTEGER(i8) :: ii

  IF ( iam == 0 ) WRITE(6,*) "...Collecting right EOM vector..."
  
  CALL r2_paeom_add_cross
  CALL mpi_barrier(mpi_comm_world,ierror)
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,r2_paeom_eqn(ch)%cval,size(r2_paeom_eqn(ch)%cval),mpi_complex16,mpi_sum,mpi_comm_world,ierror)
     CALL check_mpi(ierror, 'populate_paeom_vec', 'allreduce')
  end DO
  
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF ( test == 4 ) then
     CALL build_r1_paeom_eqn_test
     CALL build_r2_paeom_eqn_test
     IF ( eom_approx > 0 ) CALL build_r3_paeom_eqn_test
  end IF
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  vecout = 0.d0

  ii = 1
  IF ( eom%my_start <= ii .and. ii <= eom%my_stop ) then
     vecout(ii) = r1_paeom_eqn(ii)
  end IF
  ii = ii + 1
  IF ( eom%my_start <= ii .and. ii <= eom%my_stop ) then
     vecout(ii) = r1_paeom_eqn(ii)
  end IF

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     DO bra = 1, bra_confs
        DO ket = 1, 2
           ii = ii + 1
           IF ( eom%my_start <= ii .and. ii <= eom%my_stop ) then
              vecout(ii) = r2_paeom_eqn(ch)%cval(bra,ket)
           end IF
        end DO
     end DO
  end DO

  
  IF ( eom_approx > 0 ) then

     ii = max(eom_ndim2, eom%my_start-1)
     DO ch3 = ch3_paeom_min, ch3_paeom_max
        ch2       = ch2_paeom3(ch3)
        ket_confs = number_2b(1,ch2)
        DO aind = climits_paeom3(ch3,1), climits_paeom3(ch3,2)
           a       = clist_paeom3(ch3)%ival2(aind,1)
           ch1     = clist_paeom3(ch3)%ival2(aind,2)
           bra_min = mapping_paeom3(ch3)%ival2(aind,1)
           bra_max = mapping_paeom3(ch3)%ival2(aind,2)           
           DO bra = bra_min, bra_max
              b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
              c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
              IF ( a >= b ) cycle
              DO ket = 1, ket_confs
                 j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
                 k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
                 ii = ii + 1
                 vecout(ii) = r3_paeom(ch3)%val1(aind)%cval(bra,ket)
              end DO
           end DO
        end DO
     end DO
     
  end IF

end SUBROUTINE populate_paeom_vec


SUBROUTINE expand_paeom_3p2h
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE parallel
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch3,ch1,ch2, bra,bra1
  INTEGER :: bra_confs,ket_confs
  INTEGER :: bra_min,bra_max, a,b,c
  INTEGER :: aind,bind,cind1, ch_b
  INTEGER, allocatable :: clist_inv(:)
  TYPE (superblock_storage) :: r3_temp
  INTEGER :: group0
  
  IF ( iam == 0 ) WRITE(6,*) "...Expanding R3..."

  ALLOCATE( clist_inv(below_ef+1:tot_orbs) )
  
  DO ch3 = ch3_paeom_min, ch3_paeom_max        
     ch2       = ch2_paeom3(ch3)
     ket_confs = number_2b(1,ch2)
     
     ! Build inverse c map
     clist_inv = 0
     DO aind = 1, climit_paeom3(ch3)
        a = clist_paeom3(ch3)%ival2(aind,1)
        clist_inv(a) = aind
     end DO
     
     ! Allocate r3_temp for all c_channels
     ALLOCATE( r3_temp%val1(climit_paeom3(ch3)) )
     DO aind = 1, climit_paeom3(ch3)
        ch1       = clist_paeom3(ch3)%ival2(aind,2)
        bra_confs = number_2b(3,ch1)
        ALLOCATE( r3_temp%val1(aind)%cval(bra_confs,ket_confs) )
        r3_temp%val1(aind)%cval = 0.d0
     end DO
     
     ! Fill r3_temp for local c_channels
     DO aind = climits_paeom3(ch3,1), climits_paeom3(ch3,2)
        bra_min = mapping_paeom3(ch3)%ival2(aind,1)
        bra_max = mapping_paeom3(ch3)%ival2(aind,2)
        r3_temp%val1(aind)%cval(bra_min:bra_max,:) = r3_paeom(ch3)%val1(aind)%cval(bra_min:bra_max,:)
     end DO
     
     ! Allreduce r3_temp for all c_channels
     group0 = group_paeom3(iam+1, ch3)
     IF ( group0 > 0 ) then
        DO aind = 1, climit_paeom3(ch3)
           CALL mpi_allreduce(mpi_in_place,r3_temp%val1(aind)%cval,size(r3_temp%val1(aind)%cval),&
                mpi_complex16,mpi_sum,subcomm_paeom3(group0),ierror)
           CALL check_mpi(ierror, 'expand_paeom_3p2h', 'allreduce')
        end DO
     end IF

     DO aind = climits_paeom3(ch3,1), climits_paeom3(ch3,2)
        a         = clist_paeom3(ch3)%ival2(aind,1)
        ch1       = clist_paeom3(ch3)%ival2(aind,2)
        bra_min   = mapping_paeom3(ch3)%ival2(aind,1)
        bra_max   = mapping_paeom3(ch3)%ival2(aind,2)
        bra_confs = number_2b(3,ch1)
        DO bra = bra_min, bra_max
           b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
           c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
           IF ( a <= b ) cycle ! a(bc)

           bind = clist_inv(b)
           ch_b = clist_paeom3(ch3)%ival2(bind,2)
           IF ( ch_b <= 0 ) cycle
           IF ( a == c ) cycle           
           IF ( a > c ) then ! b(ca)
              bra1 = pp_config_2b%ival2(c,a)
              r3_paeom(ch3)%val1(aind)%cval(bra,:) = r3_paeom(ch3)%val1(aind)%cval(bra,:) &
                   + r3_temp%val1(bind)%cval(bra1,:)
           ELSE ! b(ac)
              bra1 = pp_config_2b%ival2(a,c)
              r3_paeom(ch3)%val1(aind)%cval(bra,:) = r3_paeom(ch3)%val1(aind)%cval(bra,:) &
                   - r3_temp%val1(bind)%cval(bra1,:)
           end IF
           
        end DO
     end DO
     
     DO cind1 = 1, climit_paeom3(ch3)
        DEALLOCATE( r3_temp%val1(cind1)%cval )
     end DO
     DEALLOCATE( r3_temp%val1 )     
  end DO
  DEALLOCATE( clist_inv )

end SUBROUTINE expand_paeom_3p2h


SUBROUTINE r2_paeom_cross_recouple
  use parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch,ch2, bra, bra2,ket2
  INTEGER :: ket_min,ket_max
  INTEGER :: a,b,i
  REAL(dp) :: phase

  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_paeom_ind_cross(2*ch2) == 0 ) cycle
     
     r2_paeom_cross(ch2)%cval = 0.d0
     
     IF ( check_my_channel_paeom_php_cross(ch2) == 0 ) cycle
     ket_min = mapping_paeom_php_cross(iam+1,ch2,2)
     ket_max = mapping_paeom_php_cross(iam+1,ch2,3)

     DO bra2 = 1, 2
        a = r2_paeom_ind_cross(2*(ch2-1)+bra2)
        DO ket2 = ket_min, ket_max
           i   = lookup_2bcross_configs(2,ch2)%ival2(1,ket2)
           b   = lookup_2bcross_configs(2,ch2)%ival2(2,ket2)

           ch = pp_channel_2b%ival2(a,b)
           IF ( ch == 0 ) cycle
           phase = 1
           bra = pp_config_2b%ival2(a,b)
           IF ( b < a ) phase = -phase
           IF ( all_orbit%sz(i) == -1 ) then
              r2_paeom_cross(ch2)%cval(bra2,ket2) = phase*r2_paeom(ch)%cval(bra,1)
           ELSE
              r2_paeom_cross(ch2)%cval(bra2,ket2) = phase*r2_paeom(ch)%cval(bra,2)
           end IF
        end DO
     end DO
  end DO  
  
  CALL mpi_barrier(mpi_comm_world,ierror)
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_paeom_ind_cross(2*ch2) == 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,r2_paeom_cross(ch2)%cval,size(r2_paeom_cross(ch2)%cval),&
          mpi_complex16,mpi_sum,mpi_comm_world,ierror)
     CALL check_mpi(ierror, 'r2_paeom_cross_recouple', 'allreduce')
  end DO
  
END SUBROUTINE r2_paeom_cross_recouple


SUBROUTINE r2_paeom_add_cross
  use parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE ang_mom_functions
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch,ch2, bra,ket, ket2
  INTEGER :: bra_min,bra_max
  INTEGER :: a,b,i
  
  CALL mpi_barrier(mpi_comm_world,ierror)
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_paeom_ind_cross(2*ch2) == 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,r2_paeom_eqn_cross(ch2)%cval,size(r2_paeom_eqn_cross(ch2)%cval),&
          mpi_complex16,mpi_sum,mpi_comm_world,ierror)
     CALL check_mpi(ierror, 'r2_paeom_add_cross', 'allreduce')
  end DO

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     
     IF ( check_my_channel_paeom_pph(ch) == 0 ) cycle
     bra_min = mapping_paeom_pph(iam+1,ch,2)
     bra_max = mapping_paeom_pph(iam+1,ch,3)

     DO bra = bra_min, bra_max
        a   = lookup_2b_configs(3,ch)%ival2(1,bra)
        b   = lookup_2b_configs(3,ch)%ival2(2,bra)
        DO ket = 1, 2
           i = r2_paeom_ind(2*(ch-1)+ket)
           ch2 = hp_channel_2bcross%ival2(i,b)
           IF ( ch2 == 0 ) cycle
           ket2 = hp_config_2bcross%ival2(i,b)
           IF ( all_orbit%sz(a) == -1 ) then
              r2_paeom_eqn(ch)%cval(bra,ket) = r2_paeom_eqn(ch)%cval(bra,ket) + r2_paeom_eqn_cross(ch2)%cval(1,ket2)
           ELSE
              r2_paeom_eqn(ch)%cval(bra,ket) = r2_paeom_eqn(ch)%cval(bra,ket) + r2_paeom_eqn_cross(ch2)%cval(2,ket2)
           end IF
        end DO
        ! P(ab)
        DO ket = 1, 2
           i = r2_paeom_ind(2*(ch-1)+ket)
           ch2 = hp_channel_2bcross%ival2(i,a)
           IF ( ch2 == 0 ) cycle
           ket2 = hp_config_2bcross%ival2(i,a)
           IF ( all_orbit%sz(b) == -1 ) then
              r2_paeom_eqn(ch)%cval(bra,ket) = r2_paeom_eqn(ch)%cval(bra,ket) - r2_paeom_eqn_cross(ch2)%cval(1,ket2)
           ELSE
              r2_paeom_eqn(ch)%cval(bra,ket) = r2_paeom_eqn(ch)%cval(bra,ket) - r2_paeom_eqn_cross(ch2)%cval(2,ket2)
           end IF
        end DO
     end DO
  end DO           

END SUBROUTINE r2_paeom_add_cross


SUBROUTINE deallocate_paeom
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch1, ch2, ch3, cind1

  CALL deallocate_hbar1b_I2
  CALL deallocate_hbar1b_I3
  CALL deallocate_hbar2b_I3pa
  CALL deallocate_hbar2b_I5pa_cross
  CALL deallocate_hbar2b_I6pa
  CALL deallocate_hbar3b_paeom
  
  DEALLOCATE( paeom_eigs )
  
  DEALLOCATE( r1_paeom_ind )
  DEALLOCATE( r1_paeom )
  DEALLOCATE( r1_paeom_eqn )
  
  DO ch1 = 1, channels_2b%number_confs
     IF ( number_2b(3,ch1) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch1) == 0 ) cycle
     DEALLOCATE( r2_paeom(ch1)%cval )
     DEALLOCATE( r2_paeom_eqn(ch1)%cval )
  end DO
  DEALLOCATE( r2_paeom )
  DEALLOCATE( r2_paeom_eqn )
  DEALLOCATE( r2_paeom_ind ) !!

  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_paeom_ind_cross(2*ch2) == 0 ) cycle
     DEALLOCATE( r2_paeom_cross(ch2)%cval )
     DEALLOCATE( r2_paeom_eqn_cross(ch2)%cval )
  end DO
  DEALLOCATE( r2_paeom_cross )
  DEALLOCATE( r2_paeom_eqn_cross )
  DEALLOCATE( r2_paeom_ind_cross ) !!

  DEALLOCATE( mapping_paeom_pph )
  DEALLOCATE( check_my_channel_paeom_pph )
  DEALLOCATE( mapping_paeom_php_cross )
  DEALLOCATE( check_my_channel_paeom_php_cross )

  DEALLOCATE( eom%all_starts )
  DEALLOCATE( eom%all_stops )

  IF ( eom_approx > 0 ) then
     DEALLOCATE( channels_paeom3%config_NxNyNz_Tz )
     DEALLOCATE( ch2_paeom3 )
     DEALLOCATE( group_paeom3 )
     DEALLOCATE( subcomm_paeom3 )
  
     DO ch3 = ch3_paeom_min, ch3_paeom_max
        DO cind1 = climits_paeom3(ch3,1), climits_paeom3(ch3,2)
           DEALLOCATE( r3_paeom(ch3)%val1(cind1)%cval )
        end DO
        DEALLOCATE( r3_paeom(ch3)%val1 )
     end DO
     DEALLOCATE( r3_paeom )

     DO ch3 = ch3_paeom_min, ch3_paeom_max
        DEALLOCATE( clist_paeom3(ch3)%ival2 )     
     end DO
     DEALLOCATE( clist_paeom3 )
     DEALLOCATE( climit_paeom3 )

     DO ch3 = ch3_paeom_min, ch3_paeom_max
        IF ( ALLOCATED(mapping_paeom3(ch3)%ival2) ) then
           DEALLOCATE( mapping_paeom3(ch3)%ival2 )
        end IF
     end DO
     DEALLOCATE( climits_paeom3 )
     DEALLOCATE( mapping_paeom3 )
  end IF
  
end SUBROUTINE deallocate_paeom
