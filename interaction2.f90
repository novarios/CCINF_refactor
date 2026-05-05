

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE setup_2b_interaction
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE chiral_potentials
  USE mem_tracker
  ! USE NNForceCartesian
  
  IMPLICIT NONE
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max, bra,ket,ket0
  INTEGER :: ch_v2, v2bra,v2ket, size1
  INTEGER :: a,b,c,d, i,j,k,l, h,p, phase, count1, count2
  COMPLEX(dpc) :: v2b, v3b, vtmp, e0_1,e0_2,e0_3, e_tmp
  COMPLEX(dpc), ALLOCATABLE :: fock_mtx0(:,:), fock_tmp(:,:)

  !
  ! Allocate and Initialize Structures
  !

  ! IF ( interaction(1:2) == "pw" ) then
  !    IF ( iam == 0 ) write(6,*) " n1max = ", n1max
  !    IF ( iam == 0 ) write(6,*) " N2max = ", N2max
  !    call vrel%init( Lbox=lx, Nmax=n1max, N2max=N2max )
  !    call vrel%SetNNCartesian( interaction, "bare", -1.d0, Jmax=12, pmax=12.d0, NMesh=200 )
  !    call mpi_barrier(mpi_comm_world, ierror)
  ! end IF

  
  ! Vacuum Energy
  e0 = 0.d0
  e0_1 = 0.d0
  e0_2 = 0.d0
  e0_3 = 0.d0

  ! Fock matrix  
  ALLOCATE( fock_mtx(all_orbit%total_orbits, all_orbit%total_orbits) )
  fock_mtx = 0.d0
  CALL mem_register('v2', REAL((all_orbit%total_orbits)**2 * 16.d0, dp))
  ALLOCATE( fock_mtx0(all_orbit%total_orbits, all_orbit%total_orbits) )
  fock_mtx0 = 0.d0
  ! For OMP parallelization
  ALLOCATE( fock_tmp(all_orbit%total_orbits, all_orbit%total_orbits) )
  fock_tmp = 0.d0
  

  ! hhhh
  ALLOCATE( v2b_hhhh(channels_2b%number_confs) )
  DO ch = 1, channels_2b%number_confs
     bra_confs = number_2b(1,ch)
     IF ( bra_confs <= 0 ) cycle
     ALLOCATE( v2b_hhhh(ch)%cval(bra_confs,bra_confs) )
     v2b_hhhh(ch)%cval = 0.d0
     CALL mem_register('v2', REAL((bra_confs)**2 * 16.d0, dp))
  end DO

  ! hphh
  ALLOCATE( v2b_hphh(channels_2b%number_confs) )
  DO ch = 1, channels_2b%number_confs
     bra_confs = number_2b(2,ch)
     ket_confs = number_2b(1,ch)
     IF ( bra_confs * ket_confs <= 0 ) cycle
     ALLOCATE( v2b_hphh(ch)%cval(bra_confs,ket_confs) )
     v2b_hphh(ch)%cval = 0.d0
     CALL mem_register('v2', REAL(bra_confs*ket_confs * 16.d0, dp))
  end DO

  ! hphp
  ALLOCATE( v2b_hphp(channels_2b%number_confs) )
  DO ch = 1, channels_2b%number_confs
     bra_confs = number_2b(2,ch)
     IF ( bra_confs <= 0 ) cycle
     ALLOCATE( v2b_hphp(ch)%cval(bra_confs,bra_confs) )
     v2b_hphp(ch)%cval = 0.d0
     CALL mem_register('v2', REAL((bra_confs)**2 * 16.d0, dp))
  end DO

  ! pphp
  ALLOCATE( v2b_pphp(channels_2b%number_confs) )
  DO ch = 1, channels_2b%number_confs
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(2,ch)
     IF ( bra_confs * ket_confs <= 0 ) cycle
     ALLOCATE( v2b_pphp(ch)%cval(bra_confs,ket_confs) )
     v2b_pphp(ch)%cval = 0.d0
     CALL mem_register('v2', REAL(bra_confs*ket_confs * 16.d0, dp))
  end DO

  ! pphh
  ALLOCATE( v2b_pphh(channels_2b%number_confs) )
  DO ch = 1, channels_2b%number_confs
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)
     IF ( bra_confs * ket_confs <= 0 ) cycle
     ALLOCATE( v2b_pphh(ch)%cval(bra_confs,ket_confs) )
     v2b_pphh(ch)%cval = 0.d0
     CALL mem_register('v2', REAL(bra_confs*ket_confs * 16.d0, dp))
  end DO

  ! pppp
  ALLOCATE( v2b_pppp(channels_2b%number_confs) )
  DO ch = 1, channels_2b%number_confs
     IF ( check_my_channel_v2b_pppp(ch) == 0 ) cycle
     bra_min = mapping_v2b_pppp(iam+1,ch,2)
     bra_max = mapping_v2b_pppp(iam+1,ch,3)
     bra_confs = number_2b(3,ch)
     ALLOCATE( v2b_pppp(ch)%cval(bra_min:bra_max,bra_confs) )
     v2b_pppp(ch)%cval = 0.d0
     CALL mem_register('v2', REAL((bra_max-bra_min+1)*bra_confs * 16.d0, dp))
  end DO

  
  !
  ! 3N contributions
  !
  IF ( tnf_approx > 0 ) then
     IF ( iam == 0 ) write(6,*) '...Normal-Ordering 3NF...'
     ! CALL v2b_hphp_cross_3b(fock_mtx0)
     CALL v2b_hhhh_3b(fock_mtx0, e0_3)
     IF ( iam == 0 ) write(6,*) 'hhhh done.'
     CALL v2b_hphh_3b(fock_mtx0)
     IF ( iam == 0 ) write(6,*) 'hphh done.'
     CALL v2b_hphp_3b(fock_mtx0)
     IF ( iam == 0 ) write(6,*) 'hphp done.'
     CALL v2b_pphp_3b
     IF ( iam == 0 ) write(6,*) 'pphp done.'
     CALL v2b_pphh_3b
     IF ( iam == 0 ) write(6,*) 'pphh done.'
     CALL v2b_pppp_3b
     IF ( iam == 0 ) write(6,*) 'pppp done.'
  end IF
    
  !
  ! 2N contributions  
  !
  IF ( iam == 0 ) write(6,*) '...Normal-Ordering 2NF...'
  ! CALL v2b_hphp_cross_2b(fock_mtx0)
  CALL v2b_hhhh_2b(fock_mtx0, e0_2)
  IF ( iam == 0 ) write(6,*) 'hhhh done.'
  CALL v2b_hphh_2b(fock_mtx0)
  IF ( iam == 0 ) write(6,*) 'hphh done.'
  CALL v2b_hphp_2b(fock_mtx0)
  IF ( iam == 0 ) write(6,*) 'hphp done.'
  CALL v2b_pphh_2b
  IF ( iam == 0 ) write(6,*) 'pphh done.'
  CALL v2b_pphp_2b
  IF ( iam == 0 ) write(6,*) 'pphp done.'
  CALL v2b_pppp_2b
  IF ( iam == 0 ) write(6,*) 'pppp done.'
  ! !
  ! IF ( interaction(1:2) == "pw" ) call vrel%fin()

  
  ! Reduce Arrays
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(2,ch) <= 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,v2b_hphp(ch)%cval,size(v2b_hphp(ch)%cval),&
          mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  end DO
  IF ( iam == 0 ) write(6,*) 'hphp reduction done.'
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) <= 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,v2b_hhhh(ch)%cval,size(v2b_hhhh(ch)%cval),&
          mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  end DO
  IF ( iam == 0 ) write(6,*) 'hhhh reduction done.'
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(2,ch) * number_2b(1,ch) <= 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,v2b_hphh(ch)%cval,size(v2b_hphh(ch)%cval),&
          mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  end DO
  IF ( iam == 0 ) write(6,*) 'hphh reduction done.'
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) * number_2b(2,ch) <= 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,v2b_pphp(ch)%cval,size(v2b_pphp(ch)%cval),&
          mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  end DO
  IF ( iam == 0 ) write(6,*) 'pphp reduction done.'
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) * number_2b(1,ch) <= 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,v2b_pphh(ch)%cval,size(v2b_pphh(ch)%cval),&
          mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  end DO
  IF ( iam == 0 ) write(6,*) 'pphh reduction done.'
  IF ( iam == 0 ) write(6,*)

  
  !
  ! Compute Fock Matrix
  !
  CALL mpi_allreduce(mpi_in_place, fock_mtx0, size(fock_mtx0), mpi_complex16, mpi_sum, mpi_comm_world, ierror)
  fock_mtx = fock_mtx + fock_mtx0
  DEALLOCATE(fock_mtx0)
  DO p = 1, tot_orbs
     fock_mtx(p,p) = fock_mtx(p,p) + all_orbit%e(p)
  end DO

  !
  ! Compute Vacuum Energy
  !
  CALL mpi_allreduce(mpi_in_place, e0_3, 1, mpi_complex16, mpi_sum, mpi_comm_world, ierror)
  CALL mpi_allreduce(mpi_in_place, e0_2, 1, mpi_complex16, mpi_sum, mpi_comm_world, ierror)
  DO h = 1, below_ef
     e0_1 = e0_1 + all_orbit%e(h)
  end DO
  e0 = e0_1 + e0_2 + e0_3
  
  IF ( iam == 0 ) write(6,*)
  IF ( iam == 0 ) write(6,*) 'E0_1b/A = ', e0_1/below_ef
  IF ( iam == 0 ) write(6,*) 'E0_2b/A = ', e0_2/below_ef
  IF ( iam == 0 ) write(6,*) 'E0_3b/A = ', e0_3/below_ef
  IF ( iam == 0 ) write(6,*) 'E0/A    = ', e0/below_ef
  IF ( iam == 0 ) write(6,*)


  ! hphp <-- v2b_hphp
  ALLOCATE( v2b_hphp_cross(channels_2bcross%number_confs) )
  DO ch = 1, channels_2bcross%number_confs
     bra_confs = number_2bcross(2,ch)
     IF ( bra_confs <= 0 ) cycle
     ALLOCATE( v2b_hphp_cross(ch)%cval(bra_confs,bra_confs) )
     v2b_hphp_cross(ch)%cval = 0.d0
     CALL mem_register('v2', REAL(bra_confs**2 * 16.d0, dp))
     IF ( check_my_channel_v2b_hphp_cross(ch) == 0 ) cycle
     bra_min = mapping_v2b_hphp_cross(iam+1,ch,2)
     bra_max = mapping_v2b_hphp_cross(iam+1,ch,3)
     ket_confs = number_2bcross(2,ch)
     
     DO bra = bra_min, bra_max
        i   = lookup_2bcross_configs(2,ch)%ival2(1,bra)
        b   = lookup_2bcross_configs(2,ch)%ival2(2,bra)
        DO ket = 1, ket_confs
           j   = lookup_2bcross_configs(2,ch)%ival2(1,ket)
           a   = lookup_2bcross_configs(2,ch)%ival2(2,ket)

           ch_v2 = hp_channel_2b%ival2(i,a)
           IF ( ch_v2 == 0 ) cycle
           IF ( ch_v2 /= hp_channel_2b%ival2(j,b) ) cycle
           
           v2bra = hp_config_2b%ival2(i,a)
           v2ket = hp_config_2b%ival2(j,b)
           v2b_hphp_cross(ch)%cval(bra,ket) = v2b_hphp(ch_v2)%cval(v2bra,v2ket)
        end DO
     end DO
  end DO
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(2,ch) <= 0 ) cycle
     DEALLOCATE( v2b_hphp(ch)%cval )
  end DO
  DEALLOCATE( v2b_hphp )
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch) <= 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,v2b_hphp_cross(ch)%cval,size(v2b_hphp_cross(ch)%cval),&
          mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  end DO
  
  ! phhp <-- v2b_pphh
  ALLOCATE( v2b_phhp_cross(channels_2bcross%number_confs) )
  DO ch = 1, channels_2bcross%number_confs
     bra_confs = number_2bcross(3,ch)
     ket_confs = number_2bcross(2,ch)
     IF ( bra_confs * ket_confs <= 0 ) cycle
     ALLOCATE( v2b_phhp_cross(ch)%cval(bra_confs,ket_confs) )
     v2b_phhp_cross(ch)%cval = 0.d0
     CALL mem_register('v2', REAL(bra_confs*ket_confs * 16.d0, dp))
     IF ( check_my_channel_v2b_hphp_cross(ch) == 0 ) cycle
     bra_min = mapping_v2b_phhp_cross(iam+1,ch,2)
     bra_max = mapping_v2b_phhp_cross(iam+1,ch,3)
     ket_confs = number_2bcross(2,ch)
     
     DO bra = bra_min, bra_max
        a   = lookup_2bcross_configs(3,ch)%ival2(1,bra)
        j   = lookup_2bcross_configs(3,ch)%ival2(2,bra)
        DO ket = 1, ket_confs
           i   = lookup_2bcross_configs(2,ch)%ival2(1,ket)
           b   = lookup_2bcross_configs(2,ch)%ival2(2,ket)

           ch_v2 = pp_channel_2b%ival2(a,b)
           IF ( ch_v2 == 0 ) cycle
           IF ( ch_v2 /= hh_channel_2b%ival2(i,j) ) cycle
           
           phase = 1
           v2bra = pp_config_2b%ival2(a,b)
           v2ket = hh_config_2b%ival2(i,j)
           IF ( b < a ) phase = -phase
           IF ( j < i ) phase = -phase
           v2b_phhp_cross(ch)%cval(bra,ket) = phase * v2b_pphh(ch_v2)%cval(v2bra,v2ket)
        end DO
     end DO
  end DO
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(3,ch) * number_2bcross(2,ch) <= 0 ) cycle
     CALL mpi_allreduce(mpi_in_place,v2b_phhp_cross(ch)%cval,size(v2b_phhp_cross(ch)%cval),&
          mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  end DO
  
end SUBROUTINE setup_2b_interaction




SUBROUTINE v2b_hphp_2b(fock_mtx0)
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE chiral_potentials

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: fock_mtx0(all_orbit%total_orbits, all_orbit%total_orbits)
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max, bra,ket,ket0
  INTEGER :: ch_v2, v2bra,v2ket, size1, count, ii
  INTEGER :: i,a,j,b, phase
  COMPLEX(dpc) :: v2b
  INTEGER, ALLOCATABLE :: loop_inds(:,:)
  COMPLEX(dpc), ALLOCATABLE :: loop_vals(:)
  
  DO ch = 1, channels_2b%number_confs
     IF ( check_my_channel_v2b_hphp(ch) == 0 ) cycle
     bra_min = mapping_v2b_hphp(iam+1,ch,2)
     bra_max = mapping_v2b_hphp(iam+1,ch,3)
     bra_confs = number_2bcross(2,ch)
     
     count = 0
     DO bra = bra_min, bra_max
        size1 = bra + bra_confs-1
        DO ket0 = bra, size1
           ket = mod(ket0-1,bra_confs) + 1
           IF ( bra_min <= ket .and. ket < bra ) cycle
           count = count + 1
        end DO
     end DO
     ALLOCATE( loop_inds(count, 2) )
     ALLOCATE( loop_vals(count) )
     loop_inds = 0
     loop_vals = 0.d0
     count = 0
     DO bra = bra_min, bra_max
        size1 = bra + bra_confs-1
        DO ket0 = bra, size1
           ket = mod(ket0-1,bra_confs) + 1
           IF ( bra_min <= ket .and. ket < bra ) cycle
           count = count + 1
           loop_inds(count,1) = bra
           loop_inds(count,2) = ket
        end DO
     end DO
     
     !$omp parallel default(shared) private(ii, bra,ket,i,a,j,b, v2b)
     !$omp do schedule(static)
     DO ii = 1, count
        bra = loop_inds(ii, 1)
        ket = loop_inds(ii, 2)
        i   = lookup_2b_configs(2,ch)%ival2(1,bra)
        a   = lookup_2b_configs(2,ch)%ival2(2,bra)
        j   = lookup_2b_configs(2,ch)%ival2(1,ket)
        b   = lookup_2b_configs(2,ch)%ival2(2,ket)
        loop_vals(ii) = v2int(i,a,j,b)
     end DO
     !$omp end do
     !$omp end parallel

     DO ii = 1, count
        bra = loop_inds(ii, 1)
        ket = loop_inds(ii, 2)
        i   = lookup_2b_configs(2,ch)%ival2(1,bra)
        a   = lookup_2b_configs(2,ch)%ival2(2,bra)
        j   = lookup_2b_configs(2,ch)%ival2(1,ket)
        b   = lookup_2b_configs(2,ch)%ival2(2,ket)
        v2b = loop_vals(ii)     
        v2b_hphp(ch)%cval(bra,ket) = v2b_hphp(ch)%cval(bra,ket) + v2b
        ! Compute fock_mtx
        IF ( i == j ) then
           fock_mtx0(a,b) = fock_mtx0(a,b) + v2b
           IF ( bra < ket .and. ket <= bra_max ) then
              fock_mtx0(b,a) = fock_mtx0(b,a) + conjg(v2b)
           end IF
        end IF
     end DO
     
     DO bra = bra_min, bra_max
        DO ket = bra_min, bra-1
           v2b_hphp(ch)%cval(bra,ket) = conjg(v2b_hphp(ch)%cval(ket,bra))
        end DO
     end DO
     DEALLOCATE( loop_inds, loop_vals )
  end DO

end SUBROUTINE v2b_hphp_2b

SUBROUTINE v2b_hhhh_2b(fock_mtx0, e0_2)
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE chiral_potentials

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: e0_2
  COMPLEX(dpc), INTENT(INOUT) :: fock_mtx0(all_orbit%total_orbits, all_orbit%total_orbits)
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max, bra,ket,ket0
  INTEGER :: ch_v2, v2bra,v2ket, size1, count, ii
  INTEGER :: i,j,k,l, phase
  COMPLEX(dpc) :: v2b
  INTEGER, ALLOCATABLE :: loop_inds(:,:)
  COMPLEX(dpc), ALLOCATABLE :: loop_vals(:)
  
  DO ch = 1, channels_2b%number_confs
     IF ( check_my_channel_v2b_hhhh(ch) == 0 ) cycle
     bra_min = mapping_v2b_hhhh(iam+1,ch,2)
     bra_max = mapping_v2b_hhhh(iam+1,ch,3)
     bra_confs = number_2b(1,ch)

     count = 0
     DO bra = bra_min, bra_max
        size1 = bra + bra_confs-1
        DO ket0 = bra, size1
           ket = mod(ket0-1,bra_confs) + 1
           IF ( bra_min <= ket .and. ket < bra ) cycle
           count = count + 1
        end DO
     end DO
     ALLOCATE( loop_inds(count, 2) )
     ALLOCATE( loop_vals(count) )
     loop_inds = 0
     loop_vals = 0.d0
     count = 0
     DO bra = bra_min, bra_max
        size1 = bra + bra_confs-1
        DO ket0 = bra, size1
           ket = mod(ket0-1,bra_confs) + 1
           IF ( bra_min <= ket .and. ket < bra ) cycle
           count = count + 1
           loop_inds(count,1) = bra
           loop_inds(count,2) = ket
        end DO
     end DO

     !$omp parallel default(shared) private(ii, bra,ket,i,j,k,l, v2b)
     !$omp do schedule(static)
     DO ii = 1, count
        bra = loop_inds(ii, 1)
        ket = loop_inds(ii, 2)
        i   = lookup_2b_configs(1,ch)%ival2(1,bra)
        j   = lookup_2b_configs(1,ch)%ival2(2,bra)
        k   = lookup_2b_configs(1,ch)%ival2(1,ket)
        l   = lookup_2b_configs(1,ch)%ival2(2,ket)
        loop_vals(ii) = v2int(i,j,k,l)
     end DO
     !$omp end do
     !$omp end parallel
     
     DO ii = 1, count
        bra = loop_inds(ii, 1)
        ket = loop_inds(ii, 2)
        i   = lookup_2b_configs(1,ch)%ival2(1,bra)
        j   = lookup_2b_configs(1,ch)%ival2(2,bra)
        k   = lookup_2b_configs(1,ch)%ival2(1,ket)
        l   = lookup_2b_configs(1,ch)%ival2(2,ket)
        v2b = loop_vals(ii)
        v2b_hhhh(ch)%cval(bra,ket) = v2b_hhhh(ch)%cval(bra,ket) + v2b
        ! Compute fock_mtx
        IF ( i == k ) fock_mtx0(j,l) = fock_mtx0(j,l) + v2b
        IF ( i == l ) fock_mtx0(j,k) = fock_mtx0(j,k) - v2b
        IF ( j == l ) fock_mtx0(i,k) = fock_mtx0(i,k) + v2b
        IF ( j == k ) fock_mtx0(i,l) = fock_mtx0(i,l) - v2b
        IF ( bra < ket .and. ket <= bra_max ) then
           IF ( i == k ) fock_mtx0(l,j) = fock_mtx0(l,j) + conjg(v2b)
           IF ( i == l ) fock_mtx0(k,j) = fock_mtx0(k,j) - conjg(v2b)
           IF ( j == l ) fock_mtx0(k,i) = fock_mtx0(k,i) + conjg(v2b)
           IF ( j == k ) fock_mtx0(l,i) = fock_mtx0(l,i) - conjg(v2b)
        end IF
        ! Compute 2b vacuum expectation value
        IF ( bra == ket ) then
           e0_2 = e0_2 + v2b
        end IF
     end DO
     
     DO bra = bra_min, bra_max
        DO ket = bra_min, bra-1
           v2b_hhhh(ch)%cval(bra,ket) = conjg(v2b_hhhh(ch)%cval(ket,bra))
        end DO
     end DO
     DEALLOCATE( loop_inds, loop_vals )
  end DO
  
end SUBROUTINE v2b_hhhh_2b

SUBROUTINE v2b_hphh_2b(fock_mtx0)
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE chiral_potentials

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: fock_mtx0(all_orbit%total_orbits, all_orbit%total_orbits)
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max, bra,ket,ket0
  INTEGER :: ch_v2, v2bra,v2ket, size1, count, ii
  INTEGER :: i,a,j,k, phase
  COMPLEX(dpc) :: v2b
  INTEGER, ALLOCATABLE :: loop_inds(:,:)
  COMPLEX(dpc), ALLOCATABLE :: loop_vals(:)
  
  ! hphh
  DO ch = 1, channels_2b%number_confs
     IF ( check_my_channel_v2b_hphh(ch) == 0 ) cycle
     bra_min = mapping_v2b_hphh(iam+1,ch,2)
     bra_max = mapping_v2b_hphh(iam+1,ch,3)
     bra_confs = number_2b(2,ch)
     ket_confs = number_2b(1,ch)

     count = 0
     DO bra = bra_min, bra_max
        DO ket = 1, ket_confs
           count = count + 1
        end DO
     end DO
     ALLOCATE( loop_inds(count, 2) )
     ALLOCATE( loop_vals(count) )
     loop_inds = 0
     loop_vals = 0.d0
     count = 0
     DO bra = bra_min, bra_max
        DO ket = 1, ket_confs
           count = count + 1
           loop_inds(count,1) = bra
           loop_inds(count,2) = ket
        end DO
     end DO

     !$omp parallel default(shared) private(ii, bra,ket,i,a,j,k, v2b)
     !$omp do schedule(static)
     DO ii = 1, count
        bra = loop_inds(ii, 1)
        ket = loop_inds(ii, 2)
        i   = lookup_2b_configs(2,ch)%ival2(1,bra)
        a   = lookup_2b_configs(2,ch)%ival2(2,bra)
        j   = lookup_2b_configs(1,ch)%ival2(1,ket)
        k   = lookup_2b_configs(1,ch)%ival2(2,ket)
        loop_vals(ii) = v2int(i,a,j,k)
     end DO
     !$omp end do
     !$omp end parallel
     
     DO ii = 1, count
        bra = loop_inds(ii, 1)
        ket = loop_inds(ii, 2)
        i   = lookup_2b_configs(2,ch)%ival2(1,bra)
        a   = lookup_2b_configs(2,ch)%ival2(2,bra)
        j   = lookup_2b_configs(1,ch)%ival2(1,ket)
        k   = lookup_2b_configs(1,ch)%ival2(2,ket)
        v2b = loop_vals(ii)
        v2b_hphh(ch)%cval(bra,ket) = v2b_hphh(ch)%cval(bra,ket) + v2b
        ! Compute fock_mtx
        IF ( i == j ) then
           fock_mtx0(k,a) = fock_mtx0(k,a) + v2b
           fock_mtx0(a,k) = fock_mtx0(a,k) + conjg(v2b)
        else IF ( i == k ) then
           fock_mtx0(j,a) = fock_mtx0(j,a) - v2b
           fock_mtx0(a,j) = fock_mtx0(a,j) - conjg(v2b)
        end IF
     end DO
     
     DEALLOCATE( loop_inds, loop_vals )
  end DO

end SUBROUTINE v2b_hphh_2b

SUBROUTINE v2b_pphp_2b
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE chiral_potentials

  IMPLICIT NONE
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max, bra,ket,ket0
  INTEGER :: ch_v2, v2bra,v2ket, size1
  INTEGER :: a,b,i,c, phase
  COMPLEX(dpc) :: v2b

  DO ch = 1, channels_2b%number_confs
     IF ( check_my_channel_v2b_pphp(ch) == 0 ) cycle
     bra_min = mapping_v2b_pphp(iam+1,ch,2)
     bra_max = mapping_v2b_pphp(iam+1,ch,3)
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(2,ch)
     !$omp parallel default(shared) private(bra,ket,ket0, a,b,i,c, v2b)
     !$omp do schedule(static)
     DO bra = bra_min, bra_max
        a   = lookup_2b_configs(3,ch)%ival2(1,bra)
        b   = lookup_2b_configs(3,ch)%ival2(2,bra)
        DO ket = 1, ket_confs
           i   = lookup_2b_configs(2,ch)%ival2(1,ket)
           c   = lookup_2b_configs(2,ch)%ival2(2,ket)
           v2b = v2int(a,b,i,c)
           v2b_pphp(ch)%cval(bra,ket) = v2b_pphp(ch)%cval(bra,ket) + v2b
        end DO
     end DO
     !$omp end do
     !$omp end parallel
  end DO
end SUBROUTINE v2b_pphp_2b

SUBROUTINE v2b_pphh_2b
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE chiral_potentials

  IMPLICIT NONE
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max, bra,ket,ket0
  INTEGER :: ch_v2, v2bra,v2ket, size1
  INTEGER :: a,b,i,j, phase
  COMPLEX(dpc) :: v2b

  DO ch = 1, channels_2b%number_confs
     IF ( check_my_channel_v2b_pphh(ch) == 0 ) cycle
     bra_min = mapping_v2b_pphh(iam+1,ch,2)
     bra_max = mapping_v2b_pphh(iam+1,ch,3)
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)
     !$omp parallel default(shared) private(bra,ket,ket0, a,b,i,j, v2b)
     !$omp do schedule(static)
     DO bra = bra_min, bra_max
        a   = lookup_2b_configs(3,ch)%ival2(1,bra)
        b   = lookup_2b_configs(3,ch)%ival2(2,bra)
        DO ket = 1, ket_confs
           i   = lookup_2b_configs(1,ch)%ival2(1,ket)
           j   = lookup_2b_configs(1,ch)%ival2(2,ket)
           v2b = v2int(a,b,i,j)
           v2b_pphh(ch)%cval(bra,ket) = v2b_pphh(ch)%cval(bra,ket) + v2b
        end DO
     end DO
     !$omp end do
     !$omp end parallel
  end DO
end SUBROUTINE v2b_pphh_2b

SUBROUTINE v2b_pppp_2b
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE chiral_potentials

  IMPLICIT NONE
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max, bra,ket,ket0
  INTEGER :: ch_v2, v2bra,v2ket, size1
  INTEGER :: a,b,c,d, phase
  COMPLEX(dpc) :: v2b

  DO ch = 1, channels_2b%number_confs
     IF ( check_my_channel_v2b_pppp(ch) == 0 ) cycle
     bra_min = mapping_v2b_pppp(iam+1,ch,2)
     bra_max = mapping_v2b_pppp(iam+1,ch,3)
     bra_confs = number_2b(3,ch)
     !$omp parallel default(shared) private(bra,ket,ket0, size1, a,b,c,d, v2b)
     !$omp do schedule(static)
     DO bra = bra_min, bra_max
        a   = lookup_2b_configs(3,ch)%ival2(1,bra)
        b   = lookup_2b_configs(3,ch)%ival2(2,bra)
        size1 = bra + bra_confs-1
        DO ket0 = bra, size1
           ket = mod(ket0-1,bra_confs) + 1
           IF ( bra_min <= ket .and. ket < bra ) cycle
           c   = lookup_2b_configs(3,ch)%ival2(1,ket)
           d   = lookup_2b_configs(3,ch)%ival2(2,ket)
           v2b = v2int(a,b,c,d)
           v2b_pppp(ch)%cval(bra,ket) = v2b_pppp(ch)%cval(bra,ket) + v2b
        end DO
     end DO
     !$omp end do
     !$omp end parallel
     DO bra = bra_min, bra_max
        DO ket = bra_min, bra-1
           v2b_pppp(ch)%cval(bra,ket) = conjg(v2b_pppp(ch)%cval(ket,bra))
        end DO
     end DO
  end DO
end SUBROUTINE v2b_pppp_2b


SUBROUTINE v2b_hphp_3b(fock_mtx0)
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE chiral_potentials

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: fock_mtx0(all_orbit%total_orbits, all_orbit%total_orbits)
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max, bra,ket,ket0
  INTEGER :: ch_v2, v2bra,v2ket, size1, count, ii
  INTEGER :: i,a,j,b, h, phase
  COMPLEX(dpc) :: v3b, vtmp
  INTEGER, ALLOCATABLE :: loop_inds(:,:)
  COMPLEX(dpc), ALLOCATABLE :: loop_vals(:)
  
  DO ch = 1, channels_2b%number_confs
     IF ( check_my_channel_v2b_hphp(ch) == 0 ) cycle
     bra_min = mapping_v2b_hphp(iam+1,ch,2)
     bra_max = mapping_v2b_hphp(iam+1,ch,3)
     bra_confs = number_2bcross(2,ch)
     
     count = 0
     DO bra = bra_min, bra_max
        size1 = bra + bra_confs-1
        DO ket0 = bra, size1
           ket = mod(ket0-1,bra_confs) + 1
           IF ( bra_min <= ket .and. ket < bra ) cycle
           count = count + 1
        end DO
     end DO
     ALLOCATE( loop_inds(count, 2) )
     ALLOCATE( loop_vals(count) )
     loop_inds = 0
     loop_vals = 0.d0
     count = 0
     DO bra = bra_min, bra_max
        size1 = bra + bra_confs-1
        DO ket0 = bra, size1
           ket = mod(ket0-1,bra_confs) + 1
           IF ( bra_min <= ket .and. ket < bra ) cycle
           count = count + 1
           loop_inds(count,1) = bra
           loop_inds(count,2) = ket
        end DO
     end DO
     
     !$omp parallel default(shared) private(ii, bra,ket,i,a,j,b,h, v3b,vtmp)
     !$omp do schedule(static)
     DO ii = 1, count
        bra = loop_inds(ii, 1)
        ket = loop_inds(ii, 2)
        i   = lookup_2b_configs(2,ch)%ival2(1,bra)
        a   = lookup_2b_configs(2,ch)%ival2(2,bra)
        j   = lookup_2b_configs(2,ch)%ival2(1,ket)
        b   = lookup_2b_configs(2,ch)%ival2(2,ket)
        vtmp = 0.d0
        DO h = 1, below_ef
           IF ( h == i .or. h == j ) cycle
           v3b = v3int(i,a,h,j,b,h)
           vtmp = vtmp + v3b
        end DO
        loop_vals(ii) = vtmp
     end DO
     !$omp end do
     !$omp end parallel
     
     DO ii = 1, count
        bra = loop_inds(ii, 1)
        ket = loop_inds(ii, 2)
        i   = lookup_2b_configs(2,ch)%ival2(1,bra)
        a   = lookup_2b_configs(2,ch)%ival2(2,bra)
        j   = lookup_2b_configs(2,ch)%ival2(1,ket)
        b   = lookup_2b_configs(2,ch)%ival2(2,ket)
        v3b = loop_vals(ii)
        v2b_hphp(ch)%cval(bra,ket) = v2b_hphp(ch)%cval(bra,ket) + v3b
        ! Compute fock_mtx
        IF ( i == j ) then
           fock_mtx0(a,b) = fock_mtx0(a,b) + v3b/2.d0
           IF ( bra < ket .and. ket <= bra_max ) then
              fock_mtx0(b,a) = fock_mtx0(b,a) + conjg(v3b)/2.d0
           end IF
        end IF
     end DO
     
     DO bra = bra_min, bra_max
        DO ket = bra_min, bra-1
           v2b_hphp(ch)%cval(bra,ket) = conjg(v2b_hphp(ch)%cval(ket,bra))
        end DO
     end DO
     DEALLOCATE( loop_inds, loop_vals )
  end DO

end SUBROUTINE v2b_hphp_3b

SUBROUTINE v2b_hhhh_3b(fock_mtx0, e0_3)
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE chiral_potentials

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: e0_3
  COMPLEX(dpc), INTENT(INOUT) :: fock_mtx0(all_orbit%total_orbits, all_orbit%total_orbits)
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max, bra,ket,ket0
  INTEGER :: ch_v2, v2bra,v2ket, size1, count, ii
  INTEGER :: i,j,k,l, h, phase
  COMPLEX(dpc) :: v3b, vtmp
  INTEGER, ALLOCATABLE :: loop_inds(:,:)
  COMPLEX(dpc), ALLOCATABLE :: loop_vals(:)
  
  DO ch = 1, channels_2b%number_confs
     IF ( check_my_channel_v2b_hhhh(ch) == 0 ) cycle
     bra_min = mapping_v2b_hhhh(iam+1,ch,2)
     bra_max = mapping_v2b_hhhh(iam+1,ch,3)
     bra_confs = number_2b(1,ch)

     count = 0
     DO bra = bra_min, bra_max
        size1 = bra + bra_confs-1
        DO ket0 = bra, size1
           ket = mod(ket0-1,bra_confs) + 1
           IF ( bra_min <= ket .and. ket < bra ) cycle
           count = count + 1
        end DO
     end DO
     ALLOCATE( loop_inds(count, 2) )
     ALLOCATE( loop_vals(count) )
     loop_inds = 0
     loop_vals = 0.d0
     count = 0
     DO bra = bra_min, bra_max
        size1 = bra + bra_confs-1
        DO ket0 = bra, size1
           ket = mod(ket0-1,bra_confs) + 1
           IF ( bra_min <= ket .and. ket < bra ) cycle
           count = count + 1
           loop_inds(count,1) = bra
           loop_inds(count,2) = ket
        end DO
     end DO
     
     !$omp parallel default(shared) private(ii, bra,ket,i,j,k,l,h, v3b,vtmp)
     !$omp do schedule(static)
     DO ii = 1, count
        bra = loop_inds(ii, 1)
        ket = loop_inds(ii, 2)
        i   = lookup_2b_configs(1,ch)%ival2(1,bra)
        j   = lookup_2b_configs(1,ch)%ival2(2,bra)
        k   = lookup_2b_configs(1,ch)%ival2(1,ket)
        l   = lookup_2b_configs(1,ch)%ival2(2,ket)
        vtmp = 0.d0
        DO h = 1, below_ef
           IF ( h == i .or. h == j .or. h == k .or. h == l ) cycle
           v3b = v3int(i,j,h,k,l,h)
           vtmp = vtmp + v3b
        end DO
        loop_vals(ii) = vtmp
     end DO
     !$omp end do
     !$omp end parallel
     
     DO ii = 1, count
        bra = loop_inds(ii, 1)
        ket = loop_inds(ii, 2)
        i   = lookup_2b_configs(1,ch)%ival2(1,bra)
        j   = lookup_2b_configs(1,ch)%ival2(2,bra)
        k   = lookup_2b_configs(1,ch)%ival2(1,ket)
        l   = lookup_2b_configs(1,ch)%ival2(2,ket)
        v3b = loop_vals(ii)
        v2b_hhhh(ch)%cval(bra,ket) = v2b_hhhh(ch)%cval(bra,ket) + v3b
        ! Compute fock_mtx
        IF ( i == k ) fock_mtx0(j,l) = fock_mtx0(j,l) + v3b/2.d0
        IF ( i == l ) fock_mtx0(j,k) = fock_mtx0(j,k) - v3b/2.d0
        IF ( j == l ) fock_mtx0(i,k) = fock_mtx0(i,k) + v3b/2.d0
        IF ( j == k ) fock_mtx0(i,l) = fock_mtx0(i,l) - v3b/2.d0
        IF ( bra < ket .and. ket <= bra_max ) then
           IF ( i == k ) fock_mtx0(l,j) = fock_mtx0(l,j) + conjg(v3b)/2.d0
           IF ( i == l ) fock_mtx0(k,j) = fock_mtx0(k,j) - conjg(v3b)/2.d0
           IF ( j == l ) fock_mtx0(k,i) = fock_mtx0(k,i) + conjg(v3b)/2.d0
           IF ( j == k ) fock_mtx0(l,i) = fock_mtx0(l,i) - conjg(v3b)/2.d0
        end IF
        ! Compute 3b vacuum expectation value
        IF ( bra == ket ) then
           e0_3 = e0_3 + v3b/3.d0
        end IF
     end DO
     
     DO bra = bra_min, bra_max
        DO ket = bra_min, bra-1
           v2b_hhhh(ch)%cval(bra,ket) = conjg(v2b_hhhh(ch)%cval(ket,bra))
        end DO
     end DO
     DEALLOCATE( loop_inds, loop_vals )
  end DO

end SUBROUTINE v2b_hhhh_3b

SUBROUTINE v2b_hphh_3b(fock_mtx0)
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE chiral_potentials

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: fock_mtx0(all_orbit%total_orbits, all_orbit%total_orbits)
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max, bra,ket,ket0
  INTEGER :: ch_v2, v2bra,v2ket, size1, count, ii
  INTEGER :: i,a,j,k, h, phase
  COMPLEX(dpc) :: v3b, vtmp
  INTEGER, ALLOCATABLE :: loop_inds(:,:)
  COMPLEX(dpc), ALLOCATABLE :: loop_vals(:)
  
  ! hphh
  DO ch = 1, channels_2b%number_confs
     IF ( check_my_channel_v2b_hphh(ch) == 0 ) cycle
     bra_min = mapping_v2b_hphh(iam+1,ch,2)
     bra_max = mapping_v2b_hphh(iam+1,ch,3)
     bra_confs = number_2b(2,ch)
     ket_confs = number_2b(1,ch)

     count = 0
     DO bra = bra_min, bra_max
        DO ket = 1, ket_confs
           count = count + 1
        end DO
     end DO
     ALLOCATE( loop_inds(count, 2) )
     ALLOCATE( loop_vals(count) )
     loop_inds = 0
     loop_vals = 0.d0
     count = 0
     DO bra = bra_min, bra_max
        DO ket = 1, ket_confs
           count = count + 1
           loop_inds(count,1) = bra
           loop_inds(count,2) = ket
        end DO
     end DO

     !$omp parallel default(shared) private(ii, bra,ket,i,a,j,k,h, v3b,vtmp)
     !$omp do schedule(static)
     DO ii = 1, count
        bra = loop_inds(ii, 1)
        ket = loop_inds(ii, 2)
        i   = lookup_2b_configs(2,ch)%ival2(1,bra)
        a   = lookup_2b_configs(2,ch)%ival2(2,bra)
        j   = lookup_2b_configs(1,ch)%ival2(1,ket)
        k   = lookup_2b_configs(1,ch)%ival2(2,ket)
        vtmp = 0.d0
        DO h = 1, below_ef
           IF ( h == i .or. h == j .or. j == k ) cycle
           v3b = v3int(i,a,h,j,k,h)
           vtmp = vtmp + v3b
        end DO
        loop_vals(ii) = vtmp
     end DO
     !$omp end do
     !$omp end parallel
     
     DO ii = 1, count
        bra = loop_inds(ii, 1)
        ket = loop_inds(ii, 2)
        i   = lookup_2b_configs(2,ch)%ival2(1,bra)
        a   = lookup_2b_configs(2,ch)%ival2(2,bra)
        j   = lookup_2b_configs(1,ch)%ival2(1,ket)
        k   = lookup_2b_configs(1,ch)%ival2(2,ket)
        v3b = loop_vals(ii)
        v2b_hphh(ch)%cval(bra,ket) = v2b_hphh(ch)%cval(bra,ket) + v3b
        ! Compute fock_mtx
        IF ( i == j ) then
           fock_mtx0(k,a) = fock_mtx0(k,a) + v3b/2.d0
           fock_mtx0(a,k) = fock_mtx0(a,k) + conjg(v3b)/2.d0
        ELSE IF ( i == k ) then
           fock_mtx0(j,a) = fock_mtx0(j,a) - v3b/2.d0
           fock_mtx0(a,j) = fock_mtx0(a,j) - conjg(v3b)/2.d0
        end IF
     end DO
     
     DEALLOCATE( loop_inds, loop_vals )
  end DO
  
end SUBROUTINE v2b_hphh_3b

SUBROUTINE v2b_pphp_3b
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE chiral_potentials

  IMPLICIT NONE
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max, bra,ket,ket0
  INTEGER :: ch_v2, v2bra,v2ket, size1
  INTEGER :: a,b,i,c, h, phase
  COMPLEX(dpc) :: v3b, vtmp

  DO ch = 1, channels_2b%number_confs
     IF ( check_my_channel_v2b_pphp(ch) == 0 ) cycle
     bra_min = mapping_v2b_pphp(iam+1,ch,2)
     bra_max = mapping_v2b_pphp(iam+1,ch,3)
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(2,ch)
     !$omp parallel default(shared) private(bra,ket,ket0, a,b,i,c,h, v3b,vtmp)
     !$omp do schedule(static)
     DO bra = bra_min, bra_max
        a   = lookup_2b_configs(3,ch)%ival2(1,bra)
        b   = lookup_2b_configs(3,ch)%ival2(2,bra)
        DO ket = 1, ket_confs
           i   = lookup_2b_configs(2,ch)%ival2(1,ket)
           c   = lookup_2b_configs(2,ch)%ival2(2,ket)
           vtmp = 0.d0
           DO h = 1, below_ef
              IF ( h == i ) cycle
              v3b = v3int(a,b,h,i,c,h)
              vtmp = vtmp + v3b
           end DO
           v2b_pphp(ch)%cval(bra,ket) = v2b_pphp(ch)%cval(bra,ket) + vtmp
        end DO
     end DO
     !$omp end do
     !$omp end parallel
  end DO
end SUBROUTINE v2b_pphp_3b

SUBROUTINE v2b_pphh_3b
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE chiral_potentials

  IMPLICIT NONE
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max, bra,ket,ket0
  INTEGER :: ch_v2, v2bra,v2ket, size1
  INTEGER :: a,b,i,j, h, phase
  COMPLEX(dpc) :: v3b, vtmp

  DO ch = 1, channels_2b%number_confs
     IF ( check_my_channel_v2b_pphh(ch) == 0 ) cycle
     bra_min = mapping_v2b_pphh(iam+1,ch,2)
     bra_max = mapping_v2b_pphh(iam+1,ch,3)
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)
     !$omp parallel default(shared) private(bra,ket,ket0, a,b,i,j,h, v3b,vtmp)
     !$omp do schedule(static)
     DO bra = bra_min, bra_max
        a   = lookup_2b_configs(3,ch)%ival2(1,bra)
        b   = lookup_2b_configs(3,ch)%ival2(2,bra)
        DO ket = 1, ket_confs
           i   = lookup_2b_configs(1,ch)%ival2(1,ket)
           j   = lookup_2b_configs(1,ch)%ival2(2,ket)
           vtmp = 0.d0
           DO h = 1, below_ef
              IF ( h == i .or. h == j ) cycle
              v3b = v3int(a,b,h,i,j,h)
              vtmp = vtmp + v3b
           end DO
           v2b_pphh(ch)%cval(bra,ket) = v2b_pphh(ch)%cval(bra,ket) + vtmp
        end DO
     end DO
     !$omp end do
     !$omp end parallel
  end DO
end SUBROUTINE v2b_pphh_3b

SUBROUTINE v2b_pppp_3b
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE chiral_potentials

  IMPLICIT NONE
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max, bra,ket,ket0
  INTEGER :: ch_v2, v2bra,v2ket, size1
  INTEGER :: a,b,c,d, h, phase
  COMPLEX(dpc) :: v2b, v3b, vtmp

  DO ch = 1, channels_2b%number_confs
     IF ( check_my_channel_v2b_pppp(ch) == 0 ) cycle
     bra_min = mapping_v2b_pppp(iam+1,ch,2)
     bra_max = mapping_v2b_pppp(iam+1,ch,3)
     bra_confs = number_2b(3,ch)
     !$omp parallel default(shared) private(bra,ket,ket0, size1, a,b,c,d,h, v3b,vtmp)
     !$omp do schedule(static)
     DO bra = bra_min, bra_max
        a   = lookup_2b_configs(3,ch)%ival2(1,bra)
        b   = lookup_2b_configs(3,ch)%ival2(2,bra)
        size1 = bra + bra_confs-1
        DO ket0 = bra, size1
           ket = mod(ket0-1,bra_confs) + 1
           IF ( bra_min <= ket .and. ket < bra ) cycle
           c   = lookup_2b_configs(3,ch)%ival2(1,ket)
           d   = lookup_2b_configs(3,ch)%ival2(2,ket)
           vtmp = 0.d0
           DO h = 1, below_ef
              v3b = v3int(a,b,h,c,d,h)
              vtmp = vtmp + v3b
           end DO
           v2b_pppp(ch)%cval(bra,ket) = v2b_pppp(ch)%cval(bra,ket) + vtmp
        end DO
     end DO
     !$omp end do
     !$omp end parallel
     DO bra = bra_min, bra_max
        DO ket = bra_min, bra-1
           v2b_pppp(ch)%cval(bra,ket) = conjg(v2b_pppp(ch)%cval(ket,bra))
        end DO
     end DO
  end DO
end SUBROUTINE v2b_pppp_3b


SUBROUTINE precalc_chp_functions
  USE PARALLEL
  USE single_particle_orbits
  USE constants
  USE chiral_potentials
  USE chiral_tables 
  USE chiral_constants
  USE mem_tracker

  implicit none 
  REAL(dp) :: k1(3),k2(3), qtrans(3)
  INTEGER :: p,q, m1,m2,m3,m4,m5,m6, t1,t2
  INTEGER :: nx1,ny1,nz1,nx2,ny2,nz2
  INTEGER :: nx3,ny3,nz3, nxmin,nxmax  
    
  IF ( .not. allocated(sigma_dot_q_tab) ) allocate( sigma_dot_q_tab(-1:1, -1:1, NNmin:NNmax, NNmin:NNmax, NNmin:NNmax) )
  sigma_dot_q_tab = 0.d0
  CALL mem_register('v2', REAL(9 * (NNmax-NNmin+1)**3 * 16.d0, dp))
  
  IF ( .not. allocated(sigma_dot_q_ope_tab) ) allocate( sigma_dot_q_ope_tab(-1:1, -1:1, NNmin:NNmax, NNmin:NNmax, NNmin:NNmax) )
  sigma_dot_q_ope_tab = 0.d0
  CALL mem_register('v2', REAL(9 * (NNmax-NNmin+1)**3 * 16.d0, dp))
  
  nxmin =  1 
  nxmax = -1
  
  DO p = 1, all_orbit%total_orbits
     
     m1 = all_orbit%sz(p) 
     t1 = all_orbit%tz(p) 
     k1(1) = all_orbit%kx(p)
     k1(2) = all_orbit%ky(p)
     k1(3) = all_orbit%kz(p)
     
     nx1 = all_orbit%nx(p)
     ny1 = all_orbit%ny(p)
     nz1 = all_orbit%nz(p)
     
     DO q = 1, all_orbit%total_orbits

        m2 = all_orbit%sz(q) 
        t2 = all_orbit%tz(q) 
        k2(1) = all_orbit%kx(q)
        k2(2) = all_orbit%ky(q)
        k2(3) = all_orbit%kz(q)
        
        qtrans = hbarc * (k2-k1)
        
        nx2 = all_orbit%nx(q)
        ny2 = all_orbit%ny(q)
        nz2 = all_orbit%nz(q)
        
        sigma_dot_q_tab(m1,m2, nx2-nx1, ny2-ny1, nz2-nz1) = chp_sigma_dot_q_mtx(m1,m2,qtrans)
        sigma_dot_q_ope_tab(m1,m2, nx2-nx1, ny2-ny1, nz2-nz1) = chp_sigma1_dot_q_ope(m1,m2,qtrans)
                
        nx3 = ny1*nz2-nz1*ny2
        IF ( nx3 > nxmax ) nxmax = nx3
        IF ( nx3 < nxmin ) nxmin = nx3
        
     end DO
  end DO
  
  IF ( .not. allocated(sigma_dot_qxq_tab) ) allocate( sigma_dot_qxq_tab(-1:1, -1:1, nxmin:nxmax, nxmin:nxmax, nxmin:nxmax) )
  sigma_dot_qxq_tab = 0.d0 
  CALL mem_register('v2', REAL(9 * (nxmax-nxmin+1)**3 * 16.d0, dp))
  
  DO p = 1, all_orbit%total_orbits
     
     m1 = all_orbit%sz(p)
     t1 = all_orbit%tz(p)
     k1(1) = all_orbit%kx(p)
     k1(2) = all_orbit%ky(p)
     k1(3) = all_orbit%kz(p)
     
     nx1 = all_orbit%nx(p)
     ny1 = all_orbit%ny(p)
     nz1 = all_orbit%nz(p)
     
     DO q = 1, all_orbit%total_orbits

        m2 = all_orbit%sz(q)
        t2 = all_orbit%tz(q)
        k2(1) = all_orbit%kx(q)
        k2(2) = all_orbit%ky(q)
        k2(3) = all_orbit%kz(q)
        
        qtrans(1) = k1(2)*k2(3)-k1(3)*k2(2)
        qtrans(2) = k1(3)*k2(1)-k1(1)*k2(3)
        qtrans(3) = k1(1)*k2(2)-k1(2)*k2(1)
        
        nx2 = all_orbit%nx(q)
        ny2 = all_orbit%ny(q)
        nz2 = all_orbit%nz(q)
        
        nx3 = ny1*nz2-nz1*ny2
        ny3 = nz1*nx2-nx1*nz2
        nz3 = nx1*ny2-ny1*nx2
        
        sigma_dot_qxq_tab(m1,m2,nx3, ny3, nz3) = chp_sigma_dot_q_mtx(m1,m2,qtrans)
        
     end DO
  end DO
  
  DO m1 = -1,1,2 
     DO m2 = -1,1,2 

        DO m3 = -1,1,2 
           DO m4 = -1,1,2 
              
              tau_dot_tau_tab(m1,m2,m3,m4) = chp_tau_dot_tau_mtx(m1,m2,m3,m4)
              DO m5 = -1,1,2 
                 DO m6 = -1,1,2 
                    
                    tau1_dot_tauXtau_tab(m1,m2,m3,m4,m5,m6) = chp_tau1_dot_tauXtau(m1,m2,m3,m4,m5,m6)
                    
                 end DO
              end DO
           end DO
        end DO
     end DO
  end DO
  sigma_dot_sigma_tab = tau_dot_tau_tab
  
  delta_tab = 0.d0
  DO m1 = -1,1,2 
     DO m2 = -1,1,2 
        IF ( m1 == m2 ) delta_tab(m1,m2) = 1.d0
     end DO
  end DO
  
end SUBROUTINE precalc_chp_functions

SUBROUTINE deallocate_chp_functions 
  USE single_particle_orbits
  USE constants
  USE chiral_potentials
  USE chiral_tables 
  USE chiral_constants

  deallocate( sigma_dot_q_tab )
  deallocate( sigma_dot_q_ope_tab )
  deallocate( sigma_dot_qxq_tab )

end SUBROUTINE deallocate_chp_functions


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
