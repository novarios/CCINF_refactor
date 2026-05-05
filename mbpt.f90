
SUBROUTINE mbpt
  USE parallel
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE diis_mod
  USE hdf5_wrapper
  USE configurations
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  REAL(dp) :: startwtime, endwtime
  COMPLEX(dpc) :: e2, e3_1, e3_2, e3_3, e3_1b, e2_3b

  e2 = 0.d0
  e2_3b = 0.d0
  e3_1 = 0.d0
  e3_2 = 0.d0
  e3_3 = 0.d0
  e3_1b = 0.d0
  
  ecorr2 = 0.d0
  ecorr3 = 0.d0
  
  ! Setup t-amplitudes
  IF ( iam == 0 ) write(6,*) '...Setting up T amplitudes...'
  CALL setup_t_amplitudes
  IF ( cc_approx > 0 .and. .not. pre_gs0 ) then
     CALL setup_t3_amplitudes(.false.)
  end IF
  IF ( .not. pre_gs0 ) CALL diis_setup
  IF ( iam == 0 ) write(6,*) '...Setting up T amplitudes done!'
  CALL mem_report('T amplitudes')

  ! Total Memory
  IF ( iam == 0 ) write(6,*) '...Solving T Amplitudes Ready!'
  CALL mem_report('total')

  
  ! Start Solving MBPT
  startwtime = MPI_WTIME()
  CALL mbpt2(e2)
  CALL mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  IF ( iam == 0 ) write(6,*) ' MBPT2 = ', real(e2)
  IF ( iam == 0 ) write(6,*) ' Total execution time MBPT2', endwtime - startwtime

  startwtime = MPI_WTIME()
  CALL mbpt3_hh(e3_1)
  CALL mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  IF ( iam == 0 ) write(6,*) ' MBPT3_hh = ', real(e3_1)
  IF ( iam == 0 ) write(6,*) ' Total execution time MBPT3_hh', endwtime - startwtime

  startwtime = MPI_WTIME()
  CALL mbpt3_pp(e3_2)
  CALL mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  IF ( iam == 0 ) write(6,*) ' MBPT3_pp = ', real(e3_2)
  IF ( iam == 0 ) write(6,*) ' Total execution time MBPT3_pp', endwtime - startwtime

  startwtime = MPI_WTIME()
  CALL mbpt3_ph(e3_3)
  CALL mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  IF ( iam == 0 ) write(6,*) ' MBPT3_ph = ', real(e3_3)
  IF ( iam == 0 ) write(6,*) ' Total execution time MBPT3_ph', endwtime - startwtime

  ! startwtime = MPI_WTIME()
  ! CALL mbpt3_1b(e3_1b)
  ! CALL mpi_barrier(mpi_comm_world,ierror)
  ! endwtime = MPI_WTIME()
  ! IF ( iam == 0 ) write(6,*) ' MBPT3_1b = ', real(e3_1b)
  ! IF ( iam == 0 ) write(6,*) ' Total execution time MBPT3_1b', endwtime - startwtime
  
  IF ( tnf_approx > 1 ) then
     startwtime = MPI_WTIME()
     CALL setup_t3_mbpt
     CALL mbpt2_3b(e2_3b)
     CALL mpi_barrier(mpi_comm_world,ierror)
     endwtime = MPI_WTIME()
     IF ( iam == 0 ) write(6,*) ' MBPT2_3b = ', real(e2_3b)
     IF ( iam == 0 ) write(6,*) 'Total execution time MBPT2_3b', endwtime - startwtime
     ecorr3 = e2_3b
  end IF
  
  CALL mpi_barrier(mpi_comm_world,ierror)
  ecorr2 = e2 + e3_1 + e3_2 + e3_3 + e3_1b 
  ecorr = ecorr2 + ecorr3
  eccdt = e0 + ecorr
  
  IF ( iam == 0 ) write(6,*) ' Energy = ', real(eccdt)
  IF ( iam == 0 ) write(6,'(A10,1x,3(f14.9,3x))') 'E/A =    ', real(e0)/below_ef, real(e0+e2+e2_3b)/below_ef, real(eccdt)/below_ef
  
END SUBROUTINE mbpt


!
!      ===========
!      /\       /\
!     /  \a    /  \b
!   --|--|-----|--|--
!    i\  /    j\  /
!      \/       \/
!      ===========
!
! (1/4).<ij|v|ab>.<ab|v|ij>/e_abij = (1/4).<ij|v|ab>.<ab|t0|ij>
!
SUBROUTINE mbpt2(energy)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: energy
  INTEGER :: ch, ket, bra_min,bra_max, ket_confs
  INTEGER :: dim1,dim2,dim3
  COMPLEX(dpc), allocatable :: VE2(:,:)

  energy = 0.d0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     ket_confs = number_2b(1,ch)

     IF ( check_my_channel_v2b_pphh(ch) == 0 ) cycle
     bra_min = mapping_v2b_pphh(iam+1,ch,2)
     bra_max = mapping_v2b_pphh(iam+1,ch,3)
     
     ! <ij|x|ij> = <ij|v|ab>.<ab|t0|ij>
     dim1 = ket_confs
     dim2 = ket_confs
     dim3 = bra_max - bra_min + 1
     ALLOCATE( VE2(ket_confs, ket_confs) )
     VE2 = 0.d0
     CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), conjg(v2b_pphh(ch)%cval(bra_min:bra_max,:)), dim3, &
          t2_ccm(ch)%cval(bra_min:bra_max,:), dim3, dcmplx(1.d0,0.d0), VE2, dim1 )
     
     ! E = <ij|x|ij>
     DO ket = 1, ket_confs
        energy = energy + VE2(ket,ket)
     end DO
     DEALLOCATE(VE2)
     
  end DO
  CALL mpi_allreduce(mpi_in_place,energy,1,mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'mbpt2', 'allreduce')
  IF ( test == 3 ) CALL mbpt2_test(energy)
    
END SUBROUTINE mbpt2


!      ===========
!      /\       /\
!    a/  \k   l/  \b
!   --|--|-----|--|--
!     |  |     |  |
!     |  =======  |
!     |  |     |  |
!   --|--|-----|--|--
!     \  /i   j\  /
!      \/       \/
!      ===========
!
! (1/8).<kl|v|ab>.<ij|v|kl>.<ab|v|ij>/(e_abij.e_abkl)
! (1/8).<kl|t0*|ab>.<ij|v|kl>.<ab|t0|ij>
! (1/8).(<kl|t0*|ab>.<ab|t0|ij>).<ij|v|kl>
!
SUBROUTINE mbpt3_hh(energy)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: energy
  INTEGER :: ch, ket, bra_min,bra_max, ket_confs
  INTEGER :: dim1,dim2,dim3
  COMPLEX(dpc), allocatable :: VE2(:,:)

  energy = 0.d0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     ket_confs = number_2b(1,ch)

     IF ( check_my_channel_v2b_pphh(ch) == 0 ) cycle
     bra_min = mapping_v2b_pphh(iam+1,ch,2)
     bra_max = mapping_v2b_pphh(iam+1,ch,3)

     ! <kl|x|ij> = <kl|t0*|ab>.<ab|t0|ij>
     dim1 = ket_confs
     dim2 = ket_confs
     dim3 = bra_max - bra_min + 1
     ALLOCATE( VE2(ket_confs, ket_confs) )
     VE2 = 0.d0
     CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), conjg(t2_ccm(ch)%cval(bra_min:bra_max,:)), dim3, &
          t2_ccm(ch)%cval(bra_min:bra_max,:), dim3, dcmplx(1.d0,0.d0), VE2, dim1 )

     ! E = <ij|v|kl>.<kl|x|ij>
     DO ket = 1, ket_confs
        energy = energy + sum(VE2(ket,:)*v2b_hhhh(ch)%cval(:,ket))
     end DO
     DEALLOCATE(VE2)
     
  end DO
  CALL mpi_allreduce(mpi_in_place,energy,1,mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'mbpt3_hh', 'allreduce')
  IF ( test == 3 ) CALL mbpt3_hh_test(energy)
    
END SUBROUTINE mbpt3_hh


!      ===========
!      /\       /\
!     /  \a   b/  \
!   --|--|-----|--|--
!     |  |     |  |
!     |  =======  |
!     |  |     |  |
!   --|--|-----|--|--
!    i\  /c   d\  /j
!      \/       \/
!      ===========
!
! (1/8).<ij|v|ab>.<ab|v|cd>.<cd|v|ij>/(e_abij.e_cdij)
! (1/8).<ij|t0*|ab>.<ab|v|cd>.<cd|t0|ij>
! (1/8).(<cd|t0|ij>.<ij|t0*|ab>).<ab|v|cd>
!
SUBROUTINE mbpt3_pp(energy)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: energy
  INTEGER :: ch, bra, ket, bra_min,bra_max, bra_confs, ket_confs
  INTEGER :: dim1,dim2,dim3
  COMPLEX(dpc), allocatable :: VE2(:,:)

  energy = 0.d0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)

     IF ( check_my_channel_v2b_pphh(ch) == 0 ) cycle
     bra_min = mapping_v2b_pphh(iam+1,ch,2)
     bra_max = mapping_v2b_pphh(iam+1,ch,3)

     ! <cd|x|ab> = <cd|t0|ij>.<ij|t0*|ab>
     dim1 = bra_confs
     dim2 = bra_max - bra_min + 1
     dim3 = ket_confs
     ALLOCATE( VE2(bra_confs, bra_min:bra_max) )
     VE2 = 0.d0
     CALL ZGEMM ( 'n', 't', dim1, dim2, dim3, dcmplx(1.d0,0.d0), t2_ccm(ch)%cval, dim1, &
          conjg(t2_ccm(ch)%cval(bra_min:bra_max,:)), dim2, dcmplx(1.d0,0.d0), VE2, dim1 )
     
     ! E = <cd|x|ab>.<ab|v|cd>
     DO bra = 1, bra_confs
        energy = energy + sum(VE2(bra,bra_min:bra_max)*v2b_pppp(ch)%cval(bra_min:bra_max,bra))
     end DO
     DEALLOCATE(VE2)
     
  end DO
  CALL mpi_allreduce(mpi_in_place,energy,1,mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'mbpt3_pp', 'allreduce')
  IF ( test == 3 ) CALL mbpt3_pp_test(energy)
  
END SUBROUTINE mbpt3_pp


!      ===========
!      /\       /\
!     /  \k   c/  \
!   --|--|-----|--|--
!     |  |     |  |
!     |  =======  |
!     |  |     |  |
!   --|--|-----|--|--
!    a\  /i   b\  /j
!      \/       \/
!      ===========
!
! -<kj|v|ac>.<ic|v|kb>.<ab|v|ij>/(e_ackj.e_abij)
!
! -<k-c|v*|a-j>.<i-b|v|k-c>.<a-j|v|i-b>/(e_ackj.e_abij)
! -<k-c|t0*|a-j>.<i-b|v|k-c>.<a-j|t0|i-b>
! -(<k-c|t0*|a-j>.<a-j|t0|i-b>).<i-b|v|k-c>
!
SUBROUTINE mbpt3_ph(energy)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: energy
  INTEGER :: ch, bra, ket, bra_min,bra_max, ket_confs
  INTEGER :: dim1,dim2,dim3
  COMPLEX(dpc), allocatable :: VE2(:,:)

  energy = 0.d0
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(3,ch) == 0 ) cycle
     IF ( number_2bcross(2,ch) == 0 ) cycle
     ket_confs = number_2bcross(2,ch)

     IF ( check_my_channel_v2b_phhp_cross(ch) == 0 ) cycle
     bra_min = mapping_v2b_phhp_cross(iam+1,ch,2)
     bra_max = mapping_v2b_phhp_cross(iam+1,ch,3)

     ! <k-c|x|i-b> = <k-c|t0*|a-j>.<a-j|t0|i-b>
     dim1 = ket_confs
     dim2 = ket_confs
     dim3 = bra_max - bra_min + 1
     ALLOCATE( VE2(ket_confs, ket_confs) )
     VE2 = 0.d0
     CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), conjg(t2_ccm_cross(ch)%cval(bra_min:bra_max,:)), dim3, &
          t2_ccm_cross(ch)%cval(bra_min:bra_max,:), dim3, dcmplx(1.d0,0.d0), VE2, dim1 )
     
     ! E = -<k-c|x|i-b>.<i-b|v|k-c>
     DO ket = 1, ket_confs
        energy = energy - sum(VE2(ket,:)*v2b_hphp_cross(ch)%cval(:,ket))
     end DO
     DEALLOCATE(VE2)
     
  end DO
  CALL mpi_allreduce(mpi_in_place,energy,1,mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'mbpt3_ph', 'allreduce')
  IF ( test == 3 ) CALL mbpt3_ph_test(energy)
  
END SUBROUTINE mbpt3_ph


!      ===========        ===========
!      /\       /\        /\       /\
!     /  \    c/  \      /  \k    /  \
!   --|--|-----|--|--  --|--|-----|--|--
!     |  |     |  |      |  |     |  |
!     |  |  ====  |      |  ====  |  |
!     |  |     |  |      |  |     |  |
!   --|--|-----|--|--  --|--|-----|--|--
!    a\  /i   b\  /j    b\  /j   a\  /i
!      \/       \/        \/       \/
!      ===========        ===========
!
! +(1/2).<ab|v|ij>.<ij|v|ac>.<c|v|b>/(e_abij.e_acij) = +(1/2).<ab|t0|ij>.<ij|t0*|ac>.<c|v|b>
! -(1/2).<ab|v|ij>.<ik|v|ab>.<j|v|k>/(e_abij.e_abik) = -(1/2).<ab|t0|ij>.<ik|t0*|ab>.<j|v|k>
!
SUBROUTINE mbpt3_1b(energy)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: energy
  INTEGER :: ch, bra, ket, bra_min,bra_max, bra_confs, ket_confs
  INTEGER :: a,b,i,j, k,c, bra0,ket0, phase
  COMPLEX(dpc) :: v2b1, v2b2

  energy = 0.d0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)

     IF ( check_my_channel_v2b_pphh(ch) == 0 ) cycle
     bra_min = mapping_v2b_pphh(iam+1,ch,2)
     bra_max = mapping_v2b_pphh(iam+1,ch,3)

     DO bra = bra_min, bra_max
        a = lookup_2b_configs(3,ch)%ival2(1,bra)
        b = lookup_2b_configs(3,ch)%ival2(2,bra)
        DO ket = 1, ket_confs
           i = lookup_2b_configs(1,ch)%ival2(1,ket)
           j = lookup_2b_configs(1,ch)%ival2(2,ket)
           ! +(1/2).<ab|t0|ij>.<ij|t0*|ac>.<c|v|b>
           v2b1 = t2_ccm(ch)%cval(bra,ket)
           DO c = below_ef+1, tot_orbs
              IF ( c == b ) cycle
              IF ( pp_channel_2b%ival2(a,c) /= ch ) cycle
              phase = 1
              bra0 = pp_config_2b%ival2(a,c)
              IF ( c < a ) phase = -phase
              v2b2 = t2_ccm(ch)%cval(bra0,ket)
              energy = energy + phase * conjg(v2b2)*v2b1*fock_mtx(c,b)
           end DO
           ! -(1/2).<ab|t0|ij>.<ik|t0*|ab>.<j|v|k>
           DO k = 1, below_ef
              IF ( k == j ) cycle
              IF ( hh_channel_2b%ival2(i,k) /= ch ) cycle
              phase = 1
              ket0 = hh_config_2b%ival2(i,k)
              IF ( k < i ) phase = -phase
              v2b2 = t2_ccm(ch)%cval(bra,ket0)
              energy = energy - phase * conjg(v2b2)*v2b1*fock_mtx(j,k)
           end DO
        end DO
     end DO
  end DO  
  CALL mpi_allreduce(mpi_in_place,energy,1,mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'mbpt3_1b', 'allreduce')
  IF ( test == 3 ) CALL mbpt3_1b_test(energy)
  
END SUBROUTINE mbpt3_1b


!
!      ====================
!      /\       /\       /\
!     /  \a    /  \b    /  \c
!   --|--|-----|--|-----|--|--
!    i\  /    j\  /    k\  /
!      \/       \/       \/
!      ====================
!
! (1/36).<ijk|w|abc>.<abc|w|ijk>/e_abcijk
!
SUBROUTINE mbpt2_3b(energy)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE chiral_potentials
  USE mpi_check
  
  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: energy
  INTEGER :: ch3, ch1,ch2, ket_confs
  INTEGER :: bra,ket, bra0,ket0, bra_min,bra_max
  INTEGER :: a,b,c,i,j,k, cind1,kind1
  COMPLEX(dpc) :: v3b, energy0, denom

  ! energy = 0.d0
  ! DO ch3 = ch3_min, ch3_max
  !    DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
  !       c       = clist_t3(ch3)%ival2(cind1,1)
  !       ch1     = clist_t3(ch3)%ival2(cind1,2)
  !       bra_min = mapping_t3(ch3)%ival2(cind1,1)
  !       bra_max = mapping_t3(ch3)%ival2(cind1,2)
  !       DO kind1 = 1, klimit_t3(ch3)
  !          k         = klist_t3(ch3)%ival2(kind1,1)
  !          ch2       = klist_t3(ch3)%ival2(kind1,2)
  !          ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
  !          IF ( ch3 == 30 .and. ch2 == 63 ) write(6,*) '3**** ', ket_confs
  !          ! IF ( ch3 == 30 .and. kind1 == 1 ) write(6,*) '3%%% ', k,ch2,':',ket_confs
  !          ! write(6,*) '$$ ', ch3,kind1, ':', ket_confs, size(hh_config_t3(ch3)%ival1(ch2)%ival1)
           
  !          energy0 = 0.d0
  !          !$omp parallel default(shared) private(bra,ket,bra0,ket0, a,b,i,j, v3b,denom) reduction(+:energy0)
  !          !$omp do schedule(static) collapse(2)
  !          DO bra = bra_min, bra_max
  !             DO ket = 1, ket_confs
  !                bra0 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
  !                a = lookup_2b_configs(3,ch1)%ival2(1,bra0)
  !                b = lookup_2b_configs(3,ch1)%ival2(2,bra0)
  !                ket0 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
  !                i = lookup_2b_configs(1,ch2)%ival2(1,ket0)
  !                j = lookup_2b_configs(1,ch2)%ival2(2,ket0)
  !                v3b = v3int(a,b,c,i,j,k)
  !                denom = fock_mtx(i,i) + fock_mtx(j,j) + fock_mtx(k,k) &
  !                     - fock_mtx(a,a) - fock_mtx(b,b) - fock_mtx(c,c)
  !                energy0 = energy0 + conjg(v3b)*v3b/denom
  !             end DO
  !          end DO
  !          !$omp end do
  !          !$omp end parallel
  !          energy = energy + energy0
           
  !       end DO
  !    end DO
  ! end DO
  energy = 0.d0
  DO ch3 = ch3_min, ch3_max
     bra_min = mapping_t3_red(ch3,1)
     bra_max = mapping_t3_red(ch3,2)
     ket_confs = number_3b_t3(ch3,2)
     
     energy0 = 0.d0
     !$omp parallel default(shared) private(bra,ket, a,b,c,i,j,k, v3b,denom) reduction(+:energy0)
     !$omp do schedule(static) collapse(2)
     DO bra = bra_min, bra_max
        DO ket = 1, ket_confs
           a = lookup_t3_configs(1,ch3)%ival2(1,bra)
           b = lookup_t3_configs(1,ch3)%ival2(2,bra)
           c = lookup_t3_configs(1,ch3)%ival2(3,bra)
           i = lookup_t3_configs(2,ch3)%ival2(1,ket)
           j = lookup_t3_configs(2,ch3)%ival2(2,ket)
           k = lookup_t3_configs(2,ch3)%ival2(3,ket)
           v3b = v3int(a,b,c,i,j,k)
           denom = fock_mtx(i,i) + fock_mtx(j,j) + fock_mtx(k,k) &
                - fock_mtx(a,a) - fock_mtx(b,b) - fock_mtx(c,c)
           energy0 = energy0 + conjg(v3b)*v3b/denom
        end DO
     end DO
     !$omp end do
     !$omp end parallel
     energy = energy + energy0
     
  end DO  
  CALL mpi_allreduce(mpi_in_place,energy,1,mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'mbpt2_3b', 'allreduce')
  IF ( test == 3 ) CALL mbpt2_3b_test(energy)
    
END SUBROUTINE mbpt2_3b
