
! Linear mixing
SUBROUTINE linear_t2
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  
  IMPLICIT NONE
  INTEGER :: ch
  REAL(dp) :: mix1, mix2

  mix1 = cc_scale
  mix2 = 1.0_dp - cc_scale
  
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     t2_ccm(ch)%cval = mix1*t2_ccm(ch)%cval + mix2*t2_ccm_eqn(ch)%cval
  end DO
  
  IF ( test > 1 ) CALL build_tamp_test
  
end SUBROUTINE linear_t2


! DIIS extrapolation (Pulay mixing)
SUBROUTINE diis_t2(count)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE diis_mod
  
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: count 
  INTEGER :: ch, bra, ket, bra_confs, ket_confs
  INTEGER :: i, j, k, i1
  INTEGER :: lwork, info
  REAL(dp) :: mix1, mix2
  COMPLEX(dpc) :: t2
  INTEGER, ALLOCATABLE :: ipiv(:)
  COMPLEX(dpc), ALLOCATABLE :: work(:)
  
  mix1 = cc_scale
  mix2 = 1.0_dp - cc_scale
  
  ! Shift history buffer if subspace is full
  IF ( nstep > diis_subspace ) THEN
     nstep = diis_subspace
     DO i = 1, diis_subspace-1 
        t2_diis(:,i) = t2_diis(:,i+1)
     end DO
  END IF
  
  ! Store mixed T2 amplitudes into DIIS history
  i1 = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)
     DO bra = 1, bra_confs
        DO ket = 1, ket_confs
           i1 = i1 + 1
           t2 = t2_ccm(ch)%cval(bra,ket)
           t2_ccm(ch)%cval(bra,ket) = mix1*t2 + mix2*t2_ccm_eqn(ch)%cval(bra,ket)
           t2_diis(i1,nstep) = t2_ccm(ch)%cval(bra,ket)
        end DO
     end DO
  end DO

  ! Perform DIIS extrapolation every diis_step iterations
  IF ( count == diis_step ) THEN
     count = 0

     ! Build error vectors: difference of successive amplitudes
     DO i = 1, nstep-1
        DO k = 1, n_diis
           t2_diisd(k,i) = t2_diis(k,i+1) - t2_diis(k,i)
        end DO
     end DO

     ! Build DIIS matrix: B(i,j) = <e_i|e_j>, with Lagrange constraint row/col
     ALLOCATE( diis_mat(nstep,nstep) )
     ALLOCATE( sol_vect(nstep) )
     diis_mat = 0.0_dp
     DO i = 1, nstep-1
        DO j = 1, nstep-1
           DO k = 1, n_diis
              diis_mat(i,j) = diis_mat(i,j) + t2_diisd(k,i)*t2_diisd(k,j)
           end DO
        end DO
     end DO
     DO i = 1, nstep-1
        diis_mat(i,nstep) = -1.0_dp
        diis_mat(nstep,i) = -1.0_dp
     end DO
     diis_mat(nstep,nstep) = 0.0_dp
     
     ! Solve for DIIS coefficients
     sol_vect = 0.0_dp
     sol_vect(nstep) = -1.0_dp
     ALLOCATE( ipiv(nstep) )
     lwork = MAX(1, 2*nstep)
     ALLOCATE( work(lwork) )
     CALL ZSYSV('U', nstep, 1, diis_mat, nstep, ipiv, sol_vect, nstep, work, lwork, info)

     IF ( iam == 0 ) write(6,'(A29)') ' ...DIIS extrapolation...    '
     
     ! Construct extrapolated amplitudes: T2 = sum_i c_i * T2_i
     ! Use t2_diisd(:,1) as workspace
     DO i1 = 1, n_diis
        t2_diisd(i1,1) = sol_vect(1)*t2_diis(i1,1)
     end DO
     DO i = 2, nstep-1
        DO i1 = 1, n_diis
           t2_diisd(i1,1) = t2_diisd(i1,1) + sol_vect(i)*t2_diis(i1,i)
        end DO
     end DO

     ! Unpack extrapolated amplitudes back into t2_ccm
     i1 = 0
     DO ch = 1, channels_2b%number_confs
        IF ( number_2b(3,ch) == 0 ) cycle
        IF ( number_2b(1,ch) == 0 ) cycle
        bra_confs = number_2b(3,ch)
        ket_confs = number_2b(1,ch)
        DO bra = 1, bra_confs
           DO ket = 1, ket_confs
              i1 = i1 + 1
              t2_ccm(ch)%cval(bra,ket) = t2_diisd(i1,1)                     
           end DO
        end DO
     end DO
     
     DEALLOCATE( diis_mat, sol_vect, ipiv, work )
  END IF

  IF ( test > 1 ) CALL build_tamp_test
  
end SUBROUTINE diis_t2


SUBROUTINE diis_setup
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE diis_mod
  USE mem_tracker

  IMPLICIT NONE
  INTEGER :: ch, i1, bra_confs, ket_confs

  i1 = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)
     i1 = i1 + bra_confs * ket_confs
  end DO
  
  n_diis  = i1
  nn_diis = diis_subspace
  ALLOCATE( t2_diis(n_diis, nn_diis+1) )
  ALLOCATE( t2_diisd(n_diis, nn_diis) )
  t2_diis  = 0.0_dp
  t2_diisd = 0.0_dp
  CALL mem_register('t2', REAL(2 * n_diis*(nn_diis+1.0_dp) * 16.0_dp, dp))

  IF ( iam == 0 ) write(6,'(A,I10,A,I4)') '  DIIS setup: n_diis =', n_diis, '  subspace =', nn_diis
  
end SUBROUTINE diis_setup


SUBROUTINE diis_take_down
  USE constants
  USE diis_mod
  USE mem_tracker

  IMPLICIT NONE  
  DEALLOCATE( t2_diis )
  DEALLOCATE( t2_diisd )
  CALL mem_release('t2', REAL(2 * n_diis*(nn_diis+1.0_dp) * 16.0_dp, dp))
end SUBROUTINE diis_take_down
