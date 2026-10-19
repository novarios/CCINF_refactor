
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
  USE mem_tracker

  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: count
  INTEGER :: ch, bra, ket, bra_confs, ket_confs
  INTEGER :: i, j, k, i1
  INTEGER :: lwork, info
  REAL(dp) :: mix1, mix2
  COMPLEX(dpc) :: t2
  INTEGER, ALLOCATABLE :: ipiv(:)
  COMPLEX(dpc), ALLOCATABLE :: work(:), buf(:)

  mix1 = cc_scale
  mix2 = 1.0_dp - cc_scale

  ! History of the last nstep steps: the mixed output y_i in t2_diis(:,i) and its residual
  ! r_i = y_i - (input T2 of step i) in t2_diisd(:,i). Drop step 1 if the subspace is full.
  ! Each rank stores only its slice diis_lo..diis_hi of the packed T2 index.
  IF ( nstep > diis_subspace ) THEN
     nstep = diis_subspace
     DO i = 1, diis_subspace-1
        t2_diis(:,i) = t2_diis(:,i+1)
        t2_diisd(:,i) = t2_diisd(:,i+1)
     end DO
  END IF

  ! Mix, then store y_nstep and r_nstep
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
           IF ( i1 < diis_lo .or. i1 > diis_hi ) cycle
           t2_diis(i1,nstep) = t2_ccm(ch)%cval(bra,ket)
           t2_diisd(i1,nstep) = t2_ccm(ch)%cval(bra,ket) - t2
        end DO
     end DO
  end DO

  ! Perform DIIS extrapolation every diis_step iterations, once there are 3 steps of history
  IF ( count >= diis_step .and. nstep >= 3 ) THEN
     count = 0

     ! Build DIIS matrix: B(i,j) = <r_i|r_j>, with Lagrange constraint row/col
     ALLOCATE( diis_mat(nstep+1,nstep+1) )
     ALLOCATE( sol_vect(nstep+1) )
     diis_mat = 0.0_dp
     DO i = 1, nstep
        DO j = 1, nstep
           DO k = diis_lo, diis_hi
              diis_mat(i,j) = diis_mat(i,j) + t2_diisd(k,i)*t2_diisd(k,j)
           end DO
        end DO
     end DO
     CALL mpi_allreduce(mpi_in_place, diis_mat, (nstep+1)**2, mpi_complex16, mpi_sum, mpi_comm_world, ierror)
     DO i = 1, nstep
        diis_mat(i,nstep+1) = -1.0_dp
        diis_mat(nstep+1,i) = -1.0_dp
     end DO
     diis_mat(nstep+1,nstep+1) = 0.0_dp

     ! Solve for DIIS coefficients; keep the mixed T2 if B is singular
     sol_vect = 0.0_dp
     sol_vect(nstep+1) = -1.0_dp
     ALLOCATE( ipiv(nstep+1) )
     lwork = MAX(1, 2*(nstep+1))
     ALLOCATE( work(lwork) )
     CALL ZSYSV('U', nstep+1, 1, diis_mat, nstep+1, ipiv, sol_vect, nstep+1, work, lwork, info)
     IF ( info /= 0 ) THEN
        IF ( iam == 0 ) write(6,'(A,I4)') ' ...DIIS skipped, ZSYSV info =', info
        DEALLOCATE( diis_mat, sol_vect, ipiv, work )
        IF ( test > 1 ) CALL build_tamp_test
        RETURN
     END IF

     IF ( iam == 0 ) write(6,'(A29)') ' ...DIIS extrapolation...    '

     ! Extrapolate T2 = sum_i c_i y_i. Each rank does its slice of a transient
     ! full-length buffer, then an allgatherv rebuilds the replicated t2_ccm.
     ALLOCATE( buf(n_diis) )
     CALL mem_register('t2', REAL(n_diis,dp) * 16.0_dp)
     DO k = diis_lo, diis_hi
        t2 = 0.0_dp
        DO j = 1, nstep
           t2 = t2 + sol_vect(j)*t2_diis(k,j)
        end DO
        buf(k) = t2
     end DO
     CALL mpi_allgatherv(mpi_in_place, 0, mpi_datatype_null, buf, diis_counts, diis_displs, &
          mpi_complex16, mpi_comm_world, ierror)
     i1 = 0
     DO ch = 1, channels_2b%number_confs
        IF ( number_2b(3,ch) == 0 ) cycle
        IF ( number_2b(1,ch) == 0 ) cycle
        bra_confs = number_2b(3,ch)
        ket_confs = number_2b(1,ch)
        DO bra = 1, bra_confs
           DO ket = 1, ket_confs
              i1 = i1 + 1
              t2_ccm(ch)%cval(bra,ket) = buf(i1)
           end DO
        end DO
     end DO

     DEALLOCATE( diis_mat, sol_vect, ipiv, work, buf )
     CALL mem_release('t2', REAL(n_diis,dp) * 16.0_dp)
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
  INTEGER :: ch, i1, bra_confs, ket_confs, r, nbase, nrem

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

  ! Split the packed T2 index into contiguous, near-equal slices, one per rank
  ALLOCATE( diis_counts(0:num_procs-1), diis_displs(0:num_procs-1) )
  nbase = n_diis / num_procs
  nrem  = MOD(n_diis, num_procs)
  DO r = 0, num_procs-1
     diis_counts(r) = nbase
     IF ( r < nrem ) diis_counts(r) = nbase + 1
     diis_displs(r) = r*nbase + MIN(r, nrem)
  end DO
  diis_lo = diis_displs(iam) + 1
  diis_hi = diis_displs(iam) + diis_counts(iam)

  ALLOCATE( t2_diis(diis_lo:diis_hi, nn_diis) )
  ALLOCATE( t2_diisd(diis_lo:diis_hi, nn_diis) )
  t2_diis  = 0.0_dp
  t2_diisd = 0.0_dp
  CALL mem_register('t2', REAL(diis_counts(iam),dp) * 2 * nn_diis * 16.0_dp)

  IF ( iam == 0 ) write(6,'(A,I10,A,I4,A,I10)') '  DIIS setup: n_diis =', n_diis, '  subspace =', nn_diis, &
       '  per rank <=', MAXVAL(diis_counts)
  
end SUBROUTINE diis_setup


SUBROUTINE diis_take_down
  USE parallel
  USE constants
  USE diis_mod
  USE mem_tracker

  IMPLICIT NONE
  CALL mem_release('t2', REAL(diis_counts(iam),dp) * 2 * nn_diis * 16.0_dp)
  DEALLOCATE( t2_diis )
  DEALLOCATE( t2_diisd )
  DEALLOCATE( diis_counts, diis_displs )
end SUBROUTINE diis_take_down
