
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

  ! History: x_1..x_nstep (mixed T2 per step) is kept as x_nstep in t2_diis(:,1) plus the
  ! error vectors e_i = x_(i+1) - x_i in t2_diisd(:,1:nstep-1). Drop e_1 if the subspace is full.
  ! Each rank stores only its slice diis_lo..diis_hi of the packed T2 index.
  IF ( nstep > diis_subspace ) THEN
     nstep = diis_subspace
     DO i = 1, diis_subspace-2
        t2_diisd(:,i) = t2_diisd(:,i+1)
     end DO
  END IF

  ! Mix, then store the new error vector and x_nstep
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
           IF ( nstep > 1 ) t2_diisd(i1,nstep-1) = t2_ccm(ch)%cval(bra,ket) - t2_diis(i1,1)
           t2_diis(i1,1) = t2_ccm(ch)%cval(bra,ket)
        end DO
     end DO
  end DO

  ! Perform DIIS extrapolation every diis_step iterations
  IF ( count == diis_step ) THEN
     count = 0

     ! Build DIIS matrix: B(i,j) = <e_i|e_j>, with Lagrange constraint row/col
     ALLOCATE( diis_mat(nstep,nstep) )
     ALLOCATE( sol_vect(nstep) )
     diis_mat = 0.0_dp
     DO i = 1, nstep-1
        DO j = 1, nstep-1
           DO k = diis_lo, diis_hi
              diis_mat(i,j) = diis_mat(i,j) + t2_diisd(k,i)*t2_diisd(k,j)
           end DO
        end DO
     end DO
     CALL mpi_allreduce(mpi_in_place, diis_mat, nstep*nstep, mpi_complex16, mpi_sum, mpi_comm_world, ierror)
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

     ! Extrapolate T2 = sum_(i<nstep) c_i x_i. With sum c_i = 1 this telescopes to
     ! x_nstep - sum_j C_j e_j, C_j = c_1 + ... + c_j. Each rank does its slice of a
     ! transient full-length buffer, then an allgatherv rebuilds the replicated t2_ccm.
     DO j = 2, nstep-1
        sol_vect(j) = sol_vect(j) + sol_vect(j-1)
     end DO
     ALLOCATE( buf(n_diis) )
     CALL mem_register('t2', REAL(n_diis,dp) * 16.0_dp)
     DO k = diis_lo, diis_hi
        t2 = t2_diis(k,1)
        DO j = 1, nstep-1
           t2 = t2 - sol_vect(j)*t2_diisd(k,j)
        end DO
        buf(k) = t2
     end DO
     ! Restart the history from the extrapolated T2, so the next error vector is its residual
     t2_diis(diis_lo:diis_hi,1) = buf(diis_lo:diis_hi)
     nstep = 1
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

  ALLOCATE( t2_diis(diis_lo:diis_hi, 1) )
  ALLOCATE( t2_diisd(diis_lo:diis_hi, nn_diis-1) )
  t2_diis  = 0.0_dp
  t2_diisd = 0.0_dp
  CALL mem_register('t2', REAL(diis_counts(iam),dp) * nn_diis * 16.0_dp)

  IF ( iam == 0 ) write(6,'(A,I10,A,I4,A,I10)') '  DIIS setup: n_diis =', n_diis, '  subspace =', nn_diis, &
       '  per rank <=', MAXVAL(diis_counts)
  
end SUBROUTINE diis_setup


SUBROUTINE diis_take_down
  USE parallel
  USE constants
  USE diis_mod
  USE mem_tracker

  IMPLICIT NONE
  CALL mem_release('t2', REAL(diis_counts(iam),dp) * nn_diis * 16.0_dp)
  DEALLOCATE( t2_diis )
  DEALLOCATE( t2_diisd )
  DEALLOCATE( diis_counts, diis_displs )
end SUBROUTINE diis_take_down
