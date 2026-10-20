
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


! Newton-Krylov helpers. Vectors are this rank's slice diis_lo..diis_hi of the packed
! T2 index, with the same packing and slicing as the DIIS history (diis_setup first).

! x <- slice of t2_ccm (mode 1), or of the residual t2_ccm_eqn - t2_ccm (mode 2)
SUBROUTINE t2_slice(x, mode)
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE diis_mod

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: mode
  COMPLEX(dpc), INTENT(OUT) :: x(diis_lo:diis_hi)
  INTEGER :: ch, bra, ket, i1

  i1 = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     DO bra = 1, number_2b(3,ch)
        DO ket = 1, number_2b(1,ch)
           i1 = i1 + 1
           IF ( i1 < diis_lo .or. i1 > diis_hi ) cycle
           IF ( mode == 1 ) then
              x(i1) = t2_ccm(ch)%cval(bra,ket)
           else
              x(i1) = t2_ccm_eqn(ch)%cval(bra,ket) - t2_ccm(ch)%cval(bra,ket)
           end IF
        end DO
     end DO
  end DO

end SUBROUTINE t2_slice


! t2_ccm <- the full vector whose slices are x on each rank; recouples t2_ccm_cross
SUBROUTINE t2_set(x)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE diis_mod
  USE mem_tracker

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(IN) :: x(diis_lo:diis_hi)
  INTEGER :: ch, bra, ket, i1
  COMPLEX(dpc), ALLOCATABLE :: buf(:)

  ALLOCATE( buf(n_diis) )
  CALL mem_register('t2', REAL(n_diis,dp) * 16.0_dp)
  buf(diis_lo:diis_hi) = x
  CALL mpi_allgatherv(mpi_in_place, 0, mpi_datatype_null, buf, diis_counts, diis_displs, &
       mpi_complex16, mpi_comm_world, ierror)
  i1 = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     DO bra = 1, number_2b(3,ch)
        DO ket = 1, number_2b(1,ch)
           i1 = i1 + 1
           t2_ccm(ch)%cval(bra,ket) = buf(i1)
        end DO
     end DO
  end DO
  DEALLOCATE( buf )
  CALL mem_release('t2', REAL(n_diis,dp) * 16.0_dp)
  CALL t2_cross_recouple

end SUBROUTINE t2_set


! d = <a|b> over the full vector (conjugate-linear in a)
SUBROUTINE slice_dot(a, b, d)
  USE parallel
  USE constants
  USE diis_mod

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(IN) :: a(diis_lo:diis_hi), b(diis_lo:diis_hi)
  COMPLEX(dpc), INTENT(OUT) :: d

  d = SUM( CONJG(a)*b )
  CALL mpi_allreduce(mpi_in_place, d, 1, mpi_complex16, mpi_sum, mpi_comm_world, ierror)

end SUBROUTINE slice_dot


! jv = J.v for F(T2) = g(T2) - T2, by a forward difference about x0 (F(x0) = f0, |x0| = xnorm,
! |v| = 1). The step is about the square root of the single-precision T3 noise. Leaves t2_ccm
! at x0 + eps.v; the caller restores it.
SUBROUTINE cc_jv(x0, f0, xnorm, v, jv)
  USE constants
  USE diis_mod

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(IN) :: x0(diis_lo:diis_hi), f0(diis_lo:diis_hi), v(diis_lo:diis_hi)
  REAL(dp), INTENT(IN) :: xnorm
  COMPLEX(dpc), INTENT(OUT) :: jv(diis_lo:diis_hi)
  REAL(dp) :: eps, rnorm

  eps = 3.d-4 * (1.d0 + xnorm)
  jv = x0 + eps*v
  CALL t2_set(jv)
  CALL cc_residual(rnorm)
  CALL t2_slice(jv, 2)
  jv = (jv - f0) / eps

end SUBROUTINE cc_jv


! GMRES for the Newton step J.delta = -f0 at x0, stopped when the least-squares residual is
! below eta.|f0| or after m_max products. nmv returns the number of J.v products used; it is
! the Newton step number, for the printout. Leaves t2_ccm at the last probe point.
SUBROUTINE cc_gmres(x0, f0, delta, eta, m_max, it, nmv)
  USE parallel
  USE constants
  USE diis_mod
  USE mem_tracker

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(IN) :: x0(diis_lo:diis_hi), f0(diis_lo:diis_hi)
  COMPLEX(dpc), INTENT(OUT) :: delta(diis_lo:diis_hi)
  REAL(dp), INTENT(IN) :: eta
  INTEGER, INTENT(IN) :: m_max, it
  INTEGER, INTENT(OUT) :: nmv
  INTEGER :: i, j, pass, info, lwork
  REAL(dp) :: xnorm, beta, hn, resid
  COMPLEX(dpc) :: d
  COMPLEX(dpc), ALLOCATABLE :: vk(:,:), h(:,:), a(:,:), b(:), work(:)

  ALLOCATE( vk(diis_lo:diis_hi, 0:m_max) )
  CALL mem_register('t2', REAL(diis_counts(iam),dp) * (m_max+1) * 16.0_dp)
  ALLOCATE( h(m_max+1, m_max), a(m_max+1, m_max), b(m_max+1) )
  lwork = 64*(m_max+1)
  ALLOCATE( work(lwork) )
  h = 0.0_dp

  CALL slice_dot(x0, x0, d)
  xnorm = SQRT(REAL(d))
  CALL slice_dot(f0, f0, d)
  beta = SQRT(REAL(d))
  vk(:,0) = -f0 / beta

  b = 0.0_dp
  nmv = 0
  DO j = 1, m_max
     CALL cc_jv(x0, f0, xnorm, vk(:,j-1), vk(:,j))
     nmv = j
     ! Arnoldi, modified Gram-Schmidt done twice
     DO pass = 1, 2
        DO i = 0, j-1
           CALL slice_dot(vk(:,i), vk(:,j), d)
           h(i+1,j) = h(i+1,j) + d
           vk(:,j) = vk(:,j) - d*vk(:,i)
        end DO
     end DO
     CALL slice_dot(vk(:,j), vk(:,j), d)
     hn = SQRT(REAL(d))
     h(j+1,j) = hn
     IF ( hn > 1.d-14*beta ) vk(:,j) = vk(:,j) / hn

     ! min | beta.e1 - H y | over the (j+1) x j Hessenberg matrix
     a(1:j+1,1:j) = h(1:j+1,1:j)
     b = 0.0_dp
     b(1) = beta
     CALL ZGELS('N', j+1, j, 1, a, m_max+1, b, m_max+1, work, lwork, info)
     resid = ABS(b(j+1))
     IF ( iam == 0 ) write(6,'(A,I4,A,I4,A,ES11.3)') ' ...GMRES        step', it, '   j', j, '   |r|/|F| =', resid/beta
     IF ( resid <= eta*beta .or. hn <= 1.d-14*beta ) EXIT
  end DO

  delta = 0.0_dp
  DO i = 1, nmv
     delta = delta + b(i)*vk(:,i-1)
  end DO

  DEALLOCATE( vk, h, a, b, work )
  CALL mem_release('t2', REAL(diis_counts(iam),dp) * (m_max+1) * 16.0_dp)

end SUBROUTINE cc_gmres
