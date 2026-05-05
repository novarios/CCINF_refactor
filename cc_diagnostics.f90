!
! MODULE cc_diagnostics
!
! Per-iteration convergence diagnostics for the CC solver. Tracks
! amplitude norms, energy changes, and iteration timing to help
! diagnose convergence failures and identify problematic channels.
!
! Usage in the CC iteration loop:
!   CALL diag_start_iteration(iter)       ! at top of iteration
!   CALL diag_t2_norm()                   ! after T2 update
!   CALL diag_t3_norm()                   ! after T3 update  
!   CALL diag_energy(e2, e3)              ! after energy evaluation
!   CALL diag_end_iteration()             ! at bottom of iteration
!
! Output (rank 0 only):
!   Prints a one-line summary per iteration with T2/T3 norms,
!   energy, delta-E, and wall time. If amplitudes are growing
!   (norm increasing), prints a warning.
!
MODULE cc_diagnostics
  USE kind_params
  IMPLICIT NONE
  PRIVATE

  INTEGER,  SAVE :: current_iter = 0
  REAL(dp), SAVE :: iter_start_time = 0.0_dp
  REAL(dp), SAVE :: prev_energy_re = 0.0_dp
  REAL(dp), SAVE :: prev_t2_norm = 0.0_dp
  REAL(dp), SAVE :: prev_t3_norm = 0.0_dp
  REAL(dp), SAVE :: last_t2_norm = 0.0_dp
  REAL(dp), SAVE :: last_t3_norm = 0.0_dp
  LOGICAL,  SAVE :: header_printed = .FALSE.

  PUBLIC :: diag_start_iteration, diag_t2_norm, diag_t3_norm
  PUBLIC :: diag_energy, diag_end_iteration

CONTAINS

  SUBROUTINE diag_start_iteration(iter)
    INCLUDE 'mpif.h'
    INTEGER, INTENT(IN) :: iter
    current_iter = iter
    iter_start_time = MPI_WTIME()
  END SUBROUTINE diag_start_iteration

  !
  ! Compute and store the Frobenius norm of all T2 amplitudes.
  !
  SUBROUTINE diag_t2_norm
    USE parallel, ONLY: iam, ierror
    USE operator_storage, ONLY: t2_ccm, number_2b
    USE configurations, ONLY: channels_2b
    INCLUDE 'mpif.h'
    INTEGER :: ch
    REAL(dp) :: local_norm, global_norm
    COMPLEX(dpc) :: z

    local_norm = 0.0_dp
    DO ch = 1, channels_2b%number_confs
       IF ( number_2b(3,ch) == 0 ) cycle
       IF ( number_2b(1,ch) == 0 ) cycle
       local_norm = local_norm + SUM(ABS(t2_ccm(ch)%cval)**2)
    END DO

    CALL MPI_ALLREDUCE(local_norm, global_norm, 1, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_COMM_WORLD, ierror)

    last_t2_norm = SQRT(global_norm)
  END SUBROUTINE diag_t2_norm

  !
  ! Compute and store the Frobenius norm of all T3 amplitudes.
  !
  SUBROUTINE diag_t3_norm
    USE parallel, ONLY: iam, ierror
    USE operator_storage, ONLY: t3_ccm, ch3_min, ch3_max, &
         climits_t3, klimit_t3, number_2b_t3, klist_t3
    INCLUDE 'mpif.h'
    INTEGER :: ch3, cind1, kind1, ch2, ket_confs
    REAL(dp) :: local_norm, global_norm

    local_norm = 0.0_dp
    DO ch3 = ch3_min, ch3_max
       DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
          DO kind1 = 1, klimit_t3(ch3)
             IF ( .NOT. ASSOCIATED(t3_ccm(ch3)%val2(cind1,kind1)%cval) ) cycle
             ch2 = klist_t3(ch3)%ival2(kind1,2)
             ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
             IF ( ket_confs <= 0 ) cycle
             local_norm = local_norm + SUM(ABS(t3_ccm(ch3)%val2(cind1,kind1)%cval)**2)
          END DO
       END DO
    END DO

    CALL MPI_ALLREDUCE(local_norm, global_norm, 1, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_COMM_WORLD, ierror)

    last_t3_norm = SQRT(global_norm)
  END SUBROUTINE diag_t3_norm

  !
  ! Record energy and print iteration summary.
  !
  SUBROUTINE diag_energy(e2, e3)
    USE parallel, ONLY: iam
    INCLUDE 'mpif.h'
    COMPLEX(dpc), INTENT(IN) :: e2, e3
    REAL(dp) :: e_re, delta_e, wtime
    LOGICAL :: t2_growing, t3_growing

    e_re = REAL(e2 + e3, dp)
    delta_e = e_re - prev_energy_re
    wtime = MPI_WTIME() - iter_start_time

    t2_growing = ( last_t2_norm > prev_t2_norm .AND. current_iter > 2 )
    t3_growing = ( last_t3_norm > prev_t3_norm .AND. current_iter > 2 &
         .AND. last_t3_norm > 0.0_dp )

    IF ( iam == 0 ) THEN
       IF ( .NOT. header_printed ) THEN
          WRITE(6,*)
          WRITE(6,'(A)') '  iter    ||T2||        ||T3||        E_corr' // &
               '           delta_E        time(s)'
          WRITE(6,'(A)') '  ----  ----------   -----------   --------' // &
               '------   -----------   --------'
          header_printed = .TRUE.
       END IF

       WRITE(6,'(I6,2x,ES12.5,1x,ES13.5,1x,F16.10,1x,ES13.5,1x,F8.2)', advance='no') &
            current_iter, last_t2_norm, last_t3_norm, e_re, delta_e, wtime

       IF ( t2_growing .OR. t3_growing ) THEN
          WRITE(6,'(A)', advance='no') '  *** WARNING:'
          IF ( t2_growing ) WRITE(6,'(A)', advance='no') ' T2 growing'
          IF ( t3_growing ) WRITE(6,'(A)', advance='no') ' T3 growing'
       END IF
       WRITE(6,*)
    END IF

    prev_energy_re = e_re
    prev_t2_norm = last_t2_norm
    prev_t3_norm = last_t3_norm
  END SUBROUTINE diag_energy

  SUBROUTINE diag_end_iteration
    ! placeholder for future per-iteration logging (e.g. to file)
  END SUBROUTINE diag_end_iteration

END MODULE cc_diagnostics
