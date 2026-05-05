!
! MODULE contracts
!
! Runtime assertion helpers for enforcing structural invariants across the
! code. The goal is to make implicit assumptions (orbit ordering, array
! allocation, build status, channel validity) explicit and checkable, so
! that a violation fails loudly at the routine that depends on it rather
! than silently propagating through the calculation.
!
! All checks are guarded by the `debug_contracts` logical flag (namelist
! controlled, default TRUE) so they can be disabled for production runs
! where the setup is trusted.
!
! Usage pattern inside a consumer routine:
!
!    CALL assert_built('t2', 'my_subroutine_name')
!    CALL assert_channel_valid(ch, 'my_subroutine_name')
!    CALL assert_holes_ordered('my_subroutine_name')
!
! On violation, contracts print a diagnostic with rank, location, and the
! specific invariant that failed, then call MPI_ABORT.
!
MODULE contracts
  USE build_status
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PUBLIC :: debug_contracts = .TRUE.
  ! namelist /args/ debug_contracts  (added by ccdt_main)

  PUBLIC :: assert_built
  PUBLIC :: assert_holes_ordered
  PUBLIC :: assert_channel_valid
  PUBLIC :: assert_channel_2b_valid
  PUBLIC :: assert_channel_t3_valid
  PUBLIC :: assert_allocated_c2
  PUBLIC :: assert_allocated_c1
  PUBLIC :: assert_positive
  PUBLIC :: contract_fail

CONTAINS

  !
  ! Abort with a formatted diagnostic. All contract failures route through
  ! here so the error format is consistent.
  !
  SUBROUTINE contract_fail(where, what)
    USE parallel, ONLY: iam, ierror
    CHARACTER(LEN=*), INTENT(IN) :: where, what
    INCLUDE 'mpif.h'
    WRITE(0,'(A,I6,A,A,A,A)') 'RANK ', iam, ' CONTRACT VIOLATION in ', &
         TRIM(where), ': ', TRIM(what)
    CALL FLUSH(0)
    CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierror)
  END SUBROUTINE contract_fail

  !
  ! Assert that a named build stage has completed.
  ! Recognized tags: 'basis', 'interactions', 'mappings', 'fock',
  !                  't2', 't3', 'hbar', 'channels_2b', 'channels_t3'
  !
  SUBROUTINE assert_built(tag, where)
    CHARACTER(LEN=*), INTENT(IN) :: tag, where
    LOGICAL :: ok
    IF ( .NOT. debug_contracts ) RETURN
    ok = .FALSE.
    SELECT CASE (TRIM(tag))
    CASE ('basis');        ok = basis_built
    CASE ('interactions'); ok = interactions_built
    CASE ('mappings');     ok = mappings_built
    CASE ('fock');         ok = fock_built
    CASE ('t2');           ok = t2_built
    CASE ('t3');           ok = t3_built
    CASE ('hbar');         ok = hbar_built
    CASE ('channels_2b');  ok = channels_2b_built
    CASE ('channels_t3');  ok = channels_t3_built
    CASE DEFAULT
       CALL contract_fail(where, 'unknown build tag: '//TRIM(tag))
    END SELECT
    IF ( .NOT. ok ) CALL contract_fail(where, 'required stage not built: '//TRIM(tag))
  END SUBROUTINE assert_built

  !
  ! Assert that hole orbits occupy indices 1..below_ef. This is the single
  ! most important structural invariant in the code; every lookup table
  ! allocation and every configuration loop depends on it.
  !
  SUBROUTINE assert_holes_ordered(where)
    USE single_particle_orbits, ONLY: all_orbit
    USE constants, ONLY: below_ef, tot_orbs
    CHARACTER(LEN=*), INTENT(IN) :: where
    INTEGER :: i
    CHARACTER(LEN=128) :: msg
    IF ( .NOT. debug_contracts ) RETURN
    DO i = 1, below_ef
       IF ( all_orbit%orbit_status(i) /= 'hole' ) THEN
          WRITE(msg,'(A,I6,A)') 'orbit ', i, ' has index <= below_ef but is not a hole'
          CALL contract_fail(where, msg)
       END IF
    END DO
    DO i = below_ef+1, tot_orbs
       IF ( all_orbit%orbit_status(i) /= 'particle' ) THEN
          WRITE(msg,'(A,I6,A)') 'orbit ', i, ' has index > below_ef but is not a particle'
          CALL contract_fail(where, msg)
       END IF
    END DO
  END SUBROUTINE assert_holes_ordered

  !
  ! Assert that a channel index is within the 2b channel range.
  !
  SUBROUTINE assert_channel_valid(ch, where)
    USE configurations, ONLY: channels_2b
    INTEGER, INTENT(IN) :: ch
    CHARACTER(LEN=*), INTENT(IN) :: where
    CHARACTER(LEN=128) :: msg
    IF ( .NOT. debug_contracts ) RETURN
    IF ( ch < 1 .OR. ch > channels_2b%number_confs ) THEN
       WRITE(msg,'(A,I8,A,I8)') 'channel out of range: ch=', ch, &
            ' max=', channels_2b%number_confs
       CALL contract_fail(where, msg)
    END IF
  END SUBROUTINE assert_channel_valid

  !
  ! Assert that a 2b channel has nonzero configs of the given struct type.
  ! struct: 1=hh, 2=hp, 3=pp
  !
  SUBROUTINE assert_channel_2b_valid(ch, struct, where)
    USE operator_storage, ONLY: number_2b
    INTEGER, INTENT(IN) :: ch, struct
    CHARACTER(LEN=*), INTENT(IN) :: where
    CHARACTER(LEN=128) :: msg
    IF ( .NOT. debug_contracts ) RETURN
    CALL assert_channel_valid(ch, where)
    IF ( struct < 1 .OR. struct > 3 ) THEN
       WRITE(msg,'(A,I4)') 'invalid struct type: ', struct
       CALL contract_fail(where, msg)
    END IF
    IF ( number_2b(struct,ch) <= 0 ) THEN
       WRITE(msg,'(A,I8,A,I2,A)') 'channel ', ch, ' has no struct=', struct, ' configs'
       CALL contract_fail(where, msg)
    END IF
  END SUBROUTINE assert_channel_2b_valid

  !
  ! Assert that a T3 channel index is within this rank's local range.
  !
  SUBROUTINE assert_channel_t3_valid(ch3, where)
    USE operator_storage, ONLY: ch3_min, ch3_max
    INTEGER, INTENT(IN) :: ch3
    CHARACTER(LEN=*), INTENT(IN) :: where
    CHARACTER(LEN=128) :: msg
    IF ( .NOT. debug_contracts ) RETURN
    IF ( ch3 < ch3_min .OR. ch3 > ch3_max ) THEN
       WRITE(msg,'(A,I8,A,I8,A,I8)') 'T3 channel out of local range: ch3=', ch3, &
            ' min=', ch3_min, ' max=', ch3_max
       CALL contract_fail(where, msg)
    END IF
  END SUBROUTINE assert_channel_t3_valid

  !
  ! Assert that a complex 2D array is allocated.
  !
  SUBROUTINE assert_allocated_c2(arr, name, where)
    COMPLEX(KIND(1.0d0)), ALLOCATABLE, INTENT(IN) :: arr(:,:)
    CHARACTER(LEN=*), INTENT(IN) :: name, where
    IF ( .NOT. debug_contracts ) RETURN
    IF ( .NOT. ALLOCATED(arr) ) &
         CALL contract_fail(where, 'array not allocated: '//TRIM(name))
  END SUBROUTINE assert_allocated_c2

  !
  ! Assert that a complex 1D array is allocated.
  !
  SUBROUTINE assert_allocated_c1(arr, name, where)
    COMPLEX(KIND(1.0d0)), ALLOCATABLE, INTENT(IN) :: arr(:)
    CHARACTER(LEN=*), INTENT(IN) :: name, where
    IF ( .NOT. debug_contracts ) RETURN
    IF ( .NOT. ALLOCATED(arr) ) &
         CALL contract_fail(where, 'array not allocated: '//TRIM(name))
  END SUBROUTINE assert_allocated_c1

  !
  ! Assert that an integer value is positive (useful for dimension checks).
  !
  SUBROUTINE assert_positive(val, name, where)
    INTEGER, INTENT(IN) :: val
    CHARACTER(LEN=*), INTENT(IN) :: name, where
    CHARACTER(LEN=128) :: msg
    IF ( .NOT. debug_contracts ) RETURN
    IF ( val <= 0 ) THEN
       WRITE(msg,'(A,A,A,I12)') 'non-positive value for ', TRIM(name), ': ', val
       CALL contract_fail(where, msg)
    END IF
  END SUBROUTINE assert_positive

END MODULE contracts
