!
! MODULE mem_tracker
!
! Tagged memory accounting to replace the scattered manual `_mem0 += ...`
! pattern. Each allocation site registers its size under a named tag; the
! module maintains per-tag running totals (per-rank) and per-tag peak
! values. Reports can be requested at any point and produce an
! MPI-reduced summary (min/max/total across ranks) for each tag.
!
! Design notes:
! - Tags are short strings, looked up linearly. With ~20 tags this is
!   faster than a hash map and easier to debug.
! - Units: all sizes are in BYTES internally. Reports convert to GB.
! - This module is ADDITIVE for v1. The old `_mem0` counters in
!   modules.f90 are left in place so the two accounting systems can be
!   compared during validation. Remove the old ones in v2.
!
! Usage:
!   CALL mem_register('t2_ccm',  nbytes)        ! on allocate
!   CALL mem_release ('t2_ccm',  nbytes)        ! on deallocate
!   CALL mem_report  ('after setup')            ! rank 0 prints summary
!
MODULE mem_tracker
  USE kind_params, ONLY: dp
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: MAX_TAGS = 32
  INTEGER, PARAMETER :: TAG_LEN  = 24

  TYPE :: mem_entry
     CHARACTER(LEN=TAG_LEN) :: tag  = ' '
     REAL(dp)               :: bytes = 0.0_dp
     REAL(dp)               :: peak  = 0.0_dp
     LOGICAL                :: used  = .FALSE.
  END TYPE mem_entry

  TYPE(mem_entry), SAVE :: entries(MAX_TAGS)
  INTEGER,         SAVE :: n_entries = 0

  PUBLIC :: mem_register, mem_release, mem_report, mem_total_local

CONTAINS

  !
  ! Find or create an entry for the given tag. Returns the index.
  !
  FUNCTION find_or_create(tag) RESULT(idx)
    CHARACTER(LEN=*), INTENT(IN) :: tag
    INTEGER :: idx, i
    DO i = 1, n_entries
       IF ( TRIM(entries(i)%tag) == TRIM(tag) ) THEN
          idx = i
          RETURN
       END IF
    END DO
    IF ( n_entries >= MAX_TAGS ) THEN
       WRITE(0,*) 'mem_tracker: MAX_TAGS exceeded, increase MAX_TAGS'
       idx = 1
       RETURN
    END IF
    n_entries = n_entries + 1
    entries(n_entries)%tag   = tag
    entries(n_entries)%bytes = 0.0_dp
    entries(n_entries)%peak  = 0.0_dp
    entries(n_entries)%used  = .TRUE.
    idx = n_entries
  END FUNCTION find_or_create

  !
  ! Register an allocation of `bytes` under the given tag.
  !
  SUBROUTINE mem_register(tag, bytes)
    CHARACTER(LEN=*), INTENT(IN) :: tag
    REAL(dp),         INTENT(IN) :: bytes
    INTEGER :: idx
    idx = find_or_create(tag)
    entries(idx)%bytes = entries(idx)%bytes + bytes
    IF ( entries(idx)%bytes > entries(idx)%peak ) &
         entries(idx)%peak = entries(idx)%bytes
  END SUBROUTINE mem_register

  !
  ! Register a deallocation of `bytes` under the given tag.
  !
  SUBROUTINE mem_release(tag, bytes)
    CHARACTER(LEN=*), INTENT(IN) :: tag
    REAL(dp),         INTENT(IN) :: bytes
    INTEGER :: idx
    idx = find_or_create(tag)
    entries(idx)%bytes = entries(idx)%bytes - bytes
    IF ( entries(idx)%bytes < 0.0_dp ) entries(idx)%bytes = 0.0_dp
  END SUBROUTINE mem_release

  !
  ! Return current total memory for this rank, in bytes.
  !
  FUNCTION mem_total_local() RESULT(total)
    REAL(dp) :: total
    INTEGER :: i
    total = 0.0_dp
    DO i = 1, n_entries
       total = total + entries(i)%bytes
    END DO
  END FUNCTION mem_total_local

  !
  ! Produce an MPI-reduced summary and print it on rank 0.
  ! Each tag shows: min / max / sum across ranks, in GB.
  !
  SUBROUTINE mem_report(section)
    USE parallel, ONLY: iam, ierror
    CHARACTER(LEN=*), INTENT(IN) :: section
    INCLUDE 'mpif.h'
    REAL(dp) :: vmin, vmax, vsum, cur, pk
    REAL(dp), PARAMETER :: GB = 1.0e9_dp
    INTEGER :: i
    IF ( iam == 0 ) THEN
       WRITE(6,*)
       WRITE(6,'(A)') '=============================================================='
       WRITE(6,'(A,A)') '  Memory report: ', TRIM(section)
       WRITE(6,'(A)') '=============================================================='
       WRITE(6,'(A24,3x,A12,3x,A12,3x,A12)') &
            'tag', 'min (GB)', 'max (GB)', 'sum (GB)'
       WRITE(6,'(A)') '--------------------------------------------------------------'
    END IF
    DO i = 1, n_entries
       cur = entries(i)%bytes
       pk  = entries(i)%peak
       vmin = cur; vmax = cur; vsum = cur
       CALL MPI_ALLREDUCE(cur, vmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
       CALL MPI_ALLREDUCE(cur, vmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
       CALL MPI_ALLREDUCE(cur, vsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
       IF ( iam == 0 .AND. vmax > 0.0_dp ) THEN
          WRITE(6,'(A24,3x,F12.4,3x,F12.4,3x,F12.4)') &
               entries(i)%tag, vmin/GB, vmax/GB, vsum/GB
       END IF
    END DO
    cur = mem_total_local()
    CALL MPI_ALLREDUCE(cur, vmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
    CALL MPI_ALLREDUCE(cur, vmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
    CALL MPI_ALLREDUCE(cur, vsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    IF ( iam == 0 ) THEN
       WRITE(6,'(A)') '--------------------------------------------------------------'
       WRITE(6,'(A24,3x,F12.4,3x,F12.4,3x,F12.4)') &
            'TOTAL', vmin/GB, vmax/GB, vsum/GB
       WRITE(6,'(A)') '=============================================================='
       WRITE(6,*)
    END IF
  END SUBROUTINE mem_report

END MODULE mem_tracker
