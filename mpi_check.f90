!
! MODULE mpi_check
!
! Lightweight MPI error checking. Wraps the most common MPI calls
! with automatic ierror checking. On failure, prints rank, location,
! operation, and MPI error string, then aborts.
!
! Usage:
!   USE mpi_check
!   CALL checked_allreduce(buf, count, datatype, op, comm, 'my_subroutine')
!
! Or for checking ierror after a raw MPI call:
!   CALL mpi_allreduce(mpi_in_place, buf, n, ..., ierror)
!   CALL check_mpi(ierror, 'my_subroutine', 'allreduce on buf')
!
MODULE mpi_check
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: check_mpi

CONTAINS

  !
  ! Check ierror from any MPI call. If nonzero, print diagnostic and abort.
  !
  SUBROUTINE check_mpi(ierr, where, what)
    USE parallel, ONLY: iam
    INCLUDE 'mpif.h'
    INTEGER, INTENT(IN) :: ierr
    CHARACTER(LEN=*), INTENT(IN) :: where, what
    CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: errstr
    INTEGER :: errlen, ierr2

    IF ( ierr == MPI_SUCCESS ) RETURN

    CALL MPI_ERROR_STRING(ierr, errstr, errlen, ierr2)
    WRITE(0,'(A,I6,A,A,A,A,A,A)') 'RANK ', iam, &
         ' MPI ERROR in ', TRIM(where), ': ', TRIM(what), &
         ' — ', errstr(1:errlen)
    CALL FLUSH(0)
    CALL MPI_ABORT(MPI_COMM_WORLD, ierr, ierr2)
  END SUBROUTINE check_mpi

END MODULE mpi_check
