!
! MODULE build_status
!
! Tracks which stages of the setup pipeline have been completed. Consumers
! (t2_eqn, t3_eqn, build_hbar_t2_iter, cc_energy, etc.) assert on these
! flags at entry via the contracts module, so a routine called before its
! prerequisites are built fails loudly with a clear message instead of
! silently reading uninitialized memory.
!
! Producers set their flag at the END of their setup routine (after all
! allreduces complete) and clear it at the START of teardown. The flags
! are all PUBLIC — they are meant to be written only by producers and read
! only via the contracts module.
!
! Adding a new stage:
!   1. Add a logical flag here, default .FALSE.
!   2. Set it .TRUE. at the end of the setup routine.
!   3. Clear it in the corresponding teardown.
!   4. Add a case to assert_built in contracts.f90.
!
MODULE build_status
  IMPLICIT NONE
  PUBLIC

  LOGICAL :: basis_built        = .FALSE.
  LOGICAL :: channels_2b_built  = .FALSE.
  LOGICAL :: channels_t3_built  = .FALSE.
  LOGICAL :: mappings_built     = .FALSE.
  LOGICAL :: interactions_built = .FALSE.
  LOGICAL :: fock_built         = .FALSE.
  LOGICAL :: t2_built           = .FALSE.
  LOGICAL :: t3_built           = .FALSE.
  LOGICAL :: hbar_built         = .FALSE.

CONTAINS

  !
  ! Clear every flag. Called on a hard restart (e.g. before a second run
  ! in the same process, if that ever becomes a use case). Not currently
  ! called from anywhere but useful for future test harnesses.
  !
  SUBROUTINE clear_all_build_flags
    basis_built        = .FALSE.
    channels_2b_built  = .FALSE.
    channels_t3_built  = .FALSE.
    mappings_built     = .FALSE.
    interactions_built = .FALSE.
    fock_built         = .FALSE.
    t2_built           = .FALSE.
    t3_built           = .FALSE.
    hbar_built         = .FALSE.
  END SUBROUTINE clear_all_build_flags

END MODULE build_status
