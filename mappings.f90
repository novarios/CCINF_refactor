
SUBROUTINE setup_proc_mappings_v2b
  use parallel
  use constants
  use configurations
  use operator_storage
  USE mem_tracker
  implicit none
  
  INTEGER(i8) :: ltot_work
  INTEGER(i8) :: work_pr_proc, curr_work, temp
  INTEGER :: curr_proc
  INTEGER :: ch, ii, bra_confs, ket_confs, bra_min,bra_max
  LOGICAL :: has_ch_been_added
  
  ! Distribute bra for all  

  ! pppp
  ltot_work = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) CYCLE
     IF ( number_2b(1,ch) == 0 ) CYCLE ! Only connects to pphh
     bra_confs = number_2b(3,ch)
     ket_confs = bra_confs

     ltot_work = ltot_work + int(bra_confs,8)*int(ket_confs,8)        
  end DO

  work_pr_proc = int( ltot_work/int(num_procs,8) )
  ALLOCATE( mapping_v2b_pppp(1:num_procs, 1:channels_2b%number_confs, 1:3) )
  mapping_v2b_pppp = 0
  
  curr_work = 0
  curr_proc = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) CYCLE
     IF ( number_2b(1,ch) == 0 ) CYCLE ! Only connects to pphh
     bra_confs = number_2b(3,ch)
     ket_confs = bra_confs

     has_ch_been_added = .FALSE.        
     ii = 0
     DO WHILE ( ii < bra_confs )
        ii = ii + 1
        temp = curr_work + ket_confs
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( has_ch_been_added ) THEN
              mapping_v2b_pppp(curr_proc+1, ch, 3) = ii
           ELSE
              has_ch_been_added = .true.
              mapping_v2b_pppp(curr_proc+1, ch, 1) = ch
              mapping_v2b_pppp(curr_proc+1, ch, 2) = ii
              mapping_v2b_pppp(curr_proc+1, ch, 3) = ii
           end IF
        ELSE
           has_ch_been_added = .false.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF
     end DO
  end DO
  ! IF ( iam == 0 ) write(6,*)
  ! IF ( iam == 0 ) write(6,'(A10)',advance="no") 'PPPP_MAP: '
  ! DO ii = 1, num_procs
  !    curr_work = 0
  !    DO ch = 1, channels_2b%number_confs
  !       IF ( mapping_v2b_pppp(ii,ch,1) == 0 ) cycle
  !       bra_min = mapping_v2b_pppp(ii,ch,2)
  !       bra_max = mapping_v2b_pppp(ii,ch,3)
  !       ket_confs = number_2b(3,ch)
  !       curr_work = curr_work + (bra_max-bra_min+1)*ket_confs
  !    end DO
  !    IF ( iam == 0 ) write(6,'(F7.5)',advance="no") real(curr_work)/real(ltot_work)
  ! end DO
  ! IF ( iam == 0 ) write(6,*)
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  ! pphp
  ltot_work = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) CYCLE
     IF ( number_2b(2,ch) == 0 ) CYCLE
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(2,ch)

     ltot_work = ltot_work + int(bra_confs,8)*int(ket_confs,8)        
  end DO

  work_pr_proc = int( ltot_work/int(num_procs,8) )
  ALLOCATE( mapping_v2b_pphp(1:num_procs, 1:channels_2b%number_confs, 1:3) )
  mapping_v2b_pphp = 0
  
  curr_work = 0
  curr_proc = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) CYCLE
     IF ( number_2b(2,ch) == 0 ) CYCLE
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(2,ch)

     has_ch_been_added = .FALSE.
     ii = 0
     DO WHILE ( ii < bra_confs )
        ii = ii + 1
        temp = curr_work + ket_confs
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( has_ch_been_added ) THEN
              mapping_v2b_pphp(curr_proc+1, ch, 3) = ii
           ELSE
              has_ch_been_added = .true.
              mapping_v2b_pphp(curr_proc+1, ch, 1) = ch
              mapping_v2b_pphp(curr_proc+1, ch, 2) = ii
              mapping_v2b_pphp(curr_proc+1, ch, 3) = ii
           end IF
        ELSE
           has_ch_been_added = .false.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF   
     end DO
  end DO
  ! IF ( iam == 0 ) write(6,*)
  ! IF ( iam == 0 ) write(6,'(A10)',advance="no") 'PPHP_MAP: '
  ! DO ii = 1, num_procs
  !    curr_work = 0
  !    DO ch = 1, channels_2b%number_confs
  !       IF ( mapping_v2b_pphp(ii,ch,1) == 0 ) cycle
  !       bra_min = mapping_v2b_pphp(ii,ch,2)
  !       bra_max = mapping_v2b_pphp(ii,ch,3)
  !       ket_confs = number_2b(2,ch)
  !       curr_work = curr_work + (bra_max-bra_min+1)*ket_confs
  !    end DO
  !    IF ( iam == 0 ) write(6,'(F7.5)',advance="no") real(curr_work)/real(ltot_work)
  ! end DO
  ! IF ( iam == 0 ) write(6,*)
  CALL mpi_barrier(mpi_comm_world,ierror)

  ! hphp
  ltot_work = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(2,ch) == 0 ) CYCLE
     bra_confs = number_2b(2,ch)
     ket_confs = bra_confs

     ltot_work = ltot_work + int(bra_confs,8)*int(ket_confs,8)
  end DO

  work_pr_proc = int( ltot_work/int(num_procs,8) )
  ALLOCATE( mapping_v2b_hphp(1:num_procs, 1:channels_2b%number_confs, 1:3) )
  mapping_v2b_hphp = 0
  
  curr_work = 0
  curr_proc = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(2,ch) == 0 ) CYCLE
     bra_confs = number_2b(2,ch)
     ket_confs = bra_confs

     has_ch_been_added = .FALSE.
     ii = 0
     DO WHILE ( ii < bra_confs )
        ii = ii + 1
        temp = curr_work + ket_confs
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( has_ch_been_added ) THEN
              mapping_v2b_hphp(curr_proc+1, ch, 3) = ii
           ELSE
              has_ch_been_added = .true.
              mapping_v2b_hphp(curr_proc+1, ch, 1) = ch
              mapping_v2b_hphp(curr_proc+1, ch, 2) = ii
              mapping_v2b_hphp(curr_proc+1, ch, 3) = ii
           end IF
        ELSE
           has_ch_been_added = .false.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF
     end DO
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  ! pphh
  ltot_work = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) CYCLE
     IF ( number_2b(1,ch) == 0 ) CYCLE
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)

     ltot_work = ltot_work + int(bra_confs,8)*int(ket_confs,8)
  end DO

  work_pr_proc = int( ltot_work/int(num_procs,8) )
  ALLOCATE( mapping_v2b_pphh(1:num_procs, 1:channels_2b%number_confs, 1:3) )
  mapping_v2b_pphh = 0
  
  curr_work = 0
  curr_proc = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) CYCLE
     IF ( number_2b(1,ch) == 0 ) CYCLE
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)

     has_ch_been_added = .FALSE.
     ii = 0
     DO WHILE ( ii < bra_confs )
        ii = ii + 1
        temp = curr_work + ket_confs
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( has_ch_been_added ) THEN
              mapping_v2b_pphh(curr_proc+1, ch, 3) = ii
           ELSE
              has_ch_been_added = .true.
              mapping_v2b_pphh(curr_proc+1, ch, 1) = ch
              mapping_v2b_pphh(curr_proc+1, ch, 2) = ii
              mapping_v2b_pphh(curr_proc+1, ch, 3) = ii
           end IF
        ELSE
           has_ch_been_added = .false.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF           
     end DO
  end DO
  ! IF ( iam == 0 ) write(6,*)
  ! IF ( iam == 0 ) write(6,'(A10)',advance="no") 'PPHH_MAP: '
  ! DO ii = 1, num_procs
  !    curr_work = 0
  !    DO ch = 1, channels_2b%number_confs
  !       IF ( mapping_v2b_pphh(ii,ch,1) == 0 ) cycle
  !       bra_min = mapping_v2b_pphh(ii,ch,2)
  !       bra_max = mapping_v2b_pphh(ii,ch,3)
  !       ket_confs = number_2b(1,ch)
  !       curr_work = curr_work + (bra_max-bra_min+1)*ket_confs
  !    end DO
  !    IF ( iam == 0 ) write(6,'(F7.5)',advance="no") real(curr_work)/real(ltot_work)
  ! end DO
  ! IF ( iam == 0 ) write(6,*)
  CALL mpi_barrier(mpi_comm_world,ierror)

  ! hphh
  ltot_work = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(2,ch) == 0 ) CYCLE
     IF ( number_2b(1,ch) == 0 ) CYCLE
     bra_confs = number_2b(2,ch)
     ket_confs = number_2b(1,ch)

     ltot_work = ltot_work + int(bra_confs,8)*int(ket_confs,8)        
  end DO

  work_pr_proc = int( ltot_work/int(num_procs,8) )
  ALLOCATE( mapping_v2b_hphh(1:num_procs, 1:channels_2b%number_confs, 1:3) )
  mapping_v2b_hphh = 0
  
  curr_work = 0
  curr_proc = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(2,ch) == 0 ) CYCLE
     IF ( number_2b(1,ch) == 0 ) CYCLE
     bra_confs = number_2b(2,ch)
     ket_confs = number_2b(1,ch)

     has_ch_been_added = .FALSE.
     ii = 0
     DO WHILE ( ii < bra_confs )
        ii = ii + 1
        temp = curr_work + ket_confs
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( has_ch_been_added ) THEN
              mapping_v2b_hphh(curr_proc+1, ch, 3) = ii
           ELSE
              has_ch_been_added = .true.
              mapping_v2b_hphh(curr_proc+1, ch, 1) = ch
              mapping_v2b_hphh(curr_proc+1, ch, 2) = ii
              mapping_v2b_hphh(curr_proc+1, ch, 3) = ii
           end IF
        ELSE
           has_ch_been_added = .false.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF   
     end DO
  end DO
  ! IF ( iam == 0 ) write(6,*)
  ! IF ( iam == 0 ) write(6,'(A10)',advance="no") 'HPHH_MAP: '
  ! DO ii = 1, num_procs
  !    curr_work = 0
  !    DO ch = 1, channels_2b%number_confs
  !       IF ( mapping_v2b_hphh(ii,ch,1) == 0 ) cycle
  !       bra_min = mapping_v2b_hphh(ii,ch,2)
  !       bra_max = mapping_v2b_hphh(ii,ch,3)
  !       ket_confs = number_2b(1,ch)
  !       curr_work = curr_work + (bra_max-bra_min+1)*ket_confs
  !    end DO
  !    IF ( iam == 0 ) write(6,'(F7.5)',advance="no") real(curr_work)/real(ltot_work)
  ! end DO
  ! IF ( iam == 0 ) write(6,*)
  CALL mpi_barrier(mpi_comm_world,ierror)

  ! hhhh
  ltot_work = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) CYCLE
     bra_confs = number_2b(1,ch)
     ket_confs = bra_confs

     ltot_work = ltot_work + int(bra_confs,8)*int(ket_confs,8)        
  end DO

  work_pr_proc = int( ltot_work/int(num_procs,8) )
  ALLOCATE( mapping_v2b_hhhh(1:num_procs, 1:channels_2b%number_confs, 1:3) )
  mapping_v2b_hhhh = 0
  
  curr_work = 0
  curr_proc = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) CYCLE
     bra_confs = number_2b(1,ch)
     ket_confs = bra_confs

     has_ch_been_added = .FALSE.        
     ii = 0
     DO WHILE ( ii < bra_confs )
        ii = ii + 1
        temp = curr_work + ket_confs
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( has_ch_been_added ) THEN
              mapping_v2b_hhhh(curr_proc+1, ch, 3) = ii
           ELSE
              has_ch_been_added = .true.
              mapping_v2b_hhhh(curr_proc+1, ch, 1) = ch
              mapping_v2b_hhhh(curr_proc+1, ch, 2) = ii
              mapping_v2b_hhhh(curr_proc+1, ch, 3) = ii
           end IF
        ELSE
           has_ch_been_added = .false.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF
     end DO
  end DO
  ! IF ( iam == 0 ) write(6,*)
  ! IF ( iam == 0 ) write(6,'(A10)',advance="no") 'HHHH_MAP: '
  ! DO ii = 1, num_procs
  !    curr_work = 0
  !    DO ch = 1, channels_2b%number_confs
  !       IF ( mapping_v2b_hhhh(ii,ch,1) == 0 ) cycle
  !       bra_min = mapping_v2b_hhhh(ii,ch,2)
  !       bra_max = mapping_v2b_hhhh(ii,ch,3)
  !       ket_confs = number_2b(1,ch)
  !       curr_work = curr_work + (bra_max-bra_min+1)*ket_confs
  !    end DO
  !    IF ( iam == 0 ) write(6,'(F7.5)',advance="no") real(curr_work)/real(ltot_work)
  ! end DO
  ! IF ( iam == 0 ) write(6,*)
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  ! pppp
  ALLOCATE( check_my_channel_v2b_pppp(1:channels_2b%number_confs) )
  check_my_channel_v2b_pppp = 0
  DO ch = 1, channels_2b%number_confs
     IF ( mapping_v2b_pppp(iam+1,ch,1) == 0 ) cycle
     check_my_channel_v2b_pppp(ch) = 1
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  ! pphp
  ALLOCATE( check_my_channel_v2b_pphp(1:channels_2b%number_confs) )
  check_my_channel_v2b_pphp = 0
  DO ch = 1, channels_2b%number_confs
     IF ( mapping_v2b_pphp(iam+1,ch,1) == 0 ) cycle
     check_my_channel_v2b_pphp(ch) = 1
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  ! hphp
  ALLOCATE( check_my_channel_v2b_hphp(1:channels_2b%number_confs) )
  check_my_channel_v2b_hphp = 0
  DO ch = 1, channels_2b%number_confs
     IF ( mapping_v2b_hphp(iam+1,ch,1) == 0 ) cycle
     check_my_channel_v2b_hphp(ch) = 1
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  ! pphh
  ALLOCATE( check_my_channel_v2b_pphh(1:channels_2b%number_confs) )
  check_my_channel_v2b_pphh = 0
  DO ch = 1, channels_2b%number_confs
     IF ( mapping_v2b_pphh(iam+1,ch,1) == 0 ) cycle
     check_my_channel_v2b_pphh(ch) = 1
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  ! hphh
  ALLOCATE( check_my_channel_v2b_hphh(1:channels_2b%number_confs) )
  check_my_channel_v2b_hphh = 0
  DO ch = 1, channels_2b%number_confs
     IF ( mapping_v2b_hphh(iam+1,ch,1) == 0 ) cycle
     check_my_channel_v2b_hphh(ch) = 1
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  ! hhhh
  ALLOCATE( check_my_channel_v2b_hhhh(1:channels_2b%number_confs) )
  check_my_channel_v2b_hhhh = 0
  DO ch = 1, channels_2b%number_confs
     IF ( mapping_v2b_hhhh(iam+1,ch,1) == 0 ) cycle
     check_my_channel_v2b_hhhh(ch) = 1
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  CALL mem_register('v2', REAL(6 * channels_2b%number_confs * num_procs * 3 * 4.d0, dp))
  CALL mem_register('v2', REAL(6 * channels_2b%number_confs * 4.d0, dp))
  IF ( iam == 0 ) write(6,*)
  
END SUBROUTINE setup_proc_mappings_v2b


SUBROUTINE setup_proc_mappings_v2b_cross
  use parallel
  use constants
  use configurations
  use operator_storage
  USE mem_tracker
  implicit none
  
  INTEGER(i8) :: ltot_work
  INTEGER(i8) :: work_pr_proc, curr_work, temp
  INTEGER :: curr_proc
  INTEGER :: ch, ii, bra_confs, ket_confs, bra_min,bra_max
  LOGICAL :: has_ch_been_added
  
  ! Distribute bra for all  

  ! hphp
  ltot_work = 0
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch) == 0 ) CYCLE
     bra_confs = number_2bcross(2,ch)
     ket_confs = bra_confs

     ltot_work = ltot_work + int(bra_confs,8)*int(ket_confs,8)
  end DO

  work_pr_proc = int( ltot_work/int(num_procs,8) )
  ALLOCATE( mapping_v2b_hphp_cross(1:num_procs, 1:channels_2bcross%number_confs, 1:3) )
  mapping_v2b_hphp_cross = 0
  
  curr_work = 0
  curr_proc = 0
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch) == 0 ) CYCLE
     bra_confs = number_2bcross(2,ch)
     ket_confs = bra_confs

     has_ch_been_added = .FALSE.
     ii = 0
     DO WHILE ( ii < bra_confs )
        ii = ii + 1
        temp = curr_work + ket_confs
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( has_ch_been_added ) THEN
              mapping_v2b_hphp_cross(curr_proc+1, ch, 3) = ii
           ELSE
              has_ch_been_added = .true.
              mapping_v2b_hphp_cross(curr_proc+1, ch, 1) = ch
              mapping_v2b_hphp_cross(curr_proc+1, ch, 2) = ii
              mapping_v2b_hphp_cross(curr_proc+1, ch, 3) = ii
           end IF
        ELSE
           has_ch_been_added = .false.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF
     end DO
  end DO
  ! IF ( iam == 0 ) write(6,*)
  ! IF ( iam == 0 ) write(6,'(A10)',advance="no") 'HPHP_MAP: '
  ! DO ii = 1, num_procs
  !    curr_work = 0
  !    DO ch = 1, channels_2bcross%number_confs
  !       IF ( mapping_v2b_hphp_cross(ii,ch,1) == 0 ) cycle
  !       bra_min = mapping_v2b_hphp_cross(ii,ch,2)
  !       bra_max = mapping_v2b_hphp_cross(ii,ch,3)
  !       ket_confs = number_2bcross(2,ch)
  !       curr_work = curr_work + (bra_max-bra_min+1)*ket_confs
  !    end DO
  !    IF ( iam == 0 ) write(6,'(F7.5)',advance="no") real(curr_work)/real(ltot_work)
  ! end DO
  ! IF ( iam == 0 ) write(6,*)
  CALL mpi_barrier(mpi_comm_world,ierror)

  ! phhp
  ltot_work = 0
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(3,ch) == 0 ) CYCLE
     IF ( number_2bcross(2,ch) == 0 ) CYCLE
     bra_confs = number_2bcross(3,ch)
     ket_confs = number_2bcross(2,ch)

     ltot_work = ltot_work + int(bra_confs,8)*int(ket_confs,8)
  end DO

  work_pr_proc = int( ltot_work/int(num_procs,8) )
  ALLOCATE( mapping_v2b_phhp_cross(1:num_procs, 1:channels_2bcross%number_confs, 1:3) )
  mapping_v2b_phhp_cross = 0
  
  curr_work = 0
  curr_proc = 0
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(3,ch) == 0 ) CYCLE
     IF ( number_2bcross(2,ch) == 0 ) CYCLE
     bra_confs = number_2bcross(3,ch)
     ket_confs = number_2bcross(2,ch)

     has_ch_been_added = .FALSE.
     ii = 0
     DO WHILE ( ii < bra_confs )
        ii = ii + 1
        temp = curr_work + ket_confs
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( has_ch_been_added ) THEN
              mapping_v2b_phhp_cross(curr_proc+1, ch, 3) = ii
           ELSE
              has_ch_been_added = .true.
              mapping_v2b_phhp_cross(curr_proc+1, ch, 1) = ch
              mapping_v2b_phhp_cross(curr_proc+1, ch, 2) = ii
              mapping_v2b_phhp_cross(curr_proc+1, ch, 3) = ii
           end IF
        ELSE
           has_ch_been_added = .false.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF
     end DO
  end DO
  ! IF ( iam == 0 ) write(6,*)
  ! IF ( iam == 0 ) write(6,'(A10)',advance="no") 'PHHP_MAP: '
  ! DO ii = 1, num_procs
  !    curr_work = 0
  !    DO ch = 1, channels_2bcross%number_confs
  !       IF ( mapping_v2b_phhp_cross(ii,ch,1) == 0 ) cycle
  !       bra_min = mapping_v2b_phhp_cross(ii,ch,2)
  !       bra_max = mapping_v2b_phhp_cross(ii,ch,3)
  !       ket_confs = number_2bcross(2,ch)
  !       curr_work = curr_work + (bra_max-bra_min+1)*ket_confs
  !    end DO
  !    IF ( iam == 0 ) write(6,'(F7.5)',advance="no") real(curr_work)/real(ltot_work)
  ! end DO
  ! IF ( iam == 0 ) write(6,*)
  CALL mpi_barrier(mpi_comm_world,ierror)

  ! hphp
  ALLOCATE( check_my_channel_v2b_hphp_cross(1:channels_2bcross%number_confs) )
  check_my_channel_v2b_hphp_cross = 0
  DO ch = 1, channels_2bcross%number_confs
     IF ( mapping_v2b_hphp_cross(iam+1,ch,1) == 0 ) cycle
     check_my_channel_v2b_hphp_cross(ch) = 1
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  ! phhp
  ALLOCATE( check_my_channel_v2b_phhp_cross(1:channels_2bcross%number_confs) )
  check_my_channel_v2b_phhp_cross = 0
  DO ch = 1, channels_2bcross%number_confs
     IF ( mapping_v2b_phhp_cross(iam+1,ch,1) == 0 ) cycle
     check_my_channel_v2b_phhp_cross(ch) = 1
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  CALL mem_register('v2', REAL(2 * channels_2bcross%number_confs * num_procs * 3 * 4.d0, dp))
  CALL mem_register('v2', REAL(2 * channels_2bcross%number_confs * 4.d0, dp))
  IF ( iam == 0 ) write(6,*)
  
END SUBROUTINE setup_proc_mappings_v2b_cross


SUBROUTINE setup_proc_mappings_v2b_setup(which)
  use parallel
  use constants
  use configurations
  use operator_storage
  implicit none

  CHARACTER(len=4), INTENT(IN) :: which
  INTEGER(i8) :: ltot_work
  INTEGER(i8) :: work_pr_proc, curr_work, temp
  INTEGER :: curr_proc
  INTEGER :: ch, ii, bra_confs, ket_confs
  LOGICAL :: has_ch_been_added

  IF ( which == 'hphp' ) then
     ! hphp
     IF ( ALLOCATED(mapping_v2b_setup) ) DEALLOCATE( mapping_v2b_setup )
     ltot_work = 0
     DO ch = 1, channels_2bcross%number_confs
        IF ( number_2bcross(2,ch) == 0 ) CYCLE
        bra_confs = number_2bcross(2,ch)
        ltot_work = ltot_work + int(bra_confs*(bra_confs+1),8)/2
     end DO
     write(6,*) 'HPHP0: ', iam, ltot_work
     work_pr_proc = int( ltot_work/int(num_procs,8) )
     ALLOCATE( mapping_v2b_setup(1:num_procs, 1:channels_2bcross%number_confs, 1:3) )
     mapping_v2b_setup = 0
     
     curr_work = 0
     curr_proc = 0
     DO ch = 1, channels_2bcross%number_confs
        IF ( number_2bcross(2,ch) == 0 ) CYCLE
        bra_confs = number_2bcross(2,ch)
        
        has_ch_been_added = .FALSE.
        ii = 0
        DO WHILE ( ii < bra_confs )
           temp = curr_work + (bra_confs-ii) ! triangular
           ii = ii + 1
           
           IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
              curr_work = temp
              IF ( has_ch_been_added ) THEN
                 mapping_v2b_setup(curr_proc+1, ch, 3) = ii
              ELSE
                 has_ch_been_added = .true.
                 mapping_v2b_setup(curr_proc+1, ch, 1) = ch
                 mapping_v2b_setup(curr_proc+1, ch, 2) = ii
                 mapping_v2b_setup(curr_proc+1, ch, 3) = ii
              end IF
           ELSE
              has_ch_been_added = .false.
              curr_work = 0
              curr_proc = curr_proc + 1
              ii = ii - 1
           end IF
        end DO
     end DO
  end IF
  
END SUBROUTINE setup_proc_mappings_v2b_setup


SUBROUTINE setup_proc_mappings0_t3
  use parallel
  use constants
  use configurations
  use operator_storage
  
  IMPLICIT NONE
  INTEGER :: ch3, curr_proc, ii
  INTEGER :: bra_confs,ket_confs
  INTEGER :: nx3,ny3,nz3,tz3
  INTEGER, allocatable :: mapping0(:,:)
  INTEGER(i8) :: ltot_work, curr_work, work_pr_proc, temp
  LOGICAL :: has_ch_been_added
  
  ! ppphhh
  ltot_work = 0
  DO ch3 = 1, channels_t3%number_confs
     nx3 = channels_t3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_t3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_t3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_t3%config_NxNyNz_Tz(ch3*4)

     bra_confs = num0_t3(1,nx3,ny3,nz3,tz3)
     IF ( bra_confs <= 0 ) cycle
     ket_confs = num0_t3(2,nx3,ny3,nz3,tz3)
     IF ( ket_confs <= 0 ) cycle
     
     ltot_work = ltot_work + int(bra_confs,8)*int(ket_confs,8)
  end DO
  work_pr_proc = int(ltot_work/int(num_procs,8),8)
  
  
  ALLOCATE( mapping0(1:num_procs, 1:channels_t3%number_confs) )
  mapping0 = 0

  curr_work = 0
  curr_proc = 0
  DO ch3 = 1, channels_t3%number_confs
     nx3 = channels_t3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_t3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_t3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_t3%config_NxNyNz_Tz(ch3*4)

     bra_confs = num0_t3(1,nx3,ny3,nz3,tz3)
     IF ( bra_confs <= 0 ) cycle
     ket_confs = num0_t3(2,nx3,ny3,nz3,tz3)
     IF ( ket_confs <= 0 ) cycle
     
     has_ch_been_added = .FALSE.
     ii = 0
     DO WHILE ( ii < bra_confs )
        ii = ii + 1
        temp = curr_work + int(ket_confs,8)
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( .not. has_ch_been_added ) THEN
              has_ch_been_added = .TRUE.
              mapping0(curr_proc+1, ch3) = 1
           end IF
        ELSE
           has_ch_been_added = .FALSE.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF
     end DO
  end DO
    
  ch3_min = channels_t3%number_confs
  ch3_max = 0
  DO ch3 = 1, channels_t3%number_confs
     IF ( mapping0(iam+1,ch3) == 0 ) cycle
     IF ( ch3 < ch3_min ) ch3_min = ch3
     IF ( ch3 > ch3_max ) ch3_max = ch3
  end DO
  DEALLOCATE( mapping0 )
  
END SUBROUTINE setup_proc_mappings0_t3

SUBROUTINE setup_proc_mappings1_t3
  use parallel
  use constants
  use configurations
  use operator_storage
  USE mem_tracker
  
  IMPLICIT NONE
  INTEGER :: curr_proc, num_chan, ii
  INTEGER :: ch3, ch1, bra0
  INTEGER :: bra_confs0,bra_confs, ket_confs
  INTEGER :: nx3,ny3,nz3,tz3
  INTEGER :: c, cind1, a,b
  INTEGER(i8) :: ltot_work, curr_work, work_pr_proc, temp
  LOGICAL :: has_ch_been_added
  
  ! ppphhh
  ltot_work = 0
  DO ch3 = 1, channels_t3%number_confs
     nx3 = channels_t3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_t3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_t3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_t3%config_NxNyNz_Tz(ch3*4)

     bra_confs = num0_t3(1,nx3,ny3,nz3,tz3)
     IF ( bra_confs <= 0 ) cycle
     ket_confs = num0_t3(2,nx3,ny3,nz3,tz3)
     IF ( ket_confs <= 0 ) cycle
     
     ltot_work = ltot_work + int(bra_confs,8)*int(ket_confs,8)
  end DO
  work_pr_proc = int(ltot_work/int(num_procs,8),8)

  
  ALLOCATE( climits_t3(ch3_min:ch3_max, 1:2) )
  CALL mem_register('t3', REAL(2 * (ch3_max-ch3_min+1) * 4.d0, dp))
  ALLOCATE( mapping_t3(ch3_min:ch3_max) )
  curr_work = 0
  curr_proc = 0
  num_chan  = 0
  DO ch3 = 1, channels_t3%number_confs
     nx3 = channels_t3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_t3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_t3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_t3%config_NxNyNz_Tz(ch3*4)

     IF ( ch3_min <= ch3 .and. ch3 <= ch3_max ) then
        climits_t3(ch3, 1) = climit_t3(ch3)+1
        climits_t3(ch3, 2) = 0
        ALLOCATE( mapping_t3(ch3)%ival2(climit_t3(ch3), 1:2) )
        mapping_t3(ch3)%ival2 = 0
        CALL mem_register('t3', REAL(2 * climit_t3(ch3) * 4.d0, dp))
     end IF
     ket_confs = num0_t3(2,nx3,ny3,nz3,tz3)
     IF ( ket_confs <= 0 ) cycle
     
     DO cind1 = 1, climit_t3(ch3)
        c    = clist_t3(ch3)%ival2(cind1,1)
        ch1  = clist_t3(ch3)%ival2(cind1,2)

        bra_confs = 0 ! count c+ab(t3_cut) configs
        bra_confs0 = number_2b(3,ch1)
        DO bra0 = 1, bra_confs0
           a = lookup_2b_configs(3,ch1)%ival2(1,bra0)
           b = lookup_2b_configs(3,ch1)%ival2(2,bra0)
           IF ( cut3b(a,b,c) ) cycle
           bra_confs = bra_confs + 1
        end DO
        
        has_ch_been_added = .FALSE.        
        ii = 0
        DO WHILE ( ii < bra_confs )
           ii = ii + 1
           temp = curr_work + int(ket_confs,8)
           IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) then
              curr_work = temp
              IF ( has_ch_been_added ) then
                 IF ( curr_proc == iam ) then
                    mapping_t3(ch3)%ival2(cind1,2) = ii
                    IF ( cind1 > climits_t3(ch3, 2) ) climits_t3(ch3, 2) = cind1
                 end IF
              ELSE
                 has_ch_been_added = .TRUE.
                 IF ( curr_proc == iam ) then
                    mapping_t3(ch3)%ival2(cind1,1) = ii
                    mapping_t3(ch3)%ival2(cind1,2) = ii
                    IF ( cind1 < climits_t3(ch3, 1) ) climits_t3(ch3, 1) = cind1
                    IF ( cind1 > climits_t3(ch3, 2) ) climits_t3(ch3, 2) = cind1
                 end IF
              end IF
           ELSE
              has_ch_been_added = .false.
              curr_work = 0
              curr_proc = curr_proc + 1
              ii = ii - 1
           end IF
        end DO
        
     end DO
  end DO
  
END SUBROUTINE setup_proc_mappings1_t3


SUBROUTINE setup_proc_mappings_t3_red
  use parallel
  use constants
  use configurations
  use operator_storage
  USE mem_tracker
  
  IMPLICIT NONE
  INTEGER :: ch3, curr_proc, ii
  INTEGER :: bra_confs,ket_confs
  INTEGER :: nx3,ny3,nz3,tz3
  INTEGER, allocatable :: mapping0(:,:)
  INTEGER(i8) :: ltot_work, curr_work, work_pr_proc, temp
  LOGICAL :: has_ch_been_added
  
  ! ppphhh
  ltot_work = 0
  DO ch3 = 1, channels_t3%number_confs
     nx3 = channels_t3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_t3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_t3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_t3%config_NxNyNz_Tz(ch3*4)

     CALL number_3b_red_confs(nx3, ny3, nz3, tz3, 2, ket_confs) ! 2 for hhh
     CALL number_3b_red_confs(nx3, ny3, nz3, tz3, 1, bra_confs) ! 1 for ppp
     
     ltot_work = ltot_work + int(bra_confs,8)*int(ket_confs,8)
  end DO
  work_pr_proc = int(ltot_work/int(num_procs,8),8)
    
  ALLOCATE( mapping0(1:num_procs, 1:channels_t3%number_confs) )
  mapping0 = 0

  curr_work = 0
  curr_proc = 0
  DO ch3 = 1, channels_t3%number_confs
     nx3 = channels_t3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_t3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_t3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_t3%config_NxNyNz_Tz(ch3*4)

     CALL number_3b_red_confs(nx3, ny3, nz3, tz3, 1, bra_confs) ! 1 for ppp
     CALL number_3b_red_confs(nx3, ny3, nz3, tz3, 2, ket_confs) ! 2 for hhh
     
     has_ch_been_added = .FALSE.
     ii = 0
     DO WHILE ( ii < bra_confs )
        ii = ii + 1
        temp = curr_work + int(ket_confs,8)
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( .not. has_ch_been_added ) THEN
              has_ch_been_added = .TRUE.
              mapping0(curr_proc+1, ch3) = 1
           end IF
        ELSE
           has_ch_been_added = .FALSE.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF
     end DO
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)
    
  ch3_min = channels_t3%number_confs
  ch3_max = 0
  DO ii = 1, num_procs
     DO ch3 = 1, channels_t3%number_confs
        IF ( mapping0(iam+1,ch3) == 0 ) cycle
        IF ( ch3 < ch3_min ) ch3_min = ch3
        IF ( ch3 > ch3_max ) ch3_max = ch3
     end DO
  end DO
  DEALLOCATE( mapping0 )

  ALLOCATE( mapping_t3_red(ch3_min:ch3_max, 1:2) )
  mapping_t3_red = 0
  CALL mem_register('t3', REAL((ch3_max-ch3_min+1) * 2 * 4.d0, dp))
  
  curr_work = 0
  curr_proc = 0
  DO ch3 = 1, channels_t3%number_confs
     bra_confs = number_3b_t3(ch3,1)
     ket_confs = number_3b_t3(ch3,2)

     IF ( ch3_min <= ch3 .and. ch3 <= ch3_max ) then
        mapping_t3_red(ch3,:) = 0
     end IF
     
     has_ch_been_added = .FALSE.        
     ii = 0
     DO WHILE ( ii < bra_confs )
        ii = ii + 1
        temp = curr_work + int(ket_confs,8)
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) then
           curr_work = temp
           IF ( has_ch_been_added ) then
              IF ( curr_proc == iam ) then
                 mapping_t3_red(ch3,2) = ii
              end IF
           ELSE
              has_ch_been_added = .true.
              IF ( curr_proc == iam ) then
                 mapping_t3_red(ch3,1) = ii
                 mapping_t3_red(ch3,2) = ii
              end IF
           end IF
        ELSE
           has_ch_been_added = .false.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF
     end DO
  end DO
  
END SUBROUTINE setup_proc_mappings_t3_red



SUBROUTINE setup_proc_mappings_paeom
  use parallel
  use configurations
  use operator_storage
  use constants
  USE mem_tracker
  implicit none
  
  INTEGER(i8) :: ltot_work
  INTEGER(i8) :: work_pr_proc, curr_work, temp
  INTEGER :: curr_proc
  INTEGER :: ch, ch2, ii, bra_confs, ket_confs
  LOGICAL :: has_ch_been_added
  
  ! pph
  ltot_work = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)  
     ltot_work = ltot_work + 2*int(bra_confs,8)
  end DO
  
  work_pr_proc = int( ltot_work/int(num_procs,8) )
  ALLOCATE( mapping_paeom_pph(1:num_procs, channels_2b%number_confs, 1:3) )
  mapping_paeom_pph = 0
  
  curr_work = 0
  curr_proc = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     has_ch_been_added = .FALSE.
     ii = 0
     DO WHILE ( ii < bra_confs )
        ii = ii + 1
        temp = curr_work + 2
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( has_ch_been_added ) THEN
              mapping_paeom_pph(curr_proc+1, ch, 3) = ii
           ELSE
              has_ch_been_added = .true.
              mapping_paeom_pph(curr_proc+1, ch, 1) = ch
              mapping_paeom_pph(curr_proc+1, ch, 2) = ii
              mapping_paeom_pph(curr_proc+1, ch, 3) = ii
           end IF
        ELSE
           has_ch_been_added = .false.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF
     end DO
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  ! pph
  ALLOCATE( check_my_channel_paeom_pph(channels_2b%number_confs) )
  check_my_channel_paeom_pph = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     IF ( mapping_paeom_pph(iam+1,ch,1) == 0 ) cycle
     check_my_channel_paeom_pph(ch) = 1
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  CALL mem_register('paeom2', REAL(channels_2b%number_confs * num_procs * 3 * 4.d0, dp))
  CALL mem_register('paeom2', REAL(channels_2b%number_confs * 4.d0, dp))


  ! php
  ltot_work = 0
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_paeom_ind_cross(2*ch2) == 0 ) cycle
     ket_confs = number_2bcross(2,ch2)
     ltot_work = ltot_work + 2*int(ket_confs,8)
  end DO
  
  work_pr_proc = int( ltot_work/int(num_procs,8) )
  ALLOCATE( mapping_paeom_php_cross(1:num_procs, channels_2bcross%number_confs, 1:3) )
  mapping_paeom_php_cross = 0
  
  curr_work = 0
  curr_proc = 0
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_paeom_ind_cross(2*ch2) == 0 ) cycle
     ket_confs = number_2bcross(2,ch2)
     has_ch_been_added = .FALSE.        
     ii = 0
     DO WHILE ( ii < ket_confs )
        ii = ii + 1
        temp = curr_work + 2
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( has_ch_been_added ) THEN
              mapping_paeom_php_cross(curr_proc+1, ch2, 3) = ii
           ELSE
              has_ch_been_added = .true.
              mapping_paeom_php_cross(curr_proc+1, ch2, 1) = ch2
              mapping_paeom_php_cross(curr_proc+1, ch2, 2) = ii
              mapping_paeom_php_cross(curr_proc+1, ch2, 3) = ii
           end IF
        ELSE
           has_ch_been_added = .false.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF
     end DO
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  ! pph
  ALLOCATE( check_my_channel_paeom_php_cross(channels_2bcross%number_confs) )
  check_my_channel_paeom_php_cross = 0
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_paeom_ind_cross(2*ch2) == 0 ) cycle
     IF ( mapping_paeom_php_cross(iam+1,ch2,1) == 0 ) cycle
     check_my_channel_paeom_php_cross(ch2) = 1
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  CALL mem_register('paeom2', REAL(channels_2bcross%number_confs * num_procs * 3 * 4.d0, dp))
  CALL mem_register('paeom2', REAL(channels_2bcross%number_confs * 4.d0, dp))
  
END SUBROUTINE setup_proc_mappings_paeom


SUBROUTINE setup_proc_mappings_paeom3
  use parallel
  use constants
  use configurations
  use operator_storage
  USE mem_tracker
  
  IMPLICIT NONE
  INTEGER(i8) :: ltot_work, ltot_work1, work_pr_proc, curr_work, temp
  INTEGER :: curr_proc, num_chan, ii
  INTEGER :: ch3, ch,ch1,ch2, bra_confs, ket_confs, bra0,bra_confs0
  INTEGER :: nx3,ny3,nz3,tz3, Nx1,Ny1,Nz1,Tz1
  INTEGER :: a,b, c, c_count, cind1, nxc,nyc,nzc,tzc
  INTEGER :: group_world, count_groups, procs_per_channel0
  INTEGER, allocatable :: mapping0(:,:)
  INTEGER, allocatable :: procs_per_channel(:)
  INTEGER, allocatable :: reduced_index(:)
  LOGICAL :: has_ch_been_added
  REAL(dp) :: mix1

  ! ppphh
  ltot_work = eom_ndim2 ! first proc has all r_p and r_pph amplitudes
  ltot_work1 = eom_ndim2 ! first proc has all r_p and r_pph amplitudes
  mix1 = 1.5d0 ! use for load balancing (0 = balanced r3, ~5 = balanced r3_vec?)
  DO ch3 = 1, channels_paeom3%number_confs
     nx3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_paeom3%config_NxNyNz_Tz(ch3*4)
     
     ch2 = ch2_paeom3(ch3)
     ket_confs = number_2b(1,ch2)
     
     CALL number_3b_confs(nx3, ny3, nz3, tz3, 1, bra_confs) ! 1 for ppp
     ltot_work = ltot_work + int(bra_confs,8)*int(ket_confs,8)
     CALL number_3b_red_confs(nx3, ny3, nz3, tz3, 1, bra_confs) ! 1 for ppp
     ltot_work1 = ltot_work1 + int(bra_confs,8)*int(ket_confs,8)
  end DO
  
  work_pr_proc = int( (ltot_work + int(mix1*ltot_work1,8)) /num_procs, 8)
  ALLOCATE( mapping0(1:num_procs, 1:channels_paeom3%number_confs) )
  mapping0 = 0
  
  curr_work = eom_ndim2 ! first proc has all r_p and r_pph amplitudes
  curr_proc = 0
  DO ch3 = 1, channels_paeom3%number_confs
     nx3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_paeom3%config_NxNyNz_Tz(ch3*4)
     
     CALL number_3b_confs(nx3, ny3, nz3, tz3, 1, bra_confs) ! 1 for ppp
     ch2 = ch2_paeom3(ch3)
     ket_confs = number_2b(1,ch2)

     ! !!!!!!!!!!!!!!!!!!!!!
     ALLOCATE( reduced_index(bra_confs) )
     reduced_index = 0
     ii = 0
     DO c = below_ef+1, tot_orbs
        nxc = all_orbit%nx(c)
        nyc = all_orbit%ny(c)
        nzc = all_orbit%nz(c)
        tzc = all_orbit%tz(c)
        DO ch1 = 1, channels_2b%number_confs
           IF ( number_2b(3,ch1) == 0 ) cycle
           Nx1 = channels_2b%config_NxNyNz_Tz(4*ch1-3)
           Ny1 = channels_2b%config_NxNyNz_Tz(4*ch1-2)
           Nz1 = channels_2b%config_NxNyNz_Tz(4*ch1-1)
           Tz1 = channels_2b%config_NxNyNz_Tz(4*ch1)
           IF ( Nx1 + nxc /= nx3 ) cycle
           IF ( Ny1 + nyc /= ny3 ) cycle
           IF ( Nz1 + nzc /= nz3 ) cycle
           IF ( 2*Tz1 + tzc /= tz3 ) cycle
           bra_confs0 = number_2b(3,ch1)
           DO bra0 = 1, bra_confs0
              a = lookup_2b_configs(3,ch1)%ival2(1,bra0)
              b = lookup_2b_configs(3,ch1)%ival2(2,bra0)
              IF ( cut3b(a,b,c) ) cycle
              ii = ii + 1
              IF ( c < a ) reduced_index(ii) = 1
           end DO
           exit
        end DO
     end DO
     ! !!!!!!!!!!!!!!!!!!!!!
     
     has_ch_been_added = .FALSE.
     ii = 0
     DO WHILE ( ii < bra_confs )
        ii = ii + 1
        temp = curr_work + int(ket_confs * (1 + mix1*reduced_index(ii)),8)
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( .not. has_ch_been_added ) THEN
              has_ch_been_added = .TRUE.
              mapping0(curr_proc+1, ch3) = 1
           end IF
        ELSE
           has_ch_been_added = .FALSE.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF
     end DO
     DEALLOCATE( reduced_index )
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE( procs_per_channel(1:channels_paeom3%number_confs) )
  procs_per_channel = 0
  DO ii = 1, num_procs
     DO ch3 = 1, channels_paeom3%number_confs
        IF ( mapping0(ii,ch3) == 0 ) cycle
        procs_per_channel(ch3) = procs_per_channel(ch3) + 1
     end DO
  end DO
  
  count_groups = 0
  DO ch3 = 1, channels_paeom3%number_confs
     IF ( procs_per_channel(ch3) > 1 ) count_groups = count_groups + 1
  end DO
  
  CALL mpi_comm_group(mpi_comm_world,group_world,ierror)
  CALL mpi_barrier(mpi_comm_world,ierror)
  ALLOCATE( group_paeom3(num_procs, channels_paeom3%number_confs) )
  ALLOCATE( subcomm_paeom3(count_groups) )
  group_paeom3 = 0
  subcomm_paeom3 = 0
  CALL mem_register('paeom3', REAL(channels_paeom3%number_confs * num_procs * 4.d0, dp))
  CALL mem_register('paeom3', REAL(count_groups * 4.d0, dp))
  ALLOCATE( group(count_groups) )
  ALLOCATE( rank(count_groups) )
  group = 0
  rank  = 0
  
  count_groups = 0
  DO ch3 = 1, channels_paeom3%number_confs
     IF ( procs_per_channel(ch3) == 1 ) cycle
     count_groups = count_groups + 1
     procs_per_channel0 = procs_per_channel(ch3)
     ALLOCATE( members(procs_per_channel0) )
     
     procs_per_channel0 = 0
     DO ii = 1, num_procs
        IF ( mapping0(ii, ch3) == 0 ) cycle
        procs_per_channel0 = procs_per_channel0 + 1
        members(procs_per_channel0) = ii-1
        group_paeom3(ii, ch3) = count_groups
     end DO
     CALL mpi_group_incl(group_world, procs_per_channel0, members(:), group(count_groups), ierror)
     CALL mpi_comm_create(mpi_comm_world, group(count_groups), subcomm_paeom3(count_groups), ierror)
     
     DEALLOCATE( members )
  end DO
  DEALLOCATE( group, rank, procs_per_channel )
  CALL mpi_barrier(mpi_comm_world,ierror)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  ch3_paeom_min = channels_paeom3%number_confs
  ch3_paeom_max = 0
  DO ii = 1, num_procs
     DO ch3 = 1, channels_paeom3%number_confs
        IF ( mapping0(iam+1,ch3) == 0 ) cycle
        IF ( ch3 < ch3_paeom_min ) ch3_paeom_min = ch3
        IF ( ch3 > ch3_paeom_max ) ch3_paeom_max = ch3
     end DO
  end DO
  DEALLOCATE( mapping0 )
  CALL mpi_barrier(mpi_comm_world,ierror)

  
  ALLOCATE( clist_paeom3(ch3_paeom_min:ch3_paeom_max) )
  ALLOCATE( climit_paeom3(ch3_paeom_min:ch3_paeom_max) ) ! max c
  CALL mem_register('paeom3', REAL((ch3_paeom_max-ch3_paeom_min+1) * 4.d0, dp))
  DO ch3 = ch3_paeom_min, ch3_paeom_max
     nx3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_paeom3%config_NxNyNz_Tz(ch3*4)
     
     c_count = 0
     DO c = below_ef+1, tot_orbs
        nxc = all_orbit%nx(c)
        nyc = all_orbit%ny(c)
        nzc = all_orbit%nz(c)
        tzc = all_orbit%tz(c)
        DO ch1 = 1, channels_2b%number_confs
           IF ( number_2b(3,ch1) == 0 ) cycle
           Nx1 = channels_2b%config_NxNyNz_Tz(4*ch1-3)
           Ny1 = channels_2b%config_NxNyNz_Tz(4*ch1-2)
           Nz1 = channels_2b%config_NxNyNz_Tz(4*ch1-1)
           Tz1 = channels_2b%config_NxNyNz_Tz(4*ch1)
           IF ( Nx1 + nxc /= nx3 ) cycle
           IF ( Ny1 + nyc /= ny3 ) cycle
           IF ( Nz1 + nzc /= nz3 ) cycle
           IF ( 2*Tz1 + tzc /= tz3 ) cycle
           c_count = c_count + 1
           EXIT
        end DO
     end DO
     climit_paeom3(ch3) = c_count
     ALLOCATE( clist_paeom3(ch3)%ival2(c_count,2) )
     clist_paeom3(ch3)%ival2 = 0
     CALL mem_register('paeom3', REAL(2 * c_count * 4.d0, dp))
     
     c_count = 0
     DO c = below_ef+1, tot_orbs
        nxc = all_orbit%nx(c)
        nyc = all_orbit%ny(c)
        nzc = all_orbit%nz(c)
        tzc = all_orbit%tz(c)
        DO ch1 = 1, channels_2b%number_confs
           IF ( number_2b(3,ch1) == 0 ) cycle
           Nx1 = channels_2b%config_NxNyNz_Tz(4*ch1-3)
           Ny1 = channels_2b%config_NxNyNz_Tz(4*ch1-2)
           Nz1 = channels_2b%config_NxNyNz_Tz(4*ch1-1)
           Tz1 = channels_2b%config_NxNyNz_Tz(4*ch1)
           IF ( Nx1 + nxc /= nx3 ) cycle
           IF ( Ny1 + nyc /= ny3 ) cycle
           IF ( Nz1 + nzc /= nz3 ) cycle
           IF ( 2*Tz1 + tzc /= tz3 ) cycle
           c_count = c_count + 1
           clist_paeom3(ch3)%ival2(c_count, 1) = c
           clist_paeom3(ch3)%ival2(c_count, 2) = ch1
           EXIT
        end DO
     end DO
     
  end DO

  ALLOCATE( climits_paeom3(ch3_paeom_min:ch3_paeom_max, 1:2) )
  CALL mem_register('paeom3', REAL(2 * (ch3_paeom_max-ch3_paeom_min+1) * 4.d0, dp))
  ALLOCATE( mapping_paeom3(ch3_paeom_min:ch3_paeom_max) )
  curr_work = eom_ndim2 ! first proc has all r_p and r_pph amplitudes
  curr_proc = 0
  num_chan  = 0
  DO ch3 = 1, channels_paeom3%number_confs
     nx3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_paeom3%config_NxNyNz_Tz(ch3*4)

     IF ( ch3_paeom_min <= ch3 .and. ch3 <= ch3_paeom_max ) then
        climits_paeom3(ch3, 1) = climit_paeom3(ch3)+1
        climits_paeom3(ch3, 2) = 0
        ALLOCATE( mapping_paeom3(ch3)%ival2(climit_paeom3(ch3), 1:2) )
        mapping_paeom3(ch3)%ival2 = 0
        CALL mem_register('paeom3', REAL(2 * climit_paeom3(ch3) * 4.d0, dp))
     end IF
     ch2 = ch2_paeom3(ch3)
     ket_confs = number_2b(1,ch2)
     
     cind1 = 0
     DO c = below_ef+1, tot_orbs
        nxc = all_orbit%nx(c)
        nyc = all_orbit%ny(c)
        nzc = all_orbit%nz(c)
        tzc = all_orbit%tz(c)

        ch1 = 0
        DO ch = 1, channels_2b%number_confs
           IF ( number_2b(3,ch) == 0 ) cycle
           Nx1 = channels_2b%config_NxNyNz_Tz(4*ch-3)
           Ny1 = channels_2b%config_NxNyNz_Tz(4*ch-2)
           Nz1 = channels_2b%config_NxNyNz_Tz(4*ch-1)
           Tz1 = channels_2b%config_NxNyNz_Tz(4*ch)
           IF ( Nx1 + nxc /= nx3 ) cycle
           IF ( Ny1 + nyc /= ny3 ) cycle
           IF ( Nz1 + nzc /= nz3 ) cycle
           IF ( 2*Tz1 + tzc /= tz3 ) cycle
           ch1 = ch
           cind1 = cind1 + 1
           EXIT
        end DO
        IF ( ch1 == 0 ) cycle
        bra_confs = number_2b(3,ch1)
        
        has_ch_been_added = .FALSE.        
        ii = 0
        DO WHILE ( ii < bra_confs )
           ii = ii + 1
           temp = curr_work + int(ket_confs,8)
           IF ( c < lookup_2b_configs(3,ch1)%ival2(1,ii) ) temp = temp + int(mix1 * ket_confs,8)
           
           IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) then
              curr_work = temp
              IF ( has_ch_been_added ) then
                 IF ( curr_proc == iam ) then
                    mapping_paeom3(ch3)%ival2(cind1,2) = ii
                    IF ( cind1 > climits_paeom3(ch3, 2) ) climits_paeom3(ch3, 2) = cind1
                 end IF
              ELSE
                 has_ch_been_added = .TRUE.
                 IF ( curr_proc == iam ) then
                    mapping_paeom3(ch3)%ival2(cind1,1) = ii
                    mapping_paeom3(ch3)%ival2(cind1,2) = ii
                    IF ( cind1 < climits_paeom3(ch3, 1) ) climits_paeom3(ch3, 1) = cind1
                    IF ( cind1 > climits_paeom3(ch3, 2) ) climits_paeom3(ch3, 2) = cind1
                 end IF
              end IF
           ELSE
              has_ch_been_added = .false.
              curr_work = 0
              curr_proc = curr_proc + 1
              ii = ii - 1
           end IF
        end DO
        
     end DO
  end DO  
  
END SUBROUTINE setup_proc_mappings_paeom3


SUBROUTINE setup_proc_mappings_paeomvec
  use parallel
  use configurations
  use operator_storage
  use constants
  USE mem_tracker
  implicit none

  INTEGER :: curr_proc
  INTEGER :: ch,ch3, ch1,ch2, bra, bra_confs,ket_confs
  INTEGER :: nx3,ny3,nz3,tz3, nxc,nyc,nzc,tzc, Nx1,Ny1,Nz1,Tz1
  INTEGER :: a,b, c,c_count, bra_min,bra_max
  INTEGER(i8) :: ltot_work, work_pr_proc, curr_work, ii
  
  ALLOCATE( eom%all_starts(num_procs) )
  ALLOCATE( eom%all_stops(num_procs) )
  eom%all_starts = eom_ndim
  eom%all_stops = 0
  CALL mem_register('paeom2', REAL(2 * num_procs * 8.d0, dp))
  
  IF ( eom_approx == 0 ) then
     
     ltot_work = eom_ndim2
     work_pr_proc = int( ltot_work/int(num_procs,8), 8 )
     
     ! r1_a
     curr_proc = 0
     curr_work = 1
     ii = 1
     eom%all_starts(curr_proc+1) = min(ii, eom%all_starts(curr_proc+1))
     eom%all_stops(curr_proc+1) = max(ii, eom%all_stops(curr_proc+1))
     curr_work = curr_work + 1
     ii = ii + 1
     eom%all_starts(curr_proc+1) = min(ii, eom%all_starts(curr_proc+1))
     eom%all_stops(curr_proc+1) = max(ii, eom%all_stops(curr_proc+1))
     
     ! r2_abi
     DO ch = 1, channels_2b%number_confs
        IF ( number_2b(3,ch) == 0 ) cycle
        IF ( r2_paeom_ind(2*ch) == 0 ) cycle
        bra_confs = number_2b(3,ch)
        
        bra = 0
        DO WHILE ( bra < bra_confs )
           curr_work = curr_work + 2
           bra = bra + 1
           ii  = ii  + 2
           IF ( curr_work <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
              eom%all_starts(curr_proc+1) = min(ii-1, eom%all_starts(curr_proc+1))
              eom%all_stops(curr_proc+1) = max(ii, eom%all_stops(curr_proc+1))
           ELSE
              curr_proc = curr_proc + 1
              curr_work = 0
              bra = bra - 1
              ii  = ii  - 2
           end IF
        end DO
     end DO
     
  ELSE

     ! first proc has all r1 and r2
     eom%all_starts(1) = 1
     eom%all_stops(1) = eom_ndim2
     ii = eom_ndim2
     
     ! CCDT STUFF     
     DO ch3 = 1, channels_paeom3%number_confs
        nx3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-3)
        ny3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-2)
        nz3 = channels_paeom3%config_NxNyNz_Tz(ch3*4-1)
        tz3 = channels_paeom3%config_NxNyNz_Tz(ch3*4)
        
        ch2 = ch2_paeom3(ch3)
        ket_confs = number_2b(1,ch2)

        c_count = 0
        DO c = below_ef+1, tot_orbs
           nxc = all_orbit%nx(c)
           nyc = all_orbit%ny(c)
           nzc = all_orbit%nz(c)
           tzc = all_orbit%tz(c)
           
           DO ch1 = 1, channels_2b%number_confs
              IF ( number_2b(3,ch1) == 0 ) cycle
              Nx1 = channels_2b%config_NxNyNz_Tz(4*ch1-3)
              Ny1 = channels_2b%config_NxNyNz_Tz(4*ch1-2)
              Nz1 = channels_2b%config_NxNyNz_Tz(4*ch1-1)
              Tz1 = channels_2b%config_NxNyNz_Tz(4*ch1)
              IF ( Nx1 + nxc /= nx3 ) cycle
              IF ( Ny1 + nyc /= ny3 ) cycle
              IF ( Nz1 + nzc /= nz3 ) cycle
              IF ( 2*Tz1 + tzc /= tz3 ) cycle
              c_count = c_count + 1
              
              bra_confs = number_2b(3,ch1)
              DO bra = 1, bra_confs
                 a   = lookup_2b_configs(3,ch1)%ival2(1,bra)
                 b   = lookup_2b_configs(3,ch1)%ival2(2,bra)
                 IF ( c >= a ) cycle
                 
                 IF ( ch3_paeom_min <= ch3 .and. ch3 <= ch3_paeom_max ) then
                    IF ( climits_paeom3(ch3, 1) <= c_count .and. c_count <= climits_paeom3(ch3, 2) ) then
                       bra_min = mapping_paeom3(ch3)%ival2(c_count,1)
                       bra_max = mapping_paeom3(ch3)%ival2(c_count,2)
                       IF ( bra_min <= bra .and. bra <= bra_max ) then
                          eom%all_starts(iam+1) = min(ii+1, eom%all_starts(iam+1))
                       end IF
                    end IF
                 end IF
                 
                 ii = ii + ket_confs
                 
                 IF ( ch3_paeom_min <= ch3 .and. ch3 <= ch3_paeom_max ) then
                    IF ( climits_paeom3(ch3, 1) <= c_count .and. c_count <= climits_paeom3(ch3, 2) ) then
                       bra_min = mapping_paeom3(ch3)%ival2(c_count,1)
                       bra_max = mapping_paeom3(ch3)%ival2(c_count,2)
                       IF ( bra_min <= bra .and. bra <= bra_max ) then
                          eom%all_stops(iam+1) = max(ii, eom%all_stops(iam+1))
                       end IF
                    end IF
                 end IF
              end DO
              EXIT
              
           end DO
        end DO
     end DO
     
  end IF
  
  eom%my_start = eom%all_starts(iam + 1)
  eom%my_stop  = eom%all_stops(iam + 1)

END SUBROUTINE setup_proc_mappings_paeomvec



SUBROUTINE setup_proc_mappings_preom
  use parallel
  use configurations
  use operator_storage
  use constants
  USE mem_tracker
  implicit none
  
  INTEGER(i8) :: ltot_work
  INTEGER(i8) :: work_pr_proc, curr_work, temp
  INTEGER :: curr_proc
  INTEGER :: ch, ch2, ii, ket_confs
  LOGICAL :: has_ch_been_added
  
  ! phh
  ltot_work = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     ket_confs = number_2b(1,ch)  
     ltot_work = ltot_work + 2*int(ket_confs,8)
  end DO
  
  work_pr_proc = int( ltot_work/int(num_procs,8) )
  ALLOCATE( mapping_preom_phh(1:num_procs, channels_2b%number_confs, 1:3) )
  mapping_preom_phh = 0
  
  curr_work = 0
  curr_proc = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     ket_confs = number_2b(1,ch)
     has_ch_been_added = .FALSE.
     ii = 0
     DO WHILE ( ii < ket_confs )
        ii = ii + 1
        temp = curr_work + 2
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( has_ch_been_added ) THEN
              mapping_preom_phh(curr_proc+1, ch, 3) = ii
           ELSE
              has_ch_been_added = .true.
              mapping_preom_phh(curr_proc+1, ch, 1) = ch
              mapping_preom_phh(curr_proc+1, ch, 2) = ii
              mapping_preom_phh(curr_proc+1, ch, 3) = ii
           end IF
        ELSE
           has_ch_been_added = .false.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF
     end DO
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  ! phh
  ALLOCATE( check_my_channel_preom_phh(channels_2b%number_confs) )
  check_my_channel_preom_phh = 0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     IF ( mapping_preom_phh(iam+1,ch,1) == 0 ) cycle
     check_my_channel_preom_phh(ch) = 1
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  CALL mem_register('preom2', REAL(channels_2b%number_confs * num_procs * 3 * 4.d0, dp))
  CALL mem_register('preom2', REAL(channels_2b%number_confs * 4.d0, dp))


  ! hhp
  ltot_work = 0
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_preom_ind_cross(2*ch2) == 0 ) cycle
     ket_confs = number_2bcross(2,ch2)
     ltot_work = ltot_work + 2*int(ket_confs,8)
  end DO
  
  work_pr_proc = int( ltot_work/int(num_procs,8) )
  ALLOCATE( mapping_preom_hhp_cross(1:num_procs, channels_2bcross%number_confs, 1:3) )
  mapping_preom_hhp_cross = 0
  
  curr_work = 0
  curr_proc = 0
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_preom_ind_cross(2*ch2) == 0 ) cycle
     ket_confs = number_2bcross(2,ch2)
     has_ch_been_added = .FALSE.        
     ii = 0
     DO WHILE ( ii < ket_confs )
        ii = ii + 1
        temp = curr_work + 2
        
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( has_ch_been_added ) THEN
              mapping_preom_hhp_cross(curr_proc+1, ch2, 3) = ii
           ELSE
              has_ch_been_added = .true.
              mapping_preom_hhp_cross(curr_proc+1, ch2, 1) = ch2
              mapping_preom_hhp_cross(curr_proc+1, ch2, 2) = ii
              mapping_preom_hhp_cross(curr_proc+1, ch2, 3) = ii
           end IF
        ELSE
           has_ch_been_added = .false.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF
     end DO
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  ! hhp
  ALLOCATE( check_my_channel_preom_hhp_cross(channels_2bcross%number_confs) )
  check_my_channel_preom_hhp_cross = 0
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_preom_ind_cross(2*ch2) == 0 ) cycle
     IF ( mapping_preom_hhp_cross(iam+1,ch2,1) == 0 ) cycle
     check_my_channel_preom_hhp_cross(ch2) = 1
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)

  CALL mem_register('preom2', REAL(channels_2bcross%number_confs * num_procs * 3 * 4.d0, dp))
  CALL mem_register('preom2', REAL(channels_2bcross%number_confs * 4.d0, dp))
  
END SUBROUTINE setup_proc_mappings_preom


SUBROUTINE setup_proc_mappings_preom3
  use parallel
  use constants
  use configurations
  use operator_storage
  USE mem_tracker
  
  IMPLICIT NONE
  INTEGER(i8) :: ltot_work,ltot_work1, work_pr_proc, curr_work, temp
  INTEGER :: curr_proc, num_chan, ii
  INTEGER :: ch3, ch,ch1,ch2, bra_confs, ket_confs, ket0,ket_confs0,i
  INTEGER :: nx3,ny3,nz3,tz3, Nx2,Ny2,Nz2,Tz2
  INTEGER :: k, k_count, kind1, nxk,nyk,nzk,tzk
  INTEGER :: group_world, count_groups, procs_per_channel0
  INTEGER, allocatable :: mapping0(:,:)
  INTEGER, allocatable :: procs_per_channel(:)
  INTEGER, allocatable :: reduced_index(:)
  LOGICAL :: has_ch_been_added
  REAL(dp) :: mix1
  
  ! pphhh
  ltot_work = eom_ndim2 ! first proc has all r_h and r_phh amplitudes
  ltot_work1 = eom_ndim2 ! first proc has all r_h and r_phh amplitudes
  mix1 = 1.5d0 ! use for load balancing (0 = balanced r3, ~5 = balanced r3_vec?)
  DO ch3 = 1, channels_preom3%number_confs
     nx3 = channels_preom3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_preom3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_preom3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_preom3%config_NxNyNz_Tz(ch3*4)
     
     ch1 = ch1_preom3(ch3)
     bra_confs = number_2b(3,ch1)
     
     CALL number_3b_confs(nx3, ny3, nz3, tz3, 2, ket_confs) ! 1 for hhh
     ltot_work = ltot_work + int(bra_confs,8)*int(ket_confs,8)
     CALL number_3b_red_confs(nx3, ny3, nz3, tz3, 2, ket_confs) ! 1 for hhh
     ltot_work1 = ltot_work1 + int(bra_confs,8)*int(ket_confs,8)
  end DO

  work_pr_proc = int( (ltot_work + int(mix1*ltot_work1,8)) /num_procs, 8)
  ALLOCATE( mapping0(1:num_procs, 1:channels_preom3%number_confs) )
  mapping0 = 0
  
  curr_work = eom_ndim2 ! first proc has all r_h and r_phh amplitudes
  curr_proc = 0
  DO ch3 = 1, channels_preom3%number_confs
     nx3 = channels_preom3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_preom3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_preom3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_preom3%config_NxNyNz_Tz(ch3*4)
     
     CALL number_3b_confs(nx3, ny3, nz3, tz3, 2, ket_confs) ! 1 for hhh
     ch1 = ch1_preom3(ch3)
     bra_confs = number_2b(3,ch1)

     ! !!!!!!!!!!!!!!!!!!!!!
     ALLOCATE( reduced_index(ket_confs) )
     reduced_index = 0
     ii = 0
     DO k = 1, below_ef
        nxk = all_orbit%nx(k)
        nyk = all_orbit%ny(k)
        nzk = all_orbit%nz(k)
        tzk = all_orbit%tz(k)
        DO ch2 = 1, channels_2b%number_confs
           IF ( number_2b(1,ch2) == 0 ) cycle
           Nx2 = channels_2b%config_NxNyNz_Tz(4*ch2-3)
           Ny2 = channels_2b%config_NxNyNz_Tz(4*ch2-2)
           Nz2 = channels_2b%config_NxNyNz_Tz(4*ch2-1)
           Tz2 = channels_2b%config_NxNyNz_Tz(4*ch2)
           IF ( Nx2 + nxk /= nx3 ) cycle
           IF ( Ny2 + nyk /= ny3 ) cycle
           IF ( Nz2 + nzk /= nz3 ) cycle
           IF ( 2*Tz2 + tzk /= tz3 ) cycle
           ket_confs0 = number_2b(1,ch2)
           DO ket0 = 1, ket_confs0
              ii = ii + 1
              i = lookup_2b_configs(1,ch2)%ival2(1,ket0)
              IF ( k < i ) reduced_index(ii) = 1
           end DO
           exit
        end DO
     end DO
     ! !!!!!!!!!!!!!!!!!!!!!
     
     has_ch_been_added = .FALSE.
     ii = 0
     DO WHILE ( ii < ket_confs )
        ii = ii + 1
        temp = curr_work + int(bra_confs * (1 + mix1*reduced_index(ii)),8)
        IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
           curr_work = temp
           IF ( .not. has_ch_been_added ) THEN
              has_ch_been_added = .TRUE.
              mapping0(curr_proc+1, ch3) = 1
           end IF
        ELSE
           has_ch_been_added = .FALSE.
           curr_work = 0
           curr_proc = curr_proc + 1
           ii = ii - 1
        end IF
     end DO
     DEALLOCATE( reduced_index )
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE( procs_per_channel(1:channels_preom3%number_confs) )
  procs_per_channel = 0
  DO ii = 1, num_procs
     DO ch3 = 1, channels_preom3%number_confs
        IF ( mapping0(ii,ch3) == 0 ) cycle
        procs_per_channel(ch3) = procs_per_channel(ch3) + 1
     end DO
  end DO
  
  count_groups = 0
  DO ch3 = 1, channels_preom3%number_confs
     IF ( procs_per_channel(ch3) > 1 ) count_groups = count_groups + 1
  end DO
  
  CALL mpi_comm_group(mpi_comm_world,group_world,ierror)
  CALL mpi_barrier(mpi_comm_world,ierror)
  ALLOCATE( group_preom3(num_procs, channels_preom3%number_confs) )
  ALLOCATE( subcomm_preom3(count_groups) )
  group_preom3 = 0
  subcomm_preom3 = 0
  CALL mem_register('preom3', REAL(channels_preom3%number_confs * num_procs * 4.d0, dp))
  CALL mem_register('preom3', REAL(count_groups * 4.d0, dp))
  ALLOCATE( group(count_groups) )
  ALLOCATE( rank(count_groups) )
  group = 0
  rank  = 0
  
  count_groups = 0
  DO ch3 = 1, channels_preom3%number_confs
     IF ( procs_per_channel(ch3) == 1 ) cycle
     count_groups = count_groups + 1
     procs_per_channel0 = procs_per_channel(ch3)
     ALLOCATE( members(procs_per_channel0) )
     
     procs_per_channel0 = 0
     DO ii = 1, num_procs
        IF ( mapping0(ii, ch3) == 0 ) cycle
        procs_per_channel0 = procs_per_channel0 + 1
        members(procs_per_channel0) = ii-1
        group_preom3(ii, ch3) = count_groups
     end DO
     CALL mpi_group_incl(group_world, procs_per_channel0, members(:), group(count_groups), ierror)
     CALL mpi_comm_create(mpi_comm_world, group(count_groups), subcomm_preom3(count_groups), ierror)
     
     DEALLOCATE( members )
  end DO
  DEALLOCATE( group, rank, procs_per_channel )
  CALL mpi_barrier(mpi_comm_world,ierror)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  ch3_preom_min = channels_preom3%number_confs
  ch3_preom_max = 0
  DO ii = 1, num_procs
     DO ch3 = 1, channels_preom3%number_confs
        IF ( mapping0(iam+1,ch3) == 0 ) cycle
        IF ( ch3 < ch3_preom_min ) ch3_preom_min = ch3
        IF ( ch3 > ch3_preom_max ) ch3_preom_max = ch3
     end DO
  end DO
  DEALLOCATE( mapping0 )
  CALL mpi_barrier(mpi_comm_world,ierror)

  
  ALLOCATE( klist_preom3(ch3_preom_min:ch3_preom_max) )
  ALLOCATE( klimit_preom3(ch3_preom_min:ch3_preom_max) ) ! max k
  CALL mem_register('preom3', REAL((ch3_preom_max-ch3_preom_min+1) * 4.d0, dp))
  DO ch3 = ch3_preom_min, ch3_preom_max
     nx3 = channels_preom3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_preom3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_preom3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_preom3%config_NxNyNz_Tz(ch3*4)
     
     k_count = 0
     DO k = 1, below_ef
        nxk = all_orbit%nx(k)
        nyk = all_orbit%ny(k)
        nzk = all_orbit%nz(k)
        tzk = all_orbit%tz(k)
        DO ch2 = 1, channels_2b%number_confs
           IF ( number_2b(1,ch2) == 0 ) cycle
           Nx2 = channels_2b%config_NxNyNz_Tz(4*ch2-3)
           Ny2 = channels_2b%config_NxNyNz_Tz(4*ch2-2)
           Nz2 = channels_2b%config_NxNyNz_Tz(4*ch2-1)
           Tz2 = channels_2b%config_NxNyNz_Tz(4*ch2)
           IF ( Nx2 + nxk /= nx3 ) cycle
           IF ( Ny2 + nyk /= ny3 ) cycle
           IF ( Nz2 + nzk /= nz3 ) cycle
           IF ( 2*Tz2 + tzk /= tz3 ) cycle
           k_count = k_count + 1
           EXIT
        end DO
     end DO
     klimit_preom3(ch3) = k_count
     ALLOCATE( klist_preom3(ch3)%ival2(k_count,2) )
     klist_preom3(ch3)%ival2 = 0
     CALL mem_register('preom3', REAL(2 * k_count * 4.d0, dp))
     
     k_count = 0
     DO k = 1, below_ef
        nxk = all_orbit%nx(k)
        nyk = all_orbit%ny(k)
        nzk = all_orbit%nz(k)
        tzk = all_orbit%tz(k)
        DO ch2 = 1, channels_2b%number_confs
           IF ( number_2b(1,ch2) == 0 ) cycle
           Nx2 = channels_2b%config_NxNyNz_Tz(4*ch2-3)
           Ny2 = channels_2b%config_NxNyNz_Tz(4*ch2-2)
           Nz2 = channels_2b%config_NxNyNz_Tz(4*ch2-1)
           Tz2 = channels_2b%config_NxNyNz_Tz(4*ch2)
           IF ( Nx2 + nxk /= nx3 ) cycle
           IF ( Ny2 + nyk /= ny3 ) cycle
           IF ( Nz2 + nzk /= nz3 ) cycle
           IF ( 2*Tz2 + tzk /= tz3 ) cycle
           k_count = k_count + 1
           klist_preom3(ch3)%ival2(k_count, 1) = k
           klist_preom3(ch3)%ival2(k_count, 2) = ch2
           EXIT
        end DO
     end DO
     
  end DO

  ALLOCATE( klimits_preom3(ch3_preom_min:ch3_preom_max, 1:2) )
  CALL mem_register('preom3', REAL(2 * (ch3_preom_max-ch3_preom_min+1) * 4.d0, dp))
  ALLOCATE( mapping_preom3(ch3_preom_min:ch3_preom_max) )
  curr_work = eom_ndim2 ! first proc has all r_h and r_phh amplitudes
  curr_proc = 0
  num_chan  = 0
  DO ch3 = 1, channels_preom3%number_confs
     nx3 = channels_preom3%config_NxNyNz_Tz(ch3*4-3)
     ny3 = channels_preom3%config_NxNyNz_Tz(ch3*4-2)
     nz3 = channels_preom3%config_NxNyNz_Tz(ch3*4-1)
     tz3 = channels_preom3%config_NxNyNz_Tz(ch3*4)

     IF ( ch3_preom_min <= ch3 .and. ch3 <= ch3_preom_max ) then
        klimits_preom3(ch3, 1) = klimit_preom3(ch3)+1
        klimits_preom3(ch3, 2) = 0
        ALLOCATE( mapping_preom3(ch3)%ival2(klimit_preom3(ch3), 1:2) )
        mapping_preom3(ch3)%ival2 = 0
        CALL mem_register('preom3', REAL(2 * klimit_preom3(ch3) * 4.d0, dp))
     end IF
     ch1 = ch1_preom3(ch3)
     bra_confs = number_2b(3,ch1)
     
     kind1 = 0
     DO k = 1, below_ef
        nxk = all_orbit%nx(k)
        nyk = all_orbit%ny(k)
        nzk = all_orbit%nz(k)
        tzk = all_orbit%tz(k)

        ch2 = 0
        DO ch = 1, channels_2b%number_confs
           IF ( number_2b(1,ch) == 0 ) cycle
           Nx2 = channels_2b%config_NxNyNz_Tz(4*ch-3)
           Ny2 = channels_2b%config_NxNyNz_Tz(4*ch-2)
           Nz2 = channels_2b%config_NxNyNz_Tz(4*ch-1)
           Tz2 = channels_2b%config_NxNyNz_Tz(4*ch)
           IF ( Nx2 + nxk /= nx3 ) cycle
           IF ( Ny2 + nyk /= ny3 ) cycle
           IF ( Nz2 + nzk /= nz3 ) cycle
           IF ( 2*Tz2 + tzk /= tz3 ) cycle
           ch2 = ch
           kind1 = kind1 + 1
           EXIT
        end DO
        IF ( ch2 == 0 ) cycle
        ket_confs = number_2b(1,ch2)
        
        has_ch_been_added = .FALSE.        
        ii = 0
        DO WHILE ( ii < ket_confs )
           ii = ii + 1
           temp = curr_work + int(bra_confs,8)
           IF ( k < lookup_2b_configs(1,ch2)%ival2(1,ii) ) temp = temp + int(mix1 * bra_confs,8)
           
           IF ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) then
              curr_work = temp
              IF ( has_ch_been_added ) then
                 IF ( curr_proc == iam ) then
                    mapping_preom3(ch3)%ival2(kind1,2) = ii
                    IF ( kind1 > klimits_preom3(ch3, 2) ) klimits_preom3(ch3, 2) = kind1
                 end IF
              ELSE
                 has_ch_been_added = .TRUE.
                 IF ( curr_proc == iam ) then
                    mapping_preom3(ch3)%ival2(kind1,1) = ii
                    mapping_preom3(ch3)%ival2(kind1,2) = ii
                    IF ( kind1 < klimits_preom3(ch3, 1) ) klimits_preom3(ch3, 1) = kind1
                    IF ( kind1 > klimits_preom3(ch3, 2) ) klimits_preom3(ch3, 2) = kind1
                 end IF
              end IF
           ELSE
              has_ch_been_added = .false.
              curr_work = 0
              curr_proc = curr_proc + 1
              ii = ii - 1
           end IF
        end DO
        
     end DO
  end DO
  
END SUBROUTINE setup_proc_mappings_preom3


SUBROUTINE setup_proc_mappings_preomvec
  use parallel
  use configurations
  use operator_storage
  use constants
  USE mem_tracker
  implicit none

  INTEGER :: curr_proc
  INTEGER :: ch,ch3, ch1,ch2, ket, bra_confs,ket_confs
  INTEGER :: nx3,ny3,nz3,tz3, nxk,nyk,nzk,tzk, Nx2,Ny2,Nz2,Tz2
  INTEGER :: i,j, k,k_count, ket_min,ket_max
  INTEGER(i8) :: ltot_work, work_pr_proc, curr_work, ii
  
  ALLOCATE( eom%all_starts(num_procs) )
  ALLOCATE( eom%all_stops(num_procs) )
  eom%all_starts = eom_ndim
  eom%all_stops = 0
  CALL mem_register('preom2', REAL(2 * num_procs * 8.d0, dp))
  
  IF ( eom_approx == 0 ) then
     
     ltot_work = eom_ndim2
     work_pr_proc = int( ltot_work/int(num_procs,8), 8 )
     
     ! r1_i
     curr_proc = 0
     curr_work = 1
     ii = 1
     eom%all_starts(curr_proc+1) = min(ii, eom%all_starts(curr_proc+1))
     eom%all_stops(curr_proc+1) = max(ii, eom%all_stops(curr_proc+1))
     curr_work = curr_work + 1
     ii = ii + 1
     eom%all_starts(curr_proc+1) = min(ii, eom%all_starts(curr_proc+1))
     eom%all_stops(curr_proc+1) = max(ii, eom%all_stops(curr_proc+1))
     
     ! r2_bij
     DO ch = 1, channels_2b%number_confs
        IF ( number_2b(1,ch) == 0 ) cycle
        IF ( r2_preom_ind(2*ch) == 0 ) cycle
        ket_confs = number_2b(1,ch)
        
        ket = 0
        DO WHILE ( ket < ket_confs )
           curr_work = curr_work + 2
           ket = ket + 1
           ii  = ii  + 2
           IF ( curr_work <= work_pr_proc .or. curr_proc == num_procs - 1 ) THEN
              eom%all_starts(curr_proc+1) = min(ii-1, eom%all_starts(curr_proc+1))
              eom%all_stops(curr_proc+1) = max(ii, eom%all_stops(curr_proc+1))
           ELSE
              curr_proc = curr_proc + 1
              curr_work = 0
              ket = ket - 1
              ii  = ii  - 2
           end IF
        end DO
     end DO
     
  ELSE

     ! first proc has all r1 and r2
     eom%all_starts(1) = 1
     eom%all_stops(1) = eom_ndim2
     ii = eom_ndim2
     
     ! CCDT STUFF
     DO ch3 = 1, channels_preom3%number_confs
        nx3 = channels_preom3%config_NxNyNz_Tz(ch3*4-3)
        ny3 = channels_preom3%config_NxNyNz_Tz(ch3*4-2)
        nz3 = channels_preom3%config_NxNyNz_Tz(ch3*4-1)
        tz3 = channels_preom3%config_NxNyNz_Tz(ch3*4)
        
        ch1 = ch1_preom3(ch3)
        bra_confs = number_2b(3,ch1)

        k_count = 0
        DO k = 1, below_ef
           nxk = all_orbit%nx(k)
           nyk = all_orbit%ny(k)
           nzk = all_orbit%nz(k)
           tzk = all_orbit%tz(k)
           
           DO ch2 = 1, channels_2b%number_confs
              IF ( number_2b(1,ch2) == 0 ) cycle
              Nx2 = channels_2b%config_NxNyNz_Tz(4*ch2-3)
              Ny2 = channels_2b%config_NxNyNz_Tz(4*ch2-2)
              Nz2 = channels_2b%config_NxNyNz_Tz(4*ch2-1)
              Tz2 = channels_2b%config_NxNyNz_Tz(4*ch2)
              IF ( Nx2 + nxk /= nx3 ) cycle
              IF ( Ny2 + nyk /= ny3 ) cycle
              IF ( Nz2 + nzk /= nz3 ) cycle
              IF ( 2*Tz2 + tzk /= tz3 ) cycle
              k_count = k_count + 1
              
              ket_confs = number_2b(1,ch2)
              DO ket = 1, ket_confs
                 i   = lookup_2b_configs(1,ch2)%ival2(1,ket)
                 j   = lookup_2b_configs(1,ch2)%ival2(2,ket)
                 IF ( k >= i ) cycle
                 
                 IF ( ch3_preom_min <= ch3 .and. ch3 <= ch3_preom_max ) then
                    IF ( klimits_preom3(ch3, 1) <= k_count .and. k_count <= klimits_preom3(ch3, 2) ) then
                       ket_min = mapping_preom3(ch3)%ival2(k_count,1)
                       ket_max = mapping_preom3(ch3)%ival2(k_count,2)
                       IF ( ket_min <= ket .and. ket <= ket_max ) then
                          eom%all_starts(iam+1) = min(ii+1, eom%all_starts(iam+1))
                       end IF
                    end IF
                 end IF
                 
                 ii = ii + bra_confs
                 
                 IF ( ch3_preom_min <= ch3 .and. ch3 <= ch3_preom_max ) then
                    IF ( klimits_preom3(ch3, 1) <= k_count .and. k_count <= klimits_preom3(ch3, 2) ) then
                       ket_min = mapping_preom3(ch3)%ival2(k_count,1)
                       ket_max = mapping_preom3(ch3)%ival2(k_count,2)
                       IF ( ket_min <= ket .and. ket <= ket_max ) then
                          eom%all_stops(iam+1) = max(ii, eom%all_stops(iam+1))
                       end IF
                    end IF
                 end IF
              end DO
              EXIT
              
           end DO
        end DO
     end DO
     
  end IF
  
  eom%my_start = eom%all_starts(iam + 1)
  eom%my_stop  = eom%all_stops(iam + 1)

END SUBROUTINE setup_proc_mappings_preomvec
