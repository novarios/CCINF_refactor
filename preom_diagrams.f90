
! <|r|i> <-- - <k|x|i>.<|r|k> - (1/2).<kl|x|id>.<d|r|kl>
SUBROUTINE preom_1h
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: ch, bra,ket, ket_min,ket_max
  INTEGER :: i,k, d,d1, iind,kind1
  INTEGER :: dim1,dim2,dim3
  COMPLEX(dpc) :: r1, v1b, x2
  COMPLEX(dpc), allocatable :: temp_mtx(:,:)

  IF ( iam == 0 ) WRITE(6,*) "...Computing Right-1h amplitudes..."

  ! <|r|i>
  r1_preom_eqn = 0.d0

  ! <|r|i> <-- - (1/2).<kl|v|id>.<d|r|kl>
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( number_2b(2,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     
     IF ( check_my_channel_preom_phh(ch) == 0 ) cycle
     ket_min = mapping_preom_phh(iam+1,ch,2)
     ket_max = mapping_preom_phh(iam+1,ch,3)

     ! <d|r|kl>.<kl|v|id>
     dim1 = 2
     dim2 = number_2b(2,ch)
     dim3 = ket_max-ket_min + 1
     ALLOCATE ( temp_mtx(dim1,dim2) )
     temp_mtx = 0.d0
     CALL ZGEMM ( 'n', 't', dim1, dim2, dim3, dcmplx(1.d0,0.d0), r2_preom(ch)%cval(:,ket_min:ket_max), dim1, &
          conjg(v2b_hphh(ch)%cval(:,ket_min:ket_max)), dim2, dcmplx(1.d0,0.d0), temp_mtx, dim1 )
     
     DO ket = 1, dim2
        i   = lookup_2b_configs(2,ch)%ival2(1,ket)
        d   = lookup_2b_configs(2,ch)%ival2(2,ket)
        DO bra = 1, 2
           d1  = r2_preom_ind(2*(ch-1)+bra)
           IF ( d1 /= d ) cycle
           x2 = temp_mtx(bra,ket)
           IF ( i == r1_preom_ind(1) ) r1_preom_eqn(1) = r1_preom_eqn(1) - x2
           IF ( i == r1_preom_ind(2) ) r1_preom_eqn(2) = r1_preom_eqn(2) - x2
        end DO
     end DO
     DEALLOCATE ( temp_mtx )
  end DO
  CALL mpi_allreduce(mpi_in_place,r1_preom_eqn,size(r1_preom_eqn),mpi_complex16,mpi_sum,mpi_comm_world,ierror)

  ! <|r|i> <-- - <k|x|i>.<|r|k>
  DO iind = 1, 2
     i    = r1_preom_ind(iind)
     DO kind1 = 1, 2
        k    = r1_preom_ind(kind1)
        v1b  = hbar1b_I3(k,i)
        r1   = r1_preom(kind1)
        r1_preom_eqn(iind) = r1_preom_eqn(iind) - v1b*r1
     end DO
  end DO
  
end SUBROUTINE preom_1h

! <b|r|ij> <-- - <kb|x|ij>.<|r|k> - P(ij).<k|x|i>.<b|r|kj> + <b|x|c>.<c|r|ij> + (1/2).<kl|x|ij>.<b|r|kl> - P(ij).<kb|x|ic>.<c|r|kj> + <|z3b|c>.<cb|t|ij>
SUBROUTINE preom_1p2h
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: ch,ch2, bra,bra2, ket_confs
  INTEGER :: ket_min,ket_max
  INTEGER :: i,b, c,k, kind1
  INTEGER :: dim1,dim2,dim3
  REAL(dp) :: phase
  COMPLEX(dpc) :: v1b

  IF ( iam == 0 ) WRITE(6,*) "...Computing Right-1p2h amplitudes..."
  
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     ! <b|r|ij>
     r2_preom_eqn(ch)%cval = 0.d0
  end DO
  
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     
     IF ( check_my_channel_preom_phh(ch) == 0 ) cycle
     ket_min = mapping_preom_phh(iam+1,ch,2)
     ket_max = mapping_preom_phh(iam+1,ch,3)
     ket_confs = number_2b(1,ch)

     ! <b|r|ij> <-- + <b|x|c>.<c|r|ij>
     DO bra  = 1, 2
        b    = r2_preom_ind(2*(ch-1)+bra)
        DO bra2 = 1, 2
           c    = r2_preom_ind(2*(ch-1)+bra2)
           v1b  = hbar1b_I2(b,c)
           r2_preom_eqn(ch)%cval(bra,ket_min:ket_max) = r2_preom_eqn(ch)%cval(bra,ket_min:ket_max) &
                + r2_preom(ch)%cval(bra2,ket_min:ket_max)*v1b
        end DO
     end DO

     ! <b|r|ij> <-- + (1/2).<kl|x|ij>.<b|r|kl>
     dim1 = 2
     dim2 = ket_max-ket_min + 1
     dim3 = ket_confs
     CALL ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), r2_preom(ch)%cval, dim1, &
          hbar2b_I4(ch)%cval(:,ket_min:ket_max), dim3, dcmplx(1.d0,0.d0), r2_preom_eqn(ch)%cval(:,ket_min:ket_max), dim1 )

     ! <b|r|ij> <-- - <kb|x|ij>.<|r|k>
     DO bra = 1, 2
        b   = r2_preom_ind(2*(ch-1)+bra)
        DO kind1 = 1, 2
           k    = r1_preom_ind(kind1)
           IF ( ch /= hp_channel_2b%ival2(k,b) ) cycle
           bra2 = hp_config_2b%ival2(k,b)
           r2_preom_eqn(ch)%cval(bra,ket_min:ket_max) = r2_preom_eqn(ch)%cval(bra,ket_min:ket_max) &
                - hbar2b_I7(ch)%cval(bra2,ket_min:ket_max)*r1_preom(kind1)
        end DO
     end DO

     ! <b|r|ij> <-- + <|z3b|c>.<cb|t|ij>
     DO bra = 1, 2
        b   = r2_preom_ind(2*(ch-1)+bra)
        DO c = below_ef+1, tot_orbs
           IF ( ch /= pp_channel_2b%ival2(c,b) ) cycle
           phase = 1.d0
           bra2 = pp_config_2b%ival2(c,b)
           IF ( b < c ) phase = -phase
           r2_preom_eqn(ch)%cval(bra,ket_min:ket_max) = r2_preom_eqn(ch)%cval(bra,ket_min:ket_max) &
                + phase * t2_ccm(ch)%cval(bra2,ket_min:ket_max)*hbar3b_preom_I2(c)
        end DO
     end DO
     
  end DO

  
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_preom_ind_cross(2*ch2) == 0 ) cycle
     ! <-i|r|j-b>
     r2_preom_eqn_cross(ch2)%cval = 0.d0
  end DO

  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_preom_ind_cross(2*ch2) == 0 ) cycle

     IF ( check_my_channel_preom_hhp_cross(ch2) == 0 ) cycle
     ket_min = mapping_preom_hhp_cross(iam+1,ch2,2)
     ket_max = mapping_preom_hhp_cross(iam+1,ch2,3)
     ket_confs = number_2bcross(2,ch2)

     ! <b|r|ij>   <-- - P(ij).<k|x|i>.<b|r|kj>
     ! <-i|r|j-b> <-- - P(ij).<-i|x|-k>.<-k|r|j-b>
     DO bra  = 1, 2
        i    = r2_preom_ind_cross(2*(ch2-1)+bra)
        DO bra2 = 1, 2
           k    = r2_preom_ind_cross(2*(ch2-1)+bra2)
           v1b  = hbar1b_I3(k,i)
           r2_preom_eqn_cross(ch2)%cval(bra,ket_min:ket_max) = r2_preom_eqn_cross(ch2)%cval(bra,ket_min:ket_max) &
                - r2_preom_cross(ch2)%cval(bra2,ket_min:ket_max)*v1b
        end DO
     end DO

     ! <b|r|ij>   <-- - P(ij).<kb|x|ic>.<c|r|kj>
     ! <b|r|ij>   <-- + P(ij).<kb|x|jc>.<c|r|ki>
     ! <b|r|ij>   <-- - P(ij).<kb|x|jc>.<c|r|ik>
     ! <-i|r|j-b> <-- - P(ij).<k-c|x|j-b>.<-i|r|k-c>
     ! <-i|r|j-b> <-- - P(ij).<-i|r|k-c>.<k-c|x|j-b>
     dim1 = 2
     dim2 = ket_max-ket_min + 1
     dim3 = ket_confs
     CALL ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(-1.d0,0.d0), r2_preom_cross(ch2)%cval, &
          dim1, hbar2b_I5_cross(ch2)%cval(:,ket_min:ket_max), dim3, dcmplx(1.d0,0.d0), &
          r2_preom_eqn_cross(ch2)%cval(:,ket_min:ket_max), dim1 )
  end DO
  ! Allreduce r2_preom_eqn in "populate_preom_vec"

end SUBROUTINE preom_1p2h



! <|r|i> <-- + (1/4).<kl|v|cd>.<cd|r|ikl>
SUBROUTINE preom_1h_2p3h
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: ch3,ch1,ch2, i,iind,iind0
  INTEGER :: ket_min,ket_max, bra_confs, bra,ket
  COMPLEX(dpc) :: r3, v2b
  COMPLEX(dpc), allocatable :: temp_mtx(:)
  
  IF ( iam == 0 ) WRITE(6,*) "...Computing R1 <- R3..."

  ! <|r|i> <-- + (1/4).<kl|v|cd>.<cd|r|ikl>
  ALLOCATE( temp_mtx(2) )
  temp_mtx = 0.d0
  DO ch3 = ch3_preom_min, ch3_preom_max
     ch1       = ch1_preom3(ch3)
     bra_confs = number_2b(3,ch1)
     DO iind = klimits_preom3(ch3,1), klimits_preom3(ch3,2)
        i       = klist_preom3(ch3)%ival2(iind,1)
        ch2     = klist_preom3(ch3)%ival2(iind,2)
        IF ( ch1 /= ch2 ) cycle
        IF ( i == r1_preom_ind(1) ) then
           iind0 = 1
        ELSE IF ( i == r1_preom_ind(2) ) then
           iind0 = 2
        ELSE
           cycle
        end IF
        ket_min = mapping_preom3(ch3)%ival2(iind,1)
        ket_max = mapping_preom3(ch3)%ival2(iind,2)

        DO bra = 1, bra_confs
           DO ket = ket_min, ket_max
              r3  = r3_preom(ch3)%val1(iind)%cval(bra,ket)
              v2b = v2b_pphh(ch1)%cval(bra,ket)
              temp_mtx(iind0) = temp_mtx(iind0) + conjg(v2b)*r3
           end DO
        end DO
     end DO
  end DO

  CALL mpi_allreduce(mpi_in_place,temp_mtx,size(temp_mtx),mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  r1_preom_eqn(:) = r1_preom_eqn(:) + temp_mtx(:)
  DEALLOCATE( temp_mtx )
  
end SUBROUTINE preom_1h_2p3h

! <b|r|ij> < -- - (1/2).P(ij).<kl|v|jc>.<bc|r|ikl> - (1/2).<kb|v|cd>.<cd|r|ijk>
SUBROUTINE preom_1p2h_2p3h
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: ch3, ch1,ch2, ch
  INTEGER :: bra_confs, ket_min,ket_max
  INTEGER :: bra,ket, bra0,ket0
  INTEGER :: iind,i, kind1,k, j,b,c
  INTEGER :: dim1,dim2,dim3
  REAL(dp) :: phase, phase0
  COMPLEX(dpc), allocatable :: temp_mtx(:,:)
  
  IF ( iam == 0 ) WRITE(6,*) "...Computing R2 <- R3..."

  ! <b|r|ij> <-- - (1/2).P(ij).<kl|v|jc>.<bc|r|ikl>
  DO ch3 = ch3_preom_min, ch3_preom_max
     ch1       = ch1_preom3(ch3)
     bra_confs = number_2b(3,ch1)
     DO iind = klimits_preom3(ch3,1), klimits_preom3(ch3,2)
        i       = klist_preom3(ch3)%ival2(iind,1)
        ch2     = klist_preom3(ch3)%ival2(iind,2)
        ket_min = mapping_preom3(ch3)%ival2(iind,1)
        ket_max = mapping_preom3(ch3)%ival2(iind,2)

        ! <bc|r|(i)kl>.<kl|v|jc>
        dim1 = bra_confs
        dim2 = number_2b(2,ch2)
        dim3 = ket_max-ket_min+1
        IF ( dim2 == 0 ) cycle
        ALLOCATE ( temp_mtx(dim1,dim2) )
        temp_mtx = dcmplx(0.d0,0.d0)
        CALL ZGEMM ( 'n', 't', dim1, dim2, dim3, dcmplx(1.d0,0.d0), r3_preom(ch3)%val1(iind)%cval(:,ket_min:ket_max), dim1, &
             conjg(v2b_hphh(ch2)%cval(:,ket_min:ket_max)), dim2, dcmplx(1.d0,0.d0), temp_mtx, dim1 )
        
        DO ket0 = 1, dim2
           j    = lookup_2b_configs(2,ch2)%ival2(1,ket0)
           c    = lookup_2b_configs(2,ch2)%ival2(2,ket0)

           ch = hh_channel_2b%ival2(i,j)
           IF ( ch == 0 ) cycle
           IF ( number_2b(1,ch) == 0 ) cycle
           IF ( r2_preom_ind(2*ch) == 0 ) cycle
           phase = 1
           ket = hh_config_2b%ival2(i,j)
           IF ( j < i ) phase = -phase ! P(ij)
           
           DO bra = 1, 2
              b   = r2_preom_ind(2*(ch-1)+bra)
              IF ( ch1 /= pp_channel_2b%ival2(b,c) ) cycle
              phase0 = 1
              bra0 = pp_config_2b%ival2(b,c)
              IF ( c < b ) phase0 = -phase0
              r2_preom_eqn(ch)%cval(bra,ket) = r2_preom_eqn(ch)%cval(bra,ket) - phase * phase0 * temp_mtx(bra0,ket0)
           end DO
        end DO
        DEALLOCATE( temp_mtx )
     end DO
  end DO

  ! <b|r|ij> <-- - (1/2).<kb|v|cd>.<cd|r|ijk>
  !              - (1/2).<kb|v|cd>.<cd|r|kij>
  DO ch3 = ch3_preom_min, ch3_preom_max
     ch1       = ch1_preom3(ch3)
     bra_confs = number_2b(3,ch1)
     DO kind1 = klimits_preom3(ch3,1), klimits_preom3(ch3,2)
        k       = klist_preom3(ch3)%ival2(kind1,1)
        ch2     = klist_preom3(ch3)%ival2(kind1,2)
        ket_min = mapping_preom3(ch3)%ival2(kind1,1)
        ket_max = mapping_preom3(ch3)%ival2(kind1,2)
        IF ( number_2b(1,ch2) == 0 ) cycle
        IF ( r2_preom_ind(2*ch2) == 0 ) cycle
        
        ! <kb|v|cd>.<cd|r|(k)ij>
        dim1 = number_2b(2,ch1)
        dim2 = ket_max-ket_min+1
        dim3 = bra_confs
        IF ( dim1 == 0 ) cycle
        ALLOCATE( temp_mtx(dim1,ket_min:ket_max) )
        temp_mtx = dcmplx(0.d0,0.d0)
        CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), conjg(v2b_pphp(ch1)%cval), dim3, &
             r3_preom(ch3)%val1(kind1)%cval(:,ket_min:ket_max), dim3, dcmplx(1.d0,0.d0), temp_mtx, dim1 )

        DO bra = 1, 2
           b   = r2_preom_ind(2*(ch2-1)+bra)
           IF ( ch1 /= hp_channel_2b%ival2(k,b) ) cycle
           bra0 = hp_config_2b%ival2(k,b)
           r2_preom_eqn(ch2)%cval(bra,ket_min:ket_max) = r2_preom_eqn(ch2)%cval(bra,ket_min:ket_max) &
                - temp_mtx(bra0,ket_min:ket_max)
        end DO
        DEALLOCATE( temp_mtx )
     end DO
  end DO
  
end SUBROUTINE preom_1p2h_2p3h

! <bc|r|ijk> < -- E(bc_ijk).<bc|r|ijk> - P(bc|i/jk).<lc|v|jk>.<b|r|il> - P(i/jk).<bc|v|id>.<d|r|jk>
SUBROUTINE preom_2p3h
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: ch3, ch1,ch2, bra,ket
  INTEGER :: bra_confs, ket_min,ket_max
  INTEGER :: iind, i,j,k,b,c
  INTEGER :: ch_jk,ch_ik,ch_ji, ket_jk,ket_ik,ket_ji
  INTEGER :: phase_jk, phase_ik, phase_ji
  COMPLEX(dpc) :: denom, r3
  
  IF ( iam == 0 ) WRITE(6,*) "...Computing R3..."

  ! <bc|r|ijk> <-- E(bc_ijk).<bc|r|ijk>
  DO ch3 = ch3_preom_min, ch3_preom_max
     ch1       = ch1_preom3(ch3)
     bra_confs = number_2b(3,ch1)
     DO iind = klimits_preom3(ch3,1), klimits_preom3(ch3,2)
        i       = klist_preom3(ch3)%ival2(iind,1)
        ch2     = klist_preom3(ch3)%ival2(iind,2)
        ket_min = mapping_preom3(ch3)%ival2(iind,1)
        ket_max = mapping_preom3(ch3)%ival2(iind,2)
        DO bra = 1, bra_confs
           b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
           c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
           DO ket = ket_min, ket_max
              j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
              k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
              
              r3 = r3_preom(ch3)%val1(iind)%cval(bra,ket)
              r3_preom(ch3)%val1(iind)%cval(bra,ket) = 0.d0
              denom = ( hbar1b_I2(b,b) + hbar1b_I2(c,c) &
                   - hbar1b_I3(i,i) - hbar1b_I3(j,j) - hbar1b_I3(k,k) )
              r3_preom(ch3)%val1(iind)%cval(bra,ket) = r3 * denom
           end DO
        end DO
     end DO
  end DO

  ! <bc|r|ijk> <-- - P(bc|i/jk).<lc|v|jk>.<b|r|il>
  DO ch3 = ch3_preom_min, ch3_preom_max
     ch1       = ch1_preom3(ch3)
     bra_confs = number_2b(3,ch1)
     DO iind = klimits_preom3(ch3,1), klimits_preom3(ch3,2)
        i       = klist_preom3(ch3)%ival2(iind,1)
        ch2     = klist_preom3(ch3)%ival2(iind,2)
        ket_min = mapping_preom3(ch3)%ival2(iind,1)
        ket_max = mapping_preom3(ch3)%ival2(iind,2)

        DO bra = 1, bra_confs
           b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
           c   = lookup_2b_configs(3,ch1)%ival2(2,bra)

           ! <bc|r|ijk> <-- - P(bc).<lc|v|jk>.<b|r|il>
           DO ket = ket_min, ket_max
              j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
              k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
              IF ( j == i .or. k == i ) cycle
              phase_jk = 1
              ket_jk = ket
              ch_jk = ch2
              
              CALL r3_diag3(ch3, iind, bra,ket, ch_jk,ket_jk,phase_jk, i,b,c, 1)
              CALL r3_diag3(ch3, iind, bra,ket, ch_jk,ket_jk,phase_jk, i,c,b, 2)
           end DO

           ! <bc|r|ijk> <-- + P(bc).<lc|v|ik>.<b|r|jl>
           DO ket = ket_min, ket_max
              j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
              k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
              IF ( j == i .or. k == i ) cycle

              ch_ik = hh_channel_2b%ival2(i,k)
              IF ( ch_ik == 0 ) cycle
              phase_ik = -1
              ket_ik = hh_config_2b%ival2(i,k)
              IF ( k < i ) phase_ik = -phase_ik
              
              CALL r3_diag3(ch3, iind, bra,ket, ch_ik,ket_ik,phase_ik, j,b,c, 1)
              CALL r3_diag3(ch3, iind, bra,ket, ch_ik,ket_ik,phase_ik, j,c,b, 2)
           end DO

           ! <bc|r|ijk> <-- - P(bc).<lc|v|ji>.<b|r|kl>
           DO ket = ket_min, ket_max
              j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
              k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
              IF ( j == i .or. k == i ) cycle

              ch_ji = hh_channel_2b%ival2(j,i)
              IF ( ch_ji == 0 ) cycle
              phase_ji = -1
              ket_ji = hh_config_2b%ival2(j,i)
              IF ( i < j ) phase_ji = -phase_ji

              CALL r3_diag3(ch3, iind, bra,ket, ch_ji,ket_ji,phase_ji, k,b,c, 1)
              CALL r3_diag3(ch3, iind, bra,ket, ch_ji,ket_ji,phase_ji, k,c,b, 2)
           end DO
           
        end DO
     end DO
  end DO

  ! <bc|r|ijk> <-- - P(i/jk).<bc|v|id>.<d|r|jk>
  DO ch3 = ch3_preom_min, ch3_preom_max
     ch1       = ch1_preom3(ch3)
     bra_confs = number_2b(3,ch1)
     DO iind = klimits_preom3(ch3,1), klimits_preom3(ch3,2)
        i       = klist_preom3(ch3)%ival2(iind,1)
        ch2     = klist_preom3(ch3)%ival2(iind,2)
        ket_min = mapping_preom3(ch3)%ival2(iind,1)
        ket_max = mapping_preom3(ch3)%ival2(iind,2)

        ! <bc|r|ijk> <-- - <bc|v|id>.<d|r|jk>
        DO ket = ket_min, ket_max
           j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
           k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
           IF ( j == i .or. k == i ) cycle
           phase_jk = 1
           ket_jk = ket
           ch_jk = ch2
           
           CALL r3_diag4(ch3,ch1, iind, ket, ch_jk,ket_jk,phase_jk, i)
        end DO

        ! <bc|r|ijk> <-- + <bc|v|jd>.<d|r|ik>
        DO ket = ket_min, ket_max
           j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
           k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
           IF ( j == i .or. k == i ) cycle

           ch_ik = hh_channel_2b%ival2(i,k)
           IF ( ch_ik == 0 ) cycle
           phase_ik = -1
           ket_ik = hh_config_2b%ival2(i,k)
           IF ( k < i ) phase_ik = -phase_ik
           
           CALL r3_diag4(ch3,ch1, iind, ket, ch_ik,ket_ik,phase_ik, j)
        end DO

        ! <bc|r|ijk> <-- + <bc|v|kd>.<d|r|ji>
        DO ket = ket_min, ket_max
           j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
           k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
           IF ( j == i .or. k == i ) cycle

           ch_ji = hh_channel_2b%ival2(j,i)
           IF ( ch_ji == 0 ) cycle
           phase_ji = -1
           ket_ji = hh_config_2b%ival2(j,i)
           IF ( i < j ) phase_ji = -phase_ji

           CALL r3_diag4(ch3,ch1, iind, ket, ch_ji,ket_ji,phase_ji, k)
        end DO

     end DO
  end DO
  
end SUBROUTINE preom_2p3h

! <bra|r|(i1)ket> <-- -<lp2|v|ket>.<p1|r|i1l>
SUBROUTINE r3_diag3(ch3, iind, bra,ket, ch2,ket2,phase2, i1,p1,p2, bc_as)
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ch3, iind, bra,ket, ch2,ket2, i1,p1,p2, bc_as, phase2
  INTEGER :: ch, p1ind, bra0,ket0, l, phase_bc

  phase_bc = 1
  IF ( bc_as == 2 ) phase_bc = -1
  
  ! l > i1
  DO l = i1+1, below_ef
     IF ( ch2 /= hp_channel_2b%ival2(l,p2) ) cycle
     bra0 = hp_config_2b%ival2(l,p2)
     ch = hh_channel_2b%ival2(i1,l)
     ket0 = hh_config_2b%ival2(i1,l)
     DO p1ind = 1, 2
        IF ( p1 /= r2_preom_ind(2*(ch-1)+p1ind) ) cycle
        r3_preom(ch3)%val1(iind)%cval(bra,ket) = r3_preom(ch3)%val1(iind)%cval(bra,ket) &
             - phase2 * phase_bc * r2_preom(ch)%cval(p1ind,ket0) * v2b_hphh(ch2)%cval(bra0,ket2)
     end DO
  end DO
  ! l < i1
  DO l = 1, i1-1
     IF ( ch2 /= hp_channel_2b%ival2(l,p2) ) cycle
     bra0 = hp_config_2b%ival2(l,p2)
     ch = hh_channel_2b%ival2(l,i1)
     ket0 = hh_config_2b%ival2(l,i1)
     DO p1ind = 1, 2
        IF ( p1 /= r2_preom_ind(2*(ch-1)+p1ind) ) cycle
        r3_preom(ch3)%val1(iind)%cval(bra,ket) = r3_preom(ch3)%val1(iind)%cval(bra,ket) &
             + phase2 * phase_bc * r2_preom(ch)%cval(p1ind,ket0) * v2b_hphh(ch2)%cval(bra0,ket2)
     end DO
  end DO
  
end SUBROUTINE r3_diag3

! <:|r|(i1)ket> <-- -<:|v|i1d>.<d|r|ket>
SUBROUTINE r3_diag4(ch3,ch1, iind, ket, ch2,ket2,phase2, i1)
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ch3,ch1, iind, ket, ch2,ket2, i1, phase2
  INTEGER :: ket0, d, dind

  DO dind = 1, 2
     d = r2_preom_ind(2*(ch2-1)+dind)
     IF ( d == 0 ) cycle
     IF ( ch1 /= hp_channel_2b%ival2(i1,d) ) cycle
     ket0 = hp_config_2b%ival2(i1,d)
     r3_preom(ch3)%val1(iind)%cval(:,ket) = r3_preom(ch3)%val1(iind)%cval(:,ket) &
          - phase2 * r2_preom(ch2)%cval(dind,ket2) * v2b_pphp(ch1)%cval(:,ket0)
  end DO
  
end SUBROUTINE r3_diag4
