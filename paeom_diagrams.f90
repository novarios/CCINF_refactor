
! <a|r|> <-- + <a|x|c>.<c|r|> - (1/2).<la|x|cd>.<cd|r|l>
SUBROUTINE paeom_1p
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: ch, bra,ket, bra_min,bra_max
  INTEGER :: a,c, l,l1, aind,cind1
  INTEGER :: dim1,dim2,dim3
  COMPLEX(dpc) :: r1, v1b, x2
  COMPLEX(dpc), allocatable :: temp_mtx(:,:)

  IF ( iam == 0 ) WRITE(6,*) "...Computing Right-1p amplitudes..."

  ! <a|r|>
  r1_paeom_eqn = 0.d0

  ! <a|r|> <-- - (1/2).<la|v|cd>.<cd|r|l>
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(2,ch) == 0 ) cycle
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     
     IF ( check_my_channel_paeom_pph(ch) == 0 ) cycle
     bra_min = mapping_paeom_pph(iam+1,ch,2)
     bra_max = mapping_paeom_pph(iam+1,ch,3)
     
     ! <la|v|cd>.<cd|r|l>
     dim1 = number_2b(2,ch)
     dim2 = 2
     dim3 = bra_max-bra_min + 1
     ALLOCATE ( temp_mtx(dim1,dim2) )
     temp_mtx = 0.d0
     CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), conjg(v2b_pphp(ch)%cval(bra_min:bra_max,:)), dim3, &
          r2_paeom(ch)%cval(bra_min:bra_max,:), dim3, dcmplx(1.d0,0.d0), temp_mtx, dim1 )
     
     DO bra = 1, dim1
        l   = lookup_2b_configs(2,ch)%ival2(1,bra)
        a   = lookup_2b_configs(2,ch)%ival2(2,bra)
        DO ket = 1, 2
           l1  = r2_paeom_ind(2*(ch-1)+ket)
           IF ( l1 /= l ) cycle
           x2 = temp_mtx(bra,ket)
           IF ( a == r1_paeom_ind(1) ) r1_paeom_eqn(1) = r1_paeom_eqn(1) - x2
           IF ( a == r1_paeom_ind(2) ) r1_paeom_eqn(2) = r1_paeom_eqn(2) - x2
        end DO
     end DO
     DEALLOCATE ( temp_mtx )
  end DO
  CALL mpi_allreduce(mpi_in_place,r1_paeom_eqn,size(r1_paeom_eqn),mpi_complex16,mpi_sum,mpi_comm_world,ierror)

  ! <a|r|> <-- + <a|x|c>.<c|r|>
  DO aind = 1, 2
     a    = r1_paeom_ind(aind)
     DO cind1 = 1, 2
        c    = r1_paeom_ind(cind1)
        v1b  = hbar1b_I2(a,c)
        r1   = r1_paeom(cind1)
        r1_paeom_eqn(aind) = r1_paeom_eqn(aind) + v1b*r1
     end DO
  end DO
  
end SUBROUTINE paeom_1p

! <ab|r|j> <-- - <ab|x|jc>.<c|r|> + P(ab).<a|x|c>.<cb|r|j> - <k|x|j>.<ab|r|k> + (1/2).<ab|x|cd>.<cd|r|j> - P(ab).<ka|x|jc>.<cb|r|k> - <k|z3b|>.<ab|t|kj>
SUBROUTINE paeom_2p1h
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: ch,ch2, bra,ket,bra2,ket2, bra_confs,ket_confs
  INTEGER :: bra_min,bra_max, ket_min,ket_max
  INTEGER :: a,j, k,c, cind1
  INTEGER :: dim1,dim2,dim3
  REAL(dp) :: phase
  COMPLEX(dpc) :: v1b

  IF ( iam == 0 ) WRITE(6,*) "...Computing Right-2p1h amplitudes..."
  
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     ! <ab|r|j>
     r2_paeom_eqn(ch)%cval = 0.d0
  end DO
  
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     
     IF ( check_my_channel_paeom_pph(ch) == 0 ) cycle
     bra_min = mapping_paeom_pph(iam+1,ch,2)
     bra_max = mapping_paeom_pph(iam+1,ch,3)
     bra_confs = number_2b(3,ch)

     ! <ab|r|j> <-- - <k|x|j>.<ab|r|k>
     DO ket  = 1, 2
        j    = r2_paeom_ind(2*(ch-1)+ket)
        DO ket2 = 1, 2
           k    = r2_paeom_ind(2*(ch-1)+ket2)
           v1b  = hbar1b_I3(k,j)
           r2_paeom_eqn(ch)%cval(bra_min:bra_max,ket) = r2_paeom_eqn(ch)%cval(bra_min:bra_max,ket) &
                - r2_paeom(ch)%cval(bra_min:bra_max,ket2)*v1b
        end DO
     end DO

     ! <ab|r|j> <-- + (1/2).<ab|x|cd>.<cd|r|j>
     dim1 = bra_max-bra_min + 1
     dim2 = 2
     dim3 = bra_confs
     CALL ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), hbar2b_I3(ch)%cval(bra_min:bra_max,:), dim1, &
          r2_paeom(ch)%cval, dim3, dcmplx(1.d0,0.d0), r2_paeom_eqn(ch)%cval(bra_min:bra_max,:), dim1 )

     ! <ab|r|j> <-- - <ab|x|jc>.<c|r|>
     DO ket = 1, 2
        j   = r2_paeom_ind(2*(ch-1)+ket)
        DO cind1 = 1, 2
           c    = r1_paeom_ind(cind1)
           IF ( ch /= hp_channel_2b%ival2(j,c) ) cycle
           ket2 = hp_config_2b%ival2(j,c)
           r2_paeom_eqn(ch)%cval(bra_min:bra_max,ket) = r2_paeom_eqn(ch)%cval(bra_min:bra_max,ket) &
                - hbar2b_I6(ch)%cval(bra_min:bra_max,ket2)*r1_paeom(cind1)
        end DO
     end DO

     ! <ab|r|j> <-- - <k|z3b|>.<ab|t|kj>
     DO ket = 1, 2
        j   = r2_paeom_ind(2*(ch-1)+ket)
        DO k = 1, below_ef
           IF ( ch /= hh_channel_2b%ival2(k,j) ) cycle           
           phase = 1
           ket2 = hh_config_2b%ival2(k,j)
           IF ( j < k ) phase = -phase
           r2_paeom_eqn(ch)%cval(bra_min:bra_max,ket) = r2_paeom_eqn(ch)%cval(bra_min:bra_max,ket) &
                - phase * t2_ccm(ch)%cval(bra_min:bra_max,ket2)*hbar3b_paeom_I3(k)
        end DO
     end DO
     
  end DO

  
  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_paeom_ind_cross(2*ch2) == 0 ) cycle
     ! <a|r|j-b>
     r2_paeom_eqn_cross(ch2)%cval = 0.d0
  end DO

  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_paeom_ind_cross(2*ch2) == 0 ) cycle

     IF ( check_my_channel_paeom_php_cross(ch2) == 0 ) cycle
     ket_min = mapping_paeom_php_cross(iam+1,ch2,2)
     ket_max = mapping_paeom_php_cross(iam+1,ch2,3)
     ket_confs = number_2bcross(2,ch2)

     ! <ab|r|j>  <-- + P(ab).<a|x|c>.<cb|r|j>
     ! <a|r|j-b> <-- + P(ab).<a|x|c>.<c|r|j-b>
     DO bra  = 1, 2
        a    = r2_paeom_ind_cross(2*(ch2-1)+bra)
        DO bra2 = 1, 2
           c    = r2_paeom_ind_cross(2*(ch2-1)+bra2)
           v1b  = hbar1b_I2(a,c)
           r2_paeom_eqn_cross(ch2)%cval(bra,ket_min:ket_max) = r2_paeom_eqn_cross(ch2)%cval(bra,ket_min:ket_max) &
                + r2_paeom_cross(ch2)%cval(bra2,ket_min:ket_max)*v1b
        end DO
     end DO

     ! <ab|r|j>  <-- - P(ab).<ka|x|jc>.<cb|r|k>
     ! <ab|r|j>  <-- + P(ab).<kb|x|jc>.<ca|r|k>
     ! <ab|r|j>  <-- - P(ab).<kb|x|jc>.<ac|r|k>
     ! <a|r|j-b> <-- - P(ab).<k-c|x|j-b>.<a|r|k-c>
     ! <a|r|j-b> <-- - P(ab).<a|r|k-c>.<k-c|x|j-b>
     dim1 = 2
     dim2 = ket_max-ket_min + 1
     dim3 = ket_confs
     CALL ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(-1.d0,0.d0), r2_paeom_cross(ch2)%cval, &
          dim1, hbar2b_I5_cross(ch2)%cval(:,ket_min:ket_max), dim3, dcmplx(1.d0,0.d0), &
          r2_paeom_eqn_cross(ch2)%cval(:,ket_min:ket_max), dim1 )
  end DO
  ! Allreduce r2_paeom_eqn in "populate_paeom_vec"

end SUBROUTINE paeom_2p1h



! <a|r|> <-- + (1/4).<kl|v|cd>.<acd|r|kl>
SUBROUTINE paeom_1p_3p2h
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: ch3,ch1,ch2, a,aind,aind0
  INTEGER :: bra_min,bra_max, ket_confs, bra,ket
  COMPLEX(dpc) :: r3, v2b
  COMPLEX(dpc), allocatable :: temp_mtx(:)
  
  IF ( iam == 0 ) WRITE(6,*) "...Computing R1 <- R3..."

  ! <a|r|> <-- + (1/4).<kl|v|cd>.<acd|r|kl>
  ALLOCATE( temp_mtx(2) )
  temp_mtx = 0.d0
  DO ch3 = ch3_paeom_min, ch3_paeom_max
     ch2       = ch2_paeom3(ch3)
     ket_confs = number_2b(1,ch2)
     DO aind = climits_paeom3(ch3,1), climits_paeom3(ch3,2)
        a       = clist_paeom3(ch3)%ival2(aind,1)
        ch1     = clist_paeom3(ch3)%ival2(aind,2)
        IF ( ch1 /= ch2 ) cycle
        IF ( a == r1_paeom_ind(1) ) then
           aind0 = 1
        ELSE IF ( a == r1_paeom_ind(2) ) then
           aind0 = 2
        ELSE
           cycle
        end IF
        bra_min = mapping_paeom3(ch3)%ival2(aind,1)
        bra_max = mapping_paeom3(ch3)%ival2(aind,2)

        DO bra = bra_min, bra_max
           DO ket = 1, ket_confs
              r3  = r3_paeom(ch3)%val1(aind)%cval(bra,ket)
              v2b = v2b_pphh(ch1)%cval(bra,ket)
              temp_mtx(aind0) = temp_mtx(aind0) + conjg(v2b)*r3
           end DO
        end DO
     end DO
  end DO
  
  CALL mpi_allreduce(mpi_in_place,temp_mtx,size(temp_mtx),mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  r1_paeom_eqn(:) = r1_paeom_eqn(:) + temp_mtx(:)
  DEALLOCATE( temp_mtx )
  
end SUBROUTINE paeom_1p_3p2h

! <ab|r|j> <-- + (1/2).P(ab).<kb|v|cd>.<acd|r|kj> - (1/2).<kl|v|jc>.<abc|r|kl>
SUBROUTINE paeom_2p1h_3p2h
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: ch3, ch1,ch2, ch
  INTEGER :: ket_confs, bra_min,bra_max
  INTEGER :: bra,ket, bra0,ket0
  INTEGER :: aind,a, cind1,c, b,j,k
  INTEGER :: dim1,dim2,dim3
  REAL(dp) :: phase, phase0
  COMPLEX(dpc), allocatable :: temp_mtx(:,:)
  
  IF ( iam == 0 ) WRITE(6,*) "...Computing R2 <- R3..."

  ! <ab|r|j> <-- + (1/2).P(ab).<kb|v|cd>.<acd|r|kj>
  DO ch3 = ch3_paeom_min, ch3_paeom_max
     ch2       = ch2_paeom3(ch3)
     ket_confs = number_2b(1,ch2)
     DO aind = climits_paeom3(ch3,1), climits_paeom3(ch3,2)
        a       = clist_paeom3(ch3)%ival2(aind,1)
        ch1     = clist_paeom3(ch3)%ival2(aind,2)
        bra_min = mapping_paeom3(ch3)%ival2(aind,1)
        bra_max = mapping_paeom3(ch3)%ival2(aind,2)
        
        ! <kb|v|cd>.<(a)cd|r|kj>
        dim1 = number_2b(2,ch1)
        dim2 = ket_confs
        dim3 = bra_max-bra_min+1
        IF ( dim1 == 0 ) cycle
        ALLOCATE ( temp_mtx(dim1,dim2) )
        temp_mtx = dcmplx(0.d0,0.d0)
        CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), conjg(v2b_pphp(ch1)%cval(bra_min:bra_max,:)), dim3, &
             r3_paeom(ch3)%val1(aind)%cval(bra_min:bra_max,:), dim3, dcmplx(1.d0,0.d0), temp_mtx, dim1 )
        
        DO bra0 = 1, dim1
           k    = lookup_2b_configs(2,ch1)%ival2(1,bra0)
           b    = lookup_2b_configs(2,ch1)%ival2(2,bra0)

           ch = pp_channel_2b%ival2(a,b)
           IF ( ch == 0 ) cycle
           IF ( number_2b(3,ch) == 0 ) cycle
           IF ( r2_paeom_ind(2*ch) == 0 ) cycle
           phase = 1
           bra = pp_config_2b%ival2(a,b)
           IF ( b < a ) phase = -phase ! P(ab)
           
           DO ket = 1, 2
              j   = r2_paeom_ind(2*(ch-1)+ket)
              IF ( ch2 /= hh_channel_2b%ival2(k,j) ) cycle
              phase0 = 1
              ket0 = hh_config_2b%ival2(k,j)
              IF ( j < k ) phase0 = -phase0
              r2_paeom_eqn(ch)%cval(bra,ket) = r2_paeom_eqn(ch)%cval(bra,ket) + phase * phase0 * temp_mtx(bra0,ket0)
           end DO
        end DO
        DEALLOCATE( temp_mtx )
     end DO
  end DO

  ! <ab|r|j> <-- - (1/2).<kl|v|jc>.<abc|r|kl>
  !              - (1/2).<kl|v|jc>.<cab|r|kl>
  DO ch3 = ch3_paeom_min, ch3_paeom_max
     ch2       = ch2_paeom3(ch3)
     ket_confs = number_2b(1,ch2)
     DO cind1 = climits_paeom3(ch3,1), climits_paeom3(ch3,2)
        c       = clist_paeom3(ch3)%ival2(cind1,1)
        ch1     = clist_paeom3(ch3)%ival2(cind1,2)
        bra_min = mapping_paeom3(ch3)%ival2(cind1,1)
        bra_max = mapping_paeom3(ch3)%ival2(cind1,2)
        IF ( number_2b(3,ch1) == 0 ) cycle
        IF ( r2_paeom_ind(2*ch1) == 0 ) cycle
        
        ! <(c)ab|r|kl>.<kl|v|jc>
        dim1 = bra_max-bra_min+1
        dim2 = number_2b(2,ch2)
        dim3 = ket_confs
        IF ( dim2 == 0 ) cycle
        ALLOCATE( temp_mtx(bra_min:bra_max, dim2) )
        temp_mtx = dcmplx(0.d0,0.d0)
        CALL ZGEMM ( 'n', 't', dim1, dim2, dim3, dcmplx(1.d0,0.d0), r3_paeom(ch3)%val1(cind1)%cval(bra_min:bra_max,:), dim1, &
             conjg(v2b_hphh(ch2)%cval), dim2, dcmplx(1.d0,0.d0), temp_mtx, dim1 )

        DO ket = 1, 2   
           j   = r2_paeom_ind(2*(ch1-1)+ket)
           IF ( ch2 /= hp_channel_2b%ival2(j,c) ) cycle
           ket0 = hp_config_2b%ival2(j,c)
           r2_paeom_eqn(ch1)%cval(bra_min:bra_max,ket) = r2_paeom_eqn(ch1)%cval(bra_min:bra_max,ket) &
                - temp_mtx(bra_min:bra_max,ket0)
        end DO
        DEALLOCATE( temp_mtx )
     end DO
  end DO
  
end SUBROUTINE paeom_2p1h_3p2h

! <abc|r|jk> <-- E(abc_jk).<abc|r|jk> - P(a/bc|jk).<bc|v|kd>.<ad|r|j> - P(a/bc).<la|v|jk>.<bc|r|l>
SUBROUTINE paeom_3p2h
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: ch3, ch1,ch2, bra,ket
  INTEGER :: ket_confs, bra_min,bra_max
  INTEGER :: aind, a,b,c,j,k
  INTEGER :: ch_bc,ch_ac,ch_ba, bra_bc,bra_ac,bra_ba
  INTEGER :: phase_bc, phase_ac, phase_ba
  COMPLEX(dpc) :: denom, r3
  
  IF ( iam == 0 ) WRITE(6,*) "...Computing R3..."

  ! <abc|r|jk> <-- E(abc_jk).<abc|r|jk>
  DO ch3 = ch3_paeom_min, ch3_paeom_max
     ch2       = ch2_paeom3(ch3)
     ket_confs = number_2b(1,ch2)
     DO aind = climits_paeom3(ch3,1), climits_paeom3(ch3,2)
        a       = clist_paeom3(ch3)%ival2(aind,1)
        ch1     = clist_paeom3(ch3)%ival2(aind,2)
        bra_min = mapping_paeom3(ch3)%ival2(aind,1)
        bra_max = mapping_paeom3(ch3)%ival2(aind,2)
        DO bra = bra_min, bra_max
           b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
           c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
           DO ket = 1, ket_confs
              j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
              k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
              
              r3 = r3_paeom(ch3)%val1(aind)%cval(bra,ket)
              r3_paeom(ch3)%val1(aind)%cval(bra,ket) = 0.d0
              denom = ( hbar1b_I2(a,a) + hbar1b_I2(b,b) + hbar1b_I2(c,c) &
                   - hbar1b_I3(j,j) - hbar1b_I3(k,k) )
              r3_paeom(ch3)%val1(aind)%cval(bra,ket) = r3 * denom
           end DO
        end DO
     end DO
  end DO

  ! <abc|r|jk> <-- - P(bc/a|jk).<bc|v|kd>.<ad|r|j>
  DO ch3 = ch3_paeom_min, ch3_paeom_max
     ch2       = ch2_paeom3(ch3)
     ket_confs = number_2b(1,ch2)
     DO aind = climits_paeom3(ch3,1), climits_paeom3(ch3,2)
        a       = clist_paeom3(ch3)%ival2(aind,1)
        ch1     = clist_paeom3(ch3)%ival2(aind,2)
        bra_min = mapping_paeom3(ch3)%ival2(aind,1)
        bra_max = mapping_paeom3(ch3)%ival2(aind,2)

        DO ket = 1, ket_confs
           j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
           k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
           
           ! <abc|r|jk> <-- - P(jk).<bc|v|kd>.<ad|r|j>
           DO bra = bra_min, bra_max
              b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
              c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
              IF ( b == a .or. c == a ) cycle
              phase_bc = 1
              bra_bc = bra
              ch_bc = ch1
              
              CALL r3_diag1(ch3, aind, bra,ket, ch_bc,bra_bc,phase_bc, a,j,k, 1)
              CALL r3_diag1(ch3, aind, bra,ket, ch_bc,bra_bc,phase_bc, a,k,j, 2)
           end DO

           ! <abc|r|jk> <-- + P(jk).<ac|v|kd>.<bd|r|j>
           DO bra = bra_min, bra_max
              b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
              c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
              IF ( b == a .or. c == a ) cycle

              ch_ac = pp_channel_2b%ival2(a,c)
              IF ( ch_ac == 0 ) cycle
              phase_ac = -1
              bra_ac = pp_config_2b%ival2(a,c)
              IF ( c < a ) phase_ac = -phase_ac
              
              CALL r3_diag1(ch3, aind, bra,ket, ch_ac,bra_ac,phase_ac, b,j,k, 1)
              CALL r3_diag1(ch3, aind, bra,ket, ch_ac,bra_ac,phase_ac, b,k,j, 2)
           end DO
        
           ! <abc|r|jk> <-- + P(jk).<ba|v|kd>.<cd|r|j>
           DO bra = bra_min, bra_max
              b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
              c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
              IF ( b == a .or. c == a ) cycle

              ch_ba = pp_channel_2b%ival2(b,a)
              IF ( ch_ba == 0 ) cycle
              phase_ba = -1
              bra_ba = pp_config_2b%ival2(b,a)
              IF ( a < b ) phase_ba = -phase_ba

              CALL r3_diag1(ch3, aind, bra,ket, ch_ba,bra_ba,phase_ba, c,j,k, 1)
              CALL r3_diag1(ch3, aind, bra,ket, ch_ba,bra_ba,phase_ba, c,k,j, 2)
           end DO
           
        end DO
     end DO
  end DO

  ! <abc|r|jk> <-- - P(bc/a).<la|x|jk>.<bc|r|l>
  DO ch3 = ch3_paeom_min, ch3_paeom_max
     ch2       = ch2_paeom3(ch3)
     ket_confs = number_2b(1,ch2)
     DO aind = climits_paeom3(ch3,1), climits_paeom3(ch3,2)
        a       = clist_paeom3(ch3)%ival2(aind,1)
        ch1     = clist_paeom3(ch3)%ival2(aind,2)
        bra_min = mapping_paeom3(ch3)%ival2(aind,1)
        bra_max = mapping_paeom3(ch3)%ival2(aind,2)

        ! <abc|r|jk> <-- - <la|x|jk>.<bc|r|l>
        DO bra = bra_min, bra_max
           b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
           c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
           IF ( b == a .or. c == a ) cycle
           phase_bc = 1
           bra_bc = bra
           ch_bc = ch1
           
           CALL r3_diag2(ch3,ch2, aind, bra, ch_bc,bra_bc,phase_bc, a)
        end DO

        ! <abc|r|jk> <-- + <lb|x|jk>.<ac|r|l>
        DO bra = bra_min, bra_max
           b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
           c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
           IF ( b == a .or. c == a ) cycle

           ch_ac = pp_channel_2b%ival2(a,c)
           IF ( ch_ac == 0 ) cycle
           phase_ac = -1
           bra_ac = pp_config_2b%ival2(a,c)
           IF ( c < a ) phase_ac = -phase_ac
           
           CALL r3_diag2(ch3,ch2, aind, bra, ch_ac,bra_ac,phase_ac, b)
        end DO
           
        ! <abc|r|jk> <-- + <lc|x|jk>.<ba|r|l>
        DO bra = bra_min, bra_max
           b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
           c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
           IF ( b == a .or. c == a ) cycle

           ch_ba = pp_channel_2b%ival2(b,a)
           IF ( ch_ba == 0 ) cycle
           phase_ba = -1
           bra_ba = pp_config_2b%ival2(b,a)
           IF ( a < b ) phase_ba = -phase_ba
           
           CALL r3_diag2(ch3,ch2, aind, bra, ch_ba,bra_ba,phase_ba, c)
        end DO

     end DO
  end DO
  
end SUBROUTINE paeom_3p2h


SUBROUTINE r3_diag1(ch3, aind, bra,ket, ch1,bra1,phase1, a1,h1,h2, jk_as)
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ch3, aind, bra,ket, ch1,bra1, a1,h1,h2, jk_as, phase1
  INTEGER :: ch, h1ind, bra0,ket0, d, phase_jk

  phase_jk = 1
  IF ( jk_as == 2 ) phase_jk = -1

  ! d > a1
  DO d = a1+1, tot_orbs
     IF ( ch1 /= hp_channel_2b%ival2(h2,d) ) cycle
     ket0 = hp_config_2b%ival2(h2,d)
     ch = pp_channel_2b%ival2(a1,d)
     bra0 = pp_config_2b%ival2(a1,d)
     DO h1ind = 1, 2
        IF ( h1 /= r2_paeom_ind(2*(ch-1)+h1ind) ) cycle
        r3_paeom(ch3)%val1(aind)%cval(bra,ket) = r3_paeom(ch3)%val1(aind)%cval(bra,ket) &
             - phase1 * phase_jk * r2_paeom(ch)%cval(bra0,h1ind) * v2b_pphp(ch1)%cval(bra1,ket0)
     end DO
  end DO
  ! d < a1
  DO d = below_ef+1, a1-1
     IF ( ch1 /= hp_channel_2b%ival2(h2,d) ) cycle
     ket0 = hp_config_2b%ival2(h2,d)
     ch = pp_channel_2b%ival2(d,a1)
     bra0 = pp_config_2b%ival2(d,a1)
     DO h1ind = 1, 2
        IF ( h1 /= r2_paeom_ind(2*(ch-1)+h1ind) ) cycle
        r3_paeom(ch3)%val1(aind)%cval(bra,ket) = r3_paeom(ch3)%val1(aind)%cval(bra,ket) &
             + phase1 * phase_jk * r2_paeom(ch)%cval(bra0,h1ind) * v2b_pphp(ch1)%cval(bra1,ket0)
     end DO
  end DO
  
end SUBROUTINE r3_diag1


SUBROUTINE r3_diag2(ch3,ch2, aind, bra, ch1,bra1,phase1, a1)
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ch3,ch2, aind, bra, ch1,bra1, a1, phase1
  INTEGER :: bra0, l, lind
  
  DO lind = 1, 2
     l = r2_paeom_ind(2*(ch1-1)+lind)
     IF ( l == 0 ) cycle
     IF ( ch2 /= hp_channel_2b%ival2(l,a1) ) cycle
     bra0 = hp_config_2b%ival2(l,a1)
     r3_paeom(ch3)%val1(aind)%cval(bra,:) = r3_paeom(ch3)%val1(aind)%cval(bra,:) &
          - phase1 * r2_paeom(ch1)%cval(bra1,lind) * v2b_hphh(ch2)%cval(bra0,:)
  end DO
  
end SUBROUTINE r3_diag2
