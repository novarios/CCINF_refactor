
MODULE t3_diagrams
  
CONTAINS

  !
  ! T3 diagram 1: <abc|t|ijk> <-- -<cd|t|ij>.<ab|v|kd>
  ! Loops over particle index d, contracts T2 with V2(pphp).
  !
  ! cval is passed explicitly (not accessed via t3_ccm module variable)
  ! so that each OpenMP thread gets its own pointer descriptor on the
  ! stack, avoiding concurrent reads of the shared pointer-remapped
  ! descriptor.
  !
  SUBROUTINE t3_diag1(ch3, ch1,ch2, cind1,kind1, c1,k, bra,bra1, phase1, cval)
    USE single_particle_orbits
    USE configurations
    USE constants
    USE operator_storage
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ch3, ch1,ch2, cind1,kind1, c1,k, bra,bra1, phase1
    COMPLEX(dpc), INTENT(INOUT) :: cval(:,:)
    INTEGER :: bra0,ket0, d, ket,ket1,ket_confs
    
    ket_confs = number_2b_t3(ch3)%ival2(2,ch2)  
    IF ( ket_confs <= 0 ) return
    
    ! d < c1
    DO d = below_ef+1, c1-1
       IF ( ch1 /= hp_channel_2b%ival2(k,d) ) cycle
       ket0 = hp_config_2b%ival2(k,d)
       IF ( ch2 /= pp_channel_2b%ival2(d,c1) ) cycle
       bra0 = pp_config_2b%ival2(d,c1)
       DO ket = 1, ket_confs
          ket1 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
          cval(bra,ket) = cval(bra,ket) &
               + phase1 * t2_ccm(ch2)%cval(bra0,ket1) * v2b_pphp(ch1)%cval(bra1,ket0)
       end DO
    end DO
    ! d > c1
    DO d = c1+1, tot_orbs
       IF ( ch1 /= hp_channel_2b%ival2(k,d) ) cycle
       ket0 = hp_config_2b%ival2(k,d)
       IF ( ch2 /= pp_channel_2b%ival2(c1,d) ) cycle
       bra0 = pp_config_2b%ival2(c1,d)
       DO ket = 1, ket_confs
          ket1 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
          cval(bra,ket) = cval(bra,ket) &
               - phase1 * t2_ccm(ch2)%cval(bra0,ket1) * v2b_pphp(ch1)%cval(bra1,ket0)
       end DO
    end DO
    
  end SUBROUTINE t3_diag1
  
  !
  ! T3 diagram 2: <abc|t|ijk> <-- +<ab|t|lk>.<lc|v|ij>
  ! Loops over hole index l, contracts T2 with V2(hphh).
  !
  SUBROUTINE t3_diag2(ch3, ch1,ch2, cind1,kind1, c1,k, bra,bra1, phase1, cval)
    USE single_particle_orbits
    USE configurations
    USE constants
    USE operator_storage
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ch3, ch1,ch2, cind1,kind1, c1,k, bra,bra1, phase1
    COMPLEX(dpc), INTENT(INOUT) :: cval(:,:)
    INTEGER :: bra0,ket0, l, ket,ket1,ket_confs
    
    ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
    IF ( ket_confs <= 0 ) return
    
    ! l < k
    DO l = 1, k-1
       IF ( ch2 /= hp_channel_2b%ival2(l,c1) ) cycle
       bra0 = hp_config_2b%ival2(l,c1)
       IF ( ch1 /= hh_channel_2b%ival2(l,k) ) cycle
       ket0 = hh_config_2b%ival2(l,k)
       DO ket = 1, ket_confs
          ket1 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
          cval(bra,ket) = cval(bra,ket) &
               + phase1 * t2_ccm(ch1)%cval(bra1,ket0) * v2b_hphh(ch2)%cval(bra0,ket1)
       end DO
    end DO
    ! l > k
    DO l = k+1, below_ef
       IF ( ch2 /= hp_channel_2b%ival2(l,c1) ) cycle
       bra0 = hp_config_2b%ival2(l,c1)
       IF ( ch1 /= hh_channel_2b%ival2(k,l) ) cycle
       ket0 = hh_config_2b%ival2(k,l)
       DO ket = 1, ket_confs
          ket1 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
          cval(bra,ket) = cval(bra,ket) &
               - phase1 * t2_ccm(ch1)%cval(bra1,ket0) * v2b_hphh(ch2)%cval(bra0,ket1)
       end DO
    end DO
    
  end SUBROUTINE t3_diag2  
  
end MODULE t3_diagrams
